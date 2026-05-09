"""Wave 2.5 — VCC-pattern hybrid ablation for neoantigen prediction.

5 tracks:
  T1: PWM (statistical prior) features
  T2: Residual / delta prediction (target = label - PWM_baseline)
  T3: Metric-aware loss (BCE + RankNet)
  T4: Pseudo-bulk per-HLA normalization (post-hoc)
  T5: Combined (PWM-feat + residual + RankNet + pseudo-bulk)

All tracks evaluated on:
  - in-domain 5-fold (training pool n=2396)
  - LOSO mean across 4 train sources
  - ITSNdb_combined (n=311)
  - ITSNdb_no_overlap (n=103)  ← headline
  - VenusVaccine TumorBinary (n=156 protein-level, top1pct_mean aggregator)
  - A*02:01 LOSO (most common allele)

Bootstrap CI: 1000 resamples for each AUROC on no_overlap.

Local CPU only.
"""
from __future__ import annotations
import argparse
import json
import time
import math
from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader, TensorDataset, WeightedRandomSampler
from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss
from sklearn.model_selection import StratifiedKFold

DEVICE = "cpu"
HIDDEN = 640
RNG = 42

AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
AA2I = {a: i for i, a in enumerate(AA_ALPHABET)}
N_AA = len(AA_ALPHABET)


# ---------------------------------------------------------------------------
# Track 1: PWM
# ---------------------------------------------------------------------------
def anchor_positions(peptide: str, max_anchor=11):
    """Return position indices: P1, P2, P3, ..., PΩ-2, PΩ-1, PΩ.
    For length 9: positions = [0,1,2,3,4, -3,-2,-1] mapped to canonical slots.
    We use a fixed 11-slot frame: 0..5 from N-term, 6..10 from C-term (last 5 incl. PΩ).
    """
    L = len(peptide)
    slots = [None] * 11
    # N-terminal anchors (P1..P6)
    for i in range(min(6, L)):
        slots[i] = peptide[i]
    # C-terminal anchors (PΩ-4..PΩ → slots 6..10), only fill if not already filled by N-side
    for j in range(min(5, L)):
        slot = 10 - j
        if slots[slot] is None:
            slots[slot] = peptide[L - 1 - j]
    return slots  # list of 11 chars or None


def build_pwm(train_df, pseudocount=1.0):
    """Build PWM[hla][slot, aa] = log((freq_pos+pc)/(freq_neg+pc)).

    Only HLAs with >=50 examples and pos_rate in (0.05, 0.95) get a PWM.
    Other HLAs get a global PWM (pooled across all sufficient HLAs).
    """
    pwms = {}
    hla_counts = train_df.groupby("HLA_norm")["label"].agg(["count", "mean"]).reset_index()
    eligible = hla_counts[(hla_counts["count"] >= 50) &
                          (hla_counts["mean"] > 0.05) &
                          (hla_counts["mean"] < 0.95)]["HLA_norm"].tolist()
    for hla in eligible:
        sub = train_df[train_df["HLA_norm"] == hla]
        pos_counts = np.zeros((11, N_AA))
        neg_counts = np.zeros((11, N_AA))
        for _, r in sub.iterrows():
            slots = anchor_positions(r["peptide"])
            for s, ch in enumerate(slots):
                if ch is None or ch not in AA2I:
                    continue
                if int(r["label"]) == 1:
                    pos_counts[s, AA2I[ch]] += 1
                else:
                    neg_counts[s, AA2I[ch]] += 1
        # log odds with pseudocount (per-slot normalized)
        pwm = np.log(pos_counts + pseudocount) - np.log(neg_counts + pseudocount)
        # subtract slot mean so PWM is mean-zero per slot (HLA-bias adjusted)
        pwm = pwm - pwm.mean(axis=1, keepdims=True)
        pwms[hla] = pwm

    # Global pooled PWM
    pos_g = np.zeros((11, N_AA))
    neg_g = np.zeros((11, N_AA))
    for _, r in train_df.iterrows():
        slots = anchor_positions(r["peptide"])
        for s, ch in enumerate(slots):
            if ch is None or ch not in AA2I:
                continue
            if int(r["label"]) == 1:
                pos_g[s, AA2I[ch]] += 1
            else:
                neg_g[s, AA2I[ch]] += 1
    pwm_g = np.log(pos_g + pseudocount) - np.log(neg_g + pseudocount)
    pwm_g = pwm_g - pwm_g.mean(axis=1, keepdims=True)
    pwms["__global__"] = pwm_g
    return pwms


def score_pwm(peptide, hla, pwms):
    pwm = pwms.get(hla, pwms["__global__"])
    slots = anchor_positions(peptide)
    s = 0.0
    n = 0
    for slot, ch in enumerate(slots):
        if ch is None or ch not in AA2I:
            continue
        s += pwm[slot, AA2I[ch]]
        n += 1
    return s / max(n, 1)


def pwm_baseline_for(df, pwms):
    """PWM raw score per row, then within-HLA rank-percentile in [0,1]."""
    raw = np.array([score_pwm(p, h, pwms) for p, h in zip(df["peptide"], df["HLA_norm"])])
    return raw


def rank_percentile_within_hla(scores, hlas):
    """Map raw scores → rank percentile in [0,1] within each HLA group."""
    out = np.zeros_like(scores, dtype=np.float64)
    hlas = np.asarray(hlas)
    for h in np.unique(hlas):
        idx = np.where(hlas == h)[0]
        if len(idx) == 1:
            out[idx] = 0.5
            continue
        ranks = np.argsort(np.argsort(scores[idx]))
        out[idx] = ranks / max(len(idx) - 1, 1)
    return out


# ---------------------------------------------------------------------------
# MLP heads
# ---------------------------------------------------------------------------
class MLPHead(nn.Module):
    def __init__(self, in_dim, hid=256, p_drop=0.3):
        super().__init__()
        self.fc1 = nn.Linear(in_dim, hid)
        self.bn1 = nn.BatchNorm1d(hid)
        self.fc2 = nn.Linear(hid, hid)
        self.bn2 = nn.BatchNorm1d(hid)
        self.head = nn.Linear(hid, 1)
        self.p = p_drop

    def forward(self, x):
        h = F.relu(self.bn1(self.fc1(x)))
        h = F.dropout(h, p=self.p, training=self.training)
        h = F.relu(self.bn2(self.fc2(h)))
        h = F.dropout(h, p=self.p, training=self.training)
        return self.head(h).squeeze(-1)


def build_X(df, pep_idx, hla_idx, pep_emb, hla_emb, extra_feats=None):
    n = len(df)
    pi = np.fromiter((pep_idx[p] for p in df["peptide"].tolist()), dtype=np.int64, count=n)
    hi = np.fromiter((hla_idx[h] for h in df["HLA_norm"].tolist()), dtype=np.int64, count=n)
    pe = pep_emb[torch.from_numpy(pi)]
    he = hla_emb[torch.from_numpy(hi)]
    X = torch.cat([pe, he], dim=-1)
    if extra_feats is not None:
        X = torch.cat([X, torch.from_numpy(extra_feats).float()], dim=-1)
    return X


def make_balanced_sampler(df):
    n_per = df["source"].value_counts().to_dict()
    w = df["source"].map(lambda s: 1.0 / n_per[s]).values
    w = w * len(df) / w.sum()
    return WeightedRandomSampler(weights=w.tolist(), num_samples=len(df), replacement=True)


def train_one(model, X, y, df, epochs=15, lr=3e-4, wd=1e-4, batch=256,
              loss_kind="bce", residual_baseline=None, lam_rank=1.0,
              n_pairs_per_batch=128, label_smooth=0.05):
    """Train MLP. loss_kind in {"bce", "residual_bce", "rank", "bce_rank"}.

    For residual_bce: target = label, but model predicts delta; final logit = baseline_logit + delta.
    For rank: BCE + RankNet pairs (pos vs neg in same batch).
    """
    model = model.to(DEVICE).train()
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=wd)
    sampler = make_balanced_sampler(df)
    if residual_baseline is None:
        residual_baseline = np.zeros(len(df), dtype=np.float32)
    base_t = torch.from_numpy(residual_baseline.astype(np.float32))
    y_t = torch.from_numpy(y.astype(np.float32))

    ds = TensorDataset(X, y_t, base_t)
    loader = DataLoader(ds, batch_size=batch, sampler=sampler, drop_last=False)

    for ep in range(epochs):
        model.train()
        for xb, yb, bb in loader:
            xb = xb.to(DEVICE)
            yb = yb.to(DEVICE)
            bb = bb.to(DEVICE)
            yb_s = yb * (1 - label_smooth) + 0.5 * label_smooth
            out = model(xb)  # raw logit "delta" or score
            if loss_kind == "bce":
                logit = out
                loss = F.binary_cross_entropy_with_logits(logit, yb_s)
            elif loss_kind == "residual_bce":
                # baseline is in [0,1] -> convert to logit, add delta, BCE
                bb_logit = torch.logit(bb.clamp(1e-3, 1 - 1e-3))
                logit = bb_logit + out
                loss = F.binary_cross_entropy_with_logits(logit, yb_s)
            elif loss_kind == "rank":
                # pure rank loss: BCE + RankNet
                logit = out
                bce = F.binary_cross_entropy_with_logits(logit, yb_s)
                pos_idx = (yb > 0.5).nonzero(as_tuple=True)[0]
                neg_idx = (yb < 0.5).nonzero(as_tuple=True)[0]
                if len(pos_idx) > 0 and len(neg_idx) > 0:
                    np_pairs = min(n_pairs_per_batch, len(pos_idx) * len(neg_idx))
                    pi = pos_idx[torch.randint(0, len(pos_idx), (np_pairs,))]
                    ni = neg_idx[torch.randint(0, len(neg_idx), (np_pairs,))]
                    rank_loss = -F.logsigmoid(logit[pi] - logit[ni]).mean()
                else:
                    rank_loss = torch.tensor(0.0, device=DEVICE)
                loss = bce + lam_rank * rank_loss
            elif loss_kind == "vcc_hybrid":
                # residual + rank
                bb_logit = torch.logit(bb.clamp(1e-3, 1 - 1e-3))
                logit = bb_logit + out
                bce = F.binary_cross_entropy_with_logits(logit, yb_s)
                pos_idx = (yb > 0.5).nonzero(as_tuple=True)[0]
                neg_idx = (yb < 0.5).nonzero(as_tuple=True)[0]
                if len(pos_idx) > 0 and len(neg_idx) > 0:
                    np_pairs = min(n_pairs_per_batch, len(pos_idx) * len(neg_idx))
                    pi = pos_idx[torch.randint(0, len(pos_idx), (np_pairs,))]
                    ni = neg_idx[torch.randint(0, len(neg_idx), (np_pairs,))]
                    rank_loss = -F.logsigmoid(logit[pi] - logit[ni]).mean()
                else:
                    rank_loss = torch.tensor(0.0, device=DEVICE)
                loss = bce + lam_rank * rank_loss
            else:
                raise ValueError(loss_kind)
            opt.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            opt.step()
    return model


@torch.no_grad()
def predict(model, X, baseline=None, loss_kind="bce", batch=2048):
    model.eval()
    n = X.size(0)
    out_all = []
    for i in range(0, n, batch):
        out_all.append(model(X[i:i+batch]).cpu())
    out = torch.cat(out_all).numpy()
    if loss_kind in ("residual_bce", "vcc_hybrid"):
        if baseline is None:
            raise ValueError("need baseline")
        bb = np.clip(baseline.astype(np.float32), 1e-3, 1 - 1e-3)
        bb_logit = np.log(bb / (1 - bb))
        return 1.0 / (1.0 + np.exp(-(bb_logit + out)))
    else:
        return 1.0 / (1.0 + np.exp(-out))


# ---------------------------------------------------------------------------
# Pseudo-bulk
# ---------------------------------------------------------------------------
def pseudo_bulk(scores, hlas):
    """Subtract per-HLA mean from scores. Returns adjusted scores in same shape."""
    scores = np.asarray(scores)
    hlas = np.asarray(hlas)
    out = scores.copy().astype(np.float64)
    for h in np.unique(hlas):
        m = hlas == h
        out[m] = scores[m] - scores[m].mean()
    return out


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------
def safe_auc(y, p):
    if len(set(y)) < 2:
        return float("nan")
    return float(roc_auc_score(y, p))


def ece(y, p, n_bins=15):
    y, p = np.asarray(y), np.asarray(p)
    bins = np.linspace(0, 1, n_bins + 1)
    val = 0.0
    n = len(y)
    for lo, hi in zip(bins[:-1], bins[1:]):
        m = (p >= lo) & (p < hi if hi < 1 else p <= hi)
        if m.sum() == 0:
            continue
        val += (m.sum() / n) * abs(p[m].mean() - y[m].mean())
    return float(val)


def bootstrap_auc_ci(y, p, n_boot=1000, seed=42):
    y = np.asarray(y); p = np.asarray(p)
    if len(set(y)) < 2:
        return (float("nan"), float("nan"))
    rng = np.random.default_rng(seed)
    aucs = []
    n = len(y)
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        if len(set(y[idx])) < 2:
            continue
        aucs.append(roc_auc_score(y[idx], p[idx]))
    if len(aucs) < 50:
        return (float("nan"), float("nan"))
    lo, hi = np.percentile(aucs, [2.5, 97.5])
    return (float(lo), float(hi))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bundle", default="project/results/p_neo_bayesian_2026_05_09/bundle.tsv")
    ap.add_argument("--embeddings", default="project/results/p_neo_bayesian_2026_05_09/embeddings.pt")
    ap.add_argument("--out_dir", default="project/results/p_neo_bayesian_2026_05_09/wave25")
    ap.add_argument("--epochs", type=int, default=15)
    ap.add_argument("--n_seeds", type=int, default=3)
    ap.add_argument("--lam_rank", type=float, default=1.0)
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    log_path = out / "wave25_run.log"
    log_f = open(log_path, "w")

    def log(msg):
        print(msg)
        log_f.write(msg + "\n")
        log_f.flush()

    t0 = time.time()
    log("[wave25] loading bundle + embeddings…")
    bundle = pd.read_csv(args.bundle, sep="\t")
    bundle["HLA_norm"] = bundle["HLA_norm"].fillna("")
    emb = torch.load(args.embeddings, map_location="cpu", weights_only=False)
    pep_idx = {k: i for i, k in enumerate(emb["pep_keys"])}
    hla_idx = {k: i for i, k in enumerate(emb["hla_keys"])}
    pep_emb = emb["pep_emb"]
    hla_emb = emb["hla_emb"]
    log(f"  pep emb: {pep_emb.shape}, hla emb: {hla_emb.shape}")

    has_emb = lambda p, h: p in pep_idx and h in hla_idx
    nv = bundle[~bundle["split"].isin(["ext_venus_test", "ext_venus_valid"])].copy()
    nv = nv[[has_emb(p, h) for p, h in zip(nv["peptide"], nv["HLA_norm"])]].copy()
    train_df = nv[nv["split"] == "train"].copy().reset_index(drop=True)
    log(f"  train n={len(train_df)}")
    itsn = nv[nv["split"].isin(["ext_itsndb_main", "ext_itsndb_val"])].copy().reset_index(drop=True)
    log(f"  ITSNdb n={len(itsn)}")

    # in_master flag for no_overlap
    full_meta = bundle.set_index(["peptide", "HLA_norm", "split"])["in_master"].to_dict()
    itsn["in_master"] = itsn.apply(
        lambda r: bool(full_meta.get((r["peptide"], r["HLA_norm"], r["split"]), False)), axis=1)
    log(f"  ITSNdb no_overlap n={(~itsn['in_master']).sum()}")

    # =========================================================================
    # TRACK 1: PWM
    # =========================================================================
    log("\n[T1] Building PWMs…")
    pwms = build_pwm(train_df, pseudocount=1.0)
    log(f"  per-HLA PWMs: {len(pwms) - 1} (+ global)")
    train_pwm_raw = pwm_baseline_for(train_df, pwms)
    itsn_pwm_raw = pwm_baseline_for(itsn, pwms)

    # rank-percentile within HLA → [0,1] sigmoid-like baseline
    train_pwm_pct = rank_percentile_within_hla(train_pwm_raw, train_df["HLA_norm"].values)
    itsn_pwm_pct = rank_percentile_within_hla(itsn_pwm_raw, itsn["HLA_norm"].values)

    # PWM-only AUROC on no_overlap (sanity)
    no_mask = ~itsn["in_master"].values
    auc_pwm_no = safe_auc(itsn["label"].values[no_mask], itsn_pwm_pct[no_mask])
    auc_pwm_comb = safe_auc(itsn["label"].values, itsn_pwm_pct)
    log(f"  PWM-alone ITSNdb_combined AUROC = {auc_pwm_comb:.4f}")
    log(f"  PWM-alone ITSNdb_no_overlap AUROC = {auc_pwm_no:.4f}")

    # save pwm features
    pd.DataFrame({
        "peptide": list(train_df["peptide"]) + list(itsn["peptide"]),
        "hla":     list(train_df["HLA_norm"]) + list(itsn["HLA_norm"]),
        "split":   ["train"] * len(train_df) + ["itsndb"] * len(itsn),
        "label":   list(train_df["label"]) + list(itsn["label"]),
        "pwm_raw": list(train_pwm_raw) + list(itsn_pwm_raw),
        "pwm_pct": list(train_pwm_pct) + list(itsn_pwm_pct),
    }).to_csv(out / "pwm_features.tsv", sep="\t", index=False)

    # Build X tensors once
    X_tr = build_X(train_df, pep_idx, hla_idx, pep_emb, hla_emb)
    X_itsn = build_X(itsn, pep_idx, hla_idx, pep_emb, hla_emb)
    # Extra-feat versions (for T5)
    train_extra = np.stack([train_pwm_raw, train_pwm_pct], axis=1).astype(np.float32)
    itsn_extra = np.stack([itsn_pwm_raw, itsn_pwm_pct], axis=1).astype(np.float32)
    X_tr_ext = build_X(train_df, pep_idx, hla_idx, pep_emb, hla_emb, extra_feats=train_extra)
    X_itsn_ext = build_X(itsn, pep_idx, hla_idx, pep_emb, hla_emb, extra_feats=itsn_extra)

    y_tr = train_df["label"].astype(int).values
    y_itsn = itsn["label"].astype(int).values

    src2id = {s: i for i, s in enumerate(sorted(train_df["source"].unique()))}

    # =========================================================================
    # Helper: train ensemble on a config, return ITSNdb predictions
    # =========================================================================
    def train_ens_and_predict(loss_kind, in_dim, X_train_full, X_eval, baseline_train, baseline_eval,
                              n_seeds=3, label="?"):
        """Train n_seeds MLPs on full train pool, average predictions."""
        ens = []
        for seed in range(n_seeds):
            torch.manual_seed(seed); np.random.seed(seed)
            model = MLPHead(in_dim=in_dim, hid=256, p_drop=0.3)
            train_one(model, X_train_full, y_tr, train_df,
                      epochs=args.epochs, batch=256,
                      loss_kind=loss_kind, residual_baseline=baseline_train,
                      lam_rank=args.lam_rank)
            p = predict(model, X_eval, baseline=baseline_eval, loss_kind=loss_kind)
            ens.append(p)
        return np.stack(ens).mean(0)

    # =========================================================================
    # In-domain 5-fold for each track
    # =========================================================================
    def kfold_eval(loss_kind, in_dim, X_full, baseline_full, n_seeds=2):
        """5-fold CV mean AUROC. Returns per-fold AUROCs + averaged predictions."""
        skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
        aucs = []
        for fold, (tr_i, te_i) in enumerate(skf.split(np.arange(len(train_df)), y_tr)):
            tr_sub = train_df.iloc[tr_i].reset_index(drop=True)
            te_sub_y = y_tr[te_i]
            X_tr_sub = X_full[tr_i]
            X_te_sub = X_full[te_i]
            base_tr = baseline_full[tr_i] if baseline_full is not None else None
            base_te = baseline_full[te_i] if baseline_full is not None else None
            ens = []
            for seed in range(n_seeds):
                torch.manual_seed(seed); np.random.seed(seed)
                model = MLPHead(in_dim=in_dim, hid=256, p_drop=0.3)
                tr_y_sub = y_tr[tr_i]
                train_one(model, X_tr_sub, tr_y_sub, tr_sub,
                          epochs=args.epochs, batch=256,
                          loss_kind=loss_kind, residual_baseline=base_tr,
                          lam_rank=args.lam_rank)
                p = predict(model, X_te_sub, baseline=base_te, loss_kind=loss_kind)
                ens.append(p)
            ens_p = np.stack(ens).mean(0)
            aucs.append(safe_auc(te_sub_y, ens_p))
            log(f"    fold{fold} {loss_kind} AUC={aucs[-1]:.4f}")
        return float(np.mean(aucs)), aucs

    # =========================================================================
    # LOSO across train sources
    # =========================================================================
    def loso_eval(loss_kind, in_dim, X_full, baseline_full, n_seeds=1):
        sources = sorted(train_df["source"].unique())
        aucs = []
        for src in sources:
            tr_mask = train_df["source"].values != src
            te_mask = ~tr_mask
            if y_tr[te_mask].sum() < 1 or y_tr[te_mask].sum() == te_mask.sum():
                continue
            tr_idx = np.where(tr_mask)[0]
            te_idx = np.where(te_mask)[0]
            tr_df_sub = train_df.iloc[tr_idx].reset_index(drop=True)
            X_tr_sub = X_full[tr_idx]
            X_te_sub = X_full[te_idx]
            base_tr = baseline_full[tr_idx] if baseline_full is not None else None
            base_te = baseline_full[te_idx] if baseline_full is not None else None
            ens = []
            for seed in range(n_seeds):
                torch.manual_seed(seed); np.random.seed(seed)
                model = MLPHead(in_dim=in_dim, hid=256, p_drop=0.3)
                train_one(model, X_tr_sub, y_tr[tr_idx], tr_df_sub,
                          epochs=args.epochs, batch=256,
                          loss_kind=loss_kind, residual_baseline=base_tr,
                          lam_rank=args.lam_rank)
                p = predict(model, X_te_sub, baseline=base_te, loss_kind=loss_kind)
                ens.append(p)
            ens_p = np.stack(ens).mean(0)
            auc = safe_auc(y_tr[te_idx], ens_p)
            aucs.append(auc)
            log(f"    LOSO {src:>22s} {loss_kind} AUC={auc:.4f}")
        return float(np.mean(aucs)), aucs

    # =========================================================================
    # A*02:01 LOSO
    # =========================================================================
    def a0201_loso(loss_kind, in_dim, X_full, baseline_full, n_seeds=1):
        a = "HLA-A*02:01"
        tr_mask = train_df["HLA_norm"].values != a
        te_mask = ~tr_mask
        tr_idx = np.where(tr_mask)[0]
        te_idx = np.where(te_mask)[0]
        if len(te_idx) < 30 or len(set(y_tr[te_idx])) < 2:
            return float("nan")
        tr_df_sub = train_df.iloc[tr_idx].reset_index(drop=True)
        X_tr_sub = X_full[tr_idx]
        X_te_sub = X_full[te_idx]
        base_tr = baseline_full[tr_idx] if baseline_full is not None else None
        base_te = baseline_full[te_idx] if baseline_full is not None else None
        ens = []
        for seed in range(n_seeds):
            torch.manual_seed(seed); np.random.seed(seed)
            model = MLPHead(in_dim=in_dim, hid=256, p_drop=0.3)
            train_one(model, X_tr_sub, y_tr[tr_idx], tr_df_sub,
                      epochs=args.epochs, batch=256,
                      loss_kind=loss_kind, residual_baseline=base_tr,
                      lam_rank=args.lam_rank)
            p = predict(model, X_te_sub, baseline=base_te, loss_kind=loss_kind)
            ens.append(p)
        ens_p = np.stack(ens).mean(0)
        return safe_auc(y_tr[te_idx], ens_p)

    # =========================================================================
    # Run each track
    # =========================================================================
    results = []  # rows for results table
    pred_records = {}  # method -> array of predictions on itsn

    # ---- Wave1 baseline (BCE, ESM2-only, 1280d) -----------------------------
    log("\n[Wave1-style baseline (BCE on ESM2 1280d)]")
    auc_kf_w1, _ = kfold_eval("bce", 2 * HIDDEN, X_tr, baseline_full=None, n_seeds=2)
    auc_loso_w1, _ = loso_eval("bce", 2 * HIDDEN, X_tr, None, n_seeds=1)
    auc_a02_w1 = a0201_loso("bce", 2 * HIDDEN, X_tr, None, n_seeds=1)
    p_itsn_w1 = train_ens_and_predict("bce", 2 * HIDDEN, X_tr, X_itsn, None, None,
                                       n_seeds=args.n_seeds, label="W1")
    pred_records["W1_baseline"] = p_itsn_w1

    # ---- Track 1: PWM concat (ESM2 + 2 PWM scalars) -----------------------
    log("\n[T1: ESM2 + PWM features]")
    in_dim_t1 = 2 * HIDDEN + 2
    auc_kf_t1, _ = kfold_eval("bce", in_dim_t1, X_tr_ext, None, n_seeds=2)
    auc_loso_t1, _ = loso_eval("bce", in_dim_t1, X_tr_ext, None, n_seeds=1)
    auc_a02_t1 = a0201_loso("bce", in_dim_t1, X_tr_ext, None, n_seeds=1)
    p_itsn_t1 = train_ens_and_predict("bce", in_dim_t1, X_tr_ext, X_itsn_ext, None, None,
                                       n_seeds=args.n_seeds)
    pred_records["T1_pwm_concat"] = p_itsn_t1

    # ---- Track 2: residual prediction --------------------------------------
    log("\n[T2: residual prediction (target = label - PWM_baseline)]")
    auc_kf_t2, _ = kfold_eval("residual_bce", 2 * HIDDEN, X_tr, train_pwm_pct, n_seeds=2)
    auc_loso_t2, _ = loso_eval("residual_bce", 2 * HIDDEN, X_tr, train_pwm_pct, n_seeds=1)
    auc_a02_t2 = a0201_loso("residual_bce", 2 * HIDDEN, X_tr, train_pwm_pct, n_seeds=1)
    p_itsn_t2 = train_ens_and_predict("residual_bce", 2 * HIDDEN, X_tr, X_itsn,
                                       train_pwm_pct, itsn_pwm_pct,
                                       n_seeds=args.n_seeds)
    pred_records["T2_residual"] = p_itsn_t2

    # ---- Track 3: BCE + RankNet --------------------------------------------
    log("\n[T3: BCE + RankNet]")
    auc_kf_t3, _ = kfold_eval("rank", 2 * HIDDEN, X_tr, None, n_seeds=2)
    auc_loso_t3, _ = loso_eval("rank", 2 * HIDDEN, X_tr, None, n_seeds=1)
    auc_a02_t3 = a0201_loso("rank", 2 * HIDDEN, X_tr, None, n_seeds=1)
    p_itsn_t3 = train_ens_and_predict("rank", 2 * HIDDEN, X_tr, X_itsn, None, None,
                                       n_seeds=args.n_seeds)
    pred_records["T3_rank"] = p_itsn_t3

    # ---- Track 4: pseudo-bulk on Wave1-baseline preds -----------------------
    log("\n[T4: pseudo-bulk post-hoc]")
    p_itsn_t4_w1 = pseudo_bulk(p_itsn_w1, itsn["HLA_norm"].values)
    # also apply on T2/T3 preds for combined assessment
    p_itsn_t4_t2 = pseudo_bulk(p_itsn_t2, itsn["HLA_norm"].values)
    p_itsn_t4_t3 = pseudo_bulk(p_itsn_t3, itsn["HLA_norm"].values)
    pred_records["T4_pb_W1"] = p_itsn_t4_w1
    pred_records["T4_pb_T2"] = p_itsn_t4_t2
    pred_records["T4_pb_T3"] = p_itsn_t4_t3

    # ---- Track 5: VCC hybrid -----------------------------------------------
    log("\n[T5: VCC hybrid (PWM-feat + residual + RankNet)]")
    in_dim_t5 = 2 * HIDDEN + 2
    auc_kf_t5, _ = kfold_eval("vcc_hybrid", in_dim_t5, X_tr_ext, train_pwm_pct, n_seeds=2)
    auc_loso_t5, _ = loso_eval("vcc_hybrid", in_dim_t5, X_tr_ext, train_pwm_pct, n_seeds=1)
    auc_a02_t5 = a0201_loso("vcc_hybrid", in_dim_t5, X_tr_ext, train_pwm_pct, n_seeds=1)
    p_itsn_t5 = train_ens_and_predict("vcc_hybrid", in_dim_t5, X_tr_ext, X_itsn_ext,
                                       train_pwm_pct, itsn_pwm_pct,
                                       n_seeds=args.n_seeds)
    # apply pseudo-bulk too
    p_itsn_t5_pb = pseudo_bulk(p_itsn_t5, itsn["HLA_norm"].values)
    pred_records["T5_vcc_hybrid"] = p_itsn_t5
    pred_records["T5_vcc_hybrid_pb"] = p_itsn_t5_pb

    # =========================================================================
    # Compute ITSNdb metrics for all methods
    # =========================================================================
    def itsndb_metrics(p):
        rows = []
        for sub_name, mask in [
            ("ITSNdb_main", itsn["split"].values == "ext_itsndb_main"),
            ("ITSNdb_combined", np.ones(len(itsn), dtype=bool)),
            ("ITSNdb_no_overlap", ~itsn["in_master"].values),
        ]:
            yy = y_itsn[mask]; pp = p[mask]
            if len(set(yy)) < 2:
                continue
            auc = safe_auc(yy, pp)
            ap = float(average_precision_score(yy, pp))
            # Brier/ECE only meaningful if scores in [0,1]
            if pp.min() >= 0 and pp.max() <= 1:
                br = float(brier_score_loss(yy, pp))
                ec = ece(yy, pp)
            else:
                br = float("nan")
                ec = float("nan")
            ci_lo, ci_hi = bootstrap_auc_ci(yy, pp, n_boot=1000, seed=42)
            rows.append({
                "subset": sub_name, "n": int(len(yy)), "n_pos": int(yy.sum()),
                "AUROC": auc, "AUPRC": ap, "Brier": br, "ECE": ec,
                "AUROC_lo95": ci_lo, "AUROC_hi95": ci_hi,
            })
        return rows

    method_results = {}
    for name, p in pred_records.items():
        method_results[name] = itsndb_metrics(p)

    # PWM-alone ITSNdb metrics
    method_results["PWM_alone"] = itsndb_metrics(itsn_pwm_pct)

    # =========================================================================
    # Aggregate into wave25_results.tsv
    # =========================================================================
    log("\n=== Wave 2.5 results ===")
    rows = []
    method_kf = {
        "W1_baseline": auc_kf_w1, "T1_pwm_concat": auc_kf_t1,
        "T2_residual": auc_kf_t2, "T3_rank": auc_kf_t3,
        "T5_vcc_hybrid": auc_kf_t5,
    }
    method_loso = {
        "W1_baseline": auc_loso_w1, "T1_pwm_concat": auc_loso_t1,
        "T2_residual": auc_loso_t2, "T3_rank": auc_loso_t3,
        "T5_vcc_hybrid": auc_loso_t5,
    }
    method_a02 = {
        "W1_baseline": auc_a02_w1, "T1_pwm_concat": auc_a02_t1,
        "T2_residual": auc_a02_t2, "T3_rank": auc_a02_t3,
        "T5_vcc_hybrid": auc_a02_t5,
    }

    for method, mlist in method_results.items():
        # extract metrics by subset
        by_subset = {r["subset"]: r for r in mlist}
        comb = by_subset.get("ITSNdb_combined", {})
        no = by_subset.get("ITSNdb_no_overlap", {})
        rows.append({
            "Method": method,
            "in_domain_5fold": method_kf.get(method, float("nan")),
            "LOSO_mean": method_loso.get(method, float("nan")),
            "ITSNdb_combined": comb.get("AUROC", float("nan")),
            "ITSNdb_combined_n": comb.get("n", 0),
            "ITSNdb_no_overlap": no.get("AUROC", float("nan")),
            "ITSNdb_no_overlap_n": no.get("n", 0),
            "ITSNdb_no_overlap_lo95": no.get("AUROC_lo95", float("nan")),
            "ITSNdb_no_overlap_hi95": no.get("AUROC_hi95", float("nan")),
            "A0201_LOSO": method_a02.get(method, float("nan")),
            "ECE_combined": comb.get("ECE", float("nan")),
        })
    res_df = pd.DataFrame(rows)
    res_df.to_csv(out / "wave25_results.tsv", sep="\t", index=False)
    log(res_df.to_string(index=False))

    # save predictions
    pred_df = pd.DataFrame({
        "peptide": itsn["peptide"], "hla": itsn["HLA_norm"], "label": y_itsn,
        "in_master": itsn["in_master"].values, "split": itsn["split"].values,
    })
    for name, p in pred_records.items():
        pred_df[f"pred_{name}"] = p
    pred_df["pred_PWM_alone"] = itsn_pwm_pct
    pred_df.to_csv(out / "predictions_all_tracks.tsv", sep="\t", index=False)
    pd.DataFrame({"peptide": itsn["peptide"], "hla": itsn["HLA_norm"], "label": y_itsn,
                  "score": pred_records["T2_residual"]}).to_csv(
        out / "predictions_residual.tsv", sep="\t", index=False)
    pd.DataFrame({"peptide": itsn["peptide"], "hla": itsn["HLA_norm"], "label": y_itsn,
                  "score": pred_records["T3_rank"]}).to_csv(
        out / "predictions_rank.tsv", sep="\t", index=False)
    pd.DataFrame({"peptide": itsn["peptide"], "hla": itsn["HLA_norm"], "label": y_itsn,
                  "score": pred_records["T5_vcc_hybrid"]}).to_csv(
        out / "predictions_vcc_hybrid.tsv", sep="\t", index=False)

    # =========================================================================
    # Figure: VCC ablation bar chart
    # =========================================================================
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(10, 5))
        order = ["PWM_alone", "W1_baseline", "T1_pwm_concat", "T2_residual", "T3_rank",
                 "T4_pb_W1", "T4_pb_T2", "T4_pb_T3", "T5_vcc_hybrid", "T5_vcc_hybrid_pb"]
        no_aucs = []
        ci_los = []
        ci_his = []
        for m in order:
            ml = method_results.get(m, [])
            no = next((r for r in ml if r["subset"] == "ITSNdb_no_overlap"), None)
            no_aucs.append(no["AUROC"] if no else 0.5)
            ci_los.append(no["AUROC_lo95"] if no else 0.5)
            ci_his.append(no["AUROC_hi95"] if no else 0.5)
        x = np.arange(len(order))
        # baseline = W1
        baseline = no_aucs[order.index("W1_baseline")]
        colors = ["gray" if m == "PWM_alone" else
                  "k" if m == "W1_baseline" else
                  "C2" if a > baseline else "C3" for m, a in zip(order, no_aucs)]
        err_lo = [a - lo for a, lo in zip(no_aucs, ci_los)]
        err_hi = [hi - a for a, hi in zip(no_aucs, ci_his)]
        ax.bar(x, no_aucs, yerr=[err_lo, err_hi], color=colors, capsize=4)
        ax.axhline(0.5, color="k", linestyle=":", alpha=0.5, label="chance")
        ax.axhline(baseline, color="C0", linestyle="--", alpha=0.5, label=f"W1 baseline {baseline:.3f}")
        ax.set_xticks(x)
        ax.set_xticklabels(order, rotation=30, ha="right")
        ax.set_ylabel("AUROC on ITSNdb_no_overlap (n=103)")
        ax.set_title("Wave 2.5 VCC-pattern ablation (95% bootstrap CI)")
        ax.set_ylim(0.2, max(no_aucs + [0.7]) + 0.05)
        ax.legend(loc="best")
        plt.tight_layout()
        plt.savefig(out / "fig_vcc_ablation.png", dpi=150)
        plt.savefig(out / "fig_vcc_ablation.pdf")
        plt.close()
        log(f"  saved fig_vcc_ablation.{{png,pdf}}")
    except Exception as e:
        log(f"  fig failed: {e}")

    log(f"\n[wave25] done in {time.time()-t0:.1f}s — outputs at {out}")
    log_f.close()


if __name__ == "__main__":
    main()
