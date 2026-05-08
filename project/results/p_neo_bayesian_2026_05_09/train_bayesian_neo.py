"""Bayesian MLP head on frozen ESM2-150M (peptide + HLA pseudo-seq) embeddings.

Training tricks:
  - Source-balanced sampler (weight = 1/n_source per row)
  - Mixup α=0.2 on peptide+HLA embeddings + soft labels
  - Label smoothing ε=0.05
  - Focal loss γ=2
  - DANN domain-adversarial branch (gradient reverse) — toggle via --no-dann
  - MC Dropout T=30 at inference
  - Deep Ensemble n_seeds=5 (per --seeds)

All eval is run inside this script so we don't move embeddings around twice.

Outputs (all under --out_dir):
  results_summary.tsv         — AUROC + ECE for each (eval-split × model variant)
  per_allele_loso_bayesian.tsv
  predictions_*.tsv
  fig_*.png/pdf
  BAYESIAN_NEO_REPORT.md
"""
from __future__ import annotations
import argparse, json, time, math, os
from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader, WeightedRandomSampler
from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss
from sklearn.model_selection import StratifiedKFold

DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
HIDDEN = 640


# ---------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------
class GradReverse(torch.autograd.Function):
    @staticmethod
    def forward(ctx, x, lam):
        ctx.lam = lam
        return x.view_as(x)

    @staticmethod
    def backward(ctx, g):
        return -ctx.lam * g, None


def grad_reverse(x, lam=1.0):
    return GradReverse.apply(x, lam)


class BayesianMLP(nn.Module):
    """[pep_640 + hla_640] -> 256 -> 1, dropout p=0.3 on every layer."""
    def __init__(self, in_dim=1280, hid=256, n_sources=4, p_drop=0.3, dann=True):
        super().__init__()
        self.dann = dann
        self.fc1 = nn.Linear(in_dim, hid)
        self.bn1 = nn.BatchNorm1d(hid)
        self.fc2 = nn.Linear(hid, hid)
        self.bn2 = nn.BatchNorm1d(hid)
        self.head = nn.Linear(hid, 1)
        self.dom = nn.Sequential(
            nn.Linear(hid, hid), nn.ReLU(), nn.Dropout(p_drop), nn.Linear(hid, n_sources)
        ) if dann else None
        self.p = p_drop

    def features(self, x):
        h = F.relu(self.bn1(self.fc1(x)))
        h = F.dropout(h, p=self.p, training=True)  # always-on for MC Dropout
        h = F.relu(self.bn2(self.fc2(h)))
        h = F.dropout(h, p=self.p, training=True)
        return h

    def forward(self, x, lam_dann=1.0):
        h = self.features(x)
        logit = self.head(h).squeeze(-1)
        if self.dom is not None:
            dom = self.dom(grad_reverse(h, lam_dann))
            return logit, dom
        return logit, None


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------
class NeoDataset(Dataset):
    def __init__(self, df, pep_idx, hla_idx, pep_emb, hla_emb, src2id):
        self.peps = df["peptide"].tolist()
        self.hlas = df["HLA_norm"].tolist()
        self.labels = df["label"].astype(int).values
        self.sources = df["source"].tolist()
        self.pep_idx = pep_idx
        self.hla_idx = hla_idx
        self.pep_emb = pep_emb
        self.hla_emb = hla_emb
        self.src2id = src2id

    def __len__(self):
        return len(self.peps)

    def __getitem__(self, i):
        pe = self.pep_emb[self.pep_idx[self.peps[i]]]
        he = self.hla_emb[self.hla_idx[self.hlas[i]]]
        x = torch.cat([pe, he], dim=-1)
        return x, float(self.labels[i]), self.src2id.get(self.sources[i], 0)


def make_loader(df, pep_idx, hla_idx, pep_emb, hla_emb, src2id, batch=256, balanced=True):
    """Vectorize: pre-build X, y, src tensors once then use TensorDataset."""
    n = len(df)
    pi = np.fromiter((pep_idx[p] for p in df["peptide"].tolist()), dtype=np.int64, count=n)
    hi = np.fromiter((hla_idx[h] for h in df["HLA_norm"].tolist()), dtype=np.int64, count=n)
    pe = pep_emb[torch.from_numpy(pi)]
    he = hla_emb[torch.from_numpy(hi)]
    X = torch.cat([pe, he], dim=-1).contiguous()
    y = torch.tensor(df["label"].astype(int).values, dtype=torch.float32)
    s = torch.tensor([src2id.get(s, 0) for s in df["source"].tolist()], dtype=torch.long)
    ds = torch.utils.data.TensorDataset(X, y, s)
    if balanced:
        n_per_src = df["source"].value_counts().to_dict()
        w = df["source"].map(lambda s: 1.0 / n_per_src[s]).values
        w = w * len(df) / w.sum()
        sampler = WeightedRandomSampler(weights=w.tolist(), num_samples=len(df), replacement=True)
        loader = DataLoader(ds, batch_size=batch, sampler=sampler, num_workers=0, drop_last=False)
    else:
        loader = DataLoader(ds, batch_size=batch, shuffle=False, num_workers=0, drop_last=False)
    return loader


# ---------------------------------------------------------------------------
# Training
# ---------------------------------------------------------------------------
def focal_bce(logit, target, gamma=2.0, label_smooth=0.05):
    """Binary focal loss with label smoothing on soft targets."""
    target = target * (1 - label_smooth) + 0.5 * label_smooth
    p = torch.sigmoid(logit)
    ce = F.binary_cross_entropy_with_logits(logit, target, reduction="none")
    pt = p * target + (1 - p) * (1 - target)
    return ((1 - pt) ** gamma * ce).mean()


def mixup_pair(x, y, alpha=0.2):
    if alpha <= 0:
        return x, y
    lam = np.random.beta(alpha, alpha)
    idx = torch.randperm(x.size(0), device=x.device)
    return lam * x + (1 - lam) * x[idx], lam * y + (1 - lam) * y[idx]


def train_one(model, loader, val_df, pep_idx, hla_idx, pep_emb, hla_emb, src2id,
              epochs=20, lr=3e-4, wd=1e-4, dann_lam=0.3, mixup_alpha=0.2,
              focal_gamma=2.0, label_smooth=0.05, log_prefix="seed0"):
    model = model.to(DEVICE)
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=wd)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=epochs)
    best_auc = 0.0
    best_state = None
    history = []
    for ep in range(epochs):
        model.train()
        losses = []
        for x, y, s in loader:
            x, y, s = x.to(DEVICE), y.float().to(DEVICE), s.long().to(DEVICE)
            x_m, y_m = mixup_pair(x, y, alpha=mixup_alpha)
            logit, dom = model(x_m, lam_dann=dann_lam)
            loss = focal_bce(logit, y_m, gamma=focal_gamma, label_smooth=label_smooth)
            if dom is not None:
                loss = loss + 0.5 * F.cross_entropy(dom, s)
            opt.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            opt.step()
            losses.append(loss.item())
        sched.step()
        # quick val AUROC (deterministic, no MC)
        if val_df is not None and len(val_df) >= 30:
            yv, sv = predict_mean(model, val_df, pep_idx, hla_idx, pep_emb, hla_emb, src2id, T=1)
            try:
                auc = roc_auc_score(val_df["label"].values, sv)
            except ValueError:
                auc = float("nan")
            history.append({"epoch": ep, "loss": float(np.mean(losses)), "val_auc": float(auc)})
            if auc > best_auc:
                best_auc = auc
                best_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}
            if (ep + 1) % 5 == 0 or ep == 0:
                print(f"  [{log_prefix}] ep{ep+1:>2d} loss={np.mean(losses):.4f} val_auc={auc:.3f}")
        else:
            history.append({"epoch": ep, "loss": float(np.mean(losses)), "val_auc": float("nan")})
    if best_state is not None:
        model.load_state_dict(best_state)
    return model, history


# ---------------------------------------------------------------------------
# Inference (MC Dropout) — vectorized tensor build
# ---------------------------------------------------------------------------
def _build_X(df, pep_idx, hla_idx, pep_emb, hla_emb):
    """Pre-build a single dense [N, 2H] tensor on CPU then push to GPU."""
    n = len(df)
    pi = np.fromiter((pep_idx[p] for p in df["peptide"].tolist()), dtype=np.int64, count=n)
    hi = np.fromiter((hla_idx[h] for h in df["HLA_norm"].tolist()), dtype=np.int64, count=n)
    pe = pep_emb[torch.from_numpy(pi)]
    he = hla_emb[torch.from_numpy(hi)]
    X = torch.cat([pe, he], dim=-1)
    return X


@torch.no_grad()
def predict_mean(model, df, pep_idx, hla_idx, pep_emb, hla_emb, src2id, T=30, batch=2048):
    X = _build_X(df, pep_idx, hla_idx, pep_emb, hla_emb).to(DEVICE)
    n = X.size(0)
    all_p = torch.zeros(T, n)
    model.eval()  # BN deterministic; F.dropout(training=True) inside model still fires
    for t in range(T):
        outs = []
        for i in range(0, n, batch):
            logit, _ = model(X[i:i+batch])
            outs.append(torch.sigmoid(logit).cpu())
        all_p[t] = torch.cat(outs)
    return df["label"].values, all_p.mean(0).numpy()


@torch.no_grad()
def predict_mc(model, df, pep_idx, hla_idx, pep_emb, hla_emb, src2id, T=30, batch=2048):
    X = _build_X(df, pep_idx, hla_idx, pep_emb, hla_emb).to(DEVICE)
    n = X.size(0)
    all_p = torch.zeros(T, n)
    model.eval()
    for t in range(T):
        outs = []
        for i in range(0, n, batch):
            logit, _ = model(X[i:i+batch])
            outs.append(torch.sigmoid(logit).cpu())
        all_p[t] = torch.cat(outs)
    return df["label"].values, all_p.mean(0).numpy(), all_p.std(0).numpy()


# ---------------------------------------------------------------------------
# Eval helpers
# ---------------------------------------------------------------------------
def ece(y, p, n_bins=15):
    """Expected Calibration Error."""
    y, p = np.asarray(y), np.asarray(p)
    bins = np.linspace(0, 1, n_bins + 1)
    ece_v = 0.0
    n = len(y)
    for lo, hi in zip(bins[:-1], bins[1:]):
        mask = (p >= lo) & (p < hi if hi < 1 else p <= hi)
        if mask.sum() == 0:
            continue
        ece_v += (mask.sum() / n) * abs(p[mask].mean() - y[mask].mean())
    return float(ece_v)


def metrics_block(y, p):
    if len(set(y)) < 2:
        return {"AUROC": None, "AUPRC": None, "Brier": None, "ECE": None,
                "n": int(len(y)), "n_pos": int(np.sum(y))}
    return {
        "AUROC": float(roc_auc_score(y, p)),
        "AUPRC": float(average_precision_score(y, p)),
        "Brier": float(brier_score_loss(y, p)),
        "ECE": ece(y, p),
        "n": int(len(y)),
        "n_pos": int(np.sum(y)),
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bundle", type=str, default="bundle.tsv")
    ap.add_argument("--embeddings", type=str, default="embeddings.pt")
    ap.add_argument("--out_dir", type=str, default=".")
    ap.add_argument("--epochs", type=int, default=20)
    ap.add_argument("--n_seeds", type=int, default=5)
    ap.add_argument("--T", type=int, default=30, help="MC Dropout samples")
    ap.add_argument("--mixup", type=float, default=0.2)
    ap.add_argument("--no_dann", action="store_true")
    ap.add_argument("--no_mixup", action="store_true")
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    print(f"loading bundle + embeddings…")
    bundle = pd.read_csv(args.bundle, sep="\t")
    bundle["HLA_norm"] = bundle["HLA_norm"].fillna("")

    emb = torch.load(args.embeddings, map_location="cpu", weights_only=False)
    pep_keys = emb["pep_keys"]
    hla_keys = emb["hla_keys"]
    pep_emb = emb["pep_emb"]
    hla_emb = emb["hla_emb"]
    pep_idx = {k: i for i, k in enumerate(pep_keys)}
    hla_idx = {k: i for i, k in enumerate(hla_keys)}
    venus_windows = emb["venus_windows"]
    print(f"  peps: {len(pep_keys)}, hlas: {len(hla_keys)}, hidden: {pep_emb.shape[1]}")

    # Filter bundle to rows whose peptide AND hla embedding exist
    def has_emb(p, h):
        return p in pep_idx and h in hla_idx

    nonvenus = bundle[~bundle["split"].isin(["ext_venus_test", "ext_venus_valid"])].copy()
    nonvenus = nonvenus[[has_emb(p, h) for p, h in zip(nonvenus["peptide"], nonvenus["HLA_norm"])]].copy()
    print(f"  rows with valid embeddings: {len(nonvenus)} / orig nonvenus")
    print(f"  per-split coverage: {nonvenus['split'].value_counts().to_dict()}")

    train_df = nonvenus[nonvenus["split"] == "train"].copy().reset_index(drop=True)
    src2id = {s: i for i, s in enumerate(sorted(train_df["source"].unique()))}
    print(f"  train sources: {src2id}")

    # ----- Cross-source LOSO + In-domain 5-fold ------------------------------
    # The two evaluations the RF baseline uses
    print("\n=== Within-source 5-fold ===")
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    fold_metrics_bayes = []
    fold_metrics_ensemble = []
    in_domain_preds = {"y": [], "p_bayes": [], "p_ens": [], "p_std": []}
    for fold, (tr, te) in enumerate(skf.split(train_df, train_df["label"].values)):
        tr_df = train_df.iloc[tr].reset_index(drop=True)
        te_df = train_df.iloc[te].reset_index(drop=True)
        loader_tr = make_loader(tr_df, pep_idx, hla_idx, pep_emb, hla_emb, src2id, balanced=True)
        # Train ensemble of n_seeds; first model also serves as "Bayesian"
        ens_probs = []
        ens_stds = []
        for seed in range(args.n_seeds):
            torch.manual_seed(seed); np.random.seed(seed)
            model = BayesianMLP(in_dim=2 * HIDDEN, hid=256,
                                n_sources=len(src2id), p_drop=0.3,
                                dann=not args.no_dann)
            model, _ = train_one(model, loader_tr, te_df, pep_idx, hla_idx,
                                 pep_emb, hla_emb, src2id,
                                 epochs=args.epochs,
                                 mixup_alpha=0.0 if args.no_mixup else args.mixup,
                                 log_prefix=f"fold{fold}-s{seed}")
            y, p_mean, p_std = predict_mc(model, te_df, pep_idx, hla_idx, pep_emb, hla_emb, src2id, T=args.T)
            ens_probs.append(p_mean)
            ens_stds.append(p_std)
            if seed == 0:
                bayes_p = p_mean
                bayes_s = p_std
        # Deep ensemble = mean of all seeds' MC means
        ens_probs = np.stack(ens_probs)
        ens_stds = np.stack(ens_stds)
        ens_mean = ens_probs.mean(0)
        ens_std = np.sqrt(ens_probs.var(0) + ens_stds.mean(0) ** 2)
        m_b = metrics_block(y, bayes_p); m_b["fold"] = fold; m_b["model"] = "BayesianMLP_seed0"
        m_e = metrics_block(y, ens_mean); m_e["fold"] = fold; m_e["model"] = "DeepEnsemble"
        fold_metrics_bayes.append(m_b)
        fold_metrics_ensemble.append(m_e)
        print(f"  fold{fold}: bayes AUROC={m_b['AUROC']:.3f}  ens AUROC={m_e['AUROC']:.3f}")
        in_domain_preds["y"].extend(y.tolist())
        in_domain_preds["p_bayes"].extend(bayes_p.tolist())
        in_domain_preds["p_ens"].extend(ens_mean.tolist())
        in_domain_preds["p_std"].extend(ens_std.tolist())

    pd.DataFrame(in_domain_preds).to_csv(out / "predictions_in_domain.tsv", sep="\t", index=False)

    # ----- Cross-source LOSO --------------------------------------------------
    print("\n=== Cross-source LOSO ===")
    loso_metrics = []
    for src in sorted(train_df["source"].unique()):
        tr_df = train_df[train_df["source"] != src].reset_index(drop=True)
        te_df = train_df[train_df["source"] == src].reset_index(drop=True)
        if te_df["label"].nunique() < 2:
            continue
        loader_tr = make_loader(tr_df, pep_idx, hla_idx, pep_emb, hla_emb, src2id, balanced=True)
        ens_probs = []
        for seed in range(args.n_seeds):
            torch.manual_seed(seed); np.random.seed(seed)
            model = BayesianMLP(in_dim=2 * HIDDEN, hid=256,
                                n_sources=len(src2id), p_drop=0.3,
                                dann=not args.no_dann)
            model, _ = train_one(model, loader_tr, None, pep_idx, hla_idx,
                                 pep_emb, hla_emb, src2id,
                                 epochs=args.epochs,
                                 mixup_alpha=0.0 if args.no_mixup else args.mixup,
                                 log_prefix=f"loso{src[:8]}-s{seed}")
            y, p_mean, p_std = predict_mc(model, te_df, pep_idx, hla_idx, pep_emb, hla_emb, src2id, T=args.T)
            ens_probs.append(p_mean)
        ens_mean = np.stack(ens_probs).mean(0)
        m = metrics_block(y, ens_mean); m["held_out_source"] = src
        loso_metrics.append(m)
        print(f"  LOSO {src:>22s}: AUROC={m['AUROC']:.3f}  AUPRC={m['AUPRC']:.3f}")

    pd.DataFrame(loso_metrics).to_csv(out / "loso_results.tsv", sep="\t", index=False)

    # ----- Train final ensemble on full pool for external evals --------------
    print("\n=== Final ensemble on full train pool ===")
    full_loader = make_loader(train_df, pep_idx, hla_idx, pep_emb, hla_emb, src2id, balanced=True)
    final_models = []
    for seed in range(args.n_seeds):
        torch.manual_seed(seed); np.random.seed(seed)
        m = BayesianMLP(in_dim=2 * HIDDEN, hid=256, n_sources=len(src2id),
                        p_drop=0.3, dann=not args.no_dann)
        m, _ = train_one(m, full_loader, None, pep_idx, hla_idx, pep_emb, hla_emb, src2id,
                         epochs=args.epochs,
                         mixup_alpha=0.0 if args.no_mixup else args.mixup,
                         log_prefix=f"final-s{seed}")
        final_models.append(m)

    def ensemble_predict(df, T=args.T):
        p_means, p_stds = [], []
        for m in final_models:
            _, pm, ps = predict_mc(m, df, pep_idx, hla_idx, pep_emb, hla_emb, src2id, T=T)
            p_means.append(pm); p_stds.append(ps)
        P = np.stack(p_means)  # [n_seeds, N]
        S = np.stack(p_stds)
        ens_mean = P.mean(0)
        ens_std = np.sqrt(P.var(0) + S.mean(0) ** 2)
        return ens_mean, ens_std

    # ----- ITSNdb -------------------------------------------------------------
    print("\n=== ITSNdb external ===")
    itsndb_results = []
    itsn = nonvenus[nonvenus["split"].isin(["ext_itsndb_main", "ext_itsndb_val"])].copy()
    if len(itsn) > 0:
        # Reload full ITSNdb to get in_master flag for no_overlap
        full_bundle = pd.read_csv(args.bundle, sep="\t")
        if "in_master" not in full_bundle.columns:
            full_bundle["in_master"] = False
        meta = full_bundle.set_index(["peptide", "HLA_norm", "split"])["in_master"].to_dict() if "in_master" in full_bundle.columns else {}
        itsn["in_master"] = itsn.apply(lambda r: bool(meta.get((r["peptide"], r["HLA_norm"], r["split"]), False)), axis=1)
        em, es = ensemble_predict(itsn)
        itsn["pred_mean"] = em; itsn["pred_std"] = es
        for sub_name, mask in [
            ("ITSNdb_main", itsn["split"] == "ext_itsndb_main"),
            ("ITSNdb_Val",  itsn["split"] == "ext_itsndb_val"),
            ("ITSNdb_combined", np.ones(len(itsn), dtype=bool)),
            ("ITSNdb_no_overlap", ~itsn["in_master"].values),
        ]:
            sub = itsn[mask]
            if len(sub) >= 5 and sub["label"].nunique() == 2:
                m = metrics_block(sub["label"].values, sub["pred_mean"].values)
                m["subset"] = sub_name
                itsndb_results.append(m)
                print(f"  {sub_name:>20s}: n={m['n']} pos={m['n_pos']} AUROC={m['AUROC']:.3f} AUPRC={m['AUPRC']:.3f} ECE={m['ECE']:.3f}")
        pd.DataFrame(itsndb_results).to_csv(out / "itsndb_bayesian_results.tsv", sep="\t", index=False)
        itsn.to_csv(out / "predictions_itsndb.tsv", sep="\t", index=False)

    # ----- VenusVaccine -------------------------------------------------------
    # Single-batch implementation: build one big DF of (window × top-10 HLA)
    # rows for ALL proteins, score once per ensemble member, then aggregate
    # back to per-protein. This avoids per-protein DataLoader overhead.
    print("\n=== VenusVaccine TumorBinary ===")
    venus_rows = bundle[bundle["split"].isin(["ext_venus_test", "ext_venus_valid"])].copy()
    venus_results = []
    if len(venus_rows) > 0:
        top_hlas = train_df["HLA_norm"].value_counts().head(10).index.tolist()
        top_hlas = [h for h in top_hlas if h in hla_idx]
        print(f"  scoring against top-10 train HLAs: {top_hlas}")
        big_rows = []
        for _, r in venus_rows.iterrows():
            pid = r["protein_id"]
            kmers = [k for k in venus_windows.get(pid, []) if k in pep_idx]
            for k in kmers:
                for h in top_hlas:
                    big_rows.append((pid, int(r["label"]), r["split"], k, h))
        big_df = pd.DataFrame(big_rows, columns=["protein_id", "label", "split", "peptide", "HLA_norm"])
        big_df["source"] = big_df["split"]
        print(f"  total (window×HLA) rows: {len(big_df):,}")
        if len(big_df) > 0:
            t_v = time.time()
            em, es = ensemble_predict(big_df, T=args.T)
            big_df["pred_mean"] = em; big_df["pred_std"] = es
            print(f"  ensemble inference done in {time.time()-t_v:.1f}s")
            # Aggregate per protein
            venus_pred_rows = []
            for (pid, label, split), g in big_df.groupby(["protein_id", "label", "split"]):
                scores = g["pred_mean"].values
                stds = g["pred_std"].values
                n_unique_w = g["peptide"].nunique()
                k = max(1, int(0.01 * len(scores)))
                venus_pred_rows.append({
                    "protein_id": pid, "label": int(label), "split": split,
                    "n_windows": int(n_unique_w),
                    "max_score": float(scores.max()),
                    "mean_score": float(scores.mean()),
                    "top10_mean": float(np.sort(scores)[::-1][:10].mean()),
                    "top1pct_mean": float(np.sort(scores)[::-1][:k].mean()),
                    "max_std": float(stds.max()),
                    "mean_std": float(stds.mean()),
                })
        vdf_pred = pd.DataFrame(venus_pred_rows)
        vdf_pred.to_csv(out / "predictions_venus.tsv", sep="\t", index=False)
        for split_name in ["ext_venus_test", "ext_venus_valid", "both"]:
            mask = (vdf_pred["split"] == split_name) if split_name != "both" else np.ones(len(vdf_pred), dtype=bool)
            sub = vdf_pred[mask]
            if len(sub) < 10 or sub["label"].nunique() < 2:
                continue
            for agg in ["max_score", "mean_score", "top10_mean", "top1pct_mean"]:
                m = metrics_block(sub["label"].values, sub[agg].values)
                m["split"] = split_name; m["aggregator"] = agg
                venus_results.append(m)
                print(f"  {split_name:>18s} {agg:>12s}: n={m['n']} AUROC={m['AUROC']:.3f}")
        pd.DataFrame(venus_results).to_csv(out / "venus_bayesian_results.tsv", sep="\t", index=False)

    # ----- Per-allele LOSO ----------------------------------------------------
    print("\n=== Per-allele LOSO ===")
    per_allele = []
    counts = train_df.groupby("HLA_norm")["label"].agg(["count", "sum", "mean"]).reset_index()
    counts.columns = ["HLA_norm", "n", "n_pos", "pos_rate"]
    eligible = counts[(counts["n"] >= 50) & (counts["pos_rate"] > 0.05) & (counts["pos_rate"] < 0.95)]
    print(f"  eligible alleles: {len(eligible)}")
    for _, row in eligible.iterrows():
        a = row["HLA_norm"]
        tr_df = train_df[train_df["HLA_norm"] != a].reset_index(drop=True)
        te_df = train_df[train_df["HLA_norm"] == a].reset_index(drop=True)
        if te_df["label"].nunique() < 2 or len(te_df) < 30:
            continue
        loader_tr = make_loader(tr_df, pep_idx, hla_idx, pep_emb, hla_emb, src2id, balanced=True)
        # Single seed (per-allele LOSO is k=11 fits → keep budget tight)
        torch.manual_seed(0); np.random.seed(0)
        model = BayesianMLP(in_dim=2 * HIDDEN, hid=256,
                            n_sources=len(src2id), p_drop=0.3,
                            dann=not args.no_dann)
        model, _ = train_one(model, loader_tr, None, pep_idx, hla_idx,
                             pep_emb, hla_emb, src2id,
                             epochs=args.epochs,
                             mixup_alpha=0.0 if args.no_mixup else args.mixup,
                             log_prefix=f"a{a[:8]}")
        y, p_mean, p_std = predict_mc(model, te_df, pep_idx, hla_idx, pep_emb, hla_emb, src2id, T=args.T)
        m = metrics_block(y, p_mean)
        m["allele"] = a
        m["n_test"] = len(te_df)
        m["n_pos_test"] = int(y.sum())
        per_allele.append(m)
        print(f"  {a:>18s} n={len(te_df)} pos={int(y.sum())} AUROC={m['AUROC']:.3f}")
        # Checkpoint after every allele
        pd.DataFrame(per_allele).to_csv(out / "per_allele_loso_bayesian.tsv", sep="\t", index=False)
    pd.DataFrame(per_allele).to_csv(out / "per_allele_loso_bayesian.tsv", sep="\t", index=False)

    # ----- Summary table ------------------------------------------------------
    summary_rows = []
    for r in fold_metrics_bayes:
        summary_rows.append({
            "eval": "in_domain_5fold", "model": "BayesianMLP_seed0",
            "fold": r["fold"], "AUROC": r["AUROC"], "AUPRC": r["AUPRC"], "ECE": r["ECE"], "n": r["n"]
        })
    for r in fold_metrics_ensemble:
        summary_rows.append({
            "eval": "in_domain_5fold", "model": "DeepEnsemble",
            "fold": r["fold"], "AUROC": r["AUROC"], "AUPRC": r["AUPRC"], "ECE": r["ECE"], "n": r["n"]
        })
    for r in loso_metrics:
        summary_rows.append({
            "eval": "cross_source_loso", "model": "DeepEnsemble",
            "held_out": r["held_out_source"], "AUROC": r["AUROC"], "AUPRC": r["AUPRC"], "ECE": r["ECE"], "n": r["n"]
        })
    for r in itsndb_results:
        summary_rows.append({
            "eval": "ITSNdb", "model": "DeepEnsemble",
            "subset": r["subset"], "AUROC": r["AUROC"], "AUPRC": r["AUPRC"], "ECE": r["ECE"],
            "n": r["n"], "n_pos": r["n_pos"],
        })
    for r in venus_results:
        summary_rows.append({
            "eval": "VenusVaccine", "model": "DeepEnsemble",
            "split": r["split"], "aggregator": r["aggregator"],
            "AUROC": r["AUROC"], "AUPRC": r["AUPRC"], "ECE": r["ECE"], "n": r["n"]
        })
    pd.DataFrame(summary_rows).to_csv(out / "results_summary.tsv", sep="\t", index=False)
    print(f"\nsaved → {out / 'results_summary.tsv'}  rows={len(summary_rows)}")


if __name__ == "__main__":
    main()
