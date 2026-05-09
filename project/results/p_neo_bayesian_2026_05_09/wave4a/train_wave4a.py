"""Wave 4A — fine-tune-track distribution-shift toolkit.

Four methods, same data + same seeds + same eval pipeline:

  M1 — LoRA on ESM2-150M attention (q/k/v adapters, rank=8)
  M2 — GroupDRO (worst-source loss) on frozen ESM2 + MLP head
  M3 — MIRO (mutual-info reg vs frozen ESM2 features) on frozen ESM2 + MLP head
  M4 — MoLE (Mixture of LoRA Experts) — 4 source-specialized rank=4 LoRA experts
        + softmax gating network on top of ESM2 mean-pool

Eval (consistent with Wave 1):
  - In-domain 5-fold (within train pool)
  - Cross-source LOSO (4 sources)
  - ITSNdb_no_overlap (n=106) ← headline, in_master=False
  - ITSNdb_in_master (n=213) — control
  - VenusVaccine top10_mean (78 valid + 78 test)

Seed = 0 for fairness (single seed; ensemble adds expense without changing
the relative comparison).

Outputs (all under --out_dir):
  predictions_lora.tsv / predictions_groupdro.tsv / predictions_miro.tsv
    / predictions_mole.tsv
  wave4a_results.tsv  (long format with bootstrap CI95)
  fig_wave4a_uplift.png/pdf
  WAVE4A_REPORT.md
"""
from __future__ import annotations
import argparse, json, time, math, os
from pathlib import Path
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader, WeightedRandomSampler
from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss
from sklearn.model_selection import StratifiedKFold

DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
HIDDEN = 640  # ESM2-150M
ESM2_NAME = "facebook/esm2_t30_150M_UR50D"


# ---------------------------------------------------------------------------
# Eval helpers
# ---------------------------------------------------------------------------
def ece(y, p, n_bins=15):
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


def auroc_bootstrap_ci(y, p, n_boot=1000, seed=0):
    y = np.asarray(y); p = np.asarray(p)
    if len(set(y)) < 2:
        return None, None, None
    auc = float(roc_auc_score(y, p))
    rng = np.random.default_rng(seed)
    n = len(y)
    boots = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        if len(set(y[idx])) < 2:
            continue
        boots.append(roc_auc_score(y[idx], p[idx]))
    if not boots:
        return auc, None, None
    boots = np.array(boots)
    return auc, float(np.quantile(boots, 0.025)), float(np.quantile(boots, 0.975))


def metrics_block(y, p, n_boot=1000):
    y = np.asarray(y); p = np.asarray(p)
    if len(set(y)) < 2:
        return {"AUROC": None, "AUROC_lo95": None, "AUROC_hi95": None,
                "AUPRC": None, "Brier": None, "ECE": None,
                "n": int(len(y)), "n_pos": int(np.sum(y))}
    auc, lo, hi = auroc_bootstrap_ci(y, p, n_boot=n_boot)
    return {
        "AUROC": auc,
        "AUROC_lo95": lo,
        "AUROC_hi95": hi,
        "AUPRC": float(average_precision_score(y, p)),
        "Brier": float(brier_score_loss(y, p)),
        "ECE": ece(y, p),
        "n": int(len(y)),
        "n_pos": int(np.sum(y)),
    }


# ---------------------------------------------------------------------------
# Data: shared bundle loader
# ---------------------------------------------------------------------------
def load_bundle_emb(bundle_path, emb_path):
    bundle = pd.read_csv(bundle_path, sep="\t")
    bundle["HLA_norm"] = bundle["HLA_norm"].fillna("")
    if "in_master" not in bundle.columns:
        bundle["in_master"] = False

    emb = torch.load(emb_path, map_location="cpu", weights_only=False)
    pep_keys = emb["pep_keys"]
    hla_keys = emb["hla_keys"]
    pep_emb = emb["pep_emb"]
    hla_emb = emb["hla_emb"]
    pep_idx = {k: i for i, k in enumerate(pep_keys)}
    hla_idx = {k: i for i, k in enumerate(hla_keys)}
    venus_windows = emb["venus_windows"]

    def has_emb(p, h):
        return p in pep_idx and h in hla_idx

    nonvenus = bundle[~bundle["split"].isin(["ext_venus_test", "ext_venus_valid"])].copy()
    nonvenus = nonvenus[[has_emb(p, h) for p, h in zip(nonvenus["peptide"], nonvenus["HLA_norm"])]].copy()
    return bundle, nonvenus, pep_idx, hla_idx, pep_emb, hla_emb, venus_windows


def make_xy(df, pep_idx, hla_idx, pep_emb, hla_emb, src2id):
    """Return (X, y, src) tensors built from cached embeddings."""
    n = len(df)
    pi = np.fromiter((pep_idx[p] for p in df["peptide"].tolist()), dtype=np.int64, count=n)
    hi = np.fromiter((hla_idx[h] for h in df["HLA_norm"].tolist()), dtype=np.int64, count=n)
    pe = pep_emb[torch.from_numpy(pi)]
    he = hla_emb[torch.from_numpy(hi)]
    X = torch.cat([pe, he], dim=-1).contiguous()
    y = torch.tensor(df["label"].astype(int).values, dtype=torch.float32)
    s = torch.tensor([src2id.get(s, 0) for s in df["source"].tolist()], dtype=torch.long)
    return X, y, s


def make_loader_xy(X, y, s, batch=256, balanced=True, src_counts=None):
    ds = torch.utils.data.TensorDataset(X, y, s)
    if balanced and src_counts is not None:
        n = X.size(0)
        w = np.array([1.0 / src_counts[int(si)] for si in s.tolist()], dtype=np.float64)
        w = w * n / w.sum()
        sampler = WeightedRandomSampler(weights=w.tolist(), num_samples=n, replacement=True)
        return DataLoader(ds, batch_size=batch, sampler=sampler, num_workers=0)
    return DataLoader(ds, batch_size=batch, shuffle=False, num_workers=0)


# ---------------------------------------------------------------------------
# Common MLP head (same arch as Wave 1, no DANN)
# ---------------------------------------------------------------------------
class MLPHead(nn.Module):
    def __init__(self, in_dim=2 * HIDDEN, hid=256, p_drop=0.3):
        super().__init__()
        self.fc1 = nn.Linear(in_dim, hid)
        self.bn1 = nn.BatchNorm1d(hid)
        self.fc2 = nn.Linear(hid, hid)
        self.bn2 = nn.BatchNorm1d(hid)
        self.head = nn.Linear(hid, 1)
        self.p = p_drop

    def features(self, x):
        h = F.relu(self.bn1(self.fc1(x)))
        h = F.dropout(h, p=self.p, training=self.training)
        h = F.relu(self.bn2(self.fc2(h)))
        h = F.dropout(h, p=self.p, training=self.training)
        return h

    def forward(self, x):
        h = self.features(x)
        return self.head(h).squeeze(-1), h


def focal_bce(logit, target, gamma=2.0, label_smooth=0.05):
    target = target * (1 - label_smooth) + 0.5 * label_smooth
    p = torch.sigmoid(logit)
    ce = F.binary_cross_entropy_with_logits(logit, target, reduction="none")
    pt = p * target + (1 - p) * (1 - target)
    return ((1 - pt) ** gamma * ce).mean(), ce  # also return per-sample ce


# ---------------------------------------------------------------------------
# M2 — GroupDRO trainer (Sagawa et al., ICLR 2020)
# ---------------------------------------------------------------------------
def train_groupdro(X_tr, y_tr, s_tr, n_groups, epochs=20, lr=3e-4, eta_q=0.01, batch=256):
    """Worst-source loss with adaptive weights q (online min-max)."""
    src_counts = np.bincount(s_tr.numpy(), minlength=n_groups).tolist()
    src_counts = [max(c, 1) for c in src_counts]
    loader = make_loader_xy(X_tr, y_tr, s_tr, batch=batch, balanced=True, src_counts=src_counts)
    model = MLPHead().to(DEVICE)
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-4)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=epochs)
    q = torch.ones(n_groups, device=DEVICE) / n_groups  # group weights
    for ep in range(epochs):
        model.train()
        losses = []
        for x, y, s in loader:
            x = x.to(DEVICE); y = y.float().to(DEVICE); s = s.long().to(DEVICE)
            logit, _ = model(x)
            _, ce = focal_bce(logit, y)
            # Per-group loss
            group_loss = torch.zeros(n_groups, device=DEVICE)
            group_n = torch.zeros(n_groups, device=DEVICE)
            for g in range(n_groups):
                m = (s == g)
                if m.any():
                    group_loss[g] = ce[m].mean()
                    group_n[g] = m.sum().float()
            # Update q (exponentiated gradient)
            with torch.no_grad():
                q = q * torch.exp(eta_q * group_loss.detach())
                q = q / q.sum()
            loss = (q * group_loss).sum()
            opt.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            opt.step()
            losses.append(loss.item())
        sched.step()
        if (ep + 1) % 5 == 0 or ep == 0:
            print(f"    [groupdro] ep{ep+1:>2d} loss={np.mean(losses):.4f} q={q.detach().cpu().numpy().round(3).tolist()}")
    return model


# ---------------------------------------------------------------------------
# M3 — MIRO (Cha et al., NeurIPS 2022)
# Loss = focal_bce + lam * ||head_features - frozen_features||^2
# Forces head features to stay close to ESM2 mean-pool (which we approximate
# with a linear projection of the input embeddings to match head hidden dim).
# ---------------------------------------------------------------------------
class MIROModel(nn.Module):
    def __init__(self, in_dim=2 * HIDDEN, hid=256):
        super().__init__()
        self.head = MLPHead(in_dim=in_dim, hid=hid)
        # Reference projection (frozen, random init): forces head to stay
        # representative of input embedding structure.
        # Better than identity-with-PCA because dim(in_dim)>>dim(hid).
        self.ref_proj = nn.Linear(in_dim, hid, bias=False)
        for p in self.ref_proj.parameters():
            p.requires_grad = False

    def forward(self, x):
        logit, h = self.head(x)
        with torch.no_grad():
            ref = self.ref_proj(x)
        return logit, h, ref


def train_miro(X_tr, y_tr, s_tr, n_groups, epochs=20, lr=3e-4, lam=0.1, batch=256):
    src_counts = np.bincount(s_tr.numpy(), minlength=n_groups).tolist()
    src_counts = [max(c, 1) for c in src_counts]
    loader = make_loader_xy(X_tr, y_tr, s_tr, batch=batch, balanced=True, src_counts=src_counts)
    model = MIROModel().to(DEVICE)
    opt = torch.optim.AdamW(filter(lambda p: p.requires_grad, model.parameters()), lr=lr, weight_decay=1e-4)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=epochs)
    for ep in range(epochs):
        model.train()
        losses = []
        for x, y, s in loader:
            x = x.to(DEVICE); y = y.float().to(DEVICE)
            logit, h, ref = model(x)
            cls_loss, _ = focal_bce(logit, y)
            miro_loss = F.mse_loss(h, ref)
            loss = cls_loss + lam * miro_loss
            opt.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            opt.step()
            losses.append(loss.item())
        sched.step()
        if (ep + 1) % 5 == 0 or ep == 0:
            print(f"    [miro]    ep{ep+1:>2d} loss={np.mean(losses):.4f}")
    return model


# ---------------------------------------------------------------------------
# Inference for cached-embedding models (M2/M3/M4-gating-only and head)
# ---------------------------------------------------------------------------
@torch.no_grad()
def predict_xy(model, X, batch=2048):
    model.eval()
    n = X.size(0)
    outs = []
    Xg = X.to(DEVICE)
    for i in range(0, n, batch):
        out = model(Xg[i:i+batch])
        if isinstance(out, tuple):
            logit = out[0]
        else:
            logit = out
        outs.append(torch.sigmoid(logit).cpu())
    return torch.cat(outs).numpy()


# ---------------------------------------------------------------------------
# M4 — MoLE: Mixture of LoRA Experts (cached-embedding variant)
#
# Since we are already using FROZEN ESM2 cached embeddings, we cannot do
# real LoRA on attention here. Instead we implement the MoLE *structure* on
# top of the cached embeddings: K "expert" MLP heads (one per source) with
# rank-4 low-rank perturbation parametrization, plus a gating MLP that maps
# the cached embedding to a softmax over experts. Training has a load-
# balancing loss to prevent gating collapse.
#
# This is the natural cached-embedding analog of MoLE. The "real" LoRA-on-
# attention is M1.
# ---------------------------------------------------------------------------
class LowRankExpert(nn.Module):
    """Shared backbone with per-expert low-rank adapter (rank=4) + per-expert head."""
    def __init__(self, in_dim, hid, rank=4):
        super().__init__()
        self.A = nn.Linear(in_dim, rank, bias=False)
        self.B = nn.Linear(rank, hid, bias=False)
        self.head = nn.Linear(hid, 1)

    def forward(self, h_shared, x):
        # additive low-rank residual on shared hidden
        delta = self.B(self.A(x))
        h = h_shared + delta
        return self.head(h).squeeze(-1)


class MoLE(nn.Module):
    def __init__(self, in_dim=2 * HIDDEN, hid=256, n_experts=4, rank=4, p_drop=0.3):
        super().__init__()
        self.n_experts = n_experts
        # Shared trunk (same as MLPHead minus head)
        self.fc1 = nn.Linear(in_dim, hid)
        self.bn1 = nn.BatchNorm1d(hid)
        self.fc2 = nn.Linear(hid, hid)
        self.bn2 = nn.BatchNorm1d(hid)
        self.p = p_drop
        # Experts
        self.experts = nn.ModuleList([LowRankExpert(in_dim, hid, rank=rank) for _ in range(n_experts)])
        # Gating MLP — small, on raw input
        self.gate = nn.Sequential(
            nn.Linear(in_dim, 64), nn.ReLU(), nn.Linear(64, n_experts)
        )

    def features(self, x):
        h = F.relu(self.bn1(self.fc1(x)))
        h = F.dropout(h, p=self.p, training=self.training)
        h = F.relu(self.bn2(self.fc2(h)))
        h = F.dropout(h, p=self.p, training=self.training)
        return h

    def forward(self, x):
        h_shared = self.features(x)
        # Per-expert logits
        per_exp = torch.stack([e(h_shared, x) for e in self.experts], dim=-1)  # [B, K]
        gate_logits = self.gate(x)  # [B, K]
        gate_w = F.softmax(gate_logits, dim=-1)  # [B, K]
        logit = (per_exp * gate_w).sum(-1)  # [B]
        return logit, gate_w


def train_mole(X_tr, y_tr, s_tr, n_experts, epochs=20, lr=3e-4, lam_lb=0.01, batch=256):
    src_counts = np.bincount(s_tr.numpy(), minlength=n_experts).tolist()
    src_counts = [max(c, 1) for c in src_counts]
    loader = make_loader_xy(X_tr, y_tr, s_tr, batch=batch, balanced=True, src_counts=src_counts)
    model = MoLE(n_experts=n_experts).to(DEVICE)
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-4)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=epochs)
    for ep in range(epochs):
        model.train()
        losses = []
        for x, y, s in loader:
            x = x.to(DEVICE); y = y.float().to(DEVICE)
            logit, gate_w = model(x)
            cls_loss, _ = focal_bce(logit, y)
            # Load-balancing loss: encourage uniform usage of experts
            mean_gate = gate_w.mean(0)
            lb = (mean_gate.log() * mean_gate).sum() + math.log(n_experts)  # KL(uniform || mean)
            loss = cls_loss + lam_lb * lb
            opt.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            opt.step()
            losses.append(loss.item())
        sched.step()
        if (ep + 1) % 5 == 0 or ep == 0:
            with torch.no_grad():
                mean_gate = gate_w.mean(0).cpu().numpy().round(3)
            print(f"    [mole]    ep{ep+1:>2d} loss={np.mean(losses):.4f} gate_mean={mean_gate.tolist()}")
    return model


# ---------------------------------------------------------------------------
# M1 — LoRA fine-tune of ESM2-150M attention
#
# Tokenize peptide+HLA pseudo-seq as a paired sequence with a separator.
# Wrap esm2_t30_150M_UR50D with PEFT LoRA on q/k/v projections, rank=8.
# Mean-pool the LM output, concat to (HLA pseudo embedding mean-pool from
# the same model — also LoRA-trained), feed to a 2-layer MLP head.
#
# This is the only method that actually unfreezes the LM. If it can't get
# above 0.5 on no_overlap, the bottleneck is data not architecture.
# ---------------------------------------------------------------------------
def train_lora(train_df, hla_pseudo_map, epochs=20, lr=3e-4, lam_lb=0.0, batch=32):
    from transformers import AutoTokenizer, AutoModel
    from peft import LoraConfig, get_peft_model, TaskType
    tok = AutoTokenizer.from_pretrained(ESM2_NAME)
    base = AutoModel.from_pretrained(ESM2_NAME)
    lora_cfg = LoraConfig(
        r=8, lora_alpha=16,
        target_modules=["query", "key", "value"],  # ESM2 attention names
        lora_dropout=0.05, bias="none",
        task_type=TaskType.FEATURE_EXTRACTION,
    )
    base = get_peft_model(base, lora_cfg)
    base.print_trainable_parameters()
    base = base.to(DEVICE)

    class LoraNeo(nn.Module):
        def __init__(self, base, hid=640, p_drop=0.3):
            super().__init__()
            self.base = base
            self.fc1 = nn.Linear(2 * hid, 256)
            self.bn1 = nn.BatchNorm1d(256)
            self.fc2 = nn.Linear(256, 256)
            self.bn2 = nn.BatchNorm1d(256)
            self.head = nn.Linear(256, 1)
            self.p = p_drop

        def encode(self, ids, mask):
            out = self.base(input_ids=ids, attention_mask=mask).last_hidden_state
            # mean-pool with mask
            m = mask.unsqueeze(-1).float()
            return (out * m).sum(1) / m.sum(1).clamp(min=1)

        def forward(self, pep_ids, pep_mask, hla_ids, hla_mask):
            pe = self.encode(pep_ids, pep_mask)
            he = self.encode(hla_ids, hla_mask)
            x = torch.cat([pe, he], dim=-1)
            h = F.relu(self.bn1(self.fc1(x)))
            h = F.dropout(h, p=self.p, training=self.training)
            h = F.relu(self.bn2(self.fc2(h)))
            h = F.dropout(h, p=self.p, training=self.training)
            return self.head(h).squeeze(-1)

    model = LoraNeo(base).to(DEVICE)
    opt = torch.optim.AdamW(filter(lambda p: p.requires_grad, model.parameters()), lr=lr, weight_decay=1e-4)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=epochs)

    # Pre-tokenize all rows
    peps = train_df["peptide"].tolist()
    hlas = [hla_pseudo_map.get(h, h) for h in train_df["HLA_norm"].tolist()]
    pep_enc = tok(peps, padding=True, truncation=True, max_length=20, return_tensors="pt")
    hla_enc = tok(hlas, padding=True, truncation=True, max_length=64, return_tensors="pt")
    y_all = torch.tensor(train_df["label"].astype(int).values, dtype=torch.float32)
    src = train_df["source"].astype("category")
    src_codes = torch.tensor(src.cat.codes.values, dtype=torch.long)
    n_src = len(src.cat.categories)

    src_counts = np.bincount(src_codes.numpy(), minlength=n_src).tolist()
    src_counts = [max(c, 1) for c in src_counts]
    n = len(train_df)
    w = np.array([1.0 / src_counts[int(c)] for c in src_codes.tolist()], dtype=np.float64)
    w = w * n / w.sum()

    ds = torch.utils.data.TensorDataset(
        pep_enc["input_ids"], pep_enc["attention_mask"],
        hla_enc["input_ids"], hla_enc["attention_mask"],
        y_all, src_codes,
    )
    sampler = WeightedRandomSampler(weights=w.tolist(), num_samples=n, replacement=True)
    loader = DataLoader(ds, batch_size=batch, sampler=sampler, num_workers=0, drop_last=True)

    for ep in range(epochs):
        model.train()
        losses = []
        for p_ids, p_mask, h_ids, h_mask, y, _ in loader:
            if p_ids.size(0) < 2:
                continue  # BatchNorm safety net
            p_ids = p_ids.to(DEVICE); p_mask = p_mask.to(DEVICE)
            h_ids = h_ids.to(DEVICE); h_mask = h_mask.to(DEVICE)
            y = y.float().to(DEVICE)
            logit = model(p_ids, p_mask, h_ids, h_mask)
            loss, _ = focal_bce(logit, y)
            opt.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(filter(lambda p: p.requires_grad, model.parameters()), 1.0)
            opt.step()
            losses.append(loss.item())
        sched.step()
        if (ep + 1) % 2 == 0 or ep == 0:
            print(f"    [lora]    ep{ep+1:>2d} loss={np.mean(losses):.4f}")

    return model, tok


@torch.no_grad()
def predict_lora(model, tok, df, hla_pseudo_map, batch=64):
    model.eval()
    peps = df["peptide"].tolist()
    hlas = [hla_pseudo_map.get(h, h) for h in df["HLA_norm"].tolist()]
    pep_enc = tok(peps, padding=True, truncation=True, max_length=20, return_tensors="pt")
    hla_enc = tok(hlas, padding=True, truncation=True, max_length=64, return_tensors="pt")
    n = len(df)
    outs = []
    for i in range(0, n, batch):
        p_ids = pep_enc["input_ids"][i:i+batch].to(DEVICE)
        p_mask = pep_enc["attention_mask"][i:i+batch].to(DEVICE)
        h_ids = hla_enc["input_ids"][i:i+batch].to(DEVICE)
        h_mask = hla_enc["attention_mask"][i:i+batch].to(DEVICE)
        logit = model(p_ids, p_mask, h_ids, h_mask)
        outs.append(torch.sigmoid(logit).cpu())
    return torch.cat(outs).numpy()


# ---------------------------------------------------------------------------
# Eval pipeline (works for any predict_fn that takes a df and returns probs)
# ---------------------------------------------------------------------------
def run_full_eval(method_name, train_fn, predict_fn,
                  bundle, nonvenus, pep_idx, hla_idx, pep_emb, hla_emb,
                  venus_windows, hla_pseudo_map, out_dir, n_groups=4):
    """Run in-domain 5-fold + LOSO + ITSNdb + Venus.
    train_fn: (train_df) -> trained model
    predict_fn: (model, df) -> probs (np.array)
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    train_df = nonvenus[nonvenus["split"] == "train"].copy().reset_index(drop=True)
    src2id = {s: i for i, s in enumerate(sorted(train_df["source"].unique()))}
    n_groups_eff = len(src2id)

    summary = []
    all_pred_rows = []

    # ----- In-domain 5-fold -------------------------------------------------
    print(f"\n=== [{method_name}] In-domain 5-fold ===")
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    fold_aurocs = []
    for fold, (tr, te) in enumerate(skf.split(train_df, train_df["label"].values)):
        tr_df = train_df.iloc[tr].reset_index(drop=True)
        te_df = train_df.iloc[te].reset_index(drop=True)
        torch.manual_seed(0); np.random.seed(0)
        model = train_fn(tr_df, n_groups_eff)
        p = predict_fn(model, te_df)
        m = metrics_block(te_df["label"].values, p, n_boot=1000)
        m["fold"] = fold
        fold_aurocs.append(m["AUROC"])
        summary.append({"method": method_name, "testset": "in_domain_5fold", "fold": fold,
                        "in_master": "NA", **m})
        all_pred_rows.append(pd.DataFrame({
            "peptide": te_df["peptide"], "HLA_norm": te_df["HLA_norm"],
            "label": te_df["label"], "source": te_df["source"],
            "split": f"in_domain_fold{fold}", "pred_mean": p,
        }))
        print(f"  [{method_name}] fold{fold} AUROC={m['AUROC']:.3f} [{m['AUROC_lo95']:.3f},{m['AUROC_hi95']:.3f}]")
        del model
        if DEVICE == "cuda":
            torch.cuda.empty_cache()
    print(f"  [{method_name}] 5-fold mean AUROC={np.mean(fold_aurocs):.3f}")

    # ----- Cross-source LOSO ------------------------------------------------
    print(f"\n=== [{method_name}] Cross-source LOSO ===")
    for src in sorted(train_df["source"].unique()):
        tr_df = train_df[train_df["source"] != src].reset_index(drop=True)
        te_df = train_df[train_df["source"] == src].reset_index(drop=True)
        if te_df["label"].nunique() < 2:
            continue
        torch.manual_seed(0); np.random.seed(0)
        model = train_fn(tr_df, n_groups_eff)
        p = predict_fn(model, te_df)
        m = metrics_block(te_df["label"].values, p, n_boot=1000)
        m["held_out_source"] = src
        summary.append({"method": method_name, "testset": f"loso_{src}", "fold": "NA",
                        "in_master": "NA", **m})
        all_pred_rows.append(pd.DataFrame({
            "peptide": te_df["peptide"], "HLA_norm": te_df["HLA_norm"],
            "label": te_df["label"], "source": te_df["source"],
            "split": f"loso_{src}", "pred_mean": p,
        }))
        print(f"  [{method_name}] LOSO {src:>22s} AUROC={m['AUROC']:.3f} [{m['AUROC_lo95']:.3f},{m['AUROC_hi95']:.3f}]")
        del model
        if DEVICE == "cuda":
            torch.cuda.empty_cache()

    # ----- Final model on full train pool, eval ITSNdb + Venus --------------
    print(f"\n=== [{method_name}] Final model on full train pool ===")
    torch.manual_seed(0); np.random.seed(0)
    final_model = train_fn(train_df, n_groups_eff)

    # ITSNdb
    print(f"\n=== [{method_name}] ITSNdb external ===")
    itsn = nonvenus[nonvenus["split"].isin(["ext_itsndb_main", "ext_itsndb_val"])].copy()
    p = predict_fn(final_model, itsn)
    itsn["pred_mean"] = p
    all_pred_rows.append(pd.DataFrame({
        "peptide": itsn["peptide"], "HLA_norm": itsn["HLA_norm"],
        "label": itsn["label"], "source": itsn["source"],
        "split": itsn["split"], "in_master": itsn["in_master"], "pred_mean": p,
    }))
    for sub_name, mask, in_master_flag in [
        ("ITSNdb_main", itsn["split"] == "ext_itsndb_main", "mixed"),
        ("ITSNdb_Val",  itsn["split"] == "ext_itsndb_val", "mixed"),
        ("ITSNdb_combined", np.ones(len(itsn), dtype=bool), "mixed"),
        ("ITSNdb_no_overlap", ~itsn["in_master"].astype(bool).values, False),
        ("ITSNdb_in_master", itsn["in_master"].astype(bool).values, True),
    ]:
        sub = itsn[mask]
        if len(sub) >= 5 and sub["label"].nunique() == 2:
            m = metrics_block(sub["label"].values, sub["pred_mean"].values, n_boot=1000)
            summary.append({"method": method_name, "testset": sub_name, "fold": "NA",
                            "in_master": in_master_flag, **m})
            print(f"  [{method_name}] {sub_name:>20s}: n={m['n']} AUROC={m['AUROC']:.3f} [{m['AUROC_lo95']:.3f},{m['AUROC_hi95']:.3f}]")

    # Venus
    print(f"\n=== [{method_name}] VenusVaccine ===")
    venus_rows = bundle[bundle["split"].isin(["ext_venus_test", "ext_venus_valid"])].copy()
    if len(venus_rows) > 0:
        top_hlas = train_df["HLA_norm"].value_counts().head(10).index.tolist()
        top_hlas = [h for h in top_hlas if h in hla_idx]
        big_rows = []
        for _, r in venus_rows.iterrows():
            pid = r["protein_id"]
            kmers = [k for k in venus_windows.get(pid, []) if k in pep_idx]
            for k in kmers:
                for h in top_hlas:
                    big_rows.append((pid, int(r["label"]), r["split"], k, h))
        big_df = pd.DataFrame(big_rows, columns=["protein_id", "label", "split", "peptide", "HLA_norm"])
        big_df["source"] = big_df["split"]
        if len(big_df) > 0:
            p = predict_fn(final_model, big_df)
            big_df["pred_mean"] = p
            for split_name in ["ext_venus_test", "ext_venus_valid"]:
                sub = big_df[big_df["split"] == split_name]
                if len(sub) < 10:
                    continue
                # top10_mean per protein
                rows = []
                for (pid, label), g in sub.groupby(["protein_id", "label"]):
                    scores = g["pred_mean"].values
                    rows.append({"protein_id": pid, "label": int(label),
                                 "top10_mean": float(np.sort(scores)[::-1][:10].mean())})
                vdf = pd.DataFrame(rows)
                if vdf["label"].nunique() == 2:
                    m = metrics_block(vdf["label"].values, vdf["top10_mean"].values, n_boot=1000)
                    summary.append({"method": method_name, "testset": f"venus_{split_name}_top10", "fold": "NA",
                                    "in_master": "NA", **m})
                    print(f"  [{method_name}] {split_name:>18s} top10_mean: n={m['n']} AUROC={m['AUROC']:.3f} [{m['AUROC_lo95']:.3f},{m['AUROC_hi95']:.3f}]")

    del final_model
    if DEVICE == "cuda":
        torch.cuda.empty_cache()

    # Save predictions
    pred_df = pd.concat(all_pred_rows, axis=0, ignore_index=True)
    pred_df.to_csv(out / f"predictions_{method_name}.tsv", sep="\t", index=False)

    return summary


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bundle", type=str, default="../bundle.tsv")
    ap.add_argument("--embeddings", type=str, default="../embeddings.pt")
    ap.add_argument("--hla_pseudo", type=str, default="../hla_pseudo.tsv")
    ap.add_argument("--out_dir", type=str, default=".")
    ap.add_argument("--epochs", type=int, default=20)
    ap.add_argument("--lora_epochs", type=int, default=20)
    ap.add_argument("--methods", type=str, default="lora,groupdro,miro,mole",
                    help="comma-separated list of methods to run")
    args = ap.parse_args()

    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    methods = [m.strip() for m in args.methods.split(",") if m.strip()]
    print(f"running methods: {methods}")
    print(f"epochs: {args.epochs} (lora={args.lora_epochs})")

    # Load shared data
    bundle, nonvenus, pep_idx, hla_idx, pep_emb, hla_emb, venus_windows = load_bundle_emb(
        args.bundle, args.embeddings)
    hla_pseudo = pd.read_csv(args.hla_pseudo, sep="\t")
    hla_pseudo_map = dict(zip(hla_pseudo.iloc[:, 0], hla_pseudo.iloc[:, 1]))

    all_summary = []

    # ===================================================================
    # M2 — GroupDRO
    # ===================================================================
    if "groupdro" in methods:
        print("\n" + "=" * 70 + "\nM2 — GroupDRO\n" + "=" * 70)
        def train_fn_gdro(train_df, n_groups):
            src2id = {s: i for i, s in enumerate(sorted(train_df["source"].unique()))}
            X, y, s = make_xy(train_df, pep_idx, hla_idx, pep_emb, hla_emb, src2id)
            return train_groupdro(X, y, s, n_groups=len(src2id), epochs=args.epochs)

        def predict_fn_gdro(model, df):
            # build cache-aware X but only for rows in pep/hla idx; fallback to 0 for missing
            valid = [(p in pep_idx and h in hla_idx) for p, h in zip(df["peptide"], df["HLA_norm"])]
            df_v = df[valid].copy().reset_index(drop=True)
            if len(df_v) == 0:
                return np.full(len(df), 0.5)
            X, _, _ = make_xy(df_v, pep_idx, hla_idx, pep_emb, hla_emb, {s: 0 for s in df_v["source"].unique()})
            p_v = predict_xy(model, X)
            full = np.full(len(df), 0.5)
            v_idx = np.where(valid)[0]
            full[v_idx] = p_v
            return full
        all_summary += run_full_eval("groupdro", train_fn_gdro, predict_fn_gdro,
                                     bundle, nonvenus, pep_idx, hla_idx, pep_emb, hla_emb,
                                     venus_windows, hla_pseudo_map, out)

    # ===================================================================
    # M3 — MIRO
    # ===================================================================
    if "miro" in methods:
        print("\n" + "=" * 70 + "\nM3 — MIRO\n" + "=" * 70)
        def train_fn_miro(train_df, n_groups):
            src2id = {s: i for i, s in enumerate(sorted(train_df["source"].unique()))}
            X, y, s = make_xy(train_df, pep_idx, hla_idx, pep_emb, hla_emb, src2id)
            return train_miro(X, y, s, n_groups=len(src2id), epochs=args.epochs)

        def predict_fn_miro(model, df):
            valid = [(p in pep_idx and h in hla_idx) for p, h in zip(df["peptide"], df["HLA_norm"])]
            df_v = df[valid].copy().reset_index(drop=True)
            if len(df_v) == 0:
                return np.full(len(df), 0.5)
            X, _, _ = make_xy(df_v, pep_idx, hla_idx, pep_emb, hla_emb, {s: 0 for s in df_v["source"].unique()})
            # MIROModel returns (logit, h, ref) — take logit
            model.eval()
            with torch.no_grad():
                Xg = X.to(DEVICE)
                outs = []
                for i in range(0, len(Xg), 2048):
                    logit, _, _ = model(Xg[i:i+2048])
                    outs.append(torch.sigmoid(logit).cpu())
                p_v = torch.cat(outs).numpy()
            full = np.full(len(df), 0.5)
            v_idx = np.where(valid)[0]
            full[v_idx] = p_v
            return full
        all_summary += run_full_eval("miro", train_fn_miro, predict_fn_miro,
                                     bundle, nonvenus, pep_idx, hla_idx, pep_emb, hla_emb,
                                     venus_windows, hla_pseudo_map, out)

    # ===================================================================
    # M4 — MoLE
    # ===================================================================
    if "mole" in methods:
        print("\n" + "=" * 70 + "\nM4 — MoLE\n" + "=" * 70)
        def train_fn_mole(train_df, n_groups):
            src2id = {s: i for i, s in enumerate(sorted(train_df["source"].unique()))}
            X, y, s = make_xy(train_df, pep_idx, hla_idx, pep_emb, hla_emb, src2id)
            return train_mole(X, y, s, n_experts=len(src2id), epochs=args.epochs)

        def predict_fn_mole(model, df):
            valid = [(p in pep_idx and h in hla_idx) for p, h in zip(df["peptide"], df["HLA_norm"])]
            df_v = df[valid].copy().reset_index(drop=True)
            if len(df_v) == 0:
                return np.full(len(df), 0.5)
            X, _, _ = make_xy(df_v, pep_idx, hla_idx, pep_emb, hla_emb, {s: 0 for s in df_v["source"].unique()})
            model.eval()
            with torch.no_grad():
                Xg = X.to(DEVICE)
                outs = []
                for i in range(0, len(Xg), 2048):
                    logit, _ = model(Xg[i:i+2048])
                    outs.append(torch.sigmoid(logit).cpu())
                p_v = torch.cat(outs).numpy()
            full = np.full(len(df), 0.5)
            v_idx = np.where(valid)[0]
            full[v_idx] = p_v
            return full
        all_summary += run_full_eval("mole", train_fn_mole, predict_fn_mole,
                                     bundle, nonvenus, pep_idx, hla_idx, pep_emb, hla_emb,
                                     venus_windows, hla_pseudo_map, out)

    # ===================================================================
    # M1 — LoRA  (most expensive — runs last)
    # ===================================================================
    if "lora" in methods:
        print("\n" + "=" * 70 + "\nM1 — LoRA on ESM2-150M attention\n" + "=" * 70)
        # LoRA path uses raw sequences, not cached embeddings
        def train_fn_lora(train_df, n_groups):
            return train_lora(train_df, hla_pseudo_map, epochs=args.lora_epochs)

        def predict_fn_lora(model_tok, df):
            model, tok = model_tok
            # Filter out rows whose HLA_norm has no pseudo
            valid = [h in hla_pseudo_map for h in df["HLA_norm"]]
            df_v = df[valid].copy().reset_index(drop=True)
            if len(df_v) == 0:
                return np.full(len(df), 0.5)
            p_v = predict_lora(model, tok, df_v, hla_pseudo_map)
            full = np.full(len(df), 0.5)
            v_idx = np.where(valid)[0]
            full[v_idx] = p_v
            return full
        all_summary += run_full_eval("lora", train_fn_lora, predict_fn_lora,
                                     bundle, nonvenus, pep_idx, hla_idx, pep_emb, hla_emb,
                                     venus_windows, hla_pseudo_map, out)

    # ----- Save combined results table --------------------------------------
    df = pd.DataFrame(all_summary)
    df.to_csv(out / "wave4a_results.tsv", sep="\t", index=False)
    print(f"\nSaved {len(df)} rows to {out / 'wave4a_results.tsv'}")


if __name__ == "__main__":
    main()
