#!/usr/bin/env python
"""Wave 5B — multi-task multi-head training.

Architecture
------------
[ESM2-150M frozen pep 640d ⊕ hla 640d → 1280d] → shared MLP 256 → 4 task heads:
  Head A  immunogenicity (binary, focal BCE), MC Dropout for Bayesian σ
  Head B  binding intensity = -log10(aff_nM)  (regression on z-scored target)
  Head C  presentation_score                  (regression on z-scored target)
  Head D  -presentation_percentile/100        (regression on z-scored target)

Loss = w_A · focal_BCE(A) + w_B · MSE(B) + w_C · MSE(C) + w_D · MSE(D)
Default weights: 1.0 / 0.3 / 0.3 / 0.3 (overridable for ablations).

Eval
----
- Head A AUROC + 1000-bootstrap 95% CI on:
    ITSNdb_combined / ITSNdb_no_overlap / ITSNdb_in_master
- Head B/C/D Pearson + Spearman vs ground-truth MHCflurry on
    ITSNdb_no_overlap (sanity that aux heads generalize).
- Per-allele LOSO via Head A (same protocol as wave1).

CLI:
  python train_multitask.py --tag full           (default 4 heads)
  python train_multitask.py --tag A_only --heads A
  python train_multitask.py --tag AB --heads A,B
  python train_multitask.py --tag AC --heads A,C
  python train_multitask.py --tag AD --heads A,D
"""
from __future__ import annotations
import argparse, json, time
from pathlib import Path
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader, TensorDataset, WeightedRandomSampler
from sklearn.metrics import roc_auc_score
from scipy.stats import pearsonr, spearmanr

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09")
W5B = ROOT / "wave5b"
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
SEED = 0


# ---------------------------------------------------------------------------
def set_seed(s):
    np.random.seed(s); torch.manual_seed(s)


def auroc_with_ci(y, s, n_boot=1000, seed=0):
    y = np.asarray(y, dtype=int); s = np.asarray(s, dtype=float)
    if len(np.unique(y)) < 2:
        return float("nan"), float("nan"), float("nan")
    auc = roc_auc_score(y, s)
    rng = np.random.default_rng(seed)
    n = len(y); aucs = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        yy = y[idx]; ss = s[idx]
        if len(np.unique(yy)) < 2:
            continue
        aucs.append(roc_auc_score(yy, ss))
    if not aucs:
        return float(auc), float("nan"), float("nan")
    return float(auc), float(np.percentile(aucs, 2.5)), float(np.percentile(aucs, 97.5))


# ---------------------------------------------------------------------------
class MultiTaskMLP(nn.Module):
    """[1280] → 256 (shared, BN+ReLU+Drop)×2 → {head_A logit, head_B/C/D scalar regression}.

    MC Dropout: F.dropout with training=True is always-on at eval time on Head A
    feature path (only used when sample_T>1).
    """
    def __init__(self, in_dim=1280, hid=256, p_drop=0.3, heads=("A", "B", "C", "D")):
        super().__init__()
        self.heads = tuple(heads)
        self.fc1 = nn.Linear(in_dim, hid)
        self.bn1 = nn.BatchNorm1d(hid)
        self.fc2 = nn.Linear(hid, hid)
        self.bn2 = nn.BatchNorm1d(hid)
        self.head_A = nn.Linear(hid, 1) if "A" in heads else None
        self.head_B = nn.Linear(hid, 1) if "B" in heads else None
        self.head_C = nn.Linear(hid, 1) if "C" in heads else None
        self.head_D = nn.Linear(hid, 1) if "D" in heads else None
        self.p = p_drop

    def features(self, x, mc=True):
        h = F.relu(self.bn1(self.fc1(x)))
        h = F.dropout(h, p=self.p, training=mc or self.training)
        h = F.relu(self.bn2(self.fc2(h)))
        h = F.dropout(h, p=self.p, training=mc or self.training)
        return h

    def forward(self, x, mc=True):
        h = self.features(x, mc=mc)
        out = {}
        if self.head_A is not None:
            out["A"] = self.head_A(h).squeeze(-1)
        if self.head_B is not None:
            out["B"] = self.head_B(h).squeeze(-1)
        if self.head_C is not None:
            out["C"] = self.head_C(h).squeeze(-1)
        if self.head_D is not None:
            out["D"] = self.head_D(h).squeeze(-1)
        return out


def focal_bce(logit, target, gamma=2.0, label_smooth=0.05):
    target = target * (1 - label_smooth) + 0.5 * label_smooth
    p = torch.sigmoid(logit)
    ce = F.binary_cross_entropy_with_logits(logit, target, reduction="none")
    pt = p * target + (1 - p) * (1 - target)
    return ((1 - pt) ** gamma * ce).mean()


# ---------------------------------------------------------------------------
def build_X(df, pep_idx, hla_idx, pep_emb, hla_emb):
    pi = np.fromiter((pep_idx[p] for p in df["peptide"].tolist()), dtype=np.int64, count=len(df))
    hi = np.fromiter((hla_idx[h] for h in df["hla"].tolist()), dtype=np.int64, count=len(df))
    pe = pep_emb[torch.from_numpy(pi)]
    he = hla_emb[torch.from_numpy(hi)]
    return torch.cat([pe, he], dim=-1).contiguous()


def make_loader(df, pep_idx, hla_idx, pep_emb, hla_emb, batch=128, balanced=True):
    X = build_X(df, pep_idx, hla_idx, pep_emb, hla_emb)
    y_A = torch.tensor(df["label_A"].astype(int).values, dtype=torch.float32)
    y_B = torch.tensor(df["label_B_z"].astype(float).values, dtype=torch.float32)
    y_C = torch.tensor(df["label_C_z"].astype(float).values, dtype=torch.float32)
    y_D = torch.tensor(df["label_D_z"].astype(float).values, dtype=torch.float32)
    ds = TensorDataset(X, y_A, y_B, y_C, y_D)
    if balanced:
        n_per_src = df["source"].value_counts().to_dict()
        w = df["source"].map(lambda s: 1.0 / n_per_src[s]).values
        w = w * len(df) / w.sum()
        sampler = WeightedRandomSampler(weights=w.tolist(), num_samples=len(df), replacement=True)
        return DataLoader(ds, batch_size=batch, sampler=sampler, num_workers=0)
    return DataLoader(ds, batch_size=batch, shuffle=False, num_workers=0)


# ---------------------------------------------------------------------------
@torch.no_grad()
def predict(model, df, pep_idx, hla_idx, pep_emb, hla_emb, T=30, batch=2048, mc=True):
    """T-sample MC Dropout prediction. Returns dict {head: [N, T] sigmoid for A, raw for B/C/D}."""
    X = build_X(df, pep_idx, hla_idx, pep_emb, hla_emb).to(DEVICE)
    n = X.size(0)
    out = {h: torch.zeros(T, n) for h in model.heads}
    model.eval()
    for t in range(T):
        for i in range(0, n, batch):
            o = model(X[i:i+batch], mc=mc)
            for h, v in o.items():
                if h == "A":
                    out[h][t, i:i+batch] = torch.sigmoid(v).cpu()
                else:
                    out[h][t, i:i+batch] = v.cpu()
    return {h: out[h].numpy() for h in out}


# ---------------------------------------------------------------------------
def train_loop(model, loader, val_df, ctx, weights, heads, epochs=25, lr=3e-4, wd=1e-4,
               focal_gamma=2.0, label_smooth=0.05, log_prefix=""):
    model = model.to(DEVICE)
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=wd)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=epochs)
    best_auc, best_state = 0.0, None
    history = []
    for ep in range(epochs):
        model.train()
        losses = []
        for x, yA, yB, yC, yD in loader:
            x = x.to(DEVICE); yA = yA.to(DEVICE)
            yB = yB.to(DEVICE); yC = yC.to(DEVICE); yD = yD.to(DEVICE)
            o = model(x)
            loss = 0.0
            if "A" in heads and "A" in o:
                loss = loss + weights["A"] * focal_bce(o["A"], yA, gamma=focal_gamma,
                                                      label_smooth=label_smooth)
            if "B" in heads and "B" in o:
                loss = loss + weights["B"] * F.mse_loss(o["B"], yB)
            if "C" in heads and "C" in o:
                loss = loss + weights["C"] * F.mse_loss(o["C"], yC)
            if "D" in heads and "D" in o:
                loss = loss + weights["D"] * F.mse_loss(o["D"], yD)
            opt.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            opt.step()
            losses.append(float(loss.detach()))
        sched.step()
        # quick val AUROC (Head A) on the held-out part of train (last 200 rows random)
        if val_df is not None and "A" in heads:
            preds = predict(model, val_df, *ctx, T=1, mc=False)
            try:
                auc = roc_auc_score(val_df["label_A"].astype(int).values, preds["A"][0])
            except ValueError:
                auc = float("nan")
            history.append(dict(epoch=ep, loss=float(np.mean(losses)), val_auc=float(auc)))
            if not np.isnan(auc) and auc > best_auc:
                best_auc = auc
                best_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}
            if (ep + 1) % 5 == 0 or ep == 0:
                print(f"  [{log_prefix}] ep{ep+1:>2d} loss={np.mean(losses):.4f} val_auc={auc:.3f}",
                      flush=True)
        else:
            history.append(dict(epoch=ep, loss=float(np.mean(losses)), val_auc=float("nan")))
    if best_state is not None:
        model.load_state_dict(best_state)
    return model, history


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="full")
    ap.add_argument("--heads", default="A,B,C,D")
    ap.add_argument("--w-A", type=float, default=1.0)
    ap.add_argument("--w-B", type=float, default=0.3)
    ap.add_argument("--w-C", type=float, default=0.3)
    ap.add_argument("--w-D", type=float, default=0.3)
    ap.add_argument("--epochs", type=int, default=25)
    ap.add_argument("--bs", type=int, default=128)
    ap.add_argument("--T", type=int, default=30)
    ap.add_argument("--seeds", default="0,1,2")
    args = ap.parse_args()

    heads = tuple([h.strip() for h in args.heads.split(",") if h.strip()])
    weights = dict(A=args.w_A, B=args.w_B, C=args.w_C, D=args.w_D)
    print(f"[cfg] tag={args.tag} heads={heads} weights={weights} epochs={args.epochs}", flush=True)

    set_seed(SEED)
    sup = pd.read_csv(W5B / "multitask_supervision.tsv", sep="\t",
                      dtype={"in_master": "boolean"})
    print(f"[load] supervision rows={len(sup)}", flush=True)
    print(sup["split"].value_counts().to_string(), flush=True)

    emb = torch.load(W5B / "embeddings.pt", weights_only=False)
    pep_keys, hla_keys = emb["pep_keys"], emb["hla_keys"]
    pep_emb, hla_emb = emb["pep_emb"], emb["hla_emb"]
    pep_idx = {k: i for i, k in enumerate(pep_keys)}
    hla_idx = {k: i for i, k in enumerate(hla_keys)}

    # filter to rows whose pep+hla both have embeddings (i.e. HLA pseudo present)
    keep = sup.apply(lambda r: r["peptide"] in pep_idx and r["hla"] in hla_idx, axis=1)
    sup = sup[keep].reset_index(drop=True)
    print(f"[filter] rows with pep+hla embed: {len(sup)}", flush=True)
    print(sup["split"].value_counts().to_string(), flush=True)

    ctx = (pep_idx, hla_idx, pep_emb, hla_emb)

    # train/val/test split
    tr = sup[sup["split"] == "train"].copy()
    its = sup[sup["split"].isin(["ext_itsndb_main", "ext_itsndb_val"])].copy()
    print(f"[split] train={len(tr)}  itsndb={len(its)}", flush=True)

    # held-out validation: stratified 10% of train
    rng = np.random.default_rng(SEED)
    vidx = rng.choice(len(tr), size=max(50, len(tr)//10), replace=False)
    vmask = np.zeros(len(tr), dtype=bool); vmask[vidx] = True
    val_df = tr[vmask].reset_index(drop=True)
    tr_df = tr[~vmask].reset_index(drop=True)
    print(f"[split] tr={len(tr_df)}  val={len(val_df)}", flush=True)

    # ensemble training (3 seeds default)
    seeds = [int(s) for s in args.seeds.split(",")]
    ensemble_pred_its = {h: [] for h in heads}
    ensemble_pred_its_std = {h: [] for h in heads}
    histories = []
    t0 = time.time()
    for sd in seeds:
        set_seed(sd)
        torch.manual_seed(sd)
        model = MultiTaskMLP(in_dim=pep_emb.shape[1] + hla_emb.shape[1], hid=256,
                             p_drop=0.3, heads=heads)
        loader = make_loader(tr_df, *ctx, batch=args.bs, balanced=True)
        model, hist = train_loop(model, loader, val_df, ctx, weights, heads,
                                 epochs=args.epochs, lr=3e-4,
                                 log_prefix=f"{args.tag}/s{sd}")
        histories.append(hist)
        # Predict on ITSNdb with MC Dropout (T samples)
        preds_its = predict(model, its, *ctx, T=args.T, mc=True)
        for h in heads:
            arr = preds_its[h]  # [T, N]
            ensemble_pred_its[h].append(arr.mean(0))
            ensemble_pred_its_std[h].append(arr.std(0))
    print(f"[train] {len(seeds)} seeds done in {time.time()-t0:.1f}s", flush=True)

    # Ensemble across seeds: mean of seed-means
    its_pred = {h: np.mean(np.stack(ensemble_pred_its[h]), axis=0) for h in heads}
    its_std = {h: np.mean(np.stack(ensemble_pred_its_std[h]), axis=0) for h in heads}

    # Per-row predictions
    pred_df = its.copy()
    for h in heads:
        pred_df[f"pred_{h}_mean"] = its_pred[h]
        pred_df[f"pred_{h}_std"] = its_std[h]
    pred_df.to_csv(W5B / f"predictions_{args.tag}.tsv", sep="\t", index=False)
    print(f"[write] predictions_{args.tag}.tsv", flush=True)

    # Save model from last seed (representative)
    if args.tag == "full":
        torch.save(model.state_dict(), W5B / "model_multitask.pt")
        # Save predictions_wave5b.tsv as canonical
        pred_df.to_csv(W5B / "predictions_wave5b.tsv", sep="\t", index=False)

    # ---- Evaluate Head A on ITSNdb stratified ----
    rows = []
    if "A" in heads:
        sA = its_pred["A"]
        # Combined
        y = pred_df["label_A"].astype(int).values
        auc, lo, hi = auroc_with_ci(y, sA)
        rows.append(dict(head="A", testset="ITSNdb_combined", n=len(y),
                         n_pos=int(y.sum()), AUROC=auc, CI_lo95=lo, CI_hi95=hi))
        # No-overlap
        m = (pred_df["in_master"] == False).values
        if m.sum() > 0:
            auc, lo, hi = auroc_with_ci(y[m], sA[m])
            rows.append(dict(head="A", testset="ITSNdb_no_overlap", n=int(m.sum()),
                             n_pos=int(y[m].sum()), AUROC=auc, CI_lo95=lo, CI_hi95=hi))
        # In-master
        m2 = (pred_df["in_master"] == True).values
        if m2.sum() > 0:
            auc, lo, hi = auroc_with_ci(y[m2], sA[m2])
            rows.append(dict(head="A", testset="ITSNdb_in_master", n=int(m2.sum()),
                             n_pos=int(y[m2].sum()), AUROC=auc, CI_lo95=lo, CI_hi95=hi))
    # Sanity: aux head correlation with their MHCflurry ground-truth on no-overlap
    m = (pred_df["in_master"] == False).values
    for h, gt_col in [("B", "label_B_z"), ("C", "label_C_z"), ("D", "label_D_z")]:
        if h in heads and m.sum() > 5:
            yhat = its_pred[h][m]
            ygt = pred_df.loc[m, gt_col].astype(float).values
            try:
                rp, _ = pearsonr(yhat, ygt)
                rs, _ = spearmanr(yhat, ygt)
            except Exception:
                rp = rs = float("nan")
            rows.append(dict(head=h, testset="ITSNdb_no_overlap_corr", n=int(m.sum()),
                             n_pos=-1, AUROC=float("nan"),
                             CI_lo95=float(rp), CI_hi95=float(rs)))
    res = pd.DataFrame(rows)
    res.to_csv(W5B / f"results_{args.tag}.tsv", sep="\t", index=False)
    print(res.to_string(index=False), flush=True)

    # Save train history
    (W5B / f"history_{args.tag}.json").write_text(json.dumps(histories, indent=2))


if __name__ == "__main__":
    main()
