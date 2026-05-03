#!/usr/bin/env python3
"""Negative controls for RAW (non-residualized) target.

For each panel: compute within-sample z-mean (no residualization) and run identical LOSO Ridge
on the cached 224 ResNet50 embeddings. Compares against random + housekeeping.
"""
from __future__ import annotations
import logging, os, time
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("neg_raw")

ROOT = Path(os.environ.get("THCA_ROOT", "/home/seungho/personal/THCA_data_analysis"))
PROCESSED = ROOT/"project/data/processed/GSE250521"
META = ROOT/"project/results/03_pathology_poc/tile_metadata_resid.tsv.gz"
EMB = ROOT/"project/results/03_pathology_poc/embeddings_resnet50_224.npz"
OUT = ROOT/"project/results/03_pathology_poc/negative_controls_raw_summary.tsv"

HK = ["ACTB","GAPDH","B2M","HPRT1","PPIA","RPL13A","RPLP0","TBP"]
RAI = ["TPO","DIO1","TSHR","PAX8","TG","FOXE1","NKX2-1","SLC5A5"]
N_RANDOM = 50
SEED = 42


def load_samples():
    out = {}
    for h in sorted(PROCESSED.rglob("*.scored.h5ad")):
        a = ad.read_h5ad(h)
        sid = a.obs["sample_id"].iloc[0]
        X = a.X
        if hasattr(X, "toarray"): X = X.toarray()
        Xl = np.log1p(X.astype(np.float32))
        qc = pd.DataFrame({
            "spot_id": a.obs.index.values,
            "sample_id": sid,
        }).reset_index(drop=True)
        out[sid] = (a.var.index.values, Xl, qc)
        log.info("loaded %s (%d×%d)", sid, *Xl.shape)
    return out


def panel_score_raw(samples_data, panel):
    rows = []
    for sid, (genes, Xl, qc) in samples_data.items():
        gset = set(panel)
        idx = [i for i, g in enumerate(genes) if g in gset]
        if len(idx) == 0: continue
        sub = Xl[:, idx]
        mu = sub.mean(0); sd = sub.std(0); sd[sd == 0] = 1
        z = (sub - mu) / sd
        score = z.mean(1)
        df = qc.copy()
        df["raw_score"] = score
        rows.append(df)
    if not rows: return pd.DataFrame()
    return pd.concat(rows, ignore_index=True)


def loso_pooled(meta, X, y_col):
    samples = sorted(meta.sample_id.unique())
    preds = np.full(len(meta), np.nan)
    for held in samples:
        tr = (meta.sample_id != held).values; te = ~tr
        y_tr = meta.loc[tr, y_col].values.astype(float)
        y_te = meta.loc[te, y_col].values.astype(float)
        keep_tr = ~np.isnan(y_tr); keep_te = ~np.isnan(y_te)
        if keep_tr.sum() < 50 or keep_te.sum() < 5: continue
        Xtr = X[tr][keep_tr]; Xte = X[te][keep_te]
        sc = StandardScaler().fit(Xtr)
        mdl = Ridge(alpha=1.0).fit(sc.transform(Xtr), y_tr[keep_tr])
        p = mdl.predict(sc.transform(Xte))
        te_idx = np.where(te)[0][keep_te]; preds[te_idx] = p
    y = meta[y_col].values.astype(float); keep = ~np.isnan(y) & ~np.isnan(preds)
    if keep.sum() < 50: return np.nan, np.nan
    sp = spearmanr(y[keep], preds[keep]).correlation
    thr = np.quantile(y[keep], 0.75)
    yhi = (y[keep] >= thr).astype(int)
    auc = roc_auc_score(yhi, preds[keep]) if 0 < yhi.sum() < len(yhi) else np.nan
    return sp, auc


def main():
    rng = np.random.default_rng(SEED)
    samples_data = load_samples()
    universe = None
    for sid, (genes, Xl, qc) in samples_data.items():
        present = set(np.array(genes)[(Xl > 0).mean(0) >= 0.30])
        universe = present if universe is None else (universe & present)
    universe = sorted(universe - set(RAI) - set(HK))
    log.info("universe: %d genes", len(universe))

    meta_full = pd.read_csv(META, sep="\t")
    meta = meta_full[meta_full.tile_size == 224].reset_index(drop=True)
    z = np.load(EMB); idx_map = {s:i for i,s in enumerate(z["spot_ids"])}
    keep = meta.spot_id.isin(idx_map).values
    meta = meta[keep].reset_index(drop=True)
    order = np.array([idx_map[s] for s in meta.spot_id.values])
    X = z["embeddings"][order]
    log.info("aligned %d tiles", len(meta))

    panels = [("housekeeping_raw", HK), ("RAI_8_real_raw", RAI)]
    for i in range(N_RANDOM):
        panels.append((f"random_raw_{i:04d}", list(rng.choice(universe, 8, replace=False))))

    rows = []
    for name, panel in panels:
        t0 = time.time()
        scored = panel_score_raw(samples_data, panel)
        m = meta.merge(scored[["spot_id","sample_id","raw_score"]], on=["spot_id","sample_id"], how="left")
        if m.raw_score.isna().sum() > 0.5*len(m): continue
        sp, auc = loso_pooled(m, X, "raw_score")
        rows.append({"panel": name, "n_genes": len(panel), "spearman_r": sp,
                     "auroc_top25": auc, "elapsed_s": time.time()-t0,
                     "panel_genes": ",".join(panel)})
        if name in ("housekeeping_raw", "RAI_8_real_raw") or name.endswith("0049"):
            log.info("[%s] sp=%.3f auc=%.3f (%.1fs)", name, sp, auc, time.time()-t0)
    df = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, sep="\t", index=False)
    log.info("wrote %s (%d panels)", OUT, len(df))
    rand = df[df.panel.str.startswith("random_raw_")]
    real = df[df.panel == "RAI_8_real_raw"]
    hk = df[df.panel == "housekeeping_raw"]
    print(f"\n=== RAW NEG CTRL SUMMARY ===")
    print(f"RAI_8 real:    sp={real.spearman_r.iloc[0]:.3f} auc={real.auroc_top25.iloc[0]:.3f}")
    print(f"housekeeping:  sp={hk.spearman_r.iloc[0]:.3f} auc={hk.auroc_top25.iloc[0]:.3f}")
    print(f"random (n=50): mean sp={rand.spearman_r.mean():.3f} max={rand.spearman_r.max():.3f}")
    print(f"               mean auc={rand.auroc_top25.mean():.3f} max={rand.auroc_top25.max():.3f}")


if __name__ == "__main__":
    main()
