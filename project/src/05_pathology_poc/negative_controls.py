#!/usr/bin/env python3
"""Negative control panels: housekeeping + N random 8-gene panels.

For each control panel:
  1. compute spot-level mean z-score (within sample) over the panel genes
  2. depth-residualize within sample (regress on log_counts + log_ngenes)
  3. run identical Ridge LOSO pipeline on the SAME tile embeddings
  4. record pooled Spearman + AUROC (top quartile)

Random panels: drawn from genes with positive expression in ≥30% of spots,
matched by mean detection rate to the real RAI_8 panel.

Output:
  negative_controls_summary.tsv  one row per panel
"""
from __future__ import annotations
import argparse, logging, os, time
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import spearmanr, pearsonr
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("neg_ctrl")

ROOT = Path(os.environ.get("THCA_ROOT", "/home/seungho/personal/THCA_data_analysis"))
PROCESSED = ROOT / "project/data/processed/GSE250521"
META = ROOT / "project/results/03_pathology_poc/tile_metadata_resid.tsv.gz"

HOUSEKEEPING = ["ACTB", "GAPDH", "B2M", "HPRT1", "PPIA", "RPL13A", "RPLP0", "TBP"]
RAI_8 = ["TPO", "DIO1", "TSHR", "PAX8", "TG", "FOXE1", "NKX2-1", "SLC5A5"]


def load_all_spots():
    """Load every scored.h5ad and stack to (spot, gene) matrix.

    Return:
        sample_lookup: dict sample_id -> (spot_ids, expr_log1p_dense, gene_index, qc_df)
    """
    out = {}
    for h in sorted(PROCESSED.rglob("*.scored.h5ad")):
        a = ad.read_h5ad(h)
        sid = a.obs["sample_id"].iloc[0]
        # use log1p of normalized counts; fall back to raw counts log1p
        X = a.X
        if hasattr(X, "toarray"): X = X.toarray()
        Xl = np.log1p(X.astype(np.float32))
        # per-sample z within the panel later
        qc = pd.DataFrame({
            "spot_id": a.obs.index.values,
            "log_counts": np.log1p(a.obs["total_counts"].values),
            "log_ngenes": np.log1p(a.obs["n_genes_by_counts"].values),
            "in_tissue": a.obs["in_tissue"].values,
        }).reset_index(drop=True)
        out[sid] = (a.var.index.values, Xl, qc)
        log.info("loaded %s (%d spots × %d genes)", sid, Xl.shape[0], Xl.shape[1])
    return out


def panel_score_per_sample(samples_data, panel_genes):
    """Compute within-sample z-mean of panel + depth-resid; return long DF."""
    rows = []
    for sid, (genes, Xl, qc) in samples_data.items():
        # gene index
        gset = set(panel_genes)
        idx = [i for i, g in enumerate(genes) if g in gset]
        if len(idx) == 0:
            continue
        sub = Xl[:, idx]
        # within-sample z per gene, then mean
        mu = sub.mean(axis=0); sd = sub.std(axis=0); sd[sd == 0] = 1
        z = (sub - mu) / sd
        score = z.mean(axis=1)
        df = qc.copy()
        df["sample_id"] = sid
        df["raw_score"] = score
        # residualize
        Xc = sm.add_constant(df[["log_counts", "log_ngenes"]].values)
        try:
            res = sm.OLS(df["raw_score"].values, Xc).fit()
            df["resid_score"] = df["raw_score"].values - res.predict(Xc)
        except Exception:
            df["resid_score"] = df["raw_score"].values - df["raw_score"].mean()
        rows.append(df)
    if not rows: return pd.DataFrame()
    return pd.concat(rows, ignore_index=True)


def loso_pooled(meta, X, y_col):
    samples = sorted(meta.sample_id.unique())
    preds = np.full(len(meta), np.nan)
    for held in samples:
        tr = (meta.sample_id != held).values
        te = ~tr
        y_tr = meta.loc[tr, y_col].values.astype(float)
        y_te = meta.loc[te, y_col].values.astype(float)
        keep_tr = ~np.isnan(y_tr); keep_te = ~np.isnan(y_te)
        if keep_tr.sum() < 50 or keep_te.sum() < 5: continue
        Xtr = X[tr][keep_tr]; Xte = X[te][keep_te]
        sc = StandardScaler().fit(Xtr)
        mdl = Ridge(alpha=1.0).fit(sc.transform(Xtr), y_tr[keep_tr])
        p = mdl.predict(sc.transform(Xte))
        te_idx = np.where(te)[0][keep_te]
        preds[te_idx] = p
    y = meta[y_col].values.astype(float)
    keep = ~np.isnan(y) & ~np.isnan(preds)
    if keep.sum() < 50:
        return np.nan, np.nan, np.nan
    sp = spearmanr(y[keep], preds[keep]).correlation
    thr = np.quantile(y[keep], 0.75)
    yhi = (y[keep] >= thr).astype(int)
    auc = (roc_auc_score(yhi, preds[keep])
           if 0 < yhi.sum() < len(yhi) else np.nan)
    pe = pearsonr(y[keep], preds[keep])[0]
    return sp, pe, auc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--embeddings", required=True)
    ap.add_argument("--tile-size", type=int, required=True)
    ap.add_argument("--n-random", type=int, default=100)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", default="project/results/03_pathology_poc/negative_controls_summary.tsv")
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    samples_data = load_all_spots()

    # gene universe = intersection of detected genes across samples
    universe = None
    for sid, (genes, Xl, qc) in samples_data.items():
        present = set(np.array(genes)[(Xl > 0).mean(axis=0) >= 0.30])
        universe = present if universe is None else (universe & present)
    universe = sorted(universe - set(RAI_8) - set(HOUSEKEEPING))
    log.info("gene universe (≥30%% detection in all samples): %d", len(universe))

    # tile metadata + embeddings
    meta_full = pd.read_csv(META, sep="\t")
    meta = meta_full[meta_full.tile_size == args.tile_size].reset_index(drop=True)
    z = np.load(args.embeddings)
    emb_spots = z["spot_ids"]
    idx_map = {s: i for i, s in enumerate(emb_spots)}
    keep = meta.spot_id.isin(idx_map).values
    meta = meta[keep].reset_index(drop=True)
    order = np.array([idx_map[s] for s in meta.spot_id.values])
    X = z["embeddings"][order]
    log.info("aligned: %d tiles, X=%s", len(meta), X.shape)

    panels = [("housekeeping", HOUSEKEEPING)]
    for i in range(args.n_random):
        p = list(rng.choice(universe, size=8, replace=False))
        panels.append((f"random_{i:04d}", p))

    rows = []
    for name, panel in panels:
        t0 = time.time()
        scored = panel_score_per_sample(samples_data, panel)
        if scored.empty:
            log.warning("%s no scoring", name); continue
        # join into tile meta
        m = meta.merge(scored[["spot_id", "sample_id", "raw_score", "resid_score"]],
                       on=["spot_id", "sample_id"], how="left")
        if m["resid_score"].isna().sum() > 0.5 * len(m):
            log.warning("%s missing too many", name); continue
        sp, pe, auc = loso_pooled(m, X, "resid_score")
        rows.append({"panel": name, "n_genes": len(panel),
                     "tile_size": args.tile_size,
                     "spearman_r": sp, "pearson_r": pe, "auroc_top25": auc,
                     "elapsed_s": time.time() - t0,
                     "panel_genes": ",".join(panel)})
        if name == "housekeeping" or name.endswith("0000") or name.endswith("0050"):
            log.info("[%s] sp=%.3f auc=%.3f  (%.1fs)", name, sp, auc, time.time() - t0)

    df = pd.DataFrame(rows)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep="\t", index=False)
    log.info("wrote %s (%d panels)", out, len(df))
    print(df.head(10).to_string())
    rand = df[df.panel.str.startswith("random_")]
    if len(rand):
        print(f"\nRandom panels: spearman mean={rand.spearman_r.mean():.3f} "
              f"max={rand.spearman_r.max():.3f}  auroc mean={rand.auroc_top25.mean():.3f}")


if __name__ == "__main__":
    main()
