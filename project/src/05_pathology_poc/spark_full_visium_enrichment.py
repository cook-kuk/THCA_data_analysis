#!/usr/bin/env python3
"""
Full-Visium SPARK enrichment — ALL spots (≈ 46k), not just 3200 tile-anchored ones.

Inputs:
  project/data/processed/GSE250521/{sid}/{sid}.scored.h5ad

For each slide and each spot:
  1. Compute per-spot module score for 6 immune/tumor modules (z-score within slide)
  2. Spatial density features (NEAREST k=6 mean distance, Moran's I-lite)
  3. Per-spot DM1_like_score / RAI_8 / TDS already in .obs

Then:
  A) Pooled mixed-model-style: per-stage spot-level mean ± 95% CI for module × DM1 quartile
  B) Full per-slide per-spot Spearman: module × DM1 (all spots, not just 200)
  C) DM1 spatial autocorrelation per slide (Moran's I-lite k=6) using ALL spots
  D) Quartile-stratified: top-25% DM1 spots vs bottom-25% — module score difference (Cohen d)

Outputs:
  spark_full_visium_per_spot.tsv.gz
  spark_full_visium_per_slide.tsv
  spark_full_visium_quartile_stratified.tsv
"""
from __future__ import annotations
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent.parent.parent.parent
GS = ROOT / "project/data/processed/GSE250521"
RES = ROOT / "project/results/03_pathology_poc"
OUT_SPOT = RES / "spark_full_visium_per_spot.tsv.gz"
OUT_SLIDE = RES / "spark_full_visium_per_slide.tsv"
OUT_QSTRAT = RES / "spark_full_visium_quartile_stratified.tsv"

MODULES = {
    "CD8_Teff":  ["CD8A", "CD8B", "GZMB", "GZMK", "PRF1", "IFNG", "NKG7", "CCL5"],
    "Bplasma":   ["CD19", "CD20", "MS4A1", "CD79A", "CD79B", "MZB1", "JCHAIN", "IGHM",
                  "IGHG1", "IGHA1", "IGKC", "IGLC2", "IGLC3"],
    "IFN_gamma": ["STAT1", "IRF1", "GBP1", "GBP4", "GBP5", "CXCL9", "CXCL10", "CXCL11",
                  "IFI44", "IFIT1", "ISG15", "MX1"],
    "M1_M2":     ["CD68", "CD163", "MRC1", "MARCO", "C1QA", "C1QB", "C1QC", "TYROBP", "AIF1"],
    "TLS":       ["CXCL13", "CCL19", "CCR7", "LTB", "LTA", "FCRL5", "FDCSP"],
    "RAI_lineage":["TPO", "DIO1", "TSHR", "PAX8", "TG", "FOXE1", "NKX2-1", "SLC5A5"],
}


def cohen_d(a, b):
    if len(a) < 5 or len(b) < 5:
        return float("nan")
    sd_pool = np.sqrt(((len(a)-1)*np.var(a, ddof=1) + (len(b)-1)*np.var(b, ddof=1)) /
                       max(len(a)+len(b)-2, 1))
    return float((np.mean(a) - np.mean(b)) / max(sd_pool, 1e-6))


def main():
    sids = sorted([d.name for d in GS.iterdir() if (d / f"{d.name}.scored.h5ad").exists()])
    print(f"slides: {len(sids)}")

    spot_rows = []
    slide_rows = []
    qstrat_rows = []

    for sid in sids:
        a = ad.read_h5ad(GS / sid / f"{sid}.scored.h5ad")
        stage = a.obs.stage.iloc[0] if "stage" in a.obs.columns else "?"
        print(f"  {sid}: {a.n_obs} spots, stage {stage}")

        # gene matrix
        X = a.X
        if hasattr(X, "toarray"):
            X = X.toarray()
        X = np.asarray(X, dtype=np.float32)
        sym_to_idx = {s: i for i, s in enumerate(a.var.index.values)}

        # module z-score per spot within slide
        mod_z = {}
        for mname, genes in MODULES.items():
            cols = [sym_to_idx[g] for g in genes if g in sym_to_idx]
            if len(cols) < 2: continue
            sub = X[:, cols]
            z = (sub - sub.mean(axis=0)) / (sub.std(axis=0) + 1e-6)
            mod_z[mname] = z.mean(axis=1)

        # spatial features full
        coords = a.obsm["spatial"] if "spatial" in a.obsm else None
        spot_density = np.full(a.n_obs, np.nan)
        morans_dm1 = float("nan")
        morans_rai8 = float("nan")
        if coords is not None and len(coords) >= 10:
            tree = cKDTree(coords)
            d, idx = tree.query(coords, k=7)
            spot_density = 1.0 / np.maximum(d[:, 1:].mean(axis=1), 1e-3)
            try:
                v = a.obs["DM1_like_score"].values
                vn = (v - v.mean()) / (v.std() + 1e-6)
                lag = vn[idx[:, 1:]].mean(axis=1)
                morans_dm1 = float(np.corrcoef(vn, lag)[0, 1])
            except Exception:
                pass
            try:
                v2 = a.obs["RAI_8_score"].values
                vn2 = (v2 - v2.mean()) / (v2.std() + 1e-6)
                lag2 = vn2[idx[:, 1:]].mean(axis=1)
                morans_rai8 = float(np.corrcoef(vn2, lag2)[0, 1])
            except Exception:
                pass

        # per-slide rollup + quartile stratified
        v_dm1 = a.obs["DM1_like_score"].values
        q = pd.qcut(v_dm1, 4, labels=False, duplicates="drop")
        slide_metrics = {"sample_id": sid, "stage": stage, "n_spots": int(a.n_obs),
                         "morans_dm1_k6": morans_dm1, "morans_rai8_k6": morans_rai8}
        for mname, mz in mod_z.items():
            keep = ~(np.isnan(mz) | np.isnan(v_dm1))
            if keep.sum() >= 30:
                r = float(spearmanr(mz[keep], v_dm1[keep]).statistic)
                slide_metrics[f"spear_{mname}_x_DM1"] = r
            # top-25 vs bottom-25 quartile Cohen d for module
            top = mz[q == 3]
            bot = mz[q == 0]
            d_val = cohen_d(top, bot)
            slide_metrics[f"d_{mname}_top25vsbot25"] = d_val
            qstrat_rows.append({"sample_id": sid, "stage": stage, "module": mname,
                                "cohen_d_top25_vs_bot25": d_val,
                                "n_top": int((q == 3).sum()), "n_bot": int((q == 0).sum())})
        slide_rows.append(slide_metrics)

        # per-spot rows (subsampled to 5k per slide max for storage)
        n = a.n_obs
        idx_save = np.arange(n) if n <= 5000 else np.random.RandomState(42).choice(n, 5000, replace=False)
        for i in idx_save:
            row = {"sample_id": sid, "stage": stage,
                   "spot_idx": int(i), "DM1": float(v_dm1[i]),
                   "RAI8": float(a.obs["RAI_8_score"].iloc[i]),
                   "spot_density": float(spot_density[i]),
                   "DM1_quartile": int(q[i]) if not np.isnan(q[i]) else -1}
            for mname, mz in mod_z.items():
                row[mname] = float(mz[i])
            spot_rows.append(row)

    spot_df = pd.DataFrame(spot_rows)
    spot_df.to_csv(OUT_SPOT, sep="\t", index=False, compression="gzip")
    print(f"\nwrote {OUT_SPOT.relative_to(ROOT)}  ({len(spot_df)} spots subsampled)")

    sl = pd.DataFrame(slide_rows)
    sl.to_csv(OUT_SLIDE, sep="\t", index=False)
    print(f"wrote {OUT_SLIDE.relative_to(ROOT)}  ({len(sl)} slides)")

    q = pd.DataFrame(qstrat_rows)
    q.to_csv(OUT_QSTRAT, sep="\t", index=False)
    print(f"wrote {OUT_QSTRAT.relative_to(ROOT)}")

    # Headline summary: per-stage mean Spearman per module
    print("\n=== Per-stage mean Spearman (module × DM1_like_score) — FULL Visium ===")
    spear_cols = [c for c in sl.columns if c.startswith("spear_")]
    g = sl.groupby("stage")[spear_cols].mean().round(3)
    print(g.T)

    print("\n=== Per-stage mean Cohen d (top25 vs bot25 DM1 quartile) ===")
    d_cols = [c for c in sl.columns if c.startswith("d_")]
    gd = sl.groupby("stage")[d_cols].mean().round(3)
    print(gd.T)

    print("\n=== Per-slide Moran's I-lite k=6 (DM1 spatial autocorrelation) ===")
    print(sl[["sample_id", "stage", "n_spots", "morans_dm1_k6", "morans_rai8_k6"]].to_string(index=False))


if __name__ == "__main__":
    main()
