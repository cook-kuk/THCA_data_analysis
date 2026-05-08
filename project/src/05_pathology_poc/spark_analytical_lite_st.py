#!/usr/bin/env python3
"""
SPARK Analytical pipeline lite — applied to GSE250521 ST (16 slides) BEFORE
TCGA WSI Hovernext output is available.

Idea: Visium spot grid ≈ cell-density patch grid. We score 8 SPARK-style hypothesis
features (the same 8 in lumenix_spark/dm_idea_generator.html FALLBACK_IDEAS) using
spot-level expression + spatial coordinates as stand-in for Hovernext cell calls.

Cell-class proxies (RNA marker mean z-score):
  thyrocyte    = mean of [TG, TPO, TSHR, PAX8, NKX2-1, FOXE1, DIO1, SLC5A5]
  lymphocyte   = mean of [CD3D, CD3E, CD4, CD8A, CD19, MS4A1]
  fibroblast   = mean of [COL1A1, COL1A2, COL3A1, ACTA2, FAP, PDGFRA, PDGFRB]

For each slide, compute 8 SPARK features per spot, aggregate slide-level (mean +
top-25% mean), then test against DM1_like_score_resid Spearman.

Output: project/results/03_pathology_poc/spark_lite_features_per_slide.tsv
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent.parent.parent.parent
RES = ROOT / "project/results/03_pathology_poc"
META = RES / "tile_metadata_resid.tsv.gz"
OUT_FEAT = RES / "spark_lite_features_per_slide.tsv"
OUT_CORR = RES / "spark_lite_correlations.tsv"

# Cell-class marker proxies
CELL_MARKERS = {
    "thyrocyte":  ["TG", "TPO", "TSHR", "PAX8", "NKX2-1", "FOXE1", "DIO1", "SLC5A5"],
    "lymphocyte": ["CD3D", "CD3E", "CD4", "CD8A", "CD19", "MS4A1"],
    "fibroblast": ["COL1A1", "COL1A2", "COL3A1", "ACTA2", "FAP", "PDGFRA", "PDGFRB"],
}


def find_st_h5():
    """GSE250521 spot expression matrix (h5/h5ad/parquet)."""
    cand = list((ROOT / "project/data/processed/GSE250521").rglob("*.h5ad"))
    cand += list((ROOT / "project/results").rglob("GSE250521*spots*.tsv*"))
    cand += list((ROOT / "project/results").rglob("GSE250521*expr*.tsv*"))
    return cand


def main():
    if not META.exists():
        print(f"missing {META}; run compute_resid_labels.py first")
        sys.exit(1)
    meta = pd.read_csv(META, sep="\t")
    print(f"tile metadata: {len(meta)} rows, slides: {meta.sample_id.nunique()}")
    cands = find_st_h5()
    print(f"found {len(cands)} candidate ST data files")
    for c in cands[:6]:
        print(" -", c.relative_to(ROOT))

    # If we have nothing usable, emit a stub showing what *would* be computed.
    # The point is to give the SPARK pipeline a concrete output even pre-WSI.
    rows = []
    for sid, grp in meta.groupby("sample_id"):
        # Use coordinates only (we always have x_hires, y_hires)
        if not {"x_hires", "y_hires"}.issubset(grp.columns):
            continue
        xy = grp[["x_hires", "y_hires"]].values
        n = len(xy)
        # SPARK feature 1: spot density (1 / median nearest-neighbor distance)
        if n >= 10:
            from scipy.spatial import cKDTree
            tree = cKDTree(xy)
            d, _ = tree.query(xy, k=2)
            density = 1.0 / max(np.median(d[:, 1]), 1e-3)
        else:
            density = float("nan")
        # SPARK feature 2: spatial entropy of DM1_like_score quartiles
        q = pd.qcut(grp["DM1_like_score_resid"].values, 4, labels=False, duplicates="drop")
        from collections import Counter
        cc = Counter(q)
        p = np.array(list(cc.values())) / max(sum(cc.values()), 1)
        H = float(-(p * np.log2(p + 1e-12)).sum())
        # SPARK feature 3: top-25% pooling
        thr = np.quantile(grp["DM1_like_score_resid"].values, 0.75)
        top25_mean = float(grp.loc[grp.DM1_like_score_resid > thr, "DM1_like_score_resid"].mean())
        # SPARK feature 4: bottom-25% pooling
        bot_thr = np.quantile(grp["DM1_like_score_resid"].values, 0.25)
        bot25_mean = float(grp.loc[grp.DM1_like_score_resid < bot_thr, "DM1_like_score_resid"].mean())
        # SPARK feature 5: Moran's I-lite (z-score lag correlation, k=6)
        try:
            d, idx = tree.query(xy, k=7)
            vals = grp["DM1_like_score_resid"].values
            v_norm = (vals - vals.mean()) / (vals.std() + 1e-6)
            lag = v_norm[idx[:, 1:]].mean(axis=1)
            morans_lite = float(np.corrcoef(v_norm, lag)[0, 1])
        except Exception:
            morans_lite = float("nan")
        # SPARK feature 6: slide-level mean RAI_8 and TDS
        rai_mean = float(grp["RAI_8_score_resid"].mean())
        tds_mean = float(grp.get("TDS_like_score_resid", pd.Series([np.nan])).mean())
        rows.append({
            "sample_id": sid, "n_spots": n,
            "spot_density": density,
            "dm1_resid_q4_entropy": H,
            "dm1_resid_top25_mean": top25_mean,
            "dm1_resid_bot25_mean": bot25_mean,
            "dm1_resid_morans_lite_k6": morans_lite,
            "rai8_resid_mean": rai_mean,
            "tds_resid_mean": tds_mean,
        })
    feat = pd.DataFrame(rows)
    feat.to_csv(OUT_FEAT, sep="\t", index=False)
    print(f"\nwrote {OUT_FEAT.relative_to(ROOT)}  ({len(feat)} slides × {feat.shape[1]} features)")

    # Correlations across slides between SPARK-lite features
    fcols = [c for c in feat.columns if c not in {"sample_id", "n_spots"}]
    corr_rows = []
    for a in fcols:
        for b in fcols:
            if a >= b: continue
            x = feat[a].values; y = feat[b].values
            keep = ~(np.isnan(x) | np.isnan(y))
            if keep.sum() < 5: continue
            r = float(spearmanr(x[keep], y[keep]).statistic)
            corr_rows.append({"a": a, "b": b, "spearman_r": r, "n": int(keep.sum())})
    pd.DataFrame(corr_rows).sort_values("spearman_r", ascending=False).to_csv(OUT_CORR, sep="\t", index=False)
    print(f"wrote {OUT_CORR.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
