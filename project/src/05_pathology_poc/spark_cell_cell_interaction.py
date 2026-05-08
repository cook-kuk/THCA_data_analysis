#!/usr/bin/env python3
"""
Cell-cell interaction (CCI) analysis on GSE250521 Visium ST.

Two approaches, both at the spot grid:

A) Spatial cross-correlation (lag correlation)
   For each pair of modules (M_a, M_b), compute Pearson(M_a at spot i, M_b at
   neighbors of spot i, k=6). Symmetric: average of (a,b)+(b,a). This is a
   classical Moran-style cross-correlation. r > 0 → spatial co-niching.

B) Ligand-Receptor activity (CellChat-lite)
   Curated 12 L-R pairs spanning tumor-immune crosstalk:
     - PD-L1 (CD274) → PD-1 (PDCD1)
     - CXCL9 / CXCL10 / CXCL11 → CXCR3
     - CCL5 → CCR5
     - IFNG → IFNGR1
     - IL10 → IL10RA
     - TGFB1 → TGFBR1
     - CXCL13 → CXCR5  (TLS marker)
     - LGALS9 → HAVCR2 (TIM-3 axis)
     - CD80 → CD28
     - CD86 → CTLA4
   Per spot: ligand expression × neighbor mean receptor expression.
   Per slide: mean L-R activity, top 25% spot region count.

Outputs:
  spark_cci_lag_correlation_per_slide.tsv  (per slide × module pair)
  spark_cci_lag_correlation_pooled.tsv     (16-slide pooled per pair, with stage)
  spark_cci_LR_activity_per_slide.tsv      (per slide × L-R pair, mean activity)
  spark_cci_LR_activity_per_stage.tsv      (per-stage rollup)
"""
from __future__ import annotations
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import pearsonr, spearmanr
from itertools import combinations

ROOT = Path(__file__).resolve().parent.parent.parent.parent
GS = ROOT / "project/data/processed/GSE250521"
RES = ROOT / "project/results/03_pathology_poc"

MODULES = {
    "CD8_Teff":  ["CD8A", "CD8B", "GZMB", "GZMK", "PRF1", "IFNG", "NKG7", "CCL5"],
    "Bplasma":   ["CD19", "CD20", "MS4A1", "CD79A", "CD79B", "MZB1", "JCHAIN", "IGHM",
                  "IGHG1", "IGHA1", "IGKC", "IGLC2", "IGLC3"],
    "IFN_gamma": ["STAT1", "IRF1", "GBP1", "GBP4", "GBP5", "CXCL9", "CXCL10", "CXCL11",
                  "IFI44", "IFIT1", "ISG15", "MX1"],
    "M1_M2":     ["CD68", "CD163", "MRC1", "MARCO", "C1QA", "C1QB", "C1QC", "TYROBP", "AIF1"],
    "TLS":       ["CXCL13", "CCL19", "CCR7", "LTB", "LTA", "FCRL5", "FDCSP"],
    "RAI_lineage":["TPO", "DIO1", "TSHR", "PAX8", "TG", "FOXE1", "NKX2-1", "SLC5A5"],
    "EMT":       ["VIM", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "TWIST1", "FN1", "CDH2"],
    "CAF":       ["FAP", "ACTA2", "PDGFRA", "PDGFRB", "COL1A1", "COL3A1", "DCN", "POSTN"],
    "Endothelial":["PECAM1", "CDH5", "VWF", "KDR", "TEK", "ENG", "FLT1"],
}

LR_PAIRS = [
    ("PD-L1→PD-1",   ["CD274"], ["PDCD1"]),
    ("CXCL9→CXCR3",  ["CXCL9"], ["CXCR3"]),
    ("CXCL10→CXCR3", ["CXCL10"], ["CXCR3"]),
    ("CXCL11→CXCR3", ["CXCL11"], ["CXCR3"]),
    ("CCL5→CCR5",    ["CCL5"],  ["CCR5"]),
    ("IFNG→IFNGR1",  ["IFNG"],  ["IFNGR1"]),
    ("IL10→IL10RA",  ["IL10"],  ["IL10RA"]),
    ("TGFB1→TGFBR1", ["TGFB1"], ["TGFBR1"]),
    ("CXCL13→CXCR5", ["CXCL13"], ["CXCR5"]),
    ("LGALS9→HAVCR2", ["LGALS9"], ["HAVCR2"]),
    ("CD80→CD28",    ["CD80"],  ["CD28"]),
    ("CD86→CTLA4",   ["CD86"],  ["CTLA4"]),
]


def module_score(X, sym_idx, genes):
    cols = [sym_idx[g] for g in genes if g in sym_idx]
    if len(cols) < 2:
        return None
    sub = X[:, cols]
    z = (sub - sub.mean(axis=0)) / (sub.std(axis=0) + 1e-6)
    return z.mean(axis=1)


def main():
    sids = sorted([d.name for d in GS.iterdir() if (d / f"{d.name}.scored.h5ad").exists()])
    print(f"slides: {len(sids)}")

    lag_rows = []
    lr_rows = []

    for sid in sids:
        a = ad.read_h5ad(GS / sid / f"{sid}.scored.h5ad")
        stage = a.obs.stage.iloc[0] if "stage" in a.obs.columns else "?"
        if "spatial" not in a.obsm or a.n_obs < 100:
            continue
        coords = a.obsm["spatial"]
        X = a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X)
        sym_idx = {s: i for i, s in enumerate(a.var.index.values)}
        tree = cKDTree(coords)
        d, idx = tree.query(coords, k=7)
        nbr_idx = idx[:, 1:]  # (n, 6)

        # module scores
        mod = {}
        for name, genes in MODULES.items():
            ms = module_score(X, sym_idx, genes)
            if ms is not None:
                mod[name] = ms

        # A) spatial cross-correlation between module pairs
        names = sorted(mod.keys())
        for ma, mb in combinations(names, 2):
            va = mod[ma]; vb = mod[mb]
            vb_lag = vb[nbr_idx].mean(axis=1)  # neighbor mean of mb
            va_lag = va[nbr_idx].mean(axis=1)
            keep = ~(np.isnan(va) | np.isnan(vb_lag))
            r1 = float(pearsonr(va[keep], vb_lag[keep])[0]) if keep.sum() > 30 else float("nan")
            keep = ~(np.isnan(vb) | np.isnan(va_lag))
            r2 = float(pearsonr(vb[keep], va_lag[keep])[0]) if keep.sum() > 30 else float("nan")
            sym = float(np.nanmean([r1, r2]))
            # also same-spot correlation for context (no lag)
            keep = ~(np.isnan(va) | np.isnan(vb))
            r0 = float(pearsonr(va[keep], vb[keep])[0]) if keep.sum() > 30 else float("nan")
            lag_rows.append({"sample_id": sid, "stage": stage,
                             "mod_a": ma, "mod_b": mb,
                             "same_spot_r": r0,
                             "lag_a_to_neighbor_b": r1,
                             "lag_b_to_neighbor_a": r2,
                             "lag_symmetric": sym})

        # B) L-R activity per spot, then aggregate
        for lr_name, ligs, recs in LR_PAIRS:
            lcols = [sym_idx[g] for g in ligs if g in sym_idx]
            rcols = [sym_idx[g] for g in recs if g in sym_idx]
            if not lcols or not rcols:
                continue
            L = X[:, lcols].mean(axis=1)
            R = X[:, rcols].mean(axis=1)
            R_nbr = R[nbr_idx].mean(axis=1)  # neighbor receptor expression
            activity = L * R_nbr  # spot-level L-R score
            mean_a = float(np.nanmean(activity))
            top25 = float(np.nanmean(activity[activity > np.nanquantile(activity, 0.75)]))
            # per DM1 quartile
            dm1 = a.obs["DM1_like_score"].values
            q = pd.qcut(dm1, 4, labels=False, duplicates="drop")
            top_dm1 = float(np.nanmean(activity[q == 3])) if (q == 3).any() else float("nan")
            bot_dm1 = float(np.nanmean(activity[q == 0])) if (q == 0).any() else float("nan")
            lr_rows.append({"sample_id": sid, "stage": stage, "lr_pair": lr_name,
                            "n_spots": int(a.n_obs),
                            "mean_activity": mean_a, "top25_activity": top25,
                            "DM1top_activity": top_dm1, "DM1bot_activity": bot_dm1,
                            "DM1_top_minus_bot": top_dm1 - bot_dm1})
        print(f"  {sid} ({stage}): {a.n_obs} spots, {len(LR_PAIRS)} L-R, {len(list(combinations(names,2)))} module pairs")

    lagdf = pd.DataFrame(lag_rows)
    lagdf.to_csv(RES / "spark_cci_lag_correlation_per_slide.tsv", sep="\t", index=False)
    pool = lagdf.groupby(["mod_a", "mod_b"]).agg(
        n_slides=("sample_id", "count"),
        median_same_spot=("same_spot_r", "median"),
        median_lag_sym=("lag_symmetric", "median"),
        q10_lag=("lag_symmetric", lambda s: float(np.quantile(s.dropna(), 0.10)) if s.notna().any() else np.nan),
        q90_lag=("lag_symmetric", lambda s: float(np.quantile(s.dropna(), 0.90)) if s.notna().any() else np.nan),
    ).reset_index().sort_values("median_lag_sym", ascending=False)
    pool.to_csv(RES / "spark_cci_lag_correlation_pooled.tsv", sep="\t", index=False)

    lrdf = pd.DataFrame(lr_rows)
    lrdf.to_csv(RES / "spark_cci_LR_activity_per_slide.tsv", sep="\t", index=False)
    lrstg = lrdf.groupby(["stage", "lr_pair"]).agg(
        n_slides=("sample_id", "count"),
        mean=("mean_activity", "mean"),
        DM1_top_minus_bot_mean=("DM1_top_minus_bot", "mean"),
    ).reset_index().sort_values("DM1_top_minus_bot_mean", ascending=False)
    lrstg.to_csv(RES / "spark_cci_LR_activity_per_stage.tsv", sep="\t", index=False)

    print(f"\nwrote {(RES / 'spark_cci_lag_correlation_pooled.tsv').relative_to(ROOT)}")
    print("=== top 12 spatially co-niched module pairs (16-slide median lag-sym Pearson) ===")
    print(pool.head(12).to_string(index=False))
    print("\n=== bottom 8 (most negative — spatial avoidance) ===")
    print(pool.tail(8).to_string(index=False))

    print(f"\nwrote {(RES / 'spark_cci_LR_activity_per_stage.tsv').relative_to(ROOT)}")
    print("=== L-R activity DM1-top minus DM1-bot (per stage mean) ===")
    pivot = lrdf.groupby(["stage", "lr_pair"])["DM1_top_minus_bot"].mean().unstack(0).round(4)
    print(pivot)


if __name__ == "__main__":
    main()
