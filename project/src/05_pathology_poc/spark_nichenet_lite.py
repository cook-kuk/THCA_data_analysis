#!/usr/bin/env python3
"""
NicheNet-lite: which CAF-niche ligands predict RAI-poor neighbor regions?

For each Visium spot in CAF-rich spots, identify ligand expression that
correlates with neighbor RAI-poor scoring. Output: top 20 candidate ligands
of CAF→RAI-thyrocyte signaling.

Approach (per slide, then pool):
  1. Define CAF-rich spots (top 25% CAF score)
  2. Extract their k=6 neighbor RAI score (lag)
  3. Spearman: ligand expression in CAF-rich spot ↔ −RAI_lag (avoidance)
  4. Top-ranked ligands = "active in CAF-rich, predict RAI-poor neighbors"
"""
from __future__ import annotations
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent.parent.parent.parent
GS = ROOT / "project/data/processed/GSE250521"
RES = ROOT / "project/results/03_pathology_poc"

# Curated thyroid-relevant ligand list (CAF/stromal/immune-secreted)
CANDIDATE_LIGANDS = [
    "TGFB1", "TGFB2", "TGFB3", "POSTN", "FN1", "COL1A1", "COL1A2", "COL3A1",
    "FGF2", "FGF7", "VEGFA", "VEGFC", "PDGFA", "PDGFB", "PDGFC",
    "WNT5A", "WNT5B", "DKK1", "DKK3", "SFRP2", "SFRP4",
    "CXCL12", "CXCL14", "CXCL13", "CCL2", "CCL19", "CCL21",
    "IL6", "IL11", "LIF", "OSM",
    "MMP2", "MMP9", "MMP11", "TIMP1", "BMP4", "BMP7", "GREM1",
    "SPP1", "TNC", "BGN", "DCN", "LUM",
]
CAF = ["FAP", "ACTA2", "PDGFRA", "PDGFRB", "COL1A1", "COL3A1", "DCN", "POSTN"]
RAI = ["TPO", "DIO1", "TSHR", "PAX8", "TG", "FOXE1", "NKX2-1", "SLC5A5"]


def module(X, sym, genes):
    cols = [sym[g] for g in genes if g in sym]
    if len(cols) < 2: return None
    sub = X[:, cols]
    return ((sub - sub.mean(0)) / (sub.std(0) + 1e-6)).mean(1)


def main():
    sids = sorted([d.name for d in GS.iterdir() if (d / f"{d.name}.scored.h5ad").exists()])
    rows = []
    for sid in sids:
        a = ad.read_h5ad(GS / sid / f"{sid}.scored.h5ad")
        if "spatial" not in a.obsm or a.n_obs < 100: continue
        coords = a.obsm["spatial"]
        X = a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X)
        sym = {s: i for i, s in enumerate(a.var.index.values)}
        caf = module(X, sym, CAF); rai = module(X, sym, RAI)
        if caf is None or rai is None: continue
        tree = cKDTree(coords)
        _, idx = tree.query(coords, k=7)
        idx = idx[:, 1:]
        # CAF-rich spots
        caf_thr = np.percentile(caf, 75)
        mask = caf > caf_thr
        if mask.sum() < 30: continue
        # neighbor RAI lag
        rai_lag_neg = -rai[idx].mean(axis=1)  # high = RAI-poor neighbor
        target = rai_lag_neg[mask]

        stage = a.obs.stage.iloc[0] if "stage" in a.obs.columns else "?"
        for lig in CANDIDATE_LIGANDS:
            if lig not in sym: continue
            lig_vals = X[:, sym[lig]][mask]
            if np.std(lig_vals) < 1e-6: continue
            try:
                r = float(spearmanr(lig_vals, target).statistic)
                rows.append({"sample_id": sid, "stage": stage, "ligand": lig,
                             "n_caf_rich_spots": int(mask.sum()),
                             "spearman": r})
            except Exception:
                continue

    df = pd.DataFrame(rows)
    df.to_csv(RES / "spark_nichenet_lite_per_slide.tsv", sep="\t", index=False)
    pool = df.groupby("ligand").agg(
        n_slides=("sample_id", "count"),
        med_r=("spearman", "median"),
        q25=("spearman", lambda s: float(np.quantile(s, 0.25))),
        q75=("spearman", lambda s: float(np.quantile(s, 0.75))),
    ).reset_index().sort_values("med_r", ascending=False)
    pool.to_csv(RES / "spark_nichenet_lite_pooled.tsv", sep="\t", index=False)

    print("=== Top ligands in CAF-rich spots predicting RAI-poor neighbors (16-slide median Spearman) ===")
    print(pool.head(15).to_string(index=False))
    # Per-stage top
    stg_pool = df.groupby(["stage", "ligand"]).agg(
        n_slides=("sample_id", "count"),
        med_r=("spearman", "median"),
    ).reset_index()
    print("\n=== Per-stage top ligand ===")
    for s in ["PT", "PTC", "LPTC", "ATC"]:
        sub = stg_pool[stg_pool.stage == s].sort_values("med_r", ascending=False).head(5)
        print(f"\n[{s}]")
        print(sub.to_string(index=False))


if __name__ == "__main__":
    main()
