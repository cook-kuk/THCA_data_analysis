#!/usr/bin/env python3
"""
SPARK feature × hallmark module enrichment per-spot Visium.

For each Visium spot we have:
  - 8 SPARK morphology features (per its tile)
  - precomputed module scores in obs (Epithelial, Proliferation, EMT, CAF_ECM, Hypoxia)
  - we add 4 immune modules computed from raw expression (CD8 T, B-cell/plasma, IFN-γ, M1/M2)

For each module × SPARK feature pair, compute pooled (across slides) Spearman
+ per-stage Spearman. Identifies which morphology features track which biology.

Output:
  project/results/03_pathology_poc/spark_module_enrichment_pooled.tsv
  project/results/03_pathology_poc/spark_module_enrichment_per_stage.tsv
"""
from __future__ import annotations
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent.parent.parent.parent
GS = ROOT / "project/data/processed/GSE250521"
SPARK_TILE = ROOT / "project/results/03_pathology_poc/spark_analytical_per_tile.tsv.gz"
OUT_POOL = ROOT / "project/results/03_pathology_poc/spark_module_enrichment_pooled.tsv"
OUT_STAGE = ROOT / "project/results/03_pathology_poc/spark_module_enrichment_per_stage.tsv"

SPARK_FEATS = [
    "B5_thyrocyte_cluster_med", "B4_tumor_stroma_ratio",
    "C2_spatial_celltype_entropy", "n_thy", "n_lym", "n_fib",
]

# Marker sets (small, conservative; ENSG-agnostic — match by symbol)
MODULES = {
    "CD8_Teff":  ["CD8A", "CD8B", "GZMB", "GZMK", "PRF1", "IFNG", "NKG7", "CCL5"],
    "Bplasma":   ["CD19", "CD20", "MS4A1", "CD79A", "CD79B", "MZB1", "JCHAIN", "IGHM",
                  "IGHG1", "IGHA1", "IGKC", "IGLC2", "IGLC3"],
    "IFN_gamma": ["STAT1", "IRF1", "GBP1", "GBP4", "GBP5", "CXCL9", "CXCL10", "CXCL11",
                  "IFI44", "IFIT1", "IFIT2", "IFIT3", "ISG15", "MX1", "OAS1"],
    "M1_M2":     ["CD68", "CD163", "MRC1", "MARCO", "C1QA", "C1QB", "C1QC", "TYROBP", "AIF1"],
    "TLS":       ["CXCL13", "CCL19", "CCR7", "LTB", "LTA", "FCRL5", "FDCSP"],
}


def main():
    spark = pd.read_csv(SPARK_TILE, sep="\t")
    sids = sorted([d.name for d in GS.iterdir() if (d / f"{d.name}.scored.h5ad").exists()])
    print(f"slides: {len(sids)}")

    pooled = []
    per_stage = []
    for sid in sids:
        a = ad.read_h5ad(GS / sid / f"{sid}.scored.h5ad")
        a.obs.index.name = "spot_id"
        obs = a.obs.reset_index()
        spark_s = spark[spark.sample_id == sid]
        common = obs.merge(spark_s, on="spot_id", how="inner")
        if len(common) < 30:
            continue

        # gene index map
        var = a.var.copy(); var["sym"] = var.index.values
        sym_to_idx = {s: i for i, s in enumerate(var.sym.values)}
        X = a.X
        if hasattr(X, "toarray"):
            X = X.toarray()
        X = np.asarray(X)

        # spot mask aligned to common
        common_idx = a.obs.index.get_indexer(common.spot_id.values)
        valid = common_idx >= 0
        if not valid.all():
            common = common[valid].reset_index(drop=True)
            common_idx = common_idx[valid]
        Xc = X[common_idx]

        # compute module scores per common spot (z-score mean over genes present)
        module_scores = {}
        for mname, genes in MODULES.items():
            cols = [sym_to_idx[g] for g in genes if g in sym_to_idx]
            if len(cols) < 2:
                continue
            sub = Xc[:, cols]
            # z-score across spots within this slide
            z = (sub - sub.mean(axis=0)) / (sub.std(axis=0) + 1e-6)
            module_scores[mname] = z.mean(axis=1)

        # also bring obs-level precomputed scores
        for c in ["DM1_like_score", "RAI_8_score", "Proliferation_score", "EMT_score",
                  "CAF_ECM_score", "Hypoxia_score", "Epithelial_score"]:
            if c in common.columns:
                module_scores[c] = common[c].values

        stage = common.stage.iloc[0] if "stage" in common.columns else "?"
        for mname, m in module_scores.items():
            for f in SPARK_FEATS:
                if f not in common.columns: continue
                v = common[f].values
                keep = ~(np.isnan(m) | np.isnan(v)) & (np.std(v) > 0)
                if keep.sum() < 20: continue
                try:
                    r = spearmanr(m[keep], v[keep]).statistic
                except Exception:
                    continue
                if not np.isnan(r):
                    pooled.append({"sample_id": sid, "stage": stage, "module": mname,
                                   "feature": f, "n": int(keep.sum()), "spearman_r": float(r)})

    df = pd.DataFrame(pooled)
    print(f"per-(slide × module × feature) rows: {len(df)}")
    df.to_csv(OUT_POOL.with_name(OUT_POOL.name.replace("_pooled","_perslide")), sep="\t", index=False)

    # Pool: median Spearman per (module, feature) across slides
    pool = df.groupby(["module", "feature"]).agg(
        n_slides=("sample_id", "count"),
        median_r=("spearman_r", "median"),
        q10=("spearman_r", lambda s: float(np.quantile(s, 0.10))),
        q90=("spearman_r", lambda s: float(np.quantile(s, 0.90))),
    ).reset_index().sort_values("median_r", ascending=False)
    pool.to_csv(OUT_POOL, sep="\t", index=False)
    print(f"\nwrote {OUT_POOL.relative_to(ROOT)}")
    print("=== top 10 module × SPARK feature (pooled median Spearman) ===")
    print(pool.head(10).to_string(index=False))
    print("\n=== bottom 10 (most negative) ===")
    print(pool.tail(10).to_string(index=False))

    # Per-stage
    stg = df.groupby(["stage", "module", "feature"]).agg(
        n=("sample_id", "count"), median_r=("spearman_r", "median"),
    ).reset_index()
    stg.to_csv(OUT_STAGE, sep="\t", index=False)
    print(f"\nwrote {OUT_STAGE.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
