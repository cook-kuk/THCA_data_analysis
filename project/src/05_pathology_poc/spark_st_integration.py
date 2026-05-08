#!/usr/bin/env python3
"""
SPARK × Spatial Transcriptomics 통합 — Visium spot × 8 SPARK features × 18k genes.

Inputs:
  project/data/processed/GSE250521/{sid}/{sid}.scored.h5ad   (Visium per-slide)
  project/results/03_pathology_poc/spark_analytical_per_tile.tsv.gz

Joins on (sample_id, spot_id).

Outputs:
  spark_st_joint_per_spot.tsv.gz             — 8521 spot × (8 SPARK feat + module scores + DM1/RAI/TDS)
  spark_st_gene_corr_top.tsv                 — top |Spearman| gene per SPARK feature (slide-level pooled)
  spark_st_feature_dm1_R2.tsv                — per-slide R² of SPARK features predicting DM1_like_score
"""
from __future__ import annotations
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import RidgeCV
from sklearn.metrics import r2_score
from sklearn.model_selection import LeaveOneGroupOut

ROOT = Path(__file__).resolve().parent.parent.parent.parent
GS = ROOT / "project/data/processed/GSE250521"
SPARK_TILE = ROOT / "project/results/03_pathology_poc/spark_analytical_per_tile.tsv.gz"
OUT_JOINT = ROOT / "project/results/03_pathology_poc/spark_st_joint_per_spot.tsv.gz"
OUT_GENECORR = ROOT / "project/results/03_pathology_poc/spark_st_gene_corr_top.tsv"
OUT_R2 = ROOT / "project/results/03_pathology_poc/spark_st_feature_dm1_R2.tsv"

SPARK_FEATS = [
    "B1_stromal_encased_thyrocyte_idx",
    "B3_nuclear_eccentricity_var",
    "B4_tumor_stroma_ratio",
    "C2_spatial_celltype_entropy",
    "B5_thyrocyte_cluster_med",
    "n_thy", "n_lym", "n_fib",
]


def main():
    spark = pd.read_csv(SPARK_TILE, sep="\t")
    print(f"SPARK tiles: {len(spark)}, slides {spark.sample_id.nunique()}")

    joint_rows = []
    feat_dm1_rows = []
    gene_corr_rows = []

    for sdir in sorted(GS.iterdir()):
        if not sdir.is_dir():
            continue
        scored = sdir / f"{sdir.name}.scored.h5ad"
        if not scored.exists():
            continue
        a = ad.read_h5ad(scored)
        sid = sdir.name
        spark_s = spark[spark.sample_id == sid]
        print(f"  {sid}: ST {a.n_obs} spots, SPARK {len(spark_s)} tiles")

        # Match Visium spot barcode to SPARK tile spot_id (same string)
        a.obs.index.name = "spot_id"
        obs = a.obs.reset_index()
        common = obs.merge(spark_s, on="spot_id", how="inner", suffixes=("", "_sp"))
        if len(common) == 0:
            # try matching on first 16 chars (Visium barcodes might have -1 suffix vs not)
            obs["spot_id_b"] = obs.spot_id.str.replace(r"-1$", "", regex=True)
            spark_s2 = spark_s.copy()
            spark_s2["spot_id_b"] = spark_s2.spot_id.str.replace(r"-1$", "", regex=True)
            common = obs.merge(spark_s2, on="spot_id_b", how="inner", suffixes=("", "_sp"))
            if len(common) == 0:
                print(f"    no spot_id overlap; skip")
                continue
        print(f"    matched {len(common)} spot×tile pairs")

        # Save joint per-spot rows
        keep_cols = ["spot_id", "stage", "DM1_like_score", "RAI_8_score", "TDS_like_score",
                     "Epithelial_score", "Proliferation_score", "EMT_score",
                     "CAF_ECM_score", "Hypoxia_score"]
        keep_cols = [c for c in keep_cols if c in common.columns]
        joint = common[["sample_id"] + keep_cols + SPARK_FEATS].copy()
        joint["sample_id"] = sid
        joint_rows.append(joint)

        # Per-slide: SPARK features → DM1_like_score residualized linear regression (R²)
        valid = common.dropna(subset=SPARK_FEATS + ["DM1_like_score"])
        if len(valid) >= 30:
            X = valid[SPARK_FEATS].values
            y = valid["DM1_like_score"].values
            mdl = RidgeCV(alphas=[0.1, 1, 10]).fit(X, y)
            pred = mdl.predict(X)
            r2 = r2_score(y, pred)
            feat_dm1_rows.append({"sample_id": sid, "n_spots": len(valid),
                                  "stage": valid.stage.iloc[0] if "stage" in valid.columns else None,
                                  "ridge_R2": r2,
                                  "spearman_pred_vs_obs": float(spearmanr(y, pred).statistic)})

        # Per-slide top-gene Spearman per SPARK feature (slow → top 500 high-variance genes only)
        if a.n_vars > 0 and a.X is not None and len(common) >= 30:
            X = a.X
            if hasattr(X, "toarray"):
                X = X.toarray()
            else:
                X = np.asarray(X)
            # filter to common spots (preserve a.obs order)
            idx = a.obs.index.isin(common.spot_id)
            Xc = X[idx]
            obs_aligned = a.obs[idx]
            common_aligned = common.set_index("spot_id").loc[obs_aligned.index, SPARK_FEATS].values
            # high-var genes
            var = Xc.var(axis=0)
            top500 = np.argsort(var)[::-1][:500]
            gene_names = a.var.index.values[top500]
            for i, feat in enumerate(SPARK_FEATS):
                y = common_aligned[:, i]
                if np.isnan(y).all() or np.var(y) < 1e-6:
                    continue
                rs = []
                for j, g in enumerate(top500):
                    x = Xc[:, g]
                    keep = ~(np.isnan(x) | np.isnan(y))
                    if keep.sum() < 20:
                        continue
                    try:
                        r = spearmanr(x[keep], y[keep]).statistic
                        if not np.isnan(r):
                            rs.append((float(r), gene_names[j]))
                    except Exception:
                        continue
                if rs:
                    rs.sort(key=lambda t: -abs(t[0]))
                    for r, g in rs[:5]:
                        gene_corr_rows.append({"sample_id": sid, "feature": feat, "gene": g, "spearman": r})

    if joint_rows:
        df = pd.concat(joint_rows, ignore_index=True)
        df.to_csv(OUT_JOINT, sep="\t", index=False, compression="gzip")
        print(f"\nwrote {OUT_JOINT.relative_to(ROOT)}  ({len(df)} spot×tile pairs)")
    if feat_dm1_rows:
        pd.DataFrame(feat_dm1_rows).to_csv(OUT_R2, sep="\t", index=False)
        print(f"wrote {OUT_R2.relative_to(ROOT)}")
    if gene_corr_rows:
        gc = pd.DataFrame(gene_corr_rows)
        gc.to_csv(OUT_GENECORR, sep="\t", index=False)
        print(f"wrote {OUT_GENECORR.relative_to(ROOT)}  ({len(gc)} rows)")


if __name__ == "__main__":
    main()
