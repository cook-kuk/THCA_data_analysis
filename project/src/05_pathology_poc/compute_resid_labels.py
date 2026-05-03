#!/usr/bin/env python3
"""Compute depth-residualized labels per sample.

Reads:
  project/results/01_spatial_score/all_spots_scored.tsv.gz   (spot_id, sample_id, scores, QC)
  project/results/03_pathology_poc/tile_metadata.tsv.gz      (tile rows)

Writes:
  project/results/03_pathology_poc/tile_metadata_resid.tsv.gz
    + DM1_like_score_resid, RAI_8_score_resid, TDS_like_score_resid
    + log_counts, log_ngenes
    + DM1_like_score_globalresid (residual on pooled depth, sensitivity)

Residualization:
  per sample_id, OLS:  score ~ log(total_counts+1) + log(n_genes_by_counts+1)
  resid = score - predicted   (within-sample)
"""
from __future__ import annotations
import logging
from pathlib import Path
import numpy as np
import pandas as pd
import os
import statsmodels.api as sm

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("resid")

ROOT = Path(os.environ.get("THCA_ROOT", "/home/seungho/personal/THCA_data_analysis"))
SPOTS = ROOT / "project/results/01_spatial_score/all_spots_scored.tsv.gz"
META  = ROOT / "project/results/03_pathology_poc/tile_metadata.tsv.gz"
OUT   = ROOT / "project/results/03_pathology_poc/tile_metadata_resid.tsv.gz"

SCORES_TO_RESID = ["DM1_like_score", "RAI_8_score", "TDS_like_score"]


def residualize_per_sample(df: pd.DataFrame, score_col: str,
                           covars: list[str]) -> np.ndarray:
    out = np.full(len(df), np.nan)
    for sid, sub in df.groupby("sample_id"):
        y = sub[score_col].values.astype(float)
        X = sub[covars].values.astype(float)
        X = sm.add_constant(X, has_constant="add")
        try:
            res = sm.OLS(y, X, missing="drop").fit()
            pred = res.predict(X)
            out[sub.index] = y - pred
        except Exception as e:
            log.warning("[%s] %s residualize failed: %s", sid, score_col, e)
            out[sub.index] = y - np.nanmean(y)
    return out


def residualize_global(df: pd.DataFrame, score_col: str,
                       covars: list[str]) -> np.ndarray:
    y = df[score_col].values.astype(float)
    X = df[covars].values.astype(float)
    X = sm.add_constant(X, has_constant="add")
    res = sm.OLS(y, X, missing="drop").fit()
    pred = res.predict(X)
    return y - pred


def main():
    spots = pd.read_csv(SPOTS, sep="\t")
    log.info("spots: %d rows, %d samples", len(spots), spots.sample_id.nunique())
    qc_cols = ["spot_id", "sample_id", "total_counts", "n_genes_by_counts"]
    spots["log_counts"] = np.log1p(spots["total_counts"])
    spots["log_ngenes"] = np.log1p(spots["n_genes_by_counts"])

    full = spots.copy()
    full = full.reset_index(drop=True)
    for s in SCORES_TO_RESID:
        if s not in full.columns:
            log.warning("score %s missing in spots table", s); continue
        full[f"{s}_resid"]       = residualize_per_sample(full, s, ["log_counts", "log_ngenes"])
        full[f"{s}_globalresid"] = residualize_global(full, s, ["log_counts", "log_ngenes"])
        log.info("residualized %s", s)

    # merge into tile metadata
    meta = pd.read_csv(META, sep="\t")
    keep_cols = ["spot_id", "sample_id", "log_counts", "log_ngenes"] + \
                [f"{s}_resid" for s in SCORES_TO_RESID if f"{s}_resid" in full.columns] + \
                [f"{s}_globalresid" for s in SCORES_TO_RESID if f"{s}_globalresid" in full.columns]
    meta = meta.merge(full[keep_cols], on=["spot_id", "sample_id"], how="left")
    n_missing = meta["DM1_like_score_resid"].isna().sum()
    log.info("merged: %d tiles, missing DM1 resid: %d", len(meta), n_missing)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    meta.to_csv(OUT, sep="\t", index=False, compression="gzip")
    log.info("wrote %s", OUT)


if __name__ == "__main__":
    main()
