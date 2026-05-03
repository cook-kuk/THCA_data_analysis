#!/usr/bin/env python3
"""Per-spot residualization of DM1_like / RAI_8 / TDS_like against log_counts + log_ngenes
within each sample, then re-run stage trend at sample-mean level.

Question: is the borderline epithelial trend a sequencing-depth artifact?
"""
from __future__ import annotations
import argparse, sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import LinearRegression

STAGE_ORDER = ["PT", "PTC", "LPTC", "ATC"]
STAGE_NUM = {s: i for i, s in enumerate(STAGE_ORDER)}
SCORES = ["RAI_8_score", "DM1_like_score", "TDS_like_score",
          "CAF_ECM_score", "EMT_score", "Hypoxia_score", "Proliferation_score"]


def residualize_within_sample(df: pd.DataFrame, target: str, predictors: list[str]) -> pd.Series:
    out = pd.Series(np.nan, index=df.index, dtype=float)
    for sid, sub in df.groupby("sample_id"):
        ok = sub[[target] + predictors].notna().all(axis=1)
        if ok.sum() < 10: continue
        X = sub.loc[ok, predictors].to_numpy()
        y = sub.loc[ok, target].to_numpy()
        m = LinearRegression().fit(X, y)
        resid = y - m.predict(X)
        out.loc[sub.index[ok]] = resid
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--combined", default="project/results/01_spatial_score/all_spots_scored.tsv.gz")
    ap.add_argument("--out", default="project/results/02_stage_trend/depth_corrected_trend.csv")
    args = ap.parse_args()

    df = pd.read_csv(args.combined, sep="\t")
    df = df[df["stage"].isin(STAGE_ORDER)].copy()
    df["log_counts"] = np.log1p(df["total_counts"])
    df["log_ngenes"] = np.log1p(df["n_genes_by_counts"])

    rows = []
    for sc in SCORES:
        if sc not in df.columns: continue
        df[f"{sc}_resid"] = residualize_within_sample(df, sc, ["log_counts", "log_ngenes"])

    epi_keep = df.groupby("sample_id", group_keys=False).apply(
        lambda g: g["Epithelial_score"] >= g["Epithelial_score"].quantile(0.5)
    )
    df["epi_keep"] = epi_keep.values

    for subset_name, mask in [("all_spots", np.ones(len(df), bool)),
                              ("epithelial_top50", df["epi_keep"].values)]:
        sub = df[mask]
        sm = sub.groupby(["sample_id", "stage"], as_index=False)[
            [f"{sc}_resid" for sc in SCORES if f"{sc}_resid" in sub.columns]
        ].mean()
        sm["stage_ord"] = sm["stage"].map(STAGE_NUM)
        for sc in SCORES:
            col = f"{sc}_resid"
            if col not in sm.columns: continue
            ok = sm[col].notna()
            if ok.sum() < 4: continue
            rho, p = spearmanr(sm.loc[ok, "stage_ord"], sm.loc[ok, col])
            rows.append({"score": sc, "subset": subset_name,
                         "test": "sample_mean_resid_spearman",
                         "n_samples": int(ok.sum()), "rho": rho, "p": p,
                         "note": "score residualized vs log_counts+log_ngenes within sample"})

    out = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)
    print(out.to_string(index=False))


if __name__ == "__main__":
    sys.exit(main())
