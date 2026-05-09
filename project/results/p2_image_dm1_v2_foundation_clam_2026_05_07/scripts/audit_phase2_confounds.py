#!/usr/bin/env python3
"""
Phase-2 audit: confounders, paper-killing slices, and provenance probes.

E1. Sex-stratified AUC (Female ≫ Male per existing report — quantify with
    bootstrap CI to know if Male AUC=0.35 is real or n=13 noise)
E2. TSS × molecular_subtype crosstab → does RAS-like concentrate in
    a few sites that share scanner/staining batch?
E3. Profile of the 3 RAS-like DM1 cases — TSS, age, stage, n_tiles,
    file size — vs the 13 RAS-like DM2 cases. Are they outliers?
E4. Within-cPTC RAS-like AUC — there are 18 (cPTC, RAS_like) cases in
    TCGA-THCA but only the 16 (FVPTC, RAS_like) made it into our 59.
    Could the model be classifying *follicular vs classical pattern*
    rather than DM1/DM2?
"""
from __future__ import annotations
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.linear_model import LogisticRegression

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/"
            "p2_image_dm1_v2_foundation_clam_2026_05_07")
META = Path("/home/seungho/personal/THCA_data_analysis/project/metadata/"
            "sample_master_v3.tsv")
DM_MASTER = Path("/home/seungho/personal/THCA_data_analysis/project/results/"
                 "dark_matter_phase1/tcga_dark_matter_master.tsv")
OUT = ROOT / "analysis_supp" / "audit_ras_auc100"
RNG = np.random.default_rng(42)
N_BOOT = 2000


def load_merged():
    pred = pd.read_csv(OUT / "merged_preds_subgroups.tsv", sep="\t")
    # add TSS (TCGA-XX-YYYY → XX is the TSS code)
    pred["tss"] = pred["submitter_id"].str.extract(r"TCGA-([^-]+)-")
    # add age, stage, sex from sample_master
    meta = pd.read_csv(META, sep="\t", low_memory=False)
    meta = meta[meta["dataset"] == "TCGA-THCA"].copy()
    meta["case_short"] = meta["sample_id"].str.extract(r"^(TCGA-[^-]+-[^-]+)")
    meta_case = (meta.sort_values("sample_id")
                     .drop_duplicates("case_short", keep="first")
                     [["case_short", "sex", "age", "ajcc_stage_group"]])
    pred = pred.merge(meta_case, left_on="submitter_id",
                      right_on="case_short", how="left")
    return pred


def bootstrap_auc(y, p, n_boot=N_BOOT):
    if len(set(y)) < 2 or len(y) < 2:
        return np.nan, np.nan, np.nan
    aucs = []
    n = len(y)
    for _ in range(n_boot):
        idx = RNG.integers(0, n, size=n)
        yb, pb = y[idx], p[idx]
        if len(set(yb)) < 2:
            continue
        aucs.append(roc_auc_score(yb, pb))
    if not aucs:
        return float(roc_auc_score(y, p)), np.nan, np.nan
    return (float(roc_auc_score(y, p)),
            float(np.percentile(aucs, 2.5)),
            float(np.percentile(aucs, 97.5)))


def e1_sex(df):
    print("\n=== E1. Sex-stratified AUC ===")
    rows = []
    for sex in ["Female", "Male"]:
        sub = df[df["sex"] == sex]
        y = sub["label"].values; p = sub["prob_DM1"].values
        auc, lo, hi = bootstrap_auc(y, p)
        rows.append({"group": f"sex={sex}", "n": len(sub),
                     "n_pos": int(y.sum()), "n_neg": int(len(y) - y.sum()),
                     "auc": auc, "ci_lo": lo, "ci_hi": hi})
        print(f"  {sex}: n={len(sub)} ({y.sum()}+/{int(len(y)-y.sum())}-) "
              f"AUC={auc:.3f} CI95=[{lo:.3f},{hi:.3f}]")
    return rows


def e2_tss_x_subtype(df):
    print("\n=== E2. TSS × molecular_subtype crosstab ===")
    ct = pd.crosstab(df["tss"], df["molecular_subtype"], margins=True,
                     margins_name="ALL")
    ct.to_csv(OUT / "e2_tss_x_subtype.tsv", sep="\t")
    # which TSS host the 16 RAS-like?
    ras = df[df["molecular_subtype"] == "RAS_like"]
    tss_ras = (ras.groupby("tss")
                  .agg(n=("label", "count"),
                       n_pos=("label", "sum"),
                       mean_prob=("prob_DM1", "mean"))
                  .reset_index())
    print("  TSS distribution among 16 RAS-like:")
    print(tss_ras.to_string(index=False))
    tss_ras.to_csv(OUT / "e2_tss_among_RASlike.tsv", sep="\t", index=False)

    # per-TSS within-RAS-like AUC (skip TSS with 0 var)
    rows = []
    for tss, sub in ras.groupby("tss"):
        if sub["label"].nunique() < 2:
            rows.append({"tss": tss, "n": len(sub),
                         "n_pos": int(sub["label"].sum()),
                         "auc": np.nan, "single_class": True})
        else:
            auc = roc_auc_score(sub["label"], sub["prob_DM1"])
            rows.append({"tss": tss, "n": len(sub),
                         "n_pos": int(sub["label"].sum()),
                         "auc": float(auc), "single_class": False})
    pd.DataFrame(rows).to_csv(OUT / "e2_per_tss_RASlike_auc.tsv",
                              sep="\t", index=False)
    return tss_ras.to_dict("records")


def e3_RASlike_DM1_profile(df):
    print("\n=== E3. Profile of 3 RAS-like DM1 vs 13 RAS-like DM2 ===")
    ras = df[df["molecular_subtype"] == "RAS_like"].copy()
    cols = ["submitter_id", "tss", "label", "prob_DM1", "n_tiles",
            "age", "ajcc_stage_group", "fold"]
    out = ras[cols].sort_values(["label", "prob_DM1"])
    out.to_csv(OUT / "e3_RASlike_full_profile.tsv", sep="\t", index=False)
    print(out.to_string(index=False))

    # are the 3 DM1 cases distinct from the 13 DM2?
    summ = (ras.groupby("label")
              .agg(n=("submitter_id", "count"),
                   mean_age=("age", "mean"),
                   tss_unique=("tss", "nunique"),
                   mean_n_tiles=("n_tiles", "mean"))
              .reset_index())
    print("\n  Group means:")
    print(summ.to_string(index=False))
    return out.to_dict("records")


def e4_within_subtype_crosshistology(df):
    """In TCGA-THCA, RAS-like and FVPTC are not 1:1 globally — there are 18
    (cPTC, RAS_like) and 11 (FVPTC, BRAF_like) cases in the full cohort,
    but our 59-slide subset doesn't contain the (cPTC, RAS_like) cells.
    Confirm and quantify the gap that prevents this disentanglement."""
    print("\n=== E4. Cross-histology × subtype gap ===")
    ct = pd.crosstab(df["histology_subtype"], df["molecular_subtype"],
                     margins=True, margins_name="ALL")
    print(ct.to_string())
    ct.to_csv(OUT / "e4_subset_histology_x_subtype.tsv", sep="\t")
    # count off-diagonal
    off = 0
    for h, m in [("cPTC", "RAS_like"), ("FVPTC", "BRAF_like")]:
        off += int(((df["histology_subtype"] == h)
                   & (df["molecular_subtype"] == m)).sum())
    print(f"  off-diagonal (cPTC×RAS_like + FVPTC×BRAF_like) = {off}")
    print("  → cannot disentangle 'histology pattern' from 'molecular subtype'")
    print("    in our 59-slide subset.")
    return {"off_diagonal_n": int(off),
            "crosstab": ct.to_dict()}


def e5_logreg_combined_baseline(df):
    """Baseline: predict label from (FVPTC indicator) AND (sex) AND (TSS dummy).
    If this gets close to CLAM's overall AUC=0.746, the image model isn't
    adding much beyond clinical covariates."""
    print("\n=== E5. Clinical-covariate-only baseline ===")
    sub = df.dropna(subset=["histology_subtype", "sex", "tss"]).copy()
    sub = sub[sub["histology_subtype"].isin(["cPTC", "FVPTC"])]
    X = pd.get_dummies(sub[["histology_subtype", "sex", "tss"]],
                       drop_first=True).astype(float).values
    y = sub["label"].values
    # k-fold to avoid in-sample inflation
    from sklearn.model_selection import KFold
    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    oof = np.full(len(y), np.nan)
    for tr, va in kf.split(X):
        lr = LogisticRegression(max_iter=1000, C=1.0)
        lr.fit(X[tr], y[tr])
        oof[va] = lr.predict_proba(X[va])[:, 1]
    auc = roc_auc_score(y, oof)
    print(f"  Clinical-only OOF LR AUC (histology+sex+TSS dummies): {auc:.3f}")
    print(f"  CLAM overall OOF AUC: 0.746")
    print(f"  CLAM image gain over clinical: +{0.746 - auc:.3f}")
    return {"clinical_only_oof_auc": float(auc),
            "clam_image_gain": float(0.746 - auc)}


def main():
    df = load_merged()
    print(f"[load] {len(df)} OOF preds with subgroup + clinical")

    summary = {}
    summary["E1_sex"] = e1_sex(df)
    summary["E2_tss"] = e2_tss_x_subtype(df)
    summary["E3_RASlike_DM1_profile"] = e3_RASlike_DM1_profile(df)
    summary["E4_subset_gap"] = e4_within_subtype_crosshistology(df)
    summary["E5_clinical_baseline"] = e5_logreg_combined_baseline(df)

    with (OUT / "PHASE2_AUDIT_SUMMARY.json").open("w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"\n[done] → {OUT}/PHASE2_AUDIT_SUMMARY.json")


if __name__ == "__main__":
    main()
