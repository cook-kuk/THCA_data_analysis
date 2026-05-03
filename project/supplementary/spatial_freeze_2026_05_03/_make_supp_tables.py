#!/usr/bin/env python3
"""Aggregate existing per-cohort outputs into manuscript-ready supp tables.
Pure aggregation — no new scoring."""
from pathlib import Path
import pandas as pd
import numpy as np

OUT = Path("project/supplementary/spatial_freeze_2026_05_03")

# ---------- SX: 28-sample sample-level summary ----------
# GSE250521 (n=16): aggregate sample means from spot-level scored table (raw within-sample z only)
g521 = pd.read_csv("project/results/01_spatial_score/all_spots_scored.tsv.gz", sep="\t")
g521_summary_rows = []
for sid, sub in g521.groupby("sample_id"):
    epi_thr_50 = sub["Epithelial_score"].quantile(0.5)
    epi_thr_25 = sub["Epithelial_score"].quantile(0.75)
    epi50 = sub[sub["Epithelial_score"] >= epi_thr_50]
    epi25 = sub[sub["Epithelial_score"] >= epi_thr_25]
    row = {"sample_id": sid, "dataset": "GSE250521",
           "condition": sub["stage"].iloc[0],  # PT/PTC/LPTC/ATC
           "n_spots": len(sub),
           "mean_RAI_8_score_raw_all":   sub["RAI_8_score"].mean(),
           "mean_RAI_8_score_raw_epi50": epi50["RAI_8_score"].mean(),
           "mean_RAI_8_score_raw_epi25": epi25["RAI_8_score"].mean(),
           "mean_DM1_like_score_raw_all":   sub["DM1_like_score"].mean(),
           "mean_DM1_like_score_raw_epi50": epi50["DM1_like_score"].mean(),
           "mean_DM1_like_score_raw_epi25": epi25["DM1_like_score"].mean(),
           "mean_TDS_overlap_score_raw_all": sub["TDS_like_score"].mean(),
           "mean_Epithelial_score_raw_all":  sub["Epithelial_score"].mean(),
           "mean_Proliferation_score_raw_all": sub["Proliferation_score"].mean(),
           "mean_THYROID_NONOVERLAP_score_raw_all":   np.nan,  # not scored in GSE250521 run
           "mean_THYROID_NONOVERLAP_score_resid_epi25": np.nan,
           "mean_DM1_like_score_resid_epi25": np.nan,  # depth-resid per-gene not stored at spot level
           "available_normalizations": "raw_within_sample_z; sample-level depth-corrected residual via depth_corrected_trend.py (epi-only, score-level)",
           "notes": "Cancer progression cohort; depth-corrected per-gene residualization not in original scoring run; refer to depth_corrected_trend.csv for sample-mean resid trend"}
    g521_summary_rows.append(row)
g521_summary = pd.DataFrame(g521_summary_rows)

# External 12 samples: rename + select columns to align with GSE250521 schema
ext = pd.read_csv("project_external_st/results/meta/sample_level_score_summary.tsv", sep="\t")
ext_keep = ["sample_id", "dataset", "condition", "n_spots",
            "mean_RAI_8_score_raw_all", "mean_RAI_8_score_raw_epi50", "mean_RAI_8_score_raw_epi25",
            "mean_DM1_like_score_raw_all", "mean_DM1_like_score_raw_epi50", "mean_DM1_like_score_raw_epi25",
            "mean_TDS_overlap_score_raw_all", "mean_Epithelial_score_raw_all",
            "mean_Proliferation_score_raw_all",
            "mean_THYROID_NONOVERLAP_score_raw_all",
            "mean_THYROID_NONOVERLAP_score_resid_epi25",
            "mean_DM1_like_score_resid_epi25"]
ext_sub = ext[[c for c in ext_keep if c in ext.columns]].copy()
ext_sub["available_normalizations"] = "raw_within_sample_z + per-gene depth-residualized z"
ext_sub["notes"] = ext_sub["dataset"].map({
    "GSE230424": "PTC + Hashimoto-overlap; per-patient HT vs PTC+HT split unavailable in GEO metadata; user spec=all PTC+HT",
    "GSE248205": "Autoimmune thyroid (no cancer); CONTROL/HT/GD; PMID 39003267"})

combined_28 = pd.concat([g521_summary, ext_sub], ignore_index=True, sort=False)
combined_28.to_csv(OUT / "SuppTable_SX_28sample_score_summary.tsv", sep="\t", index=False)
print("SX → 28 samples, %d cols → %s" % (combined_28.shape[1], OUT / "SuppTable_SX_28sample_score_summary.tsv"))
print(combined_28[["sample_id","dataset","condition","n_spots"]].to_string(index=False))

# ---------- SX+1: independent lineage + technical confounding (combined) ----------
iv = pd.read_csv("project_external_st/results/meta/independent_lineage_validation.tsv", sep="\t")
iv["table"] = "independent_lineage_validation"
tc = pd.read_csv("project_external_st/results/meta/technical_confounding_summary.tsv", sep="\t")
tc["table"] = "technical_confounding"
# rename tc cols to align: spearman_rho_with_log_counts → metric_value, p → p_value
tc_aligned = tc.rename(columns={"spearman_rho_with_log_counts": "spearman_rho",
                                "p": "spearman_p"})
tc_aligned["scope"] = "all_spots_vs_log_total_counts"
tc_aligned["version"] = tc_aligned["score"].apply(lambda s: "raw" if s.endswith("_raw") else "resid")
tc_aligned["interpretation"] = "Score ~ log_total_counts confound check; |ρ|>0.30 flagged"
tc_aligned["pearson_r"] = np.nan; tc_aligned["pearson_p"] = np.nan
keep_cols = ["table", "dataset", "score", "scope", "version", "n", "spearman_rho", "spearman_p",
             "pearson_r", "pearson_p", "flag", "interpretation"]
iv["score"] = "DM1_like vs THYROID_NONOVERLAP"
iv["flag"] = ""
combined = pd.concat([iv[[c for c in keep_cols if c in iv.columns]],
                      tc_aligned[[c for c in keep_cols if c in tc_aligned.columns]]],
                     ignore_index=True, sort=False)
combined.to_csv(OUT / "SuppTable_SX1_independent_lineage_and_technical_confounding.tsv",
                sep="\t", index=False)
print("SX+1 → %d rows, %d cols → %s" % (combined.shape[0], combined.shape[1],
                                          OUT / "SuppTable_SX1_independent_lineage_and_technical_confounding.tsv"))
