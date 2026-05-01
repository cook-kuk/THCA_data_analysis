"""E-prompt: MSK enrichment bias documentation.

MSK-IMPACT thyroid cohort vs TCGA-THCA: histology, age, stage, mutation freq.
Output to project/results/audit_2026_04_29/msk_bias/.
"""
from __future__ import annotations
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT / "project/results/audit_2026_04_29/msk_bias"
OUT.mkdir(parents=True, exist_ok=True)

msk_sample = pd.read_csv(
    ROOT / "project/results/v17_tert_recovery/v3/v3_thyroid_mskcc_2016_clinical_sample.tsv",
    sep="\t", comment="#"
)
msk_patient = pd.read_csv(
    ROOT / "project/results/v17_tert_recovery/v3/v3_thyroid_mskcc_2016_clinical_patient.tsv",
    sep="\t", comment="#"
)
msk = msk_sample.merge(msk_patient, on="PATIENT_ID", how="left")

tcga = pd.read_csv(
    ROOT / "project/results/v17_tert_recovery/v2/sample_master_v17_tert_v2.tsv",
    sep="\t", low_memory=False
)
tcga = tcga[(tcga["dataset"] == "TCGA-THCA") & (tcga["normal_vs_tumor"] == "tumor")].copy()

# MSK histology (CANCER_TYPE_DETAILED): ATC, PDTC, etc.
msk_hist = msk["CANCER_TYPE_DETAILED"].fillna("Unknown").value_counts()
tcga_hist = tcga["histology_subtype"].fillna("Unknown").value_counts()

# Build harmonized 3-class histology buckets (PTC / PDTC / ATC / Other)
def harmonize_msk(s):
    s = str(s).lower()
    if "anaplastic" in s: return "ATC"
    if "poorly" in s or "pdtc" in s: return "PDTC"
    if "papillary" in s: return "PTC"
    return "Other"
def harmonize_tcga(s):
    s = str(s).lower()
    if "anaplastic" in s: return "ATC"
    if "poorly" in s or "pdtc" in s: return "PDTC"
    if "ptc" in s or "papillary" in s: return "PTC"
    return "Other"

msk["hist_h"] = msk["CANCER_TYPE_DETAILED"].apply(harmonize_msk)
tcga["hist_h"] = tcga["histology_subtype"].apply(harmonize_tcga)

ct = pd.DataFrame({
    "MSK-IMPACT (n)": msk["hist_h"].value_counts(),
    "TCGA-THCA (n)": tcga["hist_h"].value_counts(),
}).fillna(0).astype(int)
ct["MSK %"] = (100 * ct["MSK-IMPACT (n)"] / ct["MSK-IMPACT (n)"].sum()).round(1)
ct["TCGA %"] = (100 * ct["TCGA-THCA (n)"] / ct["TCGA-THCA (n)"].sum()).round(1)
ct = ct.reindex(["PTC", "PDTC", "ATC", "Other"]).dropna(how="all")
ct.to_csv(OUT / "histology_comparison.tsv", sep="\t")

# chi-square on PTC vs (PDTC+ATC)
agg = pd.DataFrame({
    "MSK": [(msk["hist_h"] == "PTC").sum(), msk["hist_h"].isin(["PDTC", "ATC"]).sum()],
    "TCGA": [(tcga["hist_h"] == "PTC").sum(), tcga["hist_h"].isin(["PDTC", "ATC"]).sum()],
}, index=["PTC", "PDTC+ATC"])
chi2, p_chi, _, _ = stats.chi2_contingency(agg.values)

# Age comparison
msk_age = pd.to_numeric(msk["AGE"], errors="coerce").dropna()
tcga_age = pd.to_numeric(tcga["age"], errors="coerce").dropna()
mw_age = stats.mannwhitneyu(msk_age, tcga_age, alternative="two-sided")

# Stage / advanced markers — use M_STAGE for MSK; clinical_stage for TCGA
msk_advanced = msk["M_STAGE"].fillna("").str.contains("M1", na=False).mean()
tcga_advanced = tcga["clinical_stage"].fillna("").str.contains("IV", case=False, na=False).mean()

summary = {
    "n_msk": int(len(msk)),
    "n_tcga": int(len(tcga)),
    "histology_pct_PTC_MSK_vs_TCGA": [float(ct.loc["PTC", "MSK %"]) if "PTC" in ct.index else None,
                                       float(ct.loc["PTC", "TCGA %"]) if "PTC" in ct.index else None],
    "histology_pct_ATC_MSK_vs_TCGA": [float(ct.loc["ATC", "MSK %"]) if "ATC" in ct.index else None,
                                       float(ct.loc["ATC", "TCGA %"]) if "ATC" in ct.index else None],
    "histology_pct_PDTC_MSK_vs_TCGA": [float(ct.loc["PDTC", "MSK %"]) if "PDTC" in ct.index else None,
                                        float(ct.loc["PDTC", "TCGA %"]) if "PDTC" in ct.index else None],
    "chi2_PTC_vs_aggressive": {"chi2": float(chi2), "p": float(p_chi)},
    "age_median_MSK_vs_TCGA": [float(msk_age.median()), float(tcga_age.median())],
    "age_mannwhitney_p": float(mw_age.pvalue),
    "M1_frac_MSK": round(float(msk_advanced), 3),
    "stageIV_frac_TCGA": round(float(tcga_advanced), 3),
}
(OUT / "msk_bias_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

# 2-panel figure
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
ax = axes[0]
x = np.arange(len(ct.index))
w = 0.35
ax.bar(x - w/2, ct["MSK %"], w, label=f"MSK-IMPACT (n={len(msk)})", color="#d62728")
ax.bar(x + w/2, ct["TCGA %"], w, label=f"TCGA-THCA (n={len(tcga)})", color="#1f77b4")
ax.set_xticks(x); ax.set_xticklabels(ct.index)
ax.set_ylabel("Histology composition (%)")
ax.set_title(f"A. Histology distribution (chi-square p={p_chi:.2e})")
ax.legend()

ax = axes[1]
ax.hist(msk_age, bins=30, alpha=0.5, label=f"MSK (median={msk_age.median():.1f})", color="#d62728", density=True)
ax.hist(tcga_age, bins=30, alpha=0.5, label=f"TCGA (median={tcga_age.median():.1f})", color="#1f77b4", density=True)
ax.set_xlabel("Age at diagnosis"); ax.set_ylabel("Density")
ax.set_title(f"B. Age distribution (Mann-Whitney p={mw_age.pvalue:.2e})")
ax.legend()

fig.suptitle("MSK-IMPACT thyroid vs TCGA-THCA: enrichment bias documentation",
             fontsize=12, fontweight="bold")
fig.tight_layout()
fig.savefig(OUT / "msk_bias_panel.pdf", bbox_inches="tight", dpi=200)
fig.savefig(OUT / "msk_bias_panel.png", bbox_inches="tight", dpi=150)
print(json.dumps(summary, indent=2))
print(f"\nFiles: {OUT}")
