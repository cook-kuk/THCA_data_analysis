"""H-prompt: FFPE QC vs Fresh Frozen in GSE213647 (within-study comparator).

Compares panel_z (8-gene signature score) and library size between FFPE (n=80,
TruSeq RNA access kit) and Fresh-Frozen (n=169, same TruSeq kit) — controls kit
effect, isolates fixation effect.

Yoo 262 FFPE NGS as referenced in 2026-04-29 AM meeting: GSE213647 contains 80
FFPE samples (not 262). The 262 number may include additional unpublished
samples not yet in GEO; only the 80 FFPE in GEO are analyzable here.
"""
from __future__ import annotations
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT / "project/results/audit_2026_04_29/ffpe_qc"
OUT.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(ROOT / "project/results/v17_korean/GSE213647_panel_score.tsv", sep="\t")
df["panel_z"] = pd.to_numeric(df["panel_z"], errors="coerce")
df["lib_size"] = pd.to_numeric(df["lib_size"], errors="coerce")
df["age"] = pd.to_numeric(df["age"], errors="coerce")

ffpe = df[(df["fixation"] == "FFPE") & (df["library_kit"].str.contains("TruSeq", na=False))].copy()
ff_truseq = df[(df["fixation"] == "Fresh Frozen") & (df["library_kit"].str.contains("TruSeq", na=False))].copy()
ff_stranded = df[(df["fixation"] == "Fresh Frozen") & (df["library_kit"].str.contains("stranded", na=False))].copy()

# Histology distribution control
hist_ffpe = ffpe["histology"].value_counts(normalize=True).to_dict()
hist_ff = ff_truseq["histology"].value_counts(normalize=True).to_dict()

# panel_z comparison (within tumor only — exclude normals if heterogeneous)
ffpe_tumor = ffpe[ffpe["tissue_type"] != "Normal"]
ff_tumor = ff_truseq[ff_truseq["tissue_type"] != "Normal"]

mw_z = stats.mannwhitneyu(ffpe_tumor["panel_z"].dropna(), ff_tumor["panel_z"].dropna(), alternative="two-sided")
ks_z = stats.ks_2samp(ffpe_tumor["panel_z"].dropna(), ff_tumor["panel_z"].dropna())
mw_lib = stats.mannwhitneyu(ffpe["lib_size"].dropna(), ff_truseq["lib_size"].dropna(), alternative="two-sided")

summary = {
    "n_FFPE_TruSeq": int(len(ffpe)),
    "n_FF_TruSeq": int(len(ff_truseq)),
    "n_FF_stranded": int(len(ff_stranded)),
    "panel_z_FFPE_tumor": {"n": int(len(ffpe_tumor.dropna(subset=["panel_z"]))), "median": float(ffpe_tumor["panel_z"].median()), "iqr": [float(ffpe_tumor["panel_z"].quantile(0.25)), float(ffpe_tumor["panel_z"].quantile(0.75))]},
    "panel_z_FF_tumor":   {"n": int(len(ff_tumor.dropna(subset=["panel_z"]))), "median": float(ff_tumor["panel_z"].median()), "iqr": [float(ff_tumor["panel_z"].quantile(0.25)), float(ff_tumor["panel_z"].quantile(0.75))]},
    "panel_z_MW_p": float(mw_z.pvalue),
    "panel_z_KS_p": float(ks_z.pvalue),
    "panel_z_KS_stat": float(ks_z.statistic),
    "lib_size_median_FFPE": float(ffpe["lib_size"].median()),
    "lib_size_median_FF": float(ff_truseq["lib_size"].median()),
    "lib_size_MW_p": float(mw_lib.pvalue),
    "histology_pct_FFPE": {k: round(100*v, 1) for k, v in hist_ffpe.items()},
    "histology_pct_FF_TruSeq": {k: round(100*v, 1) for k, v in hist_ff.items()},
    "interpretation": (
        "If panel_z KS p > 0.05 between FFPE-TruSeq and FF-TruSeq (same kit, "
        "different fixation), the 8-gene signature is robust to formalin "
        "fixation degradation. If KS p < 0.05, FFPE introduces a detectable "
        "shift that needs ComBat-Seq adjustment before pooled analysis."
    ),
}
(OUT / "ffpe_qc_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

fig, axes = plt.subplots(1, 3, figsize=(16, 5))
ax = axes[0]
ax.boxplot([ffpe_tumor["panel_z"].dropna(), ff_tumor["panel_z"].dropna()], labels=[f"FFPE\nTruSeq\nn={len(ffpe_tumor)}", f"Fresh Frozen\nTruSeq\nn={len(ff_tumor)}"])
ax.set_ylabel("panel_z (8-gene signature score)")
ax.set_title(f"A. Tumor panel_z by fixation (KS p={ks_z.pvalue:.3f})")
ax.axhline(0, ls="--", color="gray", alpha=0.5)

ax = axes[1]
for grp, lbl, c in [(ffpe, f"FFPE-TruSeq (n={len(ffpe)})", "#d62728"), (ff_truseq, f"FF-TruSeq (n={len(ff_truseq)})", "#1f77b4"), (ff_stranded, f"FF-stranded (n={len(ff_stranded)})", "#7f7f7f")]:
    ls = grp["lib_size"].dropna() / 1e6
    ax.hist(ls, bins=30, alpha=0.5, label=lbl, color=c, density=True)
ax.set_xlabel("Library size (millions of reads)")
ax.set_ylabel("Density")
ax.set_title("B. Library size distribution by fixation × kit")
ax.legend(fontsize=8)

ax = axes[2]
all_hists = sorted(set(list(hist_ffpe.keys()) + list(hist_ff.keys())))
ffpe_pct = [100*hist_ffpe.get(h, 0) for h in all_hists]
ff_pct = [100*hist_ff.get(h, 0) for h in all_hists]
x = np.arange(len(all_hists)); w = 0.4
ax.bar(x - w/2, ffpe_pct, w, label="FFPE-TruSeq", color="#d62728")
ax.bar(x + w/2, ff_pct, w, label="FF-TruSeq", color="#1f77b4")
ax.set_xticks(x); ax.set_xticklabels(all_hists, rotation=30, ha="right")
ax.set_ylabel("Histology fraction (%)")
ax.set_title("C. Histology composition (must be similar to interpret panel_z)")
ax.legend()

fig.suptitle("FFPE-vs-Fresh-Frozen QC within GSE213647 (TruSeq RNA access kit)", fontsize=12, fontweight="bold")
fig.tight_layout()
fig.savefig(OUT / "figure_ffpe_qc.pdf", bbox_inches="tight", dpi=200)
fig.savefig(OUT / "figure_ffpe_qc.png", bbox_inches="tight", dpi=150)
print(json.dumps(summary, indent=2))
