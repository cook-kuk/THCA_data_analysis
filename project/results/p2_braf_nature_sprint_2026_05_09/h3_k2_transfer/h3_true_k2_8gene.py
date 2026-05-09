"""H3-companion: True K2 (PRJEB11591) — only 8-gene quant available; check
within-cohort bimodality and HT-axis cannot be tested (HT-13 not quantified).
But verify thyroid-differentiation 8-gene split + relate to existing TCGA-trained
calls."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture
from scipy.stats import mannwhitneyu

OUT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p2_braf_nature_sprint_2026_05_09/h3_k2_transfer")

mat = pd.read_csv(
    "/data/thca/repo_results/v17_korean/K2_8gene_tpm_matrix_v4.tsv", sep="\t"
).set_index("run")
preds = pd.read_csv(
    "/data/thca/repo_results/v17_korean/K2_korean_predictions_v4.tsv", sep="\t"
).set_index("run")

# K1A run metadata
runs = pd.read_csv(
    "/data/thca/repo_results/v17_korean/K1A_prjeb11591_runs.tsv", sep="\t"
).set_index("run_accession")

# tag tissue type from sample_alias
def tag(name: str) -> str:
    if not isinstance(name, str):
        return "unknown"
    if "FA" in name and "-N" not in name:
        return "FA"
    if "FT" in name and "-N" not in name:
        return "FTC"
    if "FV" in name and "-N" not in name:
        return "FVPTC"
    if "PT" in name and "-N" not in name:
        return "PTC"
    if "-N" in name:
        return "Normal"
    return "unknown"

runs["tissue_tag"] = runs["sample_alias"].map(tag)
df = mat.join(preds[["p_DM2", "DM_call"]]).join(runs[["sample_alias", "tissue_tag"]])
df.to_csv(OUT / "h3_true_k2_8gene_summary.tsv", sep="\t")

# log2(tpm+1) within-cohort z then mean
genes = ["SLC5A5", "TPO", "TG", "TSHR", "PAX8", "NKX2-1", "FOXE1", "DIO1"]
log = np.log2(df[genes] + 1)
z = (log - log.mean()) / log.std(ddof=0)
df["panel_z_mean_8g"] = z.mean(axis=1)

# Bimodality on tumors only (PTC + FVPTC + FTC + FA, drop normal)
tumor_mask = df["tissue_tag"].isin(["PTC", "FVPTC", "FTC", "FA"])
xt = df.loc[tumor_mask, "panel_z_mean_8g"].dropna().to_numpy().reshape(-1, 1)
gm1 = GaussianMixture(1, random_state=0).fit(xt)
gm2 = GaussianMixture(2, random_state=0).fit(xt)
b1, b2 = gm1.bic(xt), gm2.bic(xt)

# Relate within-cohort 8g panel z to existing TCGA-trained DM_call
out: dict[str, object] = {
    "n_K2_total": int(len(df)),
    "tissue_breakdown": df["tissue_tag"].value_counts().to_dict(),
    "n_tumor": int(tumor_mask.sum()),
    "panel_z_mean_8g_bic_1comp": float(b1),
    "panel_z_mean_8g_bic_2comp": float(b2),
    "panel_z_mean_8g_bic_supports_bimodal_tumor": bool(b2 < b1),
    "panel_z_DM1_mean": float(df.loc[df["DM_call"] == "DM1", "panel_z_mean_8g"].mean()),
    "panel_z_DM2_mean": float(df.loc[df["DM_call"] == "DM2", "panel_z_mean_8g"].mean()),
    "TCGA_call_dist": df["DM_call"].value_counts().to_dict(),
    "limitation": "HT-13 panel (HLA-DR/DP/DQ, CD79A/B, MS4A1, AICDA, CXCL13, CCR6, IFNG) NOT quantified in this K2 mini-index. True K2 transfer requires full-transcriptome re-quant of 33 GB FASTQ.",
}

# Mann-Whitney within-cohort 8g vs DM1/DM2 (sanity)
a = df.loc[df["DM_call"] == "DM1", "panel_z_mean_8g"].dropna()
b = df.loc[df["DM_call"] == "DM2", "panel_z_mean_8g"].dropna()
if len(a) > 0 and len(b) > 0:
    u, p = mannwhitneyu(a, b, alternative="two-sided")
    out["MW_panel_z_DM1_vs_DM2_p"] = float(p)
    out["n_DM1"] = int(len(a))
    out["n_DM2"] = int(len(b))

with open(OUT / "h3_true_k2_8gene_metrics.json", "w") as fh:
    json.dump(out, fh, indent=2, default=str)
print(json.dumps(out, indent=2, default=str))
