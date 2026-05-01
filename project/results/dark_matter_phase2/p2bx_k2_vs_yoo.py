"""K2 cluster prediction (DM1/DM2) cross-validated against Yoo 2016 mutation/subtype classification.
Tests: does our 8-gene cluster prediction recover NBNR/RAS-like/BRAF-like from Yoo 2016?"""
from pathlib import Path
import re
import json

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import cohen_kappa_score, adjusted_rand_score

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project")
OUT = ROOT / "results/dark_matter_phase2"

# 1. Load K2 predictions + run-to-sample mapping
preds = pd.read_csv(ROOT / "results/v17_korean/K2_korean_predictions_v4.tsv", sep="\t")
runs = pd.read_csv(ROOT / "results/v17_korean/K1A_prjeb11591_runs.tsv", sep="\t")
yoo = pd.read_csv(OUT / "k2_yoo2016_mutations_parsed.tsv", sep="\t")

print(f"Predictions: {len(preds)} runs")
print(f"Run metadata: {len(runs)} entries")
print(f"Yoo 2016 mutation table: {len(yoo)} patients")

# 2. Map ERR/SRR run → SampleID
# sample_alias format: "SNU-GMI-FA02", "SNU-GMI-FA03-N" (N = matched normal)
def extract_sample_id(alias):
    if pd.isna(alias):
        return None, None
    m = re.match(r"SNU-GMI-([A-Z0-9]+)(-N)?$", str(alias))
    if not m:
        return None, None
    return m.group(1), bool(m.group(2))  # SampleID, is_normal

runs["yoo_sample_id"], runs["is_normal"] = zip(*runs["sample_alias"].map(extract_sample_id))
print(f"\nSample alias parse — examples:")
print(runs[["run_accession", "sample_alias", "yoo_sample_id", "is_normal"]].head(5).to_string(index=False))

# 3. Merge: predictions × runs × Yoo
preds = preds.merge(runs[["run_accession", "yoo_sample_id", "is_normal"]],
                    left_on="run", right_on="run_accession", how="left")
preds_tumor = preds[~preds["is_normal"].fillna(False)].copy()
print(f"\nTumor-only predictions: {len(preds_tumor)} runs ({preds_tumor['yoo_sample_id'].nunique()} unique samples)")

# Some samples have multiple runs — average the prediction prob
agg = preds_tumor.groupby("yoo_sample_id", observed=True).agg(
    p_DM2_mean=("p_DM2", "mean"),
    n_runs=("run", "count"),
).reset_index()
agg["DM_call"] = (agg["p_DM2_mean"] > 0.5).map({True: "DM2", False: "DM1"})

# Merge with Yoo mutations
merged = agg.merge(yoo[["SampleID", "Driver mutation", "Molecular subtype", "mol_subtype_label",
                        "Pathology", "is_dark_matter", "has_braf_v600e", "has_ras",
                        "has_dicer1", "has_eif1ax", "has_fusion"]],
                   left_on="yoo_sample_id", right_on="SampleID", how="inner")
print(f"\nMerged: {len(merged)} samples with both K2 prediction + Yoo mutation")

# 4. Cross-tab: K2 DM_call × Yoo molecular subtype
print("\n=== K2 DM call × Yoo molecular subtype ===")
ct = pd.crosstab(merged["DM_call"], merged["mol_subtype_label"])
print(ct)
print(ct.div(ct.sum(axis=0), axis=1).round(3))

# 5. Cross-tab: K2 DM_call × Yoo Pathology
print("\n=== K2 DM call × Pathology ===")
print(pd.crosstab(merged["DM_call"], merged["Pathology"]))

# 6. Cross-tab: K2 DM_call × Dark Matter status
print("\n=== K2 DM call × is_dark_matter ===")
print(pd.crosstab(merged["DM_call"], merged["is_dark_matter"]))

# 7. Within Dark Matter: K2 DM_call × DICER1/EIF1AX
dm_only = merged[merged["is_dark_matter"]].copy()
print(f"\n=== Within Dark Matter (K2 n={len(dm_only)}): K2 DM_call × DICER1/EIF1AX ===")
dm_only["alt_driver"] = dm_only["has_dicer1"] | dm_only["has_eif1ax"]
ct_alt = pd.crosstab(dm_only["DM_call"], dm_only["alt_driver"])
print(ct_alt)
if ct_alt.shape == (2, 2):
    odds, p = stats.fisher_exact(ct_alt.values.tolist())
    print(f"Fisher OR={odds:.3f}, p={p:.4f}")

# 8. Within Dark Matter: cluster vs pathology (cPTC vs FVPTC validation)
print(f"\n=== Within Dark Matter: K2 DM_call × Pathology (cPTC/FVPTC validation) ===")
dm_only_ptc = dm_only[dm_only["Pathology"].isin(["cPTC", "fvPTC"])].copy()
ct_path = pd.crosstab(dm_only_ptc["DM_call"], dm_only_ptc["Pathology"])
print(ct_path)
if ct_path.shape == (2, 2):
    odds, p = stats.fisher_exact(ct_path.values.tolist())
    print(f"Fisher OR={odds:.3f}, p={p:.4f}  (DM1 should be cPTC-enriched, DM2 fvPTC-enriched)")
    # Concordance: how often cPTC = DM1 AND fvPTC = DM2?
    correct = ct_path.loc["DM1", "cPTC"] + ct_path.loc["DM2", "fvPTC"] if "DM1" in ct_path.index and "DM2" in ct_path.index else 0
    total = ct_path.values.sum()
    print(f"Concordance (DM1=cPTC, DM2=fvPTC): {correct}/{total} = {correct/total*100:.1f}%")

# 9. NBNR concordance: K2 DM call × NBNR subtype
nbnr_in_dm1 = ((merged["mol_subtype_label"] == "NBNR") & (merged["DM_call"] == "DM1")).sum()
nbnr_in_dm2 = ((merged["mol_subtype_label"] == "NBNR") & (merged["DM_call"] == "DM2")).sum()
nbnr_total = (merged["mol_subtype_label"] == "NBNR").sum()
print(f"\n=== NBNR samples × K2 cluster ===")
print(f"NBNR total in merged: {nbnr_total}")
print(f"NBNR in DM1: {nbnr_in_dm1} | NBNR in DM2: {nbnr_in_dm2}")

# 10. Save
merged.to_csv(OUT / "p2bx_k2_vs_yoo_merged.tsv", sep="\t", index=False)

result = {
    "n_merged_samples": int(len(merged)),
    "K2_call_distribution": merged["DM_call"].value_counts().to_dict(),
    "K2call_x_Yoo_molSubtype": ct.to_dict(),
    "K2call_x_Yoo_pathology": pd.crosstab(merged["DM_call"], merged["Pathology"]).to_dict(),
    "K2call_x_DarkMatter": pd.crosstab(merged["DM_call"], merged["is_dark_matter"]).to_dict(),
    "within_DM_K2call_x_alt_driver": ct_alt.to_dict() if ct_alt.shape == (2,2) else None,
    "within_DM_K2call_x_pathology": ct_path.to_dict() if ct_path.shape == (2,2) else None,
    "DM1_cPTC_DM2_fvPTC_concordance": (correct, total, round(correct/total*100, 1)) if ct_path.shape == (2,2) else None,
}
(OUT / "p2bx_summary.json").write_text(json.dumps(result, indent=2, default=str))

# 11. Figure: cross-tab heatmaps
fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
import seaborn as sns
sns.heatmap(ct, annot=True, fmt="d", cmap="Blues", ax=axes[0], cbar=False)
axes[0].set_title("K2 DM call × Yoo molecular subtype\n(NBNR/RAS-like/BRAF-like)")
sns.heatmap(pd.crosstab(merged["DM_call"], merged["Pathology"]),
            annot=True, fmt="d", cmap="Greens", ax=axes[1], cbar=False)
axes[1].set_title("K2 DM call × Pathology")
if ct_path.shape == (2, 2):
    sns.heatmap(ct_path, annot=True, fmt="d", cmap="Oranges", ax=axes[2], cbar=False)
    axes[2].set_title(f"Within DM: K2 cluster × cPTC/fvPTC\nconcordance = {correct/total*100:.1f}%")
fig.suptitle(f"P2-Bx K2 cluster prediction × Yoo 2016 cross-validation (n={len(merged)})", y=1.02)
fig.tight_layout()
fig.savefig(OUT / "fig_p2a/p2bx_k2_vs_yoo.png", dpi=200, bbox_inches="tight")
print(f"\nSaved p2bx_summary.json + figure")
