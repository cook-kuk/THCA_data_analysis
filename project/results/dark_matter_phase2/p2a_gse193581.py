"""P2-A: External sc validation on GSE193581 (Lu 2023 JCI, 23 samples / 6 PTC + 7 NORM + 10 ATC).

Goal: replicate Phase 1 finding r=0.905 (8-gene ↔ FVPTC) at sc level in PTC malignant cells.
Decision: pooled r > 0.7 + median patient r > 0.7 → Nat Comm reach.
"""
from pathlib import Path
import json

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path("/home/seungho/personal/THCA_data_analysis/project/results/dark_matter_phase2")
FIG = OUT / "fig_p2a"
FIG.mkdir(parents=True, exist_ok=True)
H5AD = "/home/seungho/personal/THCA_data_analysis/project/results/v17_lu2023/GSE193581_hvg_adata.h5ad"

a = ad.read_h5ad(H5AD)
print(f"Loaded: {a.shape}")

# Sample-level histology mapping
print("Per-sample histology:")
print(a.obs[["sample", "histology"]].drop_duplicates().sort_values("histology"))

# Score FVPTC and cPTC signatures (genes available in HVG)
FVPTC_GENES = ["TG", "TPO", "TSHR", "DIO2"]  # subset present
CPTC_GENES = ["KRT19", "TIMP1", "FN1", "CITED1"]
sc.tl.score_genes(a, gene_list=FVPTC_GENES, score_name="score_fvptc")
sc.tl.score_genes(a, gene_list=CPTC_GENES, score_name="score_cptc")
print(f"\nScored: FVPTC ({FVPTC_GENES}), cPTC ({CPTC_GENES})")

# Note: DM_score is precomputed by v17 work using full 8-gene panel from raw data
# Validate it's available
assert "DM_score" in a.obs.columns
print(f"DM_score range: {a.obs['DM_score'].min():.3f} to {a.obs['DM_score'].max():.3f}")

# Filter to PTC samples × Malignant cells (paper-equivalent of "tumor thyrocytes")
ptc_mask = (a.obs["histology"] == "PTC") & (a.obs["author_celltype"] == "Malignant cell")
ptc = a[ptc_mask].copy()
print(f"\nPTC malignant cells: {ptc.shape}")
print(f"Patients with PTC malignant cells:")
print(ptc.obs.groupby("sample", observed=True).size())

# Pooled correlation
r_pool, p_pool = stats.pearsonr(ptc.obs["DM_score"], ptc.obs["score_fvptc"])
print(f"\n=== POOLED Pearson r (DM_score, FVPTC) on PTC malignant cells ===")
print(f"  n = {len(ptc)}, r = {r_pool:.3f}, p = {p_pool:.2e}")

# Per-patient correlation
per_pt = []
for s, sub in ptc.obs.groupby("sample", observed=True):
    if len(sub) >= 30:
        r, p = stats.pearsonr(sub["DM_score"], sub["score_fvptc"])
        per_pt.append({"sample": s, "n_cells": len(sub), "r": r, "p": p})
pp = pd.DataFrame(per_pt).sort_values("r", ascending=False)
print(f"\n=== PER-PATIENT correlations (≥30 cells) ===")
print(pp.to_string(index=False))
print(f"\nMedian r: {pp['r'].median():.3f}")
print(f"IQR: [{pp['r'].quantile(0.25):.3f}, {pp['r'].quantile(0.75):.3f}]")
print(f"# patients with r > 0.7: {(pp['r'] > 0.7).sum()} / {len(pp)}")
print(f"# patients with r > 0.5: {(pp['r'] > 0.5).sum()} / {len(pp)}")

# Bootstrap 95% CI for pooled r
n_boot = 1000
rng = np.random.default_rng(42)
boot_r = []
n_p = len(ptc)
x = ptc.obs["DM_score"].values
y = ptc.obs["score_fvptc"].values
for _ in range(n_boot):
    idx = rng.integers(0, n_p, n_p)
    boot_r.append(stats.pearsonr(x[idx], y[idx])[0])
ci_low, ci_high = np.percentile(boot_r, [2.5, 97.5])
print(f"\nBootstrap 95% CI for pooled r: [{ci_low:.3f}, {ci_high:.3f}]")

# Decision rule
print("\n=== DECISION RULE ===")
median_r = pp["r"].median()
n_high = int((pp["r"] > 0.8).sum())
verdict = "FAIL"
if r_pool > 0.7 and median_r > 0.7 and n_high >= 3:
    verdict = "PASS — Nat Comm reach"
elif r_pool > 0.5 and median_r > 0.5:
    verdict = "PARTIAL — Cell Rep Med / JCI Insight"
print(f"  pooled r = {r_pool:.3f}  (need >0.7 for PASS)")
print(f"  median patient r = {median_r:.3f}  (need >0.7 for PASS)")
print(f"  patients with r > 0.8 = {n_high}  (need >=3 for PASS)")
print(f"  VERDICT: {verdict}")

# Random null comparison — random gene sets matched to FVPTC size
sc.tl.score_genes(a, gene_list=list(np.random.default_rng(42).choice(a.var_names, size=len(FVPTC_GENES), replace=False)), score_name="score_random")
ptc.obs["score_random"] = a[ptc_mask].obs["score_random"].values
r_null, _ = stats.pearsonr(ptc.obs["DM_score"], ptc.obs["score_random"])
print(f"\nNull-set r (random {len(FVPTC_GENES)} genes vs DM_score) = {r_null:.3f}")

# Save outputs
pp.to_csv(OUT / "p2a_gse193581_per_patient_r.tsv", sep="\t", index=False)

# Figure: per-patient r forest
fig, ax = plt.subplots(figsize=(7, max(3, len(pp)*0.4)))
y_pos = np.arange(len(pp))
ax.barh(y_pos, pp["r"], xerr=None, height=0.6, color="steelblue", alpha=0.7)
ax.axvline(0.7, color="green", linestyle="--", lw=1, alpha=0.5, label="PASS threshold (0.7)")
ax.axvline(0.5, color="orange", linestyle="--", lw=1, alpha=0.5, label="PARTIAL threshold (0.5)")
ax.axvline(r_pool, color="red", lw=1.5, label=f"Pooled r = {r_pool:.3f}")
ax.set_yticks(y_pos)
ax.set_yticklabels([f"{s} (n={n})" for s, n in zip(pp["sample"], pp["n_cells"])])
ax.set_xlabel("Pearson r (DM_score vs FVPTC signature)")
ax.set_xlim(min(-0.1, pp['r'].min()-0.05), 1.0)
ax.set_title(f"GSE193581 PTC malignant cells — per-patient 8-gene↔FVPTC correlation\n"
             f"Pooled r = {r_pool:.3f} [95%CI {ci_low:.3f}-{ci_high:.3f}], median patient r = {median_r:.3f}")
ax.legend(loc="lower right", fontsize=8)
ax.grid(alpha=0.2)
fig.tight_layout()
fig.savefig(FIG / "p2a_per_patient_r_forest.png", dpi=200, bbox_inches="tight")
fig.savefig(FIG / "p2a_per_patient_r_forest.pdf", bbox_inches="tight")
print(f"\nSaved figure to {FIG}/p2a_per_patient_r_forest.png/.pdf")

# Pooled scatter
fig2, ax2 = plt.subplots(figsize=(7, 6))
ax2.scatter(ptc.obs["DM_score"], ptc.obs["score_fvptc"], c=ptc.obs["score_cptc"],
            cmap="RdBu_r", s=2, alpha=0.4, vmin=ptc.obs["score_cptc"].quantile(0.02), vmax=ptc.obs["score_cptc"].quantile(0.98))
m, b = np.polyfit(ptc.obs["DM_score"], ptc.obs["score_fvptc"], 1)
xs = np.linspace(ptc.obs["DM_score"].min(), ptc.obs["DM_score"].max(), 100)
ax2.plot(xs, m*xs + b, "k-", lw=1.5)
ax2.set_xlabel("8-gene DM score"); ax2.set_ylabel("FVPTC signature")
ax2.set_title(f"GSE193581 PTC malignant cells (n={len(ptc):,}) — pooled r = {r_pool:.3f}")
fig2.tight_layout()
fig2.savefig(FIG / "p2a_pooled_scatter.png", dpi=180, bbox_inches="tight")

# Save summary
result = {
    "cohort": "GSE193581 (Lu 2023 JCI)",
    "n_PTC_samples": int(len(pp)),
    "n_PTC_malignant_cells": int(len(ptc)),
    "pooled_r": float(r_pool),
    "pooled_p": float(p_pool),
    "bootstrap_95CI": [float(ci_low), float(ci_high)],
    "median_patient_r": float(median_r),
    "iqr_patient_r": [float(pp["r"].quantile(0.25)), float(pp["r"].quantile(0.75))],
    "n_patients_r_gt_0.7": int(n_high),
    "n_patients_r_gt_0.5": int((pp["r"] > 0.5).sum()),
    "null_random_genes_r": float(r_null),
    "phase1_GSE241184_r": 0.905,
    "verdict": verdict,
    "FVPTC_genes_used": FVPTC_GENES,
    "FVPTC_genes_missing_from_HVG": ["DIO1", "SLC5A5", "FOXE1"],
}
(OUT / "p2a_gse193581_summary.json").write_text(json.dumps(result, indent=2))
print(f"\n=== SAVED p2a_gse193581_summary.json ===")
print(json.dumps(result, indent=2))
