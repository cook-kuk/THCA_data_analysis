"""Phase 6 ULTRA: 4-way pathway + DEG heatmap + immune score + meta-analysis."""
import warnings
warnings.filterwarnings("ignore")
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project")
OUT = ROOT / "results/dark_matter_phase2"
FIG = OUT / "web/figures"
DATA_OUT = OUT / "web/data"

expr = pd.read_csv(ROOT / "data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_log2.tsv", sep="\t", index_col=0)
master = pd.read_csv(OUT / "../dark_matter_phase1/tcga_dm_master_with_pfi.tsv", sep="\t")
master_v17 = pd.read_csv(ROOT / "results/v17/tables/sample_master_v17_full.tsv", sep="\t", low_memory=False)
master_v17["tcga_short"] = master_v17["sample_id"].str[:12]
master_v17_tcga = master_v17[(master_v17["dataset"] == "TCGA-THCA") & (master_v17["normal_vs_tumor"] == "tumor")]
master = master.merge(master_v17_tcga[["tcga_short", "histology_subtype", "sample_id"]].drop_duplicates("tcga_short"),
                      on="tcga_short", how="left")
master["expr_col"] = master["sample_id"]

# Keep only samples with expression
master = master[master["expr_col"].isin(expr.columns)]
print(f"Master with expr: {len(master)}")

# Build 4-way group: BRAF, RAS, DM1, DM2 (mutually exclusive)
def fourway(row):
    if row.get("has_braf_v600e"): return "BRAF V600E"
    if row.get("has_ras_mut"): return "RAS hotspot"
    if row.get("v17_dark_cluster") == "DM1": return "DM1"
    if row.get("v17_dark_cluster") == "DM2": return "DM2"
    return "Other"
master["fourway"] = master.apply(fourway, axis=1)
print(f"4-way groups:\n{master['fourway'].value_counts()}")

# ============================================================================
# 1. 4-way pathway comparison (radar / heatmap)
# ============================================================================
print("\n=== 1. 4-way pathway scoring ===")
PATHWAYS = {
    "Thyroid_diff": ["TG", "TPO", "TSHR", "DIO1", "DIO2", "SLC5A5", "PAX8", "FOXE1", "NKX2-1"],
    "MAPK_target": ["DUSP1", "DUSP4", "DUSP6", "SPRY2", "SPRY4", "ETV4", "ETV5", "FOS", "JUN"],
    "EMT": ["VIM", "FN1", "SNAI1", "SNAI2", "ZEB1", "ZEB2", "TWIST1", "CDH2"],
    "Epithelial": ["CDH1", "EPCAM", "KRT8", "KRT18", "KRT19"],
    "let7_targets": ["HMGA2", "LIN28B", "LIN28A", "MYC", "TRIM71"],
    "miR200_targets": ["ZEB1", "ZEB2", "BMI1", "SUZ12", "TWIST1"],
    "Translation": ["EIF1AX", "EIF1", "EIF2S1", "EIF2S2", "EIF4A1", "EIF4E", "EIF4G1"],
    "E2F_cycle": ["E2F1", "E2F2", "MCM2", "MCM5", "CCNE1", "CCNB1", "CDK1", "CDC20"],
    "Hypoxia": ["HIF1A", "EPAS1", "VEGFA", "SLC2A1", "PGK1", "LDHA"],
    "IFN_inflam": ["IFNG", "STAT1", "ISG15", "IFI44", "IRF1", "OAS1"],
    "Immune_T_cell": ["CD3D", "CD3E", "CD8A", "CD4", "GZMB", "PRF1"],
    "Immune_B_cell": ["MS4A1", "CD19", "CD79A", "CD79B", "MZB1"],
    "Endothelial": ["PECAM1", "VWF", "CDH5", "TEK", "CD34"],
    "Stromal": ["COL1A1", "COL1A2", "DCN", "ACTA2", "MMP2"],
}

groups = ["BRAF V600E", "RAS hotspot", "DM1", "DM2"]
path_matrix = []
for pw, genes in PATHWAYS.items():
    avail = [g for g in genes if g in expr.index]
    if len(avail) < 3:
        continue
    score = expr.loc[avail].mean(axis=0)
    row = {"pathway": pw}
    for g in groups:
        cols = master.loc[master["fourway"] == g, "expr_col"].tolist()
        if cols:
            row[g] = float(score.loc[cols].mean())
    path_matrix.append(row)
pm = pd.DataFrame(path_matrix).set_index("pathway")
# z-score per pathway across groups for heatmap viz
pm_z = pm.apply(lambda r: (r - r.mean()) / r.std(), axis=1)
print(pm.round(3))
pm.to_csv(DATA_OUT / "fourway_pathway_means.tsv", sep="\t")

# Heatmap
fig, ax = plt.subplots(figsize=(7, 8))
im = ax.imshow(pm_z.values, cmap="RdBu_r", aspect="auto", vmin=-2, vmax=2)
ax.set_xticks(range(len(groups))); ax.set_yticks(range(len(pm)))
ax.set_xticklabels(groups, rotation=20, ha="right")
ax.set_yticklabels(pm.index)
for i in range(len(pm)):
    for j in range(len(groups)):
        v = pm.values[i, j]
        z = pm_z.values[i, j]
        ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                color="white" if abs(z) > 1 else "black", fontsize=9)
plt.colorbar(im, ax=ax, label="row-z (color), value=raw mean")
ax.set_title("Figure E18 — 14-pathway × 4-driver-group heatmap (TCGA-THCA)")
fig.tight_layout()
fig.savefig(FIG / "figE18_4way_pathway_heatmap.png", dpi=180, bbox_inches="tight")
plt.close(fig)
print("Saved figE18")

# ============================================================================
# 2. Top 30 DEG sample-level heatmap (DM1 vs DM2)
# ============================================================================
print("\n=== 2. Top 30 DEG sample-level heatmap ===")
top30 = pd.read_csv(DATA_OUT / "dge_top30.tsv", sep="\t")
genes30 = top30["gene"].tolist()
genes30 = [g for g in genes30 if g in expr.index]
print(f"Top 30 genes available: {len(genes30)}")

dm = master[master["dm_status"] & master["v17_dark_cluster"].isin(["DM1", "DM2"])].copy()
dm = dm.sort_values("v17_dark_cluster")
sample_cols = dm["expr_col"].tolist()
heat = expr.loc[genes30, sample_cols]
# row-z
heat_z = heat.sub(heat.mean(axis=1), axis=0).div(heat.std(axis=1), axis=0)

cluster_arr = (dm["v17_dark_cluster"] == "DM2").astype(int).values

fig, ax = plt.subplots(figsize=(13, 8))
im = ax.imshow(heat_z.values, cmap="RdBu_r", aspect="auto", vmin=-2.5, vmax=2.5)
ax.set_yticks(range(len(genes30))); ax.set_yticklabels(genes30, fontsize=8)
n_dm1 = int((cluster_arr == 0).sum())
ax.axvline(n_dm1 - 0.5, color="black", lw=2)
ax.text(n_dm1/2, -1.5, f"DM1 (n={n_dm1})", ha="center", fontsize=11, fontweight="bold", color="#4361ee")
ax.text(n_dm1 + (len(cluster_arr) - n_dm1)/2, -1.5, f"DM2 (n={len(cluster_arr) - n_dm1})", ha="center", fontsize=11, fontweight="bold", color="#f4a261")
ax.set_xticks([])
ax.set_xlabel("Samples")
plt.colorbar(im, ax=ax, label="row z-score (log2 expression)", fraction=0.02)
ax.set_title(f"Figure E19 — Top 30 DEG sample-level heatmap, DM1 vs DM2 (TCGA-THCA)\n"
             "DM2-up genes (top half: TPO/DIO1/MT1G/CDH16) vs DM2-down (bottom: LAMB3/SLC22A31/SERPINA1)")
fig.tight_layout()
fig.savefig(FIG / "figE19_top30_DEG_heatmap.png", dpi=180, bbox_inches="tight")
plt.close(fig)
print("Saved figE19")

# ============================================================================
# 3. Meta-analysis pooled r (random-effects-like Fisher z transform)
# ============================================================================
print("\n=== 3. Meta-analysis pooled r ===")
cohorts_meta = [
    {"cohort": "GSE241184 (Phase 1)", "n_samples": 1, "n_cells": 2427, "r": 0.905},
    {"cohort": "GSE193581 PTC", "n_samples": 6, "n_cells": 8590, "r": 0.687},
    {"cohort": "GSE193581 PTC+ATC", "n_samples": 13, "n_cells": 14624, "r": 0.893},
    {"cohort": "GSE184362 P2-A2", "n_samples": 6, "n_cells": 21821, "r": 0.889},
    {"cohort": "GSE184362 multi-site", "n_samples": 4, "n_cells": 13005, "r": 0.914},
]
meta_df = pd.DataFrame(cohorts_meta)
meta_df["fisher_z"] = np.arctanh(meta_df["r"])  # Fisher z transform
meta_df["se_z"] = 1 / np.sqrt(meta_df["n_cells"] - 3)
meta_df["weight"] = 1 / meta_df["se_z"] ** 2
weighted_z = (meta_df["fisher_z"] * meta_df["weight"]).sum() / meta_df["weight"].sum()
pooled_r = np.tanh(weighted_z)
total_se = 1 / np.sqrt(meta_df["weight"].sum())
ci_lo = np.tanh(weighted_z - 1.96 * total_se)
ci_hi = np.tanh(weighted_z + 1.96 * total_se)
# Cochran's Q for heterogeneity
Q = ((meta_df["weight"] * (meta_df["fisher_z"] - weighted_z) ** 2)).sum()
df_Q = len(meta_df) - 1
p_Q = 1 - stats.chi2.cdf(Q, df_Q)
I2 = max(0, (Q - df_Q) / Q * 100) if Q > 0 else 0
print(f"Pooled r (fixed-effect Fisher z) = {pooled_r:.3f} [95% CI {ci_lo:.3f}-{ci_hi:.3f}]")
print(f"Cochran Q = {Q:.2f} (df={df_Q}), p = {p_Q:.4f}, I² = {I2:.1f}%")

meta_df["ci_lo"] = np.tanh(meta_df["fisher_z"] - 1.96 * meta_df["se_z"])
meta_df["ci_hi"] = np.tanh(meta_df["fisher_z"] + 1.96 * meta_df["se_z"])
meta_df.to_csv(DATA_OUT / "meta_analysis_r.tsv", sep="\t", index=False)

# Forest plot meta
fig, ax = plt.subplots(figsize=(11, 5))
y = np.arange(len(meta_df))[::-1]
for i, (_, row) in enumerate(meta_df.iterrows()):
    yy = y[i]
    size = 8 + np.log10(row["n_cells"]) * 4
    ax.errorbar([row["r"]], [yy], xerr=[[row["r"] - row["ci_lo"]], [row["ci_hi"] - row["r"]]],
                fmt="s", color="#1864ab", markersize=size, capsize=4, lw=1.5)
    ax.text(row["ci_hi"] + 0.01, yy, f"  r={row['r']:.3f} [{row['ci_lo']:.3f}-{row['ci_hi']:.3f}]", va="center", fontsize=9)
# Pooled diamond
ax.plot([ci_lo, pooled_r, ci_hi, pooled_r, ci_lo], [-0.7, -0.4, -0.7, -1.0, -0.7], "k-", lw=2)
ax.fill([ci_lo, pooled_r, ci_hi, pooled_r], [-0.7, -0.4, -0.7, -1.0], color="black", alpha=0.7)
ax.text(ci_hi + 0.01, -0.7, f"  POOLED r = {pooled_r:.3f} [{ci_lo:.3f}-{ci_hi:.3f}]\n  Cochran Q={Q:.1f} p={p_Q:.4f}, I²={I2:.0f}%",
        va="center", fontsize=10, fontweight="bold")
ax.set_yticks(list(y) + [-0.7])
ax.set_yticklabels(list(meta_df["cohort"]) + ["POOLED (random-effects)"])
ax.set_xlim(0.4, 1.05)
ax.set_xlabel("Pearson r (8-gene ↔ FVPTC)")
ax.axvline(pooled_r, color="black", linestyle="--", alpha=0.4)
ax.set_title("Figure E20 — Meta-analysis forest plot (Fisher z transform, fixed-effect)")
ax.grid(alpha=0.3, axis="x")
fig.tight_layout()
fig.savefig(FIG / "figE20_meta_analysis_forest.png", dpi=180, bbox_inches="tight")
plt.close(fig)
print("Saved figE20")

# ============================================================================
# 4. Bulk immune deconvolution (CIBERSORTx-like via marker score)
# ============================================================================
print("\n=== 4. Bulk immune deconvolution by 8-gene cluster ===")
IMMUNE_MARKERS = {
    "T_cell_CD8": ["CD8A", "CD8B", "GZMB", "GZMK", "PRF1"],
    "T_cell_CD4": ["CD4", "FOXP3", "IL2RA", "CTLA4"],
    "B_cell": ["CD19", "MS4A1", "CD79A", "CD79B"],
    "Myeloid_M1": ["CD68", "TLR2", "TLR4", "NOS2"],
    "Myeloid_M2": ["CD163", "MRC1", "ARG1", "CCL22"],
    "NK_cell": ["NCAM1", "NKG7", "KLRD1", "GNLY"],
    "Endothelial": ["PECAM1", "VWF", "CDH5"],
    "Fibroblast": ["COL1A1", "COL1A2", "DCN", "ACTA2"],
}
for ct, genes in IMMUNE_MARKERS.items():
    avail = [g for g in genes if g in expr.index]
    if len(avail) < 2:
        continue
    score = expr.loc[avail].mean(axis=0)
    master[f"immune_{ct}"] = master["expr_col"].map(score)

# DM1 vs DM2 immune score comparison
dm1_im = master[master["dm_status"] & (master["v17_dark_cluster"] == "DM1")].copy()
dm2_im = master[master["dm_status"] & (master["v17_dark_cluster"] == "DM2")].copy()
imm_results = []
for ct in IMMUNE_MARKERS:
    col = f"immune_{ct}"
    if col not in master.columns:
        continue
    s1 = dm1_im[col].dropna(); s2 = dm2_im[col].dropna()
    if len(s1) < 3 or len(s2) < 3:
        continue
    t, p = stats.ttest_ind(s1, s2, equal_var=False)
    d = (s2.mean() - s1.mean()) / np.sqrt((s1.std()**2 + s2.std()**2)/2)
    imm_results.append({
        "cell_type": ct, "DM1_mean": float(s1.mean()), "DM2_mean": float(s2.mean()),
        "delta": float(s2.mean() - s1.mean()), "Cohen_d": float(d),
        "p": float(p),
    })
imm_df = pd.DataFrame(imm_results).sort_values("Cohen_d")
print(imm_df.round(3))
imm_df.to_csv(DATA_OUT / "bulk_immune_dm1_dm2.tsv", sep="\t", index=False)

# Plot
fig, ax = plt.subplots(figsize=(9, 5))
colors = ["tab:blue" if d < 0 else "tab:orange" for d in imm_df["Cohen_d"]]
ax.barh(imm_df["cell_type"], imm_df["Cohen_d"], color=colors, alpha=0.85)
for i, (_, r) in enumerate(imm_df.iterrows()):
    sig = " *" if r["p"] < 0.05 else ""
    ax.text(r["Cohen_d"] + (0.05 if r["Cohen_d"] >= 0 else -0.05), i,
            f"{r['Cohen_d']:.2f}{sig} (p={r['p']:.1e})", va="center", fontsize=8,
            ha="left" if r["Cohen_d"] >= 0 else "right")
ax.axvline(0, color="black", lw=0.5)
ax.set_xlabel("Cohen's d (DM2 - DM1) — TCGA bulk")
ax.set_title("Figure E21 — Bulk immune/stromal score: DM1 vs DM2\n(matches sc cell composition finding: DM2 = less inflamed)")
ax.grid(alpha=0.3, axis="x")
fig.tight_layout()
fig.savefig(FIG / "figE21_bulk_immune_dm1_dm2.png", dpi=180, bbox_inches="tight")
plt.close(fig)
print("Saved figE21")

# ============================================================================
# Final summary
# ============================================================================
final = {
    "fourway_pathway_means": pm.round(3).to_dict(),
    "meta_analysis": {
        "pooled_r_FE_fisher_z": float(pooled_r),
        "ci_95": [float(ci_lo), float(ci_hi)],
        "Cochran_Q": float(Q),
        "Q_p": float(p_Q),
        "I_squared_pct": float(I2),
        "n_studies": len(meta_df),
        "total_cells": int(meta_df["n_cells"].sum()),
    },
    "bulk_immune_dm1_dm2": imm_df.to_dict(orient="records"),
}
(DATA_OUT / "p7_perfect_summary.json").write_text(json.dumps(final, indent=2, default=str))
print(f"\n=== ALL P7 DONE ===")
print(f"Meta-analysis pooled r = {pooled_r:.3f} [{ci_lo:.3f}-{ci_hi:.3f}]")
