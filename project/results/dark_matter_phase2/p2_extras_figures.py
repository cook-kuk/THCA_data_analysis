"""Additional analytical figures for v2 dashboard — 7 new figures."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project")
OUT = ROOT / "results/dark_matter_phase2"
FIG = OUT / "web/figures"
FIG.mkdir(parents=True, exist_ok=True)

# Master data
df = pd.read_csv(ROOT / "results/dark_matter_phase1/tcga_dm_master_with_pfi.tsv", sep="\t")
master_v17 = pd.read_csv(ROOT / "results/v17/tables/sample_master_v17_full.tsv", sep="\t", low_memory=False)
master_v17["tcga_short"] = master_v17["sample_id"].str[:12]
master_tcga = master_v17[(master_v17["dataset"] == "TCGA-THCA") & (master_v17["normal_vs_tumor"] == "tumor")]
df = df.merge(master_tcga[["tcga_short", "histology_subtype"]].drop_duplicates("tcga_short"),
              on="tcga_short", how="left")

def driver_class(row):
    if row.get("has_braf_v600e"): return "BRAF V600E"
    if row.get("has_ras_mut"): return "RAS hotspot"
    da = str(row.get("driver_anchor_v17", ""))
    if da == "DICER1_EIF1AX_PPM1D": return "DICER1/EIF1AX"
    if row.get("tert_pos"): return "TERT-only"
    return "True driver-neg"
df["driver_class"] = df.apply(driver_class, axis=1)

# ============================================================================
# FIG E1: Cohort sizes bar chart
# ============================================================================
cohorts = pd.DataFrame([
    {"cohort": "TCGA-THCA bulk", "n": 482, "modality": "Bulk RNA-seq", "phase": "1+2"},
    {"cohort": "Yoo 2016 K2 bulk", "n": 180, "modality": "Bulk RNA-seq", "phase": "2"},
    {"cohort": "GSE241184 sc (Phase 1)", "n": 3, "modality": "scRNA-seq", "phase": "1"},
    {"cohort": "GSE193581 sc", "n": 23, "modality": "scRNA-seq", "phase": "2"},
    {"cohort": "GSE184362 sc", "n": 23, "modality": "scRNA-seq", "phase": "2"},
])
fig, ax = plt.subplots(figsize=(10, 4.5))
colors = ["#1864ab" if "Bulk" in m else "#2a9d8f" for m in cohorts["modality"]]
bars = ax.barh(cohorts["cohort"], cohorts["n"], color=colors, alpha=0.85)
for b, n in zip(bars, cohorts["n"]):
    ax.text(n + 5, b.get_y() + b.get_height()/2, f"n={n}", va="center", fontsize=10, fontweight="bold")
ax.set_xlabel("Sample size (n)")
ax.set_title("Figure E1 — Cohort sizes used across Phase 1 + Phase 2")
ax.set_xscale("log")
from matplotlib.patches import Patch
ax.legend(handles=[Patch(facecolor="#1864ab", label="Bulk RNA-seq"),
                   Patch(facecolor="#2a9d8f", label="scRNA-seq")], loc="lower right")
ax.grid(alpha=0.3, axis="x")
fig.tight_layout()
fig.savefig(FIG / "figE1_cohort_sizes.png", dpi=180, bbox_inches="tight")
plt.close(fig)

# ============================================================================
# FIG E2: TCGA driver class pie chart + n labels
# ============================================================================
class_counts = df["driver_class"].value_counts()
fig, ax = plt.subplots(figsize=(8, 6))
colors_pie = {"BRAF V600E": "#e63946", "RAS hotspot": "#f4a261",
              "True driver-neg": "#2a9d8f", "DICER1/EIF1AX": "#9d4edd",
              "TERT-only": "#1a1a2e"}
explode = [0.05 if c == "True driver-neg" else 0 for c in class_counts.index]
wedges, texts, autotexts = ax.pie(
    class_counts.values, labels=[f"{l}\n(n={n})" for l, n in zip(class_counts.index, class_counts.values)],
    colors=[colors_pie[c] for c in class_counts.index],
    autopct=lambda p: f"{p:.1f}%", startangle=90, explode=explode,
    textprops={"fontsize": 11}, wedgeprops={"edgecolor": "white", "linewidth": 2}
)
for t in autotexts:
    t.set_color("white"); t.set_fontweight("bold"); t.set_fontsize(11)
ax.set_title(f"Figure E2 — TCGA-THCA driver class distribution (n={len(df)})\n"
             f"True driver-neg = 28% (the 'Dark Matter' subgroup)")
fig.tight_layout()
fig.savefig(FIG / "figE2_driver_pie.png", dpi=180, bbox_inches="tight")
plt.close(fig)

# ============================================================================
# FIG E3: Cross-cohort sc r forest comparison
# ============================================================================
sc_summary = pd.DataFrame([
    {"cohort": "GSE241184 (Phase 1)", "n_patient": 1, "n_cells": 2427, "r": 0.905, "ci_lo": 0.895, "ci_hi": 0.913, "anchor": "PTC tumor thyrocyte"},
    {"cohort": "GSE193581 PTC (P2-A1)", "n_patient": 6, "n_cells": 8590, "r": 0.687, "ci_lo": 0.660, "ci_hi": 0.713, "anchor": "PTC malignant"},
    {"cohort": "GSE193581 PTC+ATC", "n_patient": 13, "n_cells": 14624, "r": 0.893, "ci_lo": 0.890, "ci_hi": 0.896, "anchor": "PTC + ATC malignant"},
    {"cohort": "GSE184362 (P2-A2)", "n_patient": 6, "n_cells": 21821, "r": 0.889, "ci_lo": 0.886, "ci_hi": 0.892, "anchor": "Adult PTC tumor thyrocytes"},
    {"cohort": "GSE184362 multi-site", "n_patient": 4, "n_cells": 13005, "r": 0.914, "ci_lo": 0.911, "ci_hi": 0.917, "anchor": "T+P+LN thyrocytes"},
])
fig, ax = plt.subplots(figsize=(11, 5))
y_pos = np.arange(len(sc_summary))
colors_f = ["#4361ee", "#f4a261", "#2a9d8f", "#e63946", "#9d4edd"]
ax.errorbar(sc_summary["r"], y_pos,
            xerr=[sc_summary["r"] - sc_summary["ci_lo"], sc_summary["ci_hi"] - sc_summary["r"]],
            fmt="o", color="#1a1a2e", markersize=10, capsize=5, lw=1.5)
for i, (idx, row) in enumerate(sc_summary.iterrows()):
    ax.scatter(row["r"], i, s=200 + np.log10(row["n_cells"])*30, c=colors_f[i], alpha=0.75, edgecolors="black", lw=1, zorder=10)
    ax.text(row["r"] + 0.02, i, f"  r={row['r']:.3f} (n={row['n_cells']:,} cells, {row['n_patient']} pts)",
            va="center", fontsize=9)
ax.axvline(0.7, color="green", linestyle="--", lw=1.5, alpha=0.7, label="PASS threshold (0.7)")
ax.axvline(0.5, color="orange", linestyle="--", lw=1.5, alpha=0.7, label="PARTIAL (0.5)")
ax.set_yticks(y_pos)
ax.set_yticklabels([f"{c}\n({row['anchor']})" for c, (_, row) in zip(sc_summary["cohort"], sc_summary.iterrows())], fontsize=9)
ax.set_xlim(0.4, 1.05)
ax.set_xlabel("Pearson r (8-gene panel ↔ FVPTC signature)")
ax.set_title("Figure E3 — Cross-cohort 8-gene↔FVPTC sc-level reproducibility\n"
             "All 4 sc cohorts above PARTIAL threshold; 4/5 above PASS threshold")
ax.legend(loc="lower right", fontsize=9)
ax.grid(alpha=0.3, axis="x")
fig.tight_layout()
fig.savefig(FIG / "figE3_cross_cohort_r_forest.png", dpi=180, bbox_inches="tight")
plt.close(fig)

# ============================================================================
# FIG E4: DM1 vs DM2 age + stage overlap
# ============================================================================
dm = df[df["dm_status"] & df["v17_dark_cluster"].isin(["DM1", "DM2"])].copy()
dm["age_num"] = pd.to_numeric(dm["age"], errors="coerce")
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
ax = axes[0]
dm1_age = dm.loc[dm["v17_dark_cluster"] == "DM1", "age_num"].dropna()
dm2_age = dm.loc[dm["v17_dark_cluster"] == "DM2", "age_num"].dropna()
ax.hist(dm1_age, bins=25, alpha=0.6, label=f"DM1 (n={len(dm1_age)}, mean={dm1_age.mean():.1f})", color="#4361ee", density=True)
ax.hist(dm2_age, bins=25, alpha=0.6, label=f"DM2 (n={len(dm2_age)}, mean={dm2_age.mean():.1f})", color="#f4a261", density=True)
ax.axvline(dm1_age.mean(), color="#4361ee", linestyle="--", lw=2)
ax.axvline(dm2_age.mean(), color="#f4a261", linestyle="--", lw=2)
t = stats.ttest_ind(dm1_age, dm2_age, equal_var=False)
d = (dm1_age.mean() - dm2_age.mean()) / np.sqrt((dm1_age.std()**2 + dm2_age.std()**2)/2)
ax.set_xlabel("Age at diagnosis (years)")
ax.set_ylabel("Density")
ax.set_title(f"Figure E4A — Age distribution DM1 vs DM2\nWelch p={t.pvalue:.2e}, Cohen d={d:.2f}, Δmean=12 years")
ax.legend(); ax.grid(alpha=0.3)

# Stage stacked bar
ax = axes[1]
stage_ct = pd.crosstab(dm["v17_dark_cluster"], dm["stage"], dropna=False).fillna(0)
stage_order = ["Stage I", "Stage II", "Stage III", "Stage IVA", "Stage IVC"]
stage_ct = stage_ct.reindex(columns=[s for s in stage_order if s in stage_ct.columns]).fillna(0)
stage_ct_pct = stage_ct.div(stage_ct.sum(axis=1), axis=0) * 100
stage_ct_pct.plot(kind="bar", stacked=True, ax=ax, color=plt.cm.RdYlGn_r(np.linspace(0.1, 0.9, stage_ct.shape[1])))
ax.set_ylabel("% of cluster")
ax.set_xlabel("")
ax.set_title("Figure E4B — Stage distribution within DM cluster")
ax.legend(title="AJCC Stage", bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
ax.grid(alpha=0.3, axis="y")
fig.tight_layout()
fig.savefig(FIG / "figE4_dm_age_stage.png", dpi=180, bbox_inches="tight")
plt.close(fig)

# ============================================================================
# FIG E5: Mutation gene frequency by DM cluster
# ============================================================================
import re
def parse_genes(s):
    if pd.isna(s) or s in ("", "[]"):
        return []
    return re.findall(r"'([A-Z0-9\-]+)'", str(s))
dm["genes"] = dm["mutation_genes"].apply(parse_genes)
genes_dm1 = pd.Series([g for genes in dm.loc[dm["v17_dark_cluster"]=="DM1", "genes"] for g in genes]).value_counts()
genes_dm2 = pd.Series([g for genes in dm.loc[dm["v17_dark_cluster"]=="DM2", "genes"] for g in genes]).value_counts()
all_genes = sorted(set(list(genes_dm1.index) + list(genes_dm2.index)))
gene_df = pd.DataFrame({"gene": all_genes,
                        "DM1_n": [genes_dm1.get(g, 0) for g in all_genes],
                        "DM2_n": [genes_dm2.get(g, 0) for g in all_genes]})
gene_df = gene_df[gene_df[["DM1_n", "DM2_n"]].sum(axis=1) > 0].sort_values(["DM2_n", "DM1_n"], ascending=False)
fig, ax = plt.subplots(figsize=(11, 4.5))
x = np.arange(len(gene_df))
ax.bar(x - 0.2, gene_df["DM1_n"], 0.4, label=f"DM1 (n={(dm['v17_dark_cluster']=='DM1').sum()} pts)", color="#4361ee", alpha=0.85)
ax.bar(x + 0.2, gene_df["DM2_n"], 0.4, label=f"DM2 (n={(dm['v17_dark_cluster']=='DM2').sum()} pts)", color="#f4a261", alpha=0.85)
ax.set_xticks(x)
ax.set_xticklabels(gene_df["gene"], rotation=45, ha="right", fontsize=9)
ax.set_ylabel("# patients with mutation")
ax.set_title("Figure E5 — Non-canonical mutation gene frequency within Dark Matter\nDICER1/TP53/PPM1D enriched in DM2; CHEK2 + RTK fusions (NTRK3/ALK/RET) in DM1")
ax.legend(); ax.grid(alpha=0.3, axis="y")
fig.tight_layout()
fig.savefig(FIG / "figE5_mutation_gene_freq.png", dpi=180, bbox_inches="tight")
plt.close(fig)

# ============================================================================
# FIG E6: Forest plot — HR by driver class
# ============================================================================
hr_data = pd.DataFrame([
    {"class": "BRAF V600E (ref)", "n": 291, "HR": 1.0, "lo": 1.0, "hi": 1.0, "p": np.nan, "events": 34},
    {"class": "RAS hotspot", "n": 54, "HR": 1.22, "lo": 0.55, "hi": 2.68, "p": 0.625, "events": 7},
    {"class": "DICER1/EIF1AX", "n": 7, "HR": 1.42, "lo": 0.20, "hi": 9.95, "p": 0.722, "events": 1},
    {"class": "TERT-only", "n": 5, "HR": 4.05, "lo": 0.97, "hi": 16.98, "p": 0.056, "events": 2},
    {"class": "True driver-neg", "n": 125, "HR": 0.49, "lo": 0.23, "hi": 1.06, "p": 0.071, "events": 6},
])
fig, ax = plt.subplots(figsize=(10, 4.5))
y = np.arange(len(hr_data))[::-1]
colors_h = ["#1a1a2e", "#f4a261", "#9d4edd", "#e63946", "#2a9d8f"]
for i, (_, row) in enumerate(hr_data.iterrows()):
    yy = y[i]
    if row["class"].endswith("(ref)"):
        ax.plot([1], [yy], "o", color=colors_h[i], markersize=12, label=row["class"])
        ax.text(1.1, yy, "(reference)", va="center", fontsize=9, color="gray")
    else:
        ax.plot([row["lo"], row["hi"]], [yy, yy], "-", color=colors_h[i], lw=2)
        ax.plot([row["HR"]], [yy], "s", color=colors_h[i], markersize=10)
        sig = "*" if row["p"] < 0.1 else ""
        ax.text(row["hi"] + 0.5, yy, f"  HR={row['HR']:.2f} [{row['lo']:.2f}-{row['hi']:.2f}], p={row['p']:.3f}{sig}", va="center", fontsize=9)
ax.axvline(1.0, color="gray", linestyle="--", alpha=0.5)
ax.set_yticks(y)
ax.set_yticklabels([f"{c} (n={n}, ev={e})" for c, n, e in zip(hr_data["class"], hr_data["n"], hr_data["events"])], fontsize=10)
ax.set_xscale("log")
ax.set_xlabel("Hazard Ratio (PFI, vs BRAF V600E)")
ax.set_xlim(0.1, 25)
ax.set_title("Figure E6 — PFI Cox HR forest plot by driver class (TCGA-THCA, Liu 2018 CDR)")
ax.grid(alpha=0.3, axis="x")
fig.tight_layout()
fig.savefig(FIG / "figE6_HR_forest.png", dpi=180, bbox_inches="tight")
plt.close(fig)

# ============================================================================
# FIG E7: Heatmap — cross-cohort × variable matrix
# ============================================================================
matrix = pd.DataFrame({
    "Expression": [3, 3, 3, 3, 3],
    "Clinical metadata": [3, 3, 1, 2, 2],
    "BRAF V600E call": [3, 3, 0, 0, 0],
    "RAS hotspot call": [3, 3, 0, 0, 0],
    "TERT promoter": [3, 0, 0, 0, 0],
    "DICER1/EIF1AX": [2, 3, 0, 0, 0],
    "Fusion calls": [2, 1, 0, 0, 0],
    "8-gene cluster": [3, 3, 2, 3, 3],
    "Survival (PFI)": [3, 0, 0, 0, 0],
}, index=["TCGA-THCA", "K2 (Yoo 2016)", "GSE241184 sc", "GSE193581 sc", "GSE184362 sc"])
fig, ax = plt.subplots(figsize=(12, 4.5))
im = ax.imshow(matrix.values, cmap="RdYlGn", aspect="auto", vmin=0, vmax=3)
ax.set_xticks(np.arange(matrix.shape[1]))
ax.set_xticklabels(matrix.columns, rotation=35, ha="right")
ax.set_yticks(np.arange(matrix.shape[0]))
ax.set_yticklabels(matrix.index)
for i in range(matrix.shape[0]):
    for j in range(matrix.shape[1]):
        v = matrix.values[i, j]
        sym = {0: "✗", 1: "⚠", 2: "◐", 3: "✓"}.get(v, "")
        ax.text(j, i, sym, ha="center", va="center", color="white" if v in (0, 3) else "black", fontsize=14)
cbar = plt.colorbar(im, ax=ax, ticks=[0, 1, 2, 3])
cbar.ax.set_yticklabels(["✗ Missing", "⚠ Sparse", "◐ Partial", "✓ Complete"])
ax.set_title("Figure E7 — Cohort × Variable availability matrix")
fig.tight_layout()
fig.savefig(FIG / "figE7_cohort_variable_matrix.png", dpi=180, bbox_inches="tight")
plt.close(fig)

print("Saved 7 new figures: figE1 .. figE7")
print(f"Figures dir now has: {len(list(FIG.glob('*.png')))} PNGs")
