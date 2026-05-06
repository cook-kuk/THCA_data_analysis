"""Step 8 figures for the external expression validation sweep."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parent.parent

SCORES = pd.read_csv(ROOT / "external_sample_scores.tsv.gz", sep="\t")
DIRECTION = pd.read_csv(ROOT / "external_direction_consistency.tsv", sep="\t")
COVERAGE = pd.read_csv(ROOT / "external_gene_coverage.tsv", sep="\t")
SPEAR = pd.read_csv(ROOT / "external_spearman.tsv", sep="\t")

GROUP_ORDER = ["normal", "PTC", "FVPTC", "FTC", "PDTC", "ATC"]
GROUP_COLORS = {
    "normal": "#7f8fa6",
    "PTC":    "#79b8e0",
    "FVPTC":  "#5f9ec7",
    "FTC":    "#3f86b0",
    "PDTC":   "#e6914d",
    "ATC":    "#c0392b",
}

DATASETS = ["GSE33630", "GSE29265", "GSE65144", "GSE53157"]

# ---------------------------------------------------------------------------
# Fig 1: RAI_8 + DM1_like + THYROID_NONOVERLAP boxplots, 4 datasets x 3 panels
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(3, 4, figsize=(15, 9), sharex=False)
for col, ds in enumerate(DATASETS):
    sub = SCORES[(SCORES["dataset"] == ds) & ~SCORES["histology_clean"].isin(["UNKNOWN", "DROP_pool"])]
    groups = [g for g in GROUP_ORDER if g in sub["histology_clean"].unique()]
    for row, score in enumerate(["RAI_8_score", "DM1_like_score", "THYROID_NONOVERLAP_score"]):
        ax = axes[row, col]
        data = [sub.loc[sub["histology_clean"] == g, score].dropna().values for g in groups]
        bp = ax.boxplot(data, positions=range(len(groups)), widths=0.6,
                        patch_artist=True, showfliers=False)
        for patch, g in zip(bp["boxes"], groups):
            patch.set_facecolor(GROUP_COLORS.get(g, "#bbbbbb"))
            patch.set_alpha(0.8)
        # overlay points
        for i, vals in enumerate(data):
            if len(vals):
                ax.scatter(np.random.normal(i, 0.06, len(vals)), vals, s=10,
                           color="black", alpha=0.45, zorder=3)
        ax.set_xticks(range(len(groups)))
        ax.set_xticklabels(groups, rotation=30, fontsize=8)
        ax.axhline(0, color="grey", lw=0.5, ls="--")
        if col == 0:
            ax.set_ylabel(score, fontsize=9)
        if row == 0:
            n_total = sum(len(v) for v in data)
            ax.set_title(f"{ds}  (n={n_total})", fontsize=10)
        ax.tick_params(axis="y", labelsize=8)
fig.suptitle("External replication — RAI_8 / DM1_like / THYROID_NONOVERLAP per histology, per dataset", fontsize=12)
plt.tight_layout()
plt.savefig(ROOT / "external_rai_lineage_boxplots.png", dpi=160, bbox_inches="tight")
plt.close()

# ---------------------------------------------------------------------------
# Fig 2: DM1_like vs THYROID_NONOVERLAP scatter, 2x2 grid
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(11, 9))
for ax, ds in zip(axes.flat, DATASETS):
    sub = SCORES[(SCORES["dataset"] == ds) & ~SCORES["histology_clean"].isin(["UNKNOWN", "DROP_pool"])]
    rho_row = SPEAR[(SPEAR["dataset"] == ds) &
                    (SPEAR["x"] == "DM1_like_score") &
                    (SPEAR["y"] == "THYROID_NONOVERLAP_score")]
    rho = rho_row["rho"].iloc[0] if len(rho_row) else float("nan")
    p = rho_row["p"].iloc[0] if len(rho_row) else float("nan")
    for g in GROUP_ORDER:
        m = sub["histology_clean"] == g
        if m.sum():
            ax.scatter(sub.loc[m, "DM1_like_score"], sub.loc[m, "THYROID_NONOVERLAP_score"],
                       color=GROUP_COLORS[g], label=f"{g} (n={m.sum()})", s=32,
                       edgecolor="black", linewidth=0.4, alpha=0.85)
    ax.axhline(0, color="grey", lw=0.5, ls="--")
    ax.axvline(0, color="grey", lw=0.5, ls="--")
    ax.set_xlabel("DM1_like_score (= -RAI_8)")
    ax.set_ylabel("THYROID_NONOVERLAP_score")
    ax.set_title(f"{ds}    Spearman ρ={rho:.2f}  p={p:.1e}", fontsize=10)
    ax.legend(fontsize=8, loc="best", frameon=True)
plt.suptitle("DM1-like axis vs orthogonal thyroid lineage panel — per dataset", fontsize=12)
plt.tight_layout()
plt.savefig(ROOT / "external_dm1_nonoverlap_scatter_grid.png", dpi=160, bbox_inches="tight")
plt.close()

# ---------------------------------------------------------------------------
# Fig 3: direction-consistency forest (advanced_vs_DTC + advanced_vs_normal)
# ---------------------------------------------------------------------------
focus_scores = ["RAI_8_score", "DM1_like_score", "THYROID_NONOVERLAP_score",
                "TDS_like_score", "STAT3_AP1_DNMT_score", "TACSTD2_z"]
forest = DIRECTION[DIRECTION["score"].isin(focus_scores)].copy()
fig, ax = plt.subplots(figsize=(10, 9))
y = 0
labels = []
for ds in DATASETS:
    for contrast in ["advanced_vs_DTC", "advanced_vs_normal"]:
        sub = forest[(forest["dataset"] == ds) & (forest["contrast"] == contrast)]
        if not len(sub):
            continue
        for score in focus_scores:
            row = sub[sub["score"] == score]
            if not len(row):
                continue
            d = float(row["cohen_d"].iloc[0])
            p = float(row["p"].iloc[0])
            n1 = int(row["n1"].iloc[0]); n2 = int(row["n2"].iloc[0])
            color = "#c0392b" if d > 0 else "#2c5e9c"
            alpha = 1.0 if p < 0.05 else 0.35
            ax.scatter(d, y, s=80, color=color, alpha=alpha, edgecolor="black", linewidth=0.6)
            ax.plot([0, d], [y, y], color=color, alpha=alpha, lw=2)
            labels.append(f"{ds} | {contrast} | {score}  (n={n1}/{n2})")
            y += 1
        y += 0.5
ax.axvline(0, color="black", lw=1)
ax.set_yticks(range(len(labels)))
ax.set_yticklabels(labels, fontsize=8)
ax.invert_yaxis()
ax.set_xlabel("Cohen's d  (positive = higher in advanced)")
ax.set_title("Direction consistency: advanced (ATC/PDTC) vs DTC and vs normal — opaque points = p<0.05",
             fontsize=11)
plt.tight_layout()
plt.savefig(ROOT / "external_direction_consistency_forest.png", dpi=160, bbox_inches="tight")
plt.close()

# ---------------------------------------------------------------------------
# Fig 4: gene-coverage QC heatmap (rows = genes; cols = datasets)
# ---------------------------------------------------------------------------
panels_for_qc = ["RAI_8", "THYROID_NONOVERLAP", "TDS_TF", "MECHANISM", "OPTIONAL_TARGETS"]
cov = COVERAGE[COVERAGE["panel"].isin(panels_for_qc)].copy()
pivot = cov.pivot_table(index=["panel", "gene"], columns="dataset",
                        values="available", aggfunc="max").fillna(0)
fig, ax = plt.subplots(figsize=(8, 0.28 * len(pivot) + 1))
im = ax.imshow(pivot.values, cmap="Greens", vmin=0, vmax=1, aspect="auto")
ax.set_xticks(range(pivot.shape[1]))
ax.set_xticklabels(pivot.columns, rotation=45, ha="right")
ax.set_yticks(range(pivot.shape[0]))
ax.set_yticklabels([f"{p} | {g}" for p, g in pivot.index], fontsize=7)
for i in range(pivot.shape[0]):
    for j in range(pivot.shape[1]):
        ax.text(j, i, "✓" if pivot.values[i, j] else "✗",
                ha="center", va="center", fontsize=7,
                color="white" if pivot.values[i, j] else "black")
ax.set_title("Gene-panel availability per dataset (✓ = present after probe→gene collapse)", fontsize=10)
plt.tight_layout()
plt.savefig(ROOT / "external_dataset_qc_heatmap.png", dpi=160, bbox_inches="tight")
plt.close()

print("[figures] written:")
for f in [
    "external_rai_lineage_boxplots.png",
    "external_dm1_nonoverlap_scatter_grid.png",
    "external_direction_consistency_forest.png",
    "external_dataset_qc_heatmap.png",
]:
    p = ROOT / f
    print(f"  {p.relative_to(ROOT.parent)}  ({p.stat().st_size:,} B)")
