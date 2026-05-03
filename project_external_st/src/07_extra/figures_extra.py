#!/usr/bin/env python3
"""Figures for extra analyses A1, A2, A3."""
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUT = Path("project_external_st/results/extra")

# ---- A1: drop-one-out Pearson r ----
a1 = pd.read_csv(OUT / "a1_drop_one_out_summary.tsv", sep="\t").sort_values("pearson_r")
fig, ax = plt.subplots(figsize=(9, 5.5))
colors = ["#962E2E" if g != "(none, full RAI_8)" else "#3C6B4F" for g in a1["gene_dropped"]]
bars = ax.barh(a1["gene_dropped"], a1["pearson_r"], color=colors, edgecolor="black")
for i, (g, r, p) in enumerate(zip(a1["gene_dropped"], a1["pearson_r"], a1["pearson_p"])):
    ax.text(r - 0.005, i, f"r={r:.3f}\np={p:.1e}", va="center", ha="right",
            fontsize=10, color="white", fontweight="bold")
ax.set_xlabel("Pearson r between sample-mean DM1_drop and THYROID_NONOVERLAP (n=12 external slides)")
ax.set_xlim(-1.0, -0.97)
ax.axvline(a1[a1["gene_dropped"] == "(none, full RAI_8)"]["pearson_r"].iloc[0],
           color="black", lw=1.5, ls="--", label="Full RAI_8 baseline")
ax.legend(loc="lower right")
ax.set_title("A1 — Drop-one-out: dropping ANY single RAI_8 gene leaves r > −0.98\n"
             "Cross-validation result is panel-design-robust",
             fontsize=11)
fig.tight_layout()
fig.savefig(OUT / "a1_drop_one_out.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("→ a1_drop_one_out.png")

# ---- A2: Moran's I by condition ----
a2 = pd.read_csv(OUT / "a2_morans_i.tsv", sep="\t")
order = ["PT","PTC","LPTC","ATC","CONTROL","HT","GD","PTC_HT"]
a2 = a2[a2["condition"].isin(order)].copy()
fig, ax = plt.subplots(figsize=(11, 5.5))
import seaborn as sns
ds_palette = {"GSE250521":"#962E2E", "GSE230424":"#B8893C", "GSE248205":"#2C5C8A"}
sns.stripplot(data=a2, x="condition", y="morans_I", order=order, ax=ax,
              hue="dataset", size=12, jitter=0.18, palette=ds_palette,
              edgecolor="black", linewidth=0.5)
means = a2.groupby("condition")["morans_I"].mean().reindex(order)
for i, m in enumerate(means.values):
    ax.scatter(i, m, marker="_", s=400, color="black", lw=3, zorder=10)
ax.axhline(0, color="grey", lw=0.5, ls=":")
ax.axhline(a2["morans_I"].mean(), color="#962E2E", lw=1, ls="--",
           label=f"Global mean = {a2['morans_I'].mean():.2f}")
ax.set_xlabel(""); ax.set_ylabel("Moran's I (DM1_like, k=6 NN)")
ax.set_title("A2 — DM1_like spatial autocorrelation per slide (Moran's I)\n"
             "ATC drops to 0.13: thyroid spatial structure dissolves in advanced cancer",
             fontsize=11)
ax.legend(loc="upper right", frameon=False)
fig.tight_layout()
fig.savefig(OUT / "a2_morans_i.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("→ a2_morans_i.png")

# ---- A3: DM1 quartile × cell type stacked bar ----
a3 = pd.read_csv(OUT / "a3_DM1_quartile_x_celltype.tsv", sep="\t", index_col=0)
a3 = a3[["Epithelial","Stromal","Immune_T","Immune_B","Macrophage","Endothelial"]]
fig, ax = plt.subplots(figsize=(10, 5.5))
ct_palette = {"Epithelial":"#3C6B4F","Stromal":"#B8893C","Immune_T":"#962E2E",
              "Immune_B":"#7B1F2A","Macrophage":"#A57215","Endothelial":"#2C5C8A"}
bottoms = np.zeros(len(a3))
for ct in a3.columns:
    ax.bar(a3.index, a3[ct], bottom=bottoms, label=ct, color=ct_palette[ct],
           edgecolor="black", linewidth=0.5)
    for i, v in enumerate(a3[ct]):
        if v > 0.05:
            ax.text(i, bottoms[i] + v / 2, f"{v*100:.0f}%", ha="center", va="center",
                    color="white", fontsize=10, fontweight="bold")
    bottoms = bottoms + a3[ct].values
ax.set_xlabel("DM1_like quartile (Q1=lowest DM1, thyroid-preserved → Q4=highest DM1, dedifferentiated)")
ax.set_ylabel("Top-1 cell-type fraction (within-sample marker assignment)")
ax.set_ylim(0, 1.05)
ax.set_title("A3 — Within-slide cell-type composition by DM1 quartile (12 external slides combined)\n"
             "DM1-high spots: epithelial 33% → 10%, immune+stromal+endothelial fill in",
             fontsize=11)
ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False)
fig.tight_layout()
fig.savefig(OUT / "a3_dm1_quartile_celltype.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("→ a3_dm1_quartile_celltype.png")
