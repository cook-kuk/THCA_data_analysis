#!/usr/bin/env python3
"""F3 — Driver mRNA neutrality multi-panel (Pillar 3)."""
from __future__ import annotations
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

PROJ = Path("/opt/thyroid-dash/project")
RES = PROJ / "results/figures/F3_driver_neutrality"
RES.mkdir(parents=True, exist_ok=True)

# Data
mut = pd.read_csv(PROJ / "results/p1_driver_mrna_audit/driver_mrna_mutation_audit.tsv", sep="\t")
auc = pd.read_csv(PROJ / "results/p1_driver_mrna_audit/driver_mrna_dm_auc.tsv", sep="\t")
top20 = pd.read_csv(PROJ / "results/p1_driver_mrna_audit/top20_by_d_full_tiera67.tsv", sep="\t")
expr = pd.read_csv("/data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_log2.tsv",
                   sep="\t", index_col=0)
muts_meta = pd.read_csv(PROJ / "results/tables/tcga_thca_mutation_groups.tsv", sep="\t").set_index("sample_id")

fig = plt.figure(figsize=(15, 11))
gs = GridSpec(2, 2, figure=fig, hspace=0.35, wspace=0.30, left=0.07, right=0.96, top=0.93, bottom=0.07)

# ---------- Panel A — BRAF mRNA × V600E (boxplot) ----------
axA = fig.add_subplot(gs[0, 0])
samps = [s for s in muts_meta.index if s in expr.columns]
m = muts_meta.loc[samps]
br_pos = expr.loc["BRAF", samps][m["has_braf_v600e"] == 1].dropna().values
br_wt = expr.loc["BRAF", samps][m["has_braf_v600e"] == 0].dropna().values

axA.boxplot([br_wt, br_pos], positions=[1, 2], widths=0.5, patch_artist=True,
             boxprops=dict(facecolor="#888888", alpha=0.7),
             medianprops=dict(color="black", lw=1.5))
axA.set_xticks([1, 2])
axA.set_xticklabels([f"WT (n={len(br_wt)})", f"V600E+ (n={len(br_pos)})"], fontsize=10)
axA.set_ylabel("BRAF log2 (expression)", fontsize=11)
axA.text(1.5, axA.get_ylim()[1]*0.97,
          "Cohen d = -0.044\nMW p = 0.567\n→ NULL", ha="center", fontsize=10,
          bbox=dict(boxstyle="round,pad=0.4", facecolor="#fff8c0", edgecolor="black"))
axA.set_title("A — BRAF mRNA × V600E mutation status (TCGA n=455)\nDriver transcript ≠ driver mutation",
               fontsize=11, loc="left", fontweight="bold")
axA.spines["top"].set_visible(False); axA.spines["right"].set_visible(False)

# ---------- Panel B — Driver single-feature AUC bar ----------
axB = fig.add_subplot(gs[0, 1])
genes_b = ["BRAF", "TERT", "KRAS", "NRAS", "HRAS"]
auc_vals = [float(auc[auc["gene"] == g]["auc"].values[0]) for g in genes_b]
colors_b = ["#cc4444" if v >= 0.7 else "#888888" for v in auc_vals]
axB.bar(np.arange(len(genes_b)), auc_vals, color=colors_b, edgecolor="black", lw=0.5)
axB.axhline(0.5, color="black", linestyle="--", lw=0.8)
axB.axhline(0.7, color="black", linestyle=":", lw=0.7, label="threshold 0.7")
axB.set_xticks(np.arange(len(genes_b))); axB.set_xticklabels(genes_b, fontsize=11)
axB.set_ylabel("Single-feature AUC for DM1 vs DM2", fontsize=11)
axB.set_ylim(0.45, 0.75)
for i, v in enumerate(auc_vals):
    axB.text(i, v + 0.005, f"{v:.3f}", ha="center", fontsize=9, fontweight="bold")
axB.set_title("B — Driver mRNA single-feature AUC (TCGA n=500)\nAll drivers < 0.7 — cannot define DM cluster",
               fontsize=11, loc="left", fontweight="bold")
axB.legend(loc="upper right", fontsize=8.5)
axB.spines["top"].set_visible(False); axB.spines["right"].set_visible(False)

# ---------- Panel C — TIERA67 Cohen d ranking ----------
axC = fig.add_subplot(gs[1, 0])
top24 = top20.head(24).copy()
top24["category"] = "other"
top24.loc[top24["in_8gene"] == True, "category"] = "8-gene"
top24.loc[top24["is_driver"] == True, "category"] = "driver"

color_map = {"8-gene": "#cc4444", "driver": "#3a4ea0", "other": "#cccccc"}
colors_c = [color_map[c] for c in top24["category"]]
axC.bar(np.arange(len(top24)), top24["abs_d"], color=colors_c, edgecolor="black", lw=0.4)
axC.set_xticks(np.arange(len(top24))); axC.set_xticklabels(top24["gene"], rotation=70, fontsize=8)
axC.set_ylabel("|Cohen d| (DM1 vs DM2)", fontsize=11)
axC.set_title("C — TIERA67 67-gene Cohen d ranking (top 24)\n8-gene panel at #4-50 (top tier); drivers at #52-67",
               fontsize=11, loc="left", fontweight="bold")
# Legend
from matplotlib.patches import Patch
axC.legend(handles=[Patch(color="#cc4444", label="8-gene panel"),
                      Patch(color="#3a4ea0", label="Driver_anchor"),
                      Patch(color="#cccccc", label="Other TIERA67")], loc="upper right", fontsize=9.5)
axC.spines["top"].set_visible(False); axC.spines["right"].set_visible(False)

# ---------- Panel D — Driver_anchor cluster ARI = 0 (Pillar 4 link) ----------
axD = fig.add_subplot(gs[1, 1])
ari_data = [
    ("Driver_anchor 12", -0.007, "#3a4ea0"),
    ("8-gene panel", 0.489, "#cc4444"),
    ("TDS_core 16", 0.467, "#cc7777"),
    ("Pan-genome top 200", 0.864, "#88aa88"),
    ("Pan-genome top 5000", 0.918, "#448844"),
    ("TIERA67 (full)", 0.903, "#226622"),
]
labels_d, vals_d, colors_d = zip(*ari_data)
axD.barh(np.arange(len(labels_d)), vals_d, color=colors_d, edgecolor="black", lw=0.5)
axD.set_yticks(np.arange(len(labels_d))); axD.set_yticklabels(labels_d, fontsize=10)
axD.set_xlabel("ARI vs original DM1/DM2 cluster", fontsize=11)
axD.axvline(0.0, color="black", lw=0.5)
for i, v in enumerate(vals_d):
    axD.text(v + 0.02 if v > 0 else v - 0.02, i, f"{v:.3f}", va="center", fontsize=9, fontweight="bold",
              ha="left" if v > 0 else "right")
axD.set_xlim(-0.1, 1.0)
axD.set_title("D — Driver-only cluster ARI ≈ 0 (random)\nDrivers cannot define DM cluster",
               fontsize=11, loc="left", fontweight="bold")
axD.spines["top"].set_visible(False); axD.spines["right"].set_visible(False)

fig.suptitle("Figure 3 — Driver mRNA neutrality (Pillar 3)", fontsize=14, fontweight="bold", y=0.99)
plt.savefig(RES / "F3_driver_neutrality.pdf")
plt.savefig(RES / "F3_driver_neutrality.png", dpi=200)
plt.close()
print(f"✓ F3 saved to {RES}")
