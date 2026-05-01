#!/usr/bin/env python3
"""F4 — Pan-genome cluster robustness (Pillar 4)."""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

PROJ = Path("/opt/thyroid-dash/project")
RES = PROJ / "results/figures/F4_pangenome_robustness"
RES.mkdir(parents=True, exist_ok=True)

ari = pd.read_csv(PROJ / "results/p4_pangenome_vs_tiera67/ari_comparison.tsv", sep="\t")
ladder = pd.read_csv(PROJ / "results/p4_pangenome_vs_tiera67/topN_coverage_ladder.tsv", sep="\t")

fig = plt.figure(figsize=(14, 6))
gs = GridSpec(1, 2, figure=fig, wspace=0.3, left=0.07, right=0.96, top=0.86, bottom=0.18)

# ---------- A — ARI bar ----------
axA = fig.add_subplot(gs[0, 0])
ari_sorted = ari.sort_values("ARI")
colors_a = []
for label in ari_sorted["label"]:
    if "Driver" in label:
        colors_a.append("#3a4ea0")
    elif "8-gene" in label or "TDS" in label:
        colors_a.append("#cc4444")
    elif "TIERA67" in label:
        colors_a.append("#226622")
    else:
        colors_a.append("#88aa88")

bars = axA.barh(np.arange(len(ari_sorted)), ari_sorted["ARI"], color=colors_a, edgecolor="black", lw=0.5)
axA.set_yticks(np.arange(len(ari_sorted))); axA.set_yticklabels(ari_sorted["label"], fontsize=10)
for i, v in enumerate(ari_sorted["ARI"]):
    axA.text(v + 0.02 if v > 0 else v - 0.02, i, f"{v:+.3f}", va="center", fontsize=9, fontweight="bold",
              ha="left" if v > 0 else "right")
axA.axvline(0.0, color="black", lw=0.5)
axA.set_xlabel("ARI vs DM1/DM2 cluster", fontsize=11)
axA.set_xlim(-0.15, 1.05)
axA.set_title("A — Cluster panel ARI comparison\nTIERA67 ≈ pan-genome top 5000; 8-gene alone modest; Driver-only ≈ 0",
               fontsize=11, loc="left", fontweight="bold")
axA.spines["top"].set_visible(False); axA.spines["right"].set_visible(False)

# ---------- B — Top-N coverage ladder ----------
axB = fig.add_subplot(gs[0, 1])
N_vals = ladder["top_N"].values
TIERA_pct = ladder["pct_of_TIERA"].values
g8_pct = ladder["pct_of_g8"].values

axB.plot(N_vals, TIERA_pct, "o-", color="#226622", lw=2, markersize=8, label="TIERA67 (67 genes)")
axB.plot(N_vals, g8_pct, "s-", color="#cc4444", lw=2, markersize=8, label="8-gene panel")
axB.set_xscale("log")
axB.set_xlabel("Pan-genome top N (univariate Cohen d)", fontsize=11)
axB.set_ylabel("% covered", fontsize=11)
axB.legend(loc="upper left", fontsize=10)
axB.set_title("B — Pan-genome top-N coverage of TIERA67 + 8-gene\nHypergeometric p (TIERA67 in top 100) = 3e-4 ★",
               fontsize=11, loc="left", fontweight="bold")
axB.grid(alpha=0.3, which="both")
axB.spines["top"].set_visible(False); axB.spines["right"].set_visible(False)

fig.suptitle("Figure 4 — Pan-genome cluster robustness (Pillar 4)", fontsize=14, fontweight="bold", y=0.97)
plt.savefig(RES / "F4_pangenome_robustness.pdf")
plt.savefig(RES / "F4_pangenome_robustness.png", dpi=200)
plt.close()
print(f"✓ F4 saved to {RES}")
