#!/usr/bin/env python3
"""Render forest plots for R1 meta-survival results."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p2_braf_nature_sprint_2026_05_09/r1_meta_survival")
forest = pd.read_csv(ROOT / "r1_forest_data.tsv", sep="\t")
print(f"loaded {len(forest)} rows; panels: {forest['forest_panel'].unique()}")


def panel_plot(ax, df, title, xlim=None):
    df = df.dropna(subset=["HR"]).copy()
    df = df.reset_index(drop=True)
    if len(df) == 0:
        ax.text(0.5, 0.5, "no data", ha="center", va="center")
        ax.set_title(title)
        return
    y = np.arange(len(df))[::-1]
    for i, row in df.iterrows():
        kind = row["kind"]
        color = "#0c63b3" if kind == "cohort" else "#9c27b0"
        marker = "s" if kind == "cohort" else "D"
        size = 80 if kind == "cohort" else 120
        ax.scatter([row["HR"]], [y[i]], color=color, s=size, marker=marker, zorder=3)
        if np.isfinite(row["ci_lo"]) and np.isfinite(row["ci_hi"]):
            ax.plot([row["ci_lo"], row["ci_hi"]], [y[i], y[i]],
                     color=color, lw=1.6, zorder=2)
    ax.axvline(1.0, color="grey", lw=0.8, ls="--")
    ax.set_xscale("log")
    if xlim is None:
        # auto wide log range
        lo = max(0.001, np.nanmin(df["ci_lo"].values) * 0.7)
        hi = max(np.nanmax(df["ci_hi"].values) * 1.3, 100)
        xlim = (lo, hi)
    ax.set_xlim(*xlim)
    labels = []
    for _, row in df.iterrows():
        nlabel = ""
        if pd.notna(row["n"]):
            nlabel = f" (n={int(row['n'])}/{int(row['events'])})"
        hrtxt = f"  HR={row['HR']:.2f}"
        if np.isfinite(row['ci_lo']) and np.isfinite(row['ci_hi']):
            hrtxt += f" [{row['ci_lo']:.2f}, {row['ci_hi']:.2f}]"
        if pd.notna(row.get("I2")):
            hrtxt += f", I²={row['I2']:.0f}%"
        labels.append(f"{row['label']}{nlabel}{hrtxt}")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8.5, family="monospace")
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("Hazard Ratio (log scale)", fontsize=9)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)


fig, axes = plt.subplots(2, 2, figsize=(14, 9))
panel_plot(axes[0, 0],
           forest[forest["forest_panel"] == "DM2_vs_notDM_PFI_TCGA_drivers"],
           "(A) DM2 vs not_DM — PFI (TCGA driver-anchor sub-cohorts)",
           xlim=(0.005, 100))
panel_plot(axes[0, 1],
           forest[forest["forest_panel"] == "DM1xTERT_OS_pooled"],
           "(B) DM1×TERT comutant — OS (TCGA full + Landa, pooled)",
           xlim=(0.5, 500))
panel_plot(axes[1, 0],
           forest[forest["forest_panel"] == "DM1_vs_DM2_PFI_TCGA_drivers"],
           "(C) DM1 vs DM2 (protective) — PFI (TCGA)",
           xlim=(0.005, 100))
panel_plot(axes[1, 1],
           forest[forest["forest_panel"] == "DM2_vs_notDM_OS_TCGA_drivers"],
           "(D) DM2 vs not_DM — OS (sensitivity)",
           xlim=(0.005, 200))

fig.suptitle("R1 — Meta-survival forest plots (Paper 1+2 BRAF Nature sprint)",
              fontsize=12, fontweight="bold")
plt.tight_layout(rect=[0, 0, 1, 0.97])
plt.savefig(ROOT / "r1_forest.png", dpi=160, bbox_inches="tight")
plt.savefig(ROOT / "r1_forest.pdf", bbox_inches="tight")
print(f"wrote r1_forest.png / .pdf")
