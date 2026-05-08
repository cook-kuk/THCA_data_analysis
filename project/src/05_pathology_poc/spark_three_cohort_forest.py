#!/usr/bin/env python3
"""
F6 update: 3-cohort consistency forest plot (TCGA + Korean + GSE250521 ST).

Show CAF/M1-M2/TLS/CD8 × RAI signature replication across 3 platforms with
n=572 + n=632 + n=16 slides, plus driver-NEG and DM1/DM2 stratified slices.
"""
from __future__ import annotations
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent.parent.parent
RES = ROOT / "project/results/03_pathology_poc"
OUT = RES / "summary_figs/F6_three_cohort_forest.png"

# Pre-computed values
DATA = {
    "CAF":   {"TCGA n=572": -0.313, "Korean n=632": -0.344, "ST n=16 (lag)": -0.275,
              "TCGA driver-NEG n=132": -0.285, "TCGA DM1 n=84": -0.183, "TCGA DM2 n=61": +0.030},
    "M1/M2": {"TCGA n=572": -0.474, "Korean n=632": -0.470, "ST n=16 (lag)": -0.119,
              "TCGA driver-NEG n=132": -0.479, "TCGA DM1 n=84": -0.193, "TCGA DM2 n=61": -0.241},
    "TLS":   {"TCGA n=572": -0.304, "Korean n=632": -0.144, "ST n=16 (lag)": -0.095,
              "TCGA driver-NEG n=132": -0.253, "TCGA DM1 n=84": -0.045, "TCGA DM2 n=61": +0.091},
    "CD8":   {"TCGA n=572": -0.154, "Korean n=632": +0.010, "ST n=16 (lag)": -0.006,
              "TCGA driver-NEG n=132": -0.267, "TCGA DM1 n=84": +0.050, "TCGA DM2 n=61": -0.012},
}

COHORT_ORDER = ["TCGA n=572", "Korean n=632", "ST n=16 (lag)",
                "TCGA driver-NEG n=132", "TCGA DM1 n=84", "TCGA DM2 n=61"]
COHORT_COLOR = {
    "TCGA n=572":            "#ff7c3e",
    "Korean n=632":          "#5eead4",
    "ST n=16 (lag)":         "#fbbf24",
    "TCGA driver-NEG n=132": "#ef4444",
    "TCGA DM1 n=84":         "#a78bfa",
    "TCGA DM2 n=61":         "#9ca3af",
}


def main():
    modules = list(DATA.keys())
    fig, axes = plt.subplots(1, len(modules), figsize=(16, 5.5), constrained_layout=True, sharey=True)
    for ax, mod in zip(axes, modules):
        cohorts = COHORT_ORDER
        rs = [DATA[mod][c] for c in cohorts]
        colors = [COHORT_COLOR[c] for c in cohorts]
        y = np.arange(len(cohorts))[::-1]  # top-to-bottom
        ax.barh(y, rs, color=colors, edgecolor="white", height=0.6)
        ax.axvline(0, color="#333", lw=0.7)
        ax.set_xlim(-0.55, 0.15)
        ax.set_title(f"{mod} × RAI lineage", fontsize=11, fontweight="bold")
        ax.set_xlabel("Spearman r")
        for yi, (c, r) in zip(y, zip(cohorts, rs)):
            sig = " ★" if r < -0.25 else ""
            ax.text(r - 0.02 if r < 0 else r + 0.02, yi, f"{r:+.3f}{sig}",
                    ha="right" if r < 0 else "left", va="center", fontsize=9)
            ax.text(-0.55, yi, c, ha="left", va="center", fontsize=8.5, color="#444")
        ax.set_yticks([])
        ax.spines[["top", "right"]].set_visible(False)

    fig.suptitle("F6 · 3-cohort × subset cross-platform replication of CCI × RAI avoidance signature\n"
                 "ALL 6 strata negative for CAF & M1/M2 except TCGA DM2 (CAF +0.030 — DM1-specific finding)",
                 fontsize=12, fontweight="bold")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
