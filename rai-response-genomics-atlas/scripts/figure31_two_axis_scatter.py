#!/usr/bin/env python3
"""Figure 31 — TCGA-THCA two-axis scatter (panel z × d4p2 HT-overlap), zone-colored.

Standalone view of the canonical two-axis decomposition that defines the
R17 4-zone partition.
"""
from __future__ import annotations
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
from _nature_style import setup_rc, panel_letter, panel_title, IVORY, INK, MUTED, BLUE, RED

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

setup_rc()

ROOT = Path(__file__).resolve().parent.parent
PER = Path("/home/seungho/personal/THCA_data_analysis/project/results/r17_tcga_panel_d4p2_reconciliation/tert_interaction/r17_tert_per_sample.tsv")
OUT_PNG = ROOT / "results" / "figures" / "figure31_two_axis_scatter.png"
OUT_PDF = ROOT / "results" / "figures" / "figure31_two_axis_scatter.pdf"

ZONE_COLOR = {"WT-like": "#888888", "RAS-like": RED, "BRAF-like": BLUE, "dark-matter": "#5a4470"}


def main():
    df = pd.read_csv(PER, sep="\t").dropna(subset=["RAI_8", "sig_score", "zone"])
    print(df["zone"].value_counts())

    fig, axes = plt.subplots(1, 2, figsize=(14, 6.0), facecolor=IVORY,
                              gridspec_kw=dict(wspace=0.30, left=0.07, right=0.96, top=0.84, bottom=0.16))
    fig.suptitle("TCGA-THCA two-axis decomposition (n = " + f"{len(df):,})",
                 fontsize=13, fontweight="bold", color=INK, y=0.96)

    # Left: scatter colored by zone
    ax = axes[0]
    for z, color in ZONE_COLOR.items():
        sub = df[df["zone"] == z]
        ax.scatter(sub["RAI_8"], sub["sig_score"], color=color, s=18, alpha=0.7,
                   edgecolor="white", lw=0.4, label=f"{z} (n = {len(sub)})")
    ax.axhline(0, color=INK, lw=0.6, ls="--", alpha=0.6)
    ax.axvline(0, color=INK, lw=0.6, ls="--", alpha=0.6)
    ax.set_xlabel("panel z (RAI_8)  →  preserved")
    ax.set_ylabel("d4p2 sig_score  →  HT-immune overlap")
    ax.legend(fontsize=8.5, loc="lower right", frameon=True)
    # Quadrant labels
    xl, xh = ax.get_xlim(); yl, yh = ax.get_ylim()
    ax.text(xh*0.95, yh*0.95, "RAS-like", ha="right", va="top", fontsize=10, color=RED, fontweight="bold", alpha=0.85)
    ax.text(xl*0.95, yh*0.95, "dark-matter", ha="left", va="top", fontsize=10, color="#5a4470", fontweight="bold", alpha=0.85)
    ax.text(xl*0.95, yl*0.95, "BRAF-like", ha="left", va="bottom", fontsize=10, color=BLUE, fontweight="bold", alpha=0.85)
    ax.text(xh*0.95, yl*0.95, "WT-like", ha="right", va="bottom", fontsize=10, color="#777777", fontweight="bold", alpha=0.85)
    panel_title(ax, "Two-axis space partitions 4 zones"); panel_letter(ax, "a")

    # Right: per-zone counts and proportions
    ax = axes[1]
    counts = df["zone"].value_counts().reindex(["WT-like", "RAS-like", "BRAF-like", "dark-matter"])
    pct = counts / counts.sum() * 100
    bars = ax.barh(counts.index, counts.values, color=[ZONE_COLOR[z] for z in counts.index],
                    edgecolor="white", alpha=0.9, lw=0.7)
    for b, c, p in zip(bars, counts, pct):
        ax.text(c + 4, b.get_y() + b.get_height()/2,
                f"{int(c):,}  ({p:.1f}%)", va="center", fontsize=10, color=INK, fontweight="bold")
    ax.set_xlabel("sample count")
    ax.set_xlim(0, counts.max()*1.30)
    panel_title(ax, "Zone composition"); panel_letter(ax, "b")

    fig.text(0.5, 0.04,
             "Two-axis decomposition (panel z × d4p2 HT-overlap) cleanly partitions TCGA-THCA into 4 zones; "
             "dark-matter (silenced + HT-overlap) is the highest-risk subgroup.",
             ha="center", fontsize=9.5, color=INK, style="italic")
    fig.savefig(OUT_PNG, dpi=220, facecolor=IVORY)
    fig.savefig(OUT_PDF, facecolor=IVORY)
    plt.close(fig)
    print(f"wrote {OUT_PNG}\nwrote {OUT_PDF}")


if __name__ == "__main__":
    main()
