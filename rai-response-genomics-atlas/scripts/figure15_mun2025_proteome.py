#!/usr/bin/env python3
"""Figure 15 — Mun 2025 proteome (n = 336) zone fraction by histology group.

Stacked bar of R17 zone fraction stratified by histology group (Normal, Benign,
PTC, ATC) using pre-computed mun2025_zone_fraction_by_group.tsv.
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
ZONE_TSV = Path("/home/seungho/personal/THCA_data_analysis/project/results/r17_tcga_panel_d4p2_reconciliation/proteome_mun2025/mun2025_zone_fraction_by_group.tsv")
TWO_AXIS = Path("/home/seungho/personal/THCA_data_analysis/project/results/r17_tcga_panel_d4p2_reconciliation/proteome_mun2025/mun2025_two_axis.tsv")
OUT_PNG = ROOT / "results" / "figures" / "figure15_mun2025_proteome_zone.png"
OUT_PDF = ROOT / "results" / "figures" / "figure15_mun2025_proteome_zone.pdf"

ZONE_ORDER = ["WT-like", "RAS-like", "BRAF-like", "dark-matter"]
ZONE_COLOR = {"WT-like": "#cdc6c0", "RAS-like": RED, "BRAF-like": BLUE, "dark-matter": "#5a4470"}


def main():
    zf = pd.read_csv(ZONE_TSV, sep="\t")
    # Detect possible column names
    print(zf.head())
    print(zf.columns.tolist())

    # Try to detect a group column and zone columns / row format
    if "group" in zf.columns and "zone" in zf.columns:
        pivot = zf.pivot_table(index="group", columns="zone", values="fraction", aggfunc="first").fillna(0)
    else:
        pivot = zf.set_index(zf.columns[0])

    # Convert percentages to fractions if max > 1.5 (sum looks like 100)
    if pivot.values.max() > 1.5:
        pivot = pivot / 100.0

    # Order rows by least → most dedifferentiated
    row_order = [r for r in ["PTC", "PDTC", "ATC"] if r in pivot.index]
    if row_order:
        pivot = pivot.loc[row_order]

    # Reorder columns
    cols = [c for c in ZONE_ORDER if c in pivot.columns]
    other = [c for c in pivot.columns if c not in cols]
    pivot = pivot[cols + other]

    fig, axes = plt.subplots(1, 2, figsize=(14, 6.2), facecolor=IVORY,
                              gridspec_kw=dict(wspace=0.28, left=0.07, right=0.97, top=0.84, bottom=0.22))
    fig.suptitle("Mun 2025 proteome — R17 zone fraction by histology group (n = 336)",
                 fontsize=13.5, fontweight="bold", color=INK, y=0.97)

    # Left: stacked bar
    ax = axes[0]
    bottom = np.zeros(len(pivot))
    x = np.arange(len(pivot))
    for c in pivot.columns:
        col = ZONE_COLOR.get(c, "#777777")
        ax.bar(x, pivot[c].values, bottom=bottom, color=col, edgecolor="white", lw=0.6,
               label=c, alpha=0.92)
        # in-bar fraction labels
        for xi, v, b in zip(x, pivot[c].values, bottom):
            if v > 0.07:
                ax.text(xi, b + v/2, f"{v*100:.0f}%", ha="center", va="center",
                        fontsize=9, color="white", fontweight="bold")
        bottom += pivot[c].values

    ax.set_xticks(x); ax.set_xticklabels(pivot.index, fontsize=10)
    ax.set_ylabel("R17 zone fraction")
    ax.set_ylim(0, 1.02)
    ax.legend(fontsize=8.5, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=4, frameon=False)
    panel_title(ax, "Zone fraction by histology group  ·  ATC dark-matter enrichment"); panel_letter(ax, "a")

    # Right: per-sample two-axis scatter
    ta = pd.read_csv(TWO_AXIS, sep="\t")
    ax = axes[1]
    panel_col = next((c for c in ta.columns if c == "panel_silencing"), None)
    ht_col    = next((c for c in ta.columns if c == "HT_overlap"), None)
    group_col = next((c for c in ta.columns if c == "group"), None)
    if panel_col and ht_col:
        if group_col:
            for g, sub in ta.groupby(group_col):
                ax.scatter(sub[panel_col], sub[ht_col], s=15, alpha=0.7,
                           label=f"{g} (n={len(sub)})",
                           edgecolor="white", lw=0.4)
        else:
            ax.scatter(ta[panel_col], ta[ht_col], s=15, alpha=0.7, color=BLUE, edgecolor="white", lw=0.4)
        ax.axhline(0, color=MUTED, lw=0.7, ls="--")
        ax.axvline(0, color=MUTED, lw=0.7, ls="--")
        ax.set_xlabel(f"protein panel silencing  ({panel_col})")
        ax.set_ylabel(f"protein HT-overlap  ({ht_col})")
        if group_col:
            ax.legend(fontsize=8, loc="best", frameon=True)
        # quadrant labels — panel_silencing > 0 = silenced (right), HT_overlap > 0 = HT immune (top)
        xl, xh = ax.get_xlim(); yl, yh = ax.get_ylim()
        ax.text(xh*0.95, yh*0.95, "dark-matter", ha="right", va="top", fontsize=9, color="#5a4470", fontweight="bold")
        ax.text(xl*0.95, yh*0.95, "RAS-like", ha="left", va="top", fontsize=9, color=RED, fontweight="bold")
        ax.text(xh*0.95, yl*0.95, "BRAF-like", ha="right", va="bottom", fontsize=9, color=BLUE, fontweight="bold")
        ax.text(xl*0.95, yl*0.95, "WT-like", ha="left", va="bottom", fontsize=9, color="#777777", fontweight="bold")
        panel_title(ax, "Per-sample protein two-axis decomposition"); panel_letter(ax, "b")
    else:
        ax.text(0.5, 0.5, "two-axis columns not found", ha="center", va="center", transform=ax.transAxes)

    fig.text(0.5, 0.025,
             "ATC group shows dark-matter zone enrichment OR = 8.54, p = 2.7 × 10⁻¹⁵ "
             "(Fisher exact, ATC vs PTC).",
             ha="center", fontsize=9.5, color=INK, style="italic")
    fig.savefig(OUT_PNG, dpi=220, facecolor=IVORY)
    fig.savefig(OUT_PDF, facecolor=IVORY)
    plt.close(fig)
    print(f"wrote {OUT_PNG}\nwrote {OUT_PDF}")


if __name__ == "__main__":
    main()
