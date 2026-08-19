#!/usr/bin/env python3
"""Figure 27 — Lu 2023 sc-RNA per-cell panel_silencing violin by histology.

Per-cell distribution of panel silencing across PTC, FVPTC, ATC at the
single-cell resolution.
"""
from __future__ import annotations
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
from _nature_style import setup_rc, panel_letter, panel_title, IVORY, INK, MUTED, BLUE, RED

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import kruskal

setup_rc()

ROOT = Path(__file__).resolve().parent.parent
LU = Path("/home/seungho/personal/THCA_data_analysis/project/results/r17_tcga_panel_d4p2_reconciliation/sc/lu2023_sc_two_axis_malignant.tsv")
OUT_PNG = ROOT / "results" / "figures" / "figure27_lu2023_per_cell_violin.png"
OUT_PDF = ROOT / "results" / "figures" / "figure27_lu2023_per_cell_violin.pdf"


def main():
    df = pd.read_csv(LU, sep="\t")
    print(df["histology"].value_counts())
    # Build panel z = - panel_silencing (so positive = preserved, matches other figs)
    df["panel_z"] = -df["panel_silencing"]
    order = ["PTC", "FVPTC", "ATC"]
    order = [g for g in order if (df["histology"] == g).sum() > 0]
    palette = {"PTC": BLUE, "FVPTC": "#7e689a", "ATC": "#7c322e"}

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.6), facecolor=IVORY,
                              gridspec_kw=dict(wspace=0.32, left=0.07, right=0.97, top=0.84, bottom=0.16))
    fig.suptitle("Lu 2023 single-cell — per-cell panel z by histology (n = 14,624 malignant cells)",
                 fontsize=13, fontweight="bold", color=INK, y=0.97)

    # Left: violin
    ax = axes[0]
    data = [df.loc[df["histology"] == g, "panel_z"].dropna().values for g in order]
    parts = ax.violinplot(data, positions=range(len(order)), showmeans=False, showmedians=True,
                           widths=0.78)
    for body, g in zip(parts["bodies"], order):
        body.set_facecolor(palette[g]); body.set_alpha(0.55); body.set_edgecolor(palette[g])
    parts["cmedians"].set_color(INK); parts["cmedians"].set_lw(1.4)
    parts["cmins"].set_color(MUTED); parts["cmaxes"].set_color(MUTED); parts["cbars"].set_color(MUTED)
    # Overlay quartile lines
    for i, vals in enumerate(data):
        if len(vals) == 0: continue
        q1, q3 = np.percentile(vals, [25, 75])
        ax.hlines([q1, q3], i-0.15, i+0.15, color=INK, lw=0.8, alpha=0.7)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([f"{g}\n(cells={int((df['histology']==g).sum()):,})" for g in order], fontsize=10)
    ax.set_ylabel("per-cell panel z  (preserved → silenced flipped)")
    ax.axhline(0, color=MUTED, ls="--", lw=0.7)
    h, p = kruskal(*data)
    ax.text(0.02, 0.04, f"Kruskal–Wallis H = {h:.1f}, p = {p:.3g}",
            transform=ax.transAxes, fontsize=10, color=INK, fontweight="bold")
    panel_title(ax, "Per-cell panel z violin"); panel_letter(ax, "a")

    # Right: zone fraction stacked bar (re-used)
    ax = axes[1]
    counts = df.groupby(["histology", "zone"]).size().unstack(fill_value=0)
    counts = counts.reindex(order)
    fracs = counts.div(counts.sum(axis=1), axis=0)
    # Order zone columns visually
    zone_short = {
        "wild-type-like (preserved, no HT)": "WT-like",
        "RAS-like zone (preserved + HT)": "RAS-like",
        "BRAF-like zone (silenced, no HT)": "BRAF-like",
        "dark-matter zone (silenced + HT)": "dark-matter",
    }
    zone_colors = {"WT-like": "#cdc6c0", "RAS-like": RED, "BRAF-like": BLUE, "dark-matter": "#5a4470"}
    zone_order_full = list(zone_short.keys())
    zone_order_full = [z for z in zone_order_full if z in fracs.columns]
    bottom = np.zeros(len(order))
    x = np.arange(len(order))
    for c in zone_order_full:
        short = zone_short[c]
        vals = fracs[c].values
        ax.bar(x, vals, bottom=bottom, color=zone_colors[short], label=short, width=0.55,
               edgecolor="white", lw=0.6)
        for xi, v, b in zip(x, vals, bottom):
            if v > 0.05:
                ax.text(xi, b + v/2, f"{v*100:.0f}%", ha="center", va="center",
                        fontsize=10, color="white", fontweight="bold")
        bottom += vals
    ax.set_xticks(x); ax.set_xticklabels(order, fontsize=10)
    ax.set_ylabel("zone fraction")
    ax.set_ylim(0, 1.02)
    ax.legend(fontsize=8.5, loc="upper center", bbox_to_anchor=(0.5, -0.10), ncol=4, frameon=False)
    panel_title(ax, "Per-cell zone fraction by histology"); panel_letter(ax, "b")

    fig.text(0.5, 0.04,
             "ATC malignant cells (n = 2,749) collapse onto silenced subset; PTC cells (n = 4,510) dominate WT-like.",
             ha="center", fontsize=9.5, color=INK, style="italic")
    fig.savefig(OUT_PNG, dpi=220, facecolor=IVORY)
    fig.savefig(OUT_PDF, facecolor=IVORY)
    plt.close(fig)
    print(f"wrote {OUT_PNG}\nwrote {OUT_PDF}")


if __name__ == "__main__":
    main()
