#!/usr/bin/env python3
"""Figure 20 — Lu 2023 sc-RNA: per-sample zone fractions stratified by histology.

Stacked bar of R17 zone fraction per Lu 2023 sample (n = 14,624 malignant cells
across n samples) grouped by histology (PTC, FVPTC, ATC). Shows ATC samples
are 99% silenced (BRAF-like + dark-matter).
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
SAMP = Path("/home/seungho/personal/THCA_data_analysis/project/results/r17_tcga_panel_d4p2_reconciliation/sc/lu2023_sc_zone_fraction_by_sample.tsv")
HIST = Path("/home/seungho/personal/THCA_data_analysis/project/results/r17_tcga_panel_d4p2_reconciliation/sc/lu2023_sc_zone_fraction_by_histology.tsv")
OUT_PNG = ROOT / "results" / "figures" / "figure20_lu2023_sc_by_histology.png"
OUT_PDF = ROOT / "results" / "figures" / "figure20_lu2023_sc_by_histology.pdf"


ZONE_ORDER_FULL = [
    "wild-type-like (preserved, no HT)",
    "RAS-like zone (preserved + HT)",
    "BRAF-like zone (silenced, no HT)",
    "dark-matter zone (silenced + HT)",
]
ZONE_COLOR = {
    "wild-type-like (preserved, no HT)": "#cdc6c0",
    "RAS-like zone (preserved + HT)": RED,
    "BRAF-like zone (silenced, no HT)": BLUE,
    "dark-matter zone (silenced + HT)": "#5a4470",
}
SHORT = {
    "wild-type-like (preserved, no HT)": "WT-like",
    "RAS-like zone (preserved + HT)": "RAS-like",
    "BRAF-like zone (silenced, no HT)": "BRAF-like",
    "dark-matter zone (silenced + HT)": "dark-matter",
}


def main():
    samp = pd.read_csv(SAMP, sep="\t")
    hist = pd.read_csv(HIST, sep="\t")
    print(samp.head()); print("---"); print(hist)

    # Make sure all 4 zone columns exist
    zone_cols = [c for c in ZONE_ORDER_FULL if c in samp.columns]
    samp = samp.dropna(subset=zone_cols, how="all").copy()
    # Sort samples by histology then by dark-matter fraction
    samp["__dm"] = samp["dark-matter zone (silenced + HT)"] if "dark-matter zone (silenced + HT)" in samp.columns else 0
    samp["__hist_order"] = samp["histology"].map({"PTC": 0, "FVPTC": 1, "ATC": 2}).fillna(3)
    samp = samp.sort_values(["__hist_order", "__dm"])

    fig, axes = plt.subplots(1, 2, figsize=(15, 5.6), facecolor=IVORY,
                              gridspec_kw=dict(wspace=0.30, left=0.06, right=0.97, top=0.84, bottom=0.18, width_ratios=[1.4, 1]))
    fig.suptitle("Lu 2023 single-cell malignant — R17 zone fraction per sample × histology",
                 fontsize=13, fontweight="bold", color=INK, y=0.97)

    # Left: per-sample stacked bar
    ax = axes[0]
    x = np.arange(len(samp))
    bottom = np.zeros(len(samp))
    for c in zone_cols:
        vals = samp[c].fillna(0).values / 100.0  # convert pct → fraction
        ax.bar(x, vals, bottom=bottom, color=ZONE_COLOR[c], label=SHORT[c], width=0.85,
               edgecolor="white", lw=0.3)
        bottom += vals
    ax.set_xticks(x); ax.set_xticklabels(samp["sample"].values, fontsize=7, rotation=70, ha="right")
    ax.set_ylabel("R17 zone fraction")
    ax.set_ylim(0, 1.02)
    # Group dividers
    hg = samp["histology"].values
    for i in range(1, len(hg)):
        if hg[i] != hg[i-1]:
            ax.axvline(i - 0.5, color=INK, lw=0.6, alpha=0.6)
    # Group labels above the chart
    for g in samp["histology"].unique():
        idx = np.where(samp["histology"].values == g)[0]
        if len(idx) == 0: continue
        cx = (idx[0] + idx[-1]) / 2
        ax.text(cx, 1.07, g, ha="center", fontsize=11, color=INK, fontweight="bold")
    ax.set_ylim(0, 1.12)
    ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=4, frameon=False)
    panel_title(ax, "Per-sample zone fractions sorted by histology"); panel_letter(ax, "a")

    # Right: histology-level summary
    ax = axes[1]
    hist_mat = hist.set_index("histology")[zone_cols]
    # Normalize to fraction if needed
    if hist_mat.values.max() > 1.5: hist_mat = hist_mat / 100.0
    hist_mat = hist_mat.reindex(["PTC", "FVPTC", "ATC"]).dropna(how="all")
    x = np.arange(len(hist_mat))
    bottom = np.zeros(len(hist_mat))
    for c in zone_cols:
        vals = hist_mat[c].values
        ax.bar(x, vals, bottom=bottom, color=ZONE_COLOR[c], label=SHORT[c], width=0.55,
               edgecolor="white", lw=0.6)
        # in-bar labels
        for xi, v, b in zip(x, vals, bottom):
            if v > 0.06:
                ax.text(xi, b + v/2, f"{v*100:.0f}%", ha="center", va="center", fontsize=10,
                        color="white", fontweight="bold")
        bottom += vals
    ax.set_xticks(x); ax.set_xticklabels(hist_mat.index, fontsize=11)
    ax.set_ylabel("R17 zone fraction")
    ax.set_ylim(0, 1.02)
    panel_title(ax, "Pooled per-histology zone fraction"); panel_letter(ax, "b")

    fig.text(0.5, 0.04,
             "ATC malignant cells are 99.3% silenced (61% BRAF-like + 38.3% dark-matter); "
             "PTC cells are predominantly WT-like (72.2%).",
             ha="center", fontsize=9.5, color=INK, style="italic")
    fig.savefig(OUT_PNG, dpi=220, facecolor=IVORY)
    fig.savefig(OUT_PDF, facecolor=IVORY)
    plt.close(fig)
    print(f"wrote {OUT_PNG}\nwrote {OUT_PDF}")


if __name__ == "__main__":
    main()
