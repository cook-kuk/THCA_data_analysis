#!/usr/bin/env python3
"""Figure 13 — Mechanistic anchor: HM450 β vs RNA z, per-gene × per-zone.

Side-by-side heatmaps showing methylation hypermethylation matches RNA silencing
across the 8 panel genes and 4 R17 zones. Reads pre-computed tables from the
parent project's TCGA reconciliation module.
"""
from __future__ import annotations
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
from _nature_style import setup_rc, panel_letter, panel_title, IVORY, INK, MUTED, BLUE, RED

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

setup_rc()

ROOT = Path(__file__).resolve().parent.parent
HM450 = Path("/home/seungho/personal/THCA_data_analysis/project/results/r17_tcga_panel_d4p2_reconciliation/hm450_beta/r17_zone_hm450_beta_means.tsv")
GENE  = Path("/home/seungho/personal/THCA_data_analysis/project/results/r17_tcga_panel_d4p2_reconciliation/gene_profile/r17_zone_gene_z_means.tsv")
OUT_PNG = ROOT / "results" / "figures" / "figure13_methylation_rna_mechanism.png"
OUT_PDF = ROOT / "results" / "figures" / "figure13_methylation_rna_mechanism.pdf"

ZONE_ORDER = ["WT-like", "RAS-like", "BRAF-like", "dark-matter"]
PANEL_GENES = ["SLC5A5", "TPO", "TG", "TSHR", "PAX8", "NKX2-1", "FOXE1", "DIO1"]


def main():
    hm = pd.read_csv(HM450, sep="\t").set_index("gene")
    rna = pd.read_csv(GENE, sep="\t", index_col=0)

    # Restrict to panel genes and order both
    hm_panel = hm.reindex(PANEL_GENES)[ZONE_ORDER]
    rna_panel = rna.reindex(PANEL_GENES)[ZONE_ORDER]

    fig, axes = plt.subplots(1, 2, figsize=(15, 6.5), facecolor=IVORY,
                              gridspec_kw=dict(wspace=0.35, left=0.10, right=0.96, top=0.83, bottom=0.18))
    fig.suptitle("HM450 methylation β and RNA z per panel gene × R17 zone (TCGA-THCA, n = 518)",
                 fontsize=13.5, fontweight="bold", color=INK, y=0.96)

    # HM450 β heatmap (blue scale, low → high)
    cmap_b = LinearSegmentedColormap.from_list("beta", ["#fafafa", "#a8c5dc", "#386a99", "#1f4d7a"], N=200)
    ax = axes[0]
    im = ax.imshow(hm_panel.values, cmap=cmap_b, aspect="auto", vmin=0, vmax=hm_panel.values.max()*1.05)
    for i in range(hm_panel.shape[0]):
        for j in range(hm_panel.shape[1]):
            val = hm_panel.iloc[i, j]
            col = "white" if val > hm_panel.values.max()*0.55 else INK
            ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=9, color=col)
    ax.set_yticks(range(len(PANEL_GENES))); ax.set_yticklabels(PANEL_GENES, fontstyle="italic")
    ax.set_xticks(range(len(ZONE_ORDER))); ax.set_xticklabels(ZONE_ORDER, rotation=15, ha="right")
    cb = fig.colorbar(im, ax=ax, shrink=0.75)
    cb.set_label("mean HM450 β", fontsize=9)
    panel_title(ax, "HM450 β  ·  promoter methylation"); panel_letter(ax, "a")

    # RNA z heatmap (red-blue diverging)
    cmap_r = LinearSegmentedColormap.from_list("rna", ["#1f4d7a", "#7ba0c5", "#fafafa", "#cb857f", "#7c322e"], N=200)
    ax = axes[1]
    vmax = max(abs(rna_panel.values).max(), 1.2)
    im = ax.imshow(rna_panel.values, cmap=cmap_r, aspect="auto", vmin=-vmax, vmax=vmax)
    for i in range(rna_panel.shape[0]):
        for j in range(rna_panel.shape[1]):
            val = rna_panel.iloc[i, j]
            col = "white" if abs(val) > vmax*0.6 else INK
            ax.text(j, i, f"{val:+.2f}", ha="center", va="center", fontsize=9, color=col)
    ax.set_yticks(range(len(PANEL_GENES))); ax.set_yticklabels(PANEL_GENES, fontstyle="italic")
    ax.set_xticks(range(len(ZONE_ORDER))); ax.set_xticklabels(ZONE_ORDER, rotation=15, ha="right")
    cb = fig.colorbar(im, ax=ax, shrink=0.75)
    cb.set_label("mean RNA z", fontsize=9)
    panel_title(ax, "RNA z  ·  expression level"); panel_letter(ax, "b")

    fig.text(0.5, 0.045,
             "Dark-matter zone shows the highest mean β and the lowest mean RNA z across all panel genes — "
             "promoter hypermethylation directly tracks the silencing seen at the transcript level.",
             ha="center", fontsize=9.5, color=INK, style="italic")
    fig.savefig(OUT_PNG, dpi=220, facecolor=IVORY)
    fig.savefig(OUT_PDF, facecolor=IVORY)
    plt.close(fig)
    print(f"wrote {OUT_PNG}\nwrote {OUT_PDF}")


if __name__ == "__main__":
    main()
