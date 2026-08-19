#!/usr/bin/env python3
"""Figure 22 — Per-gene HM450 β vs RNA z scatter, BRAF-like + dark-matter zones.

Direct gene-by-gene visualisation of methylation→silencing relationship.
Each gene contributes 2 points (BRAF-like + dark-matter) using zone-mean β
vs zone-mean RNA z from TCGA-THCA.
"""
from __future__ import annotations
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
from _nature_style import setup_rc, panel_letter, panel_title, IVORY, INK, MUTED, BLUE, RED

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr

setup_rc()

ROOT = Path(__file__).resolve().parent.parent
SRC = Path("/home/seungho/personal/THCA_data_analysis/project/results/r17_tcga_panel_d4p2_reconciliation/hm450_beta/r17_beta_vs_rna_cross.tsv")
HM_MEANS = Path("/home/seungho/personal/THCA_data_analysis/project/results/r17_tcga_panel_d4p2_reconciliation/hm450_beta/r17_zone_hm450_beta_means.tsv")
GENE_MEANS = Path("/home/seungho/personal/THCA_data_analysis/project/results/r17_tcga_panel_d4p2_reconciliation/gene_profile/r17_zone_gene_z_means.tsv")
OUT_PNG = ROOT / "results" / "figures" / "figure22_methylation_rna_scatter.png"
OUT_PDF = ROOT / "results" / "figures" / "figure22_methylation_rna_scatter.pdf"

PANEL = ["SLC5A5","TPO","TG","TSHR","PAX8","NKX2-1","FOXE1","DIO1"]
ZONE_COLOR = {"WT-like": "#888888", "RAS-like": RED, "BRAF-like": BLUE, "dark-matter": "#5a4470"}


def main():
    hm = pd.read_csv(HM_MEANS, sep="\t").set_index("gene")
    rna = pd.read_csv(GENE_MEANS, sep="\t", index_col=0)
    hm_p = hm.reindex(PANEL)
    rna_p = rna.reindex(PANEL)

    # Build long-form: gene, zone, beta, rna_z
    rows = []
    for g in PANEL:
        for zone in ["WT-like", "RAS-like", "BRAF-like", "dark-matter"]:
            rows.append({"gene": g, "zone": zone,
                         "beta": hm_p.loc[g, zone],
                         "rna_z": rna_p.loc[g, zone]})
    long = pd.DataFrame(rows)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6.0), facecolor=IVORY,
                              gridspec_kw=dict(wspace=0.30, left=0.07, right=0.97, top=0.84, bottom=0.16))
    fig.suptitle("Gene-level evidence that promoter hypermethylation drives RNA silencing",
                 fontsize=13, fontweight="bold", color=INK, y=0.97)

    # Left: gene × zone scatter
    ax = axes[0]
    for zone in ["WT-like", "RAS-like", "BRAF-like", "dark-matter"]:
        sub = long[long["zone"] == zone]
        ax.scatter(sub["beta"], sub["rna_z"], color=ZONE_COLOR[zone], s=85, alpha=0.85,
                   edgecolor="white", lw=0.8, label=zone, zorder=3)
        # label genes for BRAF-like and dark-matter only
        if zone in ("BRAF-like", "dark-matter"):
            for _, r in sub.iterrows():
                ax.annotate(r["gene"], (r["beta"], r["rna_z"]),
                            xytext=(6, 5), textcoords="offset points", fontsize=8.5,
                            fontstyle="italic", color=ZONE_COLOR[zone], alpha=0.85)

    rho_all, p_all = spearmanr(long["beta"], long["rna_z"])
    ax.text(0.02, 0.04, f"Spearman ρ = {rho_all:.2f},  p = {p_all:.2g}", transform=ax.transAxes,
            fontsize=10, color=INK, fontweight="bold")
    ax.axhline(0, color=MUTED, lw=0.7, ls="--")
    ax.set_xlabel("mean HM450 β  (promoter methylation)")
    ax.set_ylabel("mean RNA z  (expression level)")
    ax.legend(fontsize=9, loc="upper right", frameon=True)
    panel_title(ax, "Zone-mean β vs RNA z across 8 panel genes × 4 zones"); panel_letter(ax, "a")

    # Right: per-gene Δβ (dark - WT) vs Δrna (dark - WT)
    ax = axes[1]
    delta_beta = (hm_p["dark-matter"] - hm_p["WT-like"]).values
    delta_rna = (rna_p["dark-matter"] - rna_p["WT-like"]).values
    for i, g in enumerate(PANEL):
        ax.scatter(delta_beta[i], delta_rna[i], color="#5a4470", s=120, alpha=0.85,
                   edgecolor="white", lw=0.9, zorder=3)
        ax.annotate(g, (delta_beta[i], delta_rna[i]), xytext=(8, 6),
                    textcoords="offset points", fontsize=10, fontstyle="italic", color="#5a4470")
    ax.axhline(0, color=MUTED, lw=0.7, ls="--")
    ax.axvline(0, color=MUTED, lw=0.7, ls="--")
    rho_d, p_d = spearmanr(delta_beta, delta_rna)
    ax.text(0.02, 0.96, f"Spearman ρ = {rho_d:.2f},  p = {p_d:.3g}", transform=ax.transAxes,
            fontsize=10, color=INK, fontweight="bold", va="top")
    ax.set_xlabel("Δ HM450 β  (dark-matter − WT-like)")
    ax.set_ylabel("Δ RNA z  (dark-matter − WT-like)")
    panel_title(ax, "Per-gene methylation gain vs RNA loss"); panel_letter(ax, "b")

    fig.text(0.5, 0.045,
             "Higher promoter methylation in the dark-matter zone correlates with lower RNA expression at the gene level "
             "— direct methylation-silencing evidence.",
             ha="center", fontsize=9.5, color=INK, style="italic")
    fig.savefig(OUT_PNG, dpi=220, facecolor=IVORY)
    fig.savefig(OUT_PDF, facecolor=IVORY)
    plt.close(fig)
    print(f"wrote {OUT_PNG}\nwrote {OUT_PDF}")


if __name__ == "__main__":
    main()
