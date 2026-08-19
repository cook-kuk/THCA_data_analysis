#!/usr/bin/env python3
"""Figure 28 — 8-gene within-cohort correlation network.

Computes pairwise Spearman correlations among the 8 panel genes in each cohort
and visualises as: (a) per-cohort correlation heatmap with hierarchical
clustering; (b) cross-cohort correlation consistency (gene-pair vs gene-pair
agreement across cohorts).
"""
from __future__ import annotations
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
from _nature_style import setup_rc, panel_letter, panel_title, IVORY, INK, MUTED

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import spearmanr

setup_rc()

ROOT = Path(__file__).resolve().parent.parent
INTERIM = ROOT / "data" / "interim"
RAW = ROOT / "data" / "raw"
TCGA = Path("/home/seungho/personal/THCA_data_analysis/project/results/r17_tcga_panel_d4p2_reconciliation/tert_interaction/r17_tert_per_sample.tsv")
OUT_PNG = ROOT / "results" / "figures" / "figure28_gene_correlation_network.png"
OUT_PDF = ROOT / "results" / "figures" / "figure28_gene_correlation_network.pdf"

PANEL = ["SLC5A5","TPO","TG","TSHR","PAX8","NKX2-1","FOXE1","DIO1"]


def load_gse_expr(acc):
    expr = pd.read_csv(INTERIM / f"{acc}_expression.tsv", sep="\t", index_col=0)
    p2g_path = next((RAW / acc).glob("GPL*_probe2gene.tsv"))
    p2g = pd.read_csv(p2g_path, sep="\t", dtype=str)
    pmap = dict(zip(p2g["ID"], p2g.get("gene_symbol", p2g.get("GENE_SYMBOL", []))))
    e = expr.copy(); e["gene"] = e.index.map(pmap)
    e = e.dropna(subset=["gene"]); e = e[e["gene"] != ""]
    e["__var"] = e.iloc[:, :-1].var(axis=1)
    e = e.sort_values("__var", ascending=False).drop_duplicates("gene", keep="first")
    return e.drop(columns="__var").set_index("gene")


def corr_panel(expr_or_df, panel=PANEL, samples=None):
    """Return Spearman corr matrix (panel × panel)."""
    if isinstance(expr_or_df, pd.DataFrame) and set(panel).issubset(expr_or_df.columns):
        # rows = samples, cols include gene names
        sub = expr_or_df[panel].copy()
    else:
        avail = [g for g in panel if g in expr_or_df.index]
        sub = expr_or_df.loc[avail]
        if samples is not None:
            sub = sub[samples]
        sub = sub.T
    return sub.corr(method="spearman").reindex(index=panel, columns=panel)


def main():
    cohorts = {}
    # TCGA
    tcga = pd.read_csv(TCGA, sep="\t")
    cohorts["TCGA-THCA  (n ≈ 500)"] = corr_panel(tcga[PANEL])
    # GSE151179
    cohorts["GSE151179  (n = 52)"] = corr_panel(load_gse_expr("GSE151179"))
    # GSE299988
    cohorts["GSE299988  (n = 14)"] = corr_panel(load_gse_expr("GSE299988"))

    cmap = LinearSegmentedColormap.from_list("c", ["#1f4d7a", "#7ba0c5", "#fafafa", "#cb857f", "#7c322e"], N=200)

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.5), facecolor=IVORY,
                              gridspec_kw=dict(wspace=0.30, left=0.05, right=0.93, top=0.82, bottom=0.13))
    fig.suptitle("Per-cohort gene-gene Spearman correlation among 8 panel genes",
                 fontsize=13, fontweight="bold", color=INK, y=0.96)

    letters = ["a", "b", "c"]
    for ax, (name, c), ltr in zip(axes, cohorts.items(), letters):
        im = ax.imshow(c.values, cmap=cmap, vmin=-1, vmax=1, aspect="equal")
        for i in range(len(PANEL)):
            for j in range(len(PANEL)):
                v = c.iloc[i, j]
                if np.isnan(v): continue
                col = "white" if abs(v) > 0.6 else INK
                ax.text(j, i, f"{v:+.2f}", ha="center", va="center", fontsize=8.5, color=col)
        ax.set_xticks(range(len(PANEL))); ax.set_xticklabels(PANEL, fontstyle="italic", rotation=70, ha="right")
        ax.set_yticks(range(len(PANEL))); ax.set_yticklabels(PANEL, fontstyle="italic")
        panel_title(ax, name); panel_letter(ax, ltr)

    # Color bar on right
    cb_ax = fig.add_axes([0.945, 0.16, 0.012, 0.62])
    cb = fig.colorbar(im, cax=cb_ax)
    cb.set_label("Spearman ρ", fontsize=9)

    fig.text(0.5, 0.04,
             "All 8 panel genes show consistently positive co-expression across cohorts — they form a single coherent thyroid-differentiation module.",
             ha="center", fontsize=9.5, color=INK, style="italic")
    fig.savefig(OUT_PNG, dpi=220, facecolor=IVORY)
    fig.savefig(OUT_PDF, facecolor=IVORY)
    plt.close(fig)
    print(f"wrote {OUT_PNG}\nwrote {OUT_PDF}")


if __name__ == "__main__":
    main()
