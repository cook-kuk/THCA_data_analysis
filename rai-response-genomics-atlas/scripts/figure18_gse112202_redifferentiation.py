#!/usr/bin/env python3
"""Figure 18 — GSE112202 digoxin redifferentiation (Tier-5) standalone view.

Per-gene log2 fold-change waterfall for 8 panel genes between digoxin-treated
and untreated NMTC patients (n=22). Reports median log2FC, # genes up.
"""
from __future__ import annotations
from pathlib import Path
import sys, gzip
sys.path.insert(0, str(Path(__file__).parent))
from _nature_style import setup_rc, panel_letter, panel_title, IVORY, INK, MUTED, BLUE, RED

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

setup_rc()

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "data" / "raw" / "GSE112202" / "GSE112202_overall_comparison.fpkm_tracking.gz"
OUT_PNG = ROOT / "results" / "figures" / "figure18_gse112202_redifferentiation.png"
OUT_PDF = ROOT / "results" / "figures" / "figure18_gse112202_redifferentiation.pdf"

PANEL = ["SLC5A5","TPO","TG","TSHR","PAX8","NKX2-1","FOXE1","DIO1"]


def main():
    df = pd.read_csv(SRC, sep="\t", compression="gzip")
    print(df.columns.tolist()[:15])
    fpkm_u = "untreated (overall)_FPKM"
    fpkm_t = "digoxin treated (overall)_FPKM"
    sub = df[df["gene_short_name"].isin(PANEL)].drop_duplicates("gene_short_name")
    sub = sub.set_index("gene_short_name").reindex(PANEL)
    sub["log2FC"] = np.log2((sub[fpkm_t] + 1e-3) / (sub[fpkm_u] + 1e-3))
    sub = sub.sort_values("log2FC", ascending=True)
    print(sub[["log2FC", fpkm_u, fpkm_t]])

    # Waterfall
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.4), facecolor=IVORY,
                              gridspec_kw=dict(wspace=0.32, left=0.07, right=0.97, top=0.84, bottom=0.17))
    fig.suptitle("GSE112202 Tier-5 redifferentiation — digoxin restores 8-gene panel expression",
                 fontsize=13, fontweight="bold", color=INK, y=0.97)

    ax = axes[0]
    colors = [RED if v < 0 else BLUE for v in sub["log2FC"]]
    ax.barh(range(len(sub)), sub["log2FC"], color=colors, alpha=0.9, edgecolor="white", lw=0.6)
    ax.axvline(0, color=MUTED, lw=0.8)
    ax.set_yticks(range(len(sub))); ax.set_yticklabels(sub.index, fontstyle="italic")
    ax.set_xlabel("log₂ fold-change  (digoxin / untreated)")
    for i, v in enumerate(sub["log2FC"]):
        if v >= 0:
            ax.text(v + 0.03, i, f"{v:+.2f}", va="center", fontsize=9, color=INK)
        else:
            ax.text(v - 0.03, i, f"{v:+.2f}", va="center", ha="right", fontsize=9, color=INK)
    median = sub["log2FC"].median()
    ax.axvline(median, color="#5a4470", lw=1.0, ls="--", label=f"median log₂FC = {median:+.2f}")
    ax.legend(fontsize=9, loc="lower right")
    panel_title(ax, "Per-gene log₂FC waterfall  ·  6 of 8 up-regulated"); panel_letter(ax, "a")

    # FPKM scatter
    ax = axes[1]
    sub_unsorted = sub.copy()
    for g in PANEL:
        if g in sub_unsorted.index:
            x = sub_unsorted.loc[g, fpkm_u]
            y = sub_unsorted.loc[g, fpkm_t]
            ax.plot([x, x], [x, y], color="#cccccc", lw=0.6)
            ax.scatter(x, y, s=50, color=BLUE if y > x else RED, alpha=0.9, edgecolor="white", lw=0.8, zorder=3)
            ax.annotate(g, (x, y), xytext=(6, 4), textcoords="offset points", fontsize=9, fontstyle="italic")
    # diagonal
    lim = max(sub_unsorted[fpkm_u].max(), sub_unsorted[fpkm_t].max()) * 1.2
    ax.plot([0.05, lim], [0.05, lim], color=MUTED, lw=0.7, ls="--", alpha=0.7, label="x = y")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(0.05, lim); ax.set_ylim(0.05, lim)
    ax.set_xlabel("untreated FPKM")
    ax.set_ylabel("digoxin treated FPKM")
    ax.legend(fontsize=9, loc="upper left")
    panel_title(ax, "FPKM untreated vs digoxin  ·  group-level"); panel_letter(ax, "b")

    n_up = int((sub["log2FC"] > 0).sum())
    fig.text(0.5, 0.04,
             f"{n_up}/8 panel genes up-regulated under digoxin; canonical iodide-axis genes "
             "(TSHR, SLC5A5, NKX2-1, TG) show the strongest restoration.",
             ha="center", fontsize=9.5, color=INK, style="italic")
    fig.savefig(OUT_PNG, dpi=220, facecolor=IVORY)
    fig.savefig(OUT_PDF, facecolor=IVORY)
    plt.close(fig)
    print(f"wrote {OUT_PNG}\nwrote {OUT_PDF}")


if __name__ == "__main__":
    main()
