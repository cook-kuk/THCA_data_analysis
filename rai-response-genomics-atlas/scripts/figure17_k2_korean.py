#!/usr/bin/env python3
"""Figure 17 — K2 PRJEB11591 Korean cohort (n = 260) panel z distribution.

Uses pre-computed K2 predictions table from the parent project. Plots panel z
distribution by histology category and by DM call.
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
K2 = Path("/home/seungho/personal/THCA_data_analysis/project/results/v17_korean/K2_korean_predictions_v4.tsv")
HIST = Path("/home/seungho/personal/THCA_data_analysis/project/results/v17_korean/K2_histology_breakdown.tsv")
OUT_PNG = ROOT / "results" / "figures" / "figure17_k2_korean.png"
OUT_PDF = ROOT / "results" / "figures" / "figure17_k2_korean.pdf"


def main():
    k2 = pd.read_csv(K2, sep="\t")
    print(k2.columns.tolist()); print(k2.head(2))
    print("n =", len(k2))
    # Compute panel z as mean of log10(FPKM/TPM)
    panel = ["SLC5A5","TPO","TG","TSHR","PAX8","NKX2-1","FOXE1","DIO1"]
    avail = [g for g in panel if g in k2.columns]
    if not avail:
        print("no panel cols")
        return
    # Use log10 + 1, then z within cohort
    x = np.log10(k2[avail].astype(float).clip(lower=1e-3) + 1.0)
    z = x.subtract(x.mean(axis=0), axis=1).div(x.std(axis=0).replace(0, np.nan), axis=1)
    panel_z = z.mean(axis=1)
    k2["panel_z"] = panel_z.values
    # Bring histology if available
    if "category_clean" not in k2.columns and HIST.exists():
        # Match via run? Most K2 tables include just run IDs; histology breakdown is at category level.
        pass

    # DM call distribution
    fig = plt.figure(figsize=(13, 5.4), facecolor=IVORY)
    gs = fig.add_gridspec(1, 2, hspace=0.3, wspace=0.30, left=0.07, right=0.97, top=0.84, bottom=0.18)
    fig.suptitle("K2 PRJEB11591 Korean cohort (n = 260) — 8-gene panel z distribution",
                 fontsize=13, fontweight="bold", color=INK, y=0.97)

    # Panel a: panel z hist
    ax = fig.add_subplot(gs[0, 0])
    ax.hist(panel_z.dropna(), bins=30, color=BLUE, alpha=0.85, edgecolor="white", lw=0.6)
    ax.axvline(panel_z.median(), color=RED, lw=1.4, ls="--", label=f"median = {panel_z.median():.2f}")
    ax.set_xlabel("8-gene panel z (within-cohort)")
    ax.set_ylabel("count")
    ax.legend(fontsize=9, loc="upper right", frameon=True)
    panel_title(ax, "Panel z distribution (n = 260)"); panel_letter(ax, "a")

    # Panel b: DM call counts
    ax = fig.add_subplot(gs[0, 1])
    if "DM_call" in k2.columns:
        counts = k2["DM_call"].value_counts().sort_index()
        colors = {"DM1": BLUE, "DM2": RED, "DM2_strong": RED, "intermediate": "#7e689a"}
        bar_colors = [colors.get(c, "#888") for c in counts.index]
        bars = ax.bar(counts.index, counts.values, color=bar_colors, alpha=0.9, edgecolor="white", lw=0.8)
        for b, v in zip(bars, counts.values):
            ax.text(b.get_x() + b.get_width()/2, v + 4, f"{v}\n({v/len(k2)*100:.1f}%)",
                    ha="center", fontsize=10, color=INK, fontweight="bold")
        ax.set_ylabel("sample count")
        ax.set_ylim(0, max(counts.values)*1.18)
        panel_title(ax, "DM (silencing) call distribution"); panel_letter(ax, "b")
    else:
        ax.text(0.5, 0.5, "DM_call column not found", ha="center", va="center", transform=ax.transAxes)

    fig.text(0.5, 0.045,
             "K2 NBNR substratum (DM-silenced) is the Korean cohort equivalent of TCGA dark-matter — "
             "see paired memory v17_D6P7_dm1_subB_NBNR.",
             ha="center", fontsize=9.5, color=INK, style="italic")
    fig.savefig(OUT_PNG, dpi=220, facecolor=IVORY)
    fig.savefig(OUT_PDF, facecolor=IVORY)
    plt.close(fig)
    print(f"wrote {OUT_PNG}\nwrote {OUT_PDF}")


if __name__ == "__main__":
    main()
