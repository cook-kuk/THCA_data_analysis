#!/usr/bin/env python3
"""Figure 21 — Lee 2024 GSE213647 Korean cohort (n = 632) panel z by histology.

Per-sample panel z stratified by histology category (Normal, PTC, PDFP, UTC/ATC)
in Korean Lee et al. 2024 cohort.
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
SRC = Path("/home/seungho/personal/THCA_data_analysis/project/results/v17_korean/GSE213647_panel_score.tsv")
OUT_PNG = ROOT / "results" / "figures" / "figure21_lee2024_korean_histology.png"
OUT_PDF = ROOT / "results" / "figures" / "figure21_lee2024_korean_histology.pdf"


def main():
    df = pd.read_csv(SRC, sep="\t")
    print(df["histology"].value_counts())
    order = ["Normal", "PTC", "PDFP", "UTC/ATC"]
    order = [g for g in order if (df["histology"] == g).sum() > 0]
    palette = {"Normal": "#888888", "PTC": BLUE, "PDFP": "#7e689a", "UTC/ATC": RED}

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.4), facecolor=IVORY,
                              gridspec_kw=dict(wspace=0.30, left=0.07, right=0.97, top=0.84, bottom=0.16))
    fig.suptitle("Lee 2024 GSE213647 Korean cohort (n = 632) — panel z by histology",
                 fontsize=13, fontweight="bold", color=INK, y=0.97)

    # Left: boxplot + jitter
    ax = axes[0]
    data = [df.loc[df["histology"] == g, "panel_z"].dropna().values for g in order]
    bp = ax.boxplot(data, positions=range(len(order)), widths=0.55,
                     patch_artist=True, medianprops=dict(color=INK, lw=1.3),
                     boxprops=dict(lw=0.7), whiskerprops=dict(lw=0.7), capprops=dict(lw=0.7),
                     flierprops=dict(marker="o", markersize=2, alpha=0.4))
    for patch, g in zip(bp["boxes"], order):
        patch.set_facecolor(palette[g]); patch.set_alpha(0.5); patch.set_edgecolor(palette[g])
    rng = np.random.default_rng(13)
    for i, vals in enumerate(data):
        jit = rng.uniform(-0.10, 0.10, len(vals))
        ax.scatter(np.full(len(vals), i)+jit, vals, color=palette[order[i]],
                   s=8, alpha=0.45, edgecolor="none", zorder=2)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([f"{g}\n(n={int((df['histology']==g).sum())})" for g in order], fontsize=10)
    ax.set_ylabel("8-gene panel z")
    ax.axhline(0, color=MUTED, ls="--", lw=0.7)
    h, p = kruskal(*data)
    ax.text(0.02, 0.04, f"Kruskal–Wallis H = {h:.1f}, p = {p:.3g}", transform=ax.transAxes,
            fontsize=10, color=INK, fontweight="bold")
    panel_title(ax, "Panel z by histology  ·  graded silencing"); panel_letter(ax, "a")

    # Right: median per histology bar
    ax = axes[1]
    medians = [np.median(v) for v in data]
    bars = ax.bar(range(len(order)), medians, color=[palette[g] for g in order],
                  alpha=0.85, edgecolor="white", lw=0.8)
    for xi, m, vals in zip(range(len(order)), medians, data):
        ax.text(xi, m + (0.05 if m >= 0 else -0.05), f"{m:+.2f}",
                ha="center", va="bottom" if m >= 0 else "top", fontsize=10, color=INK, fontweight="bold")
    ax.axhline(0, color=MUTED, ls="--", lw=0.7)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, fontsize=10)
    ax.set_ylabel("median panel z")
    ax.set_ylim(min(medians)*1.3 if min(medians) < 0 else -0.2, max(0.2, max(medians)*1.4))
    panel_title(ax, "Median panel z by histology"); panel_letter(ax, "b")

    fig.text(0.5, 0.04,
             "Korean cohort recapitulates the graded silencing axis: Normal (+0.56) → PTC (-0.24) → PDFP (-0.46) → UTC/ATC (-1.42).",
             ha="center", fontsize=9.5, color=INK, style="italic")
    fig.savefig(OUT_PNG, dpi=220, facecolor=IVORY)
    fig.savefig(OUT_PDF, facecolor=IVORY)
    plt.close(fig)
    print(f"wrote {OUT_PNG}\nwrote {OUT_PDF}")


if __name__ == "__main__":
    main()
