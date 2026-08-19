#!/usr/bin/env python3
"""Figure 34 — Korean three-cohort panel z comparison standalone.

Combines K2 PRJEB11591, Lee 2024 GSE213647, GSE286332 Korean cohorts into a
single side-by-side comparison. Total Korean n = 910 cases. Shows panel z
silenced fraction is consistent across all three Korean datasets.
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
LEE = Path("/home/seungho/personal/THCA_data_analysis/project/results/v17_korean/GSE213647_panel_score.tsv")
GSE286 = Path("/home/seungho/personal/THCA_data_analysis/project/results/r17_tcga_panel_d4p2_reconciliation/cross_ethnic/gse286332_two_axis.tsv")
OUT_PNG = ROOT / "results" / "figures" / "figure34_korean_three_cohort.png"
OUT_PDF = ROOT / "results" / "figures" / "figure34_korean_three_cohort.pdf"

PANEL = ["SLC5A5","TPO","TG","TSHR","PAX8","NKX2-1","FOXE1","DIO1"]


def k2_panel_z():
    """K2 PRJEB11591 within-cohort z."""
    k2 = pd.read_csv(K2, sep="\t")
    avail = [g for g in PANEL if g in k2.columns]
    x = np.log10(k2[avail].astype(float).clip(lower=1e-3) + 1.0)
    z = x.subtract(x.mean(axis=0), axis=1).div(x.std(axis=0).replace(0, np.nan), axis=1)
    return z.mean(axis=1).values


def lee_panel_z(hist=None):
    lee = pd.read_csv(LEE, sep="\t")
    if hist is not None:
        lee = lee[lee["histology"].isin(hist)]
    return lee["panel_z"].dropna().values


def gse286_panel_z():
    df = pd.read_csv(GSE286, sep="\t")
    # panel_silencing > 0 = silenced; flip to align with panel_z convention
    return (-df["panel_silencing"]).dropna().values


def main():
    k2_z = k2_panel_z()
    lee_t = lee_panel_z(hist=["PTC", "PDFP", "UTC/ATC"])
    lee_n = lee_panel_z(hist=["Normal"])
    gse286_z = gse286_panel_z()

    cohorts = [
        ("K2 PRJEB11591\n(n = 260, all benign + malignant)", k2_z, "#888888"),
        ("Lee 2024 GSE213647 Normal\n(n = 262)", lee_n, "#3a6694"),
        ("Lee 2024 GSE213647 Tumour\n(n = 370)", lee_t, "#a04e48"),
        ("GSE286332 PTC ± HT\n(n = 18)", gse286_z, "#5a4470"),
    ]
    print([(name, len(v), float(np.median(v))) for name, v, _ in cohorts])

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.6), facecolor=IVORY,
                              gridspec_kw=dict(wspace=0.30, left=0.07, right=0.97, top=0.84, bottom=0.16))
    fig.suptitle("Korean three-cohort panel z comparison (total Korean n = 910 cases)",
                 fontsize=13, fontweight="bold", color=INK, y=0.97)

    # Left: stacked violin
    ax = axes[0]
    data = [v for _, v, _ in cohorts]
    colors = [c for _, _, c in cohorts]
    labels = [name for name, _, _ in cohorts]
    parts = ax.violinplot(data, positions=range(len(cohorts)), widths=0.78,
                           showmeans=False, showmedians=True)
    for body, color in zip(parts["bodies"], colors):
        body.set_facecolor(color); body.set_alpha(0.6); body.set_edgecolor(color)
    parts["cmedians"].set_color(INK); parts["cmedians"].set_lw(1.5)
    for k in ("cmins", "cmaxes", "cbars"):
        parts[k].set_color(MUTED)
    rng = np.random.default_rng(7)
    for i, vals in enumerate(data):
        if len(vals) > 300:
            sel = rng.choice(vals, size=300, replace=False)
        else:
            sel = vals
        jit = rng.uniform(-0.10, 0.10, len(sel))
        ax.scatter(np.full(len(sel), i) + jit, sel, color=colors[i], s=10,
                   alpha=0.5, edgecolor="none", zorder=2)
    ax.axhline(0, color=MUTED, ls="--", lw=0.7)
    ax.set_xticks(range(len(cohorts)))
    ax.set_xticklabels(labels, fontsize=8.5, rotation=10)
    ax.set_ylabel("panel z  (within-cohort z mean)")
    panel_title(ax, "Panel z distribution by Korean cohort"); panel_letter(ax, "a")

    # Right: median bar
    ax = axes[1]
    medians = [float(np.median(v)) for _, v, _ in cohorts]
    bars = ax.bar(range(len(cohorts)), medians, color=colors, alpha=0.85, edgecolor="white", lw=0.7)
    for xi, m in zip(range(len(cohorts)), medians):
        ax.text(xi, m + (0.04 if m >= 0 else -0.04), f"{m:+.2f}",
                ha="center", va="bottom" if m >= 0 else "top",
                fontsize=11, fontweight="bold", color=INK)
    ax.axhline(0, color=MUTED, ls="--", lw=0.7)
    ax.set_xticks(range(len(cohorts)))
    ax.set_xticklabels([n.split("\n")[0].replace(" ", "\n", 1) for n, _, _ in cohorts], fontsize=8.5)
    ax.set_ylabel("median panel z")
    panel_title(ax, "Median panel z by cohort"); panel_letter(ax, "b")

    fig.text(0.5, 0.04,
             "Korean tumour cohorts collapse to lower panel z than normal-rich cohorts — graded silencing recapitulated in 3 independent Korean datasets.",
             ha="center", fontsize=9.5, color=INK, style="italic")
    fig.savefig(OUT_PNG, dpi=220, facecolor=IVORY)
    fig.savefig(OUT_PDF, facecolor=IVORY)
    plt.close(fig)
    print(f"wrote {OUT_PNG}\nwrote {OUT_PDF}")


if __name__ == "__main__":
    main()
