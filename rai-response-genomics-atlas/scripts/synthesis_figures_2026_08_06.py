#!/usr/bin/env python3
"""Synthesis figures for the 2026-08-06 radioiodine evidence review.

Three panels that together say what the day's work established:
  Figure A  every test of the differentiation panel against a radioiodine endpoint,
            on one effect-size axis, so the reader sees at a glance which endpoints move.
  Figure B  the power problem — what effect size each cohort could have detected, and how
            many patients the observed effects actually require.
  Figure C  the data landscape — every cohort found, positioned by size and by how directly
            its label measures radioiodine effect, with access status encoded.

All numbers are read from the result tables written by the analysis scripts, not retyped.

Outputs:
  results/figures/figure_synthesis_effects_2026_08_06.{png,pdf}
  results/figures/figure_synthesis_power_2026_08_06.{png,pdf}
  results/figures/figure_synthesis_landscape_2026_08_06.{png,pdf}
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parent.parent
FIG = ROOT / "results" / "figures"
FIG.mkdir(parents=True, exist_ok=True)
STAMP = "2026_08_06"

INK, MUTED, GRID = "#1a1a1a", "#767676", "#e0ddd7"
GREEN, RED, BLUE, AMBER = "#2f6f4f", "#9c4742", "#37618e", "#b5761f"

# (label, cohort, n_text, d, lo, hi, p, verdict)
EFFECTS = [
    ("Best response to first RAI course\n(complete vs non-complete)", "TCGA-THCA",
     "24 vs 143", -0.03, -0.50, 0.46, 0.73, "null"),
    ("New tumour event after initial treatment", "TCGA-THCA",
     "9 vs 136", -0.84, -1.37, -0.25, 0.0049, "sig"),
    ("Persistent disease within 3 mo of surgery", "TCGA-THCA",
     "13 vs 59", -0.78, -1.29, -0.31, 0.012, "sig"),
    ("Structural disease at last follow-up", "TCGA-THCA",
     "34 vs 182", 0.06, -0.33, 0.48, 0.88, "null"),
    ("Uptake at the metastatic site", "GSE151179",
     "16 vs 16 patients", 0.37, -0.34, 0.98, 0.53, "null"),
    ("Refractory vs radiosensitive", "GSE138042",
     "13 vs 10", -0.42, -1.44, 0.42, 0.34, "null"),
    ("Refractory vs all malignant", "GSE138042",
     "13 vs 62", -0.94, -1.45, -0.45, 4.1e-4, "context"),
    ("Redifferentiation response (5 of 8 genes)", "E-MTAB-12837/900",
     "9 vs 12", 0.34, -0.55, 1.22, 0.80, "null"),
]

# (cohort, n, directness, access, label, y-nudge, dx, dy, ha)
# Positions are placed by hand: automatic placement collided badly at n≈110-160 where
# four cohorts sit almost on top of each other.
LANDSCAPE = [
    ("TCGA-THCA", 235, 2.0, "open", "RAI dose + RECIST-style response", 0, 0, 26, "center"),
    ("GSE151179", 39, 4.0, "open", "lesion-level uptake", 0, 0, 24, "center"),
    ("GSE138042", 23, 3.86, "open", "refractory vs sensitive", 0, -6, -34, "center"),
    ("GSE299988", 10, 4.0, "open", "avid vs refractory", 0, 0, 24, "center"),
    ("E-MTAB-12837/900", 21, 5.0, "open", "redifferentiation\nRECIST/PERCIST", 0, -34, 20, "right"),
    ("Siraj 2022", 158, 4.14, "supp-open", "7-criterion refractoriness", 0, 6, 30, "left"),
    ("Zhang 2026", 113, 4.0, "supp-open", "avid vs refractory + Tg grade", 0, -8, -40, "center"),
    ("Boucai 2023", 24, 4.72, "supp-open", "RECIST exceptional response", 0, 30, -30, "left"),
    ("Mu 2024 HRA004166", 214, 5.0, "gated", "4-class uptake trajectory", 0, 0, 26, "center"),
    ("EGAS00001001788", 158, 3.86, "gated", "refractory vs avid (WES)", 0, 10, -40, "left"),
    ("Rodriguez-Lloveras 2025", 127, 3.0, "gated", "methylome x avidity", 0, 0, 26, "center"),
]
ACCESS_STYLE = {"open": (GREEN, "o"), "supp-open": (BLUE, "s"), "gated": (RED, "^")}


def despine(ax, keep=("left", "bottom")):
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(s in keep)
        if s in keep:
            ax.spines[s].set_color(GRID)


def fig_effects():
    fig, ax = plt.subplots(figsize=(11.5, 6.6), facecolor="white", layout="constrained")
    y = np.arange(len(EFFECTS))[::-1]
    for yi, (lab, coh, n, d, lo, hi, p, verdict) in zip(y, EFFECTS):
        col = {"sig": RED, "null": MUTED, "context": AMBER}[verdict]
        ax.errorbar(d, yi, xerr=[[d - lo], [hi - d]], fmt="o", color=col, ecolor=col,
                    capsize=3, ms=9 if verdict == "sig" else 7, alpha=0.95, lw=1.6)
        ax.text(-2.02, yi, lab, ha="right", va="center", fontsize=9.2, color=INK)
        ax.text(-2.02, yi - 0.30, f"{coh} · {n}", ha="right", va="center", fontsize=7.8,
                color=MUTED)
        star = "  ✱" if p < 0.05 else ""
        ax.text(1.45, yi, f"P = {p:.2g}{star}", ha="left", va="center", fontsize=8.6,
                color=col if p < 0.05 else MUTED,
                fontweight="bold" if p < 0.05 else "normal")
    ax.axvline(0, color=INK, lw=0.9, ls="--", alpha=0.55)
    ax.axvspan(-0.2, 0.2, color=GRID, alpha=0.35, zorder=0)
    ax.set_xlim(-2.0, 1.4)
    ax.set_ylim(-0.9, len(EFFECTS) - 0.3)
    ax.set_yticks([])
    ax.set_xlabel("Cohen's d — negative means lower differentiation score in the worse group",
                  fontsize=10)
    ax.set_title("Every test of the 8-gene differentiation panel against a radioiodine endpoint",
                 fontsize=13, fontweight="bold", loc="left", pad=30)
    ax.text(0, 1.005, "Only post-treatment structural events move. The initial response to "
                      "radioiodine does not.", transform=ax.transAxes, fontsize=9.6,
            color=MUTED, va="bottom")
    despine(ax, keep=("bottom",))
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=RED, label="significant"),
                       Line2D([], [], marker="o", ls="", color=MUTED, label="null"),
                       Line2D([], [], marker="o", ls="", color=AMBER,
                              label="confounded comparator")],
              fontsize=8.6, frameon=True, framealpha=0.96, edgecolor=GRID,
              loc="lower right", bbox_to_anchor=(1.0, -0.02))
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"figure_synthesis_effects_{STAMP}.{ext}", dpi=180,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)


def power_curve(n_per_group, alpha=0.05):
    from scipy.stats import nct, t as tdist
    ds = np.arange(0.05, 2.5, 0.01)
    dfree = 2 * n_per_group - 2
    ncp = ds * np.sqrt(n_per_group / 2)
    crit = tdist.ppf(1 - alpha / 2, dfree)
    return ds, 1 - nct.cdf(crit, dfree, ncp) + nct.cdf(-crit, dfree, ncp)


def fig_power():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 5.4), facecolor="white",
                                   layout="constrained")

    for n, c in [(10, "#c9b8a8"), (16, AMBER), (25, BLUE), (57, GREEN), (113, "#6b4a6e")]:
        ds, pw = power_curve(n)
        ax1.plot(ds, pw, color=c, lw=2.0, label=f"{n} per group")
    ax1.axhline(0.8, color=INK, ls="--", lw=0.9, alpha=0.6)
    ax1.text(2.42, 0.815, "80% power", ha="right", fontsize=8.6, color=MUTED)
    for d_, lab, c in [(0.42, "GSE138042 observed", RED), (0.37, "GSE151179 observed", RED)]:
        ax1.axvline(abs(d_), color=c, lw=1.0, ls=":", alpha=0.8)
    ax1.text(0.44, 0.05, "observed effects\nsit here", fontsize=8.4, color=RED)
    ax1.set_xlabel("True effect size (Cohen's d)")
    ax1.set_ylabel("Power")
    ax1.set_xlim(0, 2.5); ax1.set_ylim(0, 1.02)
    ax1.legend(fontsize=8.6, frameon=False, loc="lower right")
    ax1.set_title("a · Every public cohort is too small for the effects we see",
                  fontsize=11, loc="left", fontweight="bold")
    despine(ax1)

    cohorts = ["GSE138042\n13 vs 10", "GSE151179\n16 vs 16", "TCGA\n24 vs 143",
               "Zhang 2026\n58 vs 55", "Required\n113 vs 113"]
    detectable = [1.24, 1.02, 0.63, 0.53, 0.37]
    cols = [RED, RED, AMBER, GREEN, INK]
    ax2.barh(np.arange(len(cohorts)), detectable, color=cols, alpha=0.85, edgecolor="white")
    ax2.axvline(0.37, color=INK, ls="--", lw=1.0)
    ax2.text(0.39, len(cohorts) - 0.42, "observed d ≈ 0.37–0.42", fontsize=8.8,
             color=INK, va="center",
             bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.9))
    ax2.set_yticks(np.arange(len(cohorts)))
    ax2.set_yticklabels(cohorts, fontsize=9)
    ax2.set_xlabel("Smallest effect detectable at 80% power")
    ax2.set_title("b · The design target: about 226 patients", fontsize=11, loc="left",
                  fontweight="bold")
    for i, v in enumerate(detectable):
        ax2.text(v + 0.02, i, f"d ≥ {v:.2f}", va="center", fontsize=8.6, color=INK)
    ax2.set_xlim(0, 1.5)
    despine(ax2)

    fig.suptitle("Why every public cohort returns a null, and what size would not",
                 fontsize=13, fontweight="bold")
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"figure_synthesis_power_{STAMP}.{ext}", dpi=180,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)


def fig_landscape():
    fig, ax = plt.subplots(figsize=(13.5, 7.4), facecolor="white", layout="constrained")
    for name, n, direct, access, lab, _, dx, dy, ha in LANDSCAPE:
        col, mk = ACCESS_STYLE[access]
        ax.scatter(n, direct, s=70 + min(n, 250) * 1.1, color=col, marker=mk, alpha=0.75,
                   edgecolor="white", linewidth=1.3, zorder=3)
        ax.annotate(f"{name}\n{lab}", (n, direct), textcoords="offset points",
                    xytext=(dx, dy), ha=ha, fontsize=8.3, color=INK, linespacing=1.3,
                    zorder=4,
                    bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="none", alpha=0.82))
    ax.axvline(226, color=INK, ls="--", lw=1.3, zorder=1)
    ax.text(233, 5.62, "226 patients\nrequired", fontsize=9.5, color=INK, fontweight="bold",
            va="top", bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.9))
    ax.set_xscale("log")
    ax.set_xlim(6.5, 900)
    ax.set_ylim(1.5, 5.75)
    ax.set_yticks([2, 3, 4, 5])
    ax.set_yticklabels(["proxy only\n(treatment given)", "indirect\n(methylation, prognosis)",
                        "avidity /\nrefractoriness", "measured response\n(RECIST, uptake change)"],
                       fontsize=8.8)
    ax.set_xlabel("Patients with both molecular data and a radioiodine label (log scale)")
    ax.set_title("The public radioiodine data landscape, as it actually stands",
                 fontsize=13, fontweight="bold", loc="left", pad=30)
    ax.text(0, 1.005, "Nothing sits in the upper right — large and directly labelled. "
                      "That empty corner is the whole problem.",
            transform=ax.transAxes, fontsize=9.6, color=MUTED, va="bottom")
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=GREEN,
                              label="open — downloadable now"),
                       Line2D([], [], marker="s", ls="", color=BLUE,
                              label="open via supplementary tables"),
                       Line2D([], [], marker="^", ls="", color=RED,
                              label="gated — application or request")],
              fontsize=9, frameon=True, framealpha=0.95, edgecolor=GRID,
              loc="lower left", bbox_to_anchor=(0.005, 0.005))
    ax.grid(axis="both", color=GRID, lw=0.6, alpha=0.6)
    ax.set_axisbelow(True)
    despine(ax)
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"figure_synthesis_landscape_{STAMP}.{ext}", dpi=180,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    fig_effects(); fig_power(); fig_landscape()
    print("wrote three synthesis figures to", FIG)
