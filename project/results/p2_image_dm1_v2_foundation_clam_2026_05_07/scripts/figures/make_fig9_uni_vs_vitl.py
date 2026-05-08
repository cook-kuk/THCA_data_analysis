#!/usr/bin/env python3
"""
Paper 2 Figure 9: UNI vs ViT-L head-to-head comparison.

Two panels:
  Left  - per-fold AUC bar chart (5 folds × 2 encoders)
  Right - bootstrap CI + boxplot distribution comparison
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(
    "/home/seungho/personal/THCA_data_analysis/project/results/"
    "p2_image_dm1_v2_foundation_clam_2026_05_07"
)
OUT_DIR = ROOT / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Per-fold AUCs
VITL_FOLDS = [0.71, 0.72, 1.00, 0.71, 1.00]
UNI_FOLDS = [1.000, 1.000, 0.964, 1.000, 0.792]

# Pooled / bootstrap
VITL_POOL = 0.746
VITL_CI = (0.611, 0.862)
UNI_POOL = 0.874
UNI_CI = (0.760, 0.967)

KILL = 0.70

VITL_COLOR = "#7a7a7a"
UNI_COLOR = "#5e3b9c"


def main():
    fig, (axL, axR) = plt.subplots(
        1, 2, figsize=(13.5, 5.6), dpi=300,
        gridspec_kw={"width_ratios": [1.15, 1.0]},
    )

    # ---------------- LEFT: per-fold bars
    folds = np.arange(1, 6)
    width = 0.36
    x = folds.astype(float)

    bars_v = axL.bar(
        x - width / 2, VITL_FOLDS, width,
        color=VITL_COLOR, edgecolor="black", lw=0.7,
        label="ViT-L (ImageNet-21k)",
    )
    bars_u = axL.bar(
        x + width / 2, UNI_FOLDS, width,
        color=UNI_COLOR, edgecolor="black", lw=0.7,
        label="UNI (Mass-100K, H&E)",
    )

    for b, v in zip(bars_v, VITL_FOLDS):
        axL.text(
            b.get_x() + b.get_width() / 2, v + 0.012,
            f"{v:.2f}", ha="center", va="bottom", fontsize=9.0,
            color="#333",
        )
    for b, v in zip(bars_u, UNI_FOLDS):
        axL.text(
            b.get_x() + b.get_width() / 2, v + 0.012,
            f"{v:.3f}", ha="center", va="bottom", fontsize=9.0,
            color=UNI_COLOR, weight="bold",
        )

    # kill-switch line
    axL.axhline(KILL, color="#c0392b", lw=1.2, ls="--", zorder=1)
    axL.text(
        5.55, KILL, " kill-switch  0.70",
        ha="left", va="center", color="#c0392b", fontsize=9, weight="bold",
    )

    axL.set_xticks(folds)
    axL.set_xticklabels([f"Fold {i}" for i in folds], fontsize=10)
    axL.set_ylabel("AUC (held-out fold)", fontsize=11)
    axL.set_ylim(0, 1.12)
    axL.set_xlim(0.4, 6.2)
    axL.set_yticks(np.arange(0, 1.05, 0.1))
    axL.grid(True, axis="y", alpha=0.25, ls=":")
    axL.set_axisbelow(True)
    axL.legend(loc="lower right", fontsize=9.5, frameon=True)
    axL.set_title("A.  Per-fold AUC — ViT-L vs UNI",
                  fontsize=12, weight="bold", loc="left")

    # ---------------- RIGHT: bootstrap CI + per-fold boxplot
    # Sub-axis: top = CI bars, bottom = box of fold AUCs
    # We'll do two stacked items by manual placement.
    # Use single axR but with two y-rows separated.

    # CI bars (horizontal-style, plotted as vertical errorbars + dot)
    encoders = ["ViT-L", "UNI"]
    pools = [VITL_POOL, UNI_POOL]
    cis = [VITL_CI, UNI_CI]
    colors = [VITL_COLOR, UNI_COLOR]

    pos_ci = [1.0, 2.2]  # x positions for CI markers
    for px, p, ci, col in zip(pos_ci, pools, cis, colors):
        lo, hi = ci
        axR.errorbar(
            px, p, yerr=[[p - lo], [hi - p]],
            fmt="o", color=col, ecolor=col,
            elinewidth=2.4, capsize=8, capthick=2.0,
            markersize=11, markeredgecolor="black", markeredgewidth=0.8,
            zorder=3,
        )
        axR.text(
            px, hi + 0.025,
            f"{p:.3f}\n[{lo:.2f}, {hi:.2f}]",
            ha="center", va="bottom", fontsize=9.0, color=col,
            weight="bold",
        )

    # boxplots offset to the right
    pos_box = [3.6, 4.8]
    bp = axR.boxplot(
        [VITL_FOLDS, UNI_FOLDS],
        positions=pos_box, widths=0.55,
        patch_artist=True, showmeans=True,
        meanprops=dict(marker="D", markerfacecolor="white",
                       markeredgecolor="black", markersize=7),
    )
    for patch, col in zip(bp["boxes"], colors):
        patch.set_facecolor(col)
        patch.set_alpha(0.55)
        patch.set_edgecolor("black")
    for med in bp["medians"]:
        med.set_color("black")
        med.set_linewidth(1.5)

    # scatter raw fold points
    rng = np.random.default_rng(0)
    for px, vals, col in zip(pos_box, [VITL_FOLDS, UNI_FOLDS], colors):
        jitter = rng.uniform(-0.10, 0.10, size=len(vals))
        axR.scatter(
            np.full(len(vals), px) + jitter, vals,
            color=col, edgecolor="black", s=42, zorder=4, alpha=0.9,
        )

    # mean ± SD annotation
    v_mean = np.mean(VITL_FOLDS)
    v_sd = np.std(VITL_FOLDS, ddof=1)
    u_mean = np.mean(UNI_FOLDS)
    u_sd = np.std(UNI_FOLDS, ddof=1)
    axR.text(
        pos_box[0], 0.07,
        f"{v_mean:.2f} ± {v_sd:.2f}",
        ha="center", va="top", fontsize=9.0, color=VITL_COLOR, weight="bold",
    )
    axR.text(
        pos_box[1], 0.07,
        f"{u_mean:.2f} ± {u_sd:.2f}",
        ha="center", va="top", fontsize=9.0, color=UNI_COLOR, weight="bold",
    )

    # kill-switch line
    axR.axhline(KILL, color="#c0392b", lw=1.0, ls="--", zorder=1, alpha=0.6)

    axR.set_xticks(pos_ci + pos_box)
    axR.set_xticklabels(
        ["ViT-L\npooled", "UNI\npooled",
         "ViT-L\nfolds", "UNI\nfolds"],
        fontsize=9.5,
    )
    axR.set_ylabel("AUC", fontsize=11)
    axR.set_ylim(0, 1.12)
    axR.set_xlim(0.3, 5.6)
    axR.set_yticks(np.arange(0, 1.05, 0.1))
    axR.grid(True, axis="y", alpha=0.25, ls=":")
    axR.set_axisbelow(True)

    # Section divider
    axR.axvline(2.95, color="#bbb", lw=0.8, ls=":")
    axR.text(1.6, 1.075, "Bootstrap 95% CI",
             ha="center", fontsize=9.5, style="italic", color="#444")
    axR.text(4.2, 1.075, "Per-fold dist.",
             ha="center", fontsize=9.5, style="italic", color="#444")

    axR.set_title(
        "B.  UNI provides Δ AUC = +0.128 over ImageNet ViT-L",
        fontsize=12, weight="bold", loc="left",
    )

    fig.suptitle(
        "Figure 9.  UNI vs ViT-L head-to-head — Phase 2 TCGA-THCA H&E → DM1",
        fontsize=13.5, weight="bold", y=1.00,
    )

    plt.tight_layout(rect=(0, 0, 1, 0.95))

    png = OUT_DIR / "fig9_uni_vs_vitl_comparison.png"
    pdf = OUT_DIR / "fig9_uni_vs_vitl_comparison.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"WROTE {png}")
    print(f"WROTE {pdf}")


if __name__ == "__main__":
    main()
