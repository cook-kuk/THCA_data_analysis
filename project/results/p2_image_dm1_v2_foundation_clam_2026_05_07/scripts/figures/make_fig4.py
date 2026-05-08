#!/usr/bin/env python
"""Fig 4 — Phase 2 ROC curves (5-fold + pooled).

One panel: ROC curve for each of the 5 CV folds plus a bold pooled-prediction
ROC curve.  Diagonal dashed = chance.  Annotations include the mean +/- SD
fold AUC, the pooled AUC, and a horizontal reference at 0.55 (closure
ResNet50 baseline) for visual comparison on the y axis.

Inputs:
  - phase2_tcga_clam_QUICK/clam_per_slide_predictions.tsv

Outputs:
  - figures/fig4_phase2_roc.png
  - figures/fig4_phase2_roc.pdf
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, roc_curve

ROOT = Path(
    "/home/seungho/personal/THCA_data_analysis/project/results/"
    "p2_image_dm1_v2_foundation_clam_2026_05_07"
)
PHASE2 = ROOT / "phase2_tcga_clam_QUICK"
FIG_DIR = ROOT / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Wong color-blind palette
FOLD_COLORS = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9"]
POOLED_COLOR = "#000000"
BASELINE_COLOR = "#999999"

CLOSURE_BASELINE = 0.55  # ResNet50 closure reference


def main() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 11,
            "savefig.dpi": 300,
        }
    )

    df = pd.read_csv(PHASE2 / "clam_per_slide_predictions.tsv", sep="\t")
    folds = sorted(df["fold"].unique())

    fig, ax = plt.subplots(figsize=(7.2, 6.5))

    fold_aucs: list[float] = []
    for i, k in enumerate(folds):
        sub = df[df["fold"] == k]
        y = sub["label"].to_numpy()
        p = sub["prob_DM1"].to_numpy()
        if len(np.unique(y)) < 2:
            auc = float("nan")
            fold_aucs.append(auc)
            continue
        fpr, tpr, _ = roc_curve(y, p)
        auc = roc_auc_score(y, p)
        fold_aucs.append(auc)
        ax.plot(
            fpr,
            tpr,
            color=FOLD_COLORS[i % len(FOLD_COLORS)],
            linewidth=1.6,
            alpha=0.9,
            label=f"Fold {k} AUC={auc:.2f}",
        )

    # Pooled (cross-fold) ROC
    y_all = df["label"].to_numpy()
    p_all = df["prob_DM1"].to_numpy()
    fpr_p, tpr_p, _ = roc_curve(y_all, p_all)
    auc_p = roc_auc_score(y_all, p_all)
    ax.plot(fpr_p, tpr_p, color=POOLED_COLOR, linewidth=3.0, label=f"Pooled AUC={auc_p:.2f}")

    # Chance diagonal
    ax.plot([0, 1], [0, 1], color=BASELINE_COLOR, linestyle="--", linewidth=1.0, label="Chance")

    # Closure baseline reference (horizontal)
    ax.axhline(
        CLOSURE_BASELINE,
        color="#882255",
        linestyle=":",
        linewidth=1.0,
        label=f"Closure ResNet50 baseline ~{CLOSURE_BASELINE:.2f}",
    )

    mean_auc = float(np.nanmean(fold_aucs))
    sd_auc = float(np.nanstd(fold_aucs, ddof=0))

    ax.set_xlim(0, 1.0)
    ax.set_ylim(0, 1.0)
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.set_title(
        "Phase 2 — TCGA-THCA CLAM (UNI) DM1 vs not-DM, 5-fold CV (QUICK)",
        fontsize=12,
        pad=8,
    )

    ax.text(
        0.98,
        0.04,
        f"Mean AUC = {mean_auc:.2f} +/- {sd_auc:.2f}\nPooled AUC = {auc_p:.2f}",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=11,
        fontweight="bold",
        bbox=dict(facecolor="white", edgecolor="black", boxstyle="round,pad=0.3"),
    )

    ax.legend(loc="lower right", fontsize=9, frameon=True)
    ax.grid(linestyle=":", linewidth=0.4, alpha=0.6)
    ax.set_aspect("equal")
    fig.tight_layout()

    out_png = FIG_DIR / "fig4_phase2_roc.png"
    out_pdf = FIG_DIR / "fig4_phase2_roc.pdf"
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)

    print(f"wrote {out_png}")
    print(f"wrote {out_pdf}")
    print(f"fold AUCs = {[f'{a:.2f}' for a in fold_aucs]}")
    print(f"mean +/- sd = {mean_auc:.3f} +/- {sd_auc:.3f}; pooled = {auc_p:.3f}")


if __name__ == "__main__":
    main()
