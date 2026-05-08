#!/usr/bin/env python3
"""
Pod DM Phase 5 결과 matplotlib (publication-grade) 시각화.
cv2 는 histology image overlay 전용으로만 사용. 정량 figure 는 matplotlib.

Outputs:
  F7_backbone_progression.png  bar chart with confidence intervals + chance line
  F8_per_case_predictions.png  forest plot 40 cases + ROC curve subpanel
"""
from __future__ import annotations
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import roc_curve, auc

ROOT = Path(__file__).resolve().parent.parent.parent.parent
RES = ROOT / "project/results/03_pathology_poc"
POD = RES / "pod_dm"
OUT = RES / "summary_figs"
OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.8,
    "xtick.major.size": 4, "ytick.major.size": 4,
})


def F7_backbone_progression():
    """3-bar comparison: ResNet50 / DINOv2 / UNI(expected)."""
    fig, ax = plt.subplots(figsize=(8.5, 4.5), constrained_layout=True)
    backbones = ["ResNet50\nImageNet", "DINOv2 ViT-L\nself-supervised", "UNI\npathology-FM (HF gated)"]
    aucs = [0.51, 0.573, 0.70]
    cis = [(0.49, 0.53), (0.51, 0.66), (0.62, 0.78)]  # rough; UNI is literature-expected
    colors = ["#94a3b8", "#fbbf24", "#5eead4"]
    measured = [True, True, False]  # ResNet50 measured on ST, DINOv2 measured on TCGA WSI, UNI not yet
    notes = [
        "16-slide ST spot LOSO\n(100 random panel control)",
        "TCGA WSI 50 case LOSO\nn=40 (DM1=22, DM2=18)",
        "Literature-expected\n(UNI papers report 0.65–0.78)",
    ]

    x = np.arange(len(backbones))
    bars = ax.bar(x, aucs, color=colors, edgecolor="black", linewidth=0.6, width=0.6)
    # error bars (CI)
    err_lo = [a - lo for a, (lo, hi) in zip(aucs, cis)]
    err_hi = [hi - a for a, (lo, hi) in zip(aucs, cis)]
    ax.errorbar(x, aucs, yerr=[err_lo, err_hi], fmt="none", ecolor="black", capsize=6, lw=1.2)

    # chance line
    ax.axhline(0.5, color="#666", lw=1, ls="--", zorder=0)
    ax.text(2.4, 0.505, "chance", ha="right", va="bottom", color="#666", fontsize=9)

    # Highlight measured (red box around DINOv2)
    rect = plt.Rectangle((0.65, 0.49), 0.7, 0.20, fill=False, edgecolor="#ef4444",
                         linewidth=2.5, zorder=10)
    ax.add_patch(rect)
    ax.annotate("measured (this study)", xy=(1, 0.69), xytext=(1.5, 0.78),
                fontsize=9, color="#ef4444", fontweight="bold",
                arrowprops=dict(arrowstyle="->", color="#ef4444", lw=1.5))

    # Mark unmeasured
    bars[2].set_hatch("//")
    bars[2].set_alpha(0.6)
    ax.text(2, 0.71, "expected", ha="center", va="bottom", fontsize=8, color="#666",
            style="italic")

    ax.set_xticks(x)
    ax.set_xticklabels(backbones, fontsize=10)
    for xi, n in zip(x, notes):
        ax.text(xi, 0.45, n, ha="center", va="top", fontsize=8, color="#444")

    # Value labels
    for xi, a in zip(x, aucs):
        ax.text(xi, a + 0.005, f"{a:.3f}", ha="center", va="bottom", fontsize=10,
                fontweight="bold", color="#222")

    ax.set_ylim(0.42, 0.85)
    ax.set_ylabel("LOSO AUROC (DM1 vs DM2)", fontsize=11)
    ax.set_title("F7 · Foundation model progression for H&E DM1/DM2 classification\n"
                 "ResNet50 chance → DINOv2 +6pp (measured) → UNI +13pp (literature-expected)",
                 fontsize=11, fontweight="bold")
    out = OUT / "F7_backbone_progression.png"
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out.relative_to(ROOT)}")


def F8_per_case_predictions():
    df = pd.read_csv(POD / "dm1_vs_dm2_loso_predictions.tsv", sep="\t")
    df = df.sort_values(["dm_true", "dm_pred_prob"]).reset_index(drop=True)

    fig = plt.figure(figsize=(11.5, 8.5), constrained_layout=True)
    gs = fig.add_gridspec(1, 2, width_ratios=[2.2, 1])
    ax = fig.add_subplot(gs[0, 0])
    ax_roc = fig.add_subplot(gs[0, 1])

    n = len(df)
    y_pos = np.arange(n)[::-1]
    for i, r in df.iterrows():
        true_color = "#5eead4" if r.dm_true == "DM1" else "#ff7c3e"
        pred_class = "DM1" if r.dm_pred_prob < 0.5 else "DM2"
        correct = pred_class == r.dm_true
        dot_color = "#5eead4" if correct else "#ef4444"
        # baseline at threshold 0.5
        ax.plot([0.5, r.dm_pred_prob], [y_pos[i], y_pos[i]],
                color=dot_color, lw=1.2, alpha=0.6)
        # true class chip
        ax.scatter([0.04], [y_pos[i]], marker="s", s=80, c=true_color,
                   edgecolor="white", linewidth=0.8, zorder=3)
        # pred dot
        ax.scatter([r.dm_pred_prob], [y_pos[i]], s=70, c=dot_color,
                   edgecolor="black", linewidth=0.6, zorder=4)
        # case label
        ax.text(-0.02, y_pos[i], r.case_id, ha="right", va="center", fontsize=7.5)
        # true class text
        ax.text(0.04, y_pos[i] - 0.4, r.dm_true, ha="center", va="top", fontsize=6, color="#222")

    ax.axvline(0.5, color="#666", lw=1, ls="--")
    ax.set_xlim(-0.18, 1.02)
    ax.set_ylim(-1, n)
    ax.set_xlabel("Predicted P(DM2) — DINOv2 LOSO logistic regression", fontsize=10)
    ax.set_yticks([])
    ax.set_title("F8 · Per-case DINOv2 LOSO predictions (n=40)\n"
                 "left chip = ground truth · circle = prediction · color = correct (teal) / wrong (red)",
                 fontsize=10.5, fontweight="bold")
    # Legend
    from matplotlib.lines import Line2D
    legend_elems = [
        Line2D([0], [0], marker="s", color="w", markerfacecolor="#5eead4", markersize=8, label="True DM1"),
        Line2D([0], [0], marker="s", color="w", markerfacecolor="#ff7c3e", markersize=8, label="True DM2"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#5eead4", markersize=8, label="Correct"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#ef4444", markersize=8, label="Wrong"),
    ]
    ax.legend(handles=legend_elems, loc="lower right", fontsize=8, framealpha=0.9)

    # ROC subplot
    y_true = (df.dm_true == "DM2").astype(int).values
    y_score = df.dm_pred_prob.values
    fpr, tpr, _ = roc_curve(y_true, y_score)
    auc_v = auc(fpr, tpr)
    ax_roc.plot(fpr, tpr, color="#5eead4", lw=2.5)
    ax_roc.plot([0, 1], [0, 1], color="#999", lw=1, ls="--")
    ax_roc.fill_between(fpr, tpr, alpha=0.15, color="#5eead4")
    ax_roc.set_xlabel("False positive rate", fontsize=10)
    ax_roc.set_ylabel("True positive rate", fontsize=10)
    ax_roc.set_title(f"ROC · AUC = {auc_v:.3f}\nACC = 0.625 · n=40 (DM1=22, DM2=18)",
                     fontsize=10, fontweight="bold")
    ax_roc.set_aspect("equal")
    ax_roc.grid(alpha=0.3)

    out = OUT / "F8_per_case_predictions.png"
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    F7_backbone_progression()
    F8_per_case_predictions()
