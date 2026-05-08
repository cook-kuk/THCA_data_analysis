#!/usr/bin/env python
"""Fig 6 — Decision-curve analysis (left) + confusion matrices (right).

Left panel: decision-curve net benefit vs threshold for the pooled per-slide
CLAM probabilities, plus the standard "treat all" and "treat none" reference
strategies.

Right panel: 2x2 confusion matrices at threshold = 0.5 and at the
Youden-optimal threshold (argmax of TPR - FPR on the pooled ROC).

Inputs:
  - phase2_tcga_clam_QUICK/clam_per_slide_predictions.tsv

Outputs:
  - figures/fig6_decision_curve_confusion.png
  - figures/fig6_decision_curve_confusion.pdf
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import confusion_matrix, roc_curve

ROOT = Path(
    "/home/seungho/personal/THCA_data_analysis/project/results/"
    "p2_image_dm1_v2_foundation_clam_2026_05_07"
)
PHASE2 = ROOT / "phase2_tcga_clam_QUICK"
FIG_DIR = ROOT / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Wong color-blind palette
COL_MODEL = "#0072B2"
COL_ALL = "#E69F00"
COL_NONE = "#999999"


def net_benefit(y: np.ndarray, p: np.ndarray, thresholds: np.ndarray) -> np.ndarray:
    """Standard Vickers decision-curve net benefit at each threshold."""
    n = len(y)
    nb = np.empty_like(thresholds, dtype=float)
    for i, t in enumerate(thresholds):
        pred = (p >= t).astype(int)
        tp = int(((pred == 1) & (y == 1)).sum())
        fp = int(((pred == 1) & (y == 0)).sum())
        # Avoid div-by-zero at t -> 1
        denom = max(1.0 - t, 1e-9)
        nb[i] = (tp / n) - (fp / n) * (t / denom)
    return nb


def net_benefit_all(y: np.ndarray, thresholds: np.ndarray) -> np.ndarray:
    n = len(y)
    pos = float((y == 1).sum())
    neg = float((y == 0).sum())
    nb = np.empty_like(thresholds, dtype=float)
    for i, t in enumerate(thresholds):
        denom = max(1.0 - t, 1e-9)
        nb[i] = (pos / n) - (neg / n) * (t / denom)
    return nb


def plot_cm(ax, cm: np.ndarray, title: str) -> None:
    im = ax.imshow(cm, cmap="Blues", vmin=0)
    ax.set_xticks([0, 1])
    ax.set_yticks([0, 1])
    ax.set_xticklabels(["pred not-DM", "pred DM1"])
    ax.set_yticklabels(["true not-DM", "true DM1"])
    ax.set_title(title, fontsize=11)
    cmax = cm.max() if cm.size else 1
    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            v = int(cm[i, j])
            color = "white" if v > 0.6 * cmax else "black"
            ax.text(j, i, str(v), ha="center", va="center", fontsize=12, color=color)


def main() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 11,
            "savefig.dpi": 300,
        }
    )

    df = pd.read_csv(PHASE2 / "clam_per_slide_predictions.tsv", sep="\t")
    y = df["label"].to_numpy().astype(int)
    p = df["prob_DM1"].to_numpy().astype(float)
    n = len(y)

    # Youden-optimal threshold from pooled ROC
    fpr, tpr, thr = roc_curve(y, p)
    j = tpr - fpr
    j_idx = int(np.argmax(j))
    t_youden = float(thr[j_idx])
    if not np.isfinite(t_youden):
        # roc_curve may emit a leading +inf; fall back to next best
        finite = np.isfinite(thr)
        j_idx = int(np.argmax(np.where(finite, j, -np.inf)))
        t_youden = float(thr[j_idx])

    # Decision-curve thresholds
    ts = np.linspace(0.01, 0.99, 99)
    nb_model = net_benefit(y, p, ts)
    nb_all = net_benefit_all(y, ts)
    nb_none = np.zeros_like(ts)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), gridspec_kw={"width_ratios": [1.25, 1.0]})

    # Left: decision curve
    ax = axes[0]
    ax.plot(ts, nb_model, color=COL_MODEL, linewidth=2.2, label="CLAM (UNI) pooled")
    ax.plot(ts, nb_all, color=COL_ALL, linewidth=1.5, linestyle="--", label="Treat all")
    ax.plot(ts, nb_none, color=COL_NONE, linewidth=1.5, linestyle=":", label="Treat none")
    ax.axvline(0.5, color="black", linewidth=0.5, alpha=0.5)
    ax.axvline(t_youden, color=COL_MODEL, linewidth=0.8, linestyle="--", alpha=0.7)
    ax.set_xlim(0, 1)
    ymin = float(min(nb_model.min(), nb_all.min(), 0)) - 0.05
    ymax = float(max(nb_model.max(), nb_all.max(), 0)) + 0.05
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel("Threshold probability")
    ax.set_ylabel("Net benefit")
    ax.set_title("Decision-curve analysis (pooled, n={n})".format(n=n), fontsize=12)
    ax.legend(loc="upper right", fontsize=9, frameon=True)
    ax.grid(linestyle=":", linewidth=0.4, alpha=0.6)
    ax.text(
        t_youden,
        ymax - 0.02,
        f"Youden t={t_youden:.2f}",
        ha="left",
        va="top",
        fontsize=9,
        color=COL_MODEL,
    )

    # Right: stacked CMs
    ax_cm = axes[1]
    ax_cm.axis("off")
    # Two sub-axes inside the right panel
    sub1 = fig.add_axes([0.585, 0.55, 0.18, 0.36])
    sub2 = fig.add_axes([0.795, 0.55, 0.18, 0.36])
    sub3 = fig.add_axes([0.585, 0.10, 0.18, 0.36])
    sub4 = fig.add_axes([0.795, 0.10, 0.18, 0.36])

    pred_05 = (p >= 0.5).astype(int)
    pred_y = (p >= t_youden).astype(int)
    cm_05 = confusion_matrix(y, pred_05, labels=[0, 1])
    cm_y = confusion_matrix(y, pred_y, labels=[0, 1])

    def metrics(cm):
        tn, fp, fn, tp = cm[0, 0], cm[0, 1], cm[1, 0], cm[1, 1]
        acc = (tp + tn) / max(cm.sum(), 1)
        sens = tp / max(tp + fn, 1)
        spec = tn / max(tn + fp, 1)
        return acc, sens, spec

    acc05, sens05, spec05 = metrics(cm_05)
    accY, sensY, specY = metrics(cm_y)

    plot_cm(sub1, cm_05, f"t = 0.50\nacc={acc05:.2f}  sens={sens05:.2f}  spec={spec05:.2f}")
    sub2.axis("off")
    sub2.text(
        0.0,
        0.5,
        (
            "Threshold = 0.50\n"
            f"  TP = {cm_05[1,1]}\n  FP = {cm_05[0,1]}\n"
            f"  FN = {cm_05[1,0]}\n  TN = {cm_05[0,0]}\n"
            f"  N  = {cm_05.sum()}"
        ),
        va="center",
        ha="left",
        fontsize=10,
        family="monospace",
    )
    plot_cm(sub3, cm_y, f"t = {t_youden:.2f} (Youden)\nacc={accY:.2f}  sens={sensY:.2f}  spec={specY:.2f}")
    sub4.axis("off")
    sub4.text(
        0.0,
        0.5,
        (
            f"Threshold = {t_youden:.2f}\n"
            f"  TP = {cm_y[1,1]}\n  FP = {cm_y[0,1]}\n"
            f"  FN = {cm_y[1,0]}\n  TN = {cm_y[0,0]}\n"
            f"  N  = {cm_y.sum()}"
        ),
        va="center",
        ha="left",
        fontsize=10,
        family="monospace",
    )

    fig.suptitle(
        "Phase 2 — Decision-curve analysis and confusion matrices (TCGA-THCA CLAM QUICK)",
        fontsize=12,
        y=0.995,
    )

    out_png = FIG_DIR / "fig6_decision_curve_confusion.png"
    out_pdf = FIG_DIR / "fig6_decision_curve_confusion.pdf"
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)

    print(f"wrote {out_png}")
    print(f"wrote {out_pdf}")
    print(f"Youden t = {t_youden:.3f}")
    print(f"t=0.50  cm = {cm_05.tolist()}  acc={acc05:.3f} sens={sens05:.3f} spec={spec05:.3f}")
    print(f"t={t_youden:.2f}  cm = {cm_y.tolist()}  acc={accY:.3f} sens={sensY:.3f} spec={specY:.3f}")


if __name__ == "__main__":
    main()
