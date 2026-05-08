#!/usr/bin/env python
"""Fig 3 — Phase 1 correlation forest plot.

Horizontal bar plot of per-slide PC1 rho_DM1 values, sorted ascending.
Color: blue if rho > 0, orange if rho < 0 (color-blind friendly).
Bars with |rho| > 0.3 are saturated; others are translucent.
Vertical dashed lines mark the +/- 0.3 kill-switch thresholds.
Annotations: max-|rho| slide name and the "n_pass / 16 PASS" tally.

Inputs:
  - phase1_gse250521/uni_dm1_correlation_per_slide.tsv

Outputs:
  - figures/fig3_phase1_correlation_forest.png
  - figures/fig3_phase1_correlation_forest.pdf
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(
    "/home/seungho/personal/THCA_data_analysis/project/results/"
    "p2_image_dm1_v2_foundation_clam_2026_05_07"
)
PHASE1 = ROOT / "phase1_gse250521"
FIG_DIR = ROOT / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

THRESH = 0.3
COLOR_POS = "#0072B2"  # Wong blue
COLOR_NEG = "#E69F00"  # Wong orange


def main() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 11,
            "savefig.dpi": 300,
        }
    )

    df = pd.read_csv(PHASE1 / "uni_dm1_correlation_per_slide.tsv", sep="\t")
    pc1 = df[df["pc"] == 1].copy()
    pc1 = pc1.sort_values("rho_DM1", ascending=True).reset_index(drop=True)

    n = len(pc1)
    n_pass = int((pc1["rho_DM1"].abs() > THRESH).sum())

    short_labels = [
        s.split("_", 1)[1] if "_" in s else s[:10] for s in pc1["slide"]
    ]

    fig, ax = plt.subplots(figsize=(8.5, 6.5))

    y = np.arange(n)
    rho = pc1["rho_DM1"].to_numpy()
    passed = np.abs(rho) > THRESH

    colors = [COLOR_POS if r >= 0 else COLOR_NEG for r in rho]
    alphas = [1.0 if p else 0.35 for p in passed]

    for yi, ri, ci, ai in zip(y, rho, colors, alphas):
        ax.barh(yi, ri, color=ci, alpha=ai, edgecolor="black", linewidth=0.4, height=0.7)

    ax.axvline(0, color="black", linewidth=0.6)
    ax.axvline(THRESH, color="#117733", linestyle="--", linewidth=1.0, label=f"+/- {THRESH} kill-switch")
    ax.axvline(-THRESH, color="#117733", linestyle="--", linewidth=1.0)

    ax.set_yticks(y)
    ax.set_yticklabels(short_labels, fontsize=10)
    ax.set_xlim(-1.0, 1.0)
    ax.set_xlabel("PC1 Spearman rho (UNI embedding vs DM1 score)", fontsize=11)
    ax.set_title(
        "Phase 1 forest plot — per-slide UNI x DM1 correlation (GSE250521, n=16)",
        fontsize=12,
        pad=8,
    )

    # Max |rho| annotation
    imax = int(np.argmax(np.abs(rho)))
    rmax = float(rho[imax])
    smax = short_labels[imax]
    ax.annotate(
        f"max |rho| = {rmax:+.2f}  ({smax})",
        xy=(rmax, imax),
        xytext=(0.55 if rmax < 0 else -0.95, imax),
        ha="left",
        va="center",
        fontsize=10,
        arrowprops=dict(arrowstyle="->", color="black", lw=0.6),
    )

    ax.text(
        0.98,
        0.04,
        f"{n_pass}/{n} PASS  (|rho| > {THRESH})",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=11,
        fontweight="bold",
        bbox=dict(facecolor="white", edgecolor="black", boxstyle="round,pad=0.3"),
    )

    # Manual legend for color encoding
    from matplotlib.patches import Patch
    handles = [
        Patch(facecolor=COLOR_POS, edgecolor="black", label="rho > 0"),
        Patch(facecolor=COLOR_NEG, edgecolor="black", label="rho < 0"),
        Patch(facecolor="white", edgecolor="#117733", linestyle="--", label=f"+/- {THRESH} threshold"),
    ]
    ax.legend(handles=handles, loc="lower left", fontsize=9, frameon=True)

    ax.grid(axis="x", linestyle=":", linewidth=0.4, alpha=0.6)
    fig.tight_layout()

    out_png = FIG_DIR / "fig3_phase1_correlation_forest.png"
    out_pdf = FIG_DIR / "fig3_phase1_correlation_forest.pdf"
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)

    print(f"wrote {out_png}")
    print(f"wrote {out_pdf}")
    print(f"PASS = {n_pass}/{n}; max |rho| = {rmax:+.3f} ({smax})")


if __name__ == "__main__":
    main()
