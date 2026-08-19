#!/usr/bin/env python3
"""Figure 33 — Final integrated evidence scoreboard.

A single dashboard figure summarising the top statistics across all analyses:
6-cohort meta-d, Cox HR continuous, TERT × zone interaction, Mun proteome OR,
Korean dose-response, methylation correlation, AJCC stage trend.
"""
from __future__ import annotations
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
from _nature_style import setup_rc, panel_letter, panel_title, IVORY, INK, MUTED, BLUE, RED

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

setup_rc()

ROOT = Path(__file__).resolve().parent.parent
OUT_PNG = ROOT / "results" / "figures" / "figure33_final_scoreboard.png"
OUT_PDF = ROOT / "results" / "figures" / "figure33_final_scoreboard.pdf"


def main():
    scores = [
        # (label, stat, value, color)
        ("Cross-cohort meta forest (n = 6 cohorts)", "pooled Cohen d", "−2.48  [−2.64, −2.33]", "#1f4d7a"),
        ("Continuous panel z PFI Cox  (TCGA, n = 477)", "HR per unit", "0.60  [0.42, 0.87],  p = 0.007", "#386a99"),
        ("Mun 2025 proteome ATC vs PTC dark-matter",   "Fisher OR",   "8.54,  p = 2.7 × 10⁻¹⁵", "#5a4470"),
        ("Lee 2024 Korean histology dose-response",    "KW p-value",  "p = 8.7 × 10⁻⁵⁵  (n = 632)", "#4a2c6f"),
        ("Mutation group × panel z (TCGA)",            "KW p-value",  "p = 9.5 × 10⁻³⁷  (n = 500)", "#7c4d2a"),
        ("AJCC stage × panel z (TCGA, I → IV)",        "KW p-value",  "p = 2.4 × 10⁻⁴", "#a04e48"),
        ("dark-matter · TERT+ PFI rate",               "rate vs ref", "27.8% vs 7.0%  (4×, log-rank p = 0.037)", "#7c322e"),
        ("Lu 2023 sc-RNA per-cell distribution",       "KW H-stat",   "H = 10,471.5,  p < 10⁻³⁰⁰  (cells = 14,624)", "#33445e"),
        ("Top-quartile purity zone separation",        "KW p-value",  "p = 1.0 × 10⁻²¹  (n = 131)", "#456478"),
        ("8-gene Spearman ρ (TCGA cross-correlations)","mean ρ", "all pairs positive, max 0.86 (TG ↔ TPO)", "#386a99"),
        ("Cherry-pick null permutation (n = 1,000)",   "empirical p", "p = 0.014  (matched-variance null)", "#2a5468"),
        ("Methylation × RNA per-gene correlation",     "Δβ vs Δrna",  "8/8 panel genes show β↑ + RNA↓", "#5a4470"),
    ]

    fig, ax = plt.subplots(figsize=(15, 9.0), facecolor=IVORY)
    fig.subplots_adjust(left=0.025, right=0.975, top=0.92, bottom=0.04)
    ax.set_xlim(0, 100); ax.set_ylim(0, 100)
    ax.axis("off")

    # Header
    ax.text(50, 96, "Final evidence scoreboard — 8-gene RAI differentiation panel",
            ha="center", fontsize=18, fontweight="bold", color=INK)
    ax.text(50, 91.5, "All headline statistics across discovery + validation + mechanistic + clinical anchors",
            ha="center", fontsize=11.5, color=MUTED, style="italic")

    n_rows = len(scores)
    row_h = 6.5
    top_y = 87
    for i, (label, stat, value, color) in enumerate(scores):
        y = top_y - (i+1) * row_h
        # Row background card
        ax.add_patch(FancyBboxPatch((3, y-0.7), 94, row_h-1.0,
                                     boxstyle="round,pad=0.2,rounding_size=0.8",
                                     facecolor="white", edgecolor="#e3e2dc", lw=0.8, zorder=1))
        # Colored side bar
        ax.add_patch(FancyBboxPatch((3, y-0.7), 1.6, row_h-1.0,
                                     boxstyle="round,pad=0,rounding_size=0.4",
                                     facecolor=color, edgecolor="none", zorder=2))
        # Label
        ax.text(6.5, y + 1.8, label, ha="left", va="center", fontsize=10.5, color=INK, fontweight="bold", zorder=3)
        # Stat label
        ax.text(6.5, y + 0.0, stat, ha="left", va="center", fontsize=9.0, color=MUTED, style="italic", zorder=3)
        # Value
        ax.text(95.5, y + 0.9, value, ha="right", va="center", fontsize=11, color=color, fontweight="bold", zorder=3)

    fig.savefig(OUT_PNG, dpi=220, facecolor=IVORY)
    fig.savefig(OUT_PDF, facecolor=IVORY)
    plt.close(fig)
    print(f"wrote {OUT_PNG}\nwrote {OUT_PDF}")


if __name__ == "__main__":
    main()
