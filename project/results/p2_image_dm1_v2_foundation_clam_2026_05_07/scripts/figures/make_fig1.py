#!/usr/bin/env python3
"""
Paper 2 Figure 1: Clinical triage schematic.

Flow chart depicting how the H&E -> DM1 image classifier slots into
BRAF/RAS-negative PTC clinical triage, gating an 8-gene reflex panel
and downstream RAI-refractory / ICI-candidate / surveillance arms.

Pure matplotlib (Rectangle / FancyBboxPatch / FancyArrowPatch).
"""
from __future__ import annotations

import os
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Polygon

OUT_DIR = Path(
    "/home/seungho/personal/THCA_data_analysis/project/results/"
    "p2_image_dm1_v2_foundation_clam_2026_05_07/figures"
)
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------------ palette
GRAY = "#cfcfcf"
LIGHT_BLUE = "#bcd9ec"
PURPLE = "#c9b6e0"
YELLOW = "#fde88c"
GREEN = "#bfe3bf"
RED_GRAY = "#e5d0d0"
TEXT = "#1a1a1a"
EDGE = "#222222"


def rounded(ax, x, y, w, h, label, fc, fontsize=9, weight="normal"):
    box = FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.02,rounding_size=0.10",
        linewidth=1.2, edgecolor=EDGE, facecolor=fc, zorder=2,
    )
    ax.add_patch(box)
    ax.text(x + w / 2, y + h / 2, label, ha="center", va="center",
            fontsize=fontsize, color=TEXT, weight=weight, wrap=True, zorder=3)
    return (x, y, w, h)


def diamond(ax, cx, cy, w, h, label, fc, fontsize=9):
    pts = [(cx, cy + h / 2), (cx + w / 2, cy),
           (cx, cy - h / 2), (cx - w / 2, cy)]
    poly = Polygon(pts, closed=True, facecolor=fc, edgecolor=EDGE,
                   linewidth=1.2, zorder=2)
    ax.add_patch(poly)
    ax.text(cx, cy, label, ha="center", va="center",
            fontsize=fontsize, color=TEXT, weight="bold", zorder=3)
    return (cx, cy, w, h)


def arrow(ax, x0, y0, x1, y1, label=None, label_offset=(0, 0.12),
          color=EDGE, lw=1.4, style="-|>"):
    ar = FancyArrowPatch(
        (x0, y0), (x1, y1),
        arrowstyle=style, mutation_scale=14,
        color=color, linewidth=lw, zorder=1,
    )
    ax.add_patch(ar)
    if label:
        mx, my = (x0 + x1) / 2 + label_offset[0], (y0 + y1) / 2 + label_offset[1]
        ax.text(mx, my, label, ha="center", va="center",
                fontsize=8.5, color=color, weight="bold",
                bbox=dict(boxstyle="round,pad=0.18", fc="white",
                          ec="none", alpha=0.85))


def main():
    fig, ax = plt.subplots(figsize=(15.5, 7.5), dpi=300)
    ax.set_xlim(0, 16)
    ax.set_ylim(0, 8)
    ax.set_aspect("equal")
    ax.axis("off")

    # -------- Row 1 (input pipeline): nodes 1..3 + diamond
    n1 = rounded(ax, 0.20, 5.4, 2.5, 1.4,
                 "PTC patient diagnosed\n(BRAF / RAS-negative)",
                 GRAY, fontsize=9, weight="bold")
    n2 = rounded(ax, 3.10, 5.4, 2.5, 1.4,
                 "H&E slide\n($0, already exists)",
                 LIGHT_BLUE, fontsize=9, weight="bold")
    n3 = rounded(ax, 6.00, 5.4, 2.7, 1.4,
                 "UNI tile encoder\n+ CLAM gated MIL",
                 PURPLE, fontsize=9, weight="bold")
    d1 = diamond(ax, 10.45, 6.10, 2.6, 1.7,
                 "DM1 score\n> threshold?",
                 YELLOW, fontsize=9)

    # arrows row 1
    arrow(ax, 2.70, 6.10, 3.10, 6.10)
    arrow(ax, 5.60, 6.10, 6.00, 6.10)
    arrow(ax, 8.70, 6.10, 9.15, 6.10)

    # -------- YES branch -> reflex panel -> confirmed -> 3 outcomes
    arrow(ax, 11.75, 6.10, 12.55, 6.10, label="YES")

    rounded(ax, 12.55, 5.40, 2.95, 1.4,
            "Reflex 8-gene\nRNA panel",
            GREEN, fontsize=9, weight="bold")
    arrow(ax, 14.02, 5.40, 14.02, 4.90)
    rounded(ax, 12.55, 3.50, 2.95, 1.3,
            "DM1 confirmed",
            GREEN, fontsize=9.5, weight="bold")

    # 3 outcome boxes (right column, stacked)
    out_x = 12.55
    rounded(ax, out_x, 1.95, 3.10, 1.20,
            "RAI-refractory likelihood ↑\n→ TKI early consideration",
            "#f0d9b5", fontsize=8.5)
    rounded(ax, out_x, 0.55, 3.10, 1.20,
            "Hashimoto-like + IFN-γ\n→ ICI candidate flag",
            "#e7c4d8", fontsize=8.5)

    # surveillance — to the LEFT of confirmed (still under YES branch)
    rounded(ax, 8.55, 1.95, 3.30, 1.20,
            "6-mo follow-up\nsurveillance",
            "#cfe6ec", fontsize=8.5)

    # arrows from "DM1 confirmed" to 3 outcomes
    arrow(ax, 14.02, 3.50, 14.02, 3.18)  # to TKI
    arrow(ax, 14.02, 1.95, 14.02, 1.78)  # to ICI (chained)
    arrow(ax, 12.55, 4.05, 11.85, 3.05, color="#666666", lw=1.1)

    # -------- NO branch -> standard PTC management
    arrow(ax, 10.45, 5.25, 10.45, 4.55, label="NO", label_offset=(-0.32, 0))
    rounded(ax, 8.95, 3.30, 3.0, 1.2,
            "Standard PTC\nmanagement",
            RED_GRAY, fontsize=9.5, weight="bold")

    # -------- cost-savings annotation (below the NO branch)
    ax.text(
        10.45, 2.55,
        "Cost saved when score < threshold:\n"
        "≈ $800–1,500 RNA panel avoided / patient",
        ha="center", va="center", fontsize=8.5,
        color="#7a2b2b", style="italic", weight="bold",
        bbox=dict(boxstyle="round,pad=0.30", fc="#fff4e0",
                  ec="#7a2b2b", lw=1.0),
    )

    # -------- title + subtitle
    fig.suptitle(
        "Figure 1.  Clinical triage schematic — H&E → DM1 image classifier as "
        "a free-of-charge gate for the 8-gene reflex panel",
        fontsize=12.5, weight="bold", y=0.985,
    )
    ax.text(
        8.0, 7.55,
        "Paper 2 use case: BRAF/RAS-negative PTC enters; existing H&E slide is "
        "scored; only score-positive cases consume the 8-gene RNA panel; "
        "DM1-confirmed cases triage into RAI-refractory, ICI-candidate, "
        "or active-surveillance arms.",
        ha="center", va="center", fontsize=9, color="#333", style="italic",
    )

    # legend (bottom-left)
    legend_items = [
        ("Patient / clinical state", GRAY),
        ("Existing data ($0 input)", LIGHT_BLUE),
        ("Model / decision", PURPLE),
        ("Decision gate", YELLOW),
        ("Reflex panel / DM1 confirmed", GREEN),
        ("Standard arm", RED_GRAY),
    ]
    lx, ly = 0.30, 3.40
    ax.text(lx, ly + 0.35, "Legend", fontsize=9, weight="bold")
    for i, (lab, col) in enumerate(legend_items):
        ry = ly - i * 0.40
        rounded(ax, lx, ry - 0.15, 0.40, 0.25, "", col, fontsize=1)
        ax.text(lx + 0.55, ry - 0.02, lab, fontsize=8.2, va="center")

    plt.tight_layout(rect=(0, 0, 1, 0.96))

    png = OUT_DIR / "fig1_triage_schematic.png"
    pdf = OUT_DIR / "fig1_triage_schematic.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"WROTE {png}")
    print(f"WROTE {pdf}")


if __name__ == "__main__":
    main()
