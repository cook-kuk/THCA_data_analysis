#!/usr/bin/env python3
"""
Paper 2 Figure 7: Method comparison table figure.

4 architectures × 6 properties, rendered as a table-as-figure
(matplotlib only, no LaTeX).
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

OUT_DIR = Path(
    "/home/seungho/personal/THCA_data_analysis/project/results/"
    "p2_image_dm1_v2_foundation_clam_2026_05_07/figures"
)
OUT_DIR.mkdir(parents=True, exist_ok=True)

HEADER_BG = "#2b3e62"
HEADER_FG = "#ffffff"

ROW_COLORS = {
    "ResNet50":  "#d6d6d6",   # gray (closure baseline)
    "ViT-L":     "#fff2bf",   # pale yellow (fallback)
    "UNI":       "#9adb9e",   # green highlight (FINAL)
    "CONCH":     "#e6e6e6",   # light gray (future)
}

COLS = [
    "Architecture",
    "Training corpus",
    "Domain",
    "AUC (Phase 2 TCGA-THCA)",
    "Foundation model?",
    "Status",
]

ROWS = [
    ["ResNet50",
     "ImageNet-1k (1.3M nat. images)",
     "Natural images",
     "0.55",
     "No",
     "Closure baseline (2026-05-04)"],
    ["ViT-L",
     "ImageNet-21k (~14M nat. images)",
     "Natural images",
     "0.746",
     "No",
     "Phase 2 fallback (N=59)"],
    ["UNI",
     "Mass-100K (100K WSIs, 100M tiles)",
     "H&E pathology",
     "0.874",
     "Yes",
     "FINAL (N=54, 5-fold mean 0.951 ± 0.081)"],
    ["CONCH",
     "Quilt-1M (1.1M image-text pairs)",
     "H&E pathology + text",
     "future",
     "Yes (vision-language)",
     "Future (alternative target, gated)"],
]


def main():
    n_cols = len(COLS)
    n_rows = len(ROWS)

    # geometry — widened to fit long strings; long fields use \n wrapping
    col_widths = [1.30, 2.55, 1.55, 1.55, 1.45, 2.85]
    total_w = sum(col_widths)
    row_h = 1.20
    header_h = 0.85
    total_h = header_h + n_rows * row_h

    fig_w = total_w + 0.6
    fig_h = total_h + 1.2

    # pre-wrap long cell strings to two lines so text fits column width
    WRAPPED = [
        ["ResNet50",
         "ImageNet-1k\n(1.3M nat. images)",
         "Natural images",
         "0.55",
         "No",
         "Closure baseline\n(2026-05-04)"],
        ["ViT-L",
         "ImageNet-21k\n(~14M nat. images)",
         "Natural images",
         "0.746",
         "No",
         "Phase 2 fallback\n(N=59)"],
        ["UNI",
         "Mass-100K\n(100K WSIs, 100M tiles)",
         "H&E pathology",
         "0.874",
         "Yes",
         "FINAL (N=54)\n5-fold mean 0.951 ± 0.081"],
        ["CONCH",
         "Quilt-1M\n(1.1M image-text pairs)",
         "H&E pathology + text",
         "future",
         "Yes\n(vision-language)",
         "Future\n(alternative target, gated)"],
    ]

    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=300)
    ax.set_xlim(0, total_w)
    ax.set_ylim(0, total_h + 0.7)
    ax.axis("off")

    # column x-positions
    x_starts = [0.0]
    for w in col_widths[:-1]:
        x_starts.append(x_starts[-1] + w)

    # ---------- header
    y_header = total_h - header_h
    ax.add_patch(Rectangle(
        (0, y_header), total_w, header_h,
        facecolor=HEADER_BG, edgecolor="black", lw=1.0, zorder=1))
    HEADERS_WRAPPED = [
        "Architecture",
        "Training corpus",
        "Domain",
        "AUC\n(Phase 2 TCGA-THCA)",
        "Foundation\nmodel?",
        "Status",
    ]
    for ci, col in enumerate(HEADERS_WRAPPED):
        ax.text(
            x_starts[ci] + col_widths[ci] / 2,
            y_header + header_h / 2,
            col, ha="center", va="center",
            color=HEADER_FG, weight="bold", fontsize=10.0, zorder=2,
        )

    # ---------- data rows
    for ri, row in enumerate(WRAPPED):
        y = y_header - (ri + 1) * row_h
        bg = ROW_COLORS[row[0]]
        ax.add_patch(Rectangle(
            (0, y), total_w, row_h,
            facecolor=bg, edgecolor="black", lw=0.8, zorder=1))

        # vertical separators
        for xs in x_starts[1:]:
            ax.plot([xs, xs], [y, y + row_h],
                    color="black", lw=0.5, zorder=2)

        for ci, val in enumerate(row):
            cx = x_starts[ci] + col_widths[ci] / 2
            cy = y + row_h / 2

            weight = "normal"
            color = "#111111"
            fontsize = 9.0

            if ci == 0:
                weight = "bold"
                fontsize = 10.5

            # AUC column highlights
            if ci == 3:
                if row[0] == "UNI":
                    weight = "bold"
                    color = "#0a5d18"
                    fontsize = 14.0
                elif row[0] == "ViT-L":
                    color = "#7a4a00"
                    fontsize = 11.5
                    weight = "bold"
                elif row[0] == "ResNet50":
                    color = "#5a5a5a"
                    fontsize = 11.5
                elif val in ("pending", "future"):
                    color = "#7a7a7a"
                    fontsize = 10.0
                    weight = "normal"

            ax.text(cx, cy, val, ha="center", va="center",
                    fontsize=fontsize, color=color, weight=weight, zorder=3)

    # outer border
    ax.add_patch(Rectangle(
        (0, y_header - n_rows * row_h), total_w, total_h,
        fill=False, edgecolor="black", lw=1.4, zorder=4))

    # ---------- title + footnote
    fig.suptitle(
        "Figure 7.  Method comparison — backbones evaluated on H&E → DM1",
        fontsize=12.5, weight="bold", y=0.985,
    )
    ax.text(
        total_w / 2, total_h + 0.30,
        "Bold green AUC = UNI Mass-100K final pooled value (Phase 2, N=54, "
        "5-fold mean 0.951 ± 0.081, Δ +0.128 over ViT-L). "
        "CONCH gated on follow-up dispatch.",
        ha="center", va="center", fontsize=8.8,
        style="italic", color="#444",
    )

    plt.tight_layout(rect=(0, 0, 1, 0.96))

    png = OUT_DIR / "fig7_method_comparison.png"
    pdf = OUT_DIR / "fig7_method_comparison.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"WROTE {png}")
    print(f"WROTE {pdf}")


if __name__ == "__main__":
    main()
