#!/usr/bin/env python3
"""Figure 1 — A multi-gate model of radioiodine failure in thyroid cancer.

Nature Medicine / Nature Biomedical Engineering visual style. Four-panel composite:
  a) clinical unmet need: DTC → RAI → two outcomes (remission vs RAIR)
  b) multi-gate biology: blood → NIS → lineage TFs → organification → retention/killing
  c) three molecular states: avid (open gates) → gray zone (partial) → refractory (collapsed)
  d) evidence ladder: discovery → panel → label-anchored validation → stratification

Pure matplotlib + matplotlib.patches to keep the figure as a clean vector with no
external icon dependencies. Muted biomedical palette, thin lines, small sans labels.

Output: results/figures/figure1_multigate_rai_model.{png,pdf}
"""
from __future__ import annotations
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Circle, Wedge, Rectangle, Polygon
from matplotlib.patches import RegularPolygon
import matplotlib.patheffects as pe
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = ROOT / "results" / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Muted biomedical palette
BLUE   = "#3a6ea8"      # RAI-avid / preserved
BLUE_L = "#a9c2dc"
RED    = "#b04a44"      # RAI-refractory / silenced
RED_L  = "#d8a8a4"
GRAY_P = "#8e7aa5"      # gray zone / dark-matter
GRAY_L = "#c4baca"
BEIGE  = "#d8c8a8"      # colloid / accent
INK    = "#222222"
MUTED  = "#5a5a5a"
GRID   = "#e6e6e6"

def style_axis(ax):
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlim(0, 10); ax.set_ylim(0, 10)
    ax.set_aspect("auto")

def panel_label(ax, text):
    # panel letter goes ABOVE the axes via figure-level placement
    ax.text(-0.02, 1.13, text, transform=ax.transAxes, fontsize=15, fontweight="bold",
            family="DejaVu Sans", va="top")

def panel_title(ax, text):
    # the small header above the axes
    ax.text(0.5, 1.10, text, transform=ax.transAxes, fontsize=11, fontweight="bold",
            ha="center", color=INK, family="DejaVu Sans", va="top")

# ------------- Panel a: clinical unmet need -------------
def panel_a(ax):
    style_axis(ax)
    panel_label(ax, "a")
    panel_title(ax, "Clinical unmet need")
    # Thyroid icon (butterfly shape from two ellipses) + arrow → RAI capsule → split outcomes
    from matplotlib.patches import Ellipse
    # thyroid lobes
    for x_off in (-0.5, 0.5):
        e = Ellipse((1.3 + x_off, 7.4), 0.8, 1.2, angle=10 if x_off > 0 else -10,
                    facecolor=BEIGE, edgecolor=INK, lw=0.8, alpha=0.85)
        ax.add_patch(e)
    ax.text(1.3, 5.6, "Differentiated\nthyroid cancer\n(post-surgery)",
            ha="center", fontsize=8.5, color=INK)

    # Arrow → I-131 capsule
    ax.add_patch(FancyArrowPatch((2.5, 7.0), (4.0, 7.0), arrowstyle="->",
                                  mutation_scale=12, color=MUTED, lw=1.0))
    # Capsule
    cap = FancyBboxPatch((4.3, 6.55), 1.4, 0.9, boxstyle="round,pad=0,rounding_size=0.4",
                          facecolor=RED_L, edgecolor=RED, lw=0.8)
    ax.add_patch(cap)
    ax.text(5.0, 7.0, "I-131", fontsize=8.5, ha="center", va="center", color=RED, fontweight="bold")
    ax.text(5.0, 5.9, "Radioactive iodine therapy", ha="center", fontsize=8.2, color=MUTED)

    # Split arrows to outcomes
    ax.add_patch(FancyArrowPatch((5.7, 7.4), (8.0, 8.7), arrowstyle="->",
                                  mutation_scale=12, color=BLUE, lw=1.0))
    ax.add_patch(FancyArrowPatch((5.7, 6.6), (8.0, 5.3), arrowstyle="->",
                                  mutation_scale=12, color=RED, lw=1.0))

    # Outcome boxes
    box1 = FancyBboxPatch((7.9, 8.2), 2.0, 1.1, boxstyle="round,pad=0.05,rounding_size=0.15",
                          facecolor=BLUE_L, edgecolor=BLUE, lw=0.8)
    ax.add_patch(box1)
    ax.text(8.9, 8.95, "RAI-avid", ha="center", fontsize=9, color=BLUE, fontweight="bold")
    ax.text(8.9, 8.45, "remission", ha="center", fontsize=8.2, color=INK)

    box2 = FancyBboxPatch((7.9, 4.6), 2.0, 1.3, boxstyle="round,pad=0.05,rounding_size=0.15",
                          facecolor=RED_L, edgecolor=RED, lw=0.8)
    ax.add_patch(box2)
    ax.text(8.9, 5.55, "RAI-refractory", ha="center", fontsize=9, color=RED, fontweight="bold")
    ax.text(8.9, 5.05, "persistent /\nmetastatic", ha="center", fontsize=8.2, color=INK)

    # Bottom: key statistic
    ax.text(5.0, 3.0,
            "10-yr survival drops > 90% (RAI-responsive) → ~10–14% (distant RAIR-DTC)",
            ha="center", fontsize=8.5, color=MUTED, fontstyle="italic")
    ax.text(5.0, 2.0, "↓ a parsimonious molecular readout to triage RAI candidates",
            ha="center", fontsize=8.2, color=MUTED)
    # subtle bottom rule
    ax.plot([1.0, 9.0], [0.7, 0.7], color=GRID, lw=0.6)

# ------------- Panel b: multi-gate biology -------------
def panel_b(ax):
    style_axis(ax)
    panel_label(ax, "b")
    panel_title(ax, "Multi-gate biology of RAI therapy")

    # Vertical schematic — left: blood, middle: cell, right: colloid
    # bloodstream
    blood = Rectangle((0.6, 2.0), 1.4, 6.4, facecolor="#dceaf5", edgecolor=BLUE, lw=0.8)
    ax.add_patch(blood)
    ax.text(1.3, 8.6, "blood", ha="center", fontsize=8.5, color=BLUE, fontweight="bold")
    # iodide droplets in blood
    for y in (3.0, 4.5, 6.0, 7.5):
        ax.add_patch(Circle((1.3, y), 0.16, facecolor=RED, edgecolor="none", alpha=0.85))

    # thyrocyte cell (rounded rectangle)
    cell = FancyBboxPatch((2.8, 2.0), 4.4, 6.4, boxstyle="round,pad=0,rounding_size=0.4",
                           facecolor="#f3eee1", edgecolor=INK, lw=0.9)
    ax.add_patch(cell)
    ax.text(5.0, 8.6, "thyroid cancer cell", ha="center", fontsize=8.5, color=INK, fontweight="bold")
    # nucleus
    nuc = Circle((4.3, 5.0), 0.7, facecolor="#ebd8c8", edgecolor=MUTED, lw=0.7)
    ax.add_patch(nuc)
    ax.text(4.3, 5.0, "nucleus", ha="center", va="center", fontsize=7.5, color=MUTED)

    # colloid lumen
    coll = FancyBboxPatch((7.6, 3.5), 1.7, 3.4, boxstyle="round,pad=0,rounding_size=0.3",
                          facecolor=BEIGE, edgecolor="#a89270", lw=0.9)
    ax.add_patch(coll)
    ax.text(8.45, 6.6, "colloid", ha="center", fontsize=8.5, color="#7a6038", fontweight="bold")
    # TG storage stippling
    for y in np.linspace(4.0, 6.0, 4):
        for x in np.linspace(7.85, 9.0, 3):
            ax.add_patch(Circle((x + 0.04 * np.sin(y * 3), y), 0.07, facecolor="#a89270", alpha=0.6))

    # Gate 1: SLC5A5 (basolateral) — arrow blood → cell, label placed in the gap
    ax.add_patch(FancyArrowPatch((2.0, 7.6), (2.8, 7.6), arrowstyle="->",
                                  mutation_scale=11, color=RED, lw=1.2))
    g1 = FancyBboxPatch((2.05, 5.5), 0.75, 1.4, boxstyle="round,pad=0,rounding_size=0.12",
                        facecolor=RED_L, edgecolor=RED, lw=0.6)
    ax.add_patch(g1)
    ax.text(2.42, 6.55, "gate 1", fontsize=7.0, ha="center", color=RED, fontweight="bold")
    ax.text(2.42, 6.15, "SLC5A5\nNIS", fontsize=7.0, ha="center", color=INK)
    ax.text(2.42, 5.65, "uptake", fontsize=6.5, ha="center", color=MUTED, fontstyle="italic")

    # Gate 2: lineage TFs (in nucleus) — pulled out to top of cell
    g2 = FancyBboxPatch((3.0, 6.0), 1.7, 1.4, boxstyle="round,pad=0,rounding_size=0.15",
                        facecolor=BLUE_L, edgecolor=BLUE, lw=0.6)
    ax.add_patch(g2)
    ax.text(3.85, 7.1, "gate 2", fontsize=7.5, ha="center", color=BLUE)
    ax.text(3.85, 6.55, "PAX8 · NKX2-1\nFOXE1 · TSHR", fontsize=7.5, ha="center", color=INK)
    ax.add_patch(FancyArrowPatch((3.85, 6.0), (4.3, 5.7), arrowstyle="-",
                                  color=BLUE, lw=0.7))

    # Gate 3: organification / colloid (TG/TPO/DUOX) — apical membrane
    g3 = FancyBboxPatch((5.3, 3.0), 1.7, 1.4, boxstyle="round,pad=0,rounding_size=0.15",
                        facecolor=BLUE_L, edgecolor=BLUE, lw=0.6)
    ax.add_patch(g3)
    ax.text(6.15, 4.1, "gate 3", fontsize=7.5, ha="center", color=BLUE)
    ax.text(6.15, 3.55, "TG · TPO\nDUOX1/2", fontsize=7.5, ha="center", color=INK)
    ax.add_patch(FancyArrowPatch((7.0, 3.8), (7.55, 4.6), arrowstyle="->",
                                  mutation_scale=10, color=BLUE, lw=1.0))

    # Iodine path through cell: blood iodide → into cell → toward colloid
    for y in (3.4, 4.2, 5.7, 6.6):
        ax.add_patch(Circle((3.6, y), 0.10, facecolor=RED, edgecolor="none", alpha=0.85))
    ax.add_patch(FancyArrowPatch((3.6, 3.6), (5.6, 3.8), arrowstyle="-",
                                  color=RED, lw=0.7, linestyle="--"))

    # Gate 4: retention + radiation killing — label below colloid, rays inside
    # decay rays inside colloid
    for ang in range(0, 360, 45):
        rad = 0.28
        x0, y0 = 8.45, 5.0
        x1 = x0 + rad * np.cos(np.deg2rad(ang))
        y1 = y0 + rad * np.sin(np.deg2rad(ang))
        ax.plot([x0, x1], [y0, y1], color=RED, lw=0.7)
    # central red dot
    ax.add_patch(Circle((8.45, 5.0), 0.10, facecolor=RED, edgecolor="none"))
    # gate 4 label box below colloid
    g4 = FancyBboxPatch((7.55, 2.0), 1.8, 1.0, boxstyle="round,pad=0,rounding_size=0.12",
                        facecolor=RED_L, edgecolor=RED, lw=0.6)
    ax.add_patch(g4)
    ax.text(8.45, 2.7, "gate 4", fontsize=7.5, ha="center", color=RED, fontweight="bold")
    ax.text(8.45, 2.25, "retention &\nradiation kill", fontsize=7.0, ha="center", color=INK)

    # Caption strip
    ax.text(5.0, 1.0,
            "NIS alone is insufficient — coordinated thyroid differentiation across all four gates is required.",
            ha="center", fontsize=8.2, color=MUTED, fontstyle="italic")

# ------------- Panel c: three molecular states -------------
def panel_c(ax):
    style_axis(ax)
    panel_label(ax, "c")
    panel_title(ax, "Three molecular states along the differentiation-silencing axis")

    # Three vertically stacked horizontal panels
    state_specs = [
        ("RAI-avid differentiated state",        "all gates open · high 8-gene score",       BLUE_L,  BLUE,  +1.6),
        ("Molecular gray zone",                  "partial gate failure · intermediate score", GRAY_L,  GRAY_P, 0.0),
        ("RAI-refractory dedifferentiated state","collapsed gates · low 8-gene score",        RED_L,   RED,   -1.6),
    ]
    base_y = 8.0
    row_h  = 1.4
    for i, (title, sub, face, edge, score) in enumerate(state_specs):
        y = base_y - i * (row_h + 0.5)
        bg = FancyBboxPatch((0.6, y - row_h/2), 5.4, row_h, boxstyle="round,pad=0,rounding_size=0.15",
                            facecolor=face, edgecolor=edge, lw=0.9)
        ax.add_patch(bg)
        ax.text(0.85, y + 0.35, title, fontsize=9.5, color=edge, fontweight="bold")
        ax.text(0.85, y - 0.05, sub, fontsize=8.3, color=INK)

        # 4 gate indicators
        for g_i in range(4):
            cx = 1.0 + g_i * 0.85
            cy = y - 0.55
            if i == 0:   # all open
                state_color = BLUE
            elif i == 1: # partial
                state_color = GRAY_P if g_i < 2 else "#e0d7da"
            else:        # collapsed
                state_color = "#e0d7da" if g_i > 0 else RED
            ax.add_patch(Circle((cx, cy), 0.16, facecolor=state_color, edgecolor=edge, lw=0.6))
        ax.text(1.0 + 4*0.85, y - 0.55, "gates 1-4", fontsize=7.2, color=MUTED, va="center")

        # Score bar on the right
        bar_x = 6.5; bar_y_top = y + 0.45; bar_y_bot = y - 0.45
        ax.plot([bar_x, bar_x], [bar_y_bot, bar_y_top], color=MUTED, lw=0.6)
        ax.plot([bar_x - 0.1, bar_x + 0.1], [bar_y_top, bar_y_top], color=MUTED, lw=0.6)
        ax.plot([bar_x - 0.1, bar_x + 0.1], [bar_y_bot, bar_y_bot], color=MUTED, lw=0.6)
        norm = (score - (-1.6)) / (1.6 - (-1.6))
        marker_y = bar_y_bot + norm * (bar_y_top - bar_y_bot)
        ax.add_patch(Circle((bar_x, marker_y), 0.13, facecolor=edge, edgecolor="white", lw=0.6))
        ax.text(bar_x + 0.25, marker_y, f"z = {score:+.1f}", fontsize=8, va="center", color=edge)
        if i == 0:
            ax.text(bar_x + 1.2, bar_y_top + 0.1, "8-gene\npanel z", fontsize=7.8, color=MUTED)

        # Driver modifier badges (smaller, on the right side)
        mods = ["BRAF", "RAS", "TERT", "TP53", "fusion"]
        for m_i, m in enumerate(mods):
            mx = 8.0 + (m_i % 3) * 0.65
            my = y + 0.30 - (m_i // 3) * 0.45
            ax.add_patch(FancyBboxPatch((mx - 0.27, my - 0.16), 0.58, 0.32,
                                         boxstyle="round,pad=0,rounding_size=0.06",
                                         facecolor="white", edgecolor=MUTED, lw=0.5))
            ax.text(mx, my, m, ha="center", va="center", fontsize=6.5, color=MUTED)
        if i == 0:
            ax.text(8.65, y + 0.85, "driver modifiers", fontsize=7.5, color=MUTED, ha="center")

    # Arrow connecting states (vertical axis)
    ax.add_patch(FancyArrowPatch((0.35, 8.0), (0.35, 4.4), arrowstyle="->", mutation_scale=12, color=MUTED, lw=0.8))
    ax.text(0.15, 6.2, "silencing →", rotation=90, fontsize=8.3, color=MUTED, va="center")

# ------------- Panel d: evidence ladder -------------
def panel_d(ax):
    style_axis(ax)
    panel_label(ax, "d")
    panel_title(ax, "Evidence ladder")

    # 4 horizontal stages with arrows
    stages = [
        ("Discovery",            "TCGA-THCA + integrated\nthyroid bulk cohorts\n(n ≈ 1,500+, 4 pillars)", BLUE_L, BLUE),
        ("Panel definition",     "8-gene differentiation-\nsilencing axis\nSLC5A5 · TPO · TG · TSHR\nPAX8 · NKX2-1 · FOXE1 · DIO1", GRAY_L, GRAY_P),
        ("Label-anchored\nvalidation", "GSE151179 · GSE299988\nBoucai 2023 (request)\nMu 2024 HRA004166", RED_L, RED),
        ("Clinical stratification","RAI-avid · gray zone ·\nRAI-refractory groups", BEIGE, "#7a6038"),
    ]
    n = len(stages)
    box_w = 2.0; gap = 0.18; total = n * box_w + (n - 1) * gap
    start_x = (10 - total) / 2
    by = 5.0; bh = 4.2
    centers = []
    for i, (title, sub, face, edge) in enumerate(stages):
        bx = start_x + i * (box_w + gap)
        b = FancyBboxPatch((bx, by - bh/2), box_w, bh, boxstyle="round,pad=0,rounding_size=0.15",
                            facecolor=face, edgecolor=edge, lw=0.9)
        ax.add_patch(b)
        ax.text(bx + box_w/2, by + bh/2 - 0.45, title, ha="center", fontsize=9, color=edge, fontweight="bold")
        ax.text(bx + box_w/2, by - 0.05, sub, ha="center", fontsize=7.8, color=INK)
        centers.append((bx + box_w/2, by))
        # tier badge
        tier_text = ["Tier 4 (proxy)", "Panel", "Tier 1/2 (anchor)", "Manuscript"][i]
        tier_color = [MUTED, MUTED, RED, MUTED][i]
        ax.text(bx + box_w/2, by - bh/2 + 0.25, tier_text, ha="center", fontsize=7.0,
                color=tier_color, fontstyle="italic")
    # arrows between centers
    for i in range(n - 1):
        x1, y1 = centers[i]; x2, y2 = centers[i + 1]
        ax.add_patch(FancyArrowPatch((x1 + box_w/2 + 0.02, by), (x2 - box_w/2 - 0.02, by),
                                       arrowstyle="->", mutation_scale=14, color=MUTED, lw=0.9))

    # foot caption
    ax.text(5.0, 1.6,
            "Tier 4 proxy = discovery context (no direct RAI labels in TCGA).\n"
            "Tier 1/2 anchor = RECIST or RAI-avidity labels in independent cohorts.",
            ha="center", fontsize=7.8, color=MUTED, fontstyle="italic")


def main():
    fig = plt.figure(figsize=(14, 13.2), facecolor="white")
    # 2x2 grid; panels a, b, c, d
    gs = fig.add_gridspec(2, 2, hspace=0.55, wspace=0.20,
                          left=0.06, right=0.96, top=0.89, bottom=0.05)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    panel_a(ax_a); panel_b(ax_b); panel_c(ax_c); panel_d(ax_d)

    # Title and a horizontal separator
    fig.suptitle("A multi-gate model of radioiodine failure in thyroid cancer",
                 fontsize=15, fontweight="bold", color=INK, y=0.97, family="DejaVu Sans")

    out_png = OUT_DIR / "figure1_multigate_rai_model.png"
    out_pdf = OUT_DIR / "figure1_multigate_rai_model.pdf"
    fig.savefig(out_png, dpi=200, bbox_inches="tight", facecolor="white")
    fig.savefig(out_pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"wrote {out_png}")
    print(f"wrote {out_pdf}")


if __name__ == "__main__":
    main()
