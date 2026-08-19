#!/usr/bin/env python3
"""Figure 1 v2 — flagship Nature Medicine-style overview.

"A multi-gate model of radioiodine failure in thyroid cancer"

4-panel composition labeled A–D:
  A. Clinical unmet need with decision-node
  B. Multi-gate biology with explicit per-gate gene labels
  C. Three molecular states + horizontal gradient bar + modifier badges
  D. Tiered evidence ladder with 4 dataset tiles

16:9 aspect, warm-ivory background, muted refined palette, vector-only matplotlib.
"""
from __future__ import annotations
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import (FancyArrowPatch, FancyBboxPatch, Circle, Ellipse,
                                 Rectangle, RegularPolygon, PathPatch)
from matplotlib.path import Path as MplPath
import matplotlib.colors as mcolors
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = ROOT / "results" / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Refined Nature-Medicine-style palette
IVORY  = "#fbf8f1"
BLUE   = "#37618e"      # RAI-avid / preserved
BLUE_L = "#bccfe3"
RED    = "#9c4742"      # RAI-refractory / silenced
RED_L  = "#dfb3ae"
GRAY_P = "#7e6e94"      # gray zone
GRAY_L = "#c2b9d0"
BEIGE  = "#cdbb96"
BEIGE_L = "#ece1c3"
INK    = "#2a2a2a"
MUTED  = "#62656b"
GRID   = "#d9d4c5"
AMBER  = "#b8862a"      # I-131 particles


def style_axis(ax, xlim=(0, 100), ylim=(0, 100)):
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlim(*xlim); ax.set_ylim(*ylim)
    ax.set_aspect("auto")
    ax.set_facecolor(IVORY)


def panel_letter(ax, letter):
    ax.text(-0.012, 1.06, letter, transform=ax.transAxes, fontsize=17, fontweight="bold",
            color=INK, family="DejaVu Sans", va="top")


def panel_title(ax, text):
    ax.text(0.5, 1.04, text, transform=ax.transAxes, fontsize=12.5, fontweight="bold",
            ha="center", color=INK, family="DejaVu Sans", va="top")


# ============================================================
# Panel A · Clinical unmet need
# ============================================================
def panel_a(ax):
    style_axis(ax)
    panel_letter(ax, "a")
    panel_title(ax, "Clinical unmet need")

    # Thyroid butterfly icon
    for x_off, ang in ((-3.0, -14), (3.0, 14)):
        ax.add_patch(Ellipse((10 + x_off, 70), 6.5, 11, angle=ang,
                              facecolor=BEIGE_L, edgecolor=INK, lw=0.7))
    ax.text(10, 56, "Differentiated\nthyroid cancer", ha="center", fontsize=8.7, color=INK,
            family="DejaVu Sans")

    # Arrow → surgery (minimal scalpel)
    ax.add_patch(FancyArrowPatch((19, 70), (29, 70), arrowstyle="->",
                                  mutation_scale=11, color=MUTED, lw=0.9))
    # Scalpel — thin diagonal line + small rectangular handle
    ax.plot([30, 36], [73, 67], color=INK, lw=1.6, solid_capstyle="round")
    ax.add_patch(Rectangle((31.5, 65.5), 5.8, 1.4, facecolor=MUTED, edgecolor=INK, lw=0.5,
                            angle=-35))
    ax.text(33.5, 56, "Surgery", ha="center", fontsize=8.7, color=INK)

    # Arrow → I-131 capsule
    ax.add_patch(FancyArrowPatch((42, 70), (52, 70), arrowstyle="->",
                                  mutation_scale=11, color=MUTED, lw=0.9))
    cap = FancyBboxPatch((52, 65), 11, 10, boxstyle="round,pad=0,rounding_size=3.5",
                          facecolor=RED_L, edgecolor=RED, lw=0.8)
    ax.add_patch(cap)
    # capsule line dividing two halves
    ax.plot([57.5, 57.5], [65.5, 74.5], color=RED, lw=0.6)
    ax.text(57.5, 70, "I-131", fontsize=9, ha="center", va="center", color=RED, fontweight="bold")
    ax.text(57.5, 56, "Radioiodine therapy", ha="center", fontsize=8.7, color=INK)

    # Decision node / question mark
    ax.add_patch(Circle((73, 70), 4.0, facecolor=IVORY, edgecolor=MUTED, lw=1.0,
                        linestyle="--"))
    ax.text(73, 70, "?", ha="center", va="center", fontsize=15, color=MUTED, fontweight="bold")
    ax.text(73, 60, "decision node\n(no robust predictor)", ha="center", fontsize=7.5,
            color=MUTED, style="italic")

    # Split arrows (less steep — outcomes closer to therapy line)
    ax.add_patch(FancyArrowPatch((76.5, 73), (87, 82), arrowstyle="->",
                                  mutation_scale=11, color=BLUE, lw=0.9))
    ax.add_patch(FancyArrowPatch((76.5, 67), (87, 44), arrowstyle="->",
                                  mutation_scale=11, color=RED, lw=0.9))

    # Outcome 1: RAI-avid remission
    b1 = FancyBboxPatch((87, 78), 12, 10, boxstyle="round,pad=0,rounding_size=1.2",
                        facecolor=BLUE_L, edgecolor=BLUE, lw=0.9)
    ax.add_patch(b1)
    ax.text(93, 85, "RAI-avid", ha="center", fontsize=9.2, color=BLUE, fontweight="bold")
    ax.text(93, 81, "remission", ha="center", fontsize=8.4, color=INK)

    # Outcome 2: RAIR persistent / metastatic
    b2 = FancyBboxPatch((87, 35), 12, 12, boxstyle="round,pad=0,rounding_size=1.2",
                        facecolor=RED_L, edgecolor=RED, lw=0.9)
    ax.add_patch(b2)
    ax.text(93, 42, "RAI-refractory", ha="center", fontsize=9.2, color=RED, fontweight="bold")
    ax.text(93, 38, "persistent /\nmetastatic", ha="center", fontsize=8.0, color=INK)

    # Bottom caption strip
    ax.text(50, 18,
            "10-yr survival drops > 90% (responsive)  →  ~10–14% (distant RAI-refractory)",
            ha="center", fontsize=8.5, color=MUTED, fontstyle="italic")
    ax.text(50, 11,
            "Unmet need: a molecular readout that triages RAI candidates before therapy",
            ha="center", fontsize=8.2, color=MUTED)


# ============================================================
# Panel B · Multi-gate biology
# ============================================================
def panel_b(ax):
    style_axis(ax)
    panel_letter(ax, "b")
    panel_title(ax, "Multi-gate biology of RAI response")

    # Capillary on the left
    cap = FancyBboxPatch((4, 18), 11, 64, boxstyle="round,pad=0,rounding_size=2.5",
                          facecolor="#e0ecf6", edgecolor=BLUE, lw=0.8)
    ax.add_patch(cap)
    ax.text(9.5, 85, "blood", ha="center", fontsize=8.5, color=BLUE, fontweight="bold")
    for y in (28, 40, 52, 64, 76):
        ax.add_patch(Circle((9.5, y), 1.2, facecolor=AMBER, edgecolor="none", alpha=0.9))
        ax.text(11.5, y, "I⁻", fontsize=6.0, color=AMBER, va="center", style="italic")

    # Thyrocyte cell
    cell = FancyBboxPatch((20, 18), 50, 64, boxstyle="round,pad=0,rounding_size=3",
                           facecolor=BEIGE_L, edgecolor=INK, lw=0.9)
    ax.add_patch(cell)
    ax.text(45, 85, "thyroid cancer cell", ha="center", fontsize=8.5, color=INK, fontweight="bold")

    # Nucleus inside cell
    nuc = Ellipse((34, 50), 14, 16, facecolor="#ecddc1", edgecolor=MUTED, lw=0.7)
    ax.add_patch(nuc)
    ax.text(34, 50, "nucleus", ha="center", va="center", fontsize=7.0, color=MUTED, style="italic")

    # Colloid lumen
    col = FancyBboxPatch((76, 18), 18, 64, boxstyle="round,pad=0,rounding_size=3",
                          facecolor=BEIGE, edgecolor="#7c6a45", lw=0.9)
    ax.add_patch(col)
    ax.text(85, 85, "colloid", ha="center", fontsize=8.5, color="#5b4d2e", fontweight="bold")
    # TG/iodinated thyroglobulin scatter
    rng = np.random.default_rng(11)
    for _ in range(22):
        x = rng.uniform(78, 92)
        y = rng.uniform(30, 75)
        ax.add_patch(Circle((x, y), 0.5, facecolor="#85703f", alpha=0.55, edgecolor="none"))

    # Iodine path: blood → cell → colloid
    ax.add_patch(Circle((25, 70), 1.0, facecolor=AMBER, edgecolor="none", alpha=0.8))
    ax.add_patch(Circle((30, 65), 0.9, facecolor=AMBER, edgecolor="none", alpha=0.8))
    ax.add_patch(Circle((42, 35), 0.9, facecolor=AMBER, edgecolor="none", alpha=0.8))
    ax.add_patch(Circle((60, 33), 1.0, facecolor=AMBER, edgecolor="none", alpha=0.8))
    ax.add_patch(Circle((78, 35), 1.1, facecolor=AMBER, edgecolor="none", alpha=0.8))

    # GATE 1: SLC5A5 / NIS — at blood-cell membrane
    g1y = 70
    ax.add_patch(FancyArrowPatch((15, g1y), (20, g1y), arrowstyle="->",
                                  mutation_scale=12, color=RED, lw=1.3))
    g1 = FancyBboxPatch((15.5, 60), 9.0, 8.0, boxstyle="round,pad=0,rounding_size=1.2",
                        facecolor=RED_L, edgecolor=RED, lw=0.7)
    ax.add_patch(g1)
    ax.text(20, 65.5, "Gate 1", fontsize=7.0, ha="center", color=RED, fontweight="bold")
    ax.text(20, 61.5, "SLC5A5 / NIS", fontsize=7.0, ha="center", color=INK)

    # GATE 2: lineage TFs in nucleus
    g2 = FancyBboxPatch((26, 32), 16, 8.5, boxstyle="round,pad=0,rounding_size=1.2",
                        facecolor=BLUE_L, edgecolor=BLUE, lw=0.7)
    ax.add_patch(g2)
    ax.text(34, 38, "Gate 2  ·  thyroid lineage TFs", fontsize=7.0, ha="center", color=BLUE,
            fontweight="bold")
    ax.text(34, 34, "PAX8 · NKX2-1 · FOXE1 · TSHR", fontsize=7.2, ha="center", color=INK)
    # arrow nucleus → cell programs
    ax.add_patch(FancyArrowPatch((34, 42), (34, 47), arrowstyle="-",
                                  mutation_scale=10, color=BLUE, lw=0.6, linestyle="--"))

    # GATE 3: organification at apical membrane
    g3 = FancyBboxPatch((50, 55), 22, 13, boxstyle="round,pad=0,rounding_size=1.5",
                        facecolor=BLUE_L, edgecolor=BLUE, lw=0.7)
    ax.add_patch(g3)
    ax.text(61, 64, "Gate 3  ·  organification & storage", fontsize=7.2, ha="center", color=BLUE,
            fontweight="bold")
    ax.text(61, 60, "TG · TPO · DIO1", fontsize=7.4, ha="center", color=INK)
    ax.text(61, 56.5, "(DUOX1/2 supporting)", fontsize=6.3, ha="center", color=MUTED, style="italic")
    # arrow cell → colloid
    ax.add_patch(FancyArrowPatch((72, 58), (76, 50), arrowstyle="->",
                                  mutation_scale=11, color=BLUE, lw=1.0))

    # GATE 4: retention + radiation killing — inside colloid + DNA damage
    # radiation rays
    cx, cy = 85, 50
    for ang in range(0, 360, 30):
        r = 4.5
        x1 = cx + r * np.cos(np.deg2rad(ang))
        y1 = cy + r * np.sin(np.deg2rad(ang))
        ax.plot([cx, x1], [cy, y1], color=RED, lw=0.6)
    ax.add_patch(Circle((cx, cy), 1.3, facecolor=RED, edgecolor="white", lw=0.6))
    g4 = FancyBboxPatch((75.5, 23), 17, 11, boxstyle="round,pad=0,rounding_size=1.5",
                        facecolor=RED_L, edgecolor=RED, lw=0.7)
    ax.add_patch(g4)
    ax.text(84, 30.5, "Gate 4  ·  retention", fontsize=7.2, ha="center", color=RED,
            fontweight="bold")
    ax.text(84, 26.5, "→ DNA damage → kill", fontsize=7.0, ha="center", color=INK)

    # bottom caption strip
    cap_box = FancyBboxPatch((10, 4), 84, 8, boxstyle="round,pad=0,rounding_size=1.5",
                              facecolor="white", edgecolor=GRID, lw=0.7)
    ax.add_patch(cap_box)
    ax.text(52, 8.0, "NIS alone is insufficient  —  coordinated differentiation across all four gates is required",
            ha="center", va="center", fontsize=9.5, color=INK, fontweight="bold", style="italic")


# ============================================================
# Panel C · Three molecular states + gradient bar
# ============================================================
def panel_c(ax):
    style_axis(ax)
    panel_letter(ax, "c")
    panel_title(ax, "Three molecular states along the eight-gene differentiation-silencing axis")

    states = [
        # (label, sub, face, edge, gate_states, score_text, x_center)
        ("RAI-avid",      "preserved differentiation",
         BLUE_L,  BLUE,   ["open"]*4, "high score",    18),
        ("Gray zone",     "partial lineage loss",
         GRAY_L,  GRAY_P, ["open", "open", "half", "closed"], "intermediate", 50),
        ("RAI-refractory","differentiation silencing",
         RED_L,   RED,    ["closed"]*4, "low score",    82),
    ]

    bx_h = 35
    by = 60
    box_w = 26
    for label, sub, face, edge, gate_states, score_text, cx in states:
        b = FancyBboxPatch((cx - box_w/2, by - bx_h/2), box_w, bx_h,
                           boxstyle="round,pad=0,rounding_size=1.5",
                           facecolor=face, edgecolor=edge, lw=0.95)
        ax.add_patch(b)
        ax.text(cx, by + 13, label, ha="center", fontsize=11.5, fontweight="bold", color=edge)
        ax.text(cx, by + 9, sub, ha="center", fontsize=8.5, color=INK, style="italic")
        # gate icons row
        gate_x = np.linspace(cx - 9, cx + 9, 4)
        for gx, gs in zip(gate_x, gate_states):
            if gs == "open":
                ax.add_patch(Circle((gx, by + 2), 1.4, facecolor=edge, edgecolor="white", lw=0.6))
            elif gs == "half":
                ax.add_patch(Circle((gx, by + 2), 1.4, facecolor="white", edgecolor=edge, lw=0.8))
                ax.add_patch(Rectangle((gx - 1.4, by + 0.6), 1.4, 2.8, facecolor=edge,
                                        edgecolor="none"))
            else:
                ax.add_patch(Circle((gx, by + 2), 1.4, facecolor="white", edgecolor="#9a9a9a",
                                     lw=0.7))
                ax.plot([gx - 0.9, gx + 0.9], [by + 1.1, by + 2.9], color="#9a9a9a", lw=0.7)
                ax.plot([gx - 0.9, gx + 0.9], [by + 2.9, by + 1.1], color="#9a9a9a", lw=0.7)
        ax.text(cx, by - 3.5, "gates  1   2   3   4", fontsize=6.8, ha="center", color=MUTED)
        ax.text(cx, by - 9, score_text, fontsize=8.5, ha="center", color=edge, fontweight="bold")

    # Horizontal gradient bar — 8-gene axis
    grad_y = 28
    grad_x0, grad_x1 = 10, 90
    n_seg = 200
    cmap = mcolors.LinearSegmentedColormap.from_list("axis",
                                                       [BLUE, GRAY_P, RED], N=n_seg)
    for i in range(n_seg):
        col = cmap(i / (n_seg - 1))
        ax.add_patch(Rectangle((grad_x0 + i * (grad_x1 - grad_x0) / n_seg, grad_y - 1.5),
                                (grad_x1 - grad_x0) / n_seg + 0.05, 3,
                                facecolor=col, edgecolor="none"))
    # bar end caps + axis ticks
    ax.add_patch(Rectangle((grad_x0, grad_y - 1.5), grad_x1 - grad_x0, 3,
                            facecolor="none", edgecolor=INK, lw=0.7))
    ax.text(grad_x0 - 1, grad_y, "preserved", ha="right", va="center", fontsize=8.5, color=BLUE,
            fontweight="bold")
    ax.text(grad_x1 + 1, grad_y, "silenced", ha="left", va="center", fontsize=8.5, color=RED,
            fontweight="bold")
    ax.text(50, grad_y + 5.5, "8-gene differentiation-silencing axis",
            ha="center", fontsize=9.2, color=INK, fontweight="bold")
    ax.text(50, grad_y - 4.5,
            "SLC5A5 · TPO · TG · TSHR · PAX8 · NKX2-1 · FOXE1 · DIO1",
            ha="center", fontsize=7.8, color=MUTED, family="DejaVu Sans")

    # Driver modifier badges below the gradient
    badges = ["BRAF", "RAS", "TERT", "TP53", "Fusion"]
    bx_y = 15
    bx_w = 9.6
    bx_h2 = 4.0
    total = len(badges) * bx_w + (len(badges) - 1) * 2
    start_x = 50 - total / 2
    for i, m in enumerate(badges):
        bx = start_x + i * (bx_w + 2)
        ax.add_patch(FancyBboxPatch((bx, bx_y - bx_h2/2), bx_w, bx_h2,
                                     boxstyle="round,pad=0,rounding_size=0.8",
                                     facecolor="white", edgecolor=MUTED, lw=0.6))
        ax.text(bx + bx_w/2, bx_y, m, ha="center", va="center", fontsize=7.8, color=MUTED)
    ax.text(50, 9.5, "genomic modifiers  (secondary, not the main axis)",
            ha="center", fontsize=7.5, color=MUTED, style="italic")


# ============================================================
# Panel D · Evidence ladder
# ============================================================
def panel_d(ax):
    style_axis(ax)
    panel_letter(ax, "d")
    panel_title(ax, "Tiered evidence ladder  ·  discovery → validation → stratification")

    # Stage boxes (left → right)
    stages = [
        ("Discovery", "TCGA-THCA + integrated\nthyroid cohorts",
         "Tier 4 molecular discovery", BLUE_L, BLUE),
        ("8-gene panel", "SLC5A5 · TPO · TG · TSHR\nPAX8 · NKX2-1 · FOXE1 · DIO1",
         "differentiation-silencing axis", GRAY_L, GRAY_P),
        ("Label-anchored\nvalidation", None, "Tier 1 / Tier 2 anchors", RED_L, RED),
        ("Clinical\nstratification", "RAI-avid  ·  Gray zone\nRAI-refractory",
         "tiered patient interpretation", BEIGE, "#7c6a45"),
    ]
    n = len(stages)
    box_w = 18
    gap = 4
    total = n * box_w + (n - 1) * gap
    start_x = 50 - total / 2
    by = 70
    bh = 22
    centers = []
    for i, (title, sub, footer, face, edge) in enumerate(stages):
        bx = start_x + i * (box_w + gap)
        b = FancyBboxPatch((bx, by - bh/2), box_w, bh,
                           boxstyle="round,pad=0,rounding_size=1.5",
                           facecolor=face, edgecolor=edge, lw=0.95)
        ax.add_patch(b)
        ax.text(bx + box_w/2, by + 7, title, ha="center", fontsize=9.5, fontweight="bold",
                color=edge)
        if sub:
            ax.text(bx + box_w/2, by + 0.5, sub, ha="center", fontsize=7.7, color=INK)
        ax.text(bx + box_w/2, by - 8.5, footer, ha="center", fontsize=7.2, color=MUTED,
                style="italic")
        centers.append((bx + box_w/2, by))

    # Arrows between stages
    for i in range(n - 1):
        x1, y1 = centers[i]; x2, y2 = centers[i + 1]
        ax.add_patch(FancyArrowPatch((x1 + box_w/2 + 0.3, by), (x2 - box_w/2 - 0.3, by),
                                       arrowstyle="->", mutation_scale=14, color=MUTED, lw=1.0))

    # Dataset tiles UNDER the "Label-anchored validation" stage — single horizontal row
    third_bx = start_x + 2 * (box_w + gap)
    third_cx = third_bx + box_w / 2
    tiles = [
        ("GSE151179",   "RAI avidity",       "Tier 2 · public"),
        ("GSE299988",   "supportive small",  "Tier 2 · public"),
        ("Boucai 2023", "RECIST exceptional","Tier 1 · request"),
        ("Mu 2024",     "4-class uptake",    "Tier 2 · HRA004166"),
    ]
    tile_w = 17.0; tile_h = 11.0; tile_gap = 2.5
    total_tile_w = len(tiles) * tile_w + (len(tiles) - 1) * tile_gap
    tile_start_x = 50 - total_tile_w / 2
    tile_y = 30
    # connecting arrow from stage 3 bottom to tile row
    ax.add_patch(FancyArrowPatch((third_cx, by - bh/2 - 0.3),
                                  (third_cx, tile_y + tile_h + 1.5),
                                  arrowstyle="-", mutation_scale=10, color=MUTED, lw=0.6,
                                  linestyle=":"))
    for i, (name, kind, access) in enumerate(tiles):
        tx = tile_start_x + i * (tile_w + tile_gap)
        t = FancyBboxPatch((tx, tile_y), tile_w, tile_h,
                           boxstyle="round,pad=0,rounding_size=1.0",
                           facecolor="white", edgecolor=RED, lw=0.7)
        ax.add_patch(t)
        ax.text(tx + tile_w/2, tile_y + tile_h - 3.0, name, ha="center", fontsize=8.5, color=RED,
                fontweight="bold")
        ax.text(tx + tile_w/2, tile_y + tile_h/2 - 0.5, kind, ha="center", fontsize=7.4, color=INK)
        ax.text(tx + tile_w/2, tile_y + 2.2, access, ha="center", fontsize=6.7, color=MUTED,
                style="italic")

    # bottom caption strip
    ax.text(50, 17,
            "Tier 4 discovery  →  parsimonious 8-gene panel  →  label-anchored validation  →  three-class stratification",
            ha="center", fontsize=8.5, color=MUTED, fontstyle="italic")
    ax.text(50, 11,
            "Framework = risk stratification, not a stand-alone clinical biomarker; prospective Tier-1 confirmation is the next step.",
            ha="center", fontsize=7.5, color=MUTED, style="italic")


# ============================================================
def main():
    fig = plt.figure(figsize=(17, 11), facecolor=IVORY)
    gs = fig.add_gridspec(2, 2, hspace=0.40, wspace=0.13,
                          left=0.045, right=0.97, top=0.90, bottom=0.045)

    fig.suptitle("A multi-gate model of radioiodine failure in thyroid cancer",
                 fontsize=16, fontweight="bold", color=INK, y=0.965,
                 family="DejaVu Sans")
    fig.text(0.5, 0.929,
             "an eight-gene differentiation / iodide-handling axis links discovery cohorts to RAI-labelled validation",
             ha="center", fontsize=10.5, color=MUTED, style="italic")

    ax_a = fig.add_subplot(gs[0, 0]); ax_a.set_facecolor(IVORY)
    ax_b = fig.add_subplot(gs[0, 1]); ax_b.set_facecolor(IVORY)
    ax_c = fig.add_subplot(gs[1, 0]); ax_c.set_facecolor(IVORY)
    ax_d = fig.add_subplot(gs[1, 1]); ax_d.set_facecolor(IVORY)

    panel_a(ax_a); panel_b(ax_b); panel_c(ax_c); panel_d(ax_d)

    out_png = OUT_DIR / "figure1_v2_flagship.png"
    out_pdf = OUT_DIR / "figure1_v2_flagship.pdf"
    fig.savefig(out_png, dpi=210, bbox_inches="tight", facecolor=IVORY)
    fig.savefig(out_pdf, bbox_inches="tight", facecolor=IVORY)
    plt.close(fig)
    print(f"wrote {out_png}")
    print(f"wrote {out_pdf}")


if __name__ == "__main__":
    main()
