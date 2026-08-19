#!/usr/bin/env python3
"""Figure 1 v3 — Nature-illustrator grade flagship overview.

Push matplotlib past default look using Path + Bezier curves for organic shapes,
layered gradients via concentric ellipses, anatomically faithful thyroid follicle,
proper NIS transporter glyph, radiation decay iconography, and a real gradient
axis bar in panel c.
"""
from __future__ import annotations
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
from _nature_style import setup_rc, panel_letter, IVORY, INK, MUTED

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import (FancyArrowPatch, FancyBboxPatch, Circle, Ellipse,
                                 Rectangle, PathPatch, Polygon, Arc, RegularPolygon)
from matplotlib.path import Path as MplPath
from matplotlib.colors import LinearSegmentedColormap, to_rgba

ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = ROOT / "results" / "figures"

# Refined Nature Medicine palette (warmer, more layered)
PAPER = "#fcf9f3"
IVORY_W = "#fbf8f1"
BLUE_D = "#244b75"
BLUE   = "#3a6694"
BLUE_M = "#7a9bc1"
BLUE_L = "#c9d8e8"
BLUE_VL = "#e6eef6"
RED_D  = "#7b3833"
RED    = "#a04e48"
RED_M  = "#c98581"
RED_L  = "#e2bcb8"
RED_VL = "#f3dfdc"
PURP_D = "#5a4470"
PURP   = "#82689e"
PURP_L = "#c7b8d4"
PURP_VL = "#e9e1f0"
BEIGE_D = "#7d6736"
BEIGE   = "#b89b5e"
BEIGE_M = "#d4bf87"
BEIGE_L = "#ecd9a8"
BEIGE_VL = "#f5e9c8"
INK_D = "#1a1a1a"
INK_  = "#2a2a2a"
MUTED_  = "#5b5e64"
LINE = "#404048"
AMBER_D = "#8c5c1e"
AMBER = "#c08b35"
AMBER_L = "#e5b96b"


def fade_radial(ax, cx, cy, r_outer, color, alpha_inner=0.55, alpha_outer=0.0, n=24):
    """Simulate a radial gradient with concentric semi-transparent circles."""
    for i in range(n, 0, -1):
        frac = i / n
        a = alpha_inner * (1 - frac) ** 1.5 + alpha_outer * frac
        ax.add_patch(Circle((cx, cy), r_outer * frac, facecolor=color,
                            edgecolor="none", alpha=a, zorder=0.5))


def follicle_section(ax, cx, cy, w, h):
    """Draw an anatomical thyroid follicle cross-section using Bezier curves.

    Returns (x_basal, x_apical) approximate membrane boundaries for arrow targeting.
    """
    # Outer basement membrane (slightly oval, smooth Bezier)
    outer = MplPath([
        (cx - w/2, cy),
        (cx - w/2, cy + h*0.55), (cx + w/2, cy + h*0.55), (cx + w/2, cy),
        (cx + w/2, cy - h*0.55), (cx - w/2, cy - h*0.55), (cx - w/2, cy),
    ], [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
        MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4])
    ax.add_patch(PathPatch(outer, facecolor=BEIGE_VL, edgecolor=BEIGE, lw=0.7, zorder=1))

    # Inner colloid lumen
    lw_, lh_ = w*0.55, h*0.35
    colloid = MplPath([
        (cx - lw_/2, cy),
        (cx - lw_/2, cy + lh_), (cx + lw_/2, cy + lh_), (cx + lw_/2, cy),
        (cx + lw_/2, cy - lh_), (cx - lw_/2, cy - lh_), (cx - lw_/2, cy),
    ], [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
        MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4])
    # subtle colloid radial fade
    fade_radial(ax, cx, cy, lw_*0.6, BEIGE_M, alpha_inner=0.45)
    ax.add_patch(PathPatch(colloid, facecolor=BEIGE_L, edgecolor=BEIGE_D, lw=0.8, zorder=2))

    # TG / iodinated storage stippling
    rng = np.random.default_rng(11)
    for _ in range(38):
        # within the colloid
        rx = rng.normal(cx, lw_*0.25)
        ry = rng.normal(cy, lh_*0.32)
        if abs(rx - cx)/lw_ * 2 + abs(ry - cy)/lh_ * 2 < 1.05:
            r_s = rng.uniform(0.16, 0.34)
            ax.add_patch(Circle((rx, ry), r_s, facecolor=BEIGE_D, alpha=rng.uniform(0.35, 0.7),
                                edgecolor="none", zorder=3))

    # Apical microvilli (small bumps facing colloid) — top arc
    n_micro = 14
    for i in range(n_micro):
        ang = np.pi/2 - 0.55 + (i / (n_micro - 1)) * 1.1
        x_o = cx + (lw_/2 + 0.05) * np.cos(ang)
        y_o = cy + (lh_/2 + 0.05) * np.sin(ang)
        x_i = cx + (lw_/2 - 0.18) * np.cos(ang)
        y_i = cy + (lh_/2 - 0.18) * np.sin(ang)
        ax.plot([x_i, x_o], [y_i, y_o], color=BEIGE_D, lw=0.5, zorder=4)

    return cx - w/2, cx + w/2, lw_, lh_


def thyrocyte_cell(ax, cx, cy, r_radial, color_membrane=BLUE_D, color_cyto=BLUE_VL,
                    nucleus_offset=(0, 0), nucleus_size=(1.0, 1.2)):
    """Draw a layered thyrocyte cell with membrane, cytoplasm, nucleus."""
    # Cell membrane (slightly oval)
    fade_radial(ax, cx, cy, r_radial*1.2, color_cyto, alpha_inner=0.55)
    ax.add_patch(Ellipse((cx, cy), r_radial*2, r_radial*2.3, facecolor=color_cyto,
                         edgecolor=color_membrane, lw=0.8, zorder=5, alpha=0.85))
    # Nucleus
    nx, ny = cx + nucleus_offset[0], cy + nucleus_offset[1]
    nw, nh = nucleus_size
    fade_radial(ax, nx, ny, nw*1.3, INK_D, alpha_inner=0.18, n=14)
    ax.add_patch(Ellipse((nx, ny), nw*2, nh*2, facecolor="#dccba8", edgecolor=BEIGE_D,
                         lw=0.6, zorder=6, alpha=0.95))
    # Lineage TFs as small dots inside nucleus
    rng = np.random.default_rng(3)
    for _ in range(11):
        dx = rng.uniform(-nw*0.7, nw*0.7); dy = rng.uniform(-nh*0.65, nh*0.65)
        if (dx/nw)**2 + (dy/nh)**2 < 0.85:
            ax.add_patch(Circle((nx+dx, ny+dy), 0.10, facecolor=BLUE_D,
                                edgecolor="none", alpha=0.75, zorder=7))
    return nx, ny


def nis_transporter(ax, cx, cy, w=0.7, h=1.2, color=RED_D):
    """A stylized NIS membrane transporter glyph (two helical bundles)."""
    # outer barrel
    for dx in (-w*0.32, w*0.32):
        ax.add_patch(FancyBboxPatch((cx + dx - 0.12, cy - h/2), 0.24, h,
                                    boxstyle="round,pad=0,rounding_size=0.10",
                                    facecolor=color, edgecolor=RED_D, lw=0.4, zorder=8))
    # channel
    ax.add_patch(Rectangle((cx - 0.10, cy - h/2 + 0.05), 0.20, h - 0.10,
                            facecolor=RED_VL, edgecolor="none", zorder=8.1))


def iodine_particle(ax, x, y, size=0.22, color=AMBER, glow=True):
    if glow:
        ax.add_patch(Circle((x, y), size*2.2, facecolor=color, edgecolor="none",
                            alpha=0.10, zorder=4.5))
    ax.add_patch(Circle((x, y), size, facecolor=color, edgecolor=AMBER_D, lw=0.3,
                        alpha=0.95, zorder=4.6))


def decay_rays(ax, cx, cy, n=16, r_inner=0.25, r_outer=0.9, color=RED_D):
    """Radiation decay symbol — alternating ray lengths around a point."""
    for i in range(n):
        ang = i * (2 * np.pi / n)
        ro = r_outer if i % 2 == 0 else r_outer * 0.75
        x0 = cx + r_inner * np.cos(ang); y0 = cy + r_inner * np.sin(ang)
        x1 = cx + ro * np.cos(ang); y1 = cy + ro * np.sin(ang)
        ax.plot([x0, x1], [y0, y1], color=color, lw=0.9, alpha=0.85, zorder=9,
                solid_capstyle="round")
    ax.add_patch(Circle((cx, cy), r_inner*0.55, facecolor=color, edgecolor="white",
                        lw=0.5, zorder=10))


# ============================================================
# Panel a — Clinical unmet need (anatomical thyroid + outcome split)
# ============================================================
def panel_a(ax):
    ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.set_facecolor(PAPER)
    for s in ax.spines.values(): s.set_visible(False)
    ax.set_xticks([]); ax.set_yticks([])
    ax.text(50, 95, "Clinical unmet need", ha="center", fontsize=11.5, fontweight="bold",
            color=INK_D)

    # 1) Thyroid icon (butterfly with two lobes + isthmus)
    cx_th, cy_th = 11, 70
    for x_off, ang, lobe_c in ((-2.5, -16, BEIGE_L), (2.5, 16, BEIGE_L)):
        ax.add_patch(Ellipse((cx_th + x_off, cy_th), 4.2, 8.0, angle=ang,
                              facecolor=lobe_c, edgecolor=BEIGE_D, lw=0.6, alpha=0.92))
    # isthmus
    ax.add_patch(Rectangle((cx_th - 1.8, cy_th - 0.6), 3.6, 2.0, facecolor=BEIGE_L,
                            edgecolor=BEIGE_D, lw=0.5))
    # vessel hint
    ax.plot([cx_th - 4, cx_th + 4], [cy_th + 4, cy_th + 4.5], color=RED_M, lw=0.7, alpha=0.55)
    ax.text(cx_th, cy_th - 9, "Differentiated\nthyroid cancer", ha="center", fontsize=8.7, color=INK_)
    ax.text(cx_th, cy_th - 13, "(post-surgery)", ha="center", fontsize=7.5, color=MUTED_, style="italic")

    # 2) Arrow → I-131 capsule
    ax.add_patch(FancyArrowPatch((19, 70), (29, 70), arrowstyle="-|>",
                                  mutation_scale=11, color=MUTED_, lw=0.9))
    # Capsule with subtle two-tone halves
    ax.add_patch(FancyBboxPatch((29.5, 65), 12, 10, boxstyle="round,pad=0,rounding_size=4.5",
                                facecolor=RED_VL, edgecolor=RED_D, lw=0.8))
    ax.add_patch(FancyBboxPatch((29.5, 65), 5.8, 10, boxstyle="round,pad=0,rounding_size=4.5",
                                facecolor=RED_L, edgecolor="none", alpha=0.85))
    ax.plot([35.3, 35.3], [65.5, 74.5], color=RED_D, lw=0.5)
    # decay tick marks on capsule
    for i in range(3):
        y = 67.5 + i * 2.5
        ax.plot([37.5, 39.5], [y, y], color=RED_D, lw=0.5)
    ax.text(35.5, 70, "¹³¹I", ha="center", va="center", fontsize=9.5, color=RED_D, fontweight="bold")
    ax.text(35.5, 56, "Radioiodine therapy", ha="center", fontsize=8.7, color=INK_)

    # 3) Decision diamond (no robust predictor)
    diamond = MplPath([(50, 76), (56, 70), (50, 64), (44, 70), (50, 76)],
                      [MplPath.MOVETO]+[MplPath.LINETO]*3+[MplPath.CLOSEPOLY])
    ax.add_patch(PathPatch(diamond, facecolor=PAPER, edgecolor=MUTED_, lw=0.9, linestyle="--"))
    ax.text(50, 70, "?", ha="center", va="center", fontsize=14, color=MUTED_, fontweight="bold")
    ax.text(50, 58, "no robust\npre-RAI predictor", ha="center", fontsize=7.4, color=MUTED_, style="italic")

    # 4) Two branches with curved arrows
    arrow_up = FancyArrowPatch((57, 73), (76, 84), arrowstyle="-|>",
                                mutation_scale=11, color=BLUE_D, lw=0.9,
                                connectionstyle="arc3,rad=-0.18")
    arrow_dn = FancyArrowPatch((57, 67), (76, 35), arrowstyle="-|>",
                                mutation_scale=11, color=RED_D, lw=0.9,
                                connectionstyle="arc3,rad=0.18")
    ax.add_patch(arrow_up); ax.add_patch(arrow_dn)

    # Outcomes — with checkmark / cross icons
    # RAI-avid (top)
    ax.add_patch(FancyBboxPatch((76, 79), 22, 12, boxstyle="round,pad=0,rounding_size=1.3",
                                 facecolor=BLUE_L, edgecolor=BLUE_D, lw=0.9))
    # check icon
    ax.plot([78.5, 80.5, 84], [85, 82.5, 88], color=BLUE_D, lw=2, solid_capstyle="round")
    ax.text(89, 87, "RAI-avid", fontsize=10, color=BLUE_D, fontweight="bold")
    ax.text(89, 83.5, "remission · ≥ 90% 10-yr DSS", fontsize=8, color=INK_)

    # RAI-refractory (bottom)
    ax.add_patch(FancyBboxPatch((76, 28), 22, 12, boxstyle="round,pad=0,rounding_size=1.3",
                                 facecolor=RED_L, edgecolor=RED_D, lw=0.9))
    # cross icon
    ax.plot([78, 82], [32, 36], color=RED_D, lw=2, solid_capstyle="round")
    ax.plot([78, 82], [36, 32], color=RED_D, lw=2, solid_capstyle="round")
    ax.text(89, 36, "RAI-refractory", fontsize=10, color=RED_D, fontweight="bold")
    ax.text(89, 32.5, "persistent / metastatic disease", fontsize=8, color=INK_)
    ax.text(89, 29.5, "~10–14% 10-yr DSS  (distant)", fontsize=8, color=RED_D, style="italic")

    # Bottom take-home
    ax.add_patch(FancyBboxPatch((4, 8), 92, 11, boxstyle="round,pad=0,rounding_size=1.2",
                                 facecolor="#fff7e2", edgecolor=BEIGE_D, lw=0.6, alpha=0.7))
    ax.text(50, 14.5, "Need: a parsimonious molecular readout to triage RAI candidates before therapy",
            ha="center", fontsize=9.5, color=INK_D, fontweight="bold")
    ax.text(50, 10.5, "Current biomarkers (BRAF / RAS / TERT) do not reliably partition the RAI-refractory population.",
            ha="center", fontsize=8.0, color=MUTED_, style="italic")


# ============================================================
# Panel b — Multi-gate biology (anatomical follicle)
# ============================================================
def panel_b(ax):
    ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.set_facecolor(PAPER)
    for s in ax.spines.values(): s.set_visible(False)
    ax.set_xticks([]); ax.set_yticks([])
    ax.text(50, 95, "Multi-gate biology of RAI response", ha="center", fontsize=11.5,
            fontweight="bold", color=INK_D)

    # Capillary on far left (curved vessel)
    cap_x = 6
    cap = MplPath([(cap_x - 3, 18), (cap_x - 3.5, 50), (cap_x - 3, 82),
                   (cap_x + 3, 82), (cap_x + 3.5, 50), (cap_x + 3, 18),
                   (cap_x - 3, 18)],
                  [MplPath.MOVETO, MplPath.CURVE3, MplPath.CURVE3,
                   MplPath.LINETO, MplPath.CURVE3, MplPath.CURVE3, MplPath.CLOSEPOLY])
    ax.add_patch(PathPatch(cap, facecolor="#e8f0f8", edgecolor=BLUE_D, lw=0.7))
    ax.text(cap_x, 87, "blood", ha="center", fontsize=8.5, color=BLUE_D, fontweight="bold")
    for y in (28, 38, 50, 62, 74):
        iodine_particle(ax, cap_x + np.random.default_rng(int(y)).uniform(-1.5, 1.5), y)

    # Thyrocyte layer — 3 layered cells around the follicle
    cell_centers = []
    for i, (cx, cy) in enumerate([(28, 78), (28, 50), (28, 22)]):
        nx, ny = thyrocyte_cell(ax, cx, cy, r_radial=4.2,
                                 color_membrane=BLUE_D, color_cyto=BLUE_VL,
                                 nucleus_offset=(0, 0), nucleus_size=(2.0, 2.4))
        cell_centers.append((cx, cy, nx, ny))

    # Follicle (centered apical → colloid lumen)
    fcx, fcy = 55, 50
    bx_l, bx_r, lw_, lh_ = follicle_section(ax, fcx, fcy, w=24, h=42)

    # NIS transporters embedded in basal membranes (between blood and cells)
    for cx, cy, _, _ in cell_centers:
        nis_transporter(ax, cx - 4.4, cy, w=0.9, h=1.4)
        # iodine moving through
        iodine_particle(ax, cx - 6, cy, size=0.20)
        iodine_particle(ax, cx - 2.5, cy, size=0.18)

    # Iodine flowing into colloid (path from cells → lumen)
    for cx, cy, _, _ in cell_centers:
        iodine_particle(ax, cx + 5, cy, size=0.20)
        iodine_particle(ax, fcx - 9 + (cx - 28) * 0.0, cy, size=0.20)

    # Colloid TG/I storage (already drawn) — plus radiation killing at right
    decay_rays(ax, fcx + 2, fcy, n=18, r_inner=0.5, r_outer=2.2)

    # Gate annotation boxes — clean labels with leader lines
    gate_anno = [
        # (gate_num, label, x_label, y_label, leader_x, leader_y, color)
        (1, "SLC5A5 / NIS\niodide uptake",     16, 76, 23.5, 74, RED_D),
        (2, "PAX8 · NKX2-1 · FOXE1\nTSHR · lineage program", 18, 35, 25, 47, BLUE_D),
        (3, "TG · TPO · DIO1\norganification & storage", 75, 78, 65, 65, BLUE_D),
        (4, "retention →\nradiation killing", 75, 28, 60, 47, RED_D),
    ]
    for gn, txt, lx, ly, tx, ty, c in gate_anno:
        # background pill behind label
        ax.add_patch(FancyBboxPatch((lx - 9.5, ly - 4), 19, 8.5,
                                     boxstyle="round,pad=0,rounding_size=1.0",
                                     facecolor="white", edgecolor=c, lw=0.7, alpha=0.95, zorder=15))
        ax.text(lx, ly + 1.4, f"Gate {gn}", ha="center", fontsize=8.5, color=c, fontweight="bold")
        ax.text(lx, ly - 1.8, txt, ha="center", fontsize=7.7, color=INK_)
        # leader line
        ax.add_patch(FancyArrowPatch((lx, ly - 4), (tx, ty), arrowstyle="-",
                                      mutation_scale=8, color=c, lw=0.5, alpha=0.6,
                                      connectionstyle="arc3,rad=0.0", zorder=14))

    # Bottom take-home banner
    ax.add_patch(FancyBboxPatch((6, 6), 88, 8, boxstyle="round,pad=0,rounding_size=1.0",
                                 facecolor=RED_VL, edgecolor=RED_D, lw=0.6, alpha=0.7))
    ax.text(50, 10, "NIS alone is insufficient — coordinated thyroid differentiation across all four gates is required",
            ha="center", fontsize=9.2, color=INK_D, fontweight="bold", style="italic")


# ============================================================
# Panel c — Three molecular states + gradient axis bar
# ============================================================
def panel_c(ax):
    ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.set_facecolor(PAPER)
    for s in ax.spines.values(): s.set_visible(False)
    ax.set_xticks([]); ax.set_yticks([])
    ax.text(50, 95, "Three molecular states along the eight-gene differentiation–silencing axis",
            ha="center", fontsize=11.5, fontweight="bold", color=INK_D)

    # State boxes with proper layered design
    states = [
        ("RAI-avid", "preserved differentiation\nall gates open  ·  high score",
         BLUE_L, BLUE_D, 18, "open"),
        ("Gray zone", "partial lineage loss\nintermediate score",
         PURP_L, PURP_D, 50, "half"),
        ("RAI-refractory", "differentiation silencing\ncollapsed gates  ·  low score",
         RED_L, RED_D, 82, "closed"),
    ]
    by = 65; bh = 26; bw = 24
    for label, sub, face, edge, cx, gate_state in states:
        # soft drop shadow
        ax.add_patch(FancyBboxPatch((cx - bw/2 + 0.4, by - bh/2 - 0.5), bw, bh,
                                     boxstyle="round,pad=0,rounding_size=1.6",
                                     facecolor="#00000010", edgecolor="none", zorder=0.5))
        # box
        ax.add_patch(FancyBboxPatch((cx - bw/2, by - bh/2), bw, bh,
                                     boxstyle="round,pad=0,rounding_size=1.6",
                                     facecolor=face, edgecolor=edge, lw=1.0, zorder=1))
        ax.text(cx, by + 9, label, ha="center", fontsize=12, color=edge, fontweight="bold")
        ax.text(cx, by + 4, sub, ha="center", fontsize=8.5, color=INK_, style="italic")

        # 4 gate indicators with subtle drop shadow
        gate_x = np.linspace(cx - 9, cx + 9, 4)
        gate_y = by - 3
        for j, gx in enumerate(gate_x):
            if gate_state == "open":
                st = "open"
            elif gate_state == "half":
                st = "open" if j < 2 else ("half" if j == 2 else "closed")
            else:
                st = "closed"
            ax.add_patch(Circle((gx + 0.15, gate_y - 0.15), 1.5, facecolor="#00000018",
                                edgecolor="none", zorder=2))
            if st == "open":
                ax.add_patch(Circle((gx, gate_y), 1.5, facecolor=edge, edgecolor="white",
                                    lw=0.6, zorder=3))
            elif st == "half":
                ax.add_patch(Circle((gx, gate_y), 1.5, facecolor="white", edgecolor=edge,
                                    lw=0.7, zorder=3))
                ax.add_patch(Rectangle((gx - 1.5, gate_y - 1.4), 1.5, 2.8, facecolor=edge,
                                        edgecolor="none", zorder=4))
            else:
                ax.add_patch(Circle((gx, gate_y), 1.5, facecolor="white", edgecolor="#aaaaaa",
                                    lw=0.7, zorder=3))
                ax.plot([gx - 0.9, gx + 0.9], [gate_y - 0.9, gate_y + 0.9], color="#888",
                        lw=0.8, zorder=4)
                ax.plot([gx - 0.9, gx + 0.9], [gate_y + 0.9, gate_y - 0.9], color="#888",
                        lw=0.8, zorder=4)
        ax.text(cx, by - 8, "gates  1   2   3   4", ha="center", fontsize=7.4, color=MUTED_)

    # Gradient axis bar — real LinearSegmentedColormap gradient
    grad_y = 28
    grad_x0, grad_x1 = 12, 88
    cmap = LinearSegmentedColormap.from_list(
        "axis", [BLUE_D, BLUE_M, PURP_L, RED_M, RED_D], N=400)
    n_seg = 400
    for i in range(n_seg):
        col = cmap(i / (n_seg - 1))
        ax.add_patch(Rectangle((grad_x0 + i * (grad_x1 - grad_x0) / n_seg, grad_y - 1.8),
                                (grad_x1 - grad_x0) / n_seg + 0.05, 3.6,
                                facecolor=col, edgecolor="none", zorder=2))
    ax.add_patch(FancyBboxPatch((grad_x0, grad_y - 1.8), grad_x1 - grad_x0, 3.6,
                                 boxstyle="round,pad=0,rounding_size=0.5",
                                 facecolor="none", edgecolor=INK_D, lw=0.7, zorder=3))
    ax.text(grad_x0 - 1.5, grad_y, "preserved", ha="right", va="center", fontsize=8.5,
            color=BLUE_D, fontweight="bold")
    ax.text(grad_x1 + 1.5, grad_y, "silenced", ha="left", va="center", fontsize=8.5,
            color=RED_D, fontweight="bold")
    ax.text(50, grad_y + 5.5, "8-gene differentiation–silencing axis",
            ha="center", fontsize=9.6, color=INK_D, fontweight="bold")
    ax.text(50, grad_y - 5,
            "SLC5A5  ·  TPO  ·  TG  ·  TSHR  ·  PAX8  ·  NKX2-1  ·  FOXE1  ·  DIO1",
            ha="center", fontsize=8, color=MUTED_, family="serif", style="italic")

    # Driver modifier badges below — secondary
    badges = ["BRAF", "RAS", "TERT", "TP53", "Fusion"]
    bx_y = 14
    bw_b = 10
    bgap = 2
    total = len(badges) * bw_b + (len(badges) - 1) * bgap
    sx = 50 - total / 2
    for i, m in enumerate(badges):
        bx = sx + i * (bw_b + bgap)
        ax.add_patch(FancyBboxPatch((bx, bx_y - 2), bw_b, 4,
                                     boxstyle="round,pad=0,rounding_size=0.8",
                                     facecolor="white", edgecolor=MUTED_, lw=0.5))
        ax.text(bx + bw_b/2, bx_y, m, ha="center", va="center", fontsize=8, color=MUTED_)
    ax.text(50, 7.5, "genomic modifiers  (secondary, not the main axis)",
            ha="center", fontsize=7.8, color=MUTED_, style="italic")


# ============================================================
# Panel d — Tiered evidence ladder
# ============================================================
def panel_d(ax):
    ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.set_facecolor(PAPER)
    for s in ax.spines.values(): s.set_visible(False)
    ax.set_xticks([]); ax.set_yticks([])
    ax.text(50, 95, "Tiered evidence ladder  ·  discovery → validation → stratification",
            ha="center", fontsize=11.5, fontweight="bold", color=INK_D)

    stages = [
        ("Discovery", "TCGA-THCA + integrated\nthyroid bulk cohorts\n(n ≈ 1,500+)",
         "Tier 4 molecular discovery", BLUE_L, BLUE_D),
        ("8-gene panel", "SLC5A5 · TPO · TG · TSHR\nPAX8 · NKX2-1 · FOXE1 · DIO1",
         "differentiation-silencing axis", PURP_L, PURP_D),
        ("Label-anchored\nvalidation", None,
         "Tier 1 / Tier 2 anchors", RED_L, RED_D),
        ("Clinical\nstratification", "RAI-avid · Gray zone\nRAI-refractory",
         "tiered patient interpretation", BEIGE_L, BEIGE_D),
    ]
    n = len(stages)
    box_w = 19; gap = 3
    total = n * box_w + (n - 1) * gap
    start_x = (100 - total) / 2
    by = 68; bh = 22
    centers = []
    for i, (title, sub, footer, face, edge) in enumerate(stages):
        bx = start_x + i * (box_w + gap)
        # shadow
        ax.add_patch(FancyBboxPatch((bx + 0.3, by - bh/2 - 0.4), box_w, bh,
                                     boxstyle="round,pad=0,rounding_size=1.4",
                                     facecolor="#00000012", edgecolor="none", zorder=0.5))
        ax.add_patch(FancyBboxPatch((bx, by - bh/2), box_w, bh,
                                     boxstyle="round,pad=0,rounding_size=1.4",
                                     facecolor=face, edgecolor=edge, lw=1.0, zorder=1))
        ax.text(bx + box_w/2, by + 7, title, ha="center", fontsize=10, color=edge, fontweight="bold")
        if sub:
            ax.text(bx + box_w/2, by, sub, ha="center", fontsize=7.7, color=INK_)
        ax.text(bx + box_w/2, by - bh/2 + 2.5, footer, ha="center", fontsize=7.2, color=MUTED_,
                style="italic")
        centers.append((bx + box_w/2, by))

    for i in range(n - 1):
        x1, y1 = centers[i]; x2, y2 = centers[i + 1]
        ax.add_patch(FancyArrowPatch((x1 + box_w/2 + 0.3, by), (x2 - box_w/2 - 0.3, by),
                                       arrowstyle="-|>", mutation_scale=13, color=MUTED_, lw=0.9))

    # Dataset tiles row beneath stage 3
    tiles = [
        ("GSE151179",   "RAI avidity",       "Tier 2"),
        ("GSE299988",   "supportive",        "Tier 2"),
        ("Boucai 2023", "RECIST exceptional","Tier 1"),
        ("Mu 2024",     "4-class uptake",    "Tier 2 · HRA004166"),
    ]
    tile_w = 16; tile_h = 11; tile_gap = 2.5
    total_tw = len(tiles) * tile_w + (len(tiles) - 1) * tile_gap
    tsx = 50 - total_tw / 2
    ty = 30
    third_cx = start_x + 2 * (box_w + gap) + box_w / 2
    ax.add_patch(FancyArrowPatch((third_cx, by - bh/2 - 0.3),
                                  (third_cx, ty + tile_h + 1.5),
                                  arrowstyle="-", mutation_scale=8, color=MUTED_, lw=0.5,
                                  linestyle=":"))
    for i, (name, kind, access) in enumerate(tiles):
        tx = tsx + i * (tile_w + tile_gap)
        ax.add_patch(FancyBboxPatch((tx + 0.2, ty - 0.3), tile_w, tile_h,
                                     boxstyle="round,pad=0,rounding_size=1.0",
                                     facecolor="#00000010", edgecolor="none", zorder=0.5))
        ax.add_patch(FancyBboxPatch((tx, ty), tile_w, tile_h,
                                     boxstyle="round,pad=0,rounding_size=1.0",
                                     facecolor="white", edgecolor=RED_D, lw=0.7))
        ax.text(tx + tile_w/2, ty + tile_h - 2.5, name, ha="center", fontsize=9, color=RED_D,
                fontweight="bold")
        ax.text(tx + tile_w/2, ty + tile_h/2 - 0.5, kind, ha="center", fontsize=7.6, color=INK_)
        ax.text(tx + tile_w/2, ty + 2.0, access, ha="center", fontsize=6.8, color=MUTED_,
                style="italic")

    ax.text(50, 18, "Tier 4 discovery → parsimonious 8-gene panel → label-anchored validation → three-class stratification",
            ha="center", fontsize=9, color=INK_D, fontweight="bold")
    ax.text(50, 13, "Risk-stratification readout, not a stand-alone clinical biomarker — prospective Tier-1 confirmation is the next step.",
            ha="center", fontsize=7.8, color=MUTED_, style="italic")


# ============================================================
def main():
    setup_rc()
    plt.rcParams["font.family"] = ["DejaVu Sans"]   # for serif italic gene names
    fig = plt.figure(figsize=(17, 14), facecolor=IVORY)
    gs = fig.add_gridspec(2, 2, hspace=0.30, wspace=0.13,
                          left=0.035, right=0.97, top=0.91, bottom=0.04)

    fig.suptitle("A multi-gate model of radioiodine failure in thyroid cancer",
                 fontsize=17, fontweight="bold", color=INK_D, y=0.965)
    fig.text(0.5, 0.927,
             "An eight-gene differentiation / iodide-handling axis linking discovery to RAI-labelled validation cohorts",
             ha="center", fontsize=11, color=MUTED_, style="italic")

    ax_a = fig.add_subplot(gs[0, 0]); panel_a(ax_a); panel_letter(ax_a, "a")
    ax_b = fig.add_subplot(gs[0, 1]); panel_b(ax_b); panel_letter(ax_b, "b")
    ax_c = fig.add_subplot(gs[1, 0]); panel_c(ax_c); panel_letter(ax_c, "c")
    ax_d = fig.add_subplot(gs[1, 1]); panel_d(ax_d); panel_letter(ax_d, "d")

    out_png = OUT_DIR / "figure1_v2_flagship.png"
    out_pdf = OUT_DIR / "figure1_v2_flagship.pdf"
    fig.savefig(out_png, dpi=220, facecolor=IVORY)
    fig.savefig(out_pdf, facecolor=IVORY)
    plt.close(fig)
    print(f"wrote {out_png}\nwrote {out_pdf}")


if __name__ == "__main__":
    main()
