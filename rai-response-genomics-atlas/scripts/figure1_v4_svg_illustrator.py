#!/usr/bin/env python3
"""Figure 1 v4 — SVG-direct Nature-illustrator grade.

Direct SVG with linearGradient / radialGradient / drop-shadow filters / Bezier
paths — converted to PNG and PDF via cairosvg. Bypasses matplotlib's
geometric-shape limitation.
"""
from __future__ import annotations
from pathlib import Path
import math
import textwrap

import cairosvg

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

# Canvas: 2000 × 1600 SVG units → high DPI raster looks crisp
W, H = 2000, 1600
PAPER = "#fbf7ee"
IVORY = "#fcf9f3"
INK   = "#1f1f24"
INK_2 = "#3a3a40"
MUTED = "#65656d"
SOFT  = "#9a9a9f"
LINE  = "#cfcdc4"

# Refined palette
BLUE_D = "#1f4d7a"
BLUE   = "#386a99"
BLUE_M = "#7ba0c5"
BLUE_L = "#c4d6e6"
BLUE_VL = "#e3ecf3"
RED_D  = "#7c322e"
RED    = "#a44d46"
RED_M  = "#cb857f"
RED_L  = "#e5beba"
RED_VL = "#f5e0dd"
PURP_D = "#5a4670"
PURP   = "#7e689a"
PURP_L = "#c2b3d2"
PURP_VL = "#ebe3f1"
BEIGE_D = "#806832"
BEIGE   = "#b89858"
BEIGE_L = "#e4ce93"
BEIGE_VL = "#f5e8c5"
AMBER_D = "#8c5a1c"
AMBER   = "#c08f33"
AMBER_L = "#e8c170"


def header():
    return f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}"
  width="{W}" height="{H}"
  font-family="Inter, 'Helvetica Neue', Helvetica, Arial, 'Liberation Sans', sans-serif">'''


def defs():
    return f'''
<defs>
  <!-- subtle drop shadow filter -->
  <filter id="drop" x="-20%" y="-20%" width="140%" height="140%">
    <feGaussianBlur in="SourceAlpha" stdDeviation="2.5"/>
    <feOffset dx="0" dy="2.5" result="offsetblur"/>
    <feComponentTransfer><feFuncA type="linear" slope="0.20"/></feComponentTransfer>
    <feMerge><feMergeNode/><feMergeNode in="SourceGraphic"/></feMerge>
  </filter>
  <filter id="drop-soft" x="-30%" y="-30%" width="160%" height="160%">
    <feGaussianBlur in="SourceAlpha" stdDeviation="4"/>
    <feOffset dx="0" dy="3" result="offsetblur"/>
    <feComponentTransfer><feFuncA type="linear" slope="0.12"/></feComponentTransfer>
    <feMerge><feMergeNode/><feMergeNode in="SourceGraphic"/></feMerge>
  </filter>
  <filter id="glow">
    <feGaussianBlur stdDeviation="5" result="b"/>
    <feMerge><feMergeNode in="b"/><feMergeNode in="SourceGraphic"/></feMerge>
  </filter>

  <!-- radial gradients for cells -->
  <radialGradient id="cellBlue" cx="0.35" cy="0.35" r="0.75">
    <stop offset="0%"   stop-color="#f0f5fb"/>
    <stop offset="55%"  stop-color="{BLUE_VL}"/>
    <stop offset="100%" stop-color="{BLUE_M}" stop-opacity="0.85"/>
  </radialGradient>
  <radialGradient id="nuc" cx="0.40" cy="0.35" r="0.70">
    <stop offset="0%"   stop-color="#f4ebd3"/>
    <stop offset="60%"  stop-color="{BEIGE_L}"/>
    <stop offset="100%" stop-color="{BEIGE}" stop-opacity="0.85"/>
  </radialGradient>
  <radialGradient id="colloid" cx="0.45" cy="0.40" r="0.70">
    <stop offset="0%"   stop-color="#fbf2d3"/>
    <stop offset="55%"  stop-color="{BEIGE_L}"/>
    <stop offset="100%" stop-color="{BEIGE_D}" stop-opacity="0.55"/>
  </radialGradient>
  <radialGradient id="iodine" cx="0.35" cy="0.35" r="0.70">
    <stop offset="0%"   stop-color="#fff0c4"/>
    <stop offset="50%"  stop-color="{AMBER_L}"/>
    <stop offset="100%" stop-color="{AMBER_D}"/>
  </radialGradient>

  <!-- capsule two-tone -->
  <linearGradient id="capsuleL" x1="0" y1="0" x2="0" y2="1">
    <stop offset="0%" stop-color="{RED_M}"/>
    <stop offset="100%" stop-color="{RED}"/>
  </linearGradient>
  <linearGradient id="capsuleR" x1="0" y1="0" x2="0" y2="1">
    <stop offset="0%" stop-color="{RED_VL}"/>
    <stop offset="100%" stop-color="{RED_L}"/>
  </linearGradient>

  <!-- vessel gradient -->
  <linearGradient id="vessel" x1="0" y1="0" x2="1" y2="0">
    <stop offset="0%" stop-color="{RED_VL}"/>
    <stop offset="50%" stop-color="{RED_L}" stop-opacity="0.85"/>
    <stop offset="100%" stop-color="{RED_VL}"/>
  </linearGradient>

  <!-- differentiation axis gradient -->
  <linearGradient id="axisBar" x1="0" y1="0" x2="1" y2="0">
    <stop offset="0%"   stop-color="{BLUE_D}"/>
    <stop offset="32%"  stop-color="{BLUE_M}"/>
    <stop offset="50%"  stop-color="{PURP_L}"/>
    <stop offset="68%"  stop-color="{RED_M}"/>
    <stop offset="100%" stop-color="{RED_D}"/>
  </linearGradient>

  <!-- decay ray gradient -->
  <linearGradient id="rayG" x1="0" y1="0" x2="1" y2="0">
    <stop offset="0%" stop-color="{RED_D}"/>
    <stop offset="100%" stop-color="{RED_D}" stop-opacity="0"/>
  </linearGradient>

  <!-- arrow markers -->
  <marker id="arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
    <path d="M0,0 L10,5 L0,10 z" fill="{MUTED}"/>
  </marker>
  <marker id="arrBlue" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">
    <path d="M0,0 L10,5 L0,10 z" fill="{BLUE_D}"/>
  </marker>
  <marker id="arrRed" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">
    <path d="M0,0 L10,5 L0,10 z" fill="{RED_D}"/>
  </marker>
  <marker id="arrInk" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
    <path d="M0,0 L10,5 L0,10 z" fill="{INK_2}"/>
  </marker>
</defs>
'''


def background():
    return f'<rect width="{W}" height="{H}" fill="{IVORY}"/>'


def title_block():
    return f'''
<g>
  <text x="{W/2}" y="55" text-anchor="middle" font-size="36" font-weight="700" fill="{INK}">
    A multi-gate model of radioiodine failure in thyroid cancer
  </text>
  <text x="{W/2}" y="92" text-anchor="middle" font-size="18" font-style="italic" fill="{MUTED}">
    An eight-gene differentiation / iodide-handling axis linking discovery to RAI-labelled validation cohorts
  </text>
</g>
'''


def panel_box(x, y, w, h, letter, title):
    return f'''
<g>
  <rect x="{x}" y="{y}" width="{w}" height="{h}" rx="14" ry="14"
        fill="{PAPER}" stroke="{LINE}" stroke-width="1.2" filter="url(#drop-soft)"/>
  <text x="{x+22}" y="{y+44}" font-size="28" font-weight="800" fill="{INK}">{letter}</text>
  <text x="{x+w/2}" y="{y+42}" text-anchor="middle" font-size="20" font-weight="700" fill="{INK}">{title}</text>
</g>
'''


# ============================================================
# Panel A — Clinical unmet need
# ============================================================
def panel_a(x0, y0, w, h):
    out = [panel_box(x0, y0, w, h, "a", "Clinical unmet need")]

    # Thyroid butterfly icon
    cx, cy = x0 + 110, y0 + 230
    out.append(f'''
<g filter="url(#drop)">
  <ellipse cx="{cx-32}" cy="{cy}" rx="38" ry="76" fill="{BEIGE_L}" stroke="{BEIGE_D}" stroke-width="1.5" transform="rotate(-12 {cx-32} {cy})"/>
  <ellipse cx="{cx+32}" cy="{cy}" rx="38" ry="76" fill="{BEIGE_L}" stroke="{BEIGE_D}" stroke-width="1.5" transform="rotate(12 {cx+32} {cy})"/>
  <rect x="{cx-22}" y="{cy-12}" width="44" height="22" rx="3" fill="{BEIGE_L}" stroke="{BEIGE_D}" stroke-width="1.2"/>
  <!-- vessel hint -->
  <path d="M {cx-44} {cy-58} Q {cx} {cy-68} {cx+44} {cy-56}" stroke="{RED_M}" stroke-width="1.6" fill="none" opacity="0.65"/>
</g>
<text x="{cx}" y="{cy+115}" text-anchor="middle" font-size="17" fill="{INK}">Differentiated thyroid cancer</text>
<text x="{cx}" y="{cy+138}" text-anchor="middle" font-size="14" fill="{MUTED}" font-style="italic">(post-surgery)</text>
''')

    # Arrow to capsule
    out.append(f'<line x1="{cx+85}" y1="{cy}" x2="{cx+185}" y2="{cy}" stroke="{MUTED}" stroke-width="2" marker-end="url(#arr)"/>')

    # I-131 capsule
    cap_x, cap_y, cap_w, cap_h = cx+195, cy-30, 130, 60
    out.append(f'''
<g filter="url(#drop)">
  <rect x="{cap_x}" y="{cap_y}" width="{cap_w/2}" height="{cap_h}" rx="30" ry="30" fill="url(#capsuleL)"/>
  <rect x="{cap_x+cap_w/2}" y="{cap_y}" width="{cap_w/2}" height="{cap_h}" rx="30" ry="30" fill="url(#capsuleR)"/>
  <line x1="{cap_x+cap_w/2}" y1="{cap_y+6}" x2="{cap_x+cap_w/2}" y2="{cap_y+cap_h-6}" stroke="{RED_D}" stroke-width="1"/>
</g>
<text x="{cap_x+cap_w/2}" y="{cap_y+cap_h/2+8}" text-anchor="middle" font-size="22" font-weight="700" fill="white">¹³¹I</text>
<text x="{cap_x+cap_w/2}" y="{cap_y+cap_h+30}" text-anchor="middle" font-size="16" fill="{INK}">Radioiodine therapy</text>
''')

    # Decay rays from capsule
    cap_cx, cap_cy = cap_x+cap_w+5, cap_y+cap_h/2
    for ang_deg in (-35, -18, 0, 18, 35):
        ang = math.radians(ang_deg)
        x1 = cap_cx + 4*math.cos(ang); y1 = cap_cy + 4*math.sin(ang)
        x2 = cap_cx + 28*math.cos(ang); y2 = cap_cy + 28*math.sin(ang)
        out.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{RED_D}" stroke-width="1.5" opacity="0.7"/>')

    # Decision diamond
    dcx, dcy = x0 + w/2 + 30, y0 + 230
    out.append(f'''
<g filter="url(#drop)">
  <polygon points="{dcx},{dcy-50} {dcx+60},{dcy} {dcx},{dcy+50} {dcx-60},{dcy}"
           fill="white" stroke="{MUTED}" stroke-width="1.8" stroke-dasharray="6,4"/>
</g>
<text x="{dcx}" y="{dcy+12}" text-anchor="middle" font-size="34" font-weight="700" fill="{MUTED}">?</text>
<text x="{dcx}" y="{dcy+90}" text-anchor="middle" font-size="14" fill="{MUTED}" font-style="italic">no robust pre-RAI predictor</text>
''')

    # Branches to outcomes
    branch_x = x0 + w - 280
    rai_avid_y = y0 + 130
    rai_ref_y = y0 + 360
    # Avid (upper)
    out.append(f'''
<path d="M {dcx+58} {dcy-15} Q {(dcx+branch_x)/2+30} {(dcy+rai_avid_y)/2-50} {branch_x-10} {rai_avid_y+30}"
      stroke="{BLUE_D}" stroke-width="2.5" fill="none" marker-end="url(#arrBlue)"/>
<g filter="url(#drop)">
  <rect x="{branch_x}" y="{rai_avid_y}" width="260" height="68" rx="10" fill="{BLUE_VL}" stroke="{BLUE_D}" stroke-width="1.4"/>
</g>
<!-- check icon -->
<path d="M {branch_x+18} {rai_avid_y+38} l 15 15 l 28 -30" stroke="{BLUE_D}" stroke-width="4" fill="none" stroke-linecap="round" stroke-linejoin="round"/>
<text x="{branch_x+78}" y="{rai_avid_y+30}" font-size="19" font-weight="700" fill="{BLUE_D}">RAI-avid</text>
<text x="{branch_x+78}" y="{rai_avid_y+54}" font-size="14" fill="{INK_2}">remission · ≥ 90% 10-yr DSS</text>
''')
    # Refractory (lower)
    out.append(f'''
<path d="M {dcx+58} {dcy+15} Q {(dcx+branch_x)/2+30} {(dcy+rai_ref_y)/2+50} {branch_x-10} {rai_ref_y+30}"
      stroke="{RED_D}" stroke-width="2.5" fill="none" marker-end="url(#arrRed)"/>
<g filter="url(#drop)">
  <rect x="{branch_x}" y="{rai_ref_y}" width="260" height="80" rx="10" fill="{RED_VL}" stroke="{RED_D}" stroke-width="1.4"/>
</g>
<!-- cross icon -->
<path d="M {branch_x+18} {rai_ref_y+22} l 26 26 M {branch_x+44} {rai_ref_y+22} l -26 26" stroke="{RED_D}" stroke-width="4" stroke-linecap="round"/>
<text x="{branch_x+78}" y="{rai_ref_y+30}" font-size="19" font-weight="700" fill="{RED_D}">RAI-refractory</text>
<text x="{branch_x+78}" y="{rai_ref_y+52}" font-size="14" fill="{INK_2}">persistent / metastatic disease</text>
<text x="{branch_x+78}" y="{rai_ref_y+70}" font-size="13" fill="{RED_D}" font-style="italic">~10–14% 10-yr DSS (distant)</text>
''')

    # Take-home banner
    bx, by, bw, bh = x0 + 30, y0 + h - 110, w - 60, 80
    out.append(f'''
<rect x="{bx}" y="{by}" width="{bw}" height="{bh}" rx="10" fill="#fdf2d8" stroke="{BEIGE_D}" stroke-width="1" opacity="0.75"/>
<text x="{x0+w/2}" y="{by+33}" text-anchor="middle" font-size="17" font-weight="700" fill="{INK}">
  Need: a parsimonious molecular readout to triage RAI candidates <tspan font-style="italic">before</tspan> therapy
</text>
<text x="{x0+w/2}" y="{by+60}" text-anchor="middle" font-size="14" fill="{MUTED}" font-style="italic">
  Current biomarkers (BRAF · RAS · TERT) do not reliably partition the RAI-refractory population.
</text>
''')

    return "\n".join(out)


# ============================================================
# Panel B — Multi-gate biology with anatomical follicle
# ============================================================
def panel_b(x0, y0, w, h):
    out = [panel_box(x0, y0, w, h, "b", "Multi-gate biology of RAI response")]

    # Capillary on left
    cap_x = x0 + 65
    out.append(f'''
<path d="M {cap_x-32} {y0+150} Q {cap_x-40} {y0+h/2} {cap_x-30} {y0+h-150}
         L {cap_x+30} {y0+h-150} Q {cap_x+40} {y0+h/2} {cap_x+32} {y0+150} Z"
      fill="url(#vessel)" stroke="{RED_D}" stroke-width="1.2" opacity="0.85" filter="url(#drop)"/>
<text x="{cap_x}" y="{y0+135}" text-anchor="middle" font-size="15" font-weight="700" fill="{RED_D}">blood vessel</text>
''')
    # Iodine particles in vessel
    import random
    rng = random.Random(13)
    for y in (y0+220, y0+290, y0+360, y0+450, y0+540):
        dx = rng.uniform(-15, 15)
        out.append(f'<circle cx="{cap_x+dx}" cy="{y}" r="14" fill="url(#iodine)" filter="url(#glow)"/>')
        out.append(f'<text x="{cap_x+dx}" y="{y+5}" text-anchor="middle" font-size="12" font-weight="700" fill="{AMBER_D}">I⁻</text>')

    # Three thyrocyte cells
    cell_x = x0 + 215
    for i, cy in enumerate((y0+225, y0+395, y0+565)):
        # cell
        out.append(f'''
<g filter="url(#drop)">
  <ellipse cx="{cell_x}" cy="{cy}" rx="78" ry="92" fill="url(#cellBlue)" stroke="{BLUE_D}" stroke-width="1.3"/>
  <ellipse cx="{cell_x-8}" cy="{cy-6}" rx="36" ry="44" fill="url(#nuc)" stroke="{BEIGE_D}" stroke-width="1.2"/>
</g>
''')
        # Lineage TF dots in nucleus
        for j in range(11):
            ang = (j / 11) * 2*math.pi
            r_d = rng.uniform(8, 28)
            tx = cell_x - 8 + r_d*math.cos(ang); ty = cy - 6 + r_d*math.sin(ang)*0.9
            out.append(f'<circle cx="{tx:.1f}" cy="{ty:.1f}" r="3.0" fill="{BLUE_D}" opacity="0.75"/>')

        # NIS transporter on basal side
        nx, ny = cell_x-74, cy
        out.append(f'''
<g>
  <rect x="{nx-10}" y="{ny-22}" width="8" height="44" rx="3" fill="{RED_D}"/>
  <rect x="{nx+2}" y="{ny-22}" width="8" height="44" rx="3" fill="{RED_D}"/>
  <rect x="{nx-2}" y="{ny-20}" width="4" height="40" rx="2" fill="{RED_VL}"/>
</g>
''')
        # iodine flowing through
        out.append(f'<circle cx="{nx-22}" cy="{ny}" r="10" fill="url(#iodine)" opacity="0.9"/>')
        out.append(f'<circle cx="{nx+22}" cy="{ny}" r="9" fill="url(#iodine)" opacity="0.85"/>')
        # to colloid
        out.append(f'<circle cx="{cell_x+78}" cy="{cy+0}" r="9" fill="url(#iodine)" opacity="0.85"/>')

    # Follicle / colloid lumen on right
    fcx, fcy = x0 + 480, y0 + 420
    fw, fh = 250, 360
    out.append(f'''
<g filter="url(#drop)">
  <ellipse cx="{fcx}" cy="{fcy}" rx="{fw/2}" ry="{fh/2}" fill="url(#colloid)" stroke="{BEIGE_D}" stroke-width="1.6"/>
</g>
<text x="{fcx}" y="{fcy-fh/2-20}" text-anchor="middle" font-size="16" font-weight="700" fill="{BEIGE_D}">colloid lumen</text>
''')
    # TG storage stippling
    for _ in range(46):
        rx_ = rng.uniform(-fw*0.32, fw*0.32)
        ry_ = rng.uniform(-fh*0.34, fh*0.34)
        if (rx_/(fw*0.42))**2 + (ry_/(fh*0.42))**2 < 0.95:
            r = rng.uniform(2.5, 5)
            out.append(f'<circle cx="{fcx+rx_:.1f}" cy="{fcy+ry_:.1f}" r="{r:.1f}" fill="{BEIGE_D}" opacity="{rng.uniform(0.4,0.75):.2f}"/>')

    # Microvilli on apical side (facing follicle) - top arc
    n_micro = 12
    for k in range(n_micro):
        ang = -math.pi/2 + (-0.5 + k/(n_micro-1)) * 1.0
        x_o = fcx + (fw/2+8)*math.cos(ang+math.pi); y_o = fcy + (fh/2)*math.sin(ang+math.pi)

    # Decay rays — radiation killing
    rcx, rcy = fcx+20, fcy
    for i in range(14):
        ang = i * (2*math.pi/14)
        x1 = rcx + 12*math.cos(ang); y1 = rcy + 12*math.sin(ang)
        x2 = rcx + 36*math.cos(ang); y2 = rcy + 36*math.sin(ang)
        out.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{RED_D}" stroke-width="2" opacity="0.7" stroke-linecap="round"/>')
    out.append(f'<circle cx="{rcx}" cy="{rcy}" r="8" fill="{RED_D}"/>')

    # Gate annotations
    gates = [
        (1, "SLC5A5 / NIS", "iodide uptake", x0+85, y0+700, cell_x-50, cy_target := y0+260, RED_D),
        (2, "PAX8 · NKX2-1 · FOXE1 · TSHR", "lineage program", x0+95, y0+775, cell_x+15, y0+395, BLUE_D),
        (3, "TG · TPO · DIO1", "organification &amp; storage", x0+800, y0+700, fcx-30, fcy-30, BLUE_D),
        (4, "retention →", "radiation killing", x0+800, y0+775, rcx+20, rcy+10, RED_D),
    ]
    for gn, title, sub, lx, ly, tx, ty, c in gates:
        out.append(f'''
<g>
  <rect x="{lx-110}" y="{ly-30}" width="220" height="56" rx="8" fill="white" stroke="{c}" stroke-width="1.2"/>
  <text x="{lx}" y="{ly-8}" text-anchor="middle" font-size="15" font-weight="700" fill="{c}">Gate {gn} · {title}</text>
  <text x="{lx}" y="{ly+12}" text-anchor="middle" font-size="13" fill="{INK_2}">{sub}</text>
</g>
<path d="M {lx} {ly-30} Q {(lx+tx)/2} {(ly+ty)/2-30} {tx} {ty}" stroke="{c}" stroke-width="0.9" fill="none" opacity="0.55"/>
''')

    # take-home banner
    bx, by, bw, bh = x0 + 30, y0 + h - 80, w - 60, 60
    out.append(f'''
<rect x="{bx}" y="{by}" width="{bw}" height="{bh}" rx="10" fill="{RED_VL}" stroke="{RED_D}" stroke-width="1" opacity="0.7"/>
<text x="{x0+w/2}" y="{by+38}" text-anchor="middle" font-size="16" font-weight="700" font-style="italic" fill="{INK}">
  NIS alone is insufficient — coordinated differentiation across all four gates is required.
</text>
''')

    return "\n".join(out)


# ============================================================
# Panel C — Three molecular states + gradient axis
# ============================================================
def panel_c(x0, y0, w, h):
    out = [panel_box(x0, y0, w, h, "c", "Three molecular states along the eight-gene differentiation–silencing axis")]

    states = [
        ("RAI-avid",      "preserved differentiation\nall gates open · high score",   BLUE_VL, BLUE_D,  "open"),
        ("Gray zone",     "partial lineage loss\nintermediate score",                  PURP_VL, PURP_D,  "half"),
        ("RAI-refractory","differentiation silencing\ncollapsed gates · low score",    RED_VL,  RED_D,   "closed"),
    ]
    bw_s, bh_s = 280, 200
    gap = 50
    total = 3*bw_s + 2*gap
    sx = x0 + (w - total)/2
    by_s = y0 + 140
    centers = []
    for i, (label, sub, face, edge, gate_state) in enumerate(states):
        bx = sx + i*(bw_s + gap)
        centers.append((bx + bw_s/2, by_s + bh_s/2))
        out.append(f'''
<g filter="url(#drop)">
  <rect x="{bx}" y="{by_s}" width="{bw_s}" height="{bh_s}" rx="14" fill="{face}" stroke="{edge}" stroke-width="1.6"/>
</g>
<text x="{bx+bw_s/2}" y="{by_s+45}" text-anchor="middle" font-size="22" font-weight="700" fill="{edge}">{label}</text>
''')
        # Multi-line subtitle
        for j, line in enumerate(sub.split("\n")):
            out.append(f'<text x="{bx+bw_s/2}" y="{by_s+85+j*22}" text-anchor="middle" font-size="14" fill="{INK_2}" font-style="italic">{line}</text>')

        # 4 gate dots
        gx_arr = [bx + bw_s/2 + (k-1.5)*40 for k in range(4)]
        gy = by_s + bh_s - 40
        for k, gx in enumerate(gx_arr):
            if gate_state == "open":
                st = "open"
            elif gate_state == "half":
                st = "open" if k < 2 else ("half" if k == 2 else "closed")
            else:
                st = "closed"
            if st == "open":
                out.append(f'<circle cx="{gx}" cy="{gy}" r="11" fill="{edge}" stroke="white" stroke-width="1.5"/>')
            elif st == "half":
                out.append(f'<circle cx="{gx}" cy="{gy}" r="11" fill="white" stroke="{edge}" stroke-width="1.5"/>')
                out.append(f'<path d="M {gx-11} {gy} A 11 11 0 0 1 {gx+11} {gy} L {gx} {gy} Z" fill="{edge}"/>')
            else:
                out.append(f'<circle cx="{gx}" cy="{gy}" r="11" fill="white" stroke="{SOFT}" stroke-width="1.5"/>')
                out.append(f'<path d="M {gx-7} {gy-7} L {gx+7} {gy+7} M {gx-7} {gy+7} L {gx+7} {gy-7}" stroke="{SOFT}" stroke-width="2"/>')
        out.append(f'<text x="{bx+bw_s/2}" y="{by_s+bh_s-12}" text-anchor="middle" font-size="12" fill="{MUTED}">gates  1   2   3   4</text>')

    # Gradient axis bar
    ax0 = x0 + 90; ax1 = x0 + w - 90
    ay = y0 + h - 200
    out.append(f'''
<rect x="{ax0}" y="{ay-22}" width="{ax1-ax0}" height="44" rx="10" fill="url(#axisBar)" stroke="{INK}" stroke-width="1.2" filter="url(#drop)"/>
<text x="{ax0-15}" y="{ay+8}" text-anchor="end" font-size="16" font-weight="700" fill="{BLUE_D}">preserved</text>
<text x="{ax1+15}" y="{ay+8}" text-anchor="start" font-size="16" font-weight="700" fill="{RED_D}">silenced</text>
<text x="{x0+w/2}" y="{ay-50}" text-anchor="middle" font-size="17" font-weight="700" fill="{INK}">8-gene differentiation–silencing axis</text>
<text x="{x0+w/2}" y="{ay+55}" text-anchor="middle" font-size="14" font-style="italic" fill="{MUTED}" font-family="'Liberation Serif',serif">
  SLC5A5 · TPO · TG · TSHR · PAX8 · NKX2-1 · FOXE1 · DIO1
</text>
''')

    # Driver modifier badges
    badges = ["BRAF", "RAS", "TERT", "TP53", "Fusion"]
    bw_b, bh_b = 90, 32
    bgap = 14
    total_b = len(badges)*bw_b + (len(badges)-1)*bgap
    sx_b = x0 + (w - total_b)/2
    bby = y0 + h - 75
    for i, m in enumerate(badges):
        bx = sx_b + i*(bw_b + bgap)
        out.append(f'<rect x="{bx}" y="{bby}" width="{bw_b}" height="{bh_b}" rx="8" fill="white" stroke="{MUTED}" stroke-width="1"/>')
        out.append(f'<text x="{bx+bw_b/2}" y="{bby+22}" text-anchor="middle" font-size="14" fill="{MUTED}">{m}</text>')
    out.append(f'<text x="{x0+w/2}" y="{bby+58}" text-anchor="middle" font-size="13" font-style="italic" fill="{MUTED}">genomic modifiers · secondary, not the main axis</text>')

    return "\n".join(out)


# ============================================================
# Panel D — Tiered evidence ladder
# ============================================================
def panel_d(x0, y0, w, h):
    out = [panel_box(x0, y0, w, h, "d", "Tiered evidence ladder · discovery → validation → stratification")]

    stages = [
        ("Discovery",     "TCGA-THCA + integrated\nthyroid bulk cohorts\n(n ≈ 1,500+)", "Tier 4 molecular discovery", BLUE_VL, BLUE_D),
        ("8-gene panel",  "SLC5A5 · TPO · TG · TSHR\nPAX8 · NKX2-1 · FOXE1 · DIO1",      "differentiation–silencing axis",  PURP_VL, PURP_D),
        ("Validation",    "label-anchored\nexternal cohorts",                            "Tier 1 / Tier 2 anchors", RED_VL,  RED_D),
        ("Stratification","RAI-avid · Gray zone\nRAI-refractory",                        "tiered patient interpretation",   BEIGE_VL, BEIGE_D),
    ]
    bw_s, bh_s = 230, 200
    gap_s = 32
    total = len(stages)*bw_s + (len(stages)-1)*gap_s
    sx = x0 + (w - total)/2
    by_s = y0 + 110
    centers = []
    for i, (title, sub, foot, face, edge) in enumerate(stages):
        bx = sx + i*(bw_s + gap_s)
        centers.append((bx + bw_s/2, by_s + bh_s/2))
        out.append(f'''
<g filter="url(#drop)">
  <rect x="{bx}" y="{by_s}" width="{bw_s}" height="{bh_s}" rx="14" fill="{face}" stroke="{edge}" stroke-width="1.5"/>
</g>
<text x="{bx+bw_s/2}" y="{by_s+38}" text-anchor="middle" font-size="18" font-weight="700" fill="{edge}">{title}</text>
''')
        for j, line in enumerate(sub.split("\n")):
            out.append(f'<text x="{bx+bw_s/2}" y="{by_s+75+j*22}" text-anchor="middle" font-size="14" fill="{INK_2}">{line}</text>')
        out.append(f'<text x="{bx+bw_s/2}" y="{by_s+bh_s-15}" text-anchor="middle" font-size="12" fill="{MUTED}" font-style="italic">{foot}</text>')

    # Arrows
    for i in range(len(centers)-1):
        x1 = centers[i][0] + bw_s/2 + 2
        x2 = centers[i+1][0] - bw_s/2 - 2
        out.append(f'<line x1="{x1}" y1="{by_s+bh_s/2}" x2="{x2}" y2="{by_s+bh_s/2}" stroke="{INK_2}" stroke-width="2" marker-end="url(#arrInk)"/>')

    # Dataset tiles under Validation column
    tiles = [
        ("GSE151179",   "RAI avidity",       "Tier 2"),
        ("GSE299988",   "supportive",        "Tier 2"),
        ("Boucai 2023", "RECIST exceptional","Tier 1"),
        ("Mu 2024",     "4-class uptake",    "Tier 2 · HRA004166"),
    ]
    tw_t, th_t = 175, 100
    gap_t = 18
    total_t = len(tiles)*tw_t + (len(tiles)-1)*gap_t
    tsx = x0 + (w - total_t)/2
    ty = by_s + bh_s + 90
    third_cx = centers[2][0]
    out.append(f'<line x1="{third_cx}" y1="{by_s+bh_s+5}" x2="{third_cx}" y2="{ty-15}" stroke="{MUTED}" stroke-width="1" stroke-dasharray="3,3"/>')

    for i, (name, kind, access) in enumerate(tiles):
        tx = tsx + i*(tw_t + gap_t)
        out.append(f'''
<g filter="url(#drop)">
  <rect x="{tx}" y="{ty}" width="{tw_t}" height="{th_t}" rx="10" fill="white" stroke="{RED_D}" stroke-width="1.2"/>
</g>
<text x="{tx+tw_t/2}" y="{ty+32}" text-anchor="middle" font-size="17" font-weight="700" fill="{RED_D}">{name}</text>
<text x="{tx+tw_t/2}" y="{ty+58}" text-anchor="middle" font-size="14" fill="{INK_2}">{kind}</text>
<text x="{tx+tw_t/2}" y="{ty+82}" text-anchor="middle" font-size="12" fill="{MUTED}" font-style="italic">{access}</text>
''')

    # bottom take-home
    by_h = y0 + h - 68
    out.append(f'''
<text x="{x0+w/2}" y="{by_h}" text-anchor="middle" font-size="15" font-weight="700" fill="{INK}">
  Tier-4 discovery → parsimonious 8-gene panel → label-anchored validation → three-class stratification
</text>
<text x="{x0+w/2}" y="{by_h+24}" text-anchor="middle" font-size="13" fill="{MUTED}" font-style="italic">
  Positioned as a risk-stratification readout — prospective Tier-1 confirmation is the next step.
</text>
''')
    return "\n".join(out)


def main():
    pw, ph = (W - 60)/2, (H - 180)/2
    gap = 20
    x_left, x_right = 30, 30 + pw + gap
    y_top, y_bot = 130, 130 + ph + gap
    # Recompute so total fits
    pw = (W - 30*2 - gap)/2
    ph = (H - 130 - 30 - gap)/2
    x_left, x_right = 30, 30 + pw + gap
    y_top, y_bot = 130, 130 + ph + gap

    svg = "\n".join([
        header(),
        defs(),
        background(),
        title_block(),
        panel_a(x_left,  y_top, pw, ph),
        panel_b(x_right, y_top, pw, ph),
        panel_c(x_left,  y_bot, pw, ph),
        panel_d(x_right, y_bot, pw, ph),
        "</svg>",
    ])
    svg_path = OUT / "figure1_v3_flagship.svg"
    svg_path.write_text(svg)

    png_path = OUT / "figure1_v3_flagship.png"
    pdf_path = OUT / "figure1_v3_flagship.pdf"
    cairosvg.svg2png(bytestring=svg.encode("utf-8"), write_to=str(png_path),
                     output_width=W*1.6, output_height=H*1.6)
    cairosvg.svg2pdf(bytestring=svg.encode("utf-8"), write_to=str(pdf_path))

    print(f"wrote {svg_path}\nwrote {png_path}\nwrote {pdf_path}")


if __name__ == "__main__":
    main()
