#!/usr/bin/env python3
"""
Paper 2 Fig 1 v2 — Clinical triage schematic (cv2-polished, 4K).

Upgrades over make_fig1.py:
  - Smooth gradient backgrounds for box fills (cv2.addWeighted of solid + radial fade).
  - cv2.arrowedLine + LINE_AA strokes.
  - cv2.putText with FONT_HERSHEY_DUPLEX/TRIPLEX (anti-aliased).
  - Per-step cost annotations.
  - Branded colour palette:
      Patient/H&E:     #006B7D deep teal
      Model:           #6A4C93 purple
      Decision:        #FFB627 amber
      DM1+ branch:     #2A9D8F green
      DM2 branch:      #5B85AA blue
      Standard care:   #6C757D gray
  - 4K (3840 x ~2160) output.

Outputs:
  figures/fig1_triage_schematic_v2.png
  figures/fig1_triage_schematic_v2.pdf
"""
from __future__ import annotations

from pathlib import Path

import cv2
import numpy as np

import matplotlib.pyplot as plt

OUT_DIR = Path(
    "/home/seungho/personal/THCA_data_analysis/project/results/"
    "p2_image_dm1_v2_foundation_clam_2026_05_07/figures"
)
OUT_DIR.mkdir(parents=True, exist_ok=True)

PNG_OUT = OUT_DIR / "fig1_triage_schematic_v2.png"
PDF_OUT = OUT_DIR / "fig1_triage_schematic_v2.pdf"

# --------------------------------------------------------------------------- palette (BGR)
def hex_to_bgr(hex_str: str) -> tuple[int, int, int]:
    h = hex_str.lstrip("#")
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    return (b, g, r)


C_TEAL = hex_to_bgr("#006B7D")
C_PURPLE = hex_to_bgr("#6A4C93")
C_AMBER = hex_to_bgr("#FFB627")
C_GREEN = hex_to_bgr("#2A9D8F")
C_BLUE = hex_to_bgr("#5B85AA")
C_GRAY = hex_to_bgr("#6C757D")
C_LIGHT = (245, 245, 248)
C_DARK = (40, 40, 40)
C_MID = (90, 90, 90)
WHITE = (255, 255, 255)


# --------------------------------------------------------------------------- helpers
def gradient_box(canvas: np.ndarray, x0: int, y0: int, x1: int, y1: int,
                 base_bgr: tuple, *, lighten: float = 0.18,
                 corner_r: int = 16, border_color: tuple | None = None,
                 border_px: int = 2):
    """Draw a rounded rectangle with a subtle vertical gradient (lighter at top).
    Lightens base_bgr at top by `lighten` (0-1) toward white."""
    h = y1 - y0
    w = x1 - x0
    box = np.empty((h, w, 3), dtype=np.uint8)
    base = np.array(base_bgr, dtype=np.float32)
    top = base + (np.array(WHITE, dtype=np.float32) - base) * lighten
    bot = base
    for i in range(h):
        t = i / max(1, h - 1)
        col = top * (1 - t) + bot * t
        box[i, :, :] = col.astype(np.uint8)

    # Build a rounded-rect mask
    mask = np.zeros((h, w), dtype=np.uint8)
    cv2.rectangle(mask, (corner_r, 0), (w - corner_r, h), 255, -1)
    cv2.rectangle(mask, (0, corner_r), (w, h - corner_r), 255, -1)
    cv2.circle(mask, (corner_r, corner_r), corner_r, 255, -1, cv2.LINE_AA)
    cv2.circle(mask, (w - corner_r, corner_r), corner_r, 255, -1, cv2.LINE_AA)
    cv2.circle(mask, (corner_r, h - corner_r), corner_r, 255, -1, cv2.LINE_AA)
    cv2.circle(mask, (w - corner_r, h - corner_r), corner_r, 255, -1, cv2.LINE_AA)

    region = canvas[y0:y1, x0:x1]
    m3 = np.dstack([mask] * 3) / 255.0
    region[:] = (region * (1 - m3) + box * m3).astype(np.uint8)
    canvas[y0:y1, x0:x1] = region

    # outline (rounded rect approximation)
    if border_color is None:
        border_color = base_bgr
    # draw rounded outline using polylines on the same arcs
    arc_pts = []
    # top edge
    arc_pts.append([(x0 + corner_r, y0), (x1 - corner_r, y0)])
    arc_pts.append([(x1, y0 + corner_r), (x1, y1 - corner_r)])
    arc_pts.append([(x0 + corner_r, y1), (x1 - corner_r, y1)])
    arc_pts.append([(x0, y0 + corner_r), (x0, y1 - corner_r)])
    for (p0, p1) in arc_pts:
        cv2.line(canvas, p0, p1, border_color, border_px, cv2.LINE_AA)
    # corners
    cv2.ellipse(canvas, (x0 + corner_r, y0 + corner_r), (corner_r, corner_r),
                180, 0, 90, border_color, border_px, cv2.LINE_AA)
    cv2.ellipse(canvas, (x1 - corner_r, y0 + corner_r), (corner_r, corner_r),
                270, 0, 90, border_color, border_px, cv2.LINE_AA)
    cv2.ellipse(canvas, (x1 - corner_r, y1 - corner_r), (corner_r, corner_r),
                0, 0, 90, border_color, border_px, cv2.LINE_AA)
    cv2.ellipse(canvas, (x0 + corner_r, y1 - corner_r), (corner_r, corner_r),
                90, 0, 90, border_color, border_px, cv2.LINE_AA)


def put_text_centered(canvas, text: str, cx: int, cy: int,
                      font=cv2.FONT_HERSHEY_DUPLEX,
                      scale: float = 0.9, color=C_DARK, thickness: int = 2):
    """Anti-aliased centered text."""
    lines = text.split("\n")
    # measure
    sizes = [cv2.getTextSize(ln, font, scale, thickness)[0] for ln in lines]
    line_h = max(s[1] for s in sizes) + 10
    total_h = line_h * len(lines)
    y0 = cy - total_h // 2 + line_h - 4
    for ln, (tw, th) in zip(lines, sizes):
        x = cx - tw // 2
        cv2.putText(canvas, ln, (x, y0), font, scale, color, thickness, cv2.LINE_AA)
        y0 += line_h


def diamond(canvas, cx: int, cy: int, hw: int, hh: int,
            label: str, fill_bgr: tuple, text_color=C_DARK):
    pts = np.array([[cx, cy - hh], [cx + hw, cy],
                    [cx, cy + hh], [cx - hw, cy]], dtype=np.int32)
    # gradient fill: just lighten center by drawing a slightly lighter inner diamond
    # solid first
    cv2.fillPoly(canvas, [pts], fill_bgr, lineType=cv2.LINE_AA)
    inner = (pts - np.array([cx, cy])) * 0.7 + np.array([cx, cy])
    inner = inner.astype(np.int32)
    light = tuple(int(c + (255 - c) * 0.18) for c in fill_bgr)
    cv2.fillPoly(canvas, [inner], light, lineType=cv2.LINE_AA)
    # outline
    cv2.polylines(canvas, [pts], True, fill_bgr, 3, cv2.LINE_AA)
    # text
    put_text_centered(canvas, label, cx, cy, font=cv2.FONT_HERSHEY_TRIPLEX,
                      scale=0.85, color=text_color, thickness=2)


def arrow(canvas, p0: tuple, p1: tuple, color=C_DARK, thickness: int = 4,
          tip_length: float = 0.04, label: str = "", label_offset=(0, -28),
          label_color=None):
    cv2.arrowedLine(canvas, p0, p1, color, thickness, cv2.LINE_AA,
                    tipLength=tip_length)
    if label:
        mx = (p0[0] + p1[0]) // 2 + label_offset[0]
        my = (p0[1] + p1[1]) // 2 + label_offset[1]
        # white pill behind
        (tw, th), _ = cv2.getTextSize(label, cv2.FONT_HERSHEY_DUPLEX, 0.7, 2)
        cv2.rectangle(canvas, (mx - tw // 2 - 12, my - th - 8),
                      (mx + tw // 2 + 12, my + 8), WHITE, -1, cv2.LINE_AA)
        cv2.rectangle(canvas, (mx - tw // 2 - 12, my - th - 8),
                      (mx + tw // 2 + 12, my + 8),
                      label_color or color, 1, cv2.LINE_AA)
        cv2.putText(canvas, label, (mx - tw // 2, my),
                    cv2.FONT_HERSHEY_DUPLEX, 0.7,
                    label_color or color, 2, cv2.LINE_AA)


# --------------------------------------------------------------------------- main
def main():
    W, H = 3840, 2160
    # off-white background
    canvas = np.full((H, W, 3), 250, dtype=np.uint8)

    # ------ title bar (teal -> purple gradient) ------
    bar_h = 160
    for x in range(W):
        t = x / W
        b = int(C_TEAL[0] * (1 - t) + C_PURPLE[0] * t)
        g = int(C_TEAL[1] * (1 - t) + C_PURPLE[1] * t)
        r = int(C_TEAL[2] * (1 - t) + C_PURPLE[2] * t)
        canvas[0:bar_h, x] = (b, g, r)
    cv2.putText(canvas,
                "Fig 1  -  H&E -> DM1 image classifier in BRAF/RAS-negative PTC triage",
                (60, 80), cv2.FONT_HERSHEY_TRIPLEX, 1.7, WHITE, 4, cv2.LINE_AA)
    cv2.putText(canvas,
                "Reflex panel decision support  |  Paper 2 (Cell Reports Medicine target)",
                (60, 130), cv2.FONT_HERSHEY_DUPLEX, 1.0,
                (235, 235, 235), 2, cv2.LINE_AA)

    # ============================================================
    # Layout — coordinates in image pixels (W=3840, H=2160)
    # ============================================================
    # Row 1 (y ~ 350): Patient -> H&E -> UNI features -> CLAM
    # Row 2 (y ~ 750): Decision diamond
    # Row 3 (y ~ 1150): DM1+ branch | DM2 branch  (split to L/R)
    # Row 4 (y ~ 1550): Reflex 8-gene panel  /  Standard surveillance
    # Row 5 (y ~ 1900): Outcome arrows

    BOX_H = 180
    SMALL_H = 130

    # ----- ROW 1 -----
    boxes_row1 = [
        (160, 320, 720, 320 + BOX_H, "BRAF/RAS-neg PTC patient\n(presentation)", C_TEAL),
        (820, 320, 1380, 320 + BOX_H, "H&E whole-slide image\n(routine pathology)", C_TEAL),
        (1480, 320, 2040, 320 + BOX_H, "UNI foundation features\n(1024-dim x 200 tiles)", C_PURPLE),
        (2140, 320, 2700, 320 + BOX_H, "CLAM gated-attention MIL\nDM1 vs DM2 classifier", C_PURPLE),
    ]
    for (x0, y0, x1, y1, lbl, col) in boxes_row1:
        gradient_box(canvas, x0, y0, x1, y1, col, border_color=col, border_px=3,
                     corner_r=20)
        put_text_centered(canvas, lbl, (x0 + x1) // 2, (y0 + y1) // 2,
                          font=cv2.FONT_HERSHEY_DUPLEX, scale=0.95,
                          color=WHITE, thickness=2)

    # arrows row 1
    arrow(canvas, (720, 320 + BOX_H // 2), (820, 320 + BOX_H // 2),
          C_DARK, 5, 0.25, "$0", (0, -22), C_GRAY)
    arrow(canvas, (1380, 320 + BOX_H // 2), (1480, 320 + BOX_H // 2),
          C_DARK, 5, 0.25, "GPU compute", (0, -22), C_PURPLE)
    arrow(canvas, (2040, 320 + BOX_H // 2), (2140, 320 + BOX_H // 2),
          C_DARK, 5, 0.25, "<1 sec", (0, -22), C_PURPLE)

    # right end of row 1: cost callout
    callout_x0, callout_y0, callout_x1, callout_y1 = 2780, 320, 3680, 320 + BOX_H
    gradient_box(canvas, callout_x0, callout_y0, callout_x1, callout_y1,
                 (250, 250, 250), border_color=C_GRAY, border_px=2, corner_r=14,
                 lighten=0.0)
    put_text_centered(canvas,
                      "Marginal cost\nper slide: ~$0\n(infra-only)",
                      (callout_x0 + callout_x1) // 2,
                      (callout_y0 + callout_y1) // 2,
                      font=cv2.FONT_HERSHEY_DUPLEX, scale=0.85,
                      color=C_DARK, thickness=2)

    # ----- ROW 2: decision diamond -----
    diam_cx, diam_cy = 2420, 820
    diamond(canvas, diam_cx, diam_cy, 360, 200,
            "DM1 prob >= 0.5 ?", C_AMBER, text_color=C_DARK)
    # arrow CLAM -> diamond
    arrow(canvas, (2420, 320 + BOX_H), (2420, 820 - 200),
          C_DARK, 5, 0.04)

    # ----- ROW 3: branches -----
    # DM1+ (left)
    dm1_cx0, dm1_cy0 = 1250, 1180
    gradient_box(canvas, dm1_cx0, dm1_cy0, dm1_cx0 + 920, dm1_cy0 + BOX_H,
                 C_GREEN, border_color=C_GREEN, border_px=3, corner_r=20)
    put_text_centered(canvas,
                      "DM1+ predicted\n(driver-negative dark matter, suspected)",
                      dm1_cx0 + 460, dm1_cy0 + BOX_H // 2,
                      cv2.FONT_HERSHEY_DUPLEX, 0.95, WHITE, 2)
    # DM2 (right)
    dm2_cx0, dm2_cy0 = 2620, 1180
    gradient_box(canvas, dm2_cx0, dm2_cy0, dm2_cx0 + 920, dm2_cy0 + BOX_H,
                 C_BLUE, border_color=C_BLUE, border_px=3, corner_r=20)
    put_text_centered(canvas,
                      "DM2 predicted\n(non-DM1 / RAI-likely-responsive)",
                      dm2_cx0 + 460, dm2_cy0 + BOX_H // 2,
                      cv2.FONT_HERSHEY_DUPLEX, 0.95, WHITE, 2)

    # diamond -> branches arrows with YES / NO labels
    arrow(canvas, (diam_cx - 200, diam_cy + 100),
          (dm1_cx0 + 460, dm1_cy0), C_GREEN, 5, 0.025,
          "YES (DM1+)", (-30, -22), C_GREEN)
    arrow(canvas, (diam_cx + 200, diam_cy + 100),
          (dm2_cx0 + 460, dm2_cy0), C_BLUE, 5, 0.025,
          "NO (DM2)", (30, -22), C_BLUE)

    # ----- ROW 4: actions -----
    # DM1+ -> reflex 8-gene panel
    act1_x0, act1_y0 = 1250, 1620
    gradient_box(canvas, act1_x0, act1_y0, act1_x0 + 920, act1_y0 + BOX_H,
                 (250, 250, 250), border_color=C_GREEN, border_px=3,
                 corner_r=20, lighten=0.0)
    put_text_centered(canvas,
                      "Reflex: 8-gene RNA panel\n(confirm DM1, ~$300)",
                      act1_x0 + 460, act1_y0 + BOX_H // 2,
                      cv2.FONT_HERSHEY_DUPLEX, 0.95, C_GREEN, 2)
    arrow(canvas, (dm1_cx0 + 460, dm1_cy0 + BOX_H),
          (act1_x0 + 460, act1_y0), C_GREEN, 5, 0.04,
          "reflex", (0, -22), C_GREEN)

    # DM2 -> standard surveillance
    act2_x0, act2_y0 = 2620, 1620
    gradient_box(canvas, act2_x0, act2_y0, act2_x0 + 920, act2_y0 + BOX_H,
                 (250, 250, 250), border_color=C_BLUE, border_px=3,
                 corner_r=20, lighten=0.0)
    put_text_centered(canvas,
                      "Standard care\n(no reflex panel; surveillance)",
                      act2_x0 + 460, act2_y0 + BOX_H // 2,
                      cv2.FONT_HERSHEY_DUPLEX, 0.95, C_GRAY, 2)
    arrow(canvas, (dm2_cx0 + 460, dm2_cy0 + BOX_H),
          (act2_x0 + 460, act2_y0), C_BLUE, 5, 0.04,
          "spare", (0, -22), C_BLUE)

    # ----- ROW 5: outcomes -----
    out1_x0, out1_y0 = 1250, 1920
    gradient_box(canvas, out1_x0, out1_y0, out1_x0 + 920, out1_y0 + SMALL_H,
                 C_GREEN, border_color=C_GREEN, border_px=2, corner_r=14,
                 lighten=0.05)
    put_text_centered(canvas,
                      "RAI-refractory monitoring  |  ICI candidacy review",
                      out1_x0 + 460, out1_y0 + SMALL_H // 2,
                      cv2.FONT_HERSHEY_DUPLEX, 0.85, WHITE, 2)
    arrow(canvas, (act1_x0 + 460, act1_y0 + BOX_H),
          (out1_x0 + 460, out1_y0), C_GREEN, 4, 0.05)

    out2_x0, out2_y0 = 2620, 1920
    gradient_box(canvas, out2_x0, out2_y0, out2_x0 + 920, out2_y0 + SMALL_H,
                 C_GRAY, border_color=C_GRAY, border_px=2, corner_r=14,
                 lighten=0.05)
    put_text_centered(canvas,
                      "Routine post-thyroidectomy follow-up",
                      out2_x0 + 460, out2_y0 + SMALL_H // 2,
                      cv2.FONT_HERSHEY_DUPLEX, 0.85, WHITE, 2)
    arrow(canvas, (act2_x0 + 460, act2_y0 + BOX_H),
          (out2_x0 + 460, out2_y0), C_GRAY, 4, 0.05)

    # ----- left side stats panel -----
    stats_x0, stats_y0 = 160, 820
    stats_x1, stats_y1 = 1080, 1480
    gradient_box(canvas, stats_x0, stats_y0, stats_x1, stats_y1,
                 (250, 250, 250), border_color=C_PURPLE, border_px=3,
                 corner_r=18, lighten=0.0)
    cv2.putText(canvas, "Phase 2 results", (stats_x0 + 28, stats_y0 + 60),
                cv2.FONT_HERSHEY_TRIPLEX, 1.1, C_PURPLE, 3, cv2.LINE_AA)
    cv2.putText(canvas, "TCGA-THCA (n=88 slides)",
                (stats_x0 + 28, stats_y0 + 110),
                cv2.FONT_HERSHEY_DUPLEX, 0.78, C_DARK, 2, cv2.LINE_AA)
    stats_lines = [
        ("Pooled AUC", "0.746", "[0.61 - 0.86]"),
        ("Mean fold AUC", "0.830", "+/- 0.139"),
        ("All folds", ">= 0.71", "(2 of 5 = 1.00)"),
        ("Closure baseline", "~0.55", "(ResNet50)"),
        ("Verdict", "PASS", "(threshold > 0.70)"),
    ]
    y = stats_y0 + 175
    for lbl, val, ext in stats_lines:
        cv2.putText(canvas, lbl, (stats_x0 + 28, y),
                    cv2.FONT_HERSHEY_DUPLEX, 0.78, C_MID, 2, cv2.LINE_AA)
        cv2.putText(canvas, val, (stats_x0 + 380, y),
                    cv2.FONT_HERSHEY_TRIPLEX, 0.95, C_TEAL, 3, cv2.LINE_AA)
        cv2.putText(canvas, ext, (stats_x0 + 560, y),
                    cv2.FONT_HERSHEY_DUPLEX, 0.65, C_GRAY, 1, cv2.LINE_AA)
        y += 80

    # ----- footer ----
    cv2.line(canvas, (0, H - 90), (W, H - 90), (220, 220, 220), 1, cv2.LINE_AA)
    cv2.putText(canvas,
                "Color key:  TEAL = patient/H&E   |   PURPLE = model   |   "
                "AMBER = decision   |   GREEN = DM1+ branch   |   "
                "BLUE = DM2 branch   |   GRAY = standard care",
                (60, H - 35),
                cv2.FONT_HERSHEY_DUPLEX, 0.78, C_MID, 2, cv2.LINE_AA)

    # ----- save -----
    cv2.imwrite(str(PNG_OUT), canvas, [cv2.IMWRITE_PNG_COMPRESSION, 4])
    print(f"[saved] {PNG_OUT}  ({canvas.shape[1]}x{canvas.shape[0]})")

    rgb = cv2.cvtColor(canvas, cv2.COLOR_BGR2RGB)
    fig_w_in = 16.0
    fig_h_in = fig_w_in * H / W
    fig = plt.figure(figsize=(fig_w_in, fig_h_in), dpi=300)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(rgb); ax.axis("off")
    fig.savefig(PDF_OUT, dpi=300, bbox_inches="tight", pad_inches=0)
    plt.close(fig)
    print(f"[saved] {PDF_OUT}")


if __name__ == "__main__":
    main()
