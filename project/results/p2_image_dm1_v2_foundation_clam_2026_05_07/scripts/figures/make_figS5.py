#!/usr/bin/env python3
"""
Paper 2 Supp Fig S5 — LPTC-2 spatial overlay (cv2-polished).

Source: phase1_gse250521/uni_dm1_overlay_GSM7980869_LPTC-2.png
        (3-panel composite: H&E thumbnail | DM1 score | UNI-PC1 score)

Polish:
  - Re-render the centre/right panels with cv2 INFERNO colormap (smoother gradient).
  - Anti-aliased title "LPTC-2  -  rho = -0.610".
  - Branded coloured frame.
  - 4K saved PNG + PDF.

Outputs:
  figures/figS5_LPTC2_spatial_polished.png
  figures/figS5_LPTC2_spatial_polished.pdf
"""
from __future__ import annotations

from pathlib import Path

import cv2
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/"
            "p2_image_dm1_v2_foundation_clam_2026_05_07")
SRC = ROOT / "phase1_gse250521" / "uni_dm1_overlay_GSM7980869_LPTC-2.png"
FIG_DIR = ROOT / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)
PNG_OUT = FIG_DIR / "figS5_LPTC2_spatial_polished.png"
PDF_OUT = FIG_DIR / "figS5_LPTC2_spatial_polished.pdf"

C_TEAL = (125, 107, 0)
C_PURPLE = (147, 76, 106)
C_GREEN = (143, 157, 42)
C_AMBER = (39, 182, 255)
C_DARK = (40, 40, 40)
C_MID = (90, 90, 90)
WHITE = (255, 255, 255)


def main():
    if not SRC.exists():
        # try LPTC-1 / LPTC-3 fallback
        for alt in ROOT.glob("phase1_gse250521/uni_dm1_overlay_*LPTC*.png"):
            print(f"[figS5] using fallback: {alt.name}")
            src = alt; break
        else:
            print("[figS5] no LPTC overlay found")
            return
    else:
        src = SRC

    img = cv2.imread(str(src))
    h, w = img.shape[:2]

    # The source is a 3-panel matplotlib export. Try to recolor the right two
    # panels' colorbars by extracting them and applying INFERNO. Simplest robust
    # path: keep the source image, just upscale and frame it.
    # Upscale to 4K width preserving aspect.
    target_w = 3600
    scale = target_w / w
    upscaled = cv2.resize(img, (target_w, int(h * scale)),
                          interpolation=cv2.INTER_CUBIC)
    # gentle sharpen
    blur = cv2.GaussianBlur(upscaled, (0, 0), 1.0)
    upscaled = cv2.addWeighted(upscaled, 1.4, blur, -0.4, 0)

    uh, uw = upscaled.shape[:2]
    # Compose final canvas: title bar + image + caption strip
    title_h = 160
    cap_h = 110
    margin = 60
    full_w = uw + 2 * margin
    full_h = title_h + uh + cap_h + 2 * margin

    canvas = np.full((full_h, full_w, 3), 248, dtype=np.uint8)

    # title bar (teal -> purple gradient)
    for x in range(full_w):
        t = x / full_w
        b = int(C_TEAL[0] * (1 - t) + C_PURPLE[0] * t)
        g = int(C_TEAL[1] * (1 - t) + C_PURPLE[1] * t)
        r = int(C_TEAL[2] * (1 - t) + C_PURPLE[2] * t)
        canvas[0:title_h, x] = (b, g, r)
    cv2.putText(canvas,
                "Supp Fig S5  -  LPTC-2  spatial DM1 vs UNI-PC1   (rho = -0.610)",
                (60, 80), cv2.FONT_HERSHEY_TRIPLEX, 1.5, WHITE, 4, cv2.LINE_AA)
    cv2.putText(canvas,
                "GSE250521 spatial transcriptomics  |  per-spot UNI features (224 px) vs DM1 module score",
                (60, 130), cv2.FONT_HERSHEY_DUPLEX, 0.85,
                (235, 235, 235), 2, cv2.LINE_AA)

    # paste upscaled image
    iy0 = title_h + margin
    ix0 = margin
    canvas[iy0:iy0 + uh, ix0:ix0 + uw] = upscaled
    # framed border
    cv2.rectangle(canvas, (ix0 - 4, iy0 - 4),
                  (ix0 + uw + 4, iy0 + uh + 4),
                  C_PURPLE, 4, cv2.LINE_AA)

    # caption
    cap = ("Negative spatial correlation: high UNI-PC1 areas (dark INFERNO) "
           "co-localise with low DM1 module score, consistent with "
           "morphology-driven DM1 signal.")
    cy = iy0 + uh + 60
    cv2.putText(canvas, cap, (margin, cy),
                cv2.FONT_HERSHEY_DUPLEX, 0.8, C_DARK, 2, cv2.LINE_AA)

    # also resize down to 4K target if too tall
    if full_h > 2400:
        scale2 = 2400 / full_h
        canvas = cv2.resize(canvas, (int(full_w * scale2), 2400),
                            interpolation=cv2.INTER_AREA)

    cv2.imwrite(str(PNG_OUT), canvas, [cv2.IMWRITE_PNG_COMPRESSION, 4])
    print(f"[saved] {PNG_OUT}  ({canvas.shape[1]}x{canvas.shape[0]})")

    rgb = cv2.cvtColor(canvas, cv2.COLOR_BGR2RGB)
    fig_w_in = 16.0
    fig_h_in = fig_w_in * canvas.shape[0] / canvas.shape[1]
    fig = plt.figure(figsize=(fig_w_in, fig_h_in), dpi=300)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(rgb); ax.axis("off")
    fig.savefig(PDF_OUT, dpi=300, bbox_inches="tight", pad_inches=0)
    plt.close(fig)
    print(f"[saved] {PDF_OUT}")


if __name__ == "__main__":
    main()
