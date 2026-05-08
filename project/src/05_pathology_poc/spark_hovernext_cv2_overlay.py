#!/usr/bin/env python3
"""
HoVer-NeXt 7-class cell call cv2 overlay on TCGA WSI thumbnails.

cv2 = histology image visualization only (per user instruction).

Inputs:
  project/data/processed/TCGA-THCA-WSI-DM/thumbnails/{sid}_thumb.jpg
  project/data/processed/TCGA-THCA-WSI-DM/hovernext_out/{sid}/pred_*.tsv

Output:
  project/data/processed/TCGA-THCA-WSI-DM/hovernext_overlays/{sid}_overlay.jpg
"""
from __future__ import annotations
from pathlib import Path
import cv2
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parent.parent.parent.parent
ROOTDIR = ROOT / "project/data/processed/TCGA-THCA-WSI-DM"
OUT = ROOTDIR / "hovernext_overlays"
OUT.mkdir(parents=True, exist_ok=True)

FONT = ImageFont.truetype("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc", 28)
F_SMALL = ImageFont.truetype("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc", 18)

# Cell class color palette (BGR for cv2)
PALETTE = {
    "epithelial-cell":         (220, 100, 30),    # orange-ish (thyrocyte)
    "connective-tissue-cell":  (50, 180, 220),    # yellow (CAF/stroma)
    "lymphocyte":              (50, 220, 50),     # green
    "plasma-cell":             (200, 100, 200),   # magenta
    "eosinophil":              (50, 50, 220),     # red
    "neutrophil":              (220, 50, 220),    # purple
    "mitosis":                 (255, 255, 0),     # cyan
    "dead":                    (100, 100, 100),   # gray
}
CELL_KOR = {
    "epithelial-cell": "갑상선세포",
    "connective-tissue-cell": "섬유아세포 (CAF)",
    "lymphocyte": "림프구",
    "plasma-cell": "형질세포",
    "eosinophil": "호산구",
    "neutrophil": "호중구",
    "mitosis": "유사분열",
}


def read_pred_tsvs(slide_dir):
    """Combine all pred_*.tsv into one DataFrame."""
    rows = []
    for f in sorted(slide_dir.glob("pred_*.tsv")):
        df = pd.read_csv(f, sep="\t")
        if len(df) == 0: continue
        cls = f.stem.replace("pred_", "")
        df["class"] = cls
        rows.append(df[["x", "y", "class"]])
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def render_overlay(thumb_path, slide_dir, sid):
    img = cv2.imread(str(thumb_path))  # BGR
    H, W = img.shape[:2]

    # Discover full-resolution dims (from class_inst.json or from largest cell coordinate)
    cells = read_pred_tsvs(slide_dir)
    if len(cells) == 0:
        return None

    x_max = float(cells.x.max())
    y_max = float(cells.y.max())
    # Scale factor from full-res to thumbnail
    sf_x = W / x_max if x_max > 0 else 1
    sf_y = H / y_max if y_max > 0 else 1
    sf = min(sf_x, sf_y) * 0.95

    print(f"  {sid}: {len(cells)} cells, full {x_max:.0f}x{y_max:.0f} → thumb {W}x{H} sf={sf:.4f}")

    # Subsample to keep visualization legible (max 50k per class)
    subsample = []
    for cls, grp in cells.groupby("class"):
        if len(grp) > 30000:
            grp = grp.sample(30000, random_state=0)
        subsample.append(grp)
    cells = pd.concat(subsample, ignore_index=True)

    # Make 3 panels: left=H&E original, middle=overlay, right=density per class
    panel_w = W
    canvas = np.zeros((H + 100, 3 * panel_w + 60, 3), dtype=np.uint8) + 12

    # Panel 1: H&E original
    canvas[100:100+H, :W] = img.copy()

    # Panel 2: H&E + cell overlay
    overlay = img.copy()
    radius = max(1, int(2 * sf * 100))  # ~2 pixel radius, scale-tuned
    radius = max(1, min(2, radius))
    for _, r in cells.iterrows():
        x = int(r.x * sf)
        y = int(r.y * sf)
        if 0 <= x < W and 0 <= y < H:
            color = PALETTE.get(r["class"], (180, 180, 180))
            cv2.circle(overlay, (x, y), radius, color, -1)
    canvas[100:100+H, panel_w + 30:2*panel_w + 30] = overlay

    # Panel 3: H&E + only epithelial+lymphocyte+CAF heatmap (3-class focus)
    h3 = img.copy()
    focus = cells[cells["class"].isin(["epithelial-cell", "lymphocyte", "connective-tissue-cell"])]
    for _, r in focus.iterrows():
        x = int(r.x * sf)
        y = int(r.y * sf)
        if 0 <= x < W and 0 <= y < H:
            color = PALETTE.get(r["class"], (180, 180, 180))
            cv2.circle(h3, (x, y), max(2, radius+1), color, -1)
    canvas[100:100+H, 2*panel_w + 60:3*panel_w + 60] = h3

    # Headers (with PIL for Korean text)
    pil = Image.fromarray(cv2.cvtColor(canvas, cv2.COLOR_BGR2RGB))
    draw = ImageDraw.Draw(pil)
    draw.rectangle([(0, 0), (canvas.shape[1], 100)], fill=(20, 30, 64))
    draw.text((20, 12), f"HoVer-NeXt 7-class cell calls — {sid[:18]}...", fill=(255, 255, 255), font=FONT)
    draw.text((20, 60), f"총 {len(cells):,} cells (subsampled for display) · L40S GPU 12 min", fill=(180, 200, 240), font=F_SMALL)

    # Panel labels
    labels = ["① H&E 원본", "② 7-class cell overlay", "③ 3-class focus (thyrocyte / lymphocyte / CAF)"]
    panel_xs = [W//2, panel_w + 30 + W//2, 2*panel_w + 60 + W//2]
    for lab, px in zip(labels, panel_xs):
        bbox = draw.textbbox((px, 105), lab, font=FONT, anchor="mt")
        draw.rectangle(bbox, fill=(8, 12, 26))
        draw.text((px, 105), lab, fill=(255, 255, 255), font=FONT, anchor="mt")
    canvas = cv2.cvtColor(np.array(pil), cv2.COLOR_RGB2BGR)

    # Legend (overlay panel) — bottom-left
    legend_x = panel_w + 50
    legend_y = 100 + H - 220
    cv2.rectangle(canvas, (legend_x, legend_y), (legend_x + 280, legend_y + 200), (15, 20, 40), -1)
    cv2.rectangle(canvas, (legend_x, legend_y), (legend_x + 280, legend_y + 200), (60, 80, 120), 1)
    pil = Image.fromarray(cv2.cvtColor(canvas, cv2.COLOR_BGR2RGB))
    draw = ImageDraw.Draw(pil)
    cls_counts = cells["class"].value_counts().to_dict()
    legend_classes = ["epithelial-cell", "connective-tissue-cell", "lymphocyte",
                      "plasma-cell", "eosinophil", "neutrophil"]
    for i, cls in enumerate(legend_classes):
        rgb = tuple(reversed(PALETTE[cls]))  # BGR→RGB
        cv = cls_counts.get(cls, 0)
        # Re-read original (full) cell counts (subsample shown above)
        draw.ellipse((legend_x + 12, legend_y + 14 + i*30, legend_x + 28, legend_y + 30 + i*30),
                     fill=rgb)
        draw.text((legend_x + 36, legend_y + 12 + i*30),
                  f"{CELL_KOR.get(cls, cls)} ({cv:,})", fill=(220, 230, 240), font=F_SMALL)

    # Total cells (full, not subsample) — read from disk
    full = read_pred_tsvs(slide_dir)  # original
    total_full = len(full)
    el_count = (full["class"] == "epithelial-cell").sum()
    ly_count = (full["class"] == "lymphocyte").sum()
    caf_count = (full["class"] == "connective-tissue-cell").sum()
    ratio = ly_count / max(el_count, 1)
    draw.text((legend_x + 12, legend_y + 175),
              f"L:E ratio = {ratio:.2f}", fill=(252, 211, 77), font=F_SMALL)
    canvas = cv2.cvtColor(np.array(pil), cv2.COLOR_RGB2BGR)

    out = OUT / f"{sid[:25]}_overlay.jpg"
    if max(canvas.shape) > 4000:
        scale = 3600 / max(canvas.shape)
        canvas = cv2.resize(canvas, (int(canvas.shape[1]*scale), int(canvas.shape[0]*scale)),
                           interpolation=cv2.INTER_AREA)
    cv2.imwrite(str(out), canvas, [cv2.IMWRITE_JPEG_QUALITY, 88])
    print(f"  wrote {out}")


def main():
    import json
    thumbs = list((ROOTDIR / "thumbnails").glob("*_thumb.jpg"))
    for tp in thumbs:
        sid = tp.name.replace("_thumb.jpg", "")
        slide_dir = ROOTDIR / "hovernext_out" / sid
        if not slide_dir.exists():
            print(f"  skip {sid}: no HoVer-NeXt out"); continue
        render_overlay(tp, slide_dir, sid)


if __name__ == "__main__":
    main()
