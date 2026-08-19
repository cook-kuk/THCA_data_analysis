#!/usr/bin/env python3
"""Build H&E-first overlay viewer assets for Paper 2."""
from __future__ import annotations

import json
import gzip
import re
import shutil
from pathlib import Path

import cv2
import numpy as np
from PIL import Image, ImageOps


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
P2 = ROOT / "project/results/p2_image_dm1_v2_foundation_clam_2026_05_07"
PATHO = ROOT / "project/results/03_pathology_poc"
GSE230_RAW = ROOT / "project/data/external/GSE230424/raw"
GSE250_TILES = ROOT / "project/data/processed/GSE250521/tiles"
OUT = P2 / "analysis_supp/paper2_he_overlay_viewer_2026_05_11"
HUB = ROOT / "project/papers_hub_2026_05_04"
LIVE = Path("/var/www/papers/papers_hub_2026_05_04")
ASSET = "paper2_he_overlay_viewer"
ASSET_LOCAL = HUB / "assets" / ASSET
ASSET_LIVE = LIVE / "assets" / ASSET

BG = (13, 17, 23)
PANEL = (23, 30, 41)
PANEL2 = (31, 40, 54)
LINE = (68, 79, 94)
INK = (235, 240, 246)
MUTED = (156, 169, 188)
GOLD = (94, 178, 214)
TEAL = (196, 213, 99)
BLUE = (255, 166, 87)
RED = (114, 123, 255)
GREEN = (157, 211, 53)
WHITE = (245, 246, 248)

Image.MAX_IMAGE_PIXELS = None


def put(img, text, xy, scale=0.55, color=INK, thick=1, width=None, gap=8):
    x, y = xy
    font = cv2.FONT_HERSHEY_SIMPLEX
    lines = []
    for raw in str(text).split("\n"):
        if width is None:
            lines.append(raw)
            continue
        line = ""
        for word in raw.split():
            trial = word if not line else f"{line} {word}"
            if cv2.getTextSize(trial, font, scale, thick)[0][0] <= width or not line:
                line = trial
            else:
                lines.append(line)
                line = word
        lines.append(line)
    step = cv2.getTextSize("Ag", font, scale, thick)[0][1] + gap
    for i, line in enumerate(lines):
        cv2.putText(img, line, (x, y + i * step), font, scale, color, thick, cv2.LINE_AA)
    return y + max(len(lines), 1) * step


def box(img, xyxy, color=PANEL, outline=LINE):
    x0, y0, x1, y1 = xyxy
    cv2.rectangle(img, (x0, y0), (x1, y1), color, -1)
    cv2.rectangle(img, (x0, y0), (x1, y1), outline, 1, cv2.LINE_AA)


def title(img, kicker, heading, sub):
    put(img, kicker.upper(), (70, 78), 0.58, GOLD, 2)
    put(img, heading, (70, 142), 1.08, INK, 2, 1700)
    put(img, sub, (70, 205), 0.56, MUTED, 1, 1560)


def read_image(path: Path) -> np.ndarray:
    img = cv2.imread(str(path))
    if img is None:
        canvas = np.full((500, 700, 3), PANEL2, dtype=np.uint8)
        put(canvas, f"Missing: {path.name}", (30, 250), 0.5, RED, 1, 620)
        return canvas
    return img


def read_raw_jpeg_gz(path: Path, max_side: int = 1400) -> np.ndarray:
    with gzip.open(path, "rb") as f:
        im = Image.open(f)
        im.draft("RGB", (max_side * 2, max_side * 2))
        im.load()
        im = im.convert("RGB")
        im = ImageOps.contain(im, (max_side, max_side))
    return cv2.cvtColor(np.array(im), cv2.COLOR_RGB2BGR)


def fit_image(path: Path, size: tuple[int, int], background=WHITE) -> np.ndarray:
    tw, th = size
    src = read_image(path)
    h, w = src.shape[:2]
    scale = min(tw / max(w, 1), th / max(h, 1))
    nw, nh = max(1, int(w * scale)), max(1, int(h * scale))
    resized = cv2.resize(src, (nw, nh), interpolation=cv2.INTER_AREA)
    canvas = np.full((th, tw, 3), background, dtype=np.uint8)
    x0 = (tw - nw) // 2
    y0 = (th - nh) // 2
    canvas[y0 : y0 + nh, x0 : x0 + nw] = resized
    return canvas


def card(img, path: Path, xyxy, head: str, caption: str, color=GOLD):
    x0, y0, x1, y1 = xyxy
    box(img, (x0, y0, x1, y1), PANEL2, color)
    put(img, head, (x0 + 22, y0 + 48), 0.58, color, 2, x1 - x0 - 44)
    put(img, caption, (x0 + 22, y0 + 88), 0.38, MUTED, 1, x1 - x0 - 44, 5)
    panel = fit_image(path, (x1 - x0 - 44, y1 - y0 - 138))
    img[y0 + 116 : y1 - 22, x0 + 22 : x1 - 22] = panel


def panel_array(img, arr: np.ndarray, xyxy, head: str, caption: str, color=GOLD):
    x0, y0, x1, y1 = xyxy
    box(img, (x0, y0, x1, y1), PANEL2, color)
    put(img, head, (x0 + 22, y0 + 48), 0.58, color, 2, x1 - x0 - 44)
    put(img, caption, (x0 + 22, y0 + 88), 0.38, MUTED, 1, x1 - x0 - 44, 5)
    h, w = arr.shape[:2]
    scale = min((x1 - x0 - 44) / max(w, 1), (y1 - y0 - 138) / max(h, 1))
    nw, nh = max(1, int(w * scale)), max(1, int(h * scale))
    resized = cv2.resize(arr, (nw, nh), interpolation=cv2.INTER_AREA)
    canvas = np.full((y1 - y0 - 138, x1 - x0 - 44, 3), WHITE, dtype=np.uint8)
    x = (canvas.shape[1] - nw) // 2
    y = (canvas.shape[0] - nh) // 2
    canvas[y : y + nh, x : x + nw] = resized
    img[y0 + 116 : y1 - 22, x0 + 22 : x1 - 22] = canvas


def source_paths() -> dict[str, Path]:
    return {
        "g230_he": HUB / "assets/paper2_cv2_summary/fig02_gse230424_cv2_he_overlay_mosaic.png",
        "g230_hotspot": HUB / "assets/paper2_cv2_summary/fig05_gse230424_cv2_hotspot_contours.png",
        "g250_mosaic": HUB / "assets/paper2_cv2_summary/fig03_gse250521_cv2_uni_slide_mosaic.png",
        "scorecard": HUB / "assets/paper2_cv2_summary/fig01_paper2_cv2_evidence_scorecard.png",
    }


def slide_paths() -> list[Path]:
    return sorted((HUB / "assets/paper2_individual_figures").glob("phase1_gse250521__uni_dm1_overlay_*.png"))


def draw_storyboard(paths: dict[str, Path]) -> Path:
    img = np.full((1650, 2300, 3), BG, dtype=np.uint8)
    title(
        img,
        "H&E-first overlay viewer",
        "Read The Paper 2 Signal Directly On Histology",
        "The old separated browser remains as an appendix; this board puts H&E overlay evidence first.",
    )
    card(
        img,
        paths["g230_he"],
        (70, 310, 1125, 960),
        "GSE230424 H&E + DM1/RAI spots",
        "True H&E background with image-derived thyroid/DM1 axis overlay.",
        TEAL,
    )
    card(
        img,
        paths["g230_hotspot"],
        (1175, 310, 2230, 960),
        "GSE230424 H&E + hotspot contours",
        "Observed/predicted top-decile contours over H&E; designed for visual review.",
        GOLD,
    )
    card(
        img,
        paths["g250_mosaic"],
        (70, 1030, 1125, 1530),
        "GSE250521 H&E + UNI-DM1 mosaic",
        "Sixteen thyroid spatial slides condensed only after individual overlays are available.",
        BLUE,
    )
    card(
        img,
        paths["scorecard"],
        (1175, 1030, 2230, 1530),
        "Evidence context",
        "Classifier, spatial recovery, external support, and caveats kept beside the overlays.",
        GREEN,
    )
    box(img, (70, 1575, 2230, 1630), (42, 34, 22), GOLD)
    put(img, "Use this first for review. Use the 81-separated browser only when you need a source panel or control figure.", (105, 1614), 0.54, GOLD, 2, 2000)
    path = OUT / "F01_he_first_overlay_storyboard.png"
    cv2.imwrite(str(path), img)
    return path


def draw_g250_contact_sheet(slides: list[Path]) -> Path:
    img = np.full((1720, 2400, 3), BG, dtype=np.uint8)
    title(
        img,
        "GSE250521 slide overlay contact sheet",
        "All 16 H&E + UNI-DM1 Slide Overlays In One Clean View",
        "Each source remains clickable in the web viewer; this contact sheet is for fast scanning.",
    )
    cell_w, cell_h = 540, 245
    x0, y0 = 70, 325
    for i, path in enumerate(slides[:16]):
        row, col = divmod(i, 4)
        x = x0 + col * 575
        y = y0 + row * 310
        label = re.sub(r"^phase1_gse250521__uni_dm1_overlay_", "", path.stem).replace("_", " ")
        box(img, (x, y, x + cell_w, y + cell_h), PANEL2, TEAL if "PTC" in label else GOLD if "ATC" in label else BLUE if "LPTC" in label else LINE)
        put(img, label, (x + 14, y + 36), 0.42, INK, 2, cell_w - 28)
        panel = fit_image(path, (cell_w - 28, cell_h - 66))
        img[y + 52 : y + cell_h - 14, x + 14 : x + cell_w - 14] = panel
    path = OUT / "F02_gse250521_he_slide_overlay_contact_sheet.png"
    cv2.imwrite(str(path), img)
    return path


def draw_g230_pair(paths: dict[str, Path]) -> Path:
    img = np.full((1450, 2300, 3), BG, dtype=np.uint8)
    title(
        img,
        "GSE230424 H&E overlay pair",
        "Signal And Hotspots On The Same Histology Substrate",
        "This is the fastest read for the external thyroid Visium support layer.",
    )
    card(
        img,
        paths["g230_he"],
        (80, 320, 1110, 1275),
        "H&E + image-derived DM1/RAI axis",
        "Sample-centered rho 0.644; residual-target rho 0.232.",
        TEAL,
    )
    card(
        img,
        paths["g230_hotspot"],
        (1190, 320, 2220, 1275),
        "H&E + observed/predicted hotspot contours",
        "Top-decile raw lift 3.32x; residual hotspot lift 1.63x.",
        GOLD,
    )
    path = OUT / "F03_gse230424_he_overlay_pair.png"
    cv2.imwrite(str(path), img)
    return path


def draw_legend() -> Path:
    img = np.full((1200, 1900, 3), BG, dtype=np.uint8)
    title(
        img,
        "Overlay reading guide",
        "How To Use The H&E-Fused Viewer",
        "The viewer is not a new analysis. It is a cleaner visual reading layer over existing Paper 2 assets.",
    )
    rows = [
        ("Start here", "H&E overlay storyboard", "fastest visual read"),
        ("Inspect external thyroid", "GSE230424 H&E pair", "signal + hotspot contours"),
        ("Inspect spatial trajectory", "GSE250521 16-slide contact sheet", "stage-by-stage overlay scan"),
        ("Deep dive", "individual slide overlays", "click a single WSI-derived overlay"),
        ("Appendix only", "81 separated figures", "controls/source panels"),
    ]
    for i, (head, body, tag) in enumerate(rows):
        y = 330 + i * 150
        color = [GOLD, TEAL, BLUE, GREEN, RED][i]
        box(img, (90, y, 1810, y + 105), PANEL if i % 2 == 0 else PANEL2, color)
        put(img, head, (125, y + 58), 0.6, color, 2, 360)
        put(img, body, (520, y + 58), 0.54, INK, 1, 620)
        put(img, tag, (1230, y + 58), 0.48, MUTED, 1, 460)
    path = OUT / "F04_he_overlay_reading_guide.png"
    cv2.imwrite(str(path), img)
    return path


def draw_internal_pathology_anchor() -> Path:
    img = np.full((1680, 2320, 3), BG, dtype=np.uint8)
    title(
        img,
        "Internal pathology anchor",
        "03_pathology_poc Connects Image-DM1 To A Tight Slide-Level Readout",
        "This board adds an internal pathology bridge without collapsing it into the H&E overlay appendix.",
    )
    panels = [
        (PATHO / "summary_figs/F1_trajectory.png", (70, 315, 1130, 660), "Trajectory", "Internal slide trajectory summary."),
        (PATHO / "summary_figs/F2_cross_platform_replicate.png", (1170, 315, 2250, 660), "Cross-platform replicate", "Replicate consistency across slide groups."),
        (PATHO / "summary_figs/F4_LR_heatmap.png", (70, 725, 1130, 1140), "LR heatmap", "Slide-level LR activity map."),
        (PATHO / "summary_figs/F8_per_case_predictions.png", (1170, 725, 2250, 1140), "Per-case predictions", "Case-level prediction spread."),
        (PATHO / "pred_vs_obs_resnet50_per_slide.png", (70, 1195, 1450, 1590), "Pred vs obs", "Observed versus predicted per slide."),
        (PATHO / "closure_battery_summary.png", (1490, 1195, 2250, 1590), "Closure battery", "Internal closure / audit summary."),
    ]
    for path, xyxy, head, caption in panels:
        card(img, path, xyxy, head, caption, GREEN if "replicate" in head.lower() else TEAL if "trajectory" in head.lower() else GOLD)
    path = OUT / "F05_internal_pathology_anchor.png"
    cv2.imwrite(str(path), img)
    return path


def draw_internal_pathology_atlas() -> Path:
    img = np.full((1720, 2460, 3), BG, dtype=np.uint8)
    title(
        img,
        "Internal pathology atlas",
        "Real Slide-Level Cases From 03_pathology_poc",
        "This replaces another synthetic summary with actual internal slide overlays spanning normal, PTC, LPTC, and ATC.",
    )
    panels = [
        (PATHO / "cci_overlays/GSM7980860_N-1.png", (70, 315, 1180, 765), "Normal N-1", "Internal normal reference case."),
        (PATHO / "cci_overlays/GSM7980864_PTC-1.png", (1240, 315, 2390, 765), "PTC-1", "Internal papillary thyroid carcinoma case."),
        (PATHO / "cci_overlays/GSM7980868_LPTC-1.png", (70, 820, 1180, 1270), "LPTC-1", "Internal low-papillary thyroid carcinoma case."),
        (PATHO / "cci_overlays/GSM7980872_ATC-1.png", (1240, 820, 2390, 1270), "ATC-1", "Internal anaplastic thyroid carcinoma case."),
        (PATHO / "spatial_overlays/GSM7980860_N-1.png", (70, 1325, 1180, 1650), "Spatial N-1", "Spatial overlay counterpart."),
        (PATHO / "spatial_overlays/GSM7980864_PTC-1.png", (1240, 1325, 2390, 1650), "Spatial PTC-1", "Spatial overlay counterpart."),
    ]
    for path, xyxy, head, caption in panels:
        card(img, path, xyxy, head, caption, GREEN if "Normal" in head else RED if "ATC" in head else TEAL if "PTC" in head else GOLD)
    path = OUT / "F11_internal_pathology_atlas.png"
    cv2.imwrite(str(path), img)
    return path


def draw_internal_pathology_full_atlas() -> Path:
    img = np.full((2060, 2580, 3), BG, dtype=np.uint8)
    title(
        img,
        "Internal pathology full atlas",
        "All 16 Internal CCI Overlays At Once",
        "This is the strongest internal slide-level view: N, PTC, LPTC, and ATC across the full case set.",
    )
    cases = [
        ("N-1", PATHO / "cci_overlays/GSM7980860_N-1.png"),
        ("N-2", PATHO / "cci_overlays/GSM7980861_N-2.png"),
        ("N-3", PATHO / "cci_overlays/GSM7980862_N-3.png"),
        ("N-4", PATHO / "cci_overlays/GSM7980863_N-4.png"),
        ("PTC-1", PATHO / "cci_overlays/GSM7980864_PTC-1.png"),
        ("PTC-2", PATHO / "cci_overlays/GSM7980865_PTC-2.png"),
        ("PTC-3", PATHO / "cci_overlays/GSM7980866_PTC-3.png"),
        ("PTC-4", PATHO / "cci_overlays/GSM7980867_PTC-4.png"),
        ("LPTC-1", PATHO / "cci_overlays/GSM7980868_LPTC-1.png"),
        ("LPTC-2", PATHO / "cci_overlays/GSM7980869_LPTC-2.png"),
        ("LPTC-3", PATHO / "cci_overlays/GSM7980870_LPTC-3.png"),
        ("LPTC-4", PATHO / "cci_overlays/GSM7980871_LPTC-4.png"),
        ("ATC-1", PATHO / "cci_overlays/GSM7980872_ATC-1.png"),
        ("ATC-2", PATHO / "cci_overlays/GSM7980873_ATC-2.png"),
        ("ATC-3", PATHO / "cci_overlays/GSM7980874_ATC-3.png"),
        ("ATC-4", PATHO / "cci_overlays/GSM7980875_ATC-4.png"),
    ]
    cell_w, cell_h = 555, 275
    x0, y0 = 70, 315
    for i, (label, path) in enumerate(cases):
        row, col = divmod(i, 4)
        x = x0 + col * 610
        y = y0 + row * 370
        color = GREEN if label.startswith("N") else TEAL if label.startswith("PTC") else BLUE if label.startswith("LPTC") else RED
        box(img, (x, y, x + cell_w, y + cell_h), PANEL2, color)
        put(img, label, (x + 14, y + 34), 0.44, color, 2, cell_w - 28)
        panel = fit_image(path, (cell_w - 28, cell_h - 58))
        img[y + 46 : y + cell_h - 12, x + 14 : x + cell_w - 14] = panel
    box(img, (70, 1795, 2510, 1965), (18, 23, 31), GOLD)
    put(img, "16 internal cases across normal, papillary, low-papillary, and anaplastic states. Use this as the internal pathology front page.", (110, 1872), 0.62, INK, 2, 2320)
    path = OUT / "F12_internal_pathology_full_atlas.png"
    cv2.imwrite(str(path), img)
    return path


def draw_raw_he_gallery() -> Path:
    slides = [
        ("GSE230424 P1", GSE230_RAW / "GSM7221915_P1_HE.jpg.gz"),
        ("GSE230424 P2", GSE230_RAW / "GSM7221916_P2_HE.jpg.gz"),
        ("GSE230424 P3", GSE230_RAW / "GSM7221917_P3_HE.jpg.gz"),
        ("GSE230424 P4", GSE230_RAW / "GSM7221918_P4_HE.jpg.gz"),
    ]
    img = np.full((2080, 2680, 3), BG, dtype=np.uint8)
    title(
        img,
        "Raw H&E",
        "GSE230424 Original H&E Slides",
        "This is the actual source H&E, downsampled for review, not a synthetic overlay.",
    )
    for i, (label, path) in enumerate(slides):
        raw = read_raw_jpeg_gz(path, max_side=1300)
        row, col = divmod(i, 2)
        x0 = 70 + col * 1280
        y0 = 315 + row * 840
        panel_array(
            img,
            raw,
            (x0, y0, x0 + 1200, y0 + 760),
            label,
            "Original slide image from the raw GSE230424 source package.",
            TEAL if col == 0 else GOLD,
        )
    box(img, (70, 1925, 2610, 2005), (18, 23, 31), GOLD)
    put(img, "Use this as the true external H&E front page. Everything else in Paper 2 is downstream of this substrate.", (110, 1973), 0.62, INK, 2, 2400)
    path = OUT / "F13_raw_he_gallery.png"
    cv2.imwrite(str(path), img)
    return path


def draw_raw_source_wall() -> Path:
    raw_slides = [
        ("P1", read_raw_jpeg_gz(GSE230_RAW / "GSM7221915_P1_HE.jpg.gz", max_side=1250)),
        ("P2", read_raw_jpeg_gz(GSE230_RAW / "GSM7221916_P2_HE.jpg.gz", max_side=1250)),
        ("P3", read_raw_jpeg_gz(GSE230_RAW / "GSM7221917_P3_HE.jpg.gz", max_side=1250)),
        ("P4", read_raw_jpeg_gz(GSE230_RAW / "GSM7221918_P4_HE.jpg.gz", max_side=1250)),
    ]
    tile_sets = [
        GSE250_TILES / "GSM7980860_N-1",
        GSE250_TILES / "GSM7980864_PTC-1",
        GSE250_TILES / "GSM7980868_LPTC-1",
        GSE250_TILES / "GSM7980872_ATC-1",
    ]
    tile_cells: list[np.ndarray] = []
    for folder in tile_sets:
        for path in sorted(folder.glob("*_size672.png"))[:4]:
            tile_cells.append(read_image(path))

    img = np.full((2360, 2720, 3), BG, dtype=np.uint8)
    title(
        img,
        "Raw source wall",
        "Actual Source H&E And Real Tile Atlas",
        "This is the strongest front page: raw GSE230424 slides on the left, actual GSE250521 tiles on the right.",
    )

    # Left: four raw source slides in a 2x2 wall.
    left_x, top_y = 70, 315
    cell_w, cell_h = 560, 700
    for i, (label, arr) in enumerate(raw_slides):
        row, col = divmod(i, 2)
        x = left_x + col * 600
        y = top_y + row * 745
        panel_array(
            img,
            arr,
            (x, y, x + cell_w, y + cell_h),
            f"Raw {label}",
            "Downsampled source H&E slide.",
            TEAL if i % 2 == 0 else GOLD,
        )

    # Right: sixteen real tiles in a 4x4 atlas.
    atlas_x = 1390
    atlas_y = 315
    for i, arr in enumerate(tile_cells[:16]):
        row, col = divmod(i, 4)
        x = atlas_x + col * 320
        y = atlas_y + row * 390
        panel_array(
            img,
            arr,
            (x, y, x + 300, y + 340),
            f"Tile {i + 1:02d}",
            "Actual slide tile crop.",
            GREEN if i < 4 else TEAL if i < 8 else BLUE if i < 12 else RED,
        )

    box(img, (70, 2130, 2650, 2200), (18, 23, 31), GOLD)
    put(img, "Raw slides + real tiles are the front page now. The overlay browser follows underneath, and the 81 separated figures stay as appendix.", (110, 2177), 0.60, INK, 2, 2480)
    path = OUT / "F15_raw_source_wall.png"
    cv2.imwrite(str(path), img)
    return path


def draw_real_tile_atlas() -> Path:
    tile_sets = [
        ("Normal", GSE250_TILES / "GSM7980860_N-1"),
        ("PTC", GSE250_TILES / "GSM7980864_PTC-1"),
        ("LPTC", GSE250_TILES / "GSM7980868_LPTC-1"),
        ("ATC", GSE250_TILES / "GSM7980872_ATC-1"),
    ]
    picks: list[tuple[str, Path]] = []
    for prefix, folder in tile_sets:
        tiles = sorted(folder.glob("*_size672.png"))[:4]
        for p in tiles:
            picks.append((prefix, p))
    img = np.full((1700, 2550, 3), BG, dtype=np.uint8)
    title(
        img,
        "Real tiles",
        "GSE250521 Actual Tile Atlas",
        "16 real tile crops across normal, PTC, LPTC, and ATC. This is the slide-level image evidence in compact form.",
    )
    for i, (prefix, path) in enumerate(picks):
        row, col = divmod(i, 4)
        x0 = 70 + col * 595
        y0 = 315 + row * 345
        box(img, (x0, y0, x0 + 555, y0 + 295), PANEL2, GREEN if prefix == "Normal" else TEAL if prefix == "PTC" else BLUE if prefix == "LPTC" else RED)
        put(img, prefix, (x0 + 12, y0 + 32), 0.42, INK, 2, 250)
        panel = fit_image(path, (525, 208), WHITE)
        img[y0 + 48 : y0 + 256, x0 + 15 : x0 + 540] = panel
    box(img, (70, 1595, 2480, 1670), (18, 23, 31), GOLD)
    put(img, "Actual slide tiles, not summaries. Use this as the compact internal image atlas for the manuscript.", (110, 1639), 0.60, INK, 2, 2300)
    path = OUT / "F14_real_tile_atlas.png"
    cv2.imwrite(str(path), img)
    return path


def draw_front_door_max_impact(paths: dict[str, Path]) -> Path:
    img = np.full((1840, 2360, 3), BG, dtype=np.uint8)
    title(
        img,
        "Front door",
        "Paper 2 Max Impact Assembly",
        "If a reviewer sees only one page, this is the page: histology, validation, impact, and internal pathology in one view.",
    )
    card(
        img,
        paths["g230_he"],
        (70, 315, 1140, 880),
        "External H&E anchor",
        "GSE230424 H&E + DM1/RAI overlay mosaic.",
        TEAL,
    )
    card(
        img,
        PATHO / "summary_figs/F1_trajectory.png",
        (1220, 315, 2290, 880),
        "Internal pathology anchor",
        "Internal pathology trajectory and slide-level readout.",
        GREEN,
    )
    card(
        img,
        PATHO / "summary_figs/F2_cross_platform_replicate.png",
        (70, 950, 750, 1490),
        "Replicate support",
        "Cross-platform replicate consistency.",
        BLUE,
    )
    card(
        img,
        HUB / "assets/paper2_moonshot_editorial_pack/F06_natcom_2026_precedent_bridge.png",
        (820, 950, 1540, 1490),
        "Nature Communications bridge",
        "Target-tier precedent and framing.",
        GOLD,
    )
    card(
        img,
        HUB / "assets/paper2_external_validation_ready/F05_pathologist_hotspot_audit_sheet.png",
        (1590, 950, 2290, 1490),
        "Validation audit",
        "Pathologist hotspot audit and validation lock.",
        RED,
    )
    card(
        img,
        HUB / "assets/paper2_max_impact_unlock/F02_impact_unlock_ladder.png",
        (70, 1565, 2290, 1788),
        "Impact ladder",
        "External validation, claim control, and manuscript unlock path.",
        GOLD,
    )
    path = OUT / "F06_front_door_max_impact.png"
    cv2.imwrite(str(path), img)
    return path


def draw_executive_decision_board(paths: dict[str, Path]) -> Path:
    img = np.full((1760, 2480, 3), BG, dtype=np.uint8)
    title(
        img,
        "Executive board",
        "Paper 2 Decision Surface",
        "This is the version you show first: a clean synthesis of external H&E, internal pathology, validation, and journal framing.",
    )
    card(
        img,
        paths["g230_he"],
        (70, 315, 1210, 930),
        "External H&E proof",
        "GSE230424 H&E + image-derived DM1/RAI signal.",
        TEAL,
    )
    card(
        img,
        PATHO / "summary_figs/F1_trajectory.png",
        (1270, 315, 2410, 930),
        "Internal pathology proof",
        "Trajectory and slide-level readout on local pathology data.",
        GREEN,
    )
    card(
        img,
        HUB / "assets/paper2_external_validation_ready/F05_pathologist_hotspot_audit_sheet.png",
        (70, 995, 800, 1470),
        "Validation proof",
        "Pathologist audit and locked validation lane.",
        RED,
    )
    card(
        img,
        HUB / "assets/paper2_moonshot_editorial_pack/F06_natcom_2026_precedent_bridge.png",
        (840, 995, 1610, 1470),
        "Journal bridge",
        "Nature Communications precedent and tier framing.",
        GOLD,
    )
    card(
        img,
        HUB / "assets/paper2_max_impact_unlock/F02_impact_unlock_ladder.png",
        (1680, 995, 2410, 1470),
        "Claim ladder",
        "The unlock path from evidence to manuscript.",
        BLUE,
    )
    box(img, (70, 1540, 2410, 1695), (18, 23, 31), GOLD)
    put(img, "Front door readout: H&E first, internal pathology second, validation and target-tier bridge in the same frame.", (110, 1608), 0.62, INK, 2, 2240)
    path = OUT / "F07_executive_decision_board.png"
    cv2.imwrite(str(path), img)
    return path


def draw_cover_synthesis_board(paths: dict[str, Path]) -> Path:
    img = np.full((1920, 2520, 3), BG, dtype=np.uint8)
    title(
        img,
        "Cover synthesis",
        "Paper 2 Reviewer-Ready Synthesis",
        "One opening board: external H&E, internal pathology, validation, and the target-journal bridge.",
    )
    card(
        img,
        paths["g230_he"],
        (70, 315, 1225, 940),
        "External histology anchor",
        "GSE230424 H&E + image-derived DM1/RAI overlay.",
        TEAL,
    )
    card(
        img,
        PATHO / "summary_figs/F1_trajectory.png",
        (1295, 315, 2450, 620),
        "Internal pathology trajectory",
        "Local slide-level signal and progression readout.",
        GREEN,
    )
    card(
        img,
        PATHO / "summary_figs/F2_cross_platform_replicate.png",
        (1295, 655, 2450, 940),
        "Internal replicate support",
        "Cross-platform replicate consistency.",
        BLUE,
    )
    card(
        img,
        HUB / "assets/paper2_external_validation_ready/F05_pathologist_hotspot_audit_sheet.png",
        (70, 1005, 805, 1490),
        "Validation audit",
        "Pathologist hotspot audit and locked validation lane.",
        RED,
    )
    card(
        img,
        HUB / "assets/paper2_moonshot_editorial_pack/F06_natcom_2026_precedent_bridge.png",
        (840, 1005, 1575, 1490),
        "Nature Communications bridge",
        "Target-tier precedent and submission framing.",
        GOLD,
    )
    card(
        img,
        HUB / "assets/paper2_max_impact_unlock/F02_impact_unlock_ladder.png",
        (1610, 1005, 2450, 1490),
        "Claim ladder",
        "Unlock path from evidence to manuscript.",
        GOLD,
    )
    card(
        img,
        paths["g230_hotspot"],
        (70, 1540, 1225, 1860),
        "Hotspot contours",
        "Observed vs predicted contours over H&E.",
        GOLD,
    )
    card(
        img,
        PATHO / "summary_figs/F8_per_case_predictions.png",
        (1295, 1540, 2450, 1860),
        "Per-case pathology readout",
        "Case-level response spread in the internal series.",
        GREEN,
    )
    path = OUT / "F08_cover_synthesis_board.png"
    cv2.imwrite(str(path), img)
    return path


def draw_cover_executive_board(paths: dict[str, Path]) -> Path:
    img = np.full((1980, 2600, 3), BG, dtype=np.uint8)
    title(
        img,
        "Cover executive",
        "Paper 2 Impact Cover",
        "External H&E, internal pathology, validation, and the journal bridge compressed into a single opening board.",
    )
    card(
        img,
        paths["g230_he"],
        (70, 315, 1260, 980),
        "External histology anchor",
        "GSE230424 H&E + DM1/RAI overlay mosaic; the outside-in visual proof.",
        TEAL,
    )
    card(
        img,
        PATHO / "summary_figs/F1_trajectory.png",
        (1340, 315, 2530, 725),
        "Internal pathology anchor",
        "Local pathology trajectory and slide-level readout.",
        GREEN,
    )
    card(
        img,
        PATHO / "summary_figs/F2_cross_platform_replicate.png",
        (1340, 760, 2530, 980),
        "Replicate support",
        "Cross-platform replicate consistency.",
        BLUE,
    )
    card(
        img,
        HUB / "assets/paper2_external_validation_ready/F05_pathologist_hotspot_audit_sheet.png",
        (70, 1045, 830, 1500),
        "Validation audit",
        "Pathologist hotspot audit and locked validation lane.",
        RED,
    )
    card(
        img,
        HUB / "assets/paper2_moonshot_editorial_pack/F06_natcom_2026_precedent_bridge.png",
        (860, 1045, 1620, 1500),
        "Nature Communications bridge",
        "Tier framing and precedent for the target journal.",
        GOLD,
    )
    card(
        img,
        HUB / "assets/paper2_max_impact_unlock/F02_impact_unlock_ladder.png",
        (1650, 1045, 2530, 1500),
        "Claim ladder",
        "Unlock path from evidence to manuscript.",
        GOLD,
    )
    card(
        img,
        paths["g230_hotspot"],
        (70, 1565, 1260, 1880),
        "Hotspot contours",
        "Observed vs predicted contours over histology.",
        GOLD,
    )
    card(
        img,
        PATHO / "summary_figs/F8_per_case_predictions.png",
        (1340, 1565, 2530, 1880),
        "Per-case pathology readout",
        "Case-level readout across the internal series.",
        GREEN,
    )
    box(img, (70, 1890, 2530, 1946), (18, 23, 31), GOLD)
    put(img, "One-screen narrative: what the H&E sees, what the pathology confirms, what the reviewer can audit, and how it bridges to the journal target.", (110, 1930), 0.58, INK, 2, 2350)
    path = OUT / "F09_cover_executive_board.png"
    cv2.imwrite(str(path), img)
    return path


def draw_grand_cover_board(paths: dict[str, Path]) -> Path:
    img = np.full((2040, 2680, 3), BG, dtype=np.uint8)
    title(
        img,
        "Grand cover",
        "Paper 2 Grand Synthesis",
        "The single strongest opening frame: external H&E, internal pathology, validation, and the journal bridge, with almost no visual clutter.",
    )
    card(
        img,
        paths["g230_he"],
        (70, 315, 1320, 1080),
        "External H&E anchor",
        "GSE230424 histology with image-derived DM1/RAI signal.",
        TEAL,
    )
    card(
        img,
        PATHO / "summary_figs/F1_trajectory.png",
        (1360, 315, 2610, 700),
        "Internal pathology trajectory",
        "Local readout across the internal pathology series.",
        GREEN,
    )
    card(
        img,
        PATHO / "summary_figs/F2_cross_platform_replicate.png",
        (1360, 740, 2610, 1080),
        "Internal replicate support",
        "Cross-platform replicate consistency.",
        BLUE,
    )
    card(
        img,
        HUB / "assets/paper2_external_validation_ready/F05_pathologist_hotspot_audit_sheet.png",
        (70, 1140, 855, 1600),
        "Validation audit",
        "Pathologist hotspot audit and locked validation lane.",
        RED,
    )
    card(
        img,
        HUB / "assets/paper2_moonshot_editorial_pack/F06_natcom_2026_precedent_bridge.png",
        (920, 1140, 1705, 1600),
        "Journal bridge",
        "Nature Communications precedent and tier framing.",
        GOLD,
    )
    card(
        img,
        HUB / "assets/paper2_max_impact_unlock/F02_impact_unlock_ladder.png",
        (1770, 1140, 2610, 1600),
        "Claim ladder",
        "Unlock path from evidence to manuscript.",
        GOLD,
    )
    card(
        img,
        paths["g230_hotspot"],
        (70, 1655, 1320, 1965),
        "Hotspot contours",
        "Observed vs predicted contours over histology.",
        GOLD,
    )
    card(
        img,
        PATHO / "summary_figs/F8_per_case_predictions.png",
        (1360, 1655, 2610, 1965),
        "Per-case readout",
        "Case-level pathology spread across the internal series.",
        GREEN,
    )
    path = OUT / "F10_grand_cover_board.png"
    cv2.imwrite(str(path), img)
    return path


def copy_for_asset(src: Path, dest_name: str) -> Path:
    dest = OUT / dest_name
    shutil.copy2(src, dest)
    return dest


def write_manifest(figs: list[Path], copied: list[tuple[Path, str, str, str]]) -> None:
    entries: list[dict[str, str]] = []
    for i, fig in enumerate(figs, start=1):
        entries.append(
            {
                "id": f"he{i:03d}",
                "group": "H&E fused boards",
                "title": fig.stem.replace("_", " "),
                "caption": "OpenCV-built H&E-first board for fast visual review.",
                "src": f"assets/{ASSET}/{fig.name}",
            }
        )
    for i, (path, group, title_text, caption) in enumerate(copied, start=len(entries) + 1):
        entries.append(
            {
                "id": f"he{i:03d}",
                "group": group,
                "title": title_text,
                "caption": caption,
                "src": f"assets/{ASSET}/{path.name}",
            }
        )
    payload = json.dumps(entries, indent=2, ensure_ascii=False)
    (OUT / "overlay_manifest.json").write_text(payload + "\n", encoding="utf-8")
    (OUT / "overlay_manifest.js").write_text("window.paper2HeOverlayFigures = " + payload + ";\n", encoding="utf-8")
    summary = {
        "created": "2026-05-11",
        "purpose": "H&E-first overlay viewer assets for Paper 2 CV2 page with added validation/editorial frames",
        "n_overlay_entries": len(entries),
        "figures": [p.name for p in figs],
        "note": "Internal hospital H&E data are not included; this uses existing public/current Paper 2 H&E/spatial overlay assets plus review-frame boards and an internal pathology anchor already in the repo.",
    }
    (OUT / "PAPER2_HE_OVERLAY_VIEWER_ASSETS.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    md = [
        "# Paper 2 H&E overlay viewer assets",
        "",
        "Purpose: make the CV2 page easier to read by putting H&E-fused overlays, review frames, and an internal pathology anchor before the 81 separated source panels.",
        "",
        "Generated boards:",
        *[f"- `{p.name}`" for p in figs],
        "",
        "Boundary: no new internal K2/Bundang H&E data are used here.",
    ]
    (OUT / "SUMMARY.md").write_text("\n".join(md) + "\n", encoding="utf-8")


def publish() -> None:
    ASSET_LOCAL.mkdir(parents=True, exist_ok=True)
    ASSET_LIVE.mkdir(parents=True, exist_ok=True)
    for path in OUT.iterdir():
        if path.is_file():
            shutil.copy2(path, ASSET_LOCAL / path.name)
            shutil.copy2(path, ASSET_LIVE / path.name)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    paths = source_paths()
    slides = slide_paths()
    generated = [
        draw_raw_source_wall(),
        draw_raw_he_gallery(),
        draw_real_tile_atlas(),
        draw_internal_pathology_full_atlas(),
        draw_internal_pathology_atlas(),
        draw_grand_cover_board(paths),
        draw_cover_executive_board(paths),
        draw_cover_synthesis_board(paths),
        draw_executive_decision_board(paths),
        draw_front_door_max_impact(paths),
        draw_storyboard(paths),
        draw_g250_contact_sheet(slides),
        draw_g230_pair(paths),
        draw_legend(),
        draw_internal_pathology_anchor(),
    ]
    copied: list[tuple[Path, str, str, str]] = []
    extra_dirs = [
        (
            P2 / "analysis_supp/paper2_submission_overview_2026_05_10",
            "Submission overview",
            [
                ("F00_paper2_whole_story_overview.png", "Paper 2 whole story overview", "Master figure-flow map for manuscript organization."),
                ("F01_paper2_data_use_map.png", "Paper 2 data use map", "Data provenance and how each dataset is used."),
                ("F02_paper2_figure_table_roadmap.png", "Paper 2 figure table roadmap", "Figure and table roadmap for the submission."),
                ("F03_paper2_representative_figure_bundle.png", "Paper 2 representative figure bundle", "Representative visual bundle for editor review."),
            ],
        ),
        (
            P2 / "analysis_supp/paper2_moonshot_editorial_pack_2026_05_10",
            "Editorial H&E frames",
            [
                ("F00_general_he_front_door_frame.png", "General H&E front door frame", "High-level H&E opener for the Paper 2 review deck."),
                ("F01_moonshot_graphical_abstract.png", "Moonshot graphical abstract", "Editorial overview frame for the image-DM1 story."),
                ("F02_moonshot_claim_ladder.png", "Moonshot claim ladder", "Claim progression frame for reviewer-facing triage."),
                ("F03_editor_scorecard.png", "Editor scorecard", "High-level scorecard for submission framing."),
                ("F04_presub_attachment_stack.png", "Pre-sub attachment stack", "Attachment map for pre-submission packaging."),
                ("F05_k2_bundang_unlock_board.png", "K2 Bundang unlock board", "Internal H&E unlock framing board."),
                ("F06_natcom_2026_precedent_bridge.png", "Nature Communications precedent bridge", "Precedent bridge for the target journal tier."),
            ],
        ),
        (
            P2 / "analysis_supp/paper2_external_validation_ready_2026_05_10",
            "External validation frames",
            [
                ("F01_external_validation_request_packet.png", "External validation request packet", "Validation request frame for the public-review workflow."),
                ("F02_sample_manifest_schema.png", "Sample manifest schema", "Manifest schema for locked sample handling."),
                ("F03_locked_analysis_workflow.png", "Locked analysis workflow", "Workflow frame for locked and auditable runs."),
                ("F04_go_nogo_threshold_board.png", "Go/no-go threshold board", "Threshold frame for go/no-go decisions."),
                ("F05_pathologist_hotspot_audit_sheet.png", "Pathologist hotspot audit sheet", "Audit sheet for hotspot review."),
                ("F06_five_day_execution_board.png", "Five day execution board", "Short-horizon execution frame for the validation lane."),
            ],
        ),
        (
            P2 / "analysis_supp/paper2_max_impact_unlock_2026_05_10",
            "Impact unlock frames",
            [
                ("F01_max_impact_control_room.png", "Max impact control room", "Top-level control room for the impact package."),
                ("F02_impact_unlock_ladder.png", "Impact unlock ladder", "Stepwise ladder for unlocking the strongest claim."),
                ("F03_locked_external_validation_protocol.png", "Locked external validation protocol", "Protocol frame for locked external validation."),
                ("F04_hest_benchmark_bridge.png", "HEST benchmark bridge", "Bridge frame to the HEST-style benchmark story."),
                ("F05_editor_risk_buydown.png", "Editor risk buydown", "Risk buydown board for submission strategy."),
                ("F06_editor_presub_packet.png", "Editor pre-sub packet", "Pre-sub packet framing for the editor-facing version."),
            ],
        ),
    ]
    for folder, group, items in extra_dirs:
        for fname, title_text, caption in items:
            src = folder / fname
            dest = copy_for_asset(src, src.name)
            copied.append((dest, group, title_text, caption))
    selected = [
        (paths["g230_he"], "External thyroid H&E overlays", "GSE230424 H&E + DM1/RAI overlay mosaic", "Existing H&E-centered GSE230424 overlay mosaic."),
        (paths["g230_hotspot"], "External thyroid H&E overlays", "GSE230424 H&E + hotspot contours", "Existing H&E-centered hotspot contour board."),
        (paths["g250_mosaic"], "GSE250521 H&E overlays", "GSE250521 16-slide H&E + UNI-DM1 mosaic", "Existing 16-slide H&E overlay mosaic."),
    ]
    for src, group, title_text, caption in selected:
        dest = copy_for_asset(src, src.name)
        copied.append((dest, group, title_text, caption))
    for src in slides:
        label = re.sub(r"^phase1_gse250521__uni_dm1_overlay_", "", src.stem).replace("_", " ")
        dest = copy_for_asset(src, src.name)
        copied.append((dest, "GSE250521 individual H&E overlays", label, "Single-slide H&E + UNI-DM1 overlay."))
    write_manifest(generated, copied)
    publish()
    print(json.dumps({"out": str(OUT), "asset": ASSET, "generated": [p.name for p in generated], "slides": len(slides), "copied": len(copied)}, indent=2))


if __name__ == "__main__":
    main()
