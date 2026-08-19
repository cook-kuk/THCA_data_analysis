#!/usr/bin/env python3
"""OpenCV visual summary pack for Paper 2 image-to-spatial-DM1 evidence."""
from __future__ import annotations

import gc
import gzip
import json
import math
from pathlib import Path

import cv2
import numpy as np
import pandas as pd


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
P2 = ROOT / "project/results/p2_image_dm1_v2_foundation_clam_2026_05_07"
SUPP = P2 / "analysis_supp"
RAW230 = ROOT / "project/data/external/GSE230424/raw"
OUT = SUPP / "paper2_cv2_visual_summary_2026_05_10"

G250 = SUPP / "path2space_inspired_reanalysis_2026_05_09"
G250_SPEC = SUPP / "path2space_gse250521_random_module_specificity_2026_05_09"
G250_STAGE = SUPP / "path2space_stage_generalization_controls_2026_05_09"
G230 = SUPP / "gse230424_pathology_thyroid_axis_2026_05_09"
G230_RESID = G230 / "residual_target_controls"
DECILE = SUPP / "path2space_decile_dose_response_2026_05_09"
DOWNSAMPLE = SUPP / "path2space_downsample_robustness_2026_05_09"
HOTSPOT = SUPP / "path2space_hotspot_concordance_2026_05_09"
BLOCKS = SUPP / "path2space_spatial_block_controls_2026_05_09"
SMOOTHNESS = G230 / "spatial_autocorr_specificity"
AITD = SUPP / "gse248205_pathology_aitd_axis_2026_05_09"
ASSETS_EXISTING = ROOT / "project/papers_hub_2026_05_04/assets/paper2_path2space"
G250_DOMAIN_RHO = 0.581


BG = (17, 22, 31)
PANEL = (27, 34, 45)
PANEL2 = (36, 45, 60)
INK = (232, 237, 243)
MUTED = (166, 176, 190)
GOLD = (94, 178, 214)
TEAL = (196, 213, 99)
BLUE = (255, 166, 87)
RED = (114, 123, 255)
GREEN = (157, 211, 53)
WHITE = (245, 246, 248)


def read_json(path: Path) -> dict:
    with path.open() as handle:
        return json.load(handle)


def put_text(
    img: np.ndarray,
    text: str,
    xy: tuple[int, int],
    scale: float = 0.55,
    color: tuple[int, int, int] = INK,
    thickness: int = 1,
    max_width: int | None = None,
    line_gap: int = 8,
) -> int:
    x, y = xy
    font = cv2.FONT_HERSHEY_SIMPLEX
    lines: list[str] = []
    if max_width is None:
        lines = text.split("\n")
    else:
        for raw_line in text.split("\n"):
            words = raw_line.split()
            line = ""
            for word in words:
                trial = word if not line else f"{line} {word}"
                width = cv2.getTextSize(trial, font, scale, thickness)[0][0]
                if width <= max_width or not line:
                    line = trial
                else:
                    lines.append(line)
                    line = word
            lines.append(line)
    _, line_h = cv2.getTextSize("Ag", font, scale, thickness)[0]
    step = line_h + line_gap
    for i, line in enumerate(lines):
        cv2.putText(img, line, (x, y + i * step), font, scale, color, thickness, cv2.LINE_AA)
    return y + max(len(lines), 1) * step


def rect(img: np.ndarray, xyxy: tuple[int, int, int, int], color: tuple[int, int, int], radius: int = 12, fill: bool = True) -> None:
    x0, y0, x1, y1 = xyxy
    if not fill:
        cv2.rectangle(img, (x0, y0), (x1, y1), color, 2, cv2.LINE_AA)
        return
    cv2.rectangle(img, (x0 + radius, y0), (x1 - radius, y1), color, -1)
    cv2.rectangle(img, (x0, y0 + radius), (x1, y1 - radius), color, -1)
    for cx, cy in [(x0 + radius, y0 + radius), (x1 - radius, y0 + radius), (x0 + radius, y1 - radius), (x1 - radius, y1 - radius)]:
        cv2.circle(img, (cx, cy), radius, color, -1, cv2.LINE_AA)


def draw_metric_card(
    img: np.ndarray,
    box: tuple[int, int, int, int],
    label: str,
    value: str,
    note: str,
    accent: tuple[int, int, int],
    status: str = "support",
) -> None:
    x0, y0, x1, y1 = box
    rect(img, box, PANEL, radius=14)
    cv2.rectangle(img, (x0, y0), (x1, y0 + 7), accent, -1)
    put_text(img, label.upper(), (x0 + 22, y0 + 34), 0.46, MUTED, 1, x1 - x0 - 44)
    put_text(img, value, (x0 + 22, y0 + 88), 1.12, accent, 2, x1 - x0 - 44)
    put_text(img, note, (x0 + 22, y0 + 126), 0.48, INK, 1, x1 - x0 - 44)
    chip_color = GREEN if status == "support" else (87, 204, 255) if status == "caveat" else RED
    cv2.circle(img, (x1 - 28, y0 + 30), 7, chip_color, -1, cv2.LINE_AA)


def fit_image(src: np.ndarray, size: tuple[int, int], bg: tuple[int, int, int] = WHITE) -> np.ndarray:
    tw, th = size
    h, w = src.shape[:2]
    scale = min(tw / max(w, 1), th / max(h, 1))
    nw, nh = max(1, int(w * scale)), max(1, int(h * scale))
    resized = cv2.resize(src, (nw, nh), interpolation=cv2.INTER_AREA)
    canvas = np.full((th, tw, 3), bg, dtype=np.uint8)
    x0 = (tw - nw) // 2
    y0 = (th - nh) // 2
    canvas[y0 : y0 + nh, x0 : x0 + nw] = resized
    return canvas


def value_color(v: float, lim: float) -> tuple[int, int, int]:
    if not np.isfinite(v) or lim <= 0:
        return (180, 180, 180)
    x = float(np.clip(v / lim, -1, 1))
    if x >= 0:
        t = x
        return (
            int((1 - t) * 245 + t * 60),
            int((1 - t) * 245 + t * 80),
            int((1 - t) * 245 + t * 235),
        )
    t = -x
    return (
        int((1 - t) * 245 + t * 230),
        int((1 - t) * 245 + t * 110),
        int((1 - t) * 245 + t * 45),
    )


def draw_colorbar(img: np.ndarray, x: int, y: int, w: int, h: int, lim: float, label: str) -> None:
    for i in range(w):
        v = (i / max(w - 1, 1) * 2 - 1) * lim
        cv2.line(img, (x + i, y), (x + i, y + h), value_color(v, lim), 1)
    cv2.rectangle(img, (x, y), (x + w, y + h), (35, 42, 52), 1, cv2.LINE_AA)
    put_text(img, f"-{lim:.1f}", (x, y + h + 24), 0.42, MUTED)
    put_text(img, "0", (x + w // 2 - 8, y + h + 24), 0.42, MUTED)
    put_text(img, f"+{lim:.1f}", (x + w - 48, y + h + 24), 0.42, MUTED)
    put_text(img, label, (x, y - 8), 0.45, MUTED)


def load_metrics() -> dict:
    uni = read_json(SUPP / "audit_uni_loto/UNI_LOTO_SUMMARY.json")
    g250 = read_json(G250 / "PATH2SPACE_REANALYSIS_SUMMARY.json")
    g250_spec = read_json(G250_SPEC / "GSE250521_RANDOM_MODULE_SPECIFICITY_SUMMARY.json")
    g250_stage = read_json(G250_STAGE / "GSE250521_STAGE_GENERALIZATION_SUMMARY.json")
    g230 = read_json(G230 / "GSE230424_PATHOLOGY_THYROID_SUMMARY.json")
    smooth = read_json(SMOOTHNESS / "GSE230424_SPATIAL_AUTOCORR_SPECIFICITY_SUMMARY.json")
    decile = read_json(DECILE / "PATH2SPACE_DECILE_DOSE_RESPONSE_SUMMARY.json")
    down = read_json(DOWNSAMPLE / "PATH2SPACE_DOWNSAMPLE_ROBUSTNESS_SUMMARY.json")
    return {
        "uni": uni,
        "g250": g250,
        "g250_spec": g250_spec,
        "g250_stage": g250_stage,
        "g230": g230,
        "smooth": smooth,
        "decile": decile,
        "down": down,
    }


def make_evidence_scorecard(metrics: dict) -> None:
    img = np.full((1320, 1900, 3), BG, dtype=np.uint8)
    put_text(img, "PAPER 2 CV2 VISUAL SUMMARY", (70, 80), 0.62, GOLD, 2)
    put_text(img, "H&E -> DM1 image classifier -> spatial transcriptomics projection", (70, 142), 1.08, INK, 2, 1260)
    put_text(
        img,
        "Generated from local TSV/JSON outputs and OpenCV overlays. Claim boundary: supports image-to-spatial-RNA Paper 2; does not rescue Paper 1 MAPK mechanism.",
        (70, 190),
        0.55,
        MUTED,
        1,
        1380,
    )

    steps = [
        ("TCGA H&E", "UNI/CLAM\nLOTO audit"),
        ("GSE250521", "Path2Space-style\nspatial RNA"),
        ("GSE230424", "external thyroid\nVisium H&E"),
        ("Controls", "decile, bootstrap,\nblocks, downsample"),
        ("Boundary", "QC/smoothness\ncaveats kept"),
    ]
    sx, sy, sw, sh, gap = 72, 255, 315, 112, 38
    for i, (head, sub) in enumerate(steps):
        x = sx + i * (sw + gap)
        rect(img, (x, sy, x + sw, sy + sh), PANEL2, radius=16)
        put_text(img, head, (x + 20, sy + 38), 0.67, INK, 2)
        put_text(img, sub, (x + 20, sy + 72), 0.47, MUTED, 1, sw - 40)
        if i < len(steps) - 1:
            cv2.arrowedLine(img, (x + sw + 6, sy + sh // 2), (x + sw + gap - 8, sy + sh // 2), GOLD, 2, cv2.LINE_AA, tipLength=0.35)

    cards = [
        (
            "Image-DM1 audit",
            f"AUC {metrics['uni']['pooled_overall_auc']:.3f}",
            f"{int(metrics['uni']['n_slides'])} slides; RAS-like AUC {metrics['uni']['pooled_ras_like_auc']:.3f}",
            GREEN,
            "support",
        ),
        (
            "GSE250521 spatial",
            f"rho {metrics['g250']['dm1_smoothed_rho']:.3f}",
            f"smoothed LOSO; slide-centered {metrics['g250_stage']['dm1_existing_slide_centered_rho']:.3f}",
            TEAL,
            "support",
        ),
        (
            "GSE250521 residual specificity",
            f"p {metrics['g250_spec']['actual_residual_empirical_p']:.3f}",
            f"coord+QC residual rho {metrics['g250_spec']['actual_residual_rho']:.3f}; raw smoothness caveat",
            GOLD,
            "caveat",
        ),
        (
            "GSE230424 external",
            f"rho {metrics['g230']['top_he_sample_centered_rho']:.3f}",
            f"{int(metrics['g230']['n_spots'])} spots / {int(metrics['g230']['n_samples'])} samples; residual {metrics['g230']['top_coord_qc_residual_rho']:.3f}",
            GREEN,
            "support",
        ),
        (
            "Decile dose-response",
            f"{metrics['decile']['gse230424_raw_decile_spearman']:.3f}",
            f"G230 raw Spearman; G250 top-bottom {metrics['decile']['gse250521_raw_top_minus_bottom']:.3f}",
            BLUE,
            "support",
        ),
        (
            "Downsample robustness",
            f"rho {metrics['down']['gse230424_raw_rho_median']:.3f}",
            "1000 repeats at 200 spots/sample; residual median rho "
            f"{metrics['down']['gse230424_resid_rho_median']:.3f}",
            TEAL,
            "support",
        ),
        (
            "Smoothness specificity",
            "FAIL",
            f"G230 raw adjusted residual {metrics['smooth']['raw_smoothness_adjusted_residual']:.3f}; disclose caveat",
            RED,
            "fail",
        ),
        (
            "Safe disposition",
            "P2 use",
            "Strong visual/spatial support with explicit QC and smoothness boundaries.",
            GOLD,
            "caveat",
        ),
    ]
    x0, y0 = 72, 430
    cw, ch, cg = 420, 210, 32
    for idx, card in enumerate(cards):
        row, col = divmod(idx, 4)
        x = x0 + col * (cw + cg)
        y = y0 + row * (ch + 34)
        draw_metric_card(img, (x, y, x + cw, y + ch), *card)

    rect(img, (70, 930, 1830, 1210), PANEL, radius=18)
    put_text(img, "Reviewer-facing message", (100, 985), 0.74, GOLD, 2)
    bullets = [
        "1. Pathology carries a measurable local DM1/RAI spatial signal.",
        "2. The signal is visible in GSE250521 and externally supported in GSE230424.",
        "3. Threshold, bootstrap, spatial-block, and spot-budget controls support robustness.",
        "4. QC and spatial smoothness are real confounds, so molecular specificity is framed conservatively.",
    ]
    y = 1035
    for bullet in bullets:
        y = put_text(img, bullet, (108, y), 0.62, INK, 1, 1650, 12) + 4

    cv2.imwrite(str(OUT / "fig01_paper2_cv2_evidence_scorecard.png"), img)


def sample_prefixes() -> dict[str, str]:
    out = {}
    for p in RAW230.glob("*_matrix.mtx.gz"):
        prefix = p.name.replace("_matrix.mtx.gz", "")
        sample = prefix.split("_")[-1]
        out[sample] = prefix
    return out


def read_positions(prefix: str) -> pd.DataFrame:
    pos = pd.read_csv(
        RAW230 / f"{prefix}_tissue_positions_list.csv.gz",
        header=None,
        names=["barcode", "in_tissue", "array_row", "array_col", "pxl_row_fullres", "pxl_col_fullres"],
    )
    pos["barcode"] = pos["barcode"].astype(str)
    return pos


def load_he_image(prefix: str) -> np.ndarray:
    with gzip.open(RAW230 / f"{prefix}_HE.jpg.gz", "rb") as handle:
        data = handle.read()
    return cv2.imdecode(np.frombuffer(data, np.uint8), cv2.IMREAD_COLOR)


def he_overlay_tile(
    image: np.ndarray,
    sub: pd.DataFrame,
    value_col: str,
    title: str,
    lim: float,
    tile_w: int = 420,
    tile_h: int = 420,
) -> np.ndarray:
    xs = sub["pxl_col_fullres"].to_numpy(dtype=float)
    ys = sub["pxl_row_fullres"].to_numpy(dtype=float)
    margin = 650
    h, w = image.shape[:2]
    x0 = max(0, int(np.nanmin(xs) - margin))
    x1 = min(w, int(np.nanmax(xs) + margin))
    y0 = max(0, int(np.nanmin(ys) - margin))
    y1 = min(h, int(np.nanmax(ys) + margin))
    crop = image[y0:y1, x0:x1].copy()
    base = cv2.resize(crop, (tile_w, tile_h), interpolation=cv2.INTER_AREA)
    overlay = base.copy()
    sx = tile_w / max(x1 - x0, 1)
    sy = tile_h / max(y1 - y0, 1)
    vals = sub[value_col].to_numpy(dtype=float)
    order = np.argsort(vals)
    for idx in order:
        x = int((xs[idx] - x0) * sx)
        y = int((ys[idx] - y0) * sy)
        color = value_color(vals[idx], lim)
        cv2.circle(overlay, (x, y), 4, color, -1, cv2.LINE_AA)
    out = cv2.addWeighted(overlay, 0.72, base, 0.28, 0)
    cv2.rectangle(out, (0, 0), (tile_w - 1, tile_h - 1), (35, 42, 52), 2, cv2.LINE_AA)
    cv2.rectangle(out, (0, 0), (tile_w, 42), (15, 18, 24), -1)
    put_text(out, title, (14, 28), 0.58, WHITE, 2, tile_w - 28)
    return out


def exact_top_boolean(values: np.ndarray, frac: float = 0.10) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    keep = np.isfinite(values)
    valid = np.where(keep)[0]
    mask = np.zeros(len(values), dtype=bool)
    if len(valid) == 0:
        return mask
    k = max(1, int(np.ceil(len(valid) * frac)))
    order = np.argsort(values[valid], kind="mergesort")
    mask[valid[order[-k:]]] = True
    return mask


def hotspot_contour_tile(
    image: np.ndarray,
    sub: pd.DataFrame,
    title: str,
    tile_w: int = 520,
    tile_h: int = 420,
) -> tuple[np.ndarray, dict[str, float]]:
    xs = sub["pxl_col_fullres"].to_numpy(dtype=float)
    ys = sub["pxl_row_fullres"].to_numpy(dtype=float)
    margin = 650
    h, w = image.shape[:2]
    x0 = max(0, int(np.nanmin(xs) - margin))
    x1 = min(w, int(np.nanmax(xs) + margin))
    y0 = max(0, int(np.nanmin(ys) - margin))
    y1 = min(h, int(np.nanmax(ys) + margin))
    crop = image[y0:y1, x0:x1].copy()
    base = cv2.resize(crop, (tile_w, tile_h), interpolation=cv2.INTER_AREA)
    overlay = base.copy()
    sx = tile_w / max(x1 - x0, 1)
    sy = tile_h / max(y1 - y0, 1)
    obs_top = exact_top_boolean(sub["obs_DM1_low_RAI_score_smooth8"].to_numpy(dtype=float))
    pred_top = exact_top_boolean(sub["pred_DM1_low_RAI_score_HE_tile_features"].to_numpy(dtype=float))
    overlap = obs_top & pred_top
    pred_only = pred_top & ~obs_top
    obs_only = obs_top & ~pred_top

    mask_obs = np.zeros((tile_h, tile_w), dtype=np.uint8)
    mask_pred = np.zeros((tile_h, tile_w), dtype=np.uint8)
    pts = []
    for i in range(len(sub)):
        x = int((xs[i] - x0) * sx)
        y = int((ys[i] - y0) * sy)
        pts.append((x, y))
        if obs_top[i]:
            cv2.circle(mask_obs, (x, y), 10, 255, -1, cv2.LINE_AA)
        if pred_top[i]:
            cv2.circle(mask_pred, (x, y), 10, 255, -1, cv2.LINE_AA)
    kernel = np.ones((13, 13), np.uint8)
    mask_obs = cv2.morphologyEx(mask_obs, cv2.MORPH_CLOSE, kernel)
    mask_pred = cv2.morphologyEx(mask_pred, cv2.MORPH_CLOSE, kernel)
    for mask, color in [(mask_obs, (70, 120, 255)), (mask_pred, (70, 230, 210))]:
        contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        cv2.drawContours(overlay, contours, -1, color, 2, cv2.LINE_AA)
    for i, (x, y) in enumerate(pts):
        if overlap[i]:
            cv2.circle(overlay, (x, y), 5, GREEN, -1, cv2.LINE_AA)
        elif pred_only[i]:
            cv2.circle(overlay, (x, y), 4, (70, 230, 210), -1, cv2.LINE_AA)
        elif obs_only[i]:
            cv2.circle(overlay, (x, y), 4, (70, 120, 255), -1, cv2.LINE_AA)
    out = cv2.addWeighted(overlay, 0.78, base, 0.22, 0)
    cv2.rectangle(out, (0, 0), (tile_w - 1, tile_h - 1), (35, 42, 52), 2, cv2.LINE_AA)
    cv2.rectangle(out, (0, 0), (tile_w, 50), (15, 18, 24), -1)
    put_text(out, title, (14, 32), 0.58, WHITE, 2, tile_w - 28)
    precision = float(overlap.sum() / max(pred_top.sum(), 1))
    lift = precision / 0.10
    metrics = {
        "precision": precision,
        "lift": lift,
        "pred_top": int(pred_top.sum()),
        "obs_top": int(obs_top.sum()),
        "overlap": int(overlap.sum()),
    }
    return out, metrics


def make_gse230424_he_overlay(metrics: dict) -> None:
    pred = pd.read_csv(G230 / "gse230424_pathology_predictions.tsv.gz", sep="\t")
    prefixes = sample_prefixes()
    dfs = []
    for sample, prefix in prefixes.items():
        pos = read_positions(prefix)
        sub = pred[pred["sample"] == sample].merge(
            pos[["barcode", "pxl_row_fullres", "pxl_col_fullres"]],
            on="barcode",
            how="inner",
        )
        dfs.append(sub)
    merged = pd.concat(dfs, ignore_index=True)
    vals = np.concatenate(
        [
            merged["obs_DM1_low_RAI_score_smooth8"].to_numpy(dtype=float),
            merged["pred_DM1_low_RAI_score_HE_tile_features"].to_numpy(dtype=float),
        ]
    )
    lim = float(np.nanpercentile(np.abs(vals[np.isfinite(vals)]), 98))
    lim = max(lim, 0.5)
    canvas = np.full((2380, 1860, 3), BG, dtype=np.uint8)
    put_text(canvas, "GSE230424 external thyroid Visium: H&E spot overlays", (60, 70), 0.95, INK, 2)
    put_text(
        canvas,
        f"Sample-centered H&E rho {metrics['g230']['top_he_sample_centered_rho']:.3f}; strict residual-target rho 0.232; downsample median rho {metrics['down']['gse230424_raw_rho_median']:.3f}.",
        (60, 118),
        0.56,
        MUTED,
        1,
        1550,
    )
    put_text(canvas, "Observed spatial RNA", (300, 182), 0.66, GOLD, 2)
    put_text(canvas, "Predicted from H&E", (1035, 182), 0.66, GOLD, 2)
    draw_colorbar(canvas, 1360, 138, 300, 18, lim, "DM1/low-RAI score")

    row_y = 225
    prefixes_sorted = sorted(prefixes.items())
    disease = {"P1": "PTC+HT", "P2": "PTC+HT", "P3": "HT", "P4": "HT"}
    for r, (sample, prefix) in enumerate(prefixes_sorted):
        image = load_he_image(prefix)
        sub = merged[merged["sample"] == sample].copy()
        obs = he_overlay_tile(image, sub, "obs_DM1_low_RAI_score_smooth8", f"{sample} observed", lim)
        pred_tile = he_overlay_tile(image, sub, "pred_DM1_low_RAI_score_HE_tile_features", f"{sample} predicted", lim)
        del image
        gc.collect()
        y = row_y + r * 520
        rect(canvas, (42, y - 18, 1808, y + 462), PANEL, radius=14)
        put_text(canvas, f"{sample}", (70, y + 62), 0.86, INK, 2)
        put_text(canvas, disease.get(sample, "NA"), (70, y + 105), 0.52, GOLD, 1)
        rho, _, _ = safe_spearman(
            center_by_group(sub["obs_DM1_low_RAI_score_smooth8"].to_numpy(dtype=float), sub["sample"].to_numpy(str)),
            center_by_group(sub["pred_DM1_low_RAI_score_HE_tile_features"].to_numpy(dtype=float), sub["sample"].to_numpy(str)),
        )
        put_text(canvas, f"n={len(sub):,}\nrho={rho:.3f}", (70, y + 154), 0.50, MUTED, 1)
        canvas[y : y + 420, 255 : 255 + 420] = obs
        canvas[y : y + 420, 980 : 980 + 420] = pred_tile
        put_text(canvas, "same H&E crop, spots recolored by RNA label or H&E prediction", (1425, y + 68), 0.47, MUTED, 1, 310)
    cv2.imwrite(str(OUT / "fig02_gse230424_cv2_he_overlay_mosaic.png"), canvas)


def make_gse230424_hotspot_contours(metrics: dict) -> None:
    pred = pd.read_csv(G230 / "gse230424_pathology_predictions.tsv.gz", sep="\t")
    prefixes = sample_prefixes()
    rows = []
    canvas = np.full((1240, 1900, 3), BG, dtype=np.uint8)
    put_text(canvas, "GSE230424 hotspot contours on true H&E", (60, 74), 0.95, INK, 2)
    put_text(
        canvas,
        "Top-decile observed DM1/low-RAI contours are orange; H&E-predicted top-decile contours are cyan; direct overlap spots are green.",
        (60, 122),
        0.55,
        MUTED,
        1,
        1500,
    )
    put_text(canvas, "observed top-10%", (1190, 86), 0.45, (70, 120, 255), 1)
    put_text(canvas, "predicted top-10%", (1190, 116), 0.45, (70, 230, 210), 1)
    put_text(canvas, "overlap", (1190, 146), 0.45, GREEN, 1)
    disease = {"P1": "PTC+HT", "P2": "PTC+HT", "P3": "HT", "P4": "HT"}
    for i, (sample, prefix) in enumerate(sorted(prefixes.items())):
        pos = read_positions(prefix)
        sub = pred[pred["sample"] == sample].merge(
            pos[["barcode", "pxl_row_fullres", "pxl_col_fullres"]],
            on="barcode",
            how="inner",
        )
        image = load_he_image(prefix)
        tile, m = hotspot_contour_tile(image, sub, f"{sample} {disease.get(sample, '')}")
        del image
        gc.collect()
        row, col = divmod(i, 2)
        x = 60 + col * 900
        y = 190 + row * 500
        rect(canvas, (x, y, x + 845, y + 450), PANEL, radius=14)
        canvas[y + 16 : y + 16 + tile.shape[0], x + 16 : x + 16 + tile.shape[1]] = tile
        put_text(canvas, f"precision {m['precision']:.3f}", (x + 560, y + 80), 0.55, INK, 1)
        put_text(canvas, f"lift {m['lift']:.2f}x", (x + 560, y + 118), 0.78, GOLD, 2)
        put_text(canvas, f"overlap {m['overlap']}/{m['pred_top']}", (x + 560, y + 162), 0.50, MUTED, 1)
        rows.append({"sample": sample, **m})
    pd.DataFrame(rows).to_csv(OUT / "paper2_cv2_gse230424_hotspot_contours.tsv", sep="\t", index=False)
    cv2.imwrite(str(OUT / "fig05_gse230424_cv2_hotspot_contours.png"), canvas)


def safe_spearman(x, y) -> tuple[float, float, int]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    if keep.sum() < 5 or np.nanstd(x[keep]) == 0 or np.nanstd(y[keep]) == 0:
        return np.nan, np.nan, int(keep.sum())
    from scipy import stats

    rho, p = stats.spearmanr(x[keep], y[keep])
    return float(rho), float(p), int(keep.sum())


def center_by_group(values: np.ndarray, groups: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    groups = np.asarray(groups).astype(str)
    out = np.full(len(values), np.nan)
    for group in pd.unique(groups):
        idx = groups == group
        out[idx] = values[idx] - np.nanmean(values[idx])
    return out


def draw_spatial_points(tile: np.ndarray, sub: pd.DataFrame, value_col: str, box: tuple[int, int, int, int], lim: float) -> None:
    x0, y0, x1, y1 = box
    xs = sub["array_col"].to_numpy(dtype=float)
    ys = sub["array_row"].to_numpy(dtype=float)
    vals = sub[value_col].to_numpy(dtype=float)
    xmin, xmax = np.nanmin(xs), np.nanmax(xs)
    ymin, ymax = np.nanmin(ys), np.nanmax(ys)
    order = np.argsort(vals)
    for idx in order:
        px = int(x0 + (xs[idx] - xmin) / max(xmax - xmin, 1) * (x1 - x0 - 1))
        py = int(y0 + (ys[idx] - ymin) / max(ymax - ymin, 1) * (y1 - y0 - 1))
        cv2.circle(tile, (px, py), 3, value_color(vals[idx], lim), -1, cv2.LINE_AA)
    cv2.rectangle(tile, (x0, y0), (x1, y1), (56, 65, 80), 1, cv2.LINE_AA)


def make_gse250521_uni_mosaic(metrics: dict) -> None:
    pred = pd.read_csv(G250 / "path2space_spot_predictions.tsv.gz", sep="\t")
    slide_sum = pd.read_csv(G250 / "path2space_loso_slide_summary.tsv", sep="\t")
    slide_rho = slide_sum[(slide_sum["target"] == "DM1_like_score") & (slide_sum["mode"] == "smooth8")].set_index("sample_id")["rho"].to_dict()
    vals = np.concatenate(
        [
            pred["obs_DM1_like_score_smooth8"].to_numpy(dtype=float),
            pred["pred_DM1_like_score_smooth8"].to_numpy(dtype=float),
        ]
    )
    lim = max(0.5, float(np.nanpercentile(np.abs(vals[np.isfinite(vals)]), 98)))
    canvas = np.full((1320, 1900, 3), BG, dtype=np.uint8)
    put_text(canvas, "GSE250521 UNI spatial mosaic: observed vs predicted DM1/RAI", (60, 72), 0.92, INK, 2)
    put_text(
        canvas,
        f"Smoothed LOSO rho {metrics['g250']['dm1_smoothed_rho']:.3f}; slide-centered rho {metrics['g250_stage']['dm1_existing_slide_centered_rho']:.3f}; domain rho {G250_DOMAIN_RHO:.3f}.",
        (60, 118),
        0.54,
        MUTED,
        1,
        1550,
    )
    draw_colorbar(canvas, 1420, 78, 320, 18, lim, "DM1/RAI score")
    order = list(pred[["sample_id", "stage"]].drop_duplicates().sort_values(["stage", "sample_id"])["sample_id"])
    stage_order = {"PT": 0, "PTC": 1, "LPTC": 2, "ATC": 3}
    order = sorted(order, key=lambda s: (stage_order.get(str(pred[pred["sample_id"] == s]["stage"].iloc[0]), 9), s))
    x0, y0 = 60, 180
    tile_w, tile_h = 430, 260
    for i, sample in enumerate(order):
        row, col = divmod(i, 4)
        x = x0 + col * 455
        y = y0 + row * 275
        sub = pred[pred["sample_id"] == sample]
        stage = str(sub["stage"].iloc[0])
        rect(canvas, (x, y, x + tile_w, y + tile_h), PANEL, radius=12)
        short = sample.replace("GSM7980", "G")
        put_text(canvas, f"{stage}  {short}", (x + 16, y + 30), 0.48, INK, 1, tile_w - 32)
        put_text(canvas, f"rho {slide_rho.get(sample, np.nan):.2f}", (x + tile_w - 95, y + 30), 0.46, GOLD, 1)
        put_text(canvas, "obs", (x + 68, y + 58), 0.42, MUTED, 1)
        put_text(canvas, "pred", (x + 273, y + 58), 0.42, MUTED, 1)
        draw_spatial_points(canvas, sub, "obs_DM1_like_score_smooth8", (x + 18, y + 68, x + 198, y + 235), lim)
        draw_spatial_points(canvas, sub, "pred_DM1_like_score_smooth8", (x + 228, y + 68, x + 408, y + 235), lim)
    cv2.imwrite(str(OUT / "fig03_gse250521_cv2_uni_slide_mosaic.png"), canvas)


def make_result_gallery(metrics: dict) -> None:
    gallery = [
        ("Classifier audit", SUPP / "audit_uni_loto/fig_H_uni_loto.png"),
        ("Path2Space reanalysis", ASSETS_EXISTING / "fig_path2space_inspired_reanalysis.png"),
        ("GSE230424 external", ASSETS_EXISTING / "fig_gse230424_pathology_thyroid_axis.png"),
        ("GSE250521 specificity", ASSETS_EXISTING / "fig_gse250521_random_module_specificity.png"),
        ("Stage generalization", ASSETS_EXISTING / "fig_gse250521_stage_generalization_controls.png"),
        ("Bootstrap/downsample", ASSETS_EXISTING / "fig_path2space_downsample_robustness.png"),
        ("Decile dose-response", ASSETS_EXISTING / "fig_path2space_decile_dose_response.png"),
        ("Hotspot concordance", ASSETS_EXISTING / "fig_path2space_hotspot_concordance.png"),
        ("Hotspot contours", OUT / "fig05_gse230424_cv2_hotspot_contours.png"),
        ("Spatial-block control", ASSETS_EXISTING / "fig_path2space_spatial_block_controls.png"),
        ("Sample matrix", OUT / "fig07_paper2_cv2_sample_matrix.png"),
        ("Smoothness caveat", ASSETS_EXISTING / "fig_gse230424_spatial_autocorr_specificity.png"),
        ("GSE248205 negative", ASSETS_EXISTING / "fig_gse248205_pathology_aitd_axis.png"),
        ("Claim ladder", OUT / "fig06_paper2_cv2_claim_ladder.png"),
        ("Decision", OUT / "fig01_paper2_cv2_evidence_scorecard.png"),
    ]
    n_rows = math.ceil(len(gallery) / 3)
    canvas = np.full((190 + n_rows * 370 + 60, 1900, 3), BG, dtype=np.uint8)
    put_text(canvas, "Paper 2 figure map: what each panel is for", (60, 74), 0.95, INK, 2)
    put_text(canvas, "Green = use; gold = support with caveat; red = disclose boundary.", (60, 122), 0.55, MUTED, 1)
    tw, th = 550, 360
    for i, (label, path) in enumerate(gallery):
        row, col = divmod(i, 3)
        x = 60 + col * 605
        y = 180 + row * 370
        rect(canvas, (x, y, x + tw, y + th), PANEL, radius=12)
        src = cv2.imread(str(path), cv2.IMREAD_COLOR)
        if src is None:
            src = np.full((300, 450, 3), (250, 250, 250), dtype=np.uint8)
            put_text(src, "missing", (30, 160), 1.0, RED, 2)
        thumb = fit_image(src, (tw - 28, th - 78), bg=WHITE)
        canvas[y + 50 : y + 50 + thumb.shape[0], x + 14 : x + 14 + thumb.shape[1]] = thumb
        accent = GREEN
        if "caveat" in label.lower() or "negative" in label.lower():
            accent = RED
        elif "specificity" in label.lower() or "stage" in label.lower() or "bootstrap" in label.lower():
            accent = GOLD
        cv2.rectangle(canvas, (x, y), (x + tw, y + 8), accent, -1)
        put_text(canvas, f"{i + 1:02d}. {label}", (x + 16, y + 34), 0.52, INK, 1, tw - 32)
    cv2.imwrite(str(OUT / "fig04_paper2_cv2_result_gallery.png"), canvas)


def make_onepage_storyboard(metrics: dict) -> None:
    canvas = np.full((1600, 2400, 3), BG, dtype=np.uint8)
    put_text(canvas, "Paper 2 one-page visual storyboard", (70, 82), 1.12, INK, 2)
    put_text(
        canvas,
        "H&E image classifier -> spatial DM1/RAI projection -> external thyroid H&E support -> caveat-bounded claim",
        (70, 135),
        0.62,
        MUTED,
        1,
        1600,
    )
    stat_cards = [
        ("UNI LOTO AUC", f"{metrics['uni']['pooled_overall_auc']:.3f}", GREEN),
        ("G250 smoothed rho", f"{metrics['g250']['dm1_smoothed_rho']:.3f}", TEAL),
        ("G250 slide-centered", f"{metrics['g250_stage']['dm1_existing_slide_centered_rho']:.3f}", TEAL),
        ("G230 H&E rho", f"{metrics['g230']['top_he_sample_centered_rho']:.3f}", GREEN),
        ("G230 downsample", f"{metrics['down']['gse230424_raw_rho_median']:.3f}", GREEN),
        ("Residual downsample", f"{metrics['down']['gse230424_resid_rho_median']:.3f}", GOLD),
        ("Smoothness residual", f"{metrics['smooth']['raw_smoothness_adjusted_residual']:.3f}", RED),
    ]
    x = 70
    for label, value, color in stat_cards:
        rect(canvas, (x, 190, x + 300, 315), PANEL, radius=14)
        cv2.rectangle(canvas, (x, 190), (x + 300, 197), color, -1)
        put_text(canvas, label.upper(), (x + 18, 225), 0.38, MUTED, 1, 260)
        put_text(canvas, value, (x + 18, 282), 0.88, color, 2)
        x += 325

    blocks = [
        ("1. Classifier anchor", SUPP / "audit_uni_loto/fig_H_uni_loto.png", (70, 380, 740, 815), GREEN),
        ("2. External H&E overlay", OUT / "fig02_gse230424_cv2_he_overlay_mosaic.png", (805, 380, 1535, 815), GREEN),
        ("3. GSE250521 spatial mosaic", OUT / "fig03_gse250521_cv2_uni_slide_mosaic.png", (1600, 380, 2330, 815), TEAL),
        ("4. Hotspot contour proof", OUT / "fig05_gse230424_cv2_hotspot_contours.png", (70, 900, 740, 1335), GREEN),
        ("5. Heterogeneity matrix", OUT / "fig07_paper2_cv2_sample_matrix.png", (805, 900, 1535, 1335), GOLD),
        ("6. Claim ladder", OUT / "fig06_paper2_cv2_claim_ladder.png", (1600, 900, 2330, 1335), GOLD),
    ]
    for title, path, box, accent in blocks:
        x0, y0, x1, y1 = box
        rect(canvas, box, PANEL, radius=16)
        cv2.rectangle(canvas, (x0, y0), (x1, y0 + 8), accent, -1)
        put_text(canvas, title, (x0 + 20, y0 + 38), 0.50, INK, 1, x1 - x0 - 40)
        src = cv2.imread(str(path), cv2.IMREAD_COLOR)
        if src is None:
            src = np.full((360, 520, 3), WHITE, dtype=np.uint8)
            put_text(src, "missing", (30, 190), 1.0, RED, 2)
        thumb = fit_image(src, (x1 - x0 - 34, y1 - y0 - 70), bg=WHITE)
        canvas[y0 + 55 : y0 + 55 + thumb.shape[0], x0 + 17 : x0 + 17 + thumb.shape[1]] = thumb

    rect(canvas, (70, 1410, 2330, 1535), PANEL2, radius=18)
    put_text(canvas, "Disposition", (105, 1462), 0.74, GOLD, 2)
    put_text(
        canvas,
        "Use for Paper 2 visual review: pathology carries a measurable local DM1/RAI spatial signal. Keep QC/smoothness caveats visible; do not use this as Paper 1 MAPK mechanism rescue.",
        (300, 1462),
        0.60,
        INK,
        1,
        1920,
        12,
    )
    cv2.imwrite(str(OUT / "fig08_paper2_cv2_onepage_storyboard.png"), canvas)


def make_reviewer_objection_grid(metrics: dict) -> None:
    objections = [
        ("Only image confounding?", "UNI LOTO AUC 0.852", "Use with RAS caveat", GREEN),
        ("Only stage?", "Stage-centered rho 0.392", "Stage means do not erase G250", GREEN),
        ("Only large N?", "G230 downsample rho 0.641", "200 spots/sample holds", GREEN),
        ("Only threshold?", "Decile Spearman 0.607", "Dose-response, not only top-k", GREEN),
        ("Only broad coordinates?", "Block rho 0.497", "Broad-domain artifact unlikely", GOLD),
        ("Specificity proven?", "Smoothness residual -0.121", "No; disclose caveat", RED),
        ("Every thyroid ST works?", "GSE248205 no-go", "No; negative control", RED),
        ("Paper 1 rescue?", "Spatial MAPK closed", "No; Paper 2 only", RED),
    ]
    canvas = np.full((1320, 1900, 3), BG, dtype=np.uint8)
    put_text(canvas, "Reviewer objection grid", (60, 74), 0.95, INK, 2)
    put_text(
        canvas,
        "Each box pairs an anticipated attack with the data layer that answers it. Red boxes are boundaries to disclose, not weaknesses to hide.",
        (60, 122),
        0.55,
        MUTED,
        1,
        1500,
    )
    x0, y0 = 70, 205
    w, h = 420, 235
    for i, (question, answer, disposition, accent) in enumerate(objections):
        row, col = divmod(i, 4)
        x = x0 + col * 455
        y = y0 + row * 285
        rect(canvas, (x, y, x + w, y + h), PANEL, radius=16)
        cv2.rectangle(canvas, (x, y), (x + w, y + 8), accent, -1)
        put_text(canvas, question, (x + 22, y + 46), 0.56, INK, 2, w - 44)
        put_text(canvas, answer, (x + 22, y + 112), 0.76, accent, 2, w - 44)
        put_text(canvas, disposition, (x + 22, y + 172), 0.50, MUTED, 1, w - 44)

    rect(canvas, (70, 870, 1830, 1180), PANEL2, radius=18)
    put_text(canvas, "How to use this board", (105, 930), 0.72, GOLD, 2)
    lines = [
        "Green: main support to lead with in Paper 2 figures.",
        "Gold: support that needs one sentence of caution.",
        "Red: boundaries that prevent overclaiming and make the story more credible.",
        "Do not convert the red boxes into manuscript claims; keep them as reviewer-defense caveats.",
    ]
    y = 985
    for line in lines:
        y = put_text(canvas, line, (115, y), 0.62, INK, 1, 1600, 12) + 2
    cv2.imwrite(str(OUT / "fig09_paper2_cv2_reviewer_objection_grid.png"), canvas)


def write_caption_scaffold(metrics: dict) -> None:
    rows = [
        ("Fig CV2-01", "Evidence scorecard", "Single-board summary of classifier, spatial projection, external support, and caveat layers.", "PAPER2_CV2_VISUAL_SUMMARY.json"),
        ("Fig CV2-02", "GSE230424 H&E overlay", "Observed and H&E-predicted DM1/low-RAI spot overlays on raw GSE230424 H&E crops.", "gse230424_pathology_predictions.tsv.gz; *_HE.jpg.gz"),
        ("Fig CV2-03", "GSE250521 UNI mosaic", "Observed versus predicted DM1/RAI spatial maps across 16 GSE250521 slides.", "path2space_spot_predictions.tsv.gz"),
        ("Fig CV2-04", "Result gallery", "Index of the Paper 2 evidence and control panels.", "hub assets and analysis_supp outputs"),
        ("Fig CV2-05", "GSE230424 hotspot contours", "Top-decile observed and predicted hotspot contours on true H&E.", "paper2_cv2_gse230424_hotspot_contours.tsv"),
        ("Fig CV2-06", "Claim ladder", "Reviewer-facing claim ladder with caveats co-located.", "summary JSON files"),
        ("Fig CV2-07", "Sample matrix", "Slide/sample-level rho heatmap for GSE250521 and GSE230424.", "path2space_loso_slide_summary.tsv; gse230424 predictions"),
        ("Fig CV2-08", "One-page storyboard", "Presentation-ready one-page Paper 2 visual board.", "CV2 generated from local outputs"),
        ("Fig CV2-09", "Reviewer objection grid", "Attack-response grid for figure planning and reviewer defense.", "CV2 generated from local outputs"),
    ]
    df = pd.DataFrame(rows, columns=["figure_id", "working_title", "caption_scaffold", "source"])
    df.to_csv(OUT / "paper2_cv2_caption_scaffold.tsv", sep="\t", index=False)
    lines = [
        "# Paper 2 CV2 caption scaffold",
        "",
        "Scaffold only. Convert to final author voice later.",
        "",
        df.to_markdown(index=False),
        "",
        "Safe claim boundary: Paper 2 image-to-spatial-RNA support; QC/smoothness caveats remain explicit; not Paper 1 MAPK mechanism rescue.",
    ]
    (OUT / "PAPER2_CV2_CAPTION_SCAFFOLD.md").write_text("\n".join(lines) + "\n")


def draw_ladder_node(
    img: np.ndarray,
    x: int,
    y: int,
    w: int,
    h: int,
    idx: int,
    title: str,
    value: str,
    detail: str,
    accent: tuple[int, int, int],
) -> None:
    rect(img, (x, y, x + w, y + h), PANEL, radius=18)
    cv2.rectangle(img, (x, y), (x + w, y + 8), accent, -1)
    cv2.circle(img, (x + 38, y + 48), 22, accent, -1, cv2.LINE_AA)
    put_text(img, str(idx), (x + 27, y + 58), 0.74, BG, 2)
    put_text(img, title, (x + 78, y + 44), 0.58, INK, 2, w - 100)
    put_text(img, value, (x + 78, y + 100), 0.95, accent, 2, w - 100)
    put_text(img, detail, (x + 28, y + 152), 0.48, MUTED, 1, w - 56)


def make_claim_ladder(metrics: dict) -> None:
    img = np.full((1200, 1900, 3), BG, dtype=np.uint8)
    put_text(img, "Paper 2 claim ladder: what is strong, what is bounded", (60, 74), 0.95, INK, 2)
    put_text(
        img,
        "This figure is designed for reviewer triage. Green nodes are use-now evidence; gold nodes are support with caveat; red nodes are explicit boundaries.",
        (60, 122),
        0.55,
        MUTED,
        1,
        1500,
    )
    nodes = [
        (
            "Image-DM1 classifier",
            f"AUC {metrics['uni']['pooled_overall_auc']:.3f}",
            "UNI leave-one-TSS-out audit anchors the image side. RAS-like remains weaker and is not hidden.",
            GREEN,
        ),
        (
            "Spatial RNA projection",
            f"rho {metrics['g250']['dm1_smoothed_rho']:.3f}",
            "GSE250521 follows Path2Space logic: smooth local labels, then evaluate H&E-to-spatial-RNA.",
            TEAL,
        ),
        (
            "External thyroid H&E",
            f"rho {metrics['g230']['top_he_sample_centered_rho']:.3f}",
            "GSE230424 supports a local thyroid DM1/low-RAI axis from true H&E tiles.",
            GREEN,
        ),
        (
            "Robustness controls",
            f"downsample {metrics['down']['gse230424_raw_rho_median']:.3f}",
            "Decile dose-response, blocks, bootstrap, hotspots, and 200-spot sampling reduce easy objections.",
            GOLD,
        ),
        (
            "Specificity boundary",
            "smoothness fail",
            f"GSE230424 adjusted residual {metrics['smooth']['raw_smoothness_adjusted_residual']:.3f}; do not claim smoothness-independent molecular specificity.",
            RED,
        ),
        (
            "Final disposition",
            "Paper 2 use",
            "Visual/spatial biomarker story is stronger; Paper 1 MAPK mechanism rescue remains closed.",
            GOLD,
        ),
    ]
    x_positions = [65, 665, 1265, 1265, 665, 65]
    y_positions = [205, 205, 205, 625, 625, 625]
    w, h = 540, 300
    for i, (title, value, detail, accent) in enumerate(nodes):
        draw_ladder_node(img, x_positions[i], y_positions[i], w, h, i + 1, title, value, detail, accent)
    arrows = [
        ((610, 355), (650, 355)),
        ((1210, 355), (1250, 355)),
        ((1535, 510), (1535, 610)),
        ((1260, 775), (1215, 775)),
        ((660, 775), (615, 775)),
    ]
    for p0, p1 in arrows:
        cv2.arrowedLine(img, p0, p1, GOLD, 3, cv2.LINE_AA, tipLength=0.35)
    rect(img, (65, 1000, 1830, 1130), PANEL2, radius=16)
    put_text(
        img,
        "One-line frame: H&E carries a measurable spatial DM1/RAI signal across image classifier and Visium cohorts; QC/smoothness caveats are co-located, not buried.",
        (95, 1062),
        0.66,
        INK,
        1,
        1690,
    )
    cv2.imwrite(str(OUT / "fig06_paper2_cv2_claim_ladder.png"), img)


def draw_heat_cell(
    img: np.ndarray,
    box: tuple[int, int, int, int],
    value: float,
    label: str,
    vmin: float = -0.1,
    vmax: float = 0.7,
) -> None:
    x0, y0, x1, y1 = box
    if np.isfinite(value):
        t = float(np.clip((value - vmin) / max(vmax - vmin, 1e-9), 0, 1))
        color = (
            int((1 - t) * 70 + t * 70),
            int((1 - t) * 70 + t * 220),
            int((1 - t) * 95 + t * 190),
        )
    else:
        color = (70, 70, 80)
    cv2.rectangle(img, (x0, y0), (x1, y1), color, -1)
    cv2.rectangle(img, (x0, y0), (x1, y1), (24, 29, 38), 1, cv2.LINE_AA)
    put_text(img, label, (x0 + 8, y0 + 24), 0.36, INK, 1, x1 - x0 - 16)
    val = "NA" if not np.isfinite(value) else f"{value:.2f}"
    put_text(img, val, (x0 + 8, y1 - 12), 0.48, WHITE, 2)


def make_sample_matrix(metrics: dict) -> None:
    g250_slide = pd.read_csv(G250 / "path2space_loso_slide_summary.tsv", sep="\t")
    g250_slide = g250_slide[(g250_slide["target"] == "DM1_like_score") & (g250_slide["mode"] == "smooth8")].copy()
    pred230 = pd.read_csv(G230 / "gse230424_pathology_predictions.tsv.gz", sep="\t")
    g230_rows = []
    for sample, sub in pred230.groupby("sample", sort=True):
        rho, _, _ = safe_spearman(
            center_by_group(sub["obs_DM1_low_RAI_score_smooth8"].to_numpy(dtype=float), sub["sample"].to_numpy(str)),
            center_by_group(sub["pred_DM1_low_RAI_score_HE_tile_features"].to_numpy(dtype=float), sub["sample"].to_numpy(str)),
        )
        g230_rows.append({"sample": sample, "stage": {"P1": "PTC+HT", "P2": "PTC+HT", "P3": "HT", "P4": "HT"}.get(sample, ""), "rho": rho})

    img = np.full((1320, 1900, 3), BG, dtype=np.uint8)
    put_text(img, "Sample-level evidence matrix", (60, 74), 0.95, INK, 2)
    put_text(
        img,
        "Each cell is a held-out slide/sample rho. This is the fast scan for heterogeneity: PTC/LPTC strong, ATC weaker, GSE230424 externally positive.",
        (60, 122),
        0.55,
        MUTED,
        1,
        1500,
    )
    put_text(img, "GSE250521: 16 thyroid trajectory Visium slides", (70, 200), 0.62, GOLD, 2)
    stage_order = {"PT": 0, "PTC": 1, "LPTC": 2, "ATC": 3}
    g250_slide = g250_slide.sort_values(["stage", "sample_id"], key=lambda s: s.map(stage_order).fillna(9) if s.name == "stage" else s)
    cell_w, cell_h = 205, 105
    x0, y0 = 70, 235
    for i, row in enumerate(g250_slide.itertuples(index=False)):
        r, c = divmod(i, 8)
        label = f"{row.stage} {str(row.sample_id).split('_')[-1]}"
        draw_heat_cell(img, (x0 + c * cell_w, y0 + r * cell_h, x0 + (c + 1) * cell_w - 8, y0 + (r + 1) * cell_h - 8), float(row.rho), label)
    put_text(img, "GSE230424: external thyroid Visium H&E samples", (70, 540), 0.62, GOLD, 2)
    for i, row in enumerate(g230_rows):
        draw_heat_cell(img, (70 + i * 260, 575, 70 + (i + 1) * 260 - 8, 700), float(row["rho"]), f"{row['sample']} {row['stage']}")

    rect(img, (70, 775, 1830, 1165), PANEL, radius=18)
    put_text(img, "What the matrix says", (105, 835), 0.72, GOLD, 2)
    notes = [
        f"GSE250521 overall smoothed rho {metrics['g250']['dm1_smoothed_rho']:.3f}; slide-centered rho {metrics['g250_stage']['dm1_existing_slide_centered_rho']:.3f}.",
        "PTC/LPTC are visually and statistically strongest; ATC is the weaker stage and should stay disclosed.",
        f"GSE230424 external sample-centered rho {metrics['g230']['top_he_sample_centered_rho']:.3f}; 200-spot downsample median {metrics['down']['gse230424_raw_rho_median']:.3f}.",
        "This is heterogeneity-aware support, not a claim that every sample is equally strong.",
    ]
    y = 895
    for note in notes:
        y = put_text(img, note, (115, y), 0.62, INK, 1, 1640, 12) + 4
    cv2.imwrite(str(OUT / "fig07_paper2_cv2_sample_matrix.png"), img)


def write_report(metrics: dict) -> None:
    summary = {
        "created": "2026-05-10",
        "output_dir": str(OUT),
        "figures": [
            "fig01_paper2_cv2_evidence_scorecard.png",
            "fig02_gse230424_cv2_he_overlay_mosaic.png",
        "fig03_gse250521_cv2_uni_slide_mosaic.png",
        "fig04_paper2_cv2_result_gallery.png",
        "fig05_gse230424_cv2_hotspot_contours.png",
        "fig06_paper2_cv2_claim_ladder.png",
        "fig07_paper2_cv2_sample_matrix.png",
        "fig08_paper2_cv2_onepage_storyboard.png",
        "fig09_paper2_cv2_reviewer_objection_grid.png",
        ],
        "headline": {
            "uni_loto_auc": metrics["uni"]["pooled_overall_auc"],
            "gse250521_smoothed_rho": metrics["g250"]["dm1_smoothed_rho"],
            "gse250521_slide_centered_rho": metrics["g250_stage"]["dm1_existing_slide_centered_rho"],
            "gse230424_sample_centered_rho": metrics["g230"]["top_he_sample_centered_rho"],
            "gse230424_downsample_rho": metrics["down"]["gse230424_raw_rho_median"],
            "gse230424_residual_downsample_rho": metrics["down"]["gse230424_resid_rho_median"],
            "gse230424_smoothness_adjusted_residual": metrics["smooth"]["raw_smoothness_adjusted_residual"],
        },
        "safe_disposition": "Paper 2 visual/spatial support with QC and smoothness caveats; not Paper 1 MAPK mechanism evidence.",
    }
    (OUT / "PAPER2_CV2_VISUAL_SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n")
    lines = [
        "# Paper 2 CV2 visual summary",
        "",
        "## Figures",
        "",
        "- `fig01_paper2_cv2_evidence_scorecard.png` — single-board claim map.",
        "- `fig02_gse230424_cv2_he_overlay_mosaic.png` — actual H&E crops with observed/predicted DM1/low-RAI spot overlays.",
        "- `fig03_gse250521_cv2_uni_slide_mosaic.png` — 16-slide observed vs predicted UNI spatial mosaic.",
        "- `fig04_paper2_cv2_result_gallery.png` — figure-map index for the current Paper 2 evidence stack.",
        "- `fig05_gse230424_cv2_hotspot_contours.png` — top-decile observed/predicted hotspot contours on true H&E.",
        "- `fig06_paper2_cv2_claim_ladder.png` — visual claim ladder with caveats co-located.",
        "- `fig07_paper2_cv2_sample_matrix.png` — slide/sample-level heterogeneity matrix.",
        "",
        "## Headline",
        "",
        f"- UNI LOTO AUC: {metrics['uni']['pooled_overall_auc']:.3f}.",
        f"- GSE250521 smoothed DM1/RAI rho: {metrics['g250']['dm1_smoothed_rho']:.3f}; slide-centered {metrics['g250_stage']['dm1_existing_slide_centered_rho']:.3f}.",
        f"- GSE230424 H&E sample-centered rho: {metrics['g230']['top_he_sample_centered_rho']:.3f}; downsample median {metrics['down']['gse230424_raw_rho_median']:.3f}.",
        f"- GSE230424 residual downsample median rho: {metrics['down']['gse230424_resid_rho_median']:.3f}.",
        "",
        "Safe framing: Paper 2 image-to-spatial-RNA support, with QC/smoothness caveats kept visible.",
    ]
    (OUT / "PAPER2_CV2_VISUAL_SUMMARY_REPORT.md").write_text("\n".join(lines) + "\n")
    (OUT / "SUMMARY.md").write_text("\n".join(lines[:18]) + "\n")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    metrics = load_metrics()
    make_evidence_scorecard(metrics)
    make_gse230424_he_overlay(metrics)
    make_gse230424_hotspot_contours(metrics)
    make_gse250521_uni_mosaic(metrics)
    make_claim_ladder(metrics)
    make_sample_matrix(metrics)
    make_result_gallery(metrics)
    make_onepage_storyboard(metrics)
    make_reviewer_objection_grid(metrics)
    write_caption_scaffold(metrics)
    write_report(metrics)
    print(json.dumps({"out": str(OUT), "status": "ok"}, indent=2))


if __name__ == "__main__":
    main()
