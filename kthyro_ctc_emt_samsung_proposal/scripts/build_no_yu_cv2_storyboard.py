#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import cv2
import numpy as np


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT / "kthyro_ctc_emt_samsung_proposal/outputs/in_silico_paper_no_yu"
FIG = OUT / "figures"
VIS = OUT / "visual_qc"

BG = (8, 13, 24)
PANEL = (16, 26, 44)
LINE = (48, 68, 99)
TEXT = (235, 242, 255)
MUTED = (170, 184, 202)
GOLD = (138, 210, 255)
CYAN = (255, 220, 91)


def read_img(path: Path) -> np.ndarray:
    img = cv2.imread(str(path), cv2.IMREAD_COLOR)
    if img is None:
        raise FileNotFoundError(path)
    return img


def fit_to_canvas(img: np.ndarray, w: int, h: int) -> np.ndarray:
    ih, iw = img.shape[:2]
    scale = min(w / iw, h / ih)
    nw, nh = int(iw * scale), int(ih * scale)
    resized = cv2.resize(img, (nw, nh), interpolation=cv2.INTER_AREA)
    canvas = np.full((h, w, 3), BG, dtype=np.uint8)
    x, y = (w - nw) // 2, (h - nh) // 2
    canvas[y : y + nh, x : x + nw] = resized
    return canvas


def card(img: np.ndarray, title: str, idx: int, w: int = 760, h: int = 520) -> np.ndarray:
    outer = np.full((h, w, 3), BG, dtype=np.uint8)
    cv2.rectangle(outer, (0, 0), (w - 1, h - 1), LINE, 2)
    cv2.rectangle(outer, (0, 0), (w - 1, 54), PANEL, -1)
    cv2.putText(outer, f"{idx:02d}", (22, 36), cv2.FONT_HERSHEY_SIMPLEX, 0.75, GOLD, 2, cv2.LINE_AA)
    cv2.putText(outer, title[:52], (78, 36), cv2.FONT_HERSHEY_SIMPLEX, 0.68, TEXT, 2, cv2.LINE_AA)
    fitted = fit_to_canvas(img, w - 28, h - 78)
    outer[66 : 66 + fitted.shape[0], 14 : 14 + fitted.shape[1]] = fitted
    return outer


def hstack_with_gap(items: list[np.ndarray], gap: int = 18) -> np.ndarray:
    h = max(x.shape[0] for x in items)
    total_w = sum(x.shape[1] for x in items) + gap * (len(items) - 1)
    out = np.full((h, total_w, 3), BG, dtype=np.uint8)
    x = 0
    for im in items:
        y = (h - im.shape[0]) // 2
        out[y : y + im.shape[0], x : x + im.shape[1]] = im
        x += im.shape[1] + gap
    return out


def vstack_with_gap(items: list[np.ndarray], gap: int = 18) -> np.ndarray:
    w = max(x.shape[1] for x in items)
    total_h = sum(x.shape[0] for x in items) + gap * (len(items) - 1)
    out = np.full((total_h, w, 3), BG, dtype=np.uint8)
    y = 0
    for im in items:
        x = (w - im.shape[1]) // 2
        out[y : y + im.shape[0], x : x + im.shape[1]] = im
        y += im.shape[0] + gap
    return out


def add_header(body: np.ndarray, title: str, subtitle: str) -> np.ndarray:
    h, w = body.shape[:2]
    header_h = 150
    out = np.full((h + header_h, w, 3), BG, dtype=np.uint8)
    out[header_h:] = body
    cv2.putText(out, title, (34, 62), cv2.FONT_HERSHEY_SIMPLEX, 1.35, TEXT, 3, cv2.LINE_AA)
    cv2.putText(out, subtitle, (36, 108), cv2.FONT_HERSHEY_SIMPLEX, 0.65, MUTED, 2, cv2.LINE_AA)
    cv2.line(out, (34, 132), (w - 34, 132), LINE, 2)
    return out


def main() -> None:
    VIS.mkdir(parents=True, exist_ok=True)
    specs = [
        ("F01_no_yu_claim_ladder.png", "Claim ladder"),
        ("F02_public_data_layer_map.png", "Public data layers"),
        ("F03_tcga_braf_state_split.png", "BRAF state split"),
        ("F04_tcga_driver_variance_boundary.png", "Driver variance"),
        ("F05_spatial_tissue_state_organization.png", "Spatial organization"),
        ("F06_scrna_marker_context.png", "scRNA marker context"),
        ("F07_proteomics_dediff_direction.png", "Proteomics direction"),
        ("F08_ctc_emt_ngs_prior_panel.png", "CTC-EMT-NGS panel"),
        ("F09_publishability_decision.png", "Publishability"),
    ]
    cards = []
    for i, (fn, title) in enumerate(specs, 1):
        cards.append(card(read_img(FIG / fn), title, i))
    rows = [hstack_with_gap(cards[i : i + 3]) for i in range(0, len(cards), 3)]
    montage = add_header(
        vstack_with_gap(rows),
        "No-Yu-data In Silico Paper - Visual Storyboard",
        "CV2 contact sheet: figure-first review before writing claims",
    )
    cv2.imwrite(str(VIS / "no_yu_in_silico_visual_storyboard_cv2.png"), montage)

    priority = [cards[0], cards[2], cards[4], cards[7], cards[8]]
    row1 = hstack_with_gap(priority[:2])
    row2 = hstack_with_gap(priority[2:])
    priority_board = add_header(
        vstack_with_gap([row1, row2]),
        "Decision Board: What Can Be Published Now",
        "Keep this board in front of the proposal/manuscript team",
    )
    cv2.imwrite(str(VIS / "no_yu_decision_board_cv2.png"), priority_board)

    print(VIS / "no_yu_in_silico_visual_storyboard_cv2.png")
    print(VIS / "no_yu_decision_board_cv2.png")


if __name__ == "__main__":
    main()
