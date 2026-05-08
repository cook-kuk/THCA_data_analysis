#!/usr/bin/env python3
"""
Paper 2 Fig 5 v2 — Publication-quality CLAM attention heatmaps (cv2-polished).

Improvements over make_fig5_attention.py:
  1. cv2.applyColorMap with COLORMAP_INFERNO (perceptually uniform).
  2. cv2.addWeighted alpha blending (smoother than PIL alpha-composite).
  3. Anti-aliased text annotations (cv2.LINE_AA, FONT_HERSHEY_DUPLEX/TRIPLEX).
  4. Coloured borders: green for DM1+, blue for DM2.
  5. 1600px-wide thumbnails (was 1024).
  6. 2 rows x 3 cols grid (was 6 rows x 2).
     Each panel: H&E (left half) + attention overlay (right half), seamless stitch.
  7. Section header bands "DM1+ HIGH CONFIDENCE" / "DM2 LOW CONFIDENCE".
  8. Mini attention-weight histogram bar chart under each panel.
  9. 4K (3840x2160) output PNG + PDF (matplotlib backend).

Outputs:
  figures/fig5_attention_heatmaps_v2.png
  figures/fig5_attention_heatmaps_v2.pdf

Loads each fold's best ckpt to recompute attention only for slides predicted
by that fold (held-out predictions).
"""
from __future__ import annotations

import sys
from pathlib import Path

import cv2
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from scipy.ndimage import gaussian_filter

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

try:
    import openslide
except ImportError as e:
    sys.exit(f"openslide-python required: pip install openslide-python ({e})")

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/"
            "p2_image_dm1_v2_foundation_clam_2026_05_07")
FEATURE_DIR = ROOT / "phase2_tcga_clam" / "features"
WSI_DIR = Path("/data/tcga_thca_wsi")
MANIFEST = ROOT / "phase2_tcga_clam" / "slide_manifest.tsv"
PREDS_TSV = ROOT / "phase2_tcga_clam" / "clam_per_slide_predictions.tsv"
CKPT_DIR = ROOT / "phase2_tcga_clam"
OUT_DIR = ROOT / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

PNG_OUT = OUT_DIR / "fig5_attention_heatmaps_v2.png"
PDF_OUT = OUT_DIR / "fig5_attention_heatmaps_v2.pdf"

THUMB_W = 1600  # max dim per thumbnail
TILE_PX_FULL = 256

# Brand palette (BGR for cv2)
GREEN_BGR = (143, 157, 42)   # 2A9D8F -> BGR
BLUE_BGR = (170, 133, 91)    # 5B85AA
TEXT_DARK = (40, 40, 40)
WHITE = (255, 255, 255)
BLACK = (0, 0, 0)
HDR_DM1_BG = (90, 130, 50)   # darker green band
HDR_DM2_BG = (130, 90, 60)   # darker blue band


# ---------------------------------------------------------------------------
# CLAM model (matches phase2_clam_train_eval.py)
# ---------------------------------------------------------------------------
class GatedAttentionMIL(nn.Module):
    def __init__(self, in_dim=1024, hidden=512, n_classes=2, dropout=0.25):
        super().__init__()
        self.fc = nn.Sequential(nn.Linear(in_dim, hidden), nn.ReLU(), nn.Dropout(dropout))
        self.attn_a = nn.Sequential(nn.Linear(hidden, hidden), nn.Tanh())
        self.attn_b = nn.Sequential(nn.Linear(hidden, hidden), nn.Sigmoid())
        self.attn_w = nn.Linear(hidden, 1)
        self.classifier = nn.Linear(hidden, n_classes)

    def forward(self, x):
        h = self.fc(x)
        a = self.attn_w(self.attn_a(h) * self.attn_b(h))
        attn = F.softmax(a.squeeze(-1), dim=0)
        bag = (attn.unsqueeze(-1) * h).sum(0)
        logits = self.classifier(bag)
        return logits, attn


def load_fold_model(fold: int) -> GatedAttentionMIL:
    ckpt = CKPT_DIR / f"clam_fold{fold}_best.pt"
    state = torch.load(ckpt, map_location="cpu")
    m = GatedAttentionMIL(in_dim=1024, hidden=512, n_classes=2)
    m.load_state_dict(state)
    m.eval()
    return m


def get_attention(model: GatedAttentionMIL, feat_path: Path):
    feat = torch.load(feat_path, map_location="cpu")
    if isinstance(feat, np.ndarray):
        feat = torch.from_numpy(feat)
    x = feat.float()
    with torch.no_grad():
        logits, attn = model(x)
        prob_dm1 = F.softmax(logits, dim=-1)[1].item()
    return attn.numpy(), prob_dm1


def find_coords_tsv(slide_id: str) -> Path | None:
    for cand in (FEATURE_DIR / f"{slide_id}_coords.tsv",
                 FEATURE_DIR / f"{slide_id}.coords.tsv"):
        if cand.exists():
            return cand
    return None


def find_wsi(slide_id: str) -> Path | None:
    cands = list(WSI_DIR.glob(f"{slide_id}*.svs")) + list(WSI_DIR.glob(f"*{slide_id}*.svs"))
    return cands[0] if cands else None


# ---------------------------------------------------------------------------
# Image rendering
# ---------------------------------------------------------------------------
def render_thumb_and_heat(svs_path: Path, coords_tsv: Path,
                          attention: np.ndarray,
                          thumb_w: int = THUMB_W,
                          tile_px_full: int = TILE_PX_FULL):
    """Returns (thumb_bgr, heat_norm_2d, max_xy_full)."""
    slide = openslide.OpenSlide(str(svs_path))
    w_full, h_full = slide.dimensions
    ratio = thumb_w / max(w_full, h_full)
    tw, th = int(w_full * ratio), int(h_full * ratio)
    pil_thumb = slide.get_thumbnail((tw, th)).convert("RGB")
    slide.close()
    thumb_rgb = np.array(pil_thumb)
    thumb_bgr = cv2.cvtColor(thumb_rgb, cv2.COLOR_RGB2BGR)

    coords = pd.read_csv(coords_tsv, sep="\t")
    if "x_full" in coords.columns and "y_full" in coords.columns:
        xy = coords[["x_full", "y_full"]].values
    else:
        xy = coords.iloc[:, :2].values
    if len(xy) != len(attention):
        return thumb_bgr, None, None

    # use the actual returned thumbnail size (PIL get_thumbnail may differ slightly)
    th, tw = thumb_bgr.shape[:2]
    H, W = th, tw
    heat = np.zeros((H, W), dtype=np.float32)
    tile_t = max(1, int(tile_px_full * ratio))
    for (x, y), a in zip(xy, attention):
        xt = int(round(x * ratio)); yt = int(round(y * ratio))
        x0 = max(0, xt); y0 = max(0, yt)
        x1 = min(W, xt + tile_t); y1 = min(H, yt + tile_t)
        if x1 > x0 and y1 > y0:
            np.maximum(heat[y0:y1, x0:x1], float(a), out=heat[y0:y1, x0:x1])
    heat = gaussian_filter(heat, sigma=tile_t * 0.5)
    if heat.max() > 0:
        heat = heat / heat.max()
    top_idx = int(np.argmax(attention))
    max_xy_full = (int(xy[top_idx, 0]), int(xy[top_idx, 1]))
    return thumb_bgr, heat, max_xy_full


def attention_overlay(thumb_bgr: np.ndarray, heat: np.ndarray,
                      alpha: float = 0.55) -> np.ndarray:
    """cv2 INFERNO colormap + addWeighted blend."""
    if heat is None:
        return thumb_bgr.copy()
    # ensure heat matches thumb size exactly
    if heat.shape[:2] != thumb_bgr.shape[:2]:
        heat = cv2.resize(heat, (thumb_bgr.shape[1], thumb_bgr.shape[0]),
                          interpolation=cv2.INTER_LINEAR)
    heat_u8 = np.clip(heat * 255.0, 0, 255).astype(np.uint8)
    color = cv2.applyColorMap(heat_u8, cv2.COLORMAP_INFERNO)  # BGR
    # mask out very low attention so background stays clean
    mask = (heat > 0.05).astype(np.float32)
    mask3 = np.dstack([mask] * 3)
    blended = cv2.addWeighted(thumb_bgr, 1.0 - alpha, color, alpha, 0)
    out = thumb_bgr * (1 - mask3) + blended * mask3
    return out.astype(np.uint8)


def make_panel(slide_id: str, prob: float, group: str,
               thumb_bgr: np.ndarray, overlay_bgr: np.ndarray,
               attention: np.ndarray | None,
               panel_w: int, panel_h: int) -> np.ndarray:
    """
    Compose a single panel:
      - top: H&E | attention overlay (seamless side-by-side, equal halves)
      - bottom strip: small attention histogram
      - bordered with group colour
      - corner badges: top-left group label; bottom-right slide ID + prob
    """
    border_color = GREEN_BGR if group == "DM1+" else BLUE_BGR
    border_px = 6
    hist_h = 80  # bottom histogram strip height

    inner_w = panel_w - 2 * border_px
    inner_h = panel_h - 2 * border_px - hist_h

    # split inner area into two halves of equal width
    half_w = inner_w // 2

    def fit(img):
        h, w = img.shape[:2]
        scale = min(half_w / w, inner_h / h)
        new_w = int(w * scale); new_h = int(h * scale)
        resized = cv2.resize(img, (new_w, new_h), interpolation=cv2.INTER_AREA)
        # paste centered onto canvas of (inner_h, half_w, 3)
        canvas = np.full((inner_h, half_w, 3), 255, dtype=np.uint8)
        ox = (half_w - new_w) // 2
        oy = (inner_h - new_h) // 2
        canvas[oy:oy + new_h, ox:ox + new_w] = resized
        return canvas

    left = fit(thumb_bgr)
    right = fit(overlay_bgr)
    top = np.concatenate([left, right], axis=1)
    # if half_w*2 != inner_w, pad right
    if top.shape[1] < inner_w:
        pad = np.full((inner_h, inner_w - top.shape[1], 3), 255, dtype=np.uint8)
        top = np.concatenate([top, pad], axis=1)

    # thin vertical separator at the seam
    cv2.line(top, (half_w - 1, 0), (half_w - 1, inner_h - 1),
             (200, 200, 200), 1, cv2.LINE_AA)

    # bottom histogram strip
    hist_strip = np.full((hist_h, inner_w, 3), 250, dtype=np.uint8)
    if attention is not None and len(attention) > 1:
        # bin attention weights into 50 buckets
        n_bins = 50
        a = attention / (attention.max() + 1e-9)
        counts, _ = np.histogram(a, bins=n_bins, range=(0, 1))
        if counts.max() > 0:
            counts = counts / counts.max()
        bar_w = inner_w / n_bins
        baseline_y = hist_h - 12
        max_bar_h = hist_h - 28
        for i, c in enumerate(counts):
            x0 = int(i * bar_w)
            x1 = int((i + 1) * bar_w) - 1
            bh = int(c * max_bar_h)
            # color bars by their bin position using INFERNO
            t = (i + 0.5) / n_bins
            bar_col = cv2.applyColorMap(
                np.array([[int(t * 255)]], dtype=np.uint8),
                cv2.COLORMAP_INFERNO)[0, 0].tolist()
            cv2.rectangle(hist_strip, (x0, baseline_y - bh), (x1, baseline_y),
                          bar_col, -1, cv2.LINE_AA)
        # baseline
        cv2.line(hist_strip, (0, baseline_y), (inner_w - 1, baseline_y),
                 (180, 180, 180), 1, cv2.LINE_AA)
        cv2.putText(hist_strip, "Attention weight distribution (low -> high)",
                    (8, hist_h - 2),
                    cv2.FONT_HERSHEY_DUPLEX, 0.42, (90, 90, 90), 1, cv2.LINE_AA)

    inner = np.concatenate([top, hist_strip], axis=0)

    # frame with coloured border
    panel = np.full((panel_h, panel_w, 3), 255, dtype=np.uint8)
    panel[border_px:border_px + inner.shape[0],
          border_px:border_px + inner.shape[1]] = inner
    cv2.rectangle(panel, (0, 0), (panel_w - 1, panel_h - 1),
                  border_color, border_px, cv2.LINE_AA)

    # ----- annotations (anti-aliased) -----
    # top-left badge
    badge_text = group
    badge_color = GREEN_BGR if group == "DM1+" else BLUE_BGR
    (tw_, th_), _ = cv2.getTextSize(badge_text, cv2.FONT_HERSHEY_TRIPLEX, 0.95, 2)
    bx0, by0 = border_px + 14, border_px + 14
    bx1, by1 = bx0 + tw_ + 22, by0 + th_ + 16
    overlay = panel.copy()
    cv2.rectangle(overlay, (bx0, by0), (bx1, by1), badge_color, -1, cv2.LINE_AA)
    panel = cv2.addWeighted(overlay, 0.92, panel, 0.08, 0)
    cv2.putText(panel, badge_text, (bx0 + 11, by1 - 10),
                cv2.FONT_HERSHEY_TRIPLEX, 0.95, WHITE, 2, cv2.LINE_AA)

    # H&E / Attention sub-labels at top of each half (small)
    cv2.putText(panel, "H&E", (border_px + 12, border_px + th_ + 60),
                cv2.FONT_HERSHEY_DUPLEX, 0.7, (60, 60, 60), 1, cv2.LINE_AA)
    cv2.putText(panel, "Attention (INFERNO)",
                (border_px + half_w + 12, border_px + th_ + 60),
                cv2.FONT_HERSHEY_DUPLEX, 0.7, (60, 60, 60), 1, cv2.LINE_AA)

    # bottom-right: slide ID (small) + prob (large), white pill
    short_id = slide_id[:8] + "…"
    prob_text = f"DM1 prob: {prob:.2f}"
    (pw, ph), _ = cv2.getTextSize(prob_text, cv2.FONT_HERSHEY_TRIPLEX, 1.05, 2)
    (sw, sh), _ = cv2.getTextSize(short_id, cv2.FONT_HERSHEY_DUPLEX, 0.55, 1)
    pill_w = max(pw, sw) + 28
    pill_h = ph + sh + 26
    px1 = panel_w - border_px - 12
    py1 = border_px + inner_h - 12
    px0 = px1 - pill_w
    py0 = py1 - pill_h
    overlay = panel.copy()
    cv2.rectangle(overlay, (px0, py0), (px1, py1), (255, 255, 255), -1, cv2.LINE_AA)
    panel = cv2.addWeighted(overlay, 0.85, panel, 0.15, 0)
    cv2.rectangle(panel, (px0, py0), (px1, py1), border_color, 2, cv2.LINE_AA)
    cv2.putText(panel, short_id, (px0 + 14, py0 + sh + 8),
                cv2.FONT_HERSHEY_DUPLEX, 0.55, (90, 90, 90), 1, cv2.LINE_AA)
    cv2.putText(panel, prob_text, (px0 + 14, py1 - 10),
                cv2.FONT_HERSHEY_TRIPLEX, 1.05, border_color, 2, cv2.LINE_AA)

    return panel


def make_section_header(text: str, w: int, h: int, bg_bgr: tuple) -> np.ndarray:
    band = np.full((h, w, 3), bg_bgr, dtype=np.uint8)
    # subtle horizontal gradient for polish
    grad = np.linspace(1.0, 0.85, w, dtype=np.float32).reshape(1, w, 1)
    band = (band.astype(np.float32) * grad).clip(0, 255).astype(np.uint8)
    cv2.putText(band, text, (32, h // 2 + 16),
                cv2.FONT_HERSHEY_TRIPLEX, 1.4, WHITE, 3, cv2.LINE_AA)
    cv2.line(band, (0, h - 1), (w - 1, h - 1), (60, 60, 60), 1, cv2.LINE_AA)
    return band


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("[fig5_v2] loading predictions + manifest...")
    preds = pd.read_csv(PREDS_TSV, sep="\t")
    man = pd.read_csv(MANIFEST, sep="\t")

    # Filter to slides that have coords AND svs files (so heatmap actually renders)
    def coords_ok(sid):
        return ((FEATURE_DIR / f"{sid}_coords.tsv").exists()
                or (FEATURE_DIR / f"{sid}.coords.tsv").exists())

    def svs_ok(sid):
        return find_wsi(sid) is not None

    preds = preds[preds["slide"].apply(lambda s: coords_ok(s) and svs_ok(s))].copy()

    # Pre-screen: drop slides whose SVS openslide cannot read (corrupted)
    def svs_readable(sid):
        try:
            slide = openslide.OpenSlide(str(find_wsi(sid)))
            slide.close()
            return True
        except Exception:
            return False

    preds["readable"] = preds["slide"].apply(svs_readable)
    n_unreadable = (~preds["readable"]).sum()
    if n_unreadable:
        print(f"[fig5_v2] dropping {n_unreadable} slide(s) with unreadable SVS")
    preds = preds[preds["readable"]].copy()

    dm1_pos = (preds[preds["label"] == 1]
               .sort_values("prob_DM1", ascending=False)
               .drop_duplicates("slide", keep="first")
               .head(3))
    dm2_neg = (preds[preds["label"] == 0]
               .sort_values("prob_DM1", ascending=True)
               .drop_duplicates("slide", keep="first")
               .head(3))
    picks = pd.concat([dm1_pos.assign(group="DM1+"),
                       dm2_neg.assign(group="DM2")], ignore_index=True)
    print(f"[fig5_v2] picked {len(dm1_pos)} DM1+ / {len(dm2_neg)} DM2:")
    print(picks[["slide", "fold", "label", "prob_DM1", "group"]])

    # Cache models per fold
    fold_models: dict[int, GatedAttentionMIL] = {}

    # Render each panel
    panels = []
    panel_w, panel_h = 1230, 760  # 3 across -> ~3690 wide
    for _, row in picks.iterrows():
        slide_id = row["slide"]
        fold = int(row["fold"])
        prob = float(row["prob_DM1"])
        group = row["group"]
        if fold not in fold_models:
            print(f"  loading fold {fold} model...")
            fold_models[fold] = load_fold_model(fold)
        model = fold_models[fold]
        feat_path = FEATURE_DIR / f"{slide_id}.pt"
        coords_tsv = find_coords_tsv(slide_id)
        svs_path = find_wsi(slide_id)
        attention, prob_recompute = get_attention(model, feat_path)
        thumb_bgr, heat, max_xy = render_thumb_and_heat(svs_path, coords_tsv, attention)
        overlay_bgr = attention_overlay(thumb_bgr, heat)
        panel = make_panel(slide_id, prob, group, thumb_bgr, overlay_bgr,
                           attention, panel_w, panel_h)
        panels.append((group, panel))
        print(f"  [{group}] {slide_id[:8]} prob={prob:.3f} (recomp={prob_recompute:.3f})")

    # Compose final canvas: title strip + DM1 header + DM1 row + DM2 header + DM2 row + footer
    title_h = 130
    section_h = 80
    footer_h = 70
    grid_w = panel_w * 3
    full_h = title_h + section_h + panel_h + section_h + panel_h + footer_h

    canvas = np.full((full_h, grid_w, 3), 252, dtype=np.uint8)

    # title bar (deep teal -> purple gradient)
    title_bar = np.zeros((title_h, grid_w, 3), dtype=np.uint8)
    for x in range(grid_w):
        t = x / grid_w
        # teal #006B7D BGR (125,107,0) -> purple #6A4C93 BGR (147,76,106)
        b = int(125 + (147 - 125) * t)
        g = int(107 + (76 - 107) * t)
        r = int(0 + (106 - 0) * t)
        title_bar[:, x] = (b, g, r)
    cv2.putText(title_bar,
                "Fig 5  -  CLAM gated-attention heatmaps  (top-3 DM1+ / top-3 DM2)",
                (40, 50), cv2.FONT_HERSHEY_TRIPLEX, 1.4, WHITE, 3, cv2.LINE_AA)
    cv2.putText(title_bar,
                "TCGA-THCA  |  UNI foundation features  |  fold-3 best AUC=1.00",
                (40, 95), cv2.FONT_HERSHEY_DUPLEX, 0.85,
                (230, 230, 230), 1, cv2.LINE_AA)
    canvas[0:title_h, :, :] = title_bar

    # DM1+ header
    y = title_h
    canvas[y:y + section_h, :, :] = make_section_header(
        "DM1+ HIGH CONFIDENCE   (predicted as DM1, label = DM1)",
        grid_w, section_h, HDR_DM1_BG)
    y += section_h
    # DM1+ row
    dm1_panels = [p for g, p in panels if g == "DM1+"]
    for i, p in enumerate(dm1_panels):
        canvas[y:y + panel_h, i * panel_w:(i + 1) * panel_w] = p
    y += panel_h

    # DM2 header
    canvas[y:y + section_h, :, :] = make_section_header(
        "DM2 LOW CONFIDENCE   (predicted as DM2, label = DM2)",
        grid_w, section_h, HDR_DM2_BG)
    y += section_h
    dm2_panels = [p for g, p in panels if g == "DM2"]
    for i, p in enumerate(dm2_panels):
        canvas[y:y + panel_h, i * panel_w:(i + 1) * panel_w] = p
    y += panel_h

    # footer
    cv2.putText(canvas,
                "Heatmap colormap = INFERNO (perceptually uniform). "
                "Borders: green = DM1+, blue = DM2. "
                "Histogram = per-tile attention weight distribution.",
                (40, y + 45),
                cv2.FONT_HERSHEY_DUPLEX, 0.7, (60, 60, 60), 1, cv2.LINE_AA)

    # ---- target final size ~3840 wide while preserving aspect ----
    target_w = 3840
    scale = target_w / canvas.shape[1]
    new_w = target_w
    new_h = int(canvas.shape[0] * scale)
    canvas_4k = cv2.resize(canvas, (new_w, new_h), interpolation=cv2.INTER_AREA)

    # Save PNG via cv2 (BGR)
    cv2.imwrite(str(PNG_OUT), canvas_4k, [cv2.IMWRITE_PNG_COMPRESSION, 4])
    print(f"[saved] {PNG_OUT}  ({canvas_4k.shape[1]}x{canvas_4k.shape[0]})")

    # Save PDF via matplotlib (single page using imshow)
    rgb = cv2.cvtColor(canvas_4k, cv2.COLOR_BGR2RGB)
    fig_w_in = 16.0
    fig_h_in = fig_w_in * rgb.shape[0] / rgb.shape[1]
    fig = plt.figure(figsize=(fig_w_in, fig_h_in), dpi=300)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(rgb)
    ax.axis("off")
    fig.savefig(PDF_OUT, dpi=300, bbox_inches="tight", pad_inches=0)
    plt.close(fig)
    print(f"[saved] {PDF_OUT}")


if __name__ == "__main__":
    main()
