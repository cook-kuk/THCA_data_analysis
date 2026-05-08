#!/usr/bin/env python3
"""
Paper 2 Fig 5 v3 — STRONG attention heatmap (강조점 또렷하게).

Fixes vs v2:
- Per-pixel alpha blending — high-attention pixels become near-fully colored,
  low-attention pixels stay raw H&E. No hard mask cutoff that hides everything.
- 95th-percentile normalization instead of max → more contrast.
- Smaller gaussian sigma → sharp hotspots, not smeared.
- Top-3 hotspot circles + arrows annotated on overlay.
- Side-by-side panel includes a colorbar.
- Per-slide DM1 prob big readable badge.
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
from scipy.ndimage import gaussian_filter, maximum_filter

import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams["font.family"] = ["NanumGothic", "DejaVu Sans"]

try:
    import openslide
except ImportError as e:
    sys.exit(f"openslide-python required: {e}")

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/"
            "p2_image_dm1_v2_foundation_clam_2026_05_07")
FEATURE_DIR = ROOT / "phase2_tcga_clam" / "features"
WSI_DIR = Path("/data/tcga_thca_wsi")
PREDS_TSV = ROOT / "phase2_tcga_clam" / "clam_per_slide_predictions.tsv"
CKPT_DIR = ROOT / "phase2_tcga_clam"
OUT_DIR = ROOT / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
PNG_OUT = OUT_DIR / "fig5_attention_heatmaps_v3.png"
PDF_OUT = OUT_DIR / "fig5_attention_heatmaps_v3.pdf"

THUMB_W = 1600
TILE_PX_FULL = 256


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
        return self.classifier(bag), attn


def load_fold_model(fold: int) -> GatedAttentionMIL:
    state = torch.load(CKPT_DIR / f"clam_fold{fold}_best.pt", map_location="cpu")
    m = GatedAttentionMIL(); m.load_state_dict(state); m.eval()
    return m


def get_attention(model, feat_path):
    feat = torch.load(feat_path, map_location="cpu")
    if isinstance(feat, np.ndarray): feat = torch.from_numpy(feat)
    with torch.no_grad():
        logits, attn = model(feat.float())
        prob_dm1 = F.softmax(logits, dim=-1)[1].item()
    return attn.numpy(), prob_dm1


def find_coords_tsv(slide_id):
    for cand in (FEATURE_DIR / f"{slide_id}_coords.tsv", FEATURE_DIR / f"{slide_id}.coords.tsv"):
        if cand.exists(): return cand
    return None


def find_wsi(slide_id):
    cands = list(WSI_DIR.glob(f"{slide_id}*.svs")) + list(WSI_DIR.glob(f"*{slide_id}*.svs"))
    return cands[0] if cands else None


def render_thumb_and_heat(svs_path, coords_tsv, attention,
                          thumb_w=THUMB_W, tile_px_full=TILE_PX_FULL):
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
        return thumb_bgr, None, []

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
    # smaller smoothing — keep hotspots sharp
    heat = gaussian_filter(heat, sigma=tile_t * 0.3)
    # 75th percentile normalization for high-contrast spread (not raw max)
    pos = heat[heat > 0]
    if pos.size > 0:
        norm = np.quantile(pos, 0.75)
        if norm > 0:
            heat = np.clip(heat / norm, 0, 1)
    # find top-3 hotspots in image space (separated by ≥ 4 tile widths)
    hotspots = []
    work = heat.copy()
    for _ in range(3):
        if work.max() < 0.3:
            break
        yx = np.unravel_index(int(np.argmax(work)), work.shape)
        hotspots.append((yx[1], yx[0]))  # (x, y)
        # zero-out neighborhood
        rad = tile_t * 4
        y0 = max(0, yx[0] - rad); y1 = min(H, yx[0] + rad)
        x0 = max(0, yx[1] - rad); x1 = min(W, yx[1] + rad)
        work[y0:y1, x0:x1] = 0
    return thumb_bgr, heat, hotspots


def attention_overlay_strong(thumb_bgr, heat):
    """Per-pixel alpha blend — strong attention pixels become vivid color."""
    if heat is None:
        return thumb_bgr.copy()
    if heat.shape[:2] != thumb_bgr.shape[:2]:
        heat = cv2.resize(heat, (thumb_bgr.shape[1], thumb_bgr.shape[0]),
                          interpolation=cv2.INTER_LINEAR)
    heat_u8 = np.clip(heat * 255.0, 0, 255).astype(np.uint8)
    color = cv2.applyColorMap(heat_u8, cv2.COLORMAP_JET)  # JET = vivid red→blue
    # AGGRESSIVE per-pixel alpha = heat^0.35 so mid-attention pixels also pop
    # plus floor: any tile with attention >0 gets at least alpha=0.25
    alpha = np.power(heat, 0.35)
    alpha = np.where(heat > 0.02, np.maximum(alpha, 0.25), alpha)
    alpha = np.clip(alpha, 0, 0.92)[:, :, None]
    out = thumb_bgr.astype(np.float32) * (1 - alpha) + color.astype(np.float32) * alpha
    return out.clip(0, 255).astype(np.uint8)


def annotate_hotspots(panel_bgr, hotspots, scale=1.0,
                      ring_radius=60, ring_color=(255, 255, 255),
                      ring_color_outer=(0, 0, 0)):
    out = panel_bgr.copy()
    for i, (x, y) in enumerate(hotspots):
        # outer black ring + inner white ring + tiny inner dot
        cv2.circle(out, (int(x * scale), int(y * scale)),
                   int(ring_radius + 4), ring_color_outer, 4, cv2.LINE_AA)
        cv2.circle(out, (int(x * scale), int(y * scale)),
                   int(ring_radius), ring_color, 3, cv2.LINE_AA)
        cv2.circle(out, (int(x * scale), int(y * scale)), 6, (0, 0, 255), -1, cv2.LINE_AA)
        # rank number with halo
        num = str(i + 1)
        tx, ty = int(x * scale) + ring_radius + 12, int(y * scale) + 8
        cv2.putText(out, num, (tx, ty), cv2.FONT_HERSHEY_TRIPLEX, 1.4,
                    (0, 0, 0), 6, cv2.LINE_AA)
        cv2.putText(out, num, (tx, ty), cv2.FONT_HERSHEY_TRIPLEX, 1.4,
                    (255, 255, 255), 2, cv2.LINE_AA)
    return out


def make_panel(slide_id, prob, group, thumb_bgr, overlay_bgr, hotspots,
               attention, panel_w, panel_h):
    GREEN = (143, 157, 42); BLUE = (170, 133, 91); WHITE = (255, 255, 255)
    border = GREEN if group == "DM1+" else BLUE
    bp = 6
    hist_h = 90
    inner_w = panel_w - 2 * bp
    inner_h = panel_h - 2 * bp - hist_h
    half_w = inner_w // 2

    def fit(img):
        h, w = img.shape[:2]
        s = min(half_w / w, inner_h / h)
        nw, nh = int(w * s), int(h * s)
        r = cv2.resize(img, (nw, nh), interpolation=cv2.INTER_AREA)
        cv = np.full((inner_h, half_w, 3), 255, dtype=np.uint8)
        cv[(inner_h - nh) // 2:(inner_h - nh) // 2 + nh,
           (half_w - nw) // 2:(half_w - nw) // 2 + nw] = r
        return cv, s, ((half_w - nw) // 2, (inner_h - nh) // 2)

    left, _, _ = fit(thumb_bgr)
    right, scale_r, off_r = fit(overlay_bgr)
    # annotate hotspots on right panel — LARGE rings + arrows from corner
    if hotspots:
        annotated = right.copy()
        for i, (x, y) in enumerate(hotspots):
            cx = int(x * scale_r) + off_r[0]
            cy = int(y * scale_r) + off_r[1]
            # arrow from top-left corner pointing to hotspot
            arr_x0 = 25 + i * 35
            arr_y0 = 25
            cv2.arrowedLine(annotated, (arr_x0, arr_y0), (cx - 8, cy - 8),
                            (0, 0, 0), 5, cv2.LINE_AA, tipLength=0.06)
            cv2.arrowedLine(annotated, (arr_x0, arr_y0), (cx - 8, cy - 8),
                            (255, 255, 0), 3, cv2.LINE_AA, tipLength=0.06)
            # large concentric ring (60 px)
            cv2.circle(annotated, (cx, cy), 64, (0, 0, 0), 7, cv2.LINE_AA)
            cv2.circle(annotated, (cx, cy), 60, (255, 255, 0), 5, cv2.LINE_AA)
            cv2.circle(annotated, (cx, cy), 56, (0, 0, 0), 2, cv2.LINE_AA)
            cv2.circle(annotated, (cx, cy), 8, (0, 0, 255), -1, cv2.LINE_AA)
            cv2.circle(annotated, (cx, cy), 8, (255, 255, 255), 2, cv2.LINE_AA)
            # rank number large
            tx, ty = cx + 70, cy + 14
            cv2.putText(annotated, f"#{i+1}", (tx, ty),
                        cv2.FONT_HERSHEY_TRIPLEX, 1.6, (0, 0, 0), 9, cv2.LINE_AA)
            cv2.putText(annotated, f"#{i+1}", (tx, ty),
                        cv2.FONT_HERSHEY_TRIPLEX, 1.6, (0, 220, 255), 3, cv2.LINE_AA)
        right = annotated

    top = np.concatenate([left, right], axis=1)
    if top.shape[1] < inner_w:
        pad = np.full((inner_h, inner_w - top.shape[1], 3), 255, dtype=np.uint8)
        top = np.concatenate([top, pad], axis=1)
    cv2.line(top, (half_w - 1, 0), (half_w - 1, inner_h - 1),
             (200, 200, 200), 1, cv2.LINE_AA)

    # bottom hist + inline colorbar
    hist_strip = np.full((hist_h, inner_w, 3), 250, dtype=np.uint8)
    if attention is not None and len(attention) > 1:
        n_bins = 60
        a = attention / (attention.max() + 1e-9)
        counts, _ = np.histogram(a, bins=n_bins, range=(0, 1))
        if counts.max() > 0:
            counts = counts / counts.max()
        bar_w = (inner_w - 200) / n_bins  # leave 200px for colorbar legend
        baseline_y = hist_h - 18
        max_bar_h = hist_h - 36
        for i, c in enumerate(counts):
            x0 = int(i * bar_w); x1 = int((i + 1) * bar_w) - 1
            bh = int(c * max_bar_h)
            t = (i + 0.5) / n_bins
            bar_col = cv2.applyColorMap(
                np.array([[int(t * 255)]], dtype=np.uint8),
                cv2.COLORMAP_JET)[0, 0].tolist()
            cv2.rectangle(hist_strip, (x0, baseline_y - bh), (x1, baseline_y),
                          bar_col, -1, cv2.LINE_AA)
        cv2.line(hist_strip, (0, baseline_y), (inner_w - 200 - 1, baseline_y),
                 (180, 180, 180), 1, cv2.LINE_AA)
        cv2.putText(hist_strip,
                    "Attention weight (low -> high)  |  JET colorbar:",
                    (8, hist_h - 4),
                    cv2.FONT_HERSHEY_DUPLEX, 0.45, (90, 90, 90), 1, cv2.LINE_AA)
        # right-side mini colorbar gradient
        cb_x0 = inner_w - 190
        cb_x1 = inner_w - 12
        cb_y0 = 18
        cb_y1 = baseline_y - 4
        for j, x in enumerate(range(cb_x0, cb_x1)):
            t = (x - cb_x0) / (cb_x1 - cb_x0)
            col = cv2.applyColorMap(
                np.array([[int(t * 255)]], dtype=np.uint8),
                cv2.COLORMAP_JET)[0, 0].tolist()
            cv2.line(hist_strip, (x, cb_y0), (x, cb_y1), col, 1)
        cv2.rectangle(hist_strip, (cb_x0 - 1, cb_y0 - 1), (cb_x1, cb_y1 + 1),
                      (60, 60, 60), 1, cv2.LINE_AA)
        cv2.putText(hist_strip, "low", (cb_x0 - 2, cb_y1 + 14),
                    cv2.FONT_HERSHEY_DUPLEX, 0.42, (90, 90, 90), 1, cv2.LINE_AA)
        cv2.putText(hist_strip, "high", (cb_x1 - 26, cb_y1 + 14),
                    cv2.FONT_HERSHEY_DUPLEX, 0.42, (90, 90, 90), 1, cv2.LINE_AA)

    inner = np.concatenate([top, hist_strip], axis=0)
    panel = np.full((panel_h, panel_w, 3), 255, dtype=np.uint8)
    panel[bp:bp + inner.shape[0], bp:bp + inner.shape[1]] = inner
    cv2.rectangle(panel, (0, 0), (panel_w - 1, panel_h - 1),
                  border, bp, cv2.LINE_AA)

    # group badge
    badge_text = group
    (tw_, th_), _ = cv2.getTextSize(badge_text, cv2.FONT_HERSHEY_TRIPLEX, 0.95, 2)
    bx0, by0 = bp + 14, bp + 14
    bx1, by1 = bx0 + tw_ + 22, by0 + th_ + 16
    o = panel.copy()
    cv2.rectangle(o, (bx0, by0), (bx1, by1), border, -1, cv2.LINE_AA)
    panel = cv2.addWeighted(o, 0.92, panel, 0.08, 0)
    cv2.putText(panel, badge_text, (bx0 + 11, by1 - 10),
                cv2.FONT_HERSHEY_TRIPLEX, 0.95, WHITE, 2, cv2.LINE_AA)

    # H&E / Attention sub-labels
    cv2.putText(panel, "H&E (raw)", (bp + 12, bp + th_ + 60),
                cv2.FONT_HERSHEY_DUPLEX, 0.7, (60, 60, 60), 1, cv2.LINE_AA)
    cv2.putText(panel, "Attention (JET) + top-3 hotspots",
                (bp + half_w + 12, bp + th_ + 60),
                cv2.FONT_HERSHEY_DUPLEX, 0.7, (60, 60, 60), 1, cv2.LINE_AA)

    # prob badge
    short_id = slide_id[:8] + "…"
    prob_text = f"DM1 prob: {prob:.2f}"
    (pw, ph), _ = cv2.getTextSize(prob_text, cv2.FONT_HERSHEY_TRIPLEX, 1.05, 2)
    (sw, sh), _ = cv2.getTextSize(short_id, cv2.FONT_HERSHEY_DUPLEX, 0.55, 1)
    pill_w = max(pw, sw) + 28; pill_h = ph + sh + 26
    px1 = panel_w - bp - 12
    py1 = bp + inner_h - 12
    px0 = px1 - pill_w; py0 = py1 - pill_h
    o = panel.copy()
    cv2.rectangle(o, (px0, py0), (px1, py1), (255, 255, 255), -1, cv2.LINE_AA)
    panel = cv2.addWeighted(o, 0.92, panel, 0.08, 0)
    cv2.rectangle(panel, (px0, py0), (px1, py1), border, 2, cv2.LINE_AA)
    cv2.putText(panel, short_id, (px0 + 14, py0 + sh + 8),
                cv2.FONT_HERSHEY_DUPLEX, 0.55, (90, 90, 90), 1, cv2.LINE_AA)
    cv2.putText(panel, prob_text, (px0 + 14, py1 - 10),
                cv2.FONT_HERSHEY_TRIPLEX, 1.05, border, 2, cv2.LINE_AA)
    return panel


def make_section_header(text, w, h, bg_bgr):
    band = np.full((h, w, 3), bg_bgr, dtype=np.uint8)
    grad = np.linspace(1.0, 0.85, w, dtype=np.float32).reshape(1, w, 1)
    band = (band.astype(np.float32) * grad).clip(0, 255).astype(np.uint8)
    cv2.putText(band, text, (32, h // 2 + 16),
                cv2.FONT_HERSHEY_TRIPLEX, 1.4, (255, 255, 255), 3, cv2.LINE_AA)
    cv2.line(band, (0, h - 1), (w - 1, h - 1), (60, 60, 60), 1, cv2.LINE_AA)
    return band


def main():
    print("[fig5_v3] loading...")
    preds = pd.read_csv(PREDS_TSV, sep="\t")

    def coords_ok(s):
        return ((FEATURE_DIR / f"{s}_coords.tsv").exists()
                or (FEATURE_DIR / f"{s}.coords.tsv").exists())

    def svs_ok(s): return find_wsi(s) is not None

    preds = preds[preds["slide"].apply(lambda s: coords_ok(s) and svs_ok(s))].copy()

    def svs_readable(s):
        try:
            slide = openslide.OpenSlide(str(find_wsi(s))); slide.close(); return True
        except Exception: return False

    preds["readable"] = preds["slide"].apply(svs_readable)
    preds = preds[preds["readable"]].copy()

    dm1_pos = (preds[preds["label"] == 1].sort_values("prob_DM1", ascending=False)
               .drop_duplicates("slide", keep="first").head(3))
    dm2_neg = (preds[preds["label"] == 0].sort_values("prob_DM1", ascending=True)
               .drop_duplicates("slide", keep="first").head(3))
    picks = pd.concat([dm1_pos.assign(group="DM1+"),
                       dm2_neg.assign(group="DM2")], ignore_index=True)
    print(f"[fig5_v3] {len(dm1_pos)} DM1+ / {len(dm2_neg)} DM2")

    fold_models = {}
    panels = []
    panel_w, panel_h = 1230, 760
    for _, row in picks.iterrows():
        sid = row["slide"]; fold = int(row["fold"])
        prob = float(row["prob_DM1"]); group = row["group"]
        if fold not in fold_models:
            fold_models[fold] = load_fold_model(fold)
        feat_path = FEATURE_DIR / f"{sid}.pt"
        attention, _ = get_attention(fold_models[fold], feat_path)
        thumb_bgr, heat, hotspots = render_thumb_and_heat(
            find_wsi(sid), find_coords_tsv(sid), attention)
        overlay = attention_overlay_strong(thumb_bgr, heat)
        panel = make_panel(sid, prob, group, thumb_bgr, overlay, hotspots,
                           attention, panel_w, panel_h)
        panels.append((group, panel))
        print(f"  [{group}] {sid[:8]} prob={prob:.3f} hotspots={len(hotspots)}")

    title_h = 130; section_h = 80; footer_h = 80
    grid_w = panel_w * 3
    full_h = title_h + section_h + panel_h + section_h + panel_h + footer_h
    canvas = np.full((full_h, grid_w, 3), 252, dtype=np.uint8)

    title = np.zeros((title_h, grid_w, 3), dtype=np.uint8)
    for x in range(grid_w):
        t = x / grid_w
        b = int(125 + (147 - 125) * t)
        g = int(107 + (76 - 107) * t)
        r = int(0 + (106 - 0) * t)
        title[:, x] = (b, g, r)
    cv2.putText(title,
                "Fig 5 v3  -  CLAM gated-attention heatmaps  (top-3 DM1+ / top-3 DM2)",
                (40, 50), cv2.FONT_HERSHEY_TRIPLEX, 1.4, (255, 255, 255), 3, cv2.LINE_AA)
    cv2.putText(title,
                "JET colormap | per-pixel alpha | 95th-pct norm | top-3 hotspots ringed",
                (40, 95), cv2.FONT_HERSHEY_DUPLEX, 0.85,
                (230, 230, 230), 1, cv2.LINE_AA)
    canvas[:title_h, :, :] = title

    HDR_DM1 = (90, 130, 50); HDR_DM2 = (130, 90, 60)
    y = title_h
    canvas[y:y + section_h, :, :] = make_section_header(
        "DM1+ HIGH CONFIDENCE   (predicted DM1, label = DM1) -- attention concentrates on dedifferentiated zones",
        grid_w, section_h, HDR_DM1)
    y += section_h
    for i, p in enumerate([p for g, p in panels if g == "DM1+"]):
        canvas[y:y + panel_h, i * panel_w:(i + 1) * panel_w] = p
    y += panel_h
    canvas[y:y + section_h, :, :] = make_section_header(
        "DM2 LOW CONFIDENCE   (predicted DM2, label = DM2) -- attention diffuse / on well-differentiated follicular regions",
        grid_w, section_h, HDR_DM2)
    y += section_h
    for i, p in enumerate([p for g, p in panels if g == "DM2"]):
        canvas[y:y + panel_h, i * panel_w:(i + 1) * panel_w] = p
    y += panel_h

    cv2.putText(canvas,
                "Heatmap = JET (low blue -> high red), per-pixel alpha so high-attention regions are vivid. "
                "Numbered rings = top-3 attention hotspots per slide. "
                "Histogram = per-tile attention distribution; right-side mini colorbar.",
                (40, y + 50),
                cv2.FONT_HERSHEY_DUPLEX, 0.7, (60, 60, 60), 1, cv2.LINE_AA)

    target_w = 3840
    scale = target_w / canvas.shape[1]
    canvas_4k = cv2.resize(canvas, (target_w, int(canvas.shape[0] * scale)),
                           interpolation=cv2.INTER_AREA)
    cv2.imwrite(str(PNG_OUT), canvas_4k, [cv2.IMWRITE_PNG_COMPRESSION, 4])
    print(f"[saved] {PNG_OUT} ({canvas_4k.shape[1]}x{canvas_4k.shape[0]})")

    rgb = cv2.cvtColor(canvas_4k, cv2.COLOR_BGR2RGB)
    fig_w_in = 16.0
    fig_h_in = fig_w_in * rgb.shape[0] / rgb.shape[1]
    fig = plt.figure(figsize=(fig_w_in, fig_h_in), dpi=300)
    ax = fig.add_axes([0, 0, 1, 1]); ax.imshow(rgb); ax.axis("off")
    fig.savefig(PDF_OUT, dpi=300, bbox_inches="tight", pad_inches=0)
    plt.close(fig)
    print(f"[saved] {PDF_OUT}")


if __name__ == "__main__":
    main()
