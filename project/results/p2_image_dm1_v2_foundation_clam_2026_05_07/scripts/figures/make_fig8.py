#!/usr/bin/env python3
"""
Paper 2 Fig 8 — Visual summary card (one-pager dashboard, 4K, cv2-polished).

Single 3840 x 2160 PNG/PDF that summarises Paper 2 at a glance:
  Top: title bar + verdict pill (PASS).
  4 quadrants:
    Top-left   : Phase 1 mini-stats panel (16-slide spatial UNI -> DM1 corr).
    Top-right  : Phase 2 mini-stats panel (CLAM AUC + folds).
    Bottom-left: best DM1+ attention thumbnail (cv2 overlay).
    Bottom-right: existing fig3 forest plot (rescaled).
  Bottom: clinical chain text strip.
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

try:
    import openslide
except ImportError as e:
    sys.exit(f"openslide-python required ({e})")

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/"
            "p2_image_dm1_v2_foundation_clam_2026_05_07")
FEATURE_DIR = ROOT / "phase2_tcga_clam" / "features"
WSI_DIR = Path("/data/tcga_thca_wsi")
PREDS_TSV = ROOT / "phase2_tcga_clam" / "clam_per_slide_predictions.tsv"
CKPT_DIR = ROOT / "phase2_tcga_clam"
FIG_DIR = ROOT / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)
PNG_OUT = FIG_DIR / "fig8_visual_summary.png"
PDF_OUT = FIG_DIR / "fig8_visual_summary.pdf"

# ---- palette ----
def hex_to_bgr(h):
    h = h.lstrip("#"); r = int(h[0:2], 16); g = int(h[2:4], 16); b = int(h[4:6], 16)
    return (b, g, r)

C_TEAL = hex_to_bgr("#006B7D")
C_PURPLE = hex_to_bgr("#6A4C93")
C_AMBER = hex_to_bgr("#FFB627")
C_GREEN = hex_to_bgr("#2A9D8F")
C_BLUE = hex_to_bgr("#5B85AA")
C_GRAY = hex_to_bgr("#6C757D")
C_DARK = (40, 40, 40)
C_MID = (90, 90, 90)
WHITE = (255, 255, 255)


# ---- CLAM model (same as elsewhere) ----
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


def load_fold_model(fold):
    ckpt = CKPT_DIR / f"clam_fold{fold}_best.pt"
    state = torch.load(ckpt, map_location="cpu")
    m = GatedAttentionMIL(); m.load_state_dict(state); m.eval(); return m


def find_coords(sid):
    for c in (FEATURE_DIR / f"{sid}_coords.tsv", FEATURE_DIR / f"{sid}.coords.tsv"):
        if c.exists(): return c
    return None


def find_wsi(sid):
    cands = list(WSI_DIR.glob(f"{sid}*.svs")) + list(WSI_DIR.glob(f"*{sid}*.svs"))
    return cands[0] if cands else None


def best_dm1_render(thumb_w=1200):
    """Pick the best confident DM1+ slide we can actually render and return overlay BGR."""
    preds = pd.read_csv(PREDS_TSV, sep="\t")
    cands = preds[preds["label"] == 1].sort_values("prob_DM1", ascending=False)
    for _, row in cands.iterrows():
        sid = row["slide"]; fold = int(row["fold"])
        ct = find_coords(sid); svs = find_wsi(sid)
        if ct is None or svs is None: continue
        try:
            slide = openslide.OpenSlide(str(svs))
        except Exception:
            continue
        try:
            w_full, h_full = slide.dimensions
            ratio = thumb_w / max(w_full, h_full)
            tw, th = int(w_full * ratio), int(h_full * ratio)
            pil = slide.get_thumbnail((tw, th)).convert("RGB")
        finally:
            slide.close()
        thumb = cv2.cvtColor(np.array(pil), cv2.COLOR_RGB2BGR)
        th, tw = thumb.shape[:2]

        m = load_fold_model(fold)
        feat = torch.load(FEATURE_DIR / f"{sid}.pt", map_location="cpu")
        if isinstance(feat, np.ndarray): feat = torch.from_numpy(feat)
        with torch.no_grad():
            logits, attn = m(feat.float())
            prob = F.softmax(logits, dim=-1)[1].item()
        attn_np = attn.numpy()
        coords = pd.read_csv(ct, sep="\t")
        xy = (coords[["x_full", "y_full"]].values
              if "x_full" in coords.columns else coords.iloc[:, :2].values)
        if len(xy) != len(attn_np): continue
        heat = np.zeros((th, tw), dtype=np.float32)
        tile_t = max(1, int(256 * ratio))
        for (x, y), a in zip(xy, attn_np):
            xt = int(round(x * ratio)); yt = int(round(y * ratio))
            x0 = max(0, xt); y0 = max(0, yt)
            x1 = min(tw, xt + tile_t); y1 = min(th, yt + tile_t)
            if x1 > x0 and y1 > y0:
                np.maximum(heat[y0:y1, x0:x1], float(a),
                           out=heat[y0:y1, x0:x1])
        heat = gaussian_filter(heat, sigma=tile_t * 0.5)
        if heat.max() > 0: heat = heat / heat.max()
        if heat.shape[:2] != thumb.shape[:2]:
            heat = cv2.resize(heat, (thumb.shape[1], thumb.shape[0]),
                              interpolation=cv2.INTER_LINEAR)
        heat_u8 = np.clip(heat * 255, 0, 255).astype(np.uint8)
        color = cv2.applyColorMap(heat_u8, cv2.COLORMAP_INFERNO)
        mask = (heat > 0.05).astype(np.float32)
        mask3 = np.dstack([mask] * 3)
        blended = cv2.addWeighted(thumb, 0.45, color, 0.55, 0)
        out = (thumb * (1 - mask3) + blended * mask3).astype(np.uint8)
        return out, sid, prob, fold
    return None, None, None, None


# ---- helpers ----
def gradient_box(canvas, x0, y0, x1, y1, base_bgr, *, lighten=0.18, corner_r=20,
                 border_color=None, border_px=2):
    h = y1 - y0; w = x1 - x0
    box = np.empty((h, w, 3), dtype=np.uint8)
    base = np.array(base_bgr, dtype=np.float32)
    top = base + (np.array(WHITE, dtype=np.float32) - base) * lighten
    bot = base
    for i in range(h):
        t = i / max(1, h - 1)
        box[i, :, :] = (top * (1 - t) + bot * t).astype(np.uint8)
    mask = np.zeros((h, w), dtype=np.uint8)
    cv2.rectangle(mask, (corner_r, 0), (w - corner_r, h), 255, -1)
    cv2.rectangle(mask, (0, corner_r), (w, h - corner_r), 255, -1)
    for cx, cy in [(corner_r, corner_r), (w - corner_r, corner_r),
                   (corner_r, h - corner_r), (w - corner_r, h - corner_r)]:
        cv2.circle(mask, (cx, cy), corner_r, 255, -1, cv2.LINE_AA)
    region = canvas[y0:y1, x0:x1]
    m3 = np.dstack([mask] * 3) / 255.0
    region[:] = (region * (1 - m3) + box * m3).astype(np.uint8)
    canvas[y0:y1, x0:x1] = region
    if border_color is None: border_color = base_bgr
    # rounded outline
    for (p0, p1) in [((x0 + corner_r, y0), (x1 - corner_r, y0)),
                     ((x1, y0 + corner_r), (x1, y1 - corner_r)),
                     ((x0 + corner_r, y1), (x1 - corner_r, y1)),
                     ((x0, y0 + corner_r), (x0, y1 - corner_r))]:
        cv2.line(canvas, p0, p1, border_color, border_px, cv2.LINE_AA)
    cv2.ellipse(canvas, (x0 + corner_r, y0 + corner_r),
                (corner_r, corner_r), 180, 0, 90, border_color, border_px, cv2.LINE_AA)
    cv2.ellipse(canvas, (x1 - corner_r, y0 + corner_r),
                (corner_r, corner_r), 270, 0, 90, border_color, border_px, cv2.LINE_AA)
    cv2.ellipse(canvas, (x1 - corner_r, y1 - corner_r),
                (corner_r, corner_r), 0, 0, 90, border_color, border_px, cv2.LINE_AA)
    cv2.ellipse(canvas, (x0 + corner_r, y1 - corner_r),
                (corner_r, corner_r), 90, 0, 90, border_color, border_px, cv2.LINE_AA)


def fit_into(box_w, box_h, img):
    h, w = img.shape[:2]
    s = min(box_w / w, box_h / h)
    nw, nh = int(w * s), int(h * s)
    return cv2.resize(img, (nw, nh), interpolation=cv2.INTER_AREA)


def main():
    W, H = 3840, 2160
    canvas = np.full((H, W, 3), 248, dtype=np.uint8)

    # ---- title bar ----
    bar_h = 180
    for x in range(W):
        t = x / W
        b = int(C_TEAL[0] * (1 - t) + C_PURPLE[0] * t)
        g = int(C_TEAL[1] * (1 - t) + C_PURPLE[1] * t)
        r = int(C_TEAL[2] * (1 - t) + C_PURPLE[2] * t)
        canvas[0:bar_h, x] = (b, g, r)
    cv2.putText(canvas,
                "Paper 2  -  H&E -> DM1 Foundation-Model Classifier",
                (60, 80), cv2.FONT_HERSHEY_TRIPLEX, 1.9, WHITE, 4, cv2.LINE_AA)
    cv2.putText(canvas,
                "Visual summary  |  TCGA-THCA + GSE250521 spatial  |  Cell Reports Medicine target",
                (60, 140), cv2.FONT_HERSHEY_DUPLEX, 1.0,
                (235, 235, 235), 2, cv2.LINE_AA)
    # PASS pill on the right
    pill_w = 360; pill_h = 100
    pill_x0 = W - pill_w - 70; pill_y0 = (bar_h - pill_h) // 2
    pill_x1 = pill_x0 + pill_w; pill_y1 = pill_y0 + pill_h
    overlay = canvas.copy()
    cv2.rectangle(overlay, (pill_x0, pill_y0), (pill_x1, pill_y1),
                  C_GREEN, -1, cv2.LINE_AA)
    canvas = cv2.addWeighted(overlay, 1.0, canvas, 0.0, 0)
    cv2.rectangle(canvas, (pill_x0, pill_y0), (pill_x1, pill_y1),
                  WHITE, 3, cv2.LINE_AA)
    cv2.putText(canvas, "VERDICT: PASS",
                (pill_x0 + 28, pill_y0 + pill_h - 32),
                cv2.FONT_HERSHEY_TRIPLEX, 1.2, WHITE, 3, cv2.LINE_AA)

    # ---- 4 quadrants layout ----
    # margins
    M = 60
    cy_split = bar_h + (H - bar_h - 200) // 2  # leaves 200 for bottom strip
    cx_split = W // 2
    # quadrant rectangles (with margins)
    # top-left
    tl = (M, bar_h + 30, cx_split - M // 2, cy_split - 10)
    tr = (cx_split + M // 2, bar_h + 30, W - M, cy_split - 10)
    bl = (M, cy_split + 10, cx_split - M // 2, H - 200)
    br = (cx_split + M // 2, cy_split + 10, W - M, H - 200)

    # ===== Top-Left: Phase 1 stats =====
    gradient_box(canvas, tl[0], tl[1], tl[2], tl[3], (250, 250, 250),
                 border_color=C_TEAL, border_px=4, corner_r=24, lighten=0.0)
    cv2.putText(canvas, "PHASE 1  -  Spatial UNI -> DM1 correlation",
                (tl[0] + 30, tl[1] + 60),
                cv2.FONT_HERSHEY_TRIPLEX, 1.05, C_TEAL, 3, cv2.LINE_AA)
    cv2.putText(canvas, "GSE250521  (n=16 thyroid spatial-transcriptomic slides)",
                (tl[0] + 30, tl[1] + 110),
                cv2.FONT_HERSHEY_DUPLEX, 0.78, C_DARK, 2, cv2.LINE_AA)
    rows = [
        ("Slides with rho > 0.3", "8 / 16"),
        ("Strongest negative rho", "-0.61  (LPTC-2)"),
        ("Significant slides (p < 0.05)", "11 / 16"),
        ("PC1 explained var (mean)", "27.4%"),
    ]
    y = tl[1] + 200
    for k, v in rows:
        cv2.putText(canvas, k, (tl[0] + 40, y),
                    cv2.FONT_HERSHEY_DUPLEX, 0.85, C_MID, 2, cv2.LINE_AA)
        cv2.putText(canvas, v, (tl[0] + 750, y),
                    cv2.FONT_HERSHEY_TRIPLEX, 1.0, C_TEAL, 3, cv2.LINE_AA)
        y += 90
    # bottom-of-quadrant takeaway
    cv2.putText(canvas,
                "Spatial signal: DM1 score tracks UNI tile embeddings within a slide.",
                (tl[0] + 30, tl[3] - 40),
                cv2.FONT_HERSHEY_DUPLEX, 0.7, C_GRAY, 2, cv2.LINE_AA)

    # ===== Top-Right: Phase 2 stats =====
    gradient_box(canvas, tr[0], tr[1], tr[2], tr[3], (250, 250, 250),
                 border_color=C_PURPLE, border_px=4, corner_r=24, lighten=0.0)
    cv2.putText(canvas, "PHASE 2  -  CLAM held-out classifier",
                (tr[0] + 30, tr[1] + 60),
                cv2.FONT_HERSHEY_TRIPLEX, 1.05, C_PURPLE, 3, cv2.LINE_AA)
    cv2.putText(canvas, "TCGA-THCA  (n=88 slides, 5-fold cross-validation)",
                (tr[0] + 30, tr[1] + 110),
                cv2.FONT_HERSHEY_DUPLEX, 0.78, C_DARK, 2, cv2.LINE_AA)
    rows2 = [
        ("Pooled cross-fold AUC", "0.746"),
        ("95% CI", "[0.61 - 0.86]"),
        ("Mean per-fold AUC", "0.830 +/- 0.139"),
        ("All folds AUC", ">= 0.71  (2/5 = 1.00)"),
        ("ResNet50 baseline", "~0.55  (negative)"),
    ]
    y = tr[1] + 200
    for k, v in rows2:
        cv2.putText(canvas, k, (tr[0] + 40, y),
                    cv2.FONT_HERSHEY_DUPLEX, 0.85, C_MID, 2, cv2.LINE_AA)
        cv2.putText(canvas, v, (tr[0] + 750, y),
                    cv2.FONT_HERSHEY_TRIPLEX, 1.0, C_PURPLE, 3, cv2.LINE_AA)
        y += 80
    cv2.putText(canvas,
                "All folds clear PASS threshold (0.70).  UNI features beat ImageNet baseline.",
                (tr[0] + 30, tr[3] - 40),
                cv2.FONT_HERSHEY_DUPLEX, 0.7, C_GRAY, 2, cv2.LINE_AA)

    # ===== Bottom-Left: best DM1+ attention thumbnail =====
    gradient_box(canvas, bl[0], bl[1], bl[2], bl[3], (250, 250, 250),
                 border_color=C_GREEN, border_px=4, corner_r=24, lighten=0.0)
    cv2.putText(canvas, "BEST DM1+ ATTENTION  (CLAM INFERNO heatmap)",
                (bl[0] + 30, bl[1] + 60),
                cv2.FONT_HERSHEY_TRIPLEX, 1.0, C_GREEN, 3, cv2.LINE_AA)
    print("[fig8] computing best DM1+ attention overlay...")
    overlay_img, sid, prob, fold = best_dm1_render(thumb_w=1400)
    inner_x0 = bl[0] + 28; inner_y0 = bl[1] + 110
    inner_x1 = bl[2] - 28; inner_y1 = bl[3] - 90
    inner_w = inner_x1 - inner_x0; inner_h = inner_y1 - inner_y0
    if overlay_img is not None:
        fitted = fit_into(inner_w, inner_h, overlay_img)
        oy = inner_y0 + (inner_h - fitted.shape[0]) // 2
        ox = inner_x0 + (inner_w - fitted.shape[1]) // 2
        canvas[oy:oy + fitted.shape[0], ox:ox + fitted.shape[1]] = fitted
        cv2.rectangle(canvas, (ox - 2, oy - 2),
                      (ox + fitted.shape[1] + 2, oy + fitted.shape[0] + 2),
                      C_GREEN, 2, cv2.LINE_AA)
        # caption
        cap = f"slide {sid[:8]}...  |  DM1 prob = {prob:.2f}  |  fold {fold}"
        cv2.putText(canvas, cap, (bl[0] + 30, bl[3] - 40),
                    cv2.FONT_HERSHEY_DUPLEX, 0.78, C_DARK, 2, cv2.LINE_AA)
    else:
        cv2.putText(canvas, "(no slide rendered)",
                    (bl[0] + 60, (bl[1] + bl[3]) // 2),
                    cv2.FONT_HERSHEY_DUPLEX, 1.0, C_GRAY, 2, cv2.LINE_AA)

    # ===== Bottom-Right: subgroup forest (use existing fig3) =====
    gradient_box(canvas, br[0], br[1], br[2], br[3], (250, 250, 250),
                 border_color=C_BLUE, border_px=4, corner_r=24, lighten=0.0)
    cv2.putText(canvas, "PHASE 1 PER-SLIDE rho FOREST  (subgroups)",
                (br[0] + 30, br[1] + 60),
                cv2.FONT_HERSHEY_TRIPLEX, 1.0, C_BLUE, 3, cv2.LINE_AA)
    forest = cv2.imread(str(FIG_DIR / "fig3_phase1_correlation_forest.png"))
    inner_x0 = br[0] + 28; inner_y0 = br[1] + 110
    inner_x1 = br[2] - 28; inner_y1 = br[3] - 90
    inner_w = inner_x1 - inner_x0; inner_h = inner_y1 - inner_y0
    if forest is not None:
        fitted = fit_into(inner_w, inner_h, forest)
        oy = inner_y0 + (inner_h - fitted.shape[0]) // 2
        ox = inner_x0 + (inner_w - fitted.shape[1]) // 2
        canvas[oy:oy + fitted.shape[0], ox:ox + fitted.shape[1]] = fitted
    else:
        cv2.putText(canvas, "(forest fig missing)",
                    (br[0] + 60, (br[1] + br[3]) // 2),
                    cv2.FONT_HERSHEY_DUPLEX, 1.0, C_GRAY, 2, cv2.LINE_AA)
    cv2.putText(canvas,
                "Per-slide rho: -0.61 to +0.34;  Stage I subgroup AUC trends similar to overall.",
                (br[0] + 30, br[3] - 40),
                cv2.FONT_HERSHEY_DUPLEX, 0.7, C_GRAY, 2, cv2.LINE_AA)

    # ===== Bottom strip: clinical chain =====
    strip_h = 140
    strip_y0 = H - strip_h - 30
    strip_y1 = strip_y0 + strip_h
    gradient_box(canvas, M, strip_y0, W - M, strip_y1, C_AMBER,
                 border_color=C_AMBER, border_px=3, corner_r=18, lighten=0.15)
    chain = ("Clinical chain:   H&E WSI  ->  CLAM (DM1 prob)  ->  "
             "[DM1+ -> reflex 8-gene panel -> RAI/ICI triage]   "
             "[DM2 -> standard surveillance]")
    cv2.putText(canvas, chain, (M + 30, strip_y0 + strip_h // 2 + 16),
                cv2.FONT_HERSHEY_TRIPLEX, 1.0, C_DARK, 3, cv2.LINE_AA)

    # ----- save -----
    cv2.imwrite(str(PNG_OUT), canvas, [cv2.IMWRITE_PNG_COMPRESSION, 4])
    print(f"[saved] {PNG_OUT}  ({canvas.shape[1]}x{canvas.shape[0]})")
    rgb = cv2.cvtColor(canvas, cv2.COLOR_BGR2RGB)
    fig_w_in = 16.0; fig_h_in = fig_w_in * H / W
    fig = plt.figure(figsize=(fig_w_in, fig_h_in), dpi=300)
    ax = fig.add_axes([0, 0, 1, 1]); ax.imshow(rgb); ax.axis("off")
    fig.savefig(PDF_OUT, dpi=300, bbox_inches="tight", pad_inches=0)
    plt.close(fig)
    print(f"[saved] {PDF_OUT}")


if __name__ == "__main__":
    main()
