"""
Fig 5 — CLAM attention heatmaps for Paper 2 (H&E -> DM1 image classifier).

For the top-3 highest-confidence DM1+ slides AND top-3 highest-confidence DM2 slides
(picked from the held-out predictions tsv produced by phase2_clam_train_eval.py),
re-runs the trained CLAM model to get per-tile gated-attention weights, and renders
side-by-side panels:
  - left:  H&E thumbnail (1024 px wide)
  - right: same thumbnail with attention heatmap overlay (red high, blue low; alpha 0.6)

Composites everything into a 6x2 grid and saves
  <out_dir>/fig5_attention_heatmaps.{png,pdf}  at 300 dpi.

Usage (on pod):
  python make_fig5_attention.py \
    --feature_dir   /workspace/.../uni_feats \
    --wsi_dir       /workspace/.../svs \
    --manifest      .../phase2_tcga_clam/slide_manifest.tsv \
    --predictions   .../phase2_tcga_clam/clam_per_slide_predictions.tsv \
    --ckpt_dir      .../phase2_tcga_clam \
    --out_dir       .../phase2_tcga_clam/figures

Notes:
  - The CLAM checkpoint must have been saved by phase2_clam_train_eval.py
    (clam_fold{N}_best.pt + clam_fold_summary.tsv). If absent, this script
    raises a clear error.
  - Tile coords TSV (per-slide, columns x_full y_full at slide-level pixels)
    is expected at <feature_dir>/<file_id>.coords.tsv ; if missing the panel
    falls back to thumbnail-only with a 'no coords' note.
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
import matplotlib.pyplot as plt
from matplotlib import cm
from PIL import Image
from scipy.ndimage import gaussian_filter

try:
    import openslide
except ImportError as e:
    sys.exit(f"openslide-python required: pip install openslide-python ({e})")


# ============================================================
# CLAM gated-attention MIL  (verbatim from phase2_clam_train_eval.py)
# ============================================================
class GatedAttentionMIL(nn.Module):
    def __init__(self, in_dim=1024, hidden=512, n_classes=2, dropout=0.25):
        super().__init__()
        self.fc = nn.Sequential(nn.Linear(in_dim, hidden), nn.ReLU(), nn.Dropout(dropout))
        self.attn_a = nn.Sequential(nn.Linear(hidden, hidden), nn.Tanh())
        self.attn_b = nn.Sequential(nn.Linear(hidden, hidden), nn.Sigmoid())
        self.attn_w = nn.Linear(hidden, 1)
        self.classifier = nn.Linear(hidden, n_classes)

    def forward(self, x):
        # x: (N_tiles, in_dim)
        h = self.fc(x)
        a = self.attn_w(self.attn_a(h) * self.attn_b(h))  # (N, 1)
        attn = F.softmax(a.squeeze(-1), dim=0)            # (N,)
        bag = (attn.unsqueeze(-1) * h).sum(0)             # (hidden,)
        logits = self.classifier(bag)                     # (n_classes,)
        return logits, attn


# ============================================================
# Helpers
# ============================================================
def load_model(ckpt_path: Path, in_dim: int = 1024, hidden: int = 512,
               device: str = "cuda") -> GatedAttentionMIL:
    if not ckpt_path.exists():
        raise FileNotFoundError(
            f"CLAM checkpoint not found at {ckpt_path}. "
            "CLAM checkpoint not saved by phase2_clam_train_eval.py — "
            "modify training script to save best fold weights, then rerun."
        )
    state = torch.load(ckpt_path, map_location="cpu")
    m = GatedAttentionMIL(in_dim=in_dim, hidden=hidden, n_classes=2)
    m.load_state_dict(state)
    m.eval()
    return m.to(device)


def pick_best_ckpt(ckpt_dir: Path) -> Path:
    """Read clam_fold_summary.tsv and return checkpoint path with highest val AUC."""
    summ = ckpt_dir / "clam_fold_summary.tsv"
    if not summ.exists():
        # fall back to fold1 if summary absent but ckpts present
        cands = sorted(ckpt_dir.glob("clam_fold*_best.pt"))
        if not cands:
            raise FileNotFoundError(
                f"Neither clam_fold_summary.tsv nor clam_fold*_best.pt found in {ckpt_dir}. "
                "CLAM checkpoint not saved by phase2_clam_train_eval.py — "
                "modify training script to save best fold weights, then rerun."
            )
        return cands[0]
    df = pd.read_csv(summ, sep="\t")
    best = df.sort_values("best_val_auc", ascending=False).iloc[0]
    return ckpt_dir / best["ckpt"]


def get_attention(model: GatedAttentionMIL, feat_path: Path, device: str = "cuda"):
    feat = torch.load(feat_path, map_location="cpu")
    if isinstance(feat, np.ndarray):
        feat = torch.from_numpy(feat)
    x = feat.float().to(device)
    with torch.no_grad():
        logits, attn = model(x)
        prob_dm1 = F.softmax(logits, dim=-1)[1].item()
    return attn.cpu().numpy(), prob_dm1


def find_coords_tsv(feature_dir: Path, file_id: str) -> Path | None:
    """Look for tile-coord TSV. Try a few common name patterns."""
    for cand in (
        feature_dir / f"{file_id}.coords.tsv",
        feature_dir / f"{file_id}_coords.tsv",
        feature_dir / "coords" / f"{file_id}.tsv",
    ):
        if cand.exists():
            return cand
    return None


def find_wsi(wsi_dir: Path, file_id: str, file_name: str | None = None) -> Path | None:
    """Locate the .svs file for a given file_id (TCGA gdc layout uses {file_id}/{file_name})."""
    cands = []
    if file_name:
        cands += [wsi_dir / file_id / file_name, wsi_dir / file_name]
    cands += list((wsi_dir / file_id).glob("*.svs")) if (wsi_dir / file_id).exists() else []
    cands += list(wsi_dir.glob(f"{file_id}*.svs"))
    cands += list(wsi_dir.glob(f"*{file_id}*.svs"))
    for c in cands:
        if c and c.exists():
            return c
    return None


def render_heatmap(svs_path: Path, coords_tsv: Path | None, attention: np.ndarray,
                   thumb_size: int = 1024, tile_px_full: int = 256):
    """
    Returns (thumb_rgb_pil, heat_2d_np or None, max_xy_full or None).
    heat_2d_np is float in [0,1] aligned to the thumbnail; None if coords missing.
    max_xy_full is (x,y) in slide-level pixels of the highest-attention tile.
    """
    slide = openslide.OpenSlide(str(svs_path))
    w_full, h_full = slide.dimensions
    ratio = thumb_size / max(w_full, h_full)
    tw, th = int(w_full * ratio), int(h_full * ratio)
    thumb = slide.get_thumbnail((tw, th)).convert("RGB")
    slide.close()

    if coords_tsv is None or not coords_tsv.exists():
        return thumb, None, None

    coords = pd.read_csv(coords_tsv, sep="\t")
    # accept either 'x_full,y_full' headers or unheaded x,y
    if "x_full" in coords.columns and "y_full" in coords.columns:
        xy = coords[["x_full", "y_full"]].values
    elif coords.shape[1] >= 2:
        xy = coords.iloc[:, :2].values
    else:
        return thumb, None, None

    if len(xy) != len(attention):
        # mismatch — bail out gracefully (probably stale coords from old run)
        return thumb, None, None

    # paint each tile footprint with its attention value (max-merge on overlap)
    H, W = th, tw
    heat = np.zeros((H, W), dtype=np.float32)
    tile_t = max(1, int(tile_px_full * ratio))
    for (x, y), a in zip(xy, attention):
        xt = int(round(x * ratio))
        yt = int(round(y * ratio))
        x0 = max(0, xt); y0 = max(0, yt)
        x1 = min(W, xt + tile_t); y1 = min(H, yt + tile_t)
        if x1 > x0 and y1 > y0:
            np.maximum(heat[y0:y1, x0:x1], float(a), out=heat[y0:y1, x0:x1])

    # smooth a bit so neighbouring tiles blend
    heat = gaussian_filter(heat, sigma=tile_t * 0.5)
    if heat.max() > 0:
        heat = heat / heat.max()

    # top-attention tile coords in slide-level pixels
    top_idx = int(np.argmax(attention))
    max_xy_full = (int(xy[top_idx, 0]), int(xy[top_idx, 1]))
    return thumb, heat, max_xy_full


def overlay_heatmap_on_thumb(thumb_pil: Image.Image, heat: np.ndarray | None,
                             alpha: float = 0.6, cmap_name: str = "jet") -> Image.Image:
    if heat is None:
        return thumb_pil.copy()
    cmap = cm.get_cmap(cmap_name)
    rgba = (cmap(heat) * 255).astype(np.uint8)         # (H, W, 4)
    heat_pil = Image.fromarray(rgba, mode="RGBA")
    base = thumb_pil.convert("RGBA")
    # only overlay where heat > small threshold so background stays clean
    mask = (heat * 255 * alpha).astype(np.uint8)
    mask = np.where(heat > 0.05, mask, 0).astype(np.uint8)
    heat_pil.putalpha(Image.fromarray(mask, mode="L"))
    out = Image.alpha_composite(base, heat_pil)
    return out.convert("RGB")


# ============================================================
# Main
# ============================================================
def main():
    p = argparse.ArgumentParser()
    p.add_argument("--feature_dir", type=Path, required=True,
                   help="Per-slide UNI feature .pt files (and optional .coords.tsv)")
    p.add_argument("--wsi_dir", type=Path, required=True,
                   help="Directory containing .svs files (gdc layout: {file_id}/{file_name}.svs)")
    p.add_argument("--manifest", type=Path, required=True)
    p.add_argument("--predictions", type=Path, required=True,
                   help="clam_per_slide_predictions.tsv produced by phase2_clam_train_eval.py")
    p.add_argument("--ckpt", type=Path, default=None,
                   help="CLAM checkpoint .pt; if not given, picks best fold from --ckpt_dir")
    p.add_argument("--ckpt_dir", type=Path, default=None,
                   help="Dir with clam_fold*_best.pt + clam_fold_summary.tsv")
    p.add_argument("--out_dir", type=Path, required=True)
    p.add_argument("--top_k", type=int, default=3)
    p.add_argument("--thumb_size", type=int, default=1024)
    p.add_argument("--tile_px_full", type=int, default=256,
                   help="Tile size in slide-level pixels (CLAM 20x default = 256)")
    p.add_argument("--device", default="cuda")
    args = p.parse_args()

    device = args.device if torch.cuda.is_available() else "cpu"
    args.out_dir.mkdir(parents=True, exist_ok=True)

    # --- resolve checkpoint ---
    if args.ckpt is None:
        if args.ckpt_dir is None:
            sys.exit("Must provide --ckpt or --ckpt_dir.")
        try:
            ckpt_path = pick_best_ckpt(args.ckpt_dir)
        except FileNotFoundError as e:
            sys.exit(str(e))
    else:
        ckpt_path = args.ckpt

    print(f"[ckpt] using {ckpt_path}")
    try:
        model = load_model(ckpt_path, device=device)
    except FileNotFoundError as e:
        sys.exit(str(e))

    # --- pick top slides ---
    preds = pd.read_csv(args.predictions, sep="\t")
    man = pd.read_csv(args.manifest, sep="\t")
    name_col = "file_id" if "file_id" in man.columns else man.columns[0]
    file_name_col = "file_name" if "file_name" in man.columns else None
    man_lookup = man.set_index(name_col)

    dm1_pos = preds[preds["label"] == 1].sort_values("prob_DM1", ascending=False).head(args.top_k)
    dm2_neg = preds[preds["label"] == 0].sort_values("prob_DM1", ascending=True).head(args.top_k)
    picks = pd.concat([dm1_pos.assign(group="DM1+"),
                       dm2_neg.assign(group="DM2")], ignore_index=True)
    print(f"[pick] {len(picks)} slides ({len(dm1_pos)} DM1+ / {len(dm2_neg)} DM2)")

    n_rows = len(picks)
    fig, axes = plt.subplots(n_rows, 2, figsize=(10, 4 * n_rows))
    if n_rows == 1:
        axes = np.array([axes])

    rendered = 0
    for i, row in picks.reset_index(drop=True).iterrows():
        slide_id = row["slide"]
        prob = float(row["prob_DM1"])
        group = row["group"]
        feat_path = args.feature_dir / f"{slide_id}.pt"
        coords_tsv = find_coords_tsv(args.feature_dir, slide_id)

        file_name = None
        if file_name_col and slide_id in man_lookup.index:
            file_name = str(man_lookup.loc[slide_id, file_name_col])
        svs_path = find_wsi(args.wsi_dir, slide_id, file_name)

        ax_l, ax_r = axes[i, 0], axes[i, 1]
        for ax in (ax_l, ax_r):
            ax.set_xticks([]); ax.set_yticks([])

        if not feat_path.exists():
            ax_l.text(0.5, 0.5, f"missing feat\n{slide_id[:12]}", ha="center", va="center")
            ax_r.text(0.5, 0.5, "no attention", ha="center", va="center")
            continue
        if svs_path is None:
            ax_l.text(0.5, 0.5, f"missing svs\n{slide_id[:12]}", ha="center", va="center")
            ax_r.text(0.5, 0.5, "no thumbnail", ha="center", va="center")
            continue

        try:
            attention, prob_recompute = get_attention(model, feat_path, device=device)
        except Exception as e:
            ax_l.text(0.5, 0.5, f"attn fail\n{e}", ha="center", va="center", fontsize=7)
            ax_r.text(0.5, 0.5, "—", ha="center", va="center")
            continue

        thumb, heat, max_xy = render_heatmap(
            svs_path, coords_tsv, attention,
            thumb_size=args.thumb_size, tile_px_full=args.tile_px_full,
        )
        overlay = overlay_heatmap_on_thumb(thumb, heat, alpha=0.6)

        ax_l.imshow(thumb)
        ax_r.imshow(overlay)

        max_str = f"max@({max_xy[0]},{max_xy[1]})" if max_xy is not None else "no coords"
        title = (f"{group}  {slide_id[:8]}…  "
                 f"DM1 prob={prob:.3f}  (recomp={prob_recompute:.3f})  {max_str}")
        ax_l.set_title(title, fontsize=8, loc="left")
        ax_r.set_title("attention overlay (red=high, blue=low)", fontsize=8, loc="left")
        rendered += 1

    fig.suptitle("Fig 5 — CLAM attention heatmaps (top 3 DM1+ / top 3 DM2)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.98])

    png_out = args.out_dir / "fig5_attention_heatmaps.png"
    pdf_out = args.out_dir / "fig5_attention_heatmaps.pdf"
    fig.savefig(png_out, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_out, bbox_inches="tight")
    plt.close(fig)

    print(f"[saved] {png_out}")
    print(f"[saved] {pdf_out}")
    print(f"[done] rendered {rendered}/{n_rows} slides")


if __name__ == "__main__":
    main()
