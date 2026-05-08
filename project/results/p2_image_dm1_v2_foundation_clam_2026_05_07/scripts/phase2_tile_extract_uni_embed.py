"""
Phase 2 missing piece: CLAM-standard tile extraction + UNI embedding for TCGA WSI.

Workflow per slide:
  1. Open SVS via openslide
  2. Generate tissue mask (Otsu on saturation channel; CLAM default)
  3. Extract non-overlapping 256x256 tiles at 20× (level chosen by mpp 0.5)
  4. Filter tiles with tissue fraction >= 0.5
  5. UNI embedding per tile (1024-d)
  6. Save (N_tiles, 1024) feature tensor as <file_id>.pt + tile coords TSV

Output: features/<file_id>.pt + features/<file_id>_coords.tsv

Resource: per slide ~500-3000 tiles; UNI on L40S/A100 ~1-3 min/slide; 90 slides ~1-3 hr.
"""
from __future__ import annotations
import argparse, gc, os, time
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from PIL import Image
import openslide
from skimage.filters import threshold_otsu
from skimage.color import rgb2hsv
from skimage.morphology import closing, opening, square
from torchvision import transforms

TILE = 256
TARGET_MPP = 0.5  # 20× equivalent
TISSUE_FRACTION_MIN = 0.5

def get_level_for_mpp(slide, target_mpp=TARGET_MPP, tol=0.15):
    """Find best level matching target microns-per-pixel."""
    base_mpp = float(slide.properties.get("openslide.mpp-x", 0.25))
    for level in range(slide.level_count):
        level_mpp = base_mpp * slide.level_downsamples[level]
        if abs(level_mpp - target_mpp) / target_mpp < tol:
            return level, level_mpp
    # fallback: closest
    best = min(range(slide.level_count),
               key=lambda l: abs(base_mpp*slide.level_downsamples[l] - target_mpp))
    return best, base_mpp*slide.level_downsamples[best]

def tissue_mask(thumb_rgb: np.ndarray) -> np.ndarray:
    """Otsu on saturation channel; CLAM default."""
    hsv = rgb2hsv(thumb_rgb)
    sat = (hsv[:, :, 1] * 255).astype(np.uint8)
    try:
        t = threshold_otsu(sat)
    except Exception:
        t = 30
    mask = sat > t
    mask = closing(mask, square(5))
    mask = opening(mask, square(3))
    return mask

def extract_tiles(slide_path: Path, tile_size=TILE):
    """Yield (tile_pil, x, y) at level matching 20× (TARGET_MPP)."""
    slide = openslide.OpenSlide(str(slide_path))
    level, mpp = get_level_for_mpp(slide)
    w, h = slide.level_dimensions[level]
    # thumbnail for tissue mask
    thumb_w = max(1, w // 32)
    thumb_h = max(1, h // 32)
    thumb = np.array(slide.get_thumbnail((thumb_w, thumb_h)).convert("RGB"))
    mask = tissue_mask(thumb)

    # downsample factor between level and thumbnail
    sx = w / mask.shape[1]
    sy = h / mask.shape[0]

    base_ds = int(slide.level_downsamples[level])
    n_tiles = 0
    for y in range(0, h - tile_size, tile_size):
        for x in range(0, w - tile_size, tile_size):
            mx0 = int(x / sx); mx1 = int((x + tile_size) / sx)
            my0 = int(y / sy); my1 = int((y + tile_size) / sy)
            sub = mask[my0:my1+1, mx0:mx1+1]
            if sub.size == 0: continue
            if sub.mean() < TISSUE_FRACTION_MIN: continue
            # absolute coordinates at level 0 = level * downsample
            tile = slide.read_region((x * base_ds, y * base_ds), level,
                                      (tile_size, tile_size)).convert("RGB")
            yield np.array(tile), x, y, mpp
            n_tiles += 1
    slide.close()
    return n_tiles

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--wsi_dir", type=Path, required=True,
                        help="Directory of TCGA SVS files")
    parser.add_argument("--manifest", type=Path, required=True,
                        help="slide_manifest.tsv with file_id, file_name, dm")
    parser.add_argument("--out_dir", type=Path, required=True,
                        help="Output directory for features")
    parser.add_argument("--batch_size", type=int, default=64)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--max_tiles_per_slide", type=int, default=4000,
                        help="cap tiles per slide to bound memory")
    parser.add_argument("--limit_slides", type=int, default=None)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    feat_dir = args.out_dir / "features"
    feat_dir.mkdir(exist_ok=True)

    # --- load UNI ---
    print("[load] UNI from MahmoodLab/UNI ...")
    import timm
    from huggingface_hub import hf_hub_download
    try:
        model = timm.create_model("vit_large_patch16_224", img_size=224, pretrained=False,
                                   num_classes=0, dynamic_img_size=True)
        weights_path = hf_hub_download(repo_id="MahmoodLab/UNI", filename="pytorch_model.bin")
        state_dict = torch.load(weights_path, map_location="cpu")
        model.load_state_dict(state_dict, strict=True)
        print("  UNI weights loaded")
    except Exception as e:
        print(f"[fallback] UNI gated/unavailable: {e}")
        print("[fallback] using timm ViT-L pretrained on ImageNet (BASELINE only)")
        model = timm.create_model("vit_large_patch16_224", pretrained=True, num_classes=0)
    model = model.eval().to(args.device)

    transform = transforms.Compose([
        transforms.Resize(224),
        transforms.CenterCrop(224),
        transforms.ToTensor(),
        transforms.Normalize(mean=(0.485, 0.456, 0.406), std=(0.229, 0.224, 0.225)),
    ])

    manifest = pd.read_csv(args.manifest, sep="\t")
    print(f"[manifest] {len(manifest)} slides")
    if args.limit_slides:
        manifest = manifest.head(args.limit_slides)

    summary_rows = []
    for i, row in manifest.iterrows():
        out_pt = feat_dir / f"{row['file_id']}.pt"
        out_coords = feat_dir / f"{row['file_id']}_coords.tsv"
        if out_pt.exists() and out_coords.exists():
            print(f"  [{i+1}/{len(manifest)}] skip exist: {row['file_id']}")
            continue

        # find SVS file
        svs_path = args.wsi_dir / f"{row['file_id']}_{row['file_name']}"
        if not svs_path.exists():
            # fallback: any file matching file_id
            cands = list(args.wsi_dir.glob(f"{row['file_id']}*"))
            svs_path = cands[0] if cands else None
        if svs_path is None or not svs_path.exists():
            print(f"  [{i+1}/{len(manifest)}] MISSING WSI: {row['file_id']}")
            continue

        t0 = time.time()
        feats, coords_list = [], []
        batch_imgs, batch_xy = [], []

        try:
            n = 0
            for tile_arr, x, y, mpp in extract_tiles(svs_path):
                if n >= args.max_tiles_per_slide: break
                pil = Image.fromarray(tile_arr)
                batch_imgs.append(transform(pil))
                batch_xy.append((x, y, mpp))
                if len(batch_imgs) >= args.batch_size:
                    with torch.no_grad():
                        x_t = torch.stack(batch_imgs).to(args.device)
                        f = model(x_t).cpu().numpy()
                    feats.append(f)
                    coords_list.extend(batch_xy)
                    batch_imgs, batch_xy = [], []
                n += 1
            if batch_imgs:
                with torch.no_grad():
                    x_t = torch.stack(batch_imgs).to(args.device)
                    f = model(x_t).cpu().numpy()
                feats.append(f); coords_list.extend(batch_xy)

            if not feats:
                print(f"  [{i+1}/{len(manifest)}] EMPTY tiles: {row['file_id']}")
                continue

            feat_tensor = torch.from_numpy(np.vstack(feats)).float()
            torch.save(feat_tensor, out_pt)
            pd.DataFrame(coords_list, columns=["x", "y", "mpp"]).to_csv(out_coords, sep="\t", index=False)
            dt = time.time() - t0
            print(f"  [{i+1}/{len(manifest)}] {row['file_id']} → {feat_tensor.shape} in {dt:.1f}s")
            summary_rows.append({"file_id": row["file_id"], "dm": row.get("dm",""),
                                 "n_tiles": int(feat_tensor.shape[0]), "dt_sec": round(dt,1)})
            del feat_tensor; torch.cuda.empty_cache(); gc.collect()
        except Exception as e:
            print(f"  [{i+1}/{len(manifest)}] FAIL: {e}")
            summary_rows.append({"file_id": row["file_id"], "dm": row.get("dm",""),
                                 "n_tiles": 0, "dt_sec": 0, "error": str(e)[:200]})

    pd.DataFrame(summary_rows).to_csv(args.out_dir / "tile_extract_summary.tsv", sep="\t", index=False)
    print(f"[done] features in {feat_dir}")


if __name__ == "__main__":
    main()
