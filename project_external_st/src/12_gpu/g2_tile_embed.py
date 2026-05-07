#!/usr/bin/env python3
"""G2 — Tile + embed TCGA-THCA SVS slides for digital pathology DM1 regression.
Stream mode: process one slide at a time, embed, save, delete SVS to save disk."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd


def tile_slide(svs_path, tile_size=224, overlap=0, otsu_threshold=True, max_tiles_per_slide=200):
    """Yield (tile_image, x, y) tuples for tissue regions of a SVS slide."""
    try:
        import openslide
    except ImportError:
        raise ImportError("Install openslide-python: pip install openslide-python")
    sl = openslide.OpenSlide(str(svs_path))
    # Use level 0 (highest mag) but at 224x224 px = ~110 µm
    W, H = sl.dimensions
    step = tile_size - overlap
    coords = []
    rng = np.random.default_rng(42)
    # sample candidate tile positions
    n_x = (W - tile_size) // step
    n_y = (H - tile_size) // step
    if n_x <= 0 or n_y <= 0:
        sl.close(); return iter([])
    candidates = [(x*step, y*step) for x in range(n_x) for y in range(n_y)]
    rng.shuffle(candidates)
    out = []
    for (x, y) in candidates:
        if len(out) >= max_tiles_per_slide: break
        try:
            tile = np.asarray(sl.read_region((x, y), 0, (tile_size, tile_size)).convert("RGB"))
        except Exception: continue
        # crude tissue filter: mean intensity in middle range
        mu = tile.mean()
        if mu < 30 or mu > 220: continue  # background or empty
        # luminance variance (avoid uniform regions)
        if tile.std() < 8: continue
        out.append((tile, x, y))
    sl.close()
    return iter(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", default="data/g2_tcga_he/manifest.tsv")
    ap.add_argument("--slides-dir", default="data/g2_tcga_he/slides")
    ap.add_argument("--out", default="results/g2_slide_embeddings.npz")
    ap.add_argument("--model", default="resnet50",
                    choices=["resnet50","uni","UNI","conch","CONCH"])
    ap.add_argument("--tile-size", type=int, default=224)
    ap.add_argument("--max-tiles", type=int, default=100)
    ap.add_argument("--device", default=None)
    ap.add_argument("--delete-after", action="store_true",
                    help="delete each SVS after embedding (stream mode for disk-constrained)")
    args = ap.parse_args()

    import torch
    if args.device is None:
        args.device = "cuda" if torch.cuda.is_available() else "cpu"
    sys_path = Path(__file__).parent
    import importlib.util
    spec = importlib.util.spec_from_file_location("g1_embed", sys_path / "g1_embed.py")
    g1m = importlib.util.module_from_spec(spec); spec.loader.exec_module(g1m)
    model, xform, dim = g1m.load_model(args.model, args.device)

    manifest = pd.read_csv(args.manifest, sep="\t")
    print(f"manifest: {len(manifest)} slides")

    slide_features = {}
    slide_meta = []
    for i, r in manifest.iterrows():
        slide_path = Path(args.slides_dir) / r["file_id"] / r["file_name"]
        if not slide_path.exists():
            print(f"[skip] {r['file_name']}: not downloaded"); continue
        print(f"[{i+1}/{len(manifest)}] {r['file_name']} (case {r['case_id']})...")
        try:
            tiles_iter = tile_slide(slide_path, args.tile_size,
                                     max_tiles_per_slide=args.max_tiles)
            tile_feats = []
            for (tile_img, x, y) in tiles_iter:
                from PIL import Image
                img = Image.fromarray(tile_img)
                t = xform(img).unsqueeze(0).to(args.device)
                with torch.no_grad():
                    if hasattr(model, "encode_image"):
                        f = model.encode_image(t)
                    else:
                        f = model(t)
                tile_feats.append(f.cpu().numpy().ravel())
            if tile_feats:
                slide_emb = np.stack(tile_feats).mean(axis=0)  # mean pooling = slide-level
                slide_features[r["file_id"]] = slide_emb
                slide_meta.append({"file_id": r["file_id"], "case_id": r["case_id"],
                                    "n_tiles_used": len(tile_feats)})
                print(f"  {len(tile_feats)} tiles embedded")
            if args.delete_after:
                slide_path.unlink(missing_ok=True)
                print(f"  deleted {slide_path}")
        except Exception as e:
            print(f"  FAIL: {e}")

    if slide_features:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        keys = list(slide_features.keys())
        feats = np.stack([slide_features[k] for k in keys])
        np.savez_compressed(args.out, features=feats, file_ids=keys, model=args.model)
        pd.DataFrame(slide_meta).to_csv(Path(args.out).with_suffix(".meta.tsv"),
                                         sep="\t", index=False)
        print(f"\n→ {args.out}: {feats.shape}")
        print(f"→ {Path(args.out).with_suffix('.meta.tsv')}")


if __name__ == "__main__":
    main()
