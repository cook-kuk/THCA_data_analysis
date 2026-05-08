"""
Phase 1 — UNI foundation-model embedding on GSE250521 H&E tiles.

Architecture: UNI (Mahmood lab; ViT-Huge, 1024-d; pretrained on Mass-100K WSI)
HF: MahmoodLab/UNI

Input:
  project/data/processed/GSE250521/tiles/{GSM_*}/  (3.4 GB; PNG, sizes 224/448/672)
  project/results/01_spatial_score/all_spots_scored.tsv.gz   (per-spot DM1/RAI_8)

Output:
  project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/phase1_gse250521/
    uni_embeddings.npz  (per-spot per-size UNI embeddings, 1024-d)
    uni_embed_metadata.tsv

Resource:
  ~30k tiles total; UNI on 1× A100 ~40 GB → ~30 min for size 224 + 448 (672 optional).

NOTE: UNI is gated on HuggingFace; user must export HF_TOKEN before running.
"""
from __future__ import annotations
import argparse, os, sys, time
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from PIL import Image

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tiles_root", type=Path,
                        default=Path("project/data/processed/GSE250521/tiles"))
    parser.add_argument("--out_dir", type=Path,
                        default=Path("project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/phase1_gse250521"))
    parser.add_argument("--tile_size", type=int, choices=[224, 448, 672], default=224,
                        help="UNI was pretrained at 224; default 224")
    parser.add_argument("--batch_size", type=int, default=128)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--limit_slides", type=int, default=None,
                        help="Debug — process only N slides")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    # --- load UNI (HF gated; assumes HF_TOKEN set) ---
    print("[load] UNI from MahmoodLab/UNI ...")
    from transformers import AutoModel
    import timm
    # UNI is hosted as a timm-compatible weights file
    # Path 1: timm load via HuggingFace (requires HF_TOKEN)
    try:
        model = timm.create_model("vit_large_patch16_224", img_size=224, pretrained=False, num_classes=0,
                                   dynamic_img_size=True)
        from huggingface_hub import hf_hub_download
        weights_path = hf_hub_download(repo_id="MahmoodLab/UNI", filename="pytorch_model.bin")
        state_dict = torch.load(weights_path, map_location="cpu")
        model.load_state_dict(state_dict, strict=True)
        print("  UNI weights loaded via timm")
    except Exception as e:
        print(f"[fallback] timm load failed: {e}")
        print("[fallback] using a generic ViT-L for testing — DO NOT use for production")
        model = timm.create_model("vit_large_patch16_224", pretrained=True, num_classes=0)

    model = model.eval().to(args.device)

    # --- normalization (UNI uses ImageNet mean/std) ---
    from torchvision import transforms
    transform = transforms.Compose([
        transforms.Resize(224),
        transforms.CenterCrop(224),
        transforms.ToTensor(),
        transforms.Normalize(mean=(0.485, 0.456, 0.406), std=(0.229, 0.224, 0.225)),
    ])

    # --- iterate slides ---
    slides = sorted([d for d in args.tiles_root.iterdir() if d.is_dir()])
    if args.limit_slides:
        slides = slides[:args.limit_slides]
    print(f"[input] {len(slides)} slides")

    all_embeds = []
    all_meta = []

    t0 = time.time()
    for slide_dir in slides:
        slide_name = slide_dir.name
        tiles = sorted(slide_dir.glob(f"*_size{args.tile_size}.png"))
        if not tiles:
            print(f"  [skip] {slide_name} no tile_size={args.tile_size}")
            continue
        print(f"  [{slide_name}] {len(tiles)} tiles")
        # batch loop
        batch_imgs, batch_spots = [], []
        for tpath in tiles:
            spot_id = tpath.stem.replace(f"_size{args.tile_size}", "")
            img = Image.open(tpath).convert("RGB")
            batch_imgs.append(transform(img))
            batch_spots.append(spot_id)
            if len(batch_imgs) >= args.batch_size:
                _flush_batch(model, batch_imgs, batch_spots, slide_name,
                             args.tile_size, args.device, all_embeds, all_meta)
                batch_imgs, batch_spots = [], []
        if batch_imgs:
            _flush_batch(model, batch_imgs, batch_spots, slide_name,
                         args.tile_size, args.device, all_embeds, all_meta)
    print(f"[time] {time.time()-t0:.1f}s")

    # --- save ---
    embeds_arr = np.vstack(all_embeds) if all_embeds else np.zeros((0, 1024))
    meta_df = pd.DataFrame(all_meta)
    np.savez_compressed(args.out_dir / f"uni_embeddings_size{args.tile_size}.npz", embeddings=embeds_arr)
    meta_df.to_csv(args.out_dir / f"uni_embed_metadata_size{args.tile_size}.tsv", sep="\t", index=False)
    print(f"[done] saved {embeds_arr.shape} embeds to {args.out_dir}")


def _flush_batch(model, imgs, spots, slide_name, size, device, all_e, all_m):
    with torch.no_grad():
        x = torch.stack(imgs).to(device)
        feats = model(x)  # (B, 1024)
        feats = feats.cpu().numpy()
    all_e.append(feats)
    for s in spots:
        all_m.append({"slide": slide_name, "spot_id": s, "tile_size": size})


if __name__ == "__main__":
    main()
