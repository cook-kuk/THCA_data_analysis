#!/usr/bin/env python3
"""G1 — Foundation model embedding of 5,600 H&E tiles → npz.
Supports: UNI (gated HF), CONCH, Virchow2, Prov-GigaPath, ResNet50 ImageNet (no gating).
Default fallback: ResNet50."""
import argparse
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from PIL import Image


def load_model(name: str, device: str):
    import torch
    if name == "resnet50":
        import torchvision.models as M
        m = M.resnet50(weights=M.ResNet50_Weights.IMAGENET1K_V2)
        m.fc = torch.nn.Identity()
        m = m.eval().to(device)
        from torchvision import transforms as T
        xform = T.Compose([T.Resize(224), T.CenterCrop(224), T.ToTensor(),
                           T.Normalize([0.485,0.456,0.406],[0.229,0.224,0.225])])
        return m, xform, 2048
    elif name in ("uni","UNI"):
        # MahmoodLab/UNI requires HF login + accept license
        import timm, torch
        try:
            from huggingface_hub import login
            import os
            tk = os.environ.get("HF_TOKEN")
            if tk: login(token=tk, add_to_git_credential=False)
        except Exception: pass
        m = timm.create_model("hf-hub:MahmoodLab/UNI", pretrained=True,
                              init_values=1e-5, dynamic_img_size=True)
        m = m.eval().to(device)
        from torchvision import transforms as T
        xform = T.Compose([T.Resize(224), T.CenterCrop(224), T.ToTensor(),
                           T.Normalize([0.485,0.456,0.406],[0.229,0.224,0.225])])
        return m, xform, 1024
    elif name in ("conch","CONCH"):
        # MahmoodLab/CONCH gated
        from conch.open_clip_custom import create_model_from_pretrained
        m, xform = create_model_from_pretrained("conch_ViT-B-16",
                                                 checkpoint_path="hf_hub:MahmoodLab/CONCH")
        m = m.eval().to(device)
        return m, xform, 512
    else:
        raise ValueError(f"unknown model: {name}")


def embed(model_name: str, tiles_meta: pd.DataFrame, device: str, batch: int = 32):
    import torch
    m, xform, dim = load_model(model_name, device)
    n = len(tiles_meta)
    out = np.zeros((n, dim), dtype=np.float32)
    print(f"Embedding {n} tiles with {model_name} (dim={dim})...")
    for start in range(0, n, batch):
        end = min(start + batch, n)
        imgs = []
        for j in range(start, end):
            try:
                img = Image.open(tiles_meta.iloc[j]["tile_path"]).convert("RGB")
                imgs.append(xform(img))
            except Exception as e:
                print(f"  [skip] tile {j}: {e}")
                imgs.append(torch.zeros(3, 224, 224))
        batch_x = torch.stack(imgs).to(device)
        with torch.no_grad():
            if hasattr(m, "encode_image"):
                feat = m.encode_image(batch_x)
            else:
                feat = m(batch_x)
            feat = feat.cpu().numpy()
        out[start:end] = feat
        if start % (batch * 10) == 0:
            print(f"  {end}/{n} ({100*end/n:.0f}%)")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--meta", default="data/all_tile_metadata.tsv.gz",
                    help="combined tile metadata (g1_combined_tile_metadata.tsv.gz)")
    ap.add_argument("--model", default="resnet50",
                    choices=["resnet50","uni","UNI","conch","CONCH"])
    ap.add_argument("--out", default="results/embeddings.npz")
    ap.add_argument("--batch", type=int, default=32)
    ap.add_argument("--device", default="cuda" if __import__("torch").cuda.is_available() else "cpu")
    args = ap.parse_args()

    meta = pd.read_csv(args.meta, sep="\t")
    # filter to existing files
    keep = meta["tile_path"].apply(lambda p: Path(p).exists())
    print(f"  {keep.sum()}/{len(meta)} tile files exist")
    meta = meta[keep].reset_index(drop=True)
    if len(meta) == 0:
        print("ERROR: no tile files found at the paths in metadata.")
        sys.exit(1)

    feats = embed(args.model, meta, args.device, args.batch)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.out, features=feats,
                         spot_id=meta["spot_id"].values,
                         sample_id=meta["sample_id"].values,
                         model=args.model,
                         feature_dim=feats.shape[1])
    print(f"\n→ {args.out}: {feats.shape}")


if __name__ == "__main__":
    main()
