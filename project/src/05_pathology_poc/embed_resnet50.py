#!/usr/bin/env python3
"""Embed H&E tiles with frozen ResNet50 (ImageNet).

Output: project/results/03_pathology_poc/embeddings_resnet50_{tile_size}.npz
  arrays:
    spot_ids   : (N,) str
    sample_ids : (N,) str
    stages     : (N,) str
    embeddings : (N, 2048) float32
"""
from __future__ import annotations
import argparse, logging, os, time
from pathlib import Path
import numpy as np
import pandas as pd
import torch
import torchvision.models as M
import torchvision.transforms as T
from PIL import Image
from torch.utils.data import Dataset, DataLoader

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("embed")
Image.MAX_IMAGE_PIXELS = None

ROOT = Path(os.environ.get("THCA_ROOT", "/home/seungho/personal/THCA_data_analysis"))
META = ROOT / "project/results/03_pathology_poc/tile_metadata_resid.tsv.gz"


class TileDS(Dataset):
    def __init__(self, paths: list[str], tfm):
        self.paths = paths
        self.tfm = tfm
    def __len__(self): return len(self.paths)
    def __getitem__(self, i):
        img = Image.open(self.paths[i]).convert("RGB")
        return self.tfm(img)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tile-size", type=int, default=224, choices=[224, 448, 672])
    ap.add_argument("--batch", type=int, default=64)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--device", default="auto")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    if args.device == "auto":
        dev = "cuda" if torch.cuda.is_available() else "cpu"
    else:
        dev = args.device
    log.info("torch=%s cuda_avail=%s device=%s", torch.__version__,
             torch.cuda.is_available(), dev)
    if dev == "cuda":
        log.info("device_name=%s", torch.cuda.get_device_name(0))

    meta = pd.read_csv(META, sep="\t")
    meta = meta[meta.tile_size == args.tile_size].reset_index(drop=True)
    if len(meta) == 0:
        log.error("no tiles with tile_size=%d in %s", args.tile_size, META)
        raise SystemExit(2)
    paths = meta.tile_path.tolist()
    log.info("tiles to embed: %d (size=%d)", len(paths), args.tile_size)

    # ImageNet normalization, no aug
    tfm = T.Compose([
        T.Resize(224),  # foundation models will diff; ResNet50 needs 224
        T.CenterCrop(224),
        T.ToTensor(),
        T.Normalize([0.485, 0.456, 0.406], [0.229, 0.224, 0.225]),
    ])

    backbone = M.resnet50(weights=M.ResNet50_Weights.IMAGENET1K_V2)
    backbone.fc = torch.nn.Identity()
    backbone.eval().to(dev)

    ds = TileDS(paths, tfm)
    dl = DataLoader(ds, batch_size=args.batch, shuffle=False,
                    num_workers=args.workers, pin_memory=(dev == "cuda"))
    feats = np.zeros((len(paths), 2048), dtype=np.float32)
    t0 = time.time()
    n = 0
    with torch.inference_mode():
        for x in dl:
            x = x.to(dev, non_blocking=True)
            with torch.amp.autocast(device_type="cuda", enabled=(dev == "cuda")):
                f = backbone(x)
            f = f.float().cpu().numpy()
            feats[n:n + len(f)] = f
            n += len(f)
            if n % (args.batch * 10) == 0 or n == len(paths):
                log.info("%d/%d embedded (%.1f tile/s)",
                         n, len(paths), n / max(time.time() - t0, 1e-6))

    out = Path(args.out or
               ROOT / f"project/results/03_pathology_poc/embeddings_resnet50_{args.tile_size}.npz")
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        out,
        spot_ids=meta.spot_id.values.astype(str),
        sample_ids=meta.sample_id.values.astype(str),
        stages=meta.stage.values.astype(str),
        embeddings=feats,
    )
    log.info("wrote %s (%.1f MB)  elapsed=%.1fs",
             out, out.stat().st_size / 1e6, time.time() - t0)


if __name__ == "__main__":
    main()
