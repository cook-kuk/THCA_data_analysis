#!/bin/bash
# Resume Pod DM after WSI download: tile, embed, classify, optional HoVer-NeXt POC, pack.
set -euo pipefail

WORK=/workspace/wsi_pathology_dm
mkdir -p "$WORK"/{tiles,embeddings,segm,geojson,results,logs,hovernext_out}

echo "=== [Phase 3] WSI tiling + tissue mask ==="
echo "SVS count: $(find "$WORK/wsi" -name '*.svs' -size +100M | wc -l)"

apt-get update -qq 2>/dev/null || true
apt-get install -y -qq openslide-tools libopenslide-dev libvips-dev python3-dev gcc g++ 2>&1 | tail -3 || true
pip install -q openslide-python==1.3.1 timm transformers huggingface_hub \
    pandas tqdm scikit-learn opencv-python-headless geojson shapely \
    scikit-image pillow torchvision 2>&1 | tail -5

python3 << 'PY'
import glob
import os

import cv2
import numpy as np
import openslide
from tqdm import tqdm

WORK = "/workspace/wsi_pathology_dm"
PATCH_SIZE = 256
TARGET_LEVEL = 1
MAX_TILES = 200

svs_paths = sorted(glob.glob(f"{WORK}/wsi/*/*.svs"))
for svs in tqdm(svs_paths, desc="tiling"):
    fid = os.path.basename(os.path.dirname(svs))
    out_dir = f"{WORK}/tiles/{fid}"
    os.makedirs(out_dir, exist_ok=True)
    if len(glob.glob(f"{out_dir}/*.png")) >= 50:
        continue
    try:
        slide = openslide.OpenSlide(svs)
    except Exception as e:
        print(f"open failed {svs}: {e}", flush=True)
        continue

    level = min(TARGET_LEVEL, slide.level_count - 1)
    width, height = slide.level_dimensions[level]
    n_x = max(1, width // PATCH_SIZE)
    n_y = max(1, height // PATCH_SIZE)
    stride_x = max(1, n_x // 30)
    stride_y = max(1, n_y // 30)
    downsample = float(slide.level_downsamples[level])
    saved = 0

    for ix in range(0, n_x, stride_x):
        for iy in range(0, n_y, stride_y):
            x = ix * PATCH_SIZE
            y = iy * PATCH_SIZE
            x0 = int(x * downsample)
            y0 = int(y * downsample)
            img = slide.read_region((x0, y0), level, (PATCH_SIZE, PATCH_SIZE)).convert("RGB")
            arr = np.array(img)
            if (arr > 200).all(axis=2).mean() > 0.75:
                continue
            hsv = cv2.cvtColor(arr, cv2.COLOR_RGB2HSV)
            if (hsv[:, :, 1] > 20).mean() < 0.08:
                continue
            img.save(f"{out_dir}/tile_{ix:04d}_{iy:04d}.png")
            saved += 1
            if saved >= MAX_TILES:
                break
        if saved >= MAX_TILES:
            break
print("tiling done", flush=True)
PY

echo "tile count: $(find "$WORK/tiles" -name '*.png' | wc -l)"

echo "=== [Phase 4] DINOv2/UNI foundation embedding ==="
python3 << 'PY'
import glob
import os
import pickle

import numpy as np
import timm
import torch
from PIL import Image
from torchvision import transforms
from tqdm import tqdm

WORK = "/workspace/wsi_pathology_dm"
device = "cuda" if torch.cuda.is_available() else "cpu"

try:
    if os.environ.get("HF_TOKEN"):
        from huggingface_hub import login
        login(token=os.environ["HF_TOKEN"], add_to_git_credential=False)
    model = timm.create_model(
        "hf-hub:MahmoodLab/uni",
        pretrained=True,
        init_values=1e-5,
        dynamic_img_size=True,
        num_classes=0,
    )
    backbone = "UNI"
    input_size = 224
except Exception as e:
    print(f"UNI unavailable ({type(e).__name__}); using DINOv2 ViT-L", flush=True)
    try:
        model = timm.create_model(
            "vit_large_patch14_dinov2.lvd142m",
            pretrained=True,
            dynamic_img_size=True,
            img_size=518,
            num_classes=0,
        )
        backbone = "DINOv2-ViT-L-518"
        input_size = 518
    except Exception as ee:
        print(f"DINOv2 ViT-L unavailable ({type(ee).__name__}); using ViT-B", flush=True)
        model = timm.create_model("vit_base_patch16_224", pretrained=True, num_classes=0)
        backbone = "ViT-B-imagenet"
        input_size = 224

model.eval().to(device)
print(f"backbone: {backbone}; input={input_size}; device={device}", flush=True)

xform = transforms.Compose([
    transforms.Resize(input_size),
    transforms.CenterCrop(input_size),
    transforms.ToTensor(),
    transforms.Normalize([0.485, 0.456, 0.406], [0.229, 0.224, 0.225]),
])

slide_emb = {}
batch_size = 8 if input_size > 300 else 32
for slide_dir in tqdm(sorted(os.listdir(f"{WORK}/tiles")), desc="embedding"):
    out_npy = f"{WORK}/embeddings/{slide_dir}_tiles.npy"
    if os.path.exists(out_npy):
        feats = np.load(out_npy)
        slide_emb[slide_dir] = feats.mean(axis=0)
        continue

    tile_paths = sorted(glob.glob(f"{WORK}/tiles/{slide_dir}/*.png"))
    if not tile_paths:
        continue
    feats = []
    with torch.no_grad():
        for i in range(0, len(tile_paths), batch_size):
            batch = torch.stack([
                xform(Image.open(tp).convert("RGB")) for tp in tile_paths[i:i + batch_size]
            ]).to(device)
            out = model(batch)
            if isinstance(out, dict):
                out = out.get("x_norm_clstoken") or out.get("features")
            feats.append(out.detach().cpu().numpy())
    feats = np.concatenate(feats, axis=0)
    np.save(out_npy, feats)
    slide_emb[slide_dir] = feats.mean(axis=0)

with open(f"{WORK}/embeddings/slide_embeddings.pkl", "wb") as f:
    pickle.dump({"backbone": backbone, "embeddings": slide_emb}, f)
print(f"slides embedded: {len(slide_emb)}", flush=True)
PY

echo "=== [Phase 5] split-level DM1 vs DM2 LOSO classifier ==="
python3 << 'PY'
import json
import pickle

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegressionCV
from sklearn.metrics import accuracy_score, roc_auc_score
from sklearn.model_selection import LeaveOneGroupOut

WORK = "/workspace/wsi_pathology_dm"
manifest = pd.read_csv(f"{WORK}/manifests/tcga_thca_wsi_dm_balanced.tsv", sep="\t")
resolved = pd.read_csv(f"{WORK}/wsi_resolved_manifest.tsv", sep="\t")
m = resolved.merge(manifest, on="case_id", how="left")

with open(f"{WORK}/embeddings/slide_embeddings.pkl", "rb") as f:
    pkg = pickle.load(f)
slide_emb = pkg["embeddings"]

X, y, groups, file_ids = [], [], [], []
for _, r in m.iterrows():
    fid = r.file_id
    if fid not in slide_emb or pd.isna(r.dm) or r.dm not in {"DM1", "DM2"}:
        continue
    X.append(slide_emb[fid])
    y.append(0 if r.dm == "DM1" else 1)
    groups.append(r.case_id)
    file_ids.append(fid)

X = np.asarray(X)
y = np.asarray(y)
groups = np.asarray(groups)
print(f"matrix: X={X.shape}, y_DM1={int((y == 0).sum())}, y_DM2={int((y == 1).sum())}", flush=True)

preds = np.full(len(y), np.nan)
auc = None
acc = None
if len(y) >= 6 and len(np.unique(y)) == 2 and len(np.unique(groups)) >= 4:
    for tr, te in LeaveOneGroupOut().split(X, y, groups):
        if len(np.unique(y[tr])) < 2:
            continue
        cv = min(3, np.bincount(y[tr]).min())
        if cv < 2:
            mdl = LogisticRegressionCV(Cs=[0.1, 1, 10], cv=2, max_iter=1000).fit(X[tr], y[tr])
        else:
            mdl = LogisticRegressionCV(Cs=[0.1, 1, 10], cv=int(cv), max_iter=1000).fit(X[tr], y[tr])
        preds[te] = mdl.predict_proba(X[te])[:, 1]
    ok = ~np.isnan(preds)
    if ok.sum() and len(np.unique(y[ok])) == 2:
        auc = float(roc_auc_score(y[ok], preds[ok]))
        acc = float(accuracy_score(y[ok], (preds[ok] > 0.5).astype(int)))
        print(f"DM1 vs DM2 split LOSO AUC={auc:.3f}, ACC={acc:.3f}", flush=True)
else:
    print("Not enough class/group diversity for split-level LOSO", flush=True)

pd.DataFrame({
    "file_id": file_ids,
    "case_id": groups,
    "dm_true": np.where(y == 0, "DM1", "DM2"),
    "dm2_pred_prob": preds,
}).to_csv(f"{WORK}/results/dm1_vs_dm2_loso_predictions.tsv", sep="\t", index=False)

json.dump({
    "backbone": pkg["backbone"],
    "auc": auc,
    "acc": acc,
    "n_DM1": int((y == 0).sum()),
    "n_DM2": int((y == 1).sum()),
    "note": "Split-level only. Merge both pod embeddings locally for final 50-WSI LOSO.",
}, open(f"{WORK}/results/dm1_vs_dm2_loso_metrics.json", "w"), indent=2)
print("metrics saved", flush=True)
PY

echo "=== [Phase 6] HoVer-NeXt POC, non-blocking ==="
set +e
cd "$WORK"
if [ ! -d hover_next_inference ]; then
    git clone --depth 1 https://github.com/digitalpathologybern/hover_next_inference 2>&1 | tail -3
fi
cd "$WORK/hover_next_inference" || exit 0
pip install -q --no-deps "segmentation-models-pytorch>=0.3.0" "timm>=0.9.6" "torch_geometric>=2.5" 2>&1 | tail -3
pip install -q pretrainedmodels efficientnet-pytorch fasteners asciitree "albumentations>=1.3" \
    geojson shapely "zarr>=2.16,<3" "numcodecs>=0.12,<0.14" tifffile h5py imagecodecs \
    networkx scikit-image staintools "spams-bin>=2.6" mahotas 2>&1 | tail -3
if [ ! -d lizard_convnextv2_large ]; then
    curl -sL "https://zenodo.org/records/10635618/files/lizard_convnextv2_large.zip?download=1" -o lcl.zip
    unzip -oq lcl.zip && rm -f lcl.zip
fi
ls "$WORK"/wsi/*/*.svs > "$WORK/wsi_paths.txt" 2>/dev/null
head -5 "$WORK/wsi_paths.txt" > "$WORK/wsi_paths_5.txt"
N_CORES=$(nproc)
python3 main.py \
    --input "$WORK/wsi_paths_5.txt" \
    --output_root "$WORK/hovernext_out/" \
    --cp lizard_convnextv2_large \
    --tta 2 \
    --inf_workers "$N_CORES" \
    --pp_tiling 8 \
    --pp_workers "$((N_CORES - 1))" 2>&1 | tail -60
set -e

echo "=== [Phase 8] Pack artifacts ==="
cd "$WORK"
du -sh wsi tiles embeddings results hovernext_out 2>/dev/null || true
tar czf "$WORK/dm_wsi_artifacts.tgz" \
    -C "$WORK" embeddings results wsi_resolved_manifest.tsv hovernext_out 2>/dev/null || \
tar czf "$WORK/dm_wsi_artifacts.tgz" \
    -C "$WORK" embeddings results wsi_resolved_manifest.tsv
ls -lh "$WORK/dm_wsi_artifacts.tgz"
echo "=== ALL DONE ==="
