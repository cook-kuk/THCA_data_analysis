#!/bin/bash
# Pod DM Phase 4 fix — resize tiles to 518 for DINOv2 ViT-L OR use dynamic img size.
# Then run Phase 5 (LOSO), Phase 6 (HoVer-NeXt), Phase 8 (pack).
set -uo pipefail
WORK=/workspace/wsi_pathology_dm

echo "=== [Phase 4 fix] DINOv2 ViT-L embedding (518 input) ==="
python3 << 'PY'
import os, glob, numpy as np, torch, timm
from PIL import Image
from torchvision import transforms
from tqdm import tqdm

WORK = "/workspace/wsi_pathology_dm"
device = "cuda" if torch.cuda.is_available() else "cpu"
backbone = None
model = None

# Try UNI
try:
    if os.environ.get("HF_TOKEN"):
        from huggingface_hub import login
        login(token=os.environ["HF_TOKEN"], add_to_git_credential=False)
    model = timm.create_model("hf-hub:MahmoodLab/uni", pretrained=True, init_values=1e-5, dynamic_img_size=True)
    backbone = "UNI"
except Exception as e:
    print(f"UNI unavailable; using DINOv2 ViT-L with dynamic_img_size", flush=True)
    try:
        model = timm.create_model("vit_large_patch14_dinov2.lvd142m", pretrained=True,
                                  dynamic_img_size=True, img_size=518, num_classes=0)
        backbone = "DINOv2-ViT-L-518"
    except Exception:
        # try img_size=224 explicit
        model = timm.create_model("vit_large_patch14_dinov2.lvd142m", pretrained=True,
                                  num_classes=0, dynamic_img_size=True)
        backbone = "DINOv2-ViT-L-dynamic"

model.eval().to(device)
# detect required size from patch_embed
PATCH_SIZE = getattr(model.patch_embed, "img_size", (518, 518))
if isinstance(PATCH_SIZE, tuple):
    INPUT = PATCH_SIZE[0]
else:
    INPUT = PATCH_SIZE
print(f"backbone: {backbone}, input size: {INPUT}", flush=True)

xform = transforms.Compose([
    transforms.Resize(INPUT),
    transforms.CenterCrop(INPUT),
    transforms.ToTensor(),
    transforms.Normalize([0.485,0.456,0.406],[0.229,0.224,0.225]),
])

slide_emb = {}
for slide_dir in tqdm(sorted(os.listdir(f"{WORK}/tiles"))):
    out_npy = f"{WORK}/embeddings/{slide_dir}_tiles.npy"
    if os.path.exists(out_npy):
        feats = np.load(out_npy)
        slide_emb[slide_dir] = feats.mean(axis=0)
        continue
    tps = sorted(glob.glob(f"{WORK}/tiles/{slide_dir}/*.png"))
    if not tps: continue
    feats = []
    with torch.no_grad():
        BS = 8 if INPUT > 300 else 32
        for i in range(0, len(tps), BS):
            batch = torch.stack([xform(Image.open(tp).convert("RGB")) for tp in tps[i:i+BS]]).to(device)
            try:
                f = model(batch)
            except Exception as e:
                print(f"  fwd fail at {slide_dir}: {e}", flush=True)
                break
            feats.append(f.cpu().numpy())
    if not feats: continue
    feats = np.concatenate(feats, axis=0)
    slide_emb[slide_dir] = feats.mean(axis=0)
    np.save(out_npy, feats)

import pickle
with open(f"{WORK}/embeddings/slide_embeddings.pkl", "wb") as f:
    pickle.dump({"backbone": backbone, "embeddings": slide_emb}, f)
print(f"slides embedded: {len(slide_emb)}", flush=True)
PY

# ------------------------------------------------------------------
# Phase 5 — DM1 vs DM2 LOSO classifier
# ------------------------------------------------------------------
echo "=== [Phase 5] DM1 vs DM2 classifier ==="
python3 << 'PY'
import pickle, pandas as pd, numpy as np, json
from sklearn.linear_model import LogisticRegressionCV
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.metrics import roc_auc_score, accuracy_score

WORK = "/workspace/wsi_pathology_dm"
manifest = pd.read_csv(f"{WORK}/manifests/tcga_thca_wsi_dm_balanced.tsv", sep="\t")
resolved = pd.read_csv(f"{WORK}/wsi_resolved_manifest.tsv", sep="\t")
m = resolved.merge(manifest, on="case_id", how="left")

with open(f"{WORK}/embeddings/slide_embeddings.pkl", "rb") as f:
    pkg = pickle.load(f)
slide_emb = pkg["embeddings"]

X, y, groups = [], [], []
for _, r in m.iterrows():
    fid = r.file_id
    if fid not in slide_emb or pd.isna(r.dm) or r.dm not in {"DM1", "DM2"}:
        continue
    X.append(slide_emb[fid])
    y.append(0 if r.dm == "DM1" else 1)
    groups.append(r.case_id)

X = np.array(X); y = np.array(y); groups = np.array(groups)
print(f"matrix: X={X.shape}, y_DM1={int((y==0).sum())}, y_DM2={int((y==1).sum())}", flush=True)
loso = LeaveOneGroupOut()
preds = np.zeros_like(y, dtype=float)
for tr, te in loso.split(X, y, groups):
    mdl = LogisticRegressionCV(Cs=[0.01,0.1,1,10], cv=3, max_iter=1000).fit(X[tr], y[tr])
    preds[te] = mdl.predict_proba(X[te])[:, 1]

auc = roc_auc_score(y, preds)
acc = accuracy_score(y, (preds > 0.5).astype(int))
print(f"DM1 vs DM2 LOSO AUC={auc:.3f}, ACC={acc:.3f}", flush=True)

pd.DataFrame({"case_id": groups, "dm_true": np.where(y == 0, "DM1", "DM2"),
              "dm_pred_prob": preds}).to_csv(
    f"{WORK}/results/dm1_vs_dm2_loso_predictions.tsv", sep="\t", index=False)
json.dump({"backbone": pkg["backbone"], "auc": float(auc), "acc": float(acc),
           "n_DM1": int((y==0).sum()), "n_DM2": int((y==1).sum())},
          open(f"{WORK}/results/dm1_vs_dm2_loso_metrics.json", "w"), indent=2)
print("metrics saved", flush=True)
PY

# ------------------------------------------------------------------
# Phase 6 — HoVer-NeXt (skip if too heavy / continue if OK)
# ------------------------------------------------------------------
echo "=== [Phase 6] HoVer-NeXt cell detection ==="
cd "$WORK"
[ ! -d hover_next_inference ] && git clone --depth 1 https://github.com/digitalpathologybern/hover_next_inference 2>&1 | tail -2
cd hover_next_inference
pip install -q --no-deps -r requirements.txt 2>&1 | tail -2
pip install -q pretrainedmodels efficientnet-pytorch fasteners asciitree mahotas==1.4.18 \
    segmentation-models-pytorch==0.3.4 zarr==2.16.1 numcodecs==0.12.1 \
    spams-bin staintools shapely==2.0.2 networkx==2.8.7 geojson==3.1.0 albumentations==1.3.1 \
    timm==0.9.6 tifffile h5py imagecodecs 2>&1 | tail -2

if [ ! -d lizard_convnextv2_large ]; then
    curl -sL "https://zenodo.org/records/10635618/files/lizard_convnextv2_large.zip?download=1" -o lcl.zip
    unzip -oq lcl.zip && rm lcl.zip
fi

mkdir -p "$WORK/hovernext_out"
ls "$WORK"/wsi/*/*.svs > "$WORK/wsi_paths.txt" 2>/dev/null
N_SVS=$(wc -l < "$WORK/wsi_paths.txt")
N_CORES=$(nproc)
echo "running HoVer-NeXt on $N_SVS WSIs..."

python3 main.py \
    --input "$WORK/wsi_paths.txt" \
    --output_root "$WORK/hovernext_out/" \
    --cp lizard_convnextv2_large \
    --tta 4 \
    --inf_workers "$N_CORES" \
    --pp_tiling 8 \
    --pp_workers "$((N_CORES - 1))" 2>&1 | tail -30 || echo "HoVer-NeXt encountered issues, continuing"

# ------------------------------------------------------------------
# Phase 8 — Pack artifacts
# ------------------------------------------------------------------
echo "=== [Phase 8] Pack artifacts ==="
cd "$WORK"
ls -la results/ embeddings/ 2>/dev/null
du -sh wsi tiles embeddings results geojson hovernext_out 2>/dev/null

tar czf "$WORK/dm_wsi_artifacts.tgz" \
    -C "$WORK" \
    embeddings results geojson segm wsi_resolved_manifest.tsv \
    hovernext_out 2>/dev/null || \
tar czf "$WORK/dm_wsi_artifacts.tgz" \
    -C "$WORK" embeddings results wsi_resolved_manifest.tsv

ls -la "$WORK/dm_wsi_artifacts.tgz"
echo "=== ALL DONE ==="
