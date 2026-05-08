#!/bin/bash
# Pod DM — TCGA-THCA WSI DM1 vs DM2 dark-matter classifier (FULL pathology run)
# Run inside RunPod A40 48GB pod after SSH connect.
#
# Strategy:
#   - 50-case DM-balanced subset (DM1 n=25, DM2 n=25; manifest pre-built locally)
#   - Foundation model embedding (UNI > CTransPath > DINOv2 fallback)
#   - Tumor/stroma segmentation (lightweight HSV/Otsu fallback if no GrandQC)
#   - Tile-level + slide-level DM1/DM2 classifier (LOSO CV)
#   - Hovernext-format GeoJSON cell detection (REQUIRED for SPARK Analytical pipeline)
#
# Required upload before launch (rsync from local):
#   project/data/manifests/tcga_thca_wsi_dm_balanced.tsv
#   project/data/manifests/tcga_thca_wsi_gdc_filter.json
#
# Expected runtime: 25-35 hr   Expected cost: $10-15 (A40 48GB @ $0.39/hr)
# Output sync target: project/data/processed/TCGA-THCA-WSI-DM/
set -euo pipefail

WORK=/workspace/wsi_pathology_dm
mkdir -p "$WORK"/{wsi,tiles,embeddings,segm,geojson,results,logs}
cd "$WORK"

UPLOAD_MANIFEST="$WORK/manifests/tcga_thca_wsi_dm_balanced.tsv"
UPLOAD_FILTER="$WORK/manifests/tcga_thca_wsi_gdc_filter.json"

# ------------------------------------------------------------------
# Phase 1 — System deps + foundation model environment
# ------------------------------------------------------------------
echo "=== [Phase 1] system deps ==="
DEBIAN_FRONTEND=noninteractive apt-get update -qq
DEBIAN_FRONTEND=noninteractive apt-get install -y -qq \
    curl wget git unzip rsync \
    openslide-tools libopenslide-dev \
    python3-pip python3-dev libvips-dev 2>&1 | tail -3

pip install -q --upgrade pip
pip install -q \
    openslide-python==1.3.1 \
    pyvips==2.2.3 \
    opencv-python-headless==4.8.1.78 \
    timm==1.0.11 \
    transformers==4.42.0 \
    huggingface_hub==0.24.6 \
    pandas==2.1.4 tqdm scikit-learn==1.3.2 \
    geojson==3.1.0 shapely==2.0.3 \
    matplotlib seaborn

python3 -c "import openslide, timm, torch; print('torch', torch.__version__, 'cuda', torch.cuda.is_available(), torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'NO_GPU')"

# ------------------------------------------------------------------
# Phase 2 — TCGA-THCA WSI download via DM-balanced GDC filter
# ------------------------------------------------------------------
echo "=== [Phase 2] WSI download (50 DM-balanced cases) ==="
cd "$WORK/wsi"

python3 << 'EOF'
import json, requests, pandas as pd, sys
flt = json.load(open("/workspace/wsi_pathology_dm/manifests/tcga_thca_wsi_gdc_filter.json"))
params = {
    "filters": json.dumps(flt),
    "fields": "file_id,file_name,cases.submitter_id,file_size,experimental_strategy",
    "format": "JSON",
    "size": "200",  # GDC may return >50 if multiple slides per case
}
r = requests.get("https://api.gdc.cancer.gov/files", params=params, timeout=120)
hits = r.json()["data"]["hits"]
rows = []
for h in hits:
    case = h["cases"][0]["submitter_id"] if h.get("cases") else "?"
    rows.append({
        "file_id": h["file_id"],
        "file_name": h["file_name"],
        "case_id": case,
        "size_GB": round(h["file_size"]/1e9, 2),
    })
df = pd.DataFrame(rows)
# pick 1 diagnostic slide per case (smallest if multiple)
df = df.sort_values(["case_id","size_GB"]).groupby("case_id").head(1)
df.to_csv("/workspace/wsi_pathology_dm/wsi_resolved_manifest.tsv", sep="\t", index=False)
print(f"resolved {len(df)} WSIs across {df.case_id.nunique()} cases, total {df.size_GB.sum():.1f} GB")
EOF

# Download with gdc-client
wget -q https://gdc.cancer.gov/files/public/file/gdc-client_2.4.0_Ubuntu_x64.zip -O /tmp/gdc.zip
unzip -oq /tmp/gdc.zip -d /usr/local/bin/
chmod +x /usr/local/bin/gdc-client

tail -n +2 "$WORK/wsi_resolved_manifest.tsv" | awk -F'\t' '{print $1}' > "$WORK/manifest_ids.txt"
gdc-client download -m "$WORK/manifest_ids.txt" --no-related-files 2>&1 | tail -10

ls "$WORK/wsi"/*/*.svs 2>/dev/null | wc -l

# ------------------------------------------------------------------
# Phase 3 — Tile WSI at 20x with HSV-Otsu tissue mask
# ------------------------------------------------------------------
echo "=== [Phase 3] WSI tiling + tissue mask ==="
mkdir -p "$WORK/tiles" "$WORK/segm"

python3 << 'EOF'
import openslide, os, glob, numpy as np, cv2
from PIL import Image
from tqdm import tqdm

WORK = "/workspace/wsi_pathology_dm"
PATCH_SIZE = 256
TARGET_LEVEL = 1   # roughly 20x in TCGA SVS (level 0 = 40x)

for svs in tqdm(sorted(glob.glob(f"{WORK}/wsi/*/*.svs"))):
    fid = os.path.basename(os.path.dirname(svs))
    out_dir = f"{WORK}/tiles/{fid}"
    seg_dir = f"{WORK}/segm/{fid}"
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(seg_dir, exist_ok=True)
    if len(os.listdir(out_dir)) > 50:
        continue

    slide = openslide.OpenSlide(svs)
    L = min(TARGET_LEVEL, slide.level_count - 1)
    w, h = slide.level_dimensions[L]
    # thumbnail tissue mask via HSV+Otsu at low-res
    thumb = np.array(slide.get_thumbnail((w//8, h//8)).convert("RGB"))
    hsv = cv2.cvtColor(thumb, cv2.COLOR_RGB2HSV)
    _, sat_mask = cv2.threshold(hsv[..., 1], 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    cv2.imwrite(f"{seg_dir}/tissue_mask_thumb.png", sat_mask)

    n_x = w // PATCH_SIZE
    n_y = h // PATCH_SIZE
    # Up to 200 tiles per slide for foundation-model embed
    stride_x = max(1, n_x // 15)
    stride_y = max(1, n_y // 15)
    saved = 0
    for ix in range(0, n_x, stride_x):
        for iy in range(0, n_y, stride_y):
            x = ix * PATCH_SIZE
            y = iy * PATCH_SIZE
            img = slide.read_region((x, y), L, (PATCH_SIZE, PATCH_SIZE)).convert("RGB")
            arr = np.array(img)
            # background filter
            if (arr > 200).all(axis=2).mean() > 0.75:
                continue
            img.save(f"{out_dir}/tile_{ix:04d}_{iy:04d}.png")
            saved += 1
            if saved >= 200:
                break
        if saved >= 200:
            break
print("tiling done")
EOF

# ------------------------------------------------------------------
# Phase 4 — Foundation model embedding (UNI > CTransPath > DINOv2)
# ------------------------------------------------------------------
echo "=== [Phase 4] foundation model embedding ==="

python3 << 'EOF'
import os, glob, numpy as np, torch, timm
from PIL import Image
from torchvision import transforms
from tqdm import tqdm

WORK = "/workspace/wsi_pathology_dm"
device = "cuda" if torch.cuda.is_available() else "cpu"

# Try UNI (gated, requires HF token in env HF_TOKEN). Fallback chain.
model = None
backbone = None
try:
    from huggingface_hub import login
    if os.environ.get("HF_TOKEN"):
        login(token=os.environ["HF_TOKEN"], add_to_git_credential=False)
    model = timm.create_model(
        "hf-hub:MahmoodLab/uni",
        pretrained=True, init_values=1e-5, dynamic_img_size=True,
    )
    backbone = "UNI"
except Exception as e:
    print(f"UNI unavailable ({e}); falling back to DINOv2 ViT-L")
    try:
        model = timm.create_model("vit_large_patch14_dinov2.lvd142m", pretrained=True)
        backbone = "DINOv2-ViT-L"
    except Exception as ee:
        print(f"DINOv2 unavailable ({ee}); using ViT-B ImageNet")
        model = timm.create_model("vit_base_patch16_224", pretrained=True)
        backbone = "ViT-B-imagenet"

model.eval().to(device)
print(f"backbone: {backbone}")

xform = transforms.Compose([
    transforms.Resize(224),
    transforms.CenterCrop(224),
    transforms.ToTensor(),
    transforms.Normalize([0.485,0.456,0.406],[0.229,0.224,0.225]),
])

records = []
slide_emb = {}
for slide_dir in tqdm(sorted(os.listdir(f"{WORK}/tiles"))):
    tps = sorted(glob.glob(f"{WORK}/tiles/{slide_dir}/*.png"))
    if not tps:
        continue
    feats = []
    with torch.no_grad():
        for tp in tps:
            x = xform(Image.open(tp).convert("RGB")).unsqueeze(0).to(device)
            f = model(x).cpu().numpy().flatten() if backbone == "UNI" else \
                model.forward_features(x).mean(dim=1).cpu().numpy().flatten()
            feats.append(f)
            records.append({"slide": slide_dir, "tile": os.path.basename(tp)})
    feats = np.array(feats)
    slide_emb[slide_dir] = feats.mean(axis=0)
    np.save(f"{WORK}/embeddings/{slide_dir}_tiles.npy", feats)

import pickle, json
with open(f"{WORK}/embeddings/slide_embeddings.pkl", "wb") as f:
    pickle.dump({"backbone": backbone, "embeddings": slide_emb}, f)
print(f"slides embedded: {len(slide_emb)}")
EOF

# ------------------------------------------------------------------
# Phase 5 — DM1 vs DM2 classifier (LOSO CV)
# ------------------------------------------------------------------
echo "=== [Phase 5] DM1 vs DM2 classifier ==="
python3 << 'EOF'
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
    if fid not in slide_emb or pd.isna(r.dm) or r.dm not in {"DM1","DM2"}:
        continue
    X.append(slide_emb[fid])
    y.append(0 if r.dm == "DM1" else 1)
    groups.append(r.case_id)

X = np.array(X); y = np.array(y); groups = np.array(groups)
print(f"matrix: X={X.shape}, y={y.shape}, n_DM1={int((y==0).sum())}, n_DM2={int((y==1).sum())}")

# Leave-one-CASE-out logistic regression
loso = LeaveOneGroupOut()
preds = np.zeros_like(y, dtype=float)
for tr, te in loso.split(X, y, groups):
    mdl = LogisticRegressionCV(Cs=[0.01,0.1,1,10], cv=3, max_iter=1000).fit(X[tr], y[tr])
    preds[te] = mdl.predict_proba(X[te])[:, 1]

auc = roc_auc_score(y, preds)
acc = accuracy_score(y, (preds > 0.5).astype(int))
print(f"DM1 vs DM2 LOSO AUC={auc:.3f}, ACC={acc:.3f}")

pd.DataFrame({
    "case_id": groups, "dm_true": np.where(y == 0, "DM1", "DM2"),
    "dm_pred_prob": preds,
}).to_csv(f"{WORK}/results/dm1_vs_dm2_loso_predictions.tsv", sep="\t", index=False)
json.dump({"backbone": pkg["backbone"], "auc": float(auc), "acc": float(acc),
           "n_DM1": int((y==0).sum()), "n_DM2": int((y==1).sum())},
          open(f"{WORK}/results/dm1_vs_dm2_loso_metrics.json", "w"), indent=2)
print("metrics saved")
EOF

# ------------------------------------------------------------------
# Phase 6 — HoVer-NeXt cell detection + classification (SPARK Analytical input)
# ------------------------------------------------------------------
# Reference: digitalpathologybern/hover_next_inference + Lizard-Mitosis ConvNeXtV2-Large weights
# Output: per-WSI .zarr instance map + .tsv class lookup with cell centroids,
#         convertible to Hovernext-format GeoJSON for SPARK MAIN_FUNCTION_CALLER.
echo "=== [Phase 6] HoVer-NeXt cell detection ==="
cd "$WORK"
[ ! -d hover_next_inference ] && git clone --depth 1 https://github.com/digitalpathologybern/hover_next_inference
cd hover_next_inference

# Install HoVer-NeXt deps inside this pod's torch+CUDA env (preserve existing torch wheel)
pip install -q --no-deps -r requirements.txt 2>&1 | tail -3
pip install -q pretrainedmodels efficientnet-pytorch fasteners asciitree mahotas==1.4.18 \
    segmentation-models-pytorch==0.3.4 zarr==2.16.1 numcodecs==0.12.1 \
    spams-bin staintools shapely==2.0.2 networkx==2.8.7 geojson==3.1.0 albumentations==1.3.1 \
    timm==0.9.6 tifffile h5py imagecodecs 2>&1 | tail -5

# Pretrained weights — Lizard convnextv2-large gives 7-class output (Tumor / Stroma / Lymphocyte
# / Plasma / Eosinophil / Neutrophil / Connective)
[ ! -d lizard_convnextv2_large ] && \
  curl -sL "https://zenodo.org/records/10635618/files/lizard_convnextv2_large.zip?download=1" -o lcl.zip && \
  unzip -oq lcl.zip && rm lcl.zip
ls -la lizard_convnextv2_large/ | head -6

# Run on 50 WSIs in parallel
mkdir -p "$WORK/hovernext_out"
ls "$WORK"/wsi/*/*.svs > "$WORK/wsi_paths.txt"
N_CORES=$(nproc)
python3 main.py \
    --input "$WORK/wsi_paths.txt" \
    --output_root "$WORK/hovernext_out/" \
    --cp lizard_convnextv2_large \
    --tta 4 \
    --inf_workers "$N_CORES" \
    --pp_tiling 8 \
    --pp_workers "$((N_CORES - 1))" 2>&1 | tail -20

# ------------------------------------------------------------------
# Phase 6b — Convert HoVer-NeXt output to Hovernext-format GeoJSON for SPARK
# ------------------------------------------------------------------
echo "=== [Phase 6b] HoVer-NeXt → SPARK GeoJSON ==="
python3 << 'EOF'
import os, glob, json, zarr, numpy as np, pandas as pd
WORK = "/workspace/wsi_pathology_dm"
OUT = f"{WORK}/geojson"
os.makedirs(OUT, exist_ok=True)
# HoVer-NeXt PanNuke/Lizard class indexes
CLASS_NAME = {0: "background", 1: "tumor", 2: "stroma", 3: "lymphocyte",
              4: "plasma", 5: "eosinophil", 6: "neutrophil", 7: "connective"}
for d in sorted(glob.glob(f"{WORK}/hovernext_out/*/")):
    sid = os.path.basename(d.rstrip("/"))
    inst_path = os.path.join(d, "pinst_pp.zip")
    cls_path  = os.path.join(d, "class_lookup.csv")
    if not os.path.exists(cls_path):
        continue
    df = pd.read_csv(cls_path)
    feats = []
    for _, r in df.iterrows():
        feats.append({
            "type": "Feature",
            "geometry": {"type": "Point", "coordinates": [float(r.get("x", 0)), float(r.get("y", 0))]},
            "properties": {
                "classification": {"name": CLASS_NAME.get(int(r.get("class", 0)), "other")},
                "instance_id": int(r.get("id", 0)),
            },
        })
    out_geo = f"{OUT}/{sid}.geojson"
    with open(out_geo, "w") as f:
        json.dump({"type": "FeatureCollection", "features": feats}, f)
    print(f"{sid}: {len(feats)} cells")
EOF

# ------------------------------------------------------------------
# Phase 7 — GigaTIME virtual mIF (Cell 2025; prov-gigatime/GigaTIME)
# ------------------------------------------------------------------
# Generates 23-channel virtual immunofluorescence (DAPI/CD8/CD4/CD20/CD3/CD138/
# PD-1/PD-L1/CK/Ki67/etc) from H&E. Cross-validates HoVer-NeXt cell-class calls
# and adds clinically-actionable IO markers (PD-L1, Tryptase, Caspase3, Ki67).
#
# Requires: HF_TOKEN env var with read access to prov-gigatime/GigaTIME.
#           Set on RunPod pod via:  export HF_TOKEN="hf_...."
echo "=== [Phase 7] GigaTIME virtual mIF (23-channel) ==="
if [ -z "${HF_TOKEN:-}" ]; then
  echo "HF_TOKEN unset — skipping GigaTIME (request token from prov-gigatime HF model card)"
else
  cd "$WORK"
  [ ! -d GigaTIME ] && git clone --depth 1 https://github.com/prov-gigatime/GigaTIME GigaTIME
  cd GigaTIME
  pip install -q numpy pandas scipy scikit-learn scikit-image albumentations==1.4.0 \
      easydict edict pyyaml huggingface_hub 2>&1 | tail -3

  python3 << 'EOF'
import os, glob, sys, numpy as np, torch
sys.path.insert(0, "scripts")
from huggingface_hub import snapshot_download
from PIL import Image

os.environ.setdefault("HUGGINGFACE_HUB_TOKEN", os.environ.get("HF_TOKEN", ""))
local = snapshot_download(repo_id="prov-gigatime/GigaTIME")
print(f"GigaTIME weights -> {local}")

# Load arch
import archs
device = "cuda" if torch.cuda.is_available() else "cpu"
model = archs.__dict__["gigatime"](23, 3, deep_supervision=False).to(device)
state = torch.load(os.path.join(local, "model.pth"), map_location=device)
if "state_dict" in state: state = state["state_dict"]
state = {k.replace("module.", ""): v for k, v in state.items()}
model.load_state_dict(state, strict=False)
model.eval()

WORK = "/workspace/wsi_pathology_dm"
mean = np.array([0.485, 0.456, 0.406]); std = np.array([0.229, 0.224, 0.225])
CHANNELS = ['DAPI','TRITC','Cy5','PD-1','CD14','CD4','T-bet','CD34','CD68','CD16','CD11c',
            'CD138','CD20','CD3','CD8','PD-L1','CK','Ki67','Tryptase','Actin-D',
            'Caspase3-D','PHH3-B','Transgelin']
out_root = f"{WORK}/gigatime_mIF"
os.makedirs(out_root, exist_ok=True)

# Run on first 30 tiles per slide for budget
n_done = 0
for slide in sorted(os.listdir(f"{WORK}/tiles")):
    out_dir = f"{out_root}/{slide}"
    os.makedirs(out_dir, exist_ok=True)
    paths = sorted(glob.glob(f"{WORK}/tiles/{slide}/*.png"))[:30]
    for tp in paths:
        out_npz = f"{out_dir}/{os.path.basename(tp)}.npz"
        if os.path.exists(out_npz): continue
        img = np.array(Image.open(tp).convert("RGB").resize((512, 512))) / 255.0
        x = ((img - mean) / std).transpose(2, 0, 1).astype(np.float32)
        x = torch.from_numpy(x).unsqueeze(0).to(device)
        with torch.no_grad():
            y = model(x).cpu().numpy().squeeze(0)  # (23, H, W)
        # Per-channel mean intensity (lightweight summary)
        means = y.mean(axis=(1, 2))
        np.savez_compressed(out_npz, channels=CHANNELS, mean=means, sample=y[:, ::4, ::4])
        n_done += 1
        if n_done % 50 == 0:
            print(f"  ...{n_done} GigaTIME inferences done")
print(f"GigaTIME total: {n_done} tile inferences -> {out_root}")
EOF
fi

# ------------------------------------------------------------------
# Phase 8 — Sync results back to local repo
# ------------------------------------------------------------------
echo "=== [Phase 8] result summary ==="
ls -la "$WORK/results/"
du -sh "$WORK"/{wsi,tiles,embeddings,geojson,results,gigatime_mIF} 2>/dev/null

# Pack lightweight artifacts (embeddings + classifier + geojson + gigatime channel summaries)
tar czf "$WORK/dm_wsi_artifacts.tgz" \
    -C "$WORK" \
    embeddings results geojson segm wsi_resolved_manifest.tsv \
    gigatime_mIF 2>/dev/null || tar czf "$WORK/dm_wsi_artifacts.tgz" -C "$WORK" \
        embeddings results geojson segm wsi_resolved_manifest.tsv

echo "DONE — sync $WORK/dm_wsi_artifacts.tgz back to project/data/processed/TCGA-THCA-WSI-DM/"
echo "       and $WORK/wsi/ if you want raw WSI cached locally"
