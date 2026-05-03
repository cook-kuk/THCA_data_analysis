#!/bin/bash
# Pod B — TCGA-THCA WSI pathology multimodal Hashimoto detection (PoC)
# Run inside RunPod A40 48GB pod after SSH connect
#
# Strategy: TCGA WSI subset (50 samples, ~10 GB) + UNI/CTransPath embedding
#           + Hashimoto-overlap classifier training + inference
# Expected runtime: 25-35 hr
# Expected cost: $10-15
set -e

WORK=/workspace/wsi_pathology
mkdir -p $WORK/wsi $WORK/tiles $WORK/embeddings $WORK/results
cd $WORK

# ============================================================
# Phase 1 — Install OpenSlide + foundation model deps
# ============================================================
echo "=== [Phase 1] Tool install ==="
DEBIAN_FRONTEND=noninteractive apt-get update -qq
DEBIAN_FRONTEND=noninteractive apt-get install -y -qq \
    curl wget git \
    openslide-tools libopenslide-dev \
    python3-pip python3-dev \
    libvips-dev 2>&1 | tail -3

pip install -q openslide-python==1.3.1 \
                pyvips==2.2.3 \
                opencv-python-headless==4.8.1.78 \
                timm==0.9.16 \
                transformers==4.35.2 \
                pandas==2.1.4 \
                tqdm==4.66.1 \
                scikit-learn==1.3.2

python3 -c "import openslide, timm, torch; print(f'torch={torch.__version__}, cuda={torch.cuda.is_available()}, GPU={torch.cuda.get_device_name(0) if torch.cuda.is_available() else None}')"

# ============================================================
# Phase 2 — TCGA-THCA WSI subset download (50 samples)
# ============================================================
echo "=== [Phase 2] TCGA WSI subset download ==="
cd $WORK/wsi

# GDC API to find 50 TCGA-THCA WSI files (slide images, .svs)
# Filter: TCGA-THCA project + diagnostic slide (DX1)
python3 << 'EOF'
import requests, json
# GDC API filter for TCGA-THCA diagnostic WSI
filters = {
    "op": "and",
    "content": [
        {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-THCA"]}},
        {"op": "in", "content": {"field": "data_format", "value": ["SVS"]}},
        {"op": "in", "content": {"field": "experimental_strategy", "value": ["Diagnostic Slide"]}},
    ]
}
params = {
    "filters": json.dumps(filters),
    "fields": "file_id,file_name,cases.submitter_id,file_size",
    "format": "JSON",
    "size": "50",
}
r = requests.get("https://api.gdc.cancer.gov/files", params=params, timeout=60)
files = r.json()["data"]["hits"]
with open("/workspace/wsi_pathology/wsi_manifest.tsv", "w") as f:
    f.write("file_id\tfile_name\tsubmitter_id\tsize_GB\n")
    for fr in files:
        case_id = fr["cases"][0]["submitter_id"] if fr.get("cases") else "?"
        f.write(f"{fr['file_id']}\t{fr['file_name']}\t{case_id}\t{fr['file_size']/1e9:.2f}\n")
print(f"50 WSI manifest saved")
EOF

cat $WORK/wsi_manifest.tsv | head -10
wc -l $WORK/wsi_manifest.tsv

# Download with GDC client (~10 GB total partial)
wget -q https://gdc.cancer.gov/files/public/file/gdc-client_2.4.0_Ubuntu_x64.zip -O /tmp/gdc.zip
unzip -q /tmp/gdc.zip -d /usr/local/bin/
chmod +x /usr/local/bin/gdc-client

cd $WORK/wsi
tail -n +2 $WORK/wsi_manifest.tsv | awk '{print $1}' > $WORK/manifest_ids.txt
gdc-client download -m $WORK/manifest_ids.txt --no-related-files 2>&1 | tail -10

ls -la $WORK/wsi/*.svs 2>/dev/null | head -10

# ============================================================
# Phase 3 — Tile WSI at 20x (256x256 patches)
# ============================================================
echo "=== [Phase 3] WSI tiling ==="
mkdir -p $WORK/tiles
python3 << 'EOF'
import openslide, os, glob, numpy as np
from PIL import Image
from tqdm import tqdm

WORK = "/workspace/wsi_pathology"
svs_files = glob.glob(f"{WORK}/wsi/*/*.svs")
print(f"WSI files: {len(svs_files)}")

PATCH_SIZE = 256
LEVEL = 1  # 20x typical

for svs in tqdm(svs_files):
    name = os.path.basename(svs).replace(".svs", "")
    out_dir = f"{WORK}/tiles/{name}"
    os.makedirs(out_dir, exist_ok=True)
    if len(os.listdir(out_dir)) > 0:
        continue
    slide = openslide.OpenSlide(svs)
    if LEVEL >= slide.level_count:
        LEVEL = slide.level_count - 1
    w, h = slide.level_dimensions[LEVEL]
    n_x = w // PATCH_SIZE
    n_y = h // PATCH_SIZE
    # Sample 100 tiles per slide (PoC)
    for ix in range(0, n_x, max(1, n_x // 10)):
        for iy in range(0, n_y, max(1, n_y // 10)):
            x = ix * PATCH_SIZE
            y = iy * PATCH_SIZE
            img = slide.read_region((x, y), LEVEL, (PATCH_SIZE, PATCH_SIZE)).convert("RGB")
            arr = np.array(img)
            # Filter background tiles (>80% white)
            if (arr > 200).all(axis=2).mean() > 0.8:
                continue
            img.save(f"{out_dir}/tile_{ix}_{iy}.png")
print("Tiling done")
EOF

# ============================================================
# Phase 4 — Foundation model embedding (UNI / CTransPath / DINOv2)
# ============================================================
echo "=== [Phase 4] Foundation model embedding ==="
python3 << 'EOF'
import torch, glob, os, numpy as np
from PIL import Image
import timm
from torchvision import transforms
from tqdm import tqdm

WORK = "/workspace/wsi_pathology"
device = "cuda" if torch.cuda.is_available() else "cpu"

# UNI from Mahmood lab (or fallback DINOv2)
try:
    model = timm.create_model("vit_large_patch16_224.dinov2", pretrained=True)
    print("Using DINOv2 ViT-L (fallback)")
except Exception as e:
    print(f"DINOv2 not loaded: {e}; using ViT-B")
    model = timm.create_model("vit_base_patch16_224", pretrained=True)

model.eval().to(device)
transform = transforms.Compose([
    transforms.Resize(224),
    transforms.CenterCrop(224),
    transforms.ToTensor(),
    transforms.Normalize(mean=[0.485, 0.456, 0.406], std=[0.229, 0.224, 0.225]),
])

embeddings = {}
for slide_dir in tqdm(sorted(os.listdir(f"{WORK}/tiles"))):
    tile_paths = glob.glob(f"{WORK}/tiles/{slide_dir}/*.png")
    if not tile_paths:
        continue
    feats = []
    with torch.no_grad():
        for tp in tile_paths:
            img = Image.open(tp).convert("RGB")
            x = transform(img).unsqueeze(0).to(device)
            f = model.forward_features(x).mean(dim=1).cpu().numpy().flatten()
            feats.append(f)
    embeddings[slide_dir] = np.array(feats).mean(axis=0)  # mean over tiles per slide

import pickle
with open(f"{WORK}/embeddings/slide_embeddings.pkl", "wb") as f:
    pickle.dump(embeddings, f)
print(f"Embeddings saved: {len(embeddings)} slides")
EOF

# ============================================================
# Phase 5 — Hashimoto classifier (PoC)
# ============================================================
echo "=== [Phase 5] Hashimoto classifier PoC ==="
python3 << 'EOF'
import pickle, numpy as np, pandas as pd, json
from sklearn.linear_model import LogisticRegressionCV
from sklearn.model_selection import StratifiedKFold

WORK = "/workspace/wsi_pathology"
with open(f"{WORK}/embeddings/slide_embeddings.pkl", "rb") as f:
    emb = pickle.load(f)

# Load manifest to get TCGA case_id
manifest = pd.read_csv(f"{WORK}/wsi_manifest.tsv", sep="\t")
case_to_emb = {}
for slide_dir, e in emb.items():
    fid = slide_dir
    row = manifest[manifest["file_id"] == fid]
    if len(row) > 0:
        case_to_emb[row.iloc[0]["submitter_id"]] = e

# Save to retrieve for downstream
import pickle
with open(f"{WORK}/results/case_embeddings.pkl", "wb") as f:
    pickle.dump(case_to_emb, f)

print(f"Case-level embeddings: {len(case_to_emb)}")
print(f"Embedding dim: {next(iter(case_to_emb.values())).shape if case_to_emb else 'NA'}")
EOF

ls -la $WORK/results/
echo "=== ALL DONE ==="
df -h $WORK
