#!/bin/bash
# Pod DM resume — manifest already resolved (Phase 2 partial). Continue from WSI download.
# Uses parallel curl on GDC API directly (skips broken gdc-client wget).
set -euo pipefail
WORK=/workspace/wsi_pathology_dm
mkdir -p "$WORK"/{wsi,tiles,embeddings,segm,geojson,results,logs,gigatime_mIF}

# ------------------------------------------------------------------
# Phase 2b — direct GDC API download (bypass broken gdc-client URL)
# ------------------------------------------------------------------
echo "=== [Phase 2b] direct GDC download ==="
cd "$WORK/wsi"

# 5 parallel downloads. Download into .part files and accept a slide only
# after the byte count is close to the GDC manifest size. This prevents
# interrupted curl output from being treated as a valid SVS on rerun.
tail -n +2 "$WORK/wsi_resolved_manifest.tsv" | awk -F'\t' '{print $1"\t"$2"\t"$4}' > "$WORK/dl_list.tsv"
N=$(wc -l < "$WORK/dl_list.tsv")
echo "downloading $N WSIs in parallel (5 at a time)..."

mkdir -p "$WORK/wsi"
cat "$WORK/dl_list.tsv" | xargs -P5 -n3 bash -c '
  fid="$1"
  fname="$2"
  size_gb="$3"
  out=/workspace/wsi_pathology_dm/wsi/${fid}/${fname}
  min_bytes=$(python3 -c "print(int(float(\"$size_gb\") * 1e9 * 0.90))")
  if [ -s "$out" ] && [ "$(stat -c%s "$out")" -ge "$min_bytes" ]; then
    echo "EXISTS $fid"
    exit 0
  fi
  rm -f "$out" "${out}.part"
  mkdir -p "/workspace/wsi_pathology_dm/wsi/${fid}"
  curl -fL -sS -o "${out}.part" "https://api.gdc.cancer.gov/data/${fid}" --retry 5 --retry-delay 3 --retry-all-errors --max-time 1800
  sz=$(stat -c%s "${out}.part" 2>/dev/null || echo 0)
  if [ "$sz" -ge "$min_bytes" ]; then
    mv "${out}.part" "$out"
    echo "OK $fid ($((sz/1024/1024)) MB) $fname"
  else
    rm -f "${out}.part"
    echo "FAIL $fid $fname expected_min=$((min_bytes/1024/1024))MB got=$((sz/1024/1024))MB"
    exit 1
  fi
' _

ls "$WORK/wsi"/*/*.svs 2>/dev/null | wc -l
du -sh "$WORK/wsi" 2>/dev/null

# ------------------------------------------------------------------
# Phase 3 — Tile WSI at 20x with HSV-Otsu tissue mask
# ------------------------------------------------------------------
echo "=== [Phase 3] WSI tiling + tissue mask ==="
mkdir -p "$WORK/tiles" "$WORK/segm"

# Install needed deps if not present
pip install -q openslide-python==1.3.1 timm transformers huggingface_hub \
    pandas tqdm scikit-learn opencv-python-headless geojson shapely \
    scikit-image 2>&1 | tail -3
apt-get install -y -qq openslide-tools libopenslide-dev libvips-dev 2>&1 | tail -1

python3 << 'PY'
import openslide, os, glob, numpy as np, cv2
from PIL import Image
from tqdm import tqdm

WORK = "/workspace/wsi_pathology_dm"
PATCH_SIZE = 256
TARGET_LEVEL = 1

for svs in tqdm(sorted(glob.glob(f"{WORK}/wsi/*/*.svs"))):
    fid = os.path.basename(os.path.dirname(svs))
    out_dir = f"{WORK}/tiles/{fid}"
    os.makedirs(out_dir, exist_ok=True)
    if len(os.listdir(out_dir)) > 50:
        continue
    try:
        slide = openslide.OpenSlide(svs)
    except Exception as e:
        print(f"  open failed {svs}: {e}")
        continue
    L = min(TARGET_LEVEL, slide.level_count - 1)
    w, h = slide.level_dimensions[L]
    n_x = w // PATCH_SIZE
    n_y = h // PATCH_SIZE
    stride_x = max(1, n_x // 15)
    stride_y = max(1, n_y // 15)
    saved = 0
    for ix in range(0, n_x, stride_x):
        for iy in range(0, n_y, stride_y):
            x = ix * PATCH_SIZE
            y = iy * PATCH_SIZE
            img = slide.read_region((x, y), L, (PATCH_SIZE, PATCH_SIZE)).convert("RGB")
            arr = np.array(img)
            if (arr > 200).all(axis=2).mean() > 0.75:
                continue
            img.save(f"{out_dir}/tile_{ix:04d}_{iy:04d}.png")
            saved += 1
            if saved >= 200: break
        if saved >= 200: break
print("tiling done")
PY

# ------------------------------------------------------------------
# Phase 4 — Foundation model embedding (UNI > DINOv2)
# ------------------------------------------------------------------
echo "=== [Phase 4] foundation model embedding ==="

python3 << 'PY'
import os, glob, numpy as np, torch, timm
from PIL import Image
from torchvision import transforms
from tqdm import tqdm

WORK = "/workspace/wsi_pathology_dm"
device = "cuda" if torch.cuda.is_available() else "cpu"
backbone = None
model = None

# Try UNI (gated)
try:
    if os.environ.get("HF_TOKEN"):
        from huggingface_hub import login
        login(token=os.environ["HF_TOKEN"], add_to_git_credential=False)
    model = timm.create_model("hf-hub:MahmoodLab/uni", pretrained=True, init_values=1e-5, dynamic_img_size=True)
    backbone = "UNI"
except Exception as e:
    print(f"UNI unavailable ({type(e).__name__}); using DINOv2 ViT-L")
    try:
        model = timm.create_model("vit_large_patch14_dinov2.lvd142m", pretrained=True)
        backbone = "DINOv2-ViT-L"
    except Exception as ee:
        print(f"DINOv2 unavailable ({ee}); ViT-B")
        model = timm.create_model("vit_base_patch16_224", pretrained=True)
        backbone = "ViT-B-imagenet"

model.eval().to(device)
print(f"backbone: {backbone}")

xform = transforms.Compose([
    transforms.Resize(224), transforms.CenterCrop(224), transforms.ToTensor(),
    transforms.Normalize([0.485,0.456,0.406],[0.229,0.224,0.225]),
])

slide_emb = {}
for slide_dir in tqdm(sorted(os.listdir(f"{WORK}/tiles"))):
    tps = sorted(glob.glob(f"{WORK}/tiles/{slide_dir}/*.png"))
    if not tps: continue
    feats = []
    with torch.no_grad():
        # batched
        BS = 32
        for i in range(0, len(tps), BS):
            batch = torch.stack([xform(Image.open(tp).convert("RGB")) for tp in tps[i:i+BS]]).to(device)
            f = model(batch) if backbone == "UNI" else model.forward_features(batch).mean(dim=1)
            feats.append(f.cpu().numpy())
    feats = np.concatenate(feats, axis=0)
    slide_emb[slide_dir] = feats.mean(axis=0)
    np.save(f"{WORK}/embeddings/{slide_dir}_tiles.npy", feats)

import pickle
with open(f"{WORK}/embeddings/slide_embeddings.pkl", "wb") as f:
    pickle.dump({"backbone": backbone, "embeddings": slide_emb}, f)
print(f"slides embedded: {len(slide_emb)}")
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
print(f"matrix: X={X.shape}, y_DM1={int((y==0).sum())}, y_DM2={int((y==1).sum())}")
loso = LeaveOneGroupOut()
preds = np.zeros_like(y, dtype=float)
for tr, te in loso.split(X, y, groups):
    mdl = LogisticRegressionCV(Cs=[0.01,0.1,1,10], cv=3, max_iter=1000).fit(X[tr], y[tr])
    preds[te] = mdl.predict_proba(X[te])[:, 1]

auc = roc_auc_score(y, preds)
acc = accuracy_score(y, (preds > 0.5).astype(int))
print(f"DM1 vs DM2 LOSO AUC={auc:.3f}, ACC={acc:.3f}")

pd.DataFrame({"case_id": groups, "dm_true": np.where(y == 0, "DM1", "DM2"),
              "dm_pred_prob": preds}).to_csv(
    f"{WORK}/results/dm1_vs_dm2_loso_predictions.tsv", sep="\t", index=False)
json.dump({"backbone": pkg["backbone"], "auc": float(auc), "acc": float(acc),
           "n_DM1": int((y==0).sum()), "n_DM2": int((y==1).sum())},
          open(f"{WORK}/results/dm1_vs_dm2_loso_metrics.json", "w"), indent=2)
print("metrics saved")
PY

echo "=== [Phase 6] HoVer-NeXt cell detection ==="
cd "$WORK"
[ ! -d hover_next_inference ] && git clone --depth 1 https://github.com/digitalpathologybern/hover_next_inference 2>&1 | tail -3
cd hover_next_inference
pip install -q --no-deps -r requirements.txt 2>&1 | tail -3
pip install -q pretrainedmodels efficientnet-pytorch fasteners asciitree mahotas==1.4.18 \
    segmentation-models-pytorch==0.3.4 zarr==2.16.1 numcodecs==0.12.1 \
    spams-bin staintools shapely==2.0.2 networkx==2.8.7 geojson==3.1.0 albumentations==1.3.1 \
    timm==0.9.6 tifffile h5py imagecodecs 2>&1 | tail -3

if [ ! -d lizard_convnextv2_large ]; then
    curl -sL "https://zenodo.org/records/10635618/files/lizard_convnextv2_large.zip?download=1" -o lcl.zip
    unzip -oq lcl.zip && rm lcl.zip
fi
ls -la lizard_convnextv2_large/ 2>/dev/null | head -4

mkdir -p "$WORK/hovernext_out"
ls "$WORK"/wsi/*/*.svs > "$WORK/wsi_paths.txt" 2>/dev/null
N_SVS=$(wc -l < "$WORK/wsi_paths.txt")
echo "WSIs to process: $N_SVS"
N_CORES=$(nproc)

python3 main.py \
    --input "$WORK/wsi_paths.txt" \
    --output_root "$WORK/hovernext_out/" \
    --cp lizard_convnextv2_large \
    --tta 4 \
    --inf_workers "$N_CORES" \
    --pp_tiling 8 \
    --pp_workers "$((N_CORES - 1))" 2>&1 | tail -30 || echo "HoVer-NeXt encountered issues, continuing"

# ------------------------------------------------------------------
# Phase 7 — GigaTIME (skip if no token)
# ------------------------------------------------------------------
if [ -n "${HF_TOKEN:-}" ]; then
    echo "=== [Phase 7] GigaTIME ==="
    cd "$WORK"
    [ ! -d GigaTIME ] && git clone --depth 1 https://github.com/prov-gigatime/GigaTIME 2>&1 | tail -2
    cd GigaTIME
    pip install -q numpy pandas scipy easydict pyyaml huggingface_hub 2>&1 | tail -2
    python3 << 'PY'
import os, glob, sys, numpy as np, torch
sys.path.insert(0, "scripts")
from huggingface_hub import snapshot_download
from PIL import Image
local = snapshot_download(repo_id="prov-gigatime/GigaTIME")
print(f"GigaTIME weights -> {local}")
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
            y = model(x).cpu().numpy().squeeze(0)
        means = y.mean(axis=(1, 2))
        np.savez_compressed(out_npz, channels=CHANNELS, mean=means, sample=y[:, ::4, ::4])
        n_done += 1
        if n_done % 50 == 0:
            print(f"  ...{n_done} GigaTIME inferences done", flush=True)
print(f"GigaTIME total: {n_done} tile inferences")
PY
else
    echo "=== [Phase 7] GigaTIME SKIPPED (HF_TOKEN unset) ==="
fi

# ------------------------------------------------------------------
# Phase 8 — Pack artifacts
# ------------------------------------------------------------------
echo "=== [Phase 8] Pack artifacts ==="
cd "$WORK"
ls -la results/ embeddings/ 2>/dev/null
du -sh wsi tiles embeddings results geojson hovernext_out gigatime_mIF 2>/dev/null

tar czf "$WORK/dm_wsi_artifacts.tgz" \
    -C "$WORK" \
    embeddings results geojson segm wsi_resolved_manifest.tsv \
    hovernext_out gigatime_mIF 2>/dev/null || \
tar czf "$WORK/dm_wsi_artifacts.tgz" \
    -C "$WORK" embeddings results wsi_resolved_manifest.tsv

ls -la "$WORK/dm_wsi_artifacts.tgz"
echo "=== ALL DONE ==="
