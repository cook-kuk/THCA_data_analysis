#!/bin/bash
# Pod DM Phase 6 fix — relax HoVer-NeXt deps, skip on failure.
set -uo pipefail
WORK=/workspace/wsi_pathology_dm

echo "=== [Phase 6 fix] HoVer-NeXt with relaxed deps ==="
cd "$WORK"
[ ! -d hover_next_inference ] && git clone --depth 1 https://github.com/digitalpathologybern/hover_next_inference 2>&1 | tail -2
cd hover_next_inference

# Use pip --no-deps for HoVer-NeXt + manually install relaxed compatible deps
apt-get install -y -qq python3-dev gcc g++ 2>&1 | tail -1

pip install -q --no-deps \
    "segmentation-models-pytorch>=0.3.0" \
    "timm>=0.9.6" \
    "torch_geometric>=2.5" 2>&1 | tail -3

pip install -q \
    pretrainedmodels efficientnet-pytorch fasteners asciitree \
    "albumentations>=1.3" \
    pylibCZIrw \
    geojson shapely \
    "zarr>=2.16,<3" "numcodecs>=0.12,<0.14" \
    tifffile h5py imagecodecs networkx scikit-image staintools \
    "spams-bin>=2.6" 2>&1 | tail -3

# Try mahotas — if fails, skip HoVer-NeXt entirely
pip install -q mahotas 2>&1 | tail -3 || pip install -q mahotas==1.4.18 2>&1 | tail -2 || echo "mahotas install failed"

# weights
if [ ! -d lizard_convnextv2_large ]; then
    echo "downloading HoVer-NeXt weights..."
    curl -sL "https://zenodo.org/records/10635618/files/lizard_convnextv2_large.zip?download=1" -o lcl.zip
    unzip -oq lcl.zip && rm lcl.zip
fi

# Test imports first
python3 -c "
import sys
sys.path.insert(0, '.')
sys.path.insert(0, 'src')
import importlib
fails = []
for m in ['mahotas', 'segmentation_models_pytorch', 'timm', 'torch_geometric',
          'pretrainedmodels', 'staintools', 'spams']:
    try:
        importlib.import_module(m)
    except Exception as e:
        fails.append((m, str(e)[:60]))
print('imports OK except:', fails)
"

mkdir -p "$WORK/hovernext_out"
ls "$WORK"/wsi/*/*.svs > "$WORK/wsi_paths.txt" 2>/dev/null
N_SVS=$(wc -l < "$WORK/wsi_paths.txt")
N_CORES=$(nproc)
echo "running HoVer-NeXt on $N_SVS WSIs..."

# Run subset of 5 WSIs as proof-of-concept first; if successful expand
head -5 "$WORK/wsi_paths.txt" > "$WORK/wsi_paths_5.txt"

python3 main.py \
    --input "$WORK/wsi_paths_5.txt" \
    --output_root "$WORK/hovernext_out/" \
    --cp lizard_convnextv2_large \
    --tta 2 \
    --inf_workers "$N_CORES" \
    --pp_tiling 8 \
    --pp_workers "$((N_CORES - 1))" 2>&1 | tail -30 || echo "HoVer-NeXt POC failed; logging and continuing"

ls "$WORK/hovernext_out/" 2>/dev/null

# ------------------------------------------------------------------
# Phase 8 — Pack artifacts (DINOv2 LOSO is the main result regardless)
# ------------------------------------------------------------------
echo "=== [Phase 8] Pack artifacts ==="
cd "$WORK"
ls -la results/ embeddings/ 2>/dev/null
du -sh wsi tiles embeddings results geojson hovernext_out 2>/dev/null

tar czf "$WORK/dm_wsi_artifacts.tgz" \
    -C "$WORK" \
    embeddings results wsi_resolved_manifest.tsv hovernext_out 2>/dev/null || \
tar czf "$WORK/dm_wsi_artifacts.tgz" \
    -C "$WORK" embeddings results wsi_resolved_manifest.tsv

ls -la "$WORK/dm_wsi_artifacts.tgz"
echo "=== ALL DONE ==="
