#!/bin/bash
# RunPod environment setup — Foundation-model + CLAM v2
set -euo pipefail

echo "[1/5] System packages"
apt-get update -qq
apt-get install -y -qq \
  libopenslide-dev libgl1 libglib2.0-0 git wget tar \
  python3-pip

echo "[2/5] Python packages"
pip install -q --upgrade pip
pip install -q \
  torch==2.4.0 torchvision \
  timm transformers huggingface_hub \
  openslide-python pillow tqdm \
  scikit-learn scipy pandas numpy \
  matplotlib seaborn \
  requests

echo "[3/5] CLAM utility (tile extraction etc.)"
# CLAM official: https://github.com/mahmoodlab/CLAM
# Use minimal pieces inline (we have phase2_clam_train_eval.py)
pip install -q openslide-python

echo "[4/5] SAM weights (Phase 1b cell segmentation, optional)"
mkdir -p /workspace/sam_weights
if [ ! -f /workspace/sam_weights/sam_vit_h_4b8939.pth ]; then
  wget -q -O /workspace/sam_weights/sam_vit_h_4b8939.pth \
    https://dl.fbaipublicfiles.com/segment_anything/sam_vit_h_4b8939.pth
fi
echo "  SAM ViT-H weights ready"

echo "[5/5] Ready"
echo ""
echo "Next: export HF_TOKEN=... ; bash run_phase1.sh"
