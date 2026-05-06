#!/bin/bash
# RunPod setup for Phase B GPU package
set -e

echo "[setup] Python deps"
pip install -q --upgrade pip
pip install -q numpy pandas scipy scikit-learn matplotlib seaborn anndata scanpy
pip install -q torch torchvision  # if not already in RunPod base
pip install -q openslide-python tiatoolbox 2>&1 || echo "  (openslide system pkg may be needed)"

echo "[setup] timm + huggingface (for UNI/CONCH if used)"
pip install -q timm huggingface_hub

echo "[setup] System packages (apt) — may need sudo"
if command -v apt-get >/dev/null 2>&1; then
  sudo apt-get install -y openslide-tools 2>&1 | tail -3 || echo "  apt-get failed (no sudo?)"
fi

echo "[setup] gdc-client for TCGA download"
if ! command -v gdc-client >/dev/null 2>&1; then
  cd /tmp
  wget -q https://gdc.cancer.gov/files/public/file/gdc-client_v1.6.1_Ubuntu_x64.zip
  unzip -q gdc-client_v1.6.1_Ubuntu_x64.zip
  sudo mv gdc-client /usr/local/bin/ 2>&1 || mv gdc-client ~/.local/bin/
  cd -
fi

echo "[setup] verify"
python3 -c "import torch; print('torch', torch.__version__, 'cuda', torch.cuda.is_available(), 'gpu:', torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'none')"
python3 -c "import scanpy, anndata; print('scanpy', scanpy.__version__)"
gdc-client --version 2>&1 || echo "  gdc-client missing"
echo "[setup] done"
