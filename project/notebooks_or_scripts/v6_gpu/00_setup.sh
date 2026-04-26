#!/usr/bin/env bash
# v6 Wave 2 setup — pin GPU stack for scGPT / Geneformer / GEARS / CellRank / CellOracle / VEGA.
# Idempotent: safe to re-run.
set -euo pipefail
LOG=${LOG:-/opt/thyroid-dash/project/logs/v6_gpu_setup.log}
mkdir -p "$(dirname "$LOG")"
exec > >(tee -a "$LOG") 2>&1

echo "=== v6 Wave 2 GPU setup begins $(date -u) ==="

if ! command -v nvidia-smi >/dev/null 2>&1; then
  echo "FATAL: nvidia-smi not found. This bundle requires a CUDA-capable GPU."
  echo "Provision a GPU instance (Colab Pro+, Paperspace, Lambda, AWS p3/g5) and re-run."
  exit 1
fi
nvidia-smi | head -20

PROJECT=${PROJECT:-/opt/thyroid-dash/project}
VENV=${VENV:-$PROJECT/.venv}
if [[ ! -d "$VENV" ]]; then
  python3 -m venv "$VENV"
fi
# shellcheck source=/dev/null
. "$VENV/bin/activate"
python3 -m pip install --upgrade pip wheel >/dev/null

echo "--- core stack ---"
pip install --quiet \
  torch==2.2.* \
  transformers==4.44.* \
  accelerate==0.33.* \
  huggingface-hub==0.24.* \
  scanpy==1.10.* anndata==0.10.* leidenalg==0.10.* \
  scvi-tools==1.1.* harmonypy==0.0.10 \
  cellrank==2.0.* scvelo==0.3.* \
  celloracle==0.18.* \
  gseapy==1.1.* decoupler==1.8.*

echo "--- scGPT (git) ---"
pip install --quiet git+https://github.com/bowang-lab/scGPT.git || \
  echo "WARN: scGPT git install failed; scripts will exit 2 if scGPT not importable"

echo "--- Geneformer (HF) ---"
pip install --quiet git+https://huggingface.co/ctheodoris/Geneformer || \
  echo "WARN: Geneformer git install failed; transformers AutoModel fallback in 20_geneformer.py"

echo "--- GEARS ---"
pip install --quiet git+https://github.com/snap-stanford/GEARS.git || \
  echo "WARN: GEARS install failed"

echo "--- VEGA ---"
pip install --quiet vega-tools 2>/dev/null || \
  pip install --quiet git+https://github.com/LucasESBS/vega.git || \
  echo "WARN: VEGA install failed"

echo "=== versions ==="
python3 -c "import sys; print('python', sys.version.split()[0])"
python3 -c "import torch; print('torch', torch.__version__, 'cuda', torch.cuda.is_available(), 'devices', torch.cuda.device_count())"
for pkg in scanpy anndata scvi cellrank scvelo celloracle transformers; do
  python3 -c "import importlib; m = importlib.import_module('$pkg'); print('$pkg', getattr(m,'__version__','?'))" || echo "$pkg missing"
done

echo "=== v6 Wave 2 GPU setup done $(date -u) ==="
