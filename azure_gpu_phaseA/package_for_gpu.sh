#!/usr/bin/env bash
# package_for_gpu.sh — build a self-contained tarball for Phase A GPU run.
#
# Output:
#   /home/seungho/personal/THCA_data_analysis/azure_gpu_phaseA/phaseA_gpu_pkg.tar.gz
#
# Layout inside the tarball (untars to a single dir "phaseA_gpu_pkg/"):
#   phaseA_gpu_pkg/
#     ├── README_AZURE_GPU_PHASEA.md
#     ├── requirements_pathology_gpu.txt
#     ├── run_phaseA_gpu.sh
#     └── thca_phaseA/                       ← THCA_ROOT on remote
#         └── project/
#             ├── src/05_pathology_poc/...   (5 .py scripts)
#             ├── data/processed/GSE250521/
#             │   ├── tiles/ (3,200 × 224 px PNG)
#             │   └── GSM*/*.scored.h5ad     (16 files, for neg ctrl)
#             ├── results/01_spatial_score/all_spots_scored.tsv.gz
#             ├── results/03_pathology_poc/tile_metadata*.tsv.gz
#             └── reports/                   (empty, populated on remote)
set -euo pipefail
SRC="/home/seungho/personal/THCA_data_analysis"
PKG_ROOT="${SRC}/azure_gpu_phaseA"
STAGE="${PKG_ROOT}/.stage/phaseA_gpu_pkg"
TGZ="${PKG_ROOT}/phaseA_gpu_pkg.tar.gz"

echo "[pkg] cleaning stage…"
rm -rf "${PKG_ROOT}/.stage" "${TGZ}"
mkdir -p "${STAGE}/thca_phaseA/project/src/05_pathology_poc"
mkdir -p "${STAGE}/thca_phaseA/project/data/processed/GSE250521/tiles"
mkdir -p "${STAGE}/thca_phaseA/project/results/01_spatial_score"
mkdir -p "${STAGE}/thca_phaseA/project/results/03_pathology_poc"
mkdir -p "${STAGE}/thca_phaseA/project/reports"

echo "[pkg] copying scripts…"
cp "${SRC}"/project/src/05_pathology_poc/{compute_resid_labels,embed_resnet50,train_loso_ridge,negative_controls,plot_pred_vs_obs,extract_tiles}.py \
   "${STAGE}/thca_phaseA/project/src/05_pathology_poc/"

echo "[pkg] copying tile metadata + spot scores…"
cp "${SRC}"/project/results/01_spatial_score/all_spots_scored.tsv.gz \
   "${STAGE}/thca_phaseA/project/results/01_spatial_score/"
cp "${SRC}"/project/results/03_pathology_poc/tile_metadata.tsv.gz \
   "${SRC}"/project/results/03_pathology_poc/tile_metadata_resid.tsv.gz \
   "${STAGE}/thca_phaseA/project/results/03_pathology_poc/"

echo "[pkg] copying scored h5ad (16 files)…"
for h in "${SRC}"/project/data/processed/GSE250521/*/GSM*.scored.h5ad; do
  sid=$(basename "$(dirname "$h")")
  mkdir -p "${STAGE}/thca_phaseA/project/data/processed/GSE250521/${sid}"
  cp "$h" "${STAGE}/thca_phaseA/project/data/processed/GSE250521/${sid}/"
done

echo "[pkg] copying tiles (224 px)…"
cp -r "${SRC}"/project/data/processed/GSE250521/tiles/* \
      "${STAGE}/thca_phaseA/project/data/processed/GSE250521/tiles/"

echo "[pkg] copying top-level files…"
cp "${PKG_ROOT}/README_AZURE_GPU_PHASEA.md" "${STAGE}/"
cp "${PKG_ROOT}/requirements_pathology_gpu.txt" "${STAGE}/"
cp "${PKG_ROOT}/run_phaseA_gpu.sh" "${STAGE}/"
chmod +x "${STAGE}/run_phaseA_gpu.sh"

echo "[pkg] tarball…"
cd "${PKG_ROOT}/.stage"
tar czf "${TGZ}" phaseA_gpu_pkg
echo
echo "[pkg] DONE"
ls -lh "${TGZ}"
echo
echo "stage layout:"
find "${STAGE}" -maxdepth 4 -type d | sed "s|${STAGE}|.|g"
echo
echo "=== Send to GPU VM ==="
echo "scp ${TGZ} azureuser@<GPU_VM_IP>:~/"
echo "ssh azureuser@<GPU_VM_IP> 'tar xzf ~/phaseA_gpu_pkg.tar.gz && cd phaseA_gpu_pkg && bash run_phaseA_gpu.sh'"
