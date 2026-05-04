#!/bin/bash
# Phase B Drug Discovery — RunPod environment setup
# Author: Seungho Cook · 2026-05-04
# Pre-req: A100 / A40 / H100 pod with CUDA 12.1+
set -euo pipefail

echo "=== [setup] Phase B drug discovery environment ==="
DEBIAN_FRONTEND=noninteractive apt-get update -qq
DEBIAN_FRONTEND=noninteractive apt-get install -y -qq \
    git wget curl build-essential cmake \
    python3-pip python3-dev python3-venv \
    libffi-dev libssl-dev libxml2-dev libxslt1-dev \
    libopenbabel7 openbabel pymol \
    autodock-vina 2>&1 | tail -3

# Python deps
pip install --upgrade pip
pip install -q \
    torch==2.3.0 torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121
pip install -q \
    biopython==1.84 \
    rdkit==2024.3.1 \
    openff-toolkit==0.16.0 \
    openff-units==0.2.0 \
    pdbfixer==1.10.0 \
    mdanalysis==2.7.0 \
    ProLIF==2.0.3 \
    py3Dmol==2.1.0

# Specialized DL packages
pip install -q \
    transformers==4.40.0 \
    chemberta-2 \
    deepchem==2.8.0 \
    chemprop==1.7.1 \
    deeppurpose==0.1.5

# DiffDock (separate install)
git clone https://github.com/gcorso/DiffDock.git /workspace/DiffDock || true
cd /workspace/DiffDock && pip install -q -r requirements.txt && cd -

# GNINA (binary)
wget -qO /usr/local/bin/gnina https://github.com/gnina/gnina/releases/download/v1.1/gnina && chmod +x /usr/local/bin/gnina

# AlphaFold3 (Anthropic-cloud or local)
echo "[setup] AlphaFold3: install via DeepMind official → https://github.com/google-deepmind/alphafold3"

# OpenFE (FEP)
pip install -q openfe==1.0.1

# PROTAC-DB
mkdir -p /workspace/protac_db
wget -qO /workspace/protac_db/protac_db.csv https://protacdb.weizmann.ac.il/api/protacs.csv || true

# SwissTargetPrediction REST API key (optional)
echo "[setup] SwissTargetPrediction: REST API ready at http://www.swisstargetprediction.ch/"

# Working directories
mkdir -p /workspace/azure_gpu_phaseB_drug/{results,logs,configs}
mkdir -p /workspace/azure_gpu_phaseB_drug/results/{b1_structures,b2_diffdock,b3_gnina,b4_kdeep,b5_chemberta,b6_protac,b7_adc_optimization,b8_off_target,b9_admet_foundation,b10_fep_md}

echo "=== [setup] complete ==="
echo "Verify: nvidia-smi · python3 -c 'import torch; print(torch.cuda.is_available())' · gnina --version"
