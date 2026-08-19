#!/usr/bin/env bash
set -euo pipefail
mkdir -p runs

# GADGVGKSAL HLA-C*08:02 6UON
python run_openmm_pilot.py \
  --input inputs/6UON_GADGVGKSAL_HLA-C0802.minimized.pdb \
  --outdir runs/6UON_GADGVGKSAL_HLA-C0802_10ns \
  --mode explicit \
  --platform CUDA \
  --precision mixed \
  --ns 10 \
  --report-steps 5000

# HMTEVVRHC HLA-A*02:01 6VRN
python run_openmm_pilot.py \
  --input inputs/6VRN_HMTEVVRHC_HLA-A0201.minimized.pdb \
  --outdir runs/6VRN_HMTEVVRHC_HLA-A0201_10ns \
  --mode explicit \
  --platform CUDA \
  --precision mixed \
  --ns 10 \
  --report-steps 5000
