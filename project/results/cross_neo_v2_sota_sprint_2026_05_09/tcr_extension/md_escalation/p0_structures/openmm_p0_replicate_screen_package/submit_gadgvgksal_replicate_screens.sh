#!/usr/bin/env bash
set -euo pipefail
mkdir -p screens

# GADGVGKSAL HLA-C*08:02 6UON chain C
python run_openmm_pilot.py \
  --input inputs/6UON_GADGVGKSAL_HLA-C0802_chainC.pdb \
  --outdir screens/6UON_GADGVGKSAL_HLA-C0802_chainC_1p0ns \
  --mode explicit \
  --platform CUDA \
  --ns 1.0 \
  --report-steps 25000 \
  --minimize-iterations 1000 \
  --timestep-fs 2.0 \
  --temperature-k 300
python analyze_openmm_pilot.py \
  --top screens/6UON_GADGVGKSAL_HLA-C0802_chainC_1p0ns/prepared_start.pdb \
  --traj screens/6UON_GADGVGKSAL_HLA-C0802_chainC_1p0ns/trajectory.dcd \
  --peptide-chain C \
  --mhc-chain A \
  --tcr-chains 'J|I' \
  --peptide-sequence GADGVGKSAL \
  --out screens/6UON_GADGVGKSAL_HLA-C0802_chainC_1p0ns/analysis.json
