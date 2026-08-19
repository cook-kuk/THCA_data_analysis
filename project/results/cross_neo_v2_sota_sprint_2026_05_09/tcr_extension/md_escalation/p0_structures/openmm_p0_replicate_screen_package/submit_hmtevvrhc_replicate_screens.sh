#!/usr/bin/env bash
set -euo pipefail
mkdir -p screens

# HMTEVVRHC HLA-A*02:01 6VRM chain P
python run_openmm_pilot.py \
  --input inputs/6VRM_HMTEVVRHC_HLA-A0201_chainP.pdb \
  --outdir screens/6VRM_HMTEVVRHC_HLA-A0201_chainP_0p5ns \
  --mode explicit \
  --platform CUDA \
  --ns 0.5 \
  --report-steps 50000 \
  --minimize-iterations 1000 \
  --timestep-fs 1.0 \
  --temperature-k 300
python analyze_openmm_pilot.py \
  --top screens/6VRM_HMTEVVRHC_HLA-A0201_chainP_0p5ns/prepared_start.pdb \
  --traj screens/6VRM_HMTEVVRHC_HLA-A0201_chainP_0p5ns/trajectory.dcd \
  --peptide-chain P \
  --mhc-chain A \
  --tcr-chains 'E|D' \
  --peptide-sequence HMTEVVRHC \
  --out screens/6VRM_HMTEVVRHC_HLA-A0201_chainP_0p5ns/analysis.json

# HMTEVVRHC HLA-A*02:01 6VQO chain P
python run_openmm_pilot.py \
  --input inputs/6VQO_HMTEVVRHC_HLA-A0201_chainP.pdb \
  --outdir screens/6VQO_HMTEVVRHC_HLA-A0201_chainP_0p5ns \
  --mode explicit \
  --platform CUDA \
  --ns 0.5 \
  --report-steps 50000 \
  --minimize-iterations 1000 \
  --timestep-fs 1.0 \
  --temperature-k 300
python analyze_openmm_pilot.py \
  --top screens/6VQO_HMTEVVRHC_HLA-A0201_chainP_0p5ns/prepared_start.pdb \
  --traj screens/6VQO_HMTEVVRHC_HLA-A0201_chainP_0p5ns/trajectory.dcd \
  --peptide-chain P \
  --mhc-chain A \
  --tcr-chains 'E|D' \
  --peptide-sequence HMTEVVRHC \
  --out screens/6VQO_HMTEVVRHC_HLA-A0201_chainP_0p5ns/analysis.json

# HMTEVVRHC HLA-A*02:01 6VQO chain Q
python run_openmm_pilot.py \
  --input inputs/6VQO_HMTEVVRHC_HLA-A0201_chainQ.pdb \
  --outdir screens/6VQO_HMTEVVRHC_HLA-A0201_chainQ_0p5ns \
  --mode explicit \
  --platform CUDA \
  --ns 0.5 \
  --report-steps 50000 \
  --minimize-iterations 1000 \
  --timestep-fs 1.0 \
  --temperature-k 300
python analyze_openmm_pilot.py \
  --top screens/6VQO_HMTEVVRHC_HLA-A0201_chainQ_0p5ns/prepared_start.pdb \
  --traj screens/6VQO_HMTEVVRHC_HLA-A0201_chainQ_0p5ns/trajectory.dcd \
  --peptide-chain Q \
  --mhc-chain F \
  --tcr-chains 'J|H' \
  --peptide-sequence HMTEVVRHC \
  --out screens/6VQO_HMTEVVRHC_HLA-A0201_chainQ_0p5ns/analysis.json

# HMTEVVRHC HLA-A*02:01 7RM4 chain H
python run_openmm_pilot.py \
  --input inputs/7RM4_HMTEVVRHC_HLA-A0201_chainH.pdb \
  --outdir screens/7RM4_HMTEVVRHC_HLA-A0201_chainH_0p5ns \
  --mode explicit \
  --platform CUDA \
  --ns 0.5 \
  --report-steps 50000 \
  --minimize-iterations 1000 \
  --timestep-fs 1.0 \
  --temperature-k 300
python analyze_openmm_pilot.py \
  --top screens/7RM4_HMTEVVRHC_HLA-A0201_chainH_0p5ns/prepared_start.pdb \
  --traj screens/7RM4_HMTEVVRHC_HLA-A0201_chainH_0p5ns/trajectory.dcd \
  --peptide-chain H \
  --mhc-chain F \
  --tcr-chains 'J|I' \
  --peptide-sequence HMTEVVRHC \
  --out screens/7RM4_HMTEVVRHC_HLA-A0201_chainH_0p5ns/analysis.json

# HMTEVVRHC HLA-A*02:01 7RM4 chain C
python run_openmm_pilot.py \
  --input inputs/7RM4_HMTEVVRHC_HLA-A0201_chainC.pdb \
  --outdir screens/7RM4_HMTEVVRHC_HLA-A0201_chainC_0p5ns \
  --mode explicit \
  --platform CUDA \
  --ns 0.5 \
  --report-steps 50000 \
  --minimize-iterations 1000 \
  --timestep-fs 1.0 \
  --temperature-k 300
python analyze_openmm_pilot.py \
  --top screens/7RM4_HMTEVVRHC_HLA-A0201_chainC_0p5ns/prepared_start.pdb \
  --traj screens/7RM4_HMTEVVRHC_HLA-A0201_chainC_0p5ns/trajectory.dcd \
  --peptide-chain C \
  --mhc-chain A \
  --tcr-chains 'E|D' \
  --peptide-sequence HMTEVVRHC \
  --out screens/7RM4_HMTEVVRHC_HLA-A0201_chainC_0p5ns/analysis.json

# HMTEVVRHC HLA-A*02:01 7RM4 chain R
python run_openmm_pilot.py \
  --input inputs/7RM4_HMTEVVRHC_HLA-A0201_chainR.pdb \
  --outdir screens/7RM4_HMTEVVRHC_HLA-A0201_chainR_0p5ns \
  --mode explicit \
  --platform CUDA \
  --ns 0.5 \
  --report-steps 50000 \
  --minimize-iterations 1000 \
  --timestep-fs 1.0 \
  --temperature-k 300
python analyze_openmm_pilot.py \
  --top screens/7RM4_HMTEVVRHC_HLA-A0201_chainR_0p5ns/prepared_start.pdb \
  --traj screens/7RM4_HMTEVVRHC_HLA-A0201_chainR_0p5ns/trajectory.dcd \
  --peptide-chain R \
  --mhc-chain P \
  --tcr-chains 'T|S' \
  --peptide-sequence HMTEVVRHC \
  --out screens/7RM4_HMTEVVRHC_HLA-A0201_chainR_0p5ns/analysis.json

# HMTEVVRHC HLA-A*02:01 7RM4 chain M
python run_openmm_pilot.py \
  --input inputs/7RM4_HMTEVVRHC_HLA-A0201_chainM.pdb \
  --outdir screens/7RM4_HMTEVVRHC_HLA-A0201_chainM_0p5ns \
  --mode explicit \
  --platform CUDA \
  --ns 0.5 \
  --report-steps 50000 \
  --minimize-iterations 1000 \
  --timestep-fs 1.0 \
  --temperature-k 300
python analyze_openmm_pilot.py \
  --top screens/7RM4_HMTEVVRHC_HLA-A0201_chainM_0p5ns/prepared_start.pdb \
  --traj screens/7RM4_HMTEVVRHC_HLA-A0201_chainM_0p5ns/trajectory.dcd \
  --peptide-chain M \
  --mhc-chain K \
  --tcr-chains 'O|N' \
  --peptide-sequence HMTEVVRHC \
  --out screens/7RM4_HMTEVVRHC_HLA-A0201_chainM_0p5ns/analysis.json
