#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"
mkdir -p runs logs

echo 'tier0 short template-diversity screens'

echo '[job] tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_6VRM_HMTEVVRHC_HLA-A0201_chainP_rep01_0p5ns'
if [ -f runs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_6VRM_HMTEVVRHC_HLA-A0201_chainP_rep01_0p5ns/run_metadata.json ]; then
  echo '[skip] tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_6VRM_HMTEVVRHC_HLA-A0201_chainP_rep01_0p5ns already has run_metadata.json'
else
  python run_openmm_pilot.py --input inputs/6VRM_HMTEVVRHC_HLA-A0201_chainP_e59660fd.pdb --outdir runs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_6VRM_HMTEVVRHC_HLA-A0201_chainP_rep01_0p5ns --mode explicit --platform CUDA --ns 0.5 --report-steps 50000 --minimize-iterations 1000 --timestep-fs 1 --temperature-k 300 --random-seed 10002 2>&1 | tee logs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_6VRM_HMTEVVRHC_HLA-A0201_chainP_rep01_0p5ns.openmm.log
fi
if [ -f runs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_6VRM_HMTEVVRHC_HLA-A0201_chainP_rep01_0p5ns/trajectory.dcd ]; then
  python analyze_openmm_pilot.py --top runs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_6VRM_HMTEVVRHC_HLA-A0201_chainP_rep01_0p5ns/prepared_start.pdb --traj runs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_6VRM_HMTEVVRHC_HLA-A0201_chainP_rep01_0p5ns/trajectory.dcd --peptide-chain P --mhc-chain A --tcr-chains 'E|D' --peptide-sequence HMTEVVRHC --out runs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_6VRM_HMTEVVRHC_HLA-A0201_chainP_rep01_0p5ns/analysis.json 2>&1 | tee logs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_6VRM_HMTEVVRHC_HLA-A0201_chainP_rep01_0p5ns.analysis.log
fi

echo '[job] tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_6VQO_HMTEVVRHC_HLA-A0201_chainP_rep01_0p5ns'
if [ -f runs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_6VQO_HMTEVVRHC_HLA-A0201_chainP_rep01_0p5ns/run_metadata.json ]; then
  echo '[skip] tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_6VQO_HMTEVVRHC_HLA-A0201_chainP_rep01_0p5ns already has run_metadata.json'
else
  python run_openmm_pilot.py --input inputs/6VQO_HMTEVVRHC_HLA-A0201_chainP_c3a323eb.pdb --outdir runs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_6VQO_HMTEVVRHC_HLA-A0201_chainP_rep01_0p5ns --mode explicit --platform CUDA --ns 0.5 --report-steps 50000 --minimize-iterations 1000 --timestep-fs 1 --temperature-k 300 --random-seed 10003 2>&1 | tee logs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_6VQO_HMTEVVRHC_HLA-A0201_chainP_rep01_0p5ns.openmm.log
fi
if [ -f runs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_6VQO_HMTEVVRHC_HLA-A0201_chainP_rep01_0p5ns/trajectory.dcd ]; then
  python analyze_openmm_pilot.py --top runs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_6VQO_HMTEVVRHC_HLA-A0201_chainP_rep01_0p5ns/prepared_start.pdb --traj runs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_6VQO_HMTEVVRHC_HLA-A0201_chainP_rep01_0p5ns/trajectory.dcd --peptide-chain P --mhc-chain A --tcr-chains 'E|D' --peptide-sequence HMTEVVRHC --out runs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_6VQO_HMTEVVRHC_HLA-A0201_chainP_rep01_0p5ns/analysis.json 2>&1 | tee logs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_6VQO_HMTEVVRHC_HLA-A0201_chainP_rep01_0p5ns.analysis.log
fi

echo '[job] tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_6VQO_HMTEVVRHC_HLA-A0201_chainQ_rep01_0p5ns'
if [ -f runs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_6VQO_HMTEVVRHC_HLA-A0201_chainQ_rep01_0p5ns/run_metadata.json ]; then
  echo '[skip] tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_6VQO_HMTEVVRHC_HLA-A0201_chainQ_rep01_0p5ns already has run_metadata.json'
else
  python run_openmm_pilot.py --input inputs/6VQO_HMTEVVRHC_HLA-A0201_chainQ_3a2c9b95.pdb --outdir runs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_6VQO_HMTEVVRHC_HLA-A0201_chainQ_rep01_0p5ns --mode explicit --platform CUDA --ns 0.5 --report-steps 50000 --minimize-iterations 1000 --timestep-fs 1 --temperature-k 300 --random-seed 10004 2>&1 | tee logs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_6VQO_HMTEVVRHC_HLA-A0201_chainQ_rep01_0p5ns.openmm.log
fi
if [ -f runs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_6VQO_HMTEVVRHC_HLA-A0201_chainQ_rep01_0p5ns/trajectory.dcd ]; then
  python analyze_openmm_pilot.py --top runs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_6VQO_HMTEVVRHC_HLA-A0201_chainQ_rep01_0p5ns/prepared_start.pdb --traj runs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_6VQO_HMTEVVRHC_HLA-A0201_chainQ_rep01_0p5ns/trajectory.dcd --peptide-chain Q --mhc-chain F --tcr-chains 'J|H' --peptide-sequence HMTEVVRHC --out runs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_6VQO_HMTEVVRHC_HLA-A0201_chainQ_rep01_0p5ns/analysis.json 2>&1 | tee logs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_6VQO_HMTEVVRHC_HLA-A0201_chainQ_rep01_0p5ns.analysis.log
fi

echo '[job] tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_7RM4_HMTEVVRHC_HLA-A0201_chainH_rep01_0p5ns'
if [ -f runs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_7RM4_HMTEVVRHC_HLA-A0201_chainH_rep01_0p5ns/run_metadata.json ]; then
  echo '[skip] tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_7RM4_HMTEVVRHC_HLA-A0201_chainH_rep01_0p5ns already has run_metadata.json'
else
  python run_openmm_pilot.py --input inputs/7RM4_HMTEVVRHC_HLA-A0201_chainH_ff5c4c67.pdb --outdir runs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_7RM4_HMTEVVRHC_HLA-A0201_chainH_rep01_0p5ns --mode explicit --platform CUDA --ns 0.5 --report-steps 50000 --minimize-iterations 1000 --timestep-fs 1 --temperature-k 300 --random-seed 10005 2>&1 | tee logs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_7RM4_HMTEVVRHC_HLA-A0201_chainH_rep01_0p5ns.openmm.log
fi
if [ -f runs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_7RM4_HMTEVVRHC_HLA-A0201_chainH_rep01_0p5ns/trajectory.dcd ]; then
  python analyze_openmm_pilot.py --top runs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_7RM4_HMTEVVRHC_HLA-A0201_chainH_rep01_0p5ns/prepared_start.pdb --traj runs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_7RM4_HMTEVVRHC_HLA-A0201_chainH_rep01_0p5ns/trajectory.dcd --peptide-chain H --mhc-chain F --tcr-chains 'J|I' --peptide-sequence HMTEVVRHC --out runs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_7RM4_HMTEVVRHC_HLA-A0201_chainH_rep01_0p5ns/analysis.json 2>&1 | tee logs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_7RM4_HMTEVVRHC_HLA-A0201_chainH_rep01_0p5ns.analysis.log
fi

echo '[job] tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_7RM4_HMTEVVRHC_HLA-A0201_chainC_rep01_0p5ns'
if [ -f runs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_7RM4_HMTEVVRHC_HLA-A0201_chainC_rep01_0p5ns/run_metadata.json ]; then
  echo '[skip] tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_7RM4_HMTEVVRHC_HLA-A0201_chainC_rep01_0p5ns already has run_metadata.json'
else
  python run_openmm_pilot.py --input inputs/7RM4_HMTEVVRHC_HLA-A0201_chainC_1d97444d.pdb --outdir runs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_7RM4_HMTEVVRHC_HLA-A0201_chainC_rep01_0p5ns --mode explicit --platform CUDA --ns 0.5 --report-steps 50000 --minimize-iterations 1000 --timestep-fs 1 --temperature-k 300 --random-seed 10006 2>&1 | tee logs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_7RM4_HMTEVVRHC_HLA-A0201_chainC_rep01_0p5ns.openmm.log
fi
if [ -f runs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_7RM4_HMTEVVRHC_HLA-A0201_chainC_rep01_0p5ns/trajectory.dcd ]; then
  python analyze_openmm_pilot.py --top runs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_7RM4_HMTEVVRHC_HLA-A0201_chainC_rep01_0p5ns/prepared_start.pdb --traj runs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_7RM4_HMTEVVRHC_HLA-A0201_chainC_rep01_0p5ns/trajectory.dcd --peptide-chain C --mhc-chain A --tcr-chains 'E|D' --peptide-sequence HMTEVVRHC --out runs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_7RM4_HMTEVVRHC_HLA-A0201_chainC_rep01_0p5ns/analysis.json 2>&1 | tee logs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_7RM4_HMTEVVRHC_HLA-A0201_chainC_rep01_0p5ns.analysis.log
fi

echo '[job] tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_7RM4_HMTEVVRHC_HLA-A0201_chainR_rep01_0p5ns'
if [ -f runs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_7RM4_HMTEVVRHC_HLA-A0201_chainR_rep01_0p5ns/run_metadata.json ]; then
  echo '[skip] tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_7RM4_HMTEVVRHC_HLA-A0201_chainR_rep01_0p5ns already has run_metadata.json'
else
  python run_openmm_pilot.py --input inputs/7RM4_HMTEVVRHC_HLA-A0201_chainR_1c90cd69.pdb --outdir runs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_7RM4_HMTEVVRHC_HLA-A0201_chainR_rep01_0p5ns --mode explicit --platform CUDA --ns 0.5 --report-steps 50000 --minimize-iterations 1000 --timestep-fs 1 --temperature-k 300 --random-seed 10007 2>&1 | tee logs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_7RM4_HMTEVVRHC_HLA-A0201_chainR_rep01_0p5ns.openmm.log
fi
if [ -f runs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_7RM4_HMTEVVRHC_HLA-A0201_chainR_rep01_0p5ns/trajectory.dcd ]; then
  python analyze_openmm_pilot.py --top runs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_7RM4_HMTEVVRHC_HLA-A0201_chainR_rep01_0p5ns/prepared_start.pdb --traj runs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_7RM4_HMTEVVRHC_HLA-A0201_chainR_rep01_0p5ns/trajectory.dcd --peptide-chain R --mhc-chain P --tcr-chains 'T|S' --peptide-sequence HMTEVVRHC --out runs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_7RM4_HMTEVVRHC_HLA-A0201_chainR_rep01_0p5ns/analysis.json 2>&1 | tee logs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_7RM4_HMTEVVRHC_HLA-A0201_chainR_rep01_0p5ns.analysis.log
fi

echo '[job] tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_7RM4_HMTEVVRHC_HLA-A0201_chainM_rep01_0p5ns'
if [ -f runs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_7RM4_HMTEVVRHC_HLA-A0201_chainM_rep01_0p5ns/run_metadata.json ]; then
  echo '[skip] tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_7RM4_HMTEVVRHC_HLA-A0201_chainM_rep01_0p5ns already has run_metadata.json'
else
  python run_openmm_pilot.py --input inputs/7RM4_HMTEVVRHC_HLA-A0201_chainM_1ad69c81.pdb --outdir runs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_7RM4_HMTEVVRHC_HLA-A0201_chainM_rep01_0p5ns --mode explicit --platform CUDA --ns 0.5 --report-steps 50000 --minimize-iterations 1000 --timestep-fs 1 --temperature-k 300 --random-seed 10008 2>&1 | tee logs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_7RM4_HMTEVVRHC_HLA-A0201_chainM_rep01_0p5ns.openmm.log
fi
if [ -f runs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_7RM4_HMTEVVRHC_HLA-A0201_chainM_rep01_0p5ns/trajectory.dcd ]; then
  python analyze_openmm_pilot.py --top runs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_7RM4_HMTEVVRHC_HLA-A0201_chainM_rep01_0p5ns/prepared_start.pdb --traj runs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_7RM4_HMTEVVRHC_HLA-A0201_chainM_rep01_0p5ns/trajectory.dcd --peptide-chain M --mhc-chain K --tcr-chains 'O|N' --peptide-sequence HMTEVVRHC --out runs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_7RM4_HMTEVVRHC_HLA-A0201_chainM_rep01_0p5ns/analysis.json 2>&1 | tee logs/tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_7RM4_HMTEVVRHC_HLA-A0201_chainM_rep01_0p5ns.analysis.log
fi

