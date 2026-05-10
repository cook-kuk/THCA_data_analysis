#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"
mkdir -p runs logs

echo 'HMTEVVRHC follow-up batch after current 10ns completion/sync'

echo '[job] tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_mutant_TCR-pMHC_6VRN_HMTEVVRHC_HLA-A0201_chainP_rep01_10p0ns'
if [ -f runs/tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_mutant_TCR-pMHC_6VRN_HMTEVVRHC_HLA-A0201_chainP_rep01_10p0ns/run_metadata.json ]; then
  echo '[skip] tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_mutant_TCR-pMHC_6VRN_HMTEVVRHC_HLA-A0201_chainP_rep01_10p0ns already has run_metadata.json'
else
  python run_openmm_pilot.py --input inputs/6VRN_HMTEVVRHC_HLA-A0201_chainP_f6de3374.pdb --outdir runs/tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_mutant_TCR-pMHC_6VRN_HMTEVVRHC_HLA-A0201_chainP_rep01_10p0ns --mode explicit --platform CUDA --ns 10 --report-steps 50000 --minimize-iterations 1000 --timestep-fs 2 --temperature-k 300 --random-seed 21031 2>&1 | tee logs/tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_mutant_TCR-pMHC_6VRN_HMTEVVRHC_HLA-A0201_chainP_rep01_10p0ns.openmm.log
fi
if [ -f runs/tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_mutant_TCR-pMHC_6VRN_HMTEVVRHC_HLA-A0201_chainP_rep01_10p0ns/trajectory.dcd ]; then
  python analyze_openmm_pilot.py --top runs/tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_mutant_TCR-pMHC_6VRN_HMTEVVRHC_HLA-A0201_chainP_rep01_10p0ns/prepared_start.pdb --traj runs/tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_mutant_TCR-pMHC_6VRN_HMTEVVRHC_HLA-A0201_chainP_rep01_10p0ns/trajectory.dcd --peptide-chain P --mhc-chain A --tcr-chains 'D|E' --peptide-sequence HMTEVVRHC --out runs/tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_mutant_TCR-pMHC_6VRN_HMTEVVRHC_HLA-A0201_chainP_rep01_10p0ns/analysis.json 2>&1 | tee logs/tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_mutant_TCR-pMHC_6VRN_HMTEVVRHC_HLA-A0201_chainP_rep01_10p0ns.analysis.log
fi

echo '[job] tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_mutant_TCR-pMHC_6VRN_HMTEVVRHC_HLA-A0201_chainP_rep02_10p0ns'
if [ -f runs/tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_mutant_TCR-pMHC_6VRN_HMTEVVRHC_HLA-A0201_chainP_rep02_10p0ns/run_metadata.json ]; then
  echo '[skip] tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_mutant_TCR-pMHC_6VRN_HMTEVVRHC_HLA-A0201_chainP_rep02_10p0ns already has run_metadata.json'
else
  python run_openmm_pilot.py --input inputs/6VRN_HMTEVVRHC_HLA-A0201_chainP_f6de3374.pdb --outdir runs/tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_mutant_TCR-pMHC_6VRN_HMTEVVRHC_HLA-A0201_chainP_rep02_10p0ns --mode explicit --platform CUDA --ns 10 --report-steps 50000 --minimize-iterations 1000 --timestep-fs 2 --temperature-k 300 --random-seed 21032 2>&1 | tee logs/tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_mutant_TCR-pMHC_6VRN_HMTEVVRHC_HLA-A0201_chainP_rep02_10p0ns.openmm.log
fi
if [ -f runs/tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_mutant_TCR-pMHC_6VRN_HMTEVVRHC_HLA-A0201_chainP_rep02_10p0ns/trajectory.dcd ]; then
  python analyze_openmm_pilot.py --top runs/tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_mutant_TCR-pMHC_6VRN_HMTEVVRHC_HLA-A0201_chainP_rep02_10p0ns/prepared_start.pdb --traj runs/tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_mutant_TCR-pMHC_6VRN_HMTEVVRHC_HLA-A0201_chainP_rep02_10p0ns/trajectory.dcd --peptide-chain P --mhc-chain A --tcr-chains 'D|E' --peptide-sequence HMTEVVRHC --out runs/tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_mutant_TCR-pMHC_6VRN_HMTEVVRHC_HLA-A0201_chainP_rep02_10p0ns/analysis.json 2>&1 | tee logs/tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_mutant_TCR-pMHC_6VRN_HMTEVVRHC_HLA-A0201_chainP_rep02_10p0ns.analysis.log
fi

echo '[job] tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_same_or_similar_hla_positive_control_TCR-pMHC_5NQK_ELAGIGILTV_HLA-A0201_chainP_rep01_10p0ns'
if [ -f runs/tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_same_or_similar_hla_positive_control_TCR-pMHC_5NQK_ELAGIGILTV_HLA-A0201_chainP_rep01_10p0ns/run_metadata.json ]; then
  echo '[skip] tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_same_or_similar_hla_positive_control_TCR-pMHC_5NQK_ELAGIGILTV_HLA-A0201_chainP_rep01_10p0ns already has run_metadata.json'
else
  python run_openmm_pilot.py --input inputs/5NQK_ELAGIGILTV_HLA-A0201_chainP_66c95c14.pdb --outdir runs/tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_same_or_similar_hla_positive_control_TCR-pMHC_5NQK_ELAGIGILTV_HLA-A0201_chainP_rep01_10p0ns --mode explicit --platform CUDA --ns 10 --report-steps 50000 --minimize-iterations 1000 --timestep-fs 2 --temperature-k 300 --random-seed 21041 2>&1 | tee logs/tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_same_or_similar_hla_positive_control_TCR-pMHC_5NQK_ELAGIGILTV_HLA-A0201_chainP_rep01_10p0ns.openmm.log
fi
if [ -f runs/tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_same_or_similar_hla_positive_control_TCR-pMHC_5NQK_ELAGIGILTV_HLA-A0201_chainP_rep01_10p0ns/trajectory.dcd ]; then
  python analyze_openmm_pilot.py --top runs/tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_same_or_similar_hla_positive_control_TCR-pMHC_5NQK_ELAGIGILTV_HLA-A0201_chainP_rep01_10p0ns/prepared_start.pdb --traj runs/tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_same_or_similar_hla_positive_control_TCR-pMHC_5NQK_ELAGIGILTV_HLA-A0201_chainP_rep01_10p0ns/trajectory.dcd --peptide-chain P --mhc-chain H --tcr-chains 'A|B' --peptide-sequence ELAGIGILTV --out runs/tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_same_or_similar_hla_positive_control_TCR-pMHC_5NQK_ELAGIGILTV_HLA-A0201_chainP_rep01_10p0ns/analysis.json 2>&1 | tee logs/tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_same_or_similar_hla_positive_control_TCR-pMHC_5NQK_ELAGIGILTV_HLA-A0201_chainP_rep01_10p0ns.analysis.log
fi

echo '[job] tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_same_or_similar_hla_positive_control_TCR-pMHC_5NQK_ELAGIGILTV_HLA-A0201_chainP_rep02_10p0ns'
if [ -f runs/tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_same_or_similar_hla_positive_control_TCR-pMHC_5NQK_ELAGIGILTV_HLA-A0201_chainP_rep02_10p0ns/run_metadata.json ]; then
  echo '[skip] tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_same_or_similar_hla_positive_control_TCR-pMHC_5NQK_ELAGIGILTV_HLA-A0201_chainP_rep02_10p0ns already has run_metadata.json'
else
  python run_openmm_pilot.py --input inputs/5NQK_ELAGIGILTV_HLA-A0201_chainP_66c95c14.pdb --outdir runs/tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_same_or_similar_hla_positive_control_TCR-pMHC_5NQK_ELAGIGILTV_HLA-A0201_chainP_rep02_10p0ns --mode explicit --platform CUDA --ns 10 --report-steps 50000 --minimize-iterations 1000 --timestep-fs 2 --temperature-k 300 --random-seed 21042 2>&1 | tee logs/tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_same_or_similar_hla_positive_control_TCR-pMHC_5NQK_ELAGIGILTV_HLA-A0201_chainP_rep02_10p0ns.openmm.log
fi
if [ -f runs/tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_same_or_similar_hla_positive_control_TCR-pMHC_5NQK_ELAGIGILTV_HLA-A0201_chainP_rep02_10p0ns/trajectory.dcd ]; then
  python analyze_openmm_pilot.py --top runs/tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_same_or_similar_hla_positive_control_TCR-pMHC_5NQK_ELAGIGILTV_HLA-A0201_chainP_rep02_10p0ns/prepared_start.pdb --traj runs/tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_same_or_similar_hla_positive_control_TCR-pMHC_5NQK_ELAGIGILTV_HLA-A0201_chainP_rep02_10p0ns/trajectory.dcd --peptide-chain P --mhc-chain H --tcr-chains 'A|B' --peptide-sequence ELAGIGILTV --out runs/tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_same_or_similar_hla_positive_control_TCR-pMHC_5NQK_ELAGIGILTV_HLA-A0201_chainP_rep02_10p0ns/analysis.json 2>&1 | tee logs/tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_same_or_similar_hla_positive_control_TCR-pMHC_5NQK_ELAGIGILTV_HLA-A0201_chainP_rep02_10p0ns.analysis.log
fi

echo '[job] tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_same_or_similar_hla_positive_control_TCR-pMHC_5NQK_ELAGIGILTV_HLA-A0201_chainP_rep03_10p0ns'
if [ -f runs/tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_same_or_similar_hla_positive_control_TCR-pMHC_5NQK_ELAGIGILTV_HLA-A0201_chainP_rep03_10p0ns/run_metadata.json ]; then
  echo '[skip] tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_same_or_similar_hla_positive_control_TCR-pMHC_5NQK_ELAGIGILTV_HLA-A0201_chainP_rep03_10p0ns already has run_metadata.json'
else
  python run_openmm_pilot.py --input inputs/5NQK_ELAGIGILTV_HLA-A0201_chainP_66c95c14.pdb --outdir runs/tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_same_or_similar_hla_positive_control_TCR-pMHC_5NQK_ELAGIGILTV_HLA-A0201_chainP_rep03_10p0ns --mode explicit --platform CUDA --ns 10 --report-steps 50000 --minimize-iterations 1000 --timestep-fs 2 --temperature-k 300 --random-seed 21043 2>&1 | tee logs/tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_same_or_similar_hla_positive_control_TCR-pMHC_5NQK_ELAGIGILTV_HLA-A0201_chainP_rep03_10p0ns.openmm.log
fi
if [ -f runs/tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_same_or_similar_hla_positive_control_TCR-pMHC_5NQK_ELAGIGILTV_HLA-A0201_chainP_rep03_10p0ns/trajectory.dcd ]; then
  python analyze_openmm_pilot.py --top runs/tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_same_or_similar_hla_positive_control_TCR-pMHC_5NQK_ELAGIGILTV_HLA-A0201_chainP_rep03_10p0ns/prepared_start.pdb --traj runs/tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_same_or_similar_hla_positive_control_TCR-pMHC_5NQK_ELAGIGILTV_HLA-A0201_chainP_rep03_10p0ns/trajectory.dcd --peptide-chain P --mhc-chain H --tcr-chains 'A|B' --peptide-sequence ELAGIGILTV --out runs/tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_same_or_similar_hla_positive_control_TCR-pMHC_5NQK_ELAGIGILTV_HLA-A0201_chainP_rep03_10p0ns/analysis.json 2>&1 | tee logs/tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_same_or_similar_hla_positive_control_TCR-pMHC_5NQK_ELAGIGILTV_HLA-A0201_chainP_rep03_10p0ns.analysis.log
fi

