#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"
mkdir -p runs logs

echo 'immediate ready diverse MD batch'

echo '[job] tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_mutant_TCR-pMHC_6UON_GADGVGKSAL_HLA-C0802_chainF_rep01_10p0ns'
if [ -f runs/tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_mutant_TCR-pMHC_6UON_GADGVGKSAL_HLA-C0802_chainF_rep01_10p0ns/run_metadata.json ]; then
  echo '[skip] tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_mutant_TCR-pMHC_6UON_GADGVGKSAL_HLA-C0802_chainF_rep01_10p0ns already has run_metadata.json'
else
  python run_openmm_pilot.py --input inputs/6UON_GADGVGKSAL_HLA-C0802_chainF_81df2e11.pdb --outdir runs/tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_mutant_TCR-pMHC_6UON_GADGVGKSAL_HLA-C0802_chainF_rep01_10p0ns --mode explicit --platform CUDA --ns 10 --report-steps 50000 --minimize-iterations 1000 --timestep-fs 2 --temperature-k 300 --random-seed 21011 2>&1 | tee logs/tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_mutant_TCR-pMHC_6UON_GADGVGKSAL_HLA-C0802_chainF_rep01_10p0ns.openmm.log
fi
if [ -f runs/tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_mutant_TCR-pMHC_6UON_GADGVGKSAL_HLA-C0802_chainF_rep01_10p0ns/trajectory.dcd ]; then
  python analyze_openmm_pilot.py --top runs/tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_mutant_TCR-pMHC_6UON_GADGVGKSAL_HLA-C0802_chainF_rep01_10p0ns/prepared_start.pdb --traj runs/tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_mutant_TCR-pMHC_6UON_GADGVGKSAL_HLA-C0802_chainF_rep01_10p0ns/trajectory.dcd --peptide-chain F --mhc-chain D --tcr-chains 'G|H' --peptide-sequence GADGVGKSAL --out runs/tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_mutant_TCR-pMHC_6UON_GADGVGKSAL_HLA-C0802_chainF_rep01_10p0ns/analysis.json 2>&1 | tee logs/tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_mutant_TCR-pMHC_6UON_GADGVGKSAL_HLA-C0802_chainF_rep01_10p0ns.analysis.log
fi

echo '[job] tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_mutant_TCR-pMHC_6UON_GADGVGKSAL_HLA-C0802_chainF_rep02_10p0ns'
if [ -f runs/tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_mutant_TCR-pMHC_6UON_GADGVGKSAL_HLA-C0802_chainF_rep02_10p0ns/run_metadata.json ]; then
  echo '[skip] tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_mutant_TCR-pMHC_6UON_GADGVGKSAL_HLA-C0802_chainF_rep02_10p0ns already has run_metadata.json'
else
  python run_openmm_pilot.py --input inputs/6UON_GADGVGKSAL_HLA-C0802_chainF_81df2e11.pdb --outdir runs/tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_mutant_TCR-pMHC_6UON_GADGVGKSAL_HLA-C0802_chainF_rep02_10p0ns --mode explicit --platform CUDA --ns 10 --report-steps 50000 --minimize-iterations 1000 --timestep-fs 2 --temperature-k 300 --random-seed 21012 2>&1 | tee logs/tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_mutant_TCR-pMHC_6UON_GADGVGKSAL_HLA-C0802_chainF_rep02_10p0ns.openmm.log
fi
if [ -f runs/tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_mutant_TCR-pMHC_6UON_GADGVGKSAL_HLA-C0802_chainF_rep02_10p0ns/trajectory.dcd ]; then
  python analyze_openmm_pilot.py --top runs/tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_mutant_TCR-pMHC_6UON_GADGVGKSAL_HLA-C0802_chainF_rep02_10p0ns/prepared_start.pdb --traj runs/tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_mutant_TCR-pMHC_6UON_GADGVGKSAL_HLA-C0802_chainF_rep02_10p0ns/trajectory.dcd --peptide-chain F --mhc-chain D --tcr-chains 'G|H' --peptide-sequence GADGVGKSAL --out runs/tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_mutant_TCR-pMHC_6UON_GADGVGKSAL_HLA-C0802_chainF_rep02_10p0ns/analysis.json 2>&1 | tee logs/tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_mutant_TCR-pMHC_6UON_GADGVGKSAL_HLA-C0802_chainF_rep02_10p0ns.analysis.log
fi

echo '[job] tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_same_or_similar_hla_positive_control_pMHC_6ULR_GADGVGKSA_HLA-C0802_chainC_rep01_10p0ns'
if [ -f runs/tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_same_or_similar_hla_positive_control_pMHC_6ULR_GADGVGKSA_HLA-C0802_chainC_rep01_10p0ns/run_metadata.json ]; then
  echo '[skip] tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_same_or_similar_hla_positive_control_pMHC_6ULR_GADGVGKSA_HLA-C0802_chainC_rep01_10p0ns already has run_metadata.json'
else
  python run_openmm_pilot.py --input inputs/6ULR_GADGVGKSA_HLA-C0802_chainC_d543151a.pdb --outdir runs/tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_same_or_similar_hla_positive_control_pMHC_6ULR_GADGVGKSA_HLA-C0802_chainC_rep01_10p0ns --mode explicit --platform CUDA --ns 10 --report-steps 50000 --minimize-iterations 1000 --timestep-fs 2 --temperature-k 300 --random-seed 21021 2>&1 | tee logs/tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_same_or_similar_hla_positive_control_pMHC_6ULR_GADGVGKSA_HLA-C0802_chainC_rep01_10p0ns.openmm.log
fi
if [ -f runs/tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_same_or_similar_hla_positive_control_pMHC_6ULR_GADGVGKSA_HLA-C0802_chainC_rep01_10p0ns/trajectory.dcd ]; then
  python analyze_openmm_pilot.py --top runs/tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_same_or_similar_hla_positive_control_pMHC_6ULR_GADGVGKSA_HLA-C0802_chainC_rep01_10p0ns/prepared_start.pdb --traj runs/tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_same_or_similar_hla_positive_control_pMHC_6ULR_GADGVGKSA_HLA-C0802_chainC_rep01_10p0ns/trajectory.dcd --peptide-chain C --mhc-chain A --tcr-chains '' --peptide-sequence GADGVGKSA --out runs/tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_same_or_similar_hla_positive_control_pMHC_6ULR_GADGVGKSA_HLA-C0802_chainC_rep01_10p0ns/analysis.json 2>&1 | tee logs/tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_same_or_similar_hla_positive_control_pMHC_6ULR_GADGVGKSA_HLA-C0802_chainC_rep01_10p0ns.analysis.log
fi

echo '[job] tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_same_or_similar_hla_positive_control_pMHC_6ULR_GADGVGKSA_HLA-C0802_chainC_rep02_10p0ns'
if [ -f runs/tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_same_or_similar_hla_positive_control_pMHC_6ULR_GADGVGKSA_HLA-C0802_chainC_rep02_10p0ns/run_metadata.json ]; then
  echo '[skip] tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_same_or_similar_hla_positive_control_pMHC_6ULR_GADGVGKSA_HLA-C0802_chainC_rep02_10p0ns already has run_metadata.json'
else
  python run_openmm_pilot.py --input inputs/6ULR_GADGVGKSA_HLA-C0802_chainC_d543151a.pdb --outdir runs/tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_same_or_similar_hla_positive_control_pMHC_6ULR_GADGVGKSA_HLA-C0802_chainC_rep02_10p0ns --mode explicit --platform CUDA --ns 10 --report-steps 50000 --minimize-iterations 1000 --timestep-fs 2 --temperature-k 300 --random-seed 21022 2>&1 | tee logs/tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_same_or_similar_hla_positive_control_pMHC_6ULR_GADGVGKSA_HLA-C0802_chainC_rep02_10p0ns.openmm.log
fi
if [ -f runs/tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_same_or_similar_hla_positive_control_pMHC_6ULR_GADGVGKSA_HLA-C0802_chainC_rep02_10p0ns/trajectory.dcd ]; then
  python analyze_openmm_pilot.py --top runs/tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_same_or_similar_hla_positive_control_pMHC_6ULR_GADGVGKSA_HLA-C0802_chainC_rep02_10p0ns/prepared_start.pdb --traj runs/tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_same_or_similar_hla_positive_control_pMHC_6ULR_GADGVGKSA_HLA-C0802_chainC_rep02_10p0ns/trajectory.dcd --peptide-chain C --mhc-chain A --tcr-chains '' --peptide-sequence GADGVGKSA --out runs/tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_same_or_similar_hla_positive_control_pMHC_6ULR_GADGVGKSA_HLA-C0802_chainC_rep02_10p0ns/analysis.json 2>&1 | tee logs/tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_same_or_similar_hla_positive_control_pMHC_6ULR_GADGVGKSA_HLA-C0802_chainC_rep02_10p0ns.analysis.log
fi

echo '[job] tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_same_or_similar_hla_positive_control_pMHC_6ULR_GADGVGKSA_HLA-C0802_chainC_rep03_10p0ns'
if [ -f runs/tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_same_or_similar_hla_positive_control_pMHC_6ULR_GADGVGKSA_HLA-C0802_chainC_rep03_10p0ns/run_metadata.json ]; then
  echo '[skip] tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_same_or_similar_hla_positive_control_pMHC_6ULR_GADGVGKSA_HLA-C0802_chainC_rep03_10p0ns already has run_metadata.json'
else
  python run_openmm_pilot.py --input inputs/6ULR_GADGVGKSA_HLA-C0802_chainC_d543151a.pdb --outdir runs/tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_same_or_similar_hla_positive_control_pMHC_6ULR_GADGVGKSA_HLA-C0802_chainC_rep03_10p0ns --mode explicit --platform CUDA --ns 10 --report-steps 50000 --minimize-iterations 1000 --timestep-fs 2 --temperature-k 300 --random-seed 21023 2>&1 | tee logs/tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_same_or_similar_hla_positive_control_pMHC_6ULR_GADGVGKSA_HLA-C0802_chainC_rep03_10p0ns.openmm.log
fi
if [ -f runs/tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_same_or_similar_hla_positive_control_pMHC_6ULR_GADGVGKSA_HLA-C0802_chainC_rep03_10p0ns/trajectory.dcd ]; then
  python analyze_openmm_pilot.py --top runs/tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_same_or_similar_hla_positive_control_pMHC_6ULR_GADGVGKSA_HLA-C0802_chainC_rep03_10p0ns/prepared_start.pdb --traj runs/tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_same_or_similar_hla_positive_control_pMHC_6ULR_GADGVGKSA_HLA-C0802_chainC_rep03_10p0ns/trajectory.dcd --peptide-chain C --mhc-chain A --tcr-chains '' --peptide-sequence GADGVGKSA --out runs/tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_same_or_similar_hla_positive_control_pMHC_6ULR_GADGVGKSA_HLA-C0802_chainC_rep03_10p0ns/analysis.json 2>&1 | tee logs/tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_same_or_similar_hla_positive_control_pMHC_6ULR_GADGVGKSA_HLA-C0802_chainC_rep03_10p0ns.analysis.log
fi

