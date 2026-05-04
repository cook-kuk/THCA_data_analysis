#!/bin/bash
# Phase B Drug Discovery — Master orchestrator (RunPod)
# Author: Seungho Cook · 2026-05-04
# Usage: bash run_phaseB_gpu.sh --module b1 [--pod $POD_ID] [--dry-run]
set -euo pipefail

MODULE=""
POD_ID=""
DRY_RUN=false
COST_CAP=500   # USD per K6 kill switch
TIER_A_TOP=10           # B10 default; --tier-a-top N to scope down
TARGETS=""              # B2 default = all 8; --targets TACSTD2 to scope to one
COMPOUNDS_TOP=""        # B2/B5 per-target compound cap; default = script-internal

# Parse args
while [[ $# -gt 0 ]]; do
  case $1 in
    --module) MODULE="$2"; shift 2 ;;
    --pod) POD_ID="$2"; shift 2 ;;
    --dry-run) DRY_RUN=true; shift ;;
    --cost-cap) COST_CAP="$2"; shift 2 ;;
    --tier-a-top) TIER_A_TOP="$2"; shift 2 ;;
    --targets) TARGETS="$2"; shift 2 ;;
    --compounds-top) COMPOUNDS_TOP="$2"; shift 2 ;;
    *) echo "Unknown: $1"; exit 1 ;;
  esac
done
export TIER_A_TOP TARGETS COMPOUNDS_TOP

if [ -z "$MODULE" ]; then
  echo "Available modules:"
  echo "  b1   Structure refinement (AlphaFold3 / ESMFold) · 4-8h · \$30-60"
  echo "  b2   DiffDock pose generation (8 targets × 50 cmpds) · 3-6h · \$20-40"
  echo "       Override with --targets TACSTD2 (1 target) → ~2h · ~\$2 for taster"
  echo "       Add --compounds-top N to scope per-target compounds (default 50)"
  echo "  b3   GNINA rescoring · 2-4h · \$10-20"
  echo "  b4   KDeep / DeepPurpose binding affinity · 4-6h · \$25-40"
  echo "  b5   ChemBERTa scaffold hopping · 6-10h · \$40-70"
  echo "  b6   PROTAC-DB matching (intracellular targets) · 2-3h · \$10-15"
  echo "  b7   ADC linker / payload optimization · 4-6h · \$25-40"
  echo "  b8   Off-target profiling · 2-3h · \$10-15"
  echo "  b9   Foundation-model ADMET · 6-8h · \$40-60"
  echo "  b10  FEP / MD relative ΔΔG (Tier-A top N, default 10) · 12-24h · \$200-400"
  echo "       Override with --tier-a-top N (e.g. N=3 → ~10h · ~\$5 for taster)"
  echo "  all  Run all 10 sequentially (~ 45-95h, \$610-1160)"
  echo "  cheap5  Run B2+B3+B6+B8+B9 (~ 18-26h, \$95-155)"
  echo ""
  echo "Hard gates G1-G6 must pass; kill switches K1-K6 monitored."
  echo "G6: explicit user authorization required (this script requires --module)."
  exit 0
fi

# G1: cost check (basic)
echo "[G1] Cost cap: \$${COST_CAP} (K6 kill if exceeded)"

# Run module
case $MODULE in
  b1) bash scripts/b1_structure_refinement.sh ;;
  b2) bash scripts/b2_diffdock_pose.sh ;;
  b3) bash scripts/b3_gnina_rescore.sh ;;
  b4) bash scripts/b4_kdeep_affinity.sh ;;
  b5) bash scripts/b5_chemberta_scaffold.sh ;;
  b6) bash scripts/b6_protac_db_match.sh ;;
  b7) bash scripts/b7_adc_optimization.sh ;;
  b8) bash scripts/b8_off_target.sh ;;
  b9) bash scripts/b9_admet_foundation.sh ;;
  b10) bash scripts/b10_fep_md.sh ;;
  all)
    for m in b1 b2 b3 b4 b5 b6 b7 b8 b9 b10; do
      echo "===== Running $m ====="
      bash scripts/${m}_*.sh
    done
    ;;
  cheap5)
    for m in b2 b3 b6 b8 b9; do
      echo "===== Running $m ====="
      bash scripts/${m}_*.sh
    done
    ;;
  *) echo "Unknown module: $MODULE"; exit 1 ;;
esac

echo "===== Phase B module $MODULE complete ====="
echo "Output: /workspace/azure_gpu_phaseB_drug/results/${MODULE}_*/"
echo "SCP back to local with:"
echo "  scp -r -P \$POD_PORT root@\$POD_IP:/workspace/azure_gpu_phaseB_drug/results/${MODULE}_*/ \\"
echo "      /home/seungho/personal/THCA_data_analysis/project/results/v14_drug_phaseB/"
