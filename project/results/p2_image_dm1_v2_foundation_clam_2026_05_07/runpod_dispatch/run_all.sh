#!/bin/bash
set +e
cd "$(dirname "$0")/.."
echo "==============================================="
echo "RUN_ALL — Phase 1 + 2 sequential"
echo "Start: $(date)"
echo "==============================================="
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
bash runpod_dispatch/run_phase1.sh 2>&1 | tee phase1_gse250521/run_phase1.log
P1_VERDICT=$(grep -oE "VERDICT.*?(PASS|MARGINAL|FAIL)" phase1_gse250521/PHASE1_REPORT.md 2>/dev/null | head -1 | grep -oE "PASS|MARGINAL|FAIL")
echo "Phase 1 verdict: $P1_VERDICT"
if [[ "$P1_VERDICT" == "FAIL" ]]; then
    echo "[KILL SWITCH] Phase 1 FAIL — stopping"
    exit 0
fi
echo ""; echo "===== PHASE 2 ====="
python scripts/phase2_gdc_download.py \
  --master ../../audit_2026_04_30/round9/r9_1_master_with_rai.tsv \
  --out_dir phase2_tcga_clam --n_per_group 30
python scripts/phase2_tile_extract_uni_embed.py \
  --wsi_dir phase2_tcga_clam/wsi --manifest phase2_tcga_clam/slide_manifest.tsv \
  --out_dir phase2_tcga_clam --batch_size 64
python scripts/phase2_clam_train_eval.py \
  --feature_dir phase2_tcga_clam/features --manifest phase2_tcga_clam/slide_manifest.tsv \
  --out_dir phase2_tcga_clam --n_folds 5
P2_VERDICT=$(grep -oE "VERDICT.*?(PASS|MARGINAL|FAIL)" phase2_tcga_clam/PHASE2_REPORT.md 2>/dev/null | head -1 | grep -oE "PASS|MARGINAL|FAIL")
echo "Phase 1: $P1_VERDICT  Phase 2: $P2_VERDICT"
echo "End: $(date)"
