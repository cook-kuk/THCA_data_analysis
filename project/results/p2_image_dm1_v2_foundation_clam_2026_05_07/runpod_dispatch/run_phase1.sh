#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")/.."

[[ -z "${HF_TOKEN:-}" ]] && { echo "ERROR: export HF_TOKEN=hf_..."; exit 1; }
huggingface-cli login --token "$HF_TOKEN" >/dev/null 2>&1 || true

echo "===== Phase 1a — UNI embeddings on GSE250521 ====="
python scripts/phase1_uni_embed.py \
  --tiles_root ../../data/processed/GSE250521/tiles \
  --out_dir phase1_gse250521 \
  --tile_size 224 \
  --batch_size 128

echo ""
echo "===== Phase 1b — UNI ↔ DM1 correlation + spatial overlay ====="
python scripts/phase1_overlay_correlation.py \
  --embed_npz phase1_gse250521/uni_embeddings_size224.npz \
  --meta_tsv phase1_gse250521/uni_embed_metadata_size224.tsv \
  --scores_tsv ../../results/01_spatial_score/all_spots_scored.tsv.gz \
  --out_dir phase1_gse250521

echo ""
echo "===== Phase 1 done ====="
echo "Verdict:"
cat phase1_gse250521/PHASE1_REPORT.md
