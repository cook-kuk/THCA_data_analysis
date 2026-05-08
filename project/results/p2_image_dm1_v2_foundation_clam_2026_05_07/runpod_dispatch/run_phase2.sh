#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")/.."

[[ -z "${HF_TOKEN:-}" ]] && { echo "ERROR: export HF_TOKEN=hf_..."; exit 1; }

echo "===== Phase 2a — TCGA WSI manifest + download ====="
# manifest first (no download)
python scripts/phase2_gdc_download.py \
  --master ../../results/audit_2026_04_30/round9/r9_1_master_with_rai.tsv \
  --out_dir phase2_tcga_clam \
  --n_per_group 30 \
  --manifest_only

echo "Manifest built. Review phase2_tcga_clam/slide_manifest.tsv before download."
read -p "Proceed with download? [y/N] " yn
[[ "$yn" != "y" ]] && exit 0

# actual download
python scripts/phase2_gdc_download.py \
  --master ../../results/audit_2026_04_30/round9/r9_1_master_with_rai.tsv \
  --out_dir phase2_tcga_clam \
  --n_per_group 30

echo ""
echo "===== Phase 2b — CLAM training (UNI features) ====="
echo "Note: phase2_tile_extract_uni_embed.py is yet to be implemented."
echo "Expected workflow:"
echo "  1. CLAM tile extraction (20×, 256-px) per slide"
echo "  2. UNI embedding per tile -> .pt features"
echo "  3. CLAM gated-attention MIL training"
echo ""
echo "For now run separately:"
echo "  python scripts/phase2_clam_train_eval.py \\"
echo "    --feature_dir phase2_tcga_clam/features/ \\"
echo "    --manifest phase2_tcga_clam/slide_manifest.tsv \\"
echo "    --out_dir phase2_tcga_clam"
