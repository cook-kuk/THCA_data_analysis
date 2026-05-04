#!/bin/bash
# B3 — GNINA rescoring (CNN-based binding affinity)
# Time: 2-4h on A40 · Cost: $10-20
set -euo pipefail
cd /workspace/azure_gpu_phaseB_drug

OUT=results/b3_gnina
mkdir -p $OUT
TARGETS=(TACSTD2 TMPRSS4 PLEKHA6 CYP1B1 LDLR GABRB2 B3GNT3 PTPRE)

POSES=results/b2_diffdock/poses
> $OUT/rescore.tsv
echo -e "target\tcompound\tpose_id\tvina_score\tcnn_score\tcnn_affinity" > $OUT/rescore.tsv

for t in "${TARGETS[@]}"; do
  PDB=results/b1_structures/refined/$t.pdb
  [ -f "$PDB" ] || PDB=/workspace/v13_drug_discovery/structures/$t.pdb
  for pose in $POSES/$t/*/rank1.sdf; do
    [ -f "$pose" ] || continue
    cmpd=$(basename $(dirname $pose))
    out=$(gnina --score_only -r $PDB -l $pose --cnn_scoring rescore 2>/dev/null | grep -E 'Affinity|CNNscore|CNNaffinity' | awk '{print $NF}' | tr '\n' '\t')
    echo -e "$t\t$cmpd\trank1\t$out" >> $OUT/rescore.tsv
  done
done

echo "[B3] complete · GNINA rescoring on $(wc -l < $OUT/rescore.tsv) poses"
