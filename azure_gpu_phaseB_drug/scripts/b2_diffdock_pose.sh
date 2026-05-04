#!/bin/bash
# B2 — DiffDock pose generation (8 targets × top 50 compounds = 400 poses)
# Time: 3-6h on A100 40GB · Cost: $20-40
set -euo pipefail
cd /workspace/azure_gpu_phaseB_drug

OUT=results/b2_diffdock
mkdir -p $OUT/poses
TARGETS=(TACSTD2 TMPRSS4 PLEKHA6 CYP1B1 LDLR GABRB2 B3GNT3 PTPRE)

# Source structures: B1 refined (preferred) or v13 baseline
STRUCT_DIR=results/b1_structures/refined
if [ ! -d "$STRUCT_DIR" ] || [ -z "$(ls $STRUCT_DIR 2>/dev/null)" ]; then
  STRUCT_DIR=/workspace/v13_drug_discovery/structures
  echo "[B2] using v13 baseline structures (B1 not yet run)"
fi

# Compound source: per_target_top10 + ranked_candidates (top 50 per target)
COMPOUNDS=/workspace/v13_drug_discovery/prioritization/ranked_candidates.tsv

# DiffDock per target × top 50 SMILES
for t in "${TARGETS[@]}"; do
  echo "[B2] DiffDock $t × top 50 compounds"
  PDB=$STRUCT_DIR/$t.pdb
  mkdir -p $OUT/poses/$t

  # Extract top 50 SMILES for target
  awk -F'\t' -v t="$t" 'NR==1{for(i=1;i<=NF;i++)c[$i]=i} $c["target"]==t || NR==1' $COMPOUNDS | head -51 > $OUT/poses/$t/compounds.tsv

  # DiffDock invocation
  python3 /workspace/DiffDock/inference.py \
    --protein_path $PDB \
    --ligand_csv $OUT/poses/$t/compounds.tsv \
    --out_dir $OUT/poses/$t \
    --inference_steps 20 \
    --samples_per_complex 10 \
    --batch_size 10 \
    --actual_steps 18 \
    --no_final_step_noise
done

# Aggregate confidence scores
python3 -c "
import os, pandas as pd
rows = []
for t in ['TACSTD2','TMPRSS4','PLEKHA6','CYP1B1','LDLR','GABRB2','B3GNT3','PTPRE']:
    d = '$OUT/poses/' + t
    if not os.path.isdir(d): continue
    for fn in sorted(os.listdir(d)):
        if fn.endswith('_confidence.csv'):
            df = pd.read_csv(os.path.join(d, fn))
            df['target'] = t
            rows.append(df)
if rows:
    pd.concat(rows).to_csv('$OUT/pose_confidence.tsv', sep='\t', index=False)
    print(f'aggregated {sum(len(r) for r in rows)} poses')
"

echo "[B2] complete · ~400 poses + pose_confidence.tsv"
