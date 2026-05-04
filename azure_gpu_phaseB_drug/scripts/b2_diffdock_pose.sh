#!/bin/bash
# B2 — DiffDock pose generation (8 targets × top 50 compounds = 400 poses)
# Full run (all 8 × 50):  3-6h on A100 40GB · Cost: $20-40
# Taster   (TACSTD2 × 50): ~2h on RTX 4090 · Cost: ~$2
# Scope controlled by env: TARGETS (CSV; default = all 8), COMPOUNDS_TOP (default 50)
set -euo pipefail
cd /workspace/azure_gpu_phaseB_drug

OUT=results/b2_diffdock
mkdir -p $OUT/poses

# Honor scope env from run_phaseB_gpu.sh (--targets / --compounds-top)
if [ -n "${TARGETS:-}" ]; then
  IFS=',' read -ra TARGETS <<< "$TARGETS"
  echo "[B2] scope override: targets=${TARGETS[*]}"
else
  TARGETS=(TACSTD2 TMPRSS4 PLEKHA6 CYP1B1 LDLR GABRB2 B3GNT3 PTPRE)
fi
COMPOUNDS_TOP="${COMPOUNDS_TOP:-50}"
HEAD_LINES=$((COMPOUNDS_TOP + 1))
echo "[B2] compounds per target: $COMPOUNDS_TOP"

# Source structures: B1 refined (preferred) or v13 baseline
STRUCT_DIR=results/b1_structures/refined
if [ ! -d "$STRUCT_DIR" ] || [ -z "$(ls $STRUCT_DIR 2>/dev/null)" ]; then
  STRUCT_DIR=/workspace/v13_drug_discovery/structures
  echo "[B2] using v13 baseline structures (B1 not yet run)"
fi

# Compound source: per_target_top10 + ranked_candidates (top 50 per target)
COMPOUNDS=/workspace/v13_drug_discovery/prioritization/ranked_candidates.tsv

# DiffDock per target × top N SMILES (N = COMPOUNDS_TOP)
for t in "${TARGETS[@]}"; do
  echo "[B2] DiffDock $t × top $COMPOUNDS_TOP compounds"
  PDB=$STRUCT_DIR/$t.pdb
  mkdir -p $OUT/poses/$t

  # Extract top N SMILES for target
  awk -F'\t' -v t="$t" 'NR==1{for(i=1;i<=NF;i++)c[$i]=i} $c["target"]==t || NR==1' $COMPOUNDS | head -"$HEAD_LINES" > $OUT/poses/$t/compounds.tsv

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

# Aggregate confidence scores (only over the targets we actually ran)
TARGETS_CSV=$(IFS=,; echo "${TARGETS[*]}")
python3 - <<PYEOF
import os, pandas as pd
rows = []
targets = "${TARGETS_CSV}".split(',')
for t in targets:
    d = '${OUT}/poses/' + t
    if not os.path.isdir(d): continue
    for fn in sorted(os.listdir(d)):
        if fn.endswith('_confidence.csv'):
            df = pd.read_csv(os.path.join(d, fn))
            df['target'] = t
            rows.append(df)
if rows:
    pd.concat(rows).to_csv('${OUT}/pose_confidence.tsv', sep='\t', index=False)
    print(f'aggregated {sum(len(r) for r in rows)} poses across {len(set(r["target"].iloc[0] for r in rows))} targets')
else:
    print('[B2] warning: no confidence CSVs found')
PYEOF

echo "[B2] complete · poses for ${#TARGETS[@]} target(s) × $COMPOUNDS_TOP compounds + pose_confidence.tsv"
