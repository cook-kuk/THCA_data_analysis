#!/usr/bin/env bash
# v6 Wave 2 master runner — sequential 10 → 60.
set -uo pipefail
INPUT=${INPUT:-/data/thca/scrna/processed/classical_baseline.h5ad}
RES=${RES:-/opt/thyroid-dash/project/results/v6_scrna}
LOGDIR=${LOGDIR:-/opt/thyroid-dash/project/logs}
CONT=${CONTINUE_ON_ERROR:-0}
mkdir -p "$LOGDIR"
T0=$(date +%s)

run() {
  local script="$1"; shift
  local name; name=$(basename "$script" .py)
  local log="$LOGDIR/v6_${name}.log"
  echo "=== $name begins $(date -u) ==="
  if python3 -u "$script" --input "$INPUT" --output_dir "$RES" "$@" 2>&1 | tee "$log"; then
    echo "=== $name ok ($((($(date +%s) - T0))) s elapsed total) ==="
  else
    echo "!!! $name FAILED — see $log"
    [[ "$CONT" == "1" ]] || exit 1
  fi
}

cd "$(dirname "$0")"
run 10_scgpt_embedding.py
run 20_geneformer.py
run 30_gears_perturb.py
run 40_cellrank.py
run 50_celloracle_tf.py
run 60_vega.py

echo
echo "=== v6 Wave 2 STAGE 99 ==="
echo "Stamps:"
ls -la "$RES"/.stamp_* 2>/dev/null
echo
echo "Outputs:"
find "$RES" -name "*.tsv" -newer "$RES/.stamp_10_scgpt_embedding" 2>/dev/null | head -30
echo
echo "Total wall: $(($(date +%s) - T0)) s"
