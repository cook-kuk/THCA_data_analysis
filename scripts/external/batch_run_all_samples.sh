#!/usr/bin/env bash
# Batch-run PantheonOS pipeline on all 16 GSE250521 samples (3 parallel).
set -u
cd "$(dirname "$0")/../.."
PY=.venv/bin/python3
MAX_PARALLEL=3

run_one() {
    local gsm=$1 stage=$2 hist_filter=$3 sample=$4
    local viz="project/data/processed/GSE250521/${gsm}_${sample}/${gsm}_${sample}.scored.h5ad"
    local out="project/results/pantheonos_demo/${sample}"
    if [ -f "$out/run_meta.json" ]; then
        echo "[skip] $sample"
        return 0
    fi
    mkdir -p "$out"
    echo "[start] $sample stage=$stage hist=$hist_filter $(date +%H:%M:%S)"
    "$PY" scripts/external/run_pantheonos_pipeline.py \
        --visium "$viz" --histology-filter "$hist_filter" --out "$out/" \
        > "$out/run.log" 2>&1
    local rc=$?
    if [ $rc -ne 0 ]; then
        echo "[FAIL] $sample rc=$rc — see $out/run.log"
    else
        echo "[done]  $sample $(date +%H:%M:%S)"
    fi
}

# (gsm stage sc_hist sample)
SAMPLES=(
    "GSM7980860 PT NORM N-1"     "GSM7980861 PT NORM N-2"
    "GSM7980862 PT NORM N-3"     "GSM7980863 PT NORM N-4"
    "GSM7980864 PTC PTC PTC-1"   "GSM7980865 PTC PTC PTC-2"
    "GSM7980866 PTC PTC PTC-3"   "GSM7980867 PTC PTC PTC-4"
    "GSM7980868 LPTC PTC LPTC-1" "GSM7980869 LPTC PTC LPTC-2"
    "GSM7980870 LPTC PTC LPTC-3" "GSM7980871 LPTC PTC LPTC-4"
    "GSM7980872 ATC ATC ATC-1"   "GSM7980873 ATC ATC ATC-2"
    "GSM7980874 ATC ATC ATC-3"   "GSM7980875 ATC ATC ATC-4"
)

mkdir -p project/results/pantheonos_demo
for entry in "${SAMPLES[@]}"; do
    while [ "$(jobs -r | wc -l)" -ge "$MAX_PARALLEL" ]; do
        sleep 5
    done
    # shellcheck disable=SC2086
    run_one $entry &
done
wait
echo "=== batch complete $(date +%H:%M:%S) ==="
