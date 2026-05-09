#!/usr/bin/env bash
# K2 PRJEB11591 full-transcriptome kallisto quant
# Runs paired-end where both R1+R2 present, single-end fallback for R1-only samples.
# Per CLAUDE.md: --tmpdir on /data, --compress not supported by kallisto (kept for arcasHLA-style runs only)
#
# n_threads_per_quant set to 4; with 8 cores -> 2 quants in parallel.
# Wall budget per sample (paired-end RNA-seq, ~50M reads): 5-12 min.
set -euo pipefail

KALLISTO="/opt/thyroid-dash/project/.venv/lib/python3.12/site-packages/kb_python/bins/linux/kallisto/kallisto"
INDEX="/data/thca/reference_kallisto/gencode.v44.kallisto.idx.NEW"
FASTQ_ROOT="/data/thca/PRJEB11591_fastq"
OUT_ROOT="/data/thca/v17_korean_k2_full_quant"
LOG="/home/seungho/personal/THCA_data_analysis/project/results/p2_braf_nature_sprint_2026_05_09/h_k2_full_quant/kallisto_quant.log"
TMPDIR="/data/thca/_tmp"
export TMPDIR

mkdir -p "$OUT_ROOT" "$TMPDIR"
: > "$LOG"

echo "[$(date -Iseconds)] kallisto $($KALLISTO version | head -1)  index=$INDEX" | tee -a "$LOG"
echo "[$(date -Iseconds)] index size: $(stat -c %s $INDEX) bytes" | tee -a "$LOG"

quant_one() {
  local sample="$1"
  local fq_dir="$FASTQ_ROOT/$sample"
  local out="$OUT_ROOT/$sample"
  local r1="$fq_dir/${sample}_1.fastq.gz"
  local r2="$fq_dir/${sample}_2.fastq.gz"
  mkdir -p "$out"
  local t0=$(date +%s)
  if [[ -s "$r1" && -s "$r2" ]]; then
    "$KALLISTO" quant -t 4 -i "$INDEX" -o "$out" "$r1" "$r2" \
      >>"$out/kallisto.stderr.log" 2>&1
    local mode="paired"
  elif [[ -s "$r1" ]]; then
    "$KALLISTO" quant -t 4 -i "$INDEX" -o "$out" --single -l 200 -s 30 "$r1" \
      >>"$out/kallisto.stderr.log" 2>&1
    local mode="single"
  else
    echo "[$(date -Iseconds)] $sample SKIP: no FASTQ" | tee -a "$LOG"
    return
  fi
  local t1=$(date +%s)
  local dur=$((t1 - t0))
  local nproc=$(grep -oP '"n_processed":\s*\K[0-9]+' "$out/run_info.json" 2>/dev/null || echo NA)
  local npseudo=$(grep -oP '"p_pseudoaligned":\s*\K[0-9.]+' "$out/run_info.json" 2>/dev/null || echo NA)
  echo "[$(date -Iseconds)] $sample $mode  ${dur}s  reads=$nproc  p_pseudoaligned=$npseudo" | tee -a "$LOG"
}
export -f quant_one
export KALLISTO INDEX FASTQ_ROOT OUT_ROOT LOG

# Enumerate samples that actually have FASTQ
SAMPLES=$(find "$FASTQ_ROOT" -maxdepth 2 -name "*_1.fastq.gz" -printf '%f\n' | sed 's/_1.fastq.gz//' | sort -u)
N=$(echo "$SAMPLES" | wc -l)
echo "[$(date -Iseconds)] N samples with FASTQ: $N" | tee -a "$LOG"

T_START=$(date +%s)
# Run 2 in parallel (kallisto -t 4 each = 8 cores total)
echo "$SAMPLES" | xargs -I{} -P2 bash -c 'quant_one "$@"' _ {}
T_END=$(date +%s)
TOTAL=$((T_END - T_START))
echo "[$(date -Iseconds)] DONE total wall ${TOTAL}s for N=$N samples" | tee -a "$LOG"
