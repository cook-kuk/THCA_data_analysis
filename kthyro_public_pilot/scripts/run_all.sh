#!/usr/bin/env bash
set -u
set -o pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
mkdir -p results/logs results/tables results/figures results/reports

run_stage() {
  local script="$1"
  shift || true
  local log="results/logs/${script%.py}.stdout_stderr.log"
  echo "[$(date '+%F %T')] running $script $*" | tee -a "$log"
  if python "scripts/$script" "$@" >>"$log" 2>&1; then
    echo "[$(date '+%F %T')] completed $script" | tee -a "$log"
  else
    echo "[$(date '+%F %T')] ERROR in $script; see $log" | tee -a "$log"
    exit 1
  fi
}

run_stage 00_audit_environment.py
run_stage 01_audit_local_data.py
run_stage 02_download_or_link_tcga_thca.py
run_stage 03_tcga_bulk_vulnerability.py
run_stage 04_download_or_link_spatial_thyroid.py
run_stage 05_spatial_vulnerability_mapping.py
run_stage 06_scrna_reference_audit.py
run_stage 07_depm_prism_drug_pilot.py
run_stage 08_integrate_evidence.py
run_stage 09_make_proposal_figures.py
run_stage 10_write_report.py

echo "[$(date '+%F %T')] K-Thyro public pilot complete."
