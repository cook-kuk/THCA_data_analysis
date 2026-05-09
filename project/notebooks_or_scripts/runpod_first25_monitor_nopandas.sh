#!/bin/bash
set -euo pipefail

WORK=/workspace/wsi_pathology_dm
FIRST=$WORK/wsi_resolved_first25.tsv
SECOND=$WORK/wsi_resolved_second25.tsv
LOG=/workspace/first25_monitor.log

printf "[%s] nopandas monitor start\n" "$(date -Is)" >> "$LOG"

while true; do
  DONE=$(python3 - <<'PY'
import csv
import os
import subprocess

work = "/workspace/wsi_pathology_dm"
first = f"{work}/wsi_resolved_first25.tsv"
ps = subprocess.check_output(["ps", "-ef"], text=True, errors="ignore")

ok = 0
with open(first, newline="") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        path = f"{work}/wsi/{row['file_id']}/{row['file_name']}"
        if not os.path.exists(path):
            continue
        size = os.path.getsize(path)
        expected = float(row["size_GB"]) * 1e9
        # Manifest sizes are rounded; require no active curl path to avoid active partials.
        if size >= expected * 0.85 and path not in ps:
            ok += 1
print(ok)
PY
)
  TOTAL=$(($(wc -l < "$FIRST") - 1))
  printf "[%s] first-half complete %s/%s\n" "$(date -Is)" "$DONE" "$TOTAL" >> "$LOG"

  if [ "$DONE" -ge "$TOTAL" ]; then
    printf "[%s] switching podDM to first25 manifest\n" "$(date -Is)" >> "$LOG"
    tmux kill-session -t podDM 2>/dev/null || true
    cp "$FIRST" "$WORK/wsi_resolved_manifest.tsv"
    tail -n +2 "$SECOND" | while IFS=$'\t' read -r fid fname case_id size_gb; do
      rm -rf "$WORK/wsi/$fid"
    done
    tmux new -d -s podDM "bash /workspace/runpod_pod_DM_resume.sh > /workspace/podDM.log 2>&1"
    printf "[%s] restarted podDM first25\n" "$(date -Is)" >> "$LOG"
    exit 0
  fi

  sleep 20
done
