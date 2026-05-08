#!/usr/bin/env bash
# Launch all 4 NComm-push public-data analyses in parallel.
# Logs: ./logs/<worker>.log; status: ./<worker>/_status.json
# Resource layout: 8 cores, 14GB free RAM -> 2 cores per worker, ~3.5GB each.
set -u
cd "$(dirname "$0")"
ROOT="$(pwd)"

WORKERS=(cbioportal_sweep depmap_thyroid cptac_thca geo_meta)
PIDS=()

mkdir -p logs
echo "Launching ${#WORKERS[@]} workers in parallel at $(date)" | tee logs/launcher.log

for w in "${WORKERS[@]}"; do
    log="logs/${w}.log"
    echo "  -> $w (log=$log)" | tee -a logs/launcher.log
    # Each worker gets a clean log; nohup so it survives session disconnect
    : > "$log"
    nohup python3 "${w}/run.py" >> "$log" 2>&1 &
    PIDS+=("$!")
    echo "     started PID $!" | tee -a logs/launcher.log
done

echo "All workers launched. PIDs: ${PIDS[*]}" | tee -a logs/launcher.log
echo "$(date +%s)" > .launch_timestamp
printf '%s\n' "${PIDS[@]}" > .pids

echo
echo "Monitor with:"
echo "  tail -f logs/*.log"
echo "  watch -n 5 'for d in cbioportal_sweep depmap_thyroid cptac_thca geo_meta; do echo \"== \$d ==\"; cat \$d/_status.json 2>/dev/null || echo no-status; done'"
