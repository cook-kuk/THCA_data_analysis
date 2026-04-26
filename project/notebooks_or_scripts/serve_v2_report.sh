#!/usr/bin/env bash
set -u

ROOT="/home/seungho/personal/THCA_data_analysis/project"
PORT="${1:-8012}"
LOG_DIR="$ROOT/../project/logs"
STAMP="$(date +%Y%m%d)"
LOG_FILE="$LOG_DIR/thca_v2_report_server_${STAMP}.log"

mkdir -p "$LOG_DIR"
cd "$ROOT" || exit 1

echo "[$(date '+%F %T')] keepalive server start on 0.0.0.0:${PORT}" >> "$LOG_FILE"

while true; do
  echo "[$(date '+%F %T')] launching python http.server on ${PORT}" >> "$LOG_FILE"
  python3 -m http.server "$PORT" --bind 0.0.0.0 >> "$LOG_FILE" 2>&1
  code=$?
  echo "[$(date '+%F %T')] server exited with code ${code}; restarting in 2s" >> "$LOG_FILE"
  sleep 2
done
