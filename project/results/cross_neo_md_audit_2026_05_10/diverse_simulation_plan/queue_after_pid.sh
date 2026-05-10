#!/usr/bin/env bash
set -euo pipefail
WAIT_PID=${1:?usage: queue_after_pid.sh WAIT_PID [SCRIPT]}
SCRIPT=${2:-run_immediate_ready_batch.sh}
POLL_SECONDS=${POLL_SECONDS:-300}
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "[queue] waiting for PID ${WAIT_PID} before ${SCRIPT}"
while kill -0 "$WAIT_PID" 2>/dev/null; do
  date -Is
  sleep "$POLL_SECONDS"
done
echo "[queue] PID ${WAIT_PID} finished; starting ${SCRIPT}"
bash "$SCRIPT"
