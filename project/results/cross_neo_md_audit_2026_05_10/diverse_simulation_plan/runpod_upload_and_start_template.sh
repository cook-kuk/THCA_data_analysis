#!/usr/bin/env bash
set -euo pipefail
POD_HOST=${POD_HOST:-135.84.176.142}
POD_PORT=${POD_PORT:-20878}
KEY=${KEY:-$HOME/.runpod/ssh/RunPod-Key-Go}
REMOTE=${REMOTE:-/runpod-volume/cross_neo_md_diverse_2026_05_10}
SCRIPT=${SCRIPT:-run_immediate_ready_batch.sh}

ssh -i "$KEY" -p "$POD_PORT" -o StrictHostKeyChecking=no root@"$POD_HOST" "mkdir -p '$REMOTE'"
rsync -az --info=progress2 -e "ssh -i $KEY -p $POD_PORT -o StrictHostKeyChecking=no" ./ root@"$POD_HOST":"$REMOTE"/
ssh -i "$KEY" -p "$POD_PORT" -o StrictHostKeyChecking=no root@"$POD_HOST" "cd '$REMOTE' && nohup bash '$SCRIPT' > '${SCRIPT%.sh}.nohup.log' 2>&1 & echo \$! > '${SCRIPT%.sh}.pid'"
