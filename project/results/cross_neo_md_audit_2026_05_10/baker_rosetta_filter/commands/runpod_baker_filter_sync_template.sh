#!/usr/bin/env bash
set -euo pipefail
# Non-destructive template. Edit POD_HOST/POD_PORT if RunPod ports change.
POD_HOST=${POD_HOST:-135.84.176.142}
POD_PORT=${POD_PORT:-20878}
KEY=${KEY:-$HOME/.runpod/ssh/RunPod-Key-Go}
REMOTE=${REMOTE:-/workspace/cross_neo_baker_filter}
ssh -i "$KEY" -p "$POD_PORT" -o StrictHostKeyChecking=no root@"$POD_HOST" "mkdir -p $REMOTE"
scp -i "$KEY" -P "$POD_PORT" -o StrictHostKeyChecking=no /home/seungho/personal/THCA_data_analysis/project/results/cross_neo_md_audit_2026_05_10/baker_rosetta_filter/baker_rosetta_filter_manifest.tsv root@"$POD_HOST":"$REMOTE/"
# Copy PDBs separately if needed; large trajectory files are intentionally not copied here.
