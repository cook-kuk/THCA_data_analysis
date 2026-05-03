#!/bin/bash
# GSE250521 download + extract (Lu et al. 2024 thyroid Visium ST, 16 samples)
# Stages encoded in sample names: N-* (PT), PTC-*, LPTC-*, ATC-*
set -euo pipefail

RAW_DIR="${1:-project/data/raw/GSE250521}"
TAR_URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE250nnn/GSE250521/suppl/GSE250521_RAW.tar"
TAR_PATH="${RAW_DIR}/GSE250521_RAW.tar"

mkdir -p "$RAW_DIR"
cd "$(dirname "$0")/../../.."

if [ ! -f "$TAR_PATH" ]; then
  echo "[download] $TAR_URL → $TAR_PATH"
  curl -L --fail --retry 3 -o "$TAR_PATH" "$TAR_URL"
else
  echo "[skip] $TAR_PATH already exists ($(du -h "$TAR_PATH" | cut -f1))"
fi

EXTRACT_DIR="${RAW_DIR}/extracted"
if [ ! -d "$EXTRACT_DIR" ] || [ -z "$(ls -A "$EXTRACT_DIR" 2>/dev/null)" ]; then
  mkdir -p "$EXTRACT_DIR"
  echo "[extract] tar → $EXTRACT_DIR"
  tar -xf "$TAR_PATH" -C "$EXTRACT_DIR"
else
  echo "[skip] already extracted: $EXTRACT_DIR"
fi

echo "[done] files: $(ls "$EXTRACT_DIR" | wc -l)"
ls "$EXTRACT_DIR" | head -10
