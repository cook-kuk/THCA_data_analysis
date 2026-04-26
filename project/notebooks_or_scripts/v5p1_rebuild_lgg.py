#!/usr/bin/env python3
"""One-off: rebuild LGG tcga_counts.tsv.gz after the partial-write incident."""
from __future__ import annotations

import sys
sys.path.insert(0, "/opt/thyroid-dash/project/notebooks_or_scripts")

from v5p1_common import log_line, LOGS, DATA_RAW_V5
from v5p1_download_real import download_tcga_rnaseq, extract_driver_labels

LOG = LOGS / "v5p1_download.log"

log_line(LOG, "=== LGG rebuild START ===")

# Remove the broken gz so rebuild starts clean
broken = DATA_RAW_V5 / "LGG" / "tcga_counts.tsv.gz"
if broken.exists():
    log_line(LOG, f"[LGG] removing broken {broken.name} ({broken.stat().st_size} bytes)")
    broken.unlink()

r = download_tcga_rnaseq("LGG")
log_line(LOG, f"[LGG rebuild rnaseq] {r}")

# Re-extract driver labels since they depend on MAF (unchanged) — keep them in sync
labels = extract_driver_labels("LGG")
if labels is not None:
    labels.to_csv(DATA_RAW_V5 / "LGG" / "tcga_driver_labels.tsv", sep="\t", index=False)
    log_line(LOG, f"[LGG rebuild labels] n={len(labels)}")

log_line(LOG, "=== LGG rebuild DONE ===")
