#!/usr/bin/env python3
"""One-off: re-harmonize LGG only, update the LGG row in v5p1_harmonization.tsv.
Leaves THCA/SKCM/LUAD/COAD rows untouched."""
from __future__ import annotations

import sys
import traceback
sys.path.insert(0, "/opt/thyroid-dash/project/notebooks_or_scripts")

import pandas as pd
from v5p1_harmonize import harmonize_cancer
from v5p1_common import log_line, LOGS, RESULTS_V5

LOG = LOGS / "v5p1_harmonize.log"
CANCER = "LGG"

log_line(LOG, f"=== {CANCER} harmonize re-run START ===")
try:
    r = harmonize_cancer(CANCER)
except Exception as e:
    log_line(LOG, f"[{CANCER}] ERROR: {traceback.format_exc()}")
    r = {"status": "error", "reason": str(e)}

log_line(LOG, f"[{CANCER} rerun result] {r}")

new_row = {
    "cancer": CANCER,
    "status": r.get("status"),
    "n_cohorts": r.get("n_cohorts", 0),
    "cohort_ids": ",".join(r.get("cohort_ids", [])),
    "n_shared_genes": r.get("n_shared_genes", 0),
    "n_samples": r.get("n_samples", 0),
    "n_class_A": r.get("n_a", 0),
    "n_class_B": r.get("n_b", 0),
    "note": r.get("reason", ""),
    "semi_synthetic": False,
}

tsv_path = RESULTS_V5 / "v5p1_harmonization.tsv"
df = pd.read_csv(tsv_path, sep="\t")
mask = df["cancer"] == CANCER
if not mask.any():
    df = pd.concat([df, pd.DataFrame([new_row])], ignore_index=True)
else:
    for col, val in new_row.items():
        df.loc[mask, col] = val

# atomic write
tmp_path = tsv_path.with_suffix(tsv_path.suffix + ".part")
df.to_csv(tmp_path, sep="\t", index=False)
tmp_path.replace(tsv_path)

log_line(LOG, f"=== {CANCER} harmonize re-run DONE ===")
print(f"UPDATED {CANCER} row:")
print(df[df["cancer"] == CANCER].to_string(index=False))
