#!/usr/bin/env python3
"""v5.1 Phase 3 — DIAL computation under the F1 leakage-safe ComBat fix.

Identical to v5p1_dial.py except:
  - per-fold ComBat fit on train, applied to held-out test cohort
    (compute_dial(..., leak_safe=True) is the new default).
  - results saved to v5p1_dial_*_F1.tsv to keep the v5.1a artifacts intact.
"""
from __future__ import annotations

import sys
import time
import traceback
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd

sys.path.insert(0, "/opt/thyroid-dash/project/notebooks_or_scripts")
from v5p1_common import (
    DATA_PROC_V5, RESULTS_V5, LOGS, CANCERS, log_line, run_one_classifier,
    get_classifier_factories,
)

LOG = LOGS / "v5p1_dial_F1.log"


def run_for_cancer(cancer: str) -> pd.DataFrame:
    out_dir = DATA_PROC_V5 / cancer
    X_path = out_dir / "X_combined.npz"
    Y_path = out_dir / "Y.tsv"
    B_path = out_dir / "B.tsv"
    if not (X_path.exists() and Y_path.exists() and B_path.exists()):
        log_line(LOG, f"[{cancer}] skip - missing harmonized data")
        return pd.DataFrame()

    facs = get_classifier_factories()
    tasks = [dict(cancer=cancer, clf_name=k,
                  X_path=str(X_path), Y_path=str(Y_path), B_path=str(B_path))
             for k in facs.keys()]

    results = []
    with ProcessPoolExecutor(max_workers=min(5, len(tasks))) as ex:
        futs = {ex.submit(run_one_classifier, t): t for t in tasks}
        for fut in as_completed(futs):
            try:
                r = fut.result()
                results.append(r)
                log_line(LOG, f"[{cancer}/{r['classifier']}] dial={r['dial']:.3f} "
                              f"dial_v2={r['dial_v2']:.3f} pflip={r['partial_flip_score']:.3f} "
                              f"auc_pre={r['auc_pre']:.3f} auc_post={r['auc_post']:.3f} "
                              f"interp={r['interpretation']} {r['seconds']}s")
            except Exception as e:
                t = futs[fut]
                log_line(LOG, f"[{cancer}/{t['clf_name']}] ERROR {e}")
    return pd.DataFrame(results)


def main():
    log_line(LOG, "=== v5p1 Phase 3 DIAL F1 (leak-safe ComBat) START ===")
    all_rows = []
    failed = []
    for cancer in CANCERS:
        try:
            df = run_for_cancer(cancer)
        except Exception as e:
            tb = traceback.format_exc()
            log_line(LOG, f"[{cancer}] FATAL {type(e).__name__}: {e}")
            fail_path = RESULTS_V5 / f"v5p1_dial_{cancer}_F1_FAILED.tsv"
            fail_path.write_text(
                "cancer\terror_type\terror_message\ttraceback\n"
                f"{cancer}\t{type(e).__name__}\t{e}\t{tb!r}\n"
            )
            print(f"[{cancer}] FATAL:\n{tb}", file=sys.stderr, flush=True)
            failed.append(cancer)
            continue
        if df.empty:
            continue
        df.to_csv(RESULTS_V5 / f"v5p1_dial_{cancer}_F1.tsv", sep="\t", index=False)
        all_rows.append(df)
        log_line(LOG, f"[{cancer}] wrote v5p1_dial_{cancer}_F1.tsv with {len(df)} rows")

    if all_rows:
        full = pd.concat(all_rows, ignore_index=True)
        if "semi_synthetic" not in full.columns:
            full["semi_synthetic"] = False
        full["semi_synthetic"] = False
        full.to_csv(RESULTS_V5 / "v5p1_dial_all_cancers_F1.tsv", sep="\t", index=False)
        log_line(LOG, f"wrote v5p1_dial_all_cancers_F1.tsv with {len(full)} rows")
    else:
        log_line(LOG, "NO DIAL results to write - PARTIAL_FAILURE")
    if failed:
        log_line(LOG, f"FAILED cancers: {','.join(failed)} - see v5p1_dial_*_F1_FAILED.tsv")
    log_line(LOG, "=== v5p1 Phase 3 F1 DONE ===")


if __name__ == "__main__":
    main()
