#!/usr/bin/env python3
"""v5 DIAL Phase 3 — Compute DIAL per (cancer × classifier).

5 classifiers per cancer: LogReg_l2, LogReg_elasticnet, RandomForest,
GradientBoosting, XGBoost. ProcessPoolExecutor(max_workers=5) per cancer.
Sequential across cancers with gc between.
"""
from __future__ import annotations

import gc
import sys
import time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from v5_dial_common import (
    COHORTS, DATA_PROC_V5, RESULTS_V5, LOGS, log_line,
    get_classifier_factories, run_one_classifier, peak_rss_gb,
)

LOGFILE = LOGS / "v5_dial_run.log"


def process_one_cancer(cancer: str, semi: bool) -> list[dict]:
    outdir = DATA_PROC_V5 / cancer
    X_path = outdir / "X_combined.npz"
    Y_path = outdir / "Y.tsv"
    B_path = outdir / "B.tsv"
    if not X_path.exists():
        log_line(LOGFILE, f"PHASE3 MISSING X for {cancer}, skip")
        return []

    classifiers = list(get_classifier_factories().keys())
    tasks = [
        dict(cancer=cancer, clf_name=clf,
             X_path=str(X_path), Y_path=str(Y_path), B_path=str(B_path),
             semi_synthetic=semi)
        for clf in classifiers
    ]
    rows = []
    t0 = time.time()
    # ProcessPoolExecutor per cancer
    with ProcessPoolExecutor(max_workers=min(5, len(classifiers))) as ex:
        futs = {ex.submit(run_one_classifier, t): t for t in tasks}
        for fut in as_completed(futs):
            t = futs[fut]
            try:
                res = fut.result()
                rows.append(res)
                log_line(LOGFILE, f"PHASE3 {cancer}/{t['clf_name']} DIAL={res['dial']:.3f} AUC_pre={res['auc_pre']:.3f} AUC_post={res['auc_post']:.3f} ident={res['batch_identifiability_post']:.3f} ({res['interpretation']}) {res['seconds']:.1f}s")
            except Exception as e:
                log_line(LOGFILE, f"PHASE3 FAIL {cancer}/{t['clf_name']}: {type(e).__name__}: {e}")
                rows.append(dict(
                    cancer=cancer, classifier=t["clf_name"],
                    auc_pre=float("nan"), auc_flip_pre=float("nan"),
                    auc_post=float("nan"), auc_flip_post=float("nan"),
                    dial=float("nan"), batch_identifiability_post=float("nan"),
                    interpretation="error", semi_synthetic=semi,
                    n_samples=0, n_genes=0, seconds=0,
                ))
    dt = time.time() - t0
    log_line(LOGFILE, f"PHASE3 {cancer} total {dt:.1f}s peak_RSS={peak_rss_gb():.2f}GB")
    gc.collect()
    return rows


def main():
    log_line(LOGFILE, "PHASE3 start — compute DIAL")
    all_rows = []
    # read harmonization status to know which are semi-synthetic
    harm = pd.read_csv(RESULTS_V5 / "v5_dial_harmonization.tsv", sep="\t")
    semi_map = dict(zip(harm["cancer"], harm["semi_synthetic"]))

    for cancer in ["THCA", "SKCM", "LGG", "LUAD", "COAD"]:
        rows = process_one_cancer(cancer, bool(semi_map.get(cancer, True)))
        if rows:
            pd.DataFrame(rows).to_csv(RESULTS_V5 / f"v5_dial_{cancer}.tsv", sep="\t", index=False)
            all_rows.extend(rows)

    df = pd.DataFrame(all_rows)
    col_order = [
        "cancer", "classifier", "auc_pre", "auc_flip_pre", "auc_post",
        "auc_flip_post", "dial", "batch_identifiability_post",
        "interpretation", "n_samples", "n_genes", "semi_synthetic", "seconds",
    ]
    df = df[[c for c in col_order if c in df.columns]]
    df.to_csv(RESULTS_V5 / "v5_dial_all_cancers.tsv", sep="\t", index=False)
    log_line(LOGFILE, f"PHASE3 wrote v5_dial_all_cancers.tsv n_rows={len(df)}")
    log_line(LOGFILE, "PHASE3 done")


if __name__ == "__main__":
    main()
