#!/usr/bin/env python3
"""v5.2 Task A1 runner — re-run DIAL for all 5 cancers with PROPER LODO ComBat.

Reads the same processed X/Y/B matrices used by v5.1 (so any change in DIAL
is attributable to the ComBat refactor, not data).

Output:
  results/v5p2_fix/v5p2_dial_proper_lodo.tsv     — 25 rows (new values)
  results/v5p2_fix/v5p2_lodo_comparison.tsv      — side-by-side with v5.1
"""
from __future__ import annotations

import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "/opt/thyroid-dash/project/notebooks_or_scripts")
from v5p1_common import (
    CANCERS, DATA_PROC_V5, LOGS, PROJECT, RESULTS_V5,
    get_classifier_factories, log_line,
)

RES = PROJECT / "results" / "v5p2_fix"
RES.mkdir(parents=True, exist_ok=True)
LOG = LOGS / "v5p2_critical_fix.log"


def run_one(task: dict) -> dict:
    """ProcessPool-safe single-classifier DIAL run with proper LODO ComBat."""
    import numpy as np
    from sklearn.model_selection import LeaveOneGroupOut
    from v5p1_common import get_classifier_factories
    from v5p2_combat_lodo import compute_dial_lodo_proper

    cancer = task["cancer"]
    clf_name = task["clf_name"]
    X = np.load(task["X_path"])["X"]
    Y = np.loadtxt(task["Y_path"], dtype=str)
    B = np.loadtxt(task["B_path"], dtype=str)

    facs = get_classifier_factories()
    factory = facs[clf_name]

    Y_bin = (Y == sorted(np.unique(Y).tolist())[0]).astype(int)
    logo = LeaveOneGroupOut()
    splits = list(logo.split(np.arange(len(Y_bin)), Y_bin, groups=B))

    # Same top-variance filter as v5.1 (3000 genes) for apples-to-apples comparison
    n_top = 3000
    n_orig = X.shape[1]
    if X.shape[1] > n_top:
        var_per_gene = X.var(axis=0)
        top_idx = np.argsort(var_per_gene)[-n_top:]
        X = X[:, top_idx]

    t0 = time.time()
    res = compute_dial_lodo_proper(X.astype(np.float32), Y, B, factory, splits)
    dt = time.time() - t0
    res.update(dict(cancer=cancer, classifier=clf_name, n_samples=len(Y),
                    n_genes=X.shape[1], semi_synthetic=False, seconds=round(dt, 2),
                    n_orig_genes=n_orig))
    return res


def main():
    log_line(LOG, "=== v5.2 proper-LODO DIAL START ===")
    facs = get_classifier_factories()
    tasks = []
    for cancer in CANCERS:
        proc_dir = DATA_PROC_V5 / cancer
        if not (proc_dir / "X_combined.npz").exists():
            log_line(LOG, f"[{cancer}] skip — missing processed data")
            continue
        for clf_name in facs.keys():
            tasks.append(dict(
                cancer=cancer, clf_name=clf_name,
                X_path=str(proc_dir / "X_combined.npz"),
                Y_path=str(proc_dir / "Y.tsv"),
                B_path=str(proc_dir / "B.tsv"),
            ))

    log_line(LOG, f"Dispatching {len(tasks)} (cancer × classifier) tasks")

    results = []
    # Run with modest parallelism — ComBat fit is OLS so each job is ~5-30s
    max_workers = min(6, max(1, len(tasks)))
    with ProcessPoolExecutor(max_workers=max_workers) as ex:
        futs = {ex.submit(run_one, t): t for t in tasks}
        for fut in as_completed(futs):
            t = futs[fut]
            try:
                r = fut.result()
                # Drop per_fold_auc_post list (emit elsewhere) — flatten for TSV
                per_fold = r.pop("per_fold_auc_post", [])
                r["per_fold_auc_post"] = ";".join(f"{h}:{a:.4f}" for h, a in per_fold)
                results.append(r)
                log_line(LOG, f"[{t['cancer']}/{t['clf_name']}] "
                              f"dial={r['dial']:.4f} auc_pre={r['auc_pre']:.3f} "
                              f"auc_post={r['auc_post']:.3f} interp={r['interpretation']} "
                              f"{r['seconds']}s")
            except Exception as e:
                tb = traceback.format_exc()
                log_line(LOG, f"[{t['cancer']}/{t['clf_name']}] ERROR {type(e).__name__}: {e}\n{tb}")

    # Cast to DataFrame
    df = pd.DataFrame(results)
    col_order = ["cancer", "classifier", "auc_pre", "auc_post", "auc_flip_pre",
                 "auc_flip_post", "dial", "interpretation",
                 "batch_identifiability_post", "n_samples", "n_genes", "n_orig_genes",
                 "seconds", "semi_synthetic", "per_fold_auc_post"]
    df = df[[c for c in col_order if c in df.columns]]
    # add is_label_flip for parity with v5.1
    df["is_label_flip"] = df["auc_post"] < 0.5

    out = RES / "v5p2_dial_proper_lodo.tsv"
    df.to_csv(out, sep="\t", index=False)
    log_line(LOG, f"wrote {out} ({len(df)} rows)")
    print(f"wrote {out}")

    # Side-by-side comparison
    v51_path = RESULTS_V5 / "v5p1_dial_all_cancers.tsv"
    if v51_path.exists():
        v51 = pd.read_csv(v51_path, sep="\t")
        merged = v51.merge(
            df[["cancer", "classifier", "auc_post", "dial", "interpretation",
                "is_label_flip"]]
            .rename(columns={
                "auc_post": "auc_post_v52",
                "dial": "dial_v52",
                "interpretation": "interp_v52",
                "is_label_flip": "is_flip_v52",
            }),
            on=["cancer", "classifier"], how="left",
        )
        merged = merged.rename(columns={
            "auc_post": "auc_post_v51",
            "dial": "dial_v51",
            "interpretation": "interp_v51",
            "is_label_flip": "is_flip_v51",
        })
        merged["delta_dial"] = merged["dial_v52"] - merged["dial_v51"]
        merged["delta_auc_post"] = merged["auc_post_v52"] - merged["auc_post_v51"]
        keep = ["cancer", "classifier",
                "auc_post_v51", "auc_post_v52", "delta_auc_post",
                "dial_v51", "dial_v52", "delta_dial",
                "interp_v51", "interp_v52",
                "is_flip_v51", "is_flip_v52"]
        keep = [c for c in keep if c in merged.columns]
        cmp = merged[keep]
        cmp_out = RES / "v5p2_lodo_comparison.tsv"
        cmp.to_csv(cmp_out, sep="\t", index=False)
        log_line(LOG, f"wrote {cmp_out}")
        print(f"wrote {cmp_out}")
        print("\n=== v5.1 vs v5.2 (proper LODO ComBat) ===")
        print(cmp.to_string(index=False))

    log_line(LOG, "=== v5.2 proper-LODO DIAL DONE ===")


if __name__ == "__main__":
    main()
