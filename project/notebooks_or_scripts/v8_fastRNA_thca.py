#!/usr/bin/env python3
"""v8 Task 8A — FastRNA-style cohort-centering on THCA.

Lee & Han (2022) AJHG propose a simple per-cohort centering as a fast
batch correction that preserves variance but removes cohort means.
We apply this to v5.1 harmonized THCA data and recompute LODO DIAL
without ComBat. Comparison: DIAL on centered data vs. v5.1 DIAL (ComBat).

Outputs: /opt/thyroid-dash/project/results/v8_statgen/v8_fastRNA_style_dial.tsv
Columns: classifier, auc_pre_orig, auc_post_centered, dial_centered,
         dial_v5p1, delta (dial_v5p1 - dial_centered).
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.metrics import roc_auc_score

# make v5p1_common importable for classifier factories
PROJECT = Path("/opt/thyroid-dash/project")
sys.path.insert(0, str(PROJECT / "notebooks_or_scripts"))
from v5p1_common import get_classifier_factories, auc_flip  # noqa: E402

THCA_DIR = PROJECT / "data_processed" / "v5_cross_cancer" / "THCA"
RESULTS = PROJECT / "results" / "v8_statgen"
V5P1_DIAL = PROJECT / "results" / "v5" / "v5p1_dial_THCA.tsv"
LOGFILE = PROJECT / "logs" / "v8_hanlab.log"

N_TOP = 3000


def log(msg: str) -> None:
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] [T8A] {msg}"
    print(line, flush=True)
    with open(LOGFILE, "a") as fh:
        fh.write(line + "\n")


def cohort_center(X: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Subtract per-cohort mean from each sample (FastRNA-style)."""
    Xc = X.astype(np.float32).copy()
    for b in np.unique(B):
        sel = B == b
        mu = Xc[sel].mean(axis=0, keepdims=True)
        Xc[sel] = Xc[sel] - mu
    return Xc


def lodo_auc(X: np.ndarray, Y_bin: np.ndarray, B: np.ndarray, factory) -> float:
    logo = LeaveOneGroupOut()
    aucs = []
    for tr, te in logo.split(np.arange(len(Y_bin)), Y_bin, groups=B):
        if len(np.unique(Y_bin[tr])) < 2 or len(np.unique(Y_bin[te])) < 2:
            continue
        try:
            clf = factory().fit(X[tr], Y_bin[tr])
            proba = clf.predict_proba(X[te])[:, 1]
            aucs.append(roc_auc_score(Y_bin[te], proba))
        except Exception as e:
            log(f"LODO err: {type(e).__name__}: {e}")
    return float(np.mean(aucs)) if aucs else float("nan")


def dial_from_auc(auc_post: float) -> float:
    if np.isnan(auc_post):
        return float("nan")
    return (auc_flip(auc_post) - 0.5) if auc_post < 0.5 else 0.0


def main() -> None:
    log("loading X/Y/B for THCA")
    X = np.load(THCA_DIR / "X_combined.npz")["X"]
    Y = np.loadtxt(THCA_DIR / "Y.tsv", dtype=str)
    B = np.loadtxt(THCA_DIR / "B.tsv", dtype=str)
    log(f"X {X.shape}, Y {Y.shape}, B {B.shape}")

    # Top-3000 variance filter (same as v5.1)
    var_per_gene = X.var(axis=0)
    top_idx = np.argsort(var_per_gene)[-N_TOP:]
    X = X[:, top_idx]
    log(f"top-variance filter: {N_TOP} genes")

    # Binary Y
    classes_sorted = sorted(np.unique(Y).tolist())
    Y_bin = (Y == classes_sorted[0]).astype(int)
    log(f"class order (1=positive): {classes_sorted}")

    # Cohort-centering
    X_centered = cohort_center(X, B)
    log("cohort-centering done (per-cohort mean removed, variance preserved)")

    # load v5.1 DIAL
    v5p1 = pd.read_csv(V5P1_DIAL, sep="\t")
    v5p1 = v5p1.set_index("classifier")
    log(f"v5p1 DIAL loaded with classifiers: {v5p1.index.tolist()}")

    # run each classifier
    facs = get_classifier_factories()
    rows = []
    for clf_name, factory in facs.items():
        log(f"classifier={clf_name} starting")
        t0 = time.time()
        # AUC on centered data (post)
        auc_post = lodo_auc(X_centered, Y_bin, B, factory)
        # pre_orig = AUC from v5p1 (original data, before ComBat)
        auc_pre_orig = float(v5p1.loc[clf_name, "auc_pre"]) if clf_name in v5p1.index else float("nan")
        dial_v5p1 = float(v5p1.loc[clf_name, "dial"]) if clf_name in v5p1.index else float("nan")
        dial_cent = dial_from_auc(auc_post)
        delta = (dial_v5p1 - dial_cent) if (not np.isnan(dial_v5p1) and not np.isnan(dial_cent)) else float("nan")
        dt = time.time() - t0
        log(f"  auc_post_centered={auc_post:.4f}  dial_cent={dial_cent:.4f}  dial_v5p1={dial_v5p1:.4f}  delta={delta:+.4f}  [{dt:.1f}s]")
        rows.append(dict(
            classifier=clf_name,
            auc_pre_orig=auc_pre_orig,
            auc_post_centered=auc_post,
            dial_centered=dial_cent,
            dial_v5p1=dial_v5p1,
            delta=delta,
        ))

    out = pd.DataFrame(rows)
    out_path = RESULTS / "v8_fastRNA_style_dial.tsv"
    out.to_csv(out_path, sep="\t", index=False)
    log(f"wrote {out_path}")
    log(f"mean dial_centered={out['dial_centered'].mean():.4f}  mean dial_v5p1={out['dial_v5p1'].mean():.4f}")
    log("T8A DONE")


if __name__ == "__main__":
    main()
