#!/usr/bin/env python3
"""v9 — playground grid precompute.

Generates DIAL values for (cancer × classifier × feature_count × cv_scheme)
so the interactive slider can do a cheap JS-side lookup.

Real computation with 5×5×5×2 = 250 runs is ~15 min. Here we run a
smaller real grid (feature counts {300, 1000, 3000}) and interpolate the
remaining counts {100, 2000} in JS. StratifiedKFold is run as a cheap
in-training shuffle rather than a full re-fit (approximation marked in
the JSON).
"""
from __future__ import annotations

import json
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

sys.path.insert(0, "/opt/thyroid-dash/project/notebooks_or_scripts")
from v5p1_common import (
    CANCERS, DATA_PROC_V5, _combat_preserve, get_classifier_factories, LOGS, log_line,
)

OUT = Path("/opt/thyroid-dash/project/reports/html/assets/v9/data/v9_grid.json")
LOG = LOGS / "v9_grid.log"

FEATURE_COUNTS = [300, 1000, 3000]
CV_SCHEMES = ["LODO", "StratifiedKFold"]


def auc_flip(a: float) -> float:
    return max(a, 1.0 - a)


def dial_for(cancer: str, clf_name: str, n_top: int, cv_scheme: str) -> dict:
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import LeaveOneGroupOut, StratifiedKFold
    from sklearn.metrics import roc_auc_score

    d = DATA_PROC_V5 / cancer
    X = np.load(d / "X_combined.npz")["X"]
    Y = np.loadtxt(d / "Y.tsv", dtype=str)
    B = np.loadtxt(d / "B.tsv", dtype=str)

    var = X.var(axis=0)
    idx = np.argsort(var)[-min(n_top, X.shape[1]):]
    X = X[:, idx].astype(np.float32)

    classes = sorted(np.unique(Y).tolist())
    Y_bin = (Y == classes[0]).astype(int)

    facs = get_classifier_factories()
    factory = facs.get(clf_name)
    if factory is None:
        return dict(error=f"missing classifier {clf_name}")

    if cv_scheme == "LODO":
        logo = LeaveOneGroupOut()
        splits = list(logo.split(np.arange(len(Y_bin)), Y_bin, groups=B))
    else:
        skf = StratifiedKFold(n_splits=3, shuffle=True, random_state=42)
        splits = list(skf.split(np.arange(len(Y_bin)), Y_bin))

    # pre
    pre = []
    for tr, te in splits:
        try:
            clf = factory().fit(X[tr], Y_bin[tr])
            proba = clf.predict_proba(X[te])[:, 1]
            pre.append(roc_auc_score(Y_bin[te], proba))
        except Exception as e:
            log_line(LOG, f"[{cancer}/{clf_name}/{n_top}/{cv_scheme}] pre err: {e}")
    auc_pre = float(np.mean(pre)) if pre else float("nan")

    # post — v5.2 protocol: per-fold ComBat fit on TRAIN only, held-out cohort
    # centered locally (no Y leak). Avoids the v5.1 pooled-X leakage.
    post = []
    for tr, te in splits:
        try:
            Xtrain = _combat_preserve(X[tr], Y[tr], B[tr]).astype(np.float32)
            Xtest = X[te].astype(np.float32)
            Xtest = Xtest - Xtest.mean(axis=0) + Xtrain.mean(axis=0)
            clf = factory().fit(Xtrain, Y_bin[tr])
            proba = clf.predict_proba(Xtest)[:, 1]
            post.append(roc_auc_score(Y_bin[te], proba))
        except Exception as e:
            log_line(LOG, f"[{cancer}/{clf_name}/{n_top}/{cv_scheme}] post err: {e}")
    auc_post = float(np.mean(post)) if post else float("nan")

    dial_val = (auc_flip(auc_post) - 0.5) if (not np.isnan(auc_post) and auc_post < 0.5) else 0.0

    if np.isnan(dial_val) or np.isnan(auc_post):
        interp = "ambiguous"
    elif dial_val > 0.3:
        interp = "batch_entangled"
    elif dial_val > 0.1:
        interp = "partial_batch"
    elif auc_post > 0.7 and dial_val < 0.05:
        interp = "true_biology"
    elif auc_post < 0.6:
        interp = "no_signal"
    else:
        interp = "ambiguous"

    return dict(
        cancer=cancer,
        classifier=clf_name,
        feature_count=n_top,
        cv_scheme=cv_scheme,
        auc_pre=auc_pre if not np.isnan(auc_pre) else None,
        auc_post=auc_post if not np.isnan(auc_post) else None,
        dial=float(dial_val),
        interpretation=interp,
    )


def main():
    LOG.parent.mkdir(parents=True, exist_ok=True)
    log_line(LOG, "=== v9 grid START ===")
    t0 = time.time()
    classifiers = list(get_classifier_factories().keys())
    rows = []
    for cancer in CANCERS:
        for clf_name in classifiers:
            for fc in FEATURE_COUNTS:
                for cv in CV_SCHEMES:
                    try:
                        r = dial_for(cancer, clf_name, fc, cv)
                    except Exception as e:
                        log_line(LOG, f"[{cancer}/{clf_name}/{fc}/{cv}] FAIL {e}")
                        r = dict(cancer=cancer, classifier=clf_name,
                                 feature_count=fc, cv_scheme=cv,
                                 dial=None, auc_pre=None, auc_post=None,
                                 interpretation="ambiguous", error=str(e))
                    rows.append(r)
                    log_line(LOG, f"[{cancer}/{clf_name}/{fc}/{cv}] dial={r.get('dial')}")

    # Index by (cancer, classifier, feature_count, cv_scheme) for fast lookup.
    index = {}
    for r in rows:
        key = f"{r['cancer']}|{r['classifier']}|{r['feature_count']}|{r['cv_scheme']}"
        index[key] = r
    payload = dict(
        meta=dict(
            generated_at=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            feature_counts=FEATURE_COUNTS,
            cv_schemes=CV_SCHEMES,
            classifiers=classifiers,
            cancers=CANCERS,
            runtime_seconds=round(time.time() - t0, 1),
            note="Other feature_counts (100, 2000) are interpolated on the client.",
        ),
        rows=rows,
        index=index,
    )
    OUT.write_text(json.dumps(payload, separators=(",", ":")))
    size_kb = OUT.stat().st_size / 1024
    log_line(LOG, f"Wrote {OUT} ({size_kb:.1f} KB)")
    print(f"OK {OUT} {size_kb:.1f}KB rows={len(rows)} runtime={time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
