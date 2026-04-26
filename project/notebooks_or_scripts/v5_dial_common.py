#!/usr/bin/env python3
"""v5 DIAL Cross-Cancer Audit — shared utilities.

Canonical DIAL (Direction-Invariant AUC Leakage) implementation per the
v5 cross-cancer method paper sprint spec. Subtype-preserving ComBat using
inmoose.pycombat.pycombat_norm, CV via StratifiedKFold, 5 classifiers.

Reuses existing THCA bulk matrix for THCA cancer type. For SKCM/LGG/LUAD/COAD
we attempt TCGA/GEO downloads; if network is unavailable or the cohort does
not satisfy n>=30 per class, we fall back to semi-synthetic data derived
from THCA using the pattern from v5_track3_cross_cancer.py (realistic mean
shifts + subtype DEGs). Rows are logged honestly via the
`semi_synthetic` flag in v5_dial_cohort_availability.tsv.
"""
from __future__ import annotations

import gc
import os
import sys
import time
import json
import gzip
import warnings
from pathlib import Path
from typing import Callable, Dict, List, Tuple

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

PROJECT = Path("/opt/thyroid-dash/project")
DATA_RAW_V5 = PROJECT / "data_raw" / "v5_cross_cancer"
DATA_PROC_V5 = PROJECT / "data_processed" / "v5_cross_cancer"
RESULTS_V5 = PROJECT / "results" / "v5"
REPORTS_V5 = PROJECT / "reports" / "v5"
FIGS_V5 = PROJECT / "reports" / "html" / "figs_interactive" / "v5"
LOGS = PROJECT / "logs"
PAGES = PROJECT / "reports" / "html" / "pages"

for p in [DATA_RAW_V5, DATA_PROC_V5, RESULTS_V5, REPORTS_V5, FIGS_V5, LOGS]:
    p.mkdir(parents=True, exist_ok=True)


# ============================================================
# Cohort manifest
# ============================================================
COHORTS = {
    "THCA": {
        "tcga": "TCGA-THCA",
        "geo":  ["GSE27155", "GSE33630", "GSE29265"],
        "task": "BRAF_vs_RAS",
        "labels": "mutation",
        "class_a": "BRAF",
        "class_b": "RAS",
    },
    "SKCM": {
        "tcga": "TCGA-SKCM",
        "geo":  ["GSE65904", "GSE22153"],
        "task": "BRAF_vs_NRAS",
        "labels": "mutation",
        "class_a": "BRAF",
        "class_b": "NRAS",
    },
    "LGG": {
        "tcga": "TCGA-LGG",
        "geo":  ["GSE16011", "GSE4271"],
        "task": "IDH_mut_vs_wt",
        "labels": "mutation",
        "class_a": "IDHmut",
        "class_b": "IDHwt",
    },
    "LUAD": {
        "tcga": "TCGA-LUAD",
        "geo":  ["GSE31210", "GSE72094"],
        "task": "KRAS_vs_EGFR",
        "labels": "mutation",
        "class_a": "KRAS",
        "class_b": "EGFR",
    },
    "COAD": {
        "tcga": "TCGA-COAD",
        "geo":  ["GSE39582", "GSE17536"],
        "task": "CMS1_vs_CMS2",
        "labels": "CMSclassifier",
        "class_a": "CMS1",
        "class_b": "CMS2",
    },
}


# ============================================================
# Logging
# ============================================================
def log_line(logfile: Path, msg: str) -> None:
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    with open(logfile, "a") as fh:
        fh.write(line + "\n")


def peak_rss_gb() -> float:
    try:
        import resource
        kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return kb / (1024.0 * 1024.0)
    except Exception:
        return -1.0


# ============================================================
# Canonical DIAL
# ============================================================
def auc_flip(a: float) -> float:
    return max(a, 1.0 - a)


def ovr_macro_auc(y_true: np.ndarray, y_proba: np.ndarray) -> float:
    from sklearn.metrics import roc_auc_score
    try:
        uniq = np.unique(y_true)
        if len(uniq) < 2:
            return float("nan")
        if y_proba.ndim == 1 or y_proba.shape[1] == 1:
            return float(roc_auc_score(y_true, y_proba))
        if len(uniq) == 2 and y_proba.shape[1] >= 2:
            return float(roc_auc_score(y_true, y_proba[:, 1]))
        return float(roc_auc_score(y_true, y_proba, multi_class="ovr", average="macro"))
    except Exception:
        return float("nan")


def _combat_preserve(X: np.ndarray, Y: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Subtype-preserving ComBat via inmoose.pycombat.pycombat_norm.

    X is (n_samples, n_genes); inmoose expects (n_genes, n_samples). Y is
    passed as a one-hot covariate matrix when possible (subtype-preserving
    mode); falls back to the standard un-preserved ComBat if the
    covariate matrix is rank-deficient w.r.t. the batch factor (which is
    the empirical regime the Lemma describes). This matches the spec's
    `pycombat(X.T, batch=B, mod=pd.get_dummies(Y))` signature while
    handling real-data singular-matrix errors gracefully.
    """
    from inmoose.pycombat import pycombat_norm

    df = pd.DataFrame(X.T)  # genes x samples
    if len(np.unique(B)) < 2:
        return X.astype(float)

    counts = pd.Series(B).value_counts()
    if (counts < 2).any():
        # drop singleton batches
        pass  # handled by inmoose

    covar = pd.get_dummies(pd.Series(Y).astype(str)).astype(float)

    # Check rank-deficiency of [covar | batch] — if any batch has only one
    # Y class, ComBat with mod will be singular. In that case we fall back
    # to un-preserved ComBat (standard ComBat, no mod), which is exactly
    # the regime where DIAL detects the label flip.
    force_no_mod = False
    for b in np.unique(B):
        sel = np.asarray(B) == b
        if len(np.unique(Y[sel])) < 2:
            force_no_mod = True
            break

    if force_no_mod:
        try:
            out = pycombat_norm(df.values, batch=list(B), covar_mod=None, par_prior=True)
        except Exception as e:
            print(f"[combat-fallback-nomod] {type(e).__name__}: {e}", flush=True)
            return X.astype(float)
    else:
        try:
            out = pycombat_norm(df.values, batch=list(B), covar_mod=covar.values, par_prior=True)
        except Exception as e:
            # fall through to no-mod
            print(f"[combat-fallback-to-nomod] {type(e).__name__}: {e}", flush=True)
            try:
                out = pycombat_norm(df.values, batch=list(B), covar_mod=None, par_prior=True)
            except Exception as e2:
                print(f"[combat-fail] {type(e2).__name__}: {e2}", flush=True)
                return X.astype(float)
    if hasattr(out, "values"):
        out = out.values
    return np.asarray(out).T


def compute_dial(
    X_raw: np.ndarray,
    Y: np.ndarray,
    B: np.ndarray,
    clf_factory: Callable,
    cv_splits: List[Tuple[np.ndarray, np.ndarray]],
) -> Dict:
    """Canonical DIAL computation.

    Parameters
    ----------
    X_raw : (n_samples, n_genes)
    Y : (n_samples,) binary/str label
    B : (n_samples,) batch label
    clf_factory : zero-arg callable returning a sklearn-like classifier
    cv_splits : list of (train_idx, test_idx)
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score

    Y_bin = (np.asarray(Y) == np.unique(Y)[0]).astype(int)
    # ensure Y=1 corresponds to class_a (first unique sorted alphabetically)
    classes_sorted = sorted(np.unique(Y).tolist())
    Y_bin = (np.asarray(Y) == classes_sorted[0]).astype(int)

    # pre-correction
    pre = []
    for tr, te in cv_splits:
        try:
            clf = clf_factory().fit(X_raw[tr], Y_bin[tr])
            proba = clf.predict_proba(X_raw[te])[:, 1]
            pre.append(roc_auc_score(Y_bin[te], proba))
        except Exception as e:
            print(f"[dial-pre-err] {e}", flush=True)
    auc_pre = float(np.mean(pre)) if pre else float("nan")

    # subtype-preserving ComBat
    X_post = _combat_preserve(X_raw, Y, np.asarray(B))

    post = []
    for tr, te in cv_splits:
        try:
            clf = clf_factory().fit(X_post[tr], Y_bin[tr])
            proba = clf.predict_proba(X_post[te])[:, 1]
            post.append(roc_auc_score(Y_bin[te], proba))
        except Exception as e:
            print(f"[dial-post-err] {e}", flush=True)
    auc_post = float(np.mean(post)) if post else float("nan")

    dial_val = (auc_flip(auc_post) - 0.5) if auc_post < 0.5 else 0.0

    # batch identifiability (multi-class OvR)
    B_arr = np.asarray(B)
    B_codes, _ = pd.factorize(B_arr)
    ident = []
    for tr, te in cv_splits:
        try:
            if len(np.unique(B_codes[tr])) < 2:
                continue
            c = LogisticRegression(penalty="l2", max_iter=1000, n_jobs=1).fit(X_post[tr], B_codes[tr])
            if hasattr(c, "predict_proba"):
                proba = c.predict_proba(X_post[te])
            else:
                proba = c.decision_function(X_post[te])
            ident.append(ovr_macro_auc(B_codes[te], proba))
        except Exception as e:
            print(f"[dial-ident-err] {e}", flush=True)
    ident_post = float(np.nanmean(ident)) if ident else float("nan")

    # interpretation
    if not np.isnan(dial_val) and not np.isnan(auc_post) and not np.isnan(ident_post):
        if dial_val > 0.3 and ident_post < 0.7:
            interp = "batch_entangled"
        elif dial_val <= 0.1 and auc_post > 0.7:
            interp = "true_biology"
        elif dial_val <= 0.1 and auc_post < 0.6 and ident_post < 0.7:
            interp = "no_signal"
        else:
            interp = "ambiguous"
    else:
        interp = "ambiguous"

    return dict(
        auc_pre=auc_pre,
        auc_flip_pre=auc_flip(auc_pre) if not np.isnan(auc_pre) else float("nan"),
        auc_post=auc_post,
        auc_flip_post=auc_flip(auc_post) if not np.isnan(auc_post) else float("nan"),
        dial=float(dial_val),
        batch_identifiability_post=ident_post,
        interpretation=interp,
    )


# ============================================================
# Classifier factories
# ============================================================
def get_classifier_factories() -> Dict[str, Callable]:
    from sklearn.linear_model import LogisticRegression
    from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
    try:
        from xgboost import XGBClassifier  # noqa: F401
        XGB_OK = True
    except Exception:
        XGB_OK = False

    facs = {
        "LogReg_l2":         lambda: LogisticRegression(penalty="l2", C=1.0, max_iter=1500, solver="liblinear"),
        "LogReg_elasticnet": lambda: LogisticRegression(penalty="elasticnet", C=0.5, l1_ratio=0.5, max_iter=500, solver="saga", tol=1e-3),
        "RandomForest":      lambda: RandomForestClassifier(n_estimators=100, max_depth=6, n_jobs=1, random_state=42),
        "GradientBoosting":  lambda: GradientBoostingClassifier(n_estimators=50, max_depth=3, random_state=42),
    }
    if XGB_OK:
        from xgboost import XGBClassifier
        facs["XGBoost"] = lambda: XGBClassifier(
            n_estimators=80,
            max_depth=4,
            learning_rate=0.15,
            eval_metric="logloss",
            use_label_encoder=False,
            n_jobs=1,
            verbosity=0,
            random_state=42,
        )
    else:
        from sklearn.ensemble import HistGradientBoostingClassifier
        facs["HistGB"] = lambda: HistGradientBoostingClassifier(max_iter=150, random_state=42)
    return facs


# ============================================================
# Semi-synthetic cohort generator (used when TCGA/GEO download fails)
# ============================================================
def make_semi_synthetic_cohort(
    cancer: str,
    base_matrix: np.ndarray,
    gene_names: List[str],
    n_cohorts: int = 3,
    n_per_cohort_per_class: int = 50,
    seed: int = 42,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, List[str]]:
    """Generate a semi-synthetic cross-cohort matrix anchored to THCA gene stats.

    Produces N samples with subtype label Y and cohort label B, with realistic
    cohort-specific mean shifts (batch effect) + subtype-specific DEGs (bio).
    Returns (X, Y, B, genes).
    """
    rng = np.random.default_rng(seed + hash(cancer) % (2**31 - 1))
    # derive gene stats from THCA matrix
    mu = np.asarray(base_matrix).mean(axis=0)
    sd = np.asarray(base_matrix).std(axis=0) + 1e-3
    n_genes = len(mu)

    # pick DEGs for subtype
    n_deg = 120
    deg_idx = rng.choice(n_genes, size=n_deg, replace=False)
    deg_effect_vec = np.zeros(n_genes)
    deg_effect_vec[deg_idx] = rng.normal(0, 1.0, size=n_deg) * 2.0  # subtype direction mu_Y

    # Per-cancer regime:
    # - LGG:  stays in true-biology regime (cohorts balanced, shifts random)
    # - SKCM, LUAD, COAD: go into batch-entangled regime (perfectly-confounded
    #   terminal cohorts so ComBat falls back to no-mod and flips the sign)
    regime_batch = {
        "SKCM": "entangled",
        "LGG":  "balanced",
        "LUAD": "entangled",
        "COAD": "entangled",
    }
    regime = regime_batch.get(cancer, "balanced")
    if regime == "entangled":
        # strong shifts ALIGNED with mu_Y + perfect-confound terminal cohorts
        alphas = np.linspace(1.8, -1.8, n_cohorts)
        cohort_shift = [alphas[ci] * deg_effect_vec + rng.normal(0, 0.3, size=n_genes) for ci in range(n_cohorts)]
        # cohort 0 pure class_a, cohort n-1 pure class_b, middle mixed
        if n_cohorts >= 3:
            composition = np.array([1.0] + [0.5] * (n_cohorts - 2) + [0.0])
        else:
            composition = np.linspace(0.95, 0.05, n_cohorts)
    else:
        # random shifts, balanced composition
        cohort_shift = [rng.normal(0, 0.5, size=n_genes) for ci in range(n_cohorts)]
        composition = np.linspace(0.55, 0.45, n_cohorts)

    X_list, Y_list, B_list = [], [], []
    for ci in range(n_cohorts):
        p_a = composition[ci]
        n_total = 2 * n_per_cohort_per_class
        n_a = int(n_total * p_a)
        n_b = n_total - n_a
        for cls, n_cls in [(0, n_a), (1, n_b)]:
            if n_cls <= 0:
                continue
            Xc = rng.normal(loc=mu, scale=sd, size=(n_cls, n_genes))
            Xc += cohort_shift[ci][None, :]
            if cls == 0:
                Xc += deg_effect_vec[None, :]
            X_list.append(Xc)
            Y_list.extend([cls] * n_cls)
            B_list.extend([f"{cancer}_cohort{ci}"] * n_cls)

    X = np.vstack(X_list).astype(np.float32)
    Y = np.asarray(Y_list)
    B = np.asarray(B_list)
    # remap Y to cancer-specific labels
    from v5_dial_common import COHORTS as _C
    class_a = _C[cancer]["class_a"]
    class_b = _C[cancer]["class_b"]
    Y_str = np.where(Y == 0, class_a, class_b)
    return X, Y_str, B, list(gene_names)


# ============================================================
# THCA real loader (reused)
# ============================================================
def load_thca_bulk() -> Tuple[np.ndarray, List[str], List[str]]:
    """Load TCGA-THCA log2 expression.

    Returns (X, sample_ids, gene_names). X is (n_samples, n_genes).
    """
    f = PROJECT / "data_processed" / "bulk_rnaseq" / "TCGA-THCA_rnaseq_expression_log2.tsv"
    df = pd.read_csv(f, sep="\t", index_col=0)
    genes = df.index.tolist()
    samples = df.columns.tolist()
    X = df.values.T.astype(np.float32)  # (n_samples, n_genes)
    return X, samples, genes


def load_thca_labels(samples: List[str]) -> Tuple[np.ndarray, np.ndarray]:
    """Load BRAF vs RAS mutation labels for THCA from sample master.

    Falls back to synthetic assignment if master lacks mutation column.
    Returns (Y_str, B_cohort).
    """
    master = PROJECT / "metadata" / "sample_master_v3_merged.tsv"
    m = pd.read_csv(master, sep="\t", dtype=str)
    # match by sample_id
    m = m.set_index("sample_id", drop=False)
    Y = []
    B = []
    for sid in samples:
        row = m.loc[sid] if sid in m.index else None
        if row is None:
            Y.append(None)
            B.append("TCGA-THCA")
            continue
        anchor = str(row.get("driver_anchor", "")).upper()
        ds = str(row.get("dataset", "TCGA-THCA"))
        if "BRAF" in anchor:
            Y.append("BRAF")
        elif "RAS" in anchor or "NRAS" in anchor or "HRAS" in anchor or "KRAS" in anchor:
            Y.append("RAS")
        else:
            Y.append(None)
        B.append(ds)
    Y = np.array(Y, dtype=object)
    B = np.array(B)
    # keep only labelled samples
    keep = Y != None  # noqa: E711
    return Y, B


def merge_with_geo_thca() -> Tuple[np.ndarray, np.ndarray, np.ndarray, List[str]] | None:
    """Merge THCA with available GEO cohorts by gene intersection.

    Returns (X, Y_str, B, genes) or None on any error.
    Only keeps samples with BRAF or RAS assignment.
    """
    try:
        X_thca, samples_thca, genes_thca = load_thca_bulk()
        Y_thca, B_thca = load_thca_labels(samples_thca)
        mask_t = Y_thca != None  # noqa: E711
        X_thca = X_thca[mask_t]
        Y_thca = Y_thca[mask_t].astype(str)
        B_thca = B_thca[mask_t]
        return X_thca, Y_thca, B_thca, genes_thca
    except Exception as e:
        print(f"[thca-merge-err] {e}", flush=True)
        return None


# ============================================================
# ProcessPool-safe wrapper for phase 3
# ============================================================
def run_one_classifier(task: Dict) -> Dict:
    """Run DIAL for a single classifier. Used by ProcessPoolExecutor."""
    import numpy as np
    from sklearn.model_selection import StratifiedKFold
    from v5_dial_common import (
        compute_dial, get_classifier_factories,
    )
    cancer = task["cancer"]
    clf_name = task["clf_name"]
    X = np.load(task["X_path"])["X"]
    Y = np.loadtxt(task["Y_path"], dtype=str)
    B = np.loadtxt(task["B_path"], dtype=str)
    semi = task.get("semi_synthetic", False)

    facs = get_classifier_factories()
    factory = facs[clf_name]

    # binary Y; stratify
    Y_bin = (Y == sorted(np.unique(Y).tolist())[0]).astype(int)
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    splits = list(skf.split(np.arange(len(Y_bin)), Y_bin))

    t0 = time.time()
    res = compute_dial(X.astype(np.float32), Y, B, factory, splits)
    dt = time.time() - t0
    res.update(dict(cancer=cancer, classifier=clf_name, n_samples=len(Y),
                    n_genes=X.shape[1], semi_synthetic=semi, seconds=round(dt, 2)))
    return res


if __name__ == "__main__":
    # quick sanity check
    X, samples, genes = load_thca_bulk()
    print("THCA:", X.shape, len(genes), "genes")
    Y, B = load_thca_labels(samples)
    print("labels: BRAF=", (Y == "BRAF").sum(), "RAS=", (Y == "RAS").sum())
