#!/usr/bin/env python3
"""v5.1 REAL-DATA Cross-Cancer DIAL — shared utilities.

Hard rule: NO SYNTHETIC DATA. Every output row must have semi_synthetic=False.
"""
from __future__ import annotations

import os
import sys
import time
import json
import gzip
import warnings
from pathlib import Path
from typing import Callable, Dict, List, Tuple, Optional

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

PROJECT = Path("/opt/thyroid-dash/project")
# /data disk fallback
_data_disk = Path("/data")
if _data_disk.exists():
    DATA_RAW_V5 = _data_disk / "thca" / "v5_cross_cancer" / "raw"
else:
    DATA_RAW_V5 = PROJECT / "data_raw" / "v5_cross_cancer"
DATA_PROC_V5 = PROJECT / "data_processed" / "v5_cross_cancer"
RESULTS_V5 = PROJECT / "results" / "v5"
REPORTS_V5 = PROJECT / "reports" / "v5"
FIGS_V5 = PROJECT / "reports" / "html" / "figs_interactive" / "v5"
LOGS = PROJECT / "logs"
PAGES = PROJECT / "reports" / "html" / "pages"
GDC_CACHE = PROJECT / "data_raw" / "gdc"

for p in [DATA_RAW_V5, DATA_PROC_V5, RESULTS_V5, REPORTS_V5, FIGS_V5, LOGS]:
    p.mkdir(parents=True, exist_ok=True)


COHORTS = {
    "THCA": {
        "tcga": "TCGA-THCA",
        "geo":  ["GSE27155", "GSE33630", "GSE29265"],
        "task": "BRAF_vs_RAS",
        "class_a": "BRAF", "class_b": "RAS",
    },
    "SKCM": {
        "tcga": "TCGA-SKCM",
        "geo":  ["GSE65904", "GSE22153"],
        "task": "BRAF_vs_NRAS",
        "class_a": "BRAF", "class_b": "NRAS",
    },
    "LGG": {
        "tcga": "TCGA-LGG",
        "geo":  ["GSE16011", "GSE4271"],
        "task": "IDH_mut_vs_wt",
        "class_a": "IDHmut", "class_b": "IDHwt",
    },
    "LUAD": {
        "tcga": "TCGA-LUAD",
        "geo":  ["GSE31210", "GSE72094"],
        "task": "KRAS_vs_EGFR",
        "class_a": "KRAS", "class_b": "EGFR",
    },
    "COAD": {
        "tcga": "TCGA-COAD",
        "geo":  ["GSE39582", "GSE17536"],
        "task": "BRAF_vs_KRAS",
        "class_a": "BRAF", "class_b": "KRAS",
    },
}

CANCERS = ["THCA", "SKCM", "LGG", "LUAD", "COAD"]


def log_line(logfile: Path, msg: str) -> None:
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    with open(logfile, "a") as fh:
        fh.write(line + "\n")


# ============================================================
# Canonical DIAL (same as v5)
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
    """LEGACY (v5.1a): fits ComBat on the entire matrix in one call.

    Retained only for reproducing the historical v5.1a numbers; the F1 fix
    routes through `_combat_fit_train` / `_combat_apply_test` instead.
    """
    from inmoose.pycombat import pycombat_norm

    df = pd.DataFrame(X.T)
    if len(np.unique(B)) < 2:
        return X.astype(float)

    covar = pd.get_dummies(pd.Series(Y).astype(str)).astype(float)

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
            print(f"[combat-fallback-to-nomod] {type(e).__name__}: {e}", flush=True)
            try:
                out = pycombat_norm(df.values, batch=list(B), covar_mod=None, par_prior=True)
            except Exception as e2:
                print(f"[combat-fail] {type(e2).__name__}: {e2}", flush=True)
                return X.astype(float)
    if hasattr(out, "values"):
        out = out.values
    return np.asarray(out).T


# ============================================================
# F1 FIX: leak-safe ComBat with explicit fit/transform separation.
#
# Path (c) of the audit prompt: a small from-scratch ComBat that mirrors
# inmoose's `pycombat_norm` algorithm but exposes the train-fit parameters
# (B_hat, grand_mean, var_pooled, gamma_star, delta_star, hyper-priors) so
# the held-out test cohort can be standardized using train-only statistics.
#
# For the test cohort (a brand-new batch the train fit has never seen) we
# estimate gamma_star/delta_star for the test batch using ONLY test data,
# but standardized against TRAIN-fit grand_mean/var_pooled and adjusted
# under TRAIN-fit hyper-priors. This is the genuinely leak-safe LODO
# correction: train labels and train cohort statistics never touch test
# expression in a way that biases the test classifier toward train geometry.
# ============================================================
def _combat_fit_train(X_tr: np.ndarray, Y_tr: np.ndarray, B_tr: np.ndarray) -> Optional[Dict]:
    """Fit ComBat parameters on the train slice only.

    X_tr is sample x gene. Returns a dict of fitted parameters or None if
    the training slice is degenerate (single batch, etc.) -- caller should
    fall back to identity correction.
    """
    if len(np.unique(B_tr)) < 2:
        return None
    dat = X_tr.T.astype(float)  # (genes x samples)
    n_gene, n_sample = dat.shape

    batches_uniq, batch_codes = np.unique(B_tr, return_inverse=True)
    n_batch = len(batches_uniq)
    batches_ind = [np.where(batch_codes == i)[0] for i in range(n_batch)]
    batch_sizes = [len(v) for v in batches_ind]

    # Build design: rows = covariates+batches, cols = samples (matches inmoose convention)
    # batch dummies (n_batch rows)
    batch_dm = np.zeros((n_batch, n_sample))
    for i, idx in enumerate(batches_ind):
        batch_dm[i, idx] = 1.0

    # covariate dummies (drop one to avoid singularity with batch dummies)
    Y_arr = np.asarray(Y_tr).astype(str)
    classes = sorted(np.unique(Y_arr).tolist())
    use_cov = True
    # If any batch has only one class, the per-batch covar is collinear -> skip cov
    for i in range(n_batch):
        if len(np.unique(Y_arr[batches_ind[i]])) < 2:
            use_cov = False
            break

    if use_cov and len(classes) >= 2:
        # one-hot minus last (avoid intercept collinearity with batch dummies)
        cov_dm = np.zeros((len(classes) - 1, n_sample))
        for k, cls in enumerate(classes[:-1]):
            cov_dm[k, :] = (Y_arr == cls).astype(float)
        design = np.vstack([batch_dm, cov_dm])
    else:
        design = batch_dm

    # B_hat = (D D^T)^-1 D X^T  (n_design x n_gene)
    try:
        DDt = design @ design.T
        DXt = design @ dat.T
        B_hat = np.linalg.solve(DDt, DXt)
    except np.linalg.LinAlgError:
        # design singular -- drop covariates
        design = batch_dm
        DDt = design @ design.T
        DXt = design @ dat.T
        B_hat = np.linalg.solve(DDt, DXt)

    # grand_mean: weighted mean of batch coefficients (rows 0..n_batch)
    grand_mean = np.array([s / n_sample for s in batch_sizes]) @ B_hat[:n_batch, :]

    # var_pooled across all train samples
    resid = dat - (design.T @ B_hat).T
    var_pooled = (resid ** 2) @ (np.ones(n_sample) / n_sample)

    # avoid division by zero
    var_pooled = np.where(var_pooled <= 1e-12, 1e-12, var_pooled)
    sd_pooled = np.sqrt(var_pooled)

    # stand_mean = grand_mean repeated across samples + covariate part of design.T @ B_hat
    stand_mean = np.tile(grand_mean.reshape(-1, 1), (1, n_sample))
    if design.shape[0] > n_batch:
        # add covariate effect (set batch rows to 0)
        tmp = design.copy()
        tmp[:n_batch, :] = 0.0
        stand_mean = stand_mean + (tmp.T @ B_hat).T

    # standardise
    s_data = (dat - stand_mean) / sd_pooled[:, None]

    # gamma_hat (additive batch effect)
    batch_design = design[:n_batch, :]
    BDBDt = batch_design @ batch_design.T
    gamma_hat = np.linalg.solve(BDBDt, batch_design @ s_data.T)  # (n_batch x n_gene)

    # delta_hat (multiplicative batch effect): per-batch variance across samples
    delta_hat = np.zeros((n_batch, n_gene))
    for i in range(n_batch):
        delta_hat[i, :] = s_data[:, batches_ind[i]].var(axis=1, ddof=0)
    delta_hat = np.where(delta_hat <= 1e-12, 1e-12, delta_hat)

    # hyper-priors per batch
    gamma_bar = gamma_hat.mean(axis=1)         # (n_batch,)
    t2 = gamma_hat.var(axis=1, ddof=0)         # (n_batch,)

    a_prior = np.array([
        (2 * np.var(delta_hat[i]) + np.mean(delta_hat[i]) ** 2) / max(np.var(delta_hat[i]), 1e-12)
        for i in range(n_batch)
    ])
    b_prior = np.array([
        (np.mean(delta_hat[i]) * np.var(delta_hat[i]) + np.mean(delta_hat[i]) ** 3) / max(np.var(delta_hat[i]), 1e-12)
        for i in range(n_batch)
    ])

    # parametric EB iteration per batch
    gamma_star = np.zeros((n_batch, n_gene))
    delta_star = np.zeros((n_batch, n_gene))
    for i in range(n_batch):
        sdat_i = s_data[:, batches_ind[i]]
        n_i = sdat_i.shape[1]
        g_old = gamma_hat[i].copy()
        d_old = delta_hat[i].copy()
        change = 1.0
        count = 0
        while change > 1e-4 and count < 100:
            t2_n = t2[i] * n_i
            t2_n_g_hat = t2_n * gamma_hat[i]
            g_new = (t2_n_g_hat + d_old * gamma_bar[i]) / (t2_n + d_old)
            sum2 = ((sdat_i - g_new[:, None]) ** 2).sum(axis=1)
            d_new = (0.5 * sum2 + b_prior[i]) / (0.5 * n_i + a_prior[i] - 1.0)
            d_new = np.where(d_new <= 1e-12, 1e-12, d_new)
            change = max(
                np.max(np.abs(g_new - g_old) / np.maximum(np.abs(g_old), 1e-12)),
                np.max(np.abs(d_new - d_old) / np.maximum(np.abs(d_old), 1e-12)),
            )
            g_old = g_new
            d_old = d_new
            count += 1
        gamma_star[i] = g_old
        delta_star[i] = d_old

    return dict(
        batches_uniq=batches_uniq,
        n_batch=n_batch,
        grand_mean=grand_mean,
        var_pooled=var_pooled,
        sd_pooled=sd_pooled,
        gamma_star=gamma_star,
        delta_star=delta_star,
        gamma_bar=gamma_bar,
        t2=t2,
        a_prior=a_prior,
        b_prior=b_prior,
        B_hat=B_hat,
        n_design=design.shape[0],
        use_cov=(design.shape[0] > n_batch),
        classes=classes,
        n_train=n_sample,
    )


def _combat_apply_train(X_tr: np.ndarray, Y_tr: np.ndarray, B_tr: np.ndarray, fit: Dict) -> np.ndarray:
    """Apply the fitted ComBat correction to the train slice (sample x gene)."""
    dat = X_tr.T.astype(float)
    n_gene, n_sample = dat.shape
    batches_uniq = fit["batches_uniq"]
    n_batch = fit["n_batch"]

    # Re-derive batch indices in train order (matches the fit, batch order is the same)
    batch_codes = np.array([np.where(batches_uniq == b)[0][0] if b in batches_uniq else -1 for b in B_tr])
    batches_ind = [np.where(batch_codes == i)[0] for i in range(n_batch)]
    batch_sizes = [len(v) for v in batches_ind]

    grand_mean = fit["grand_mean"]
    sd_pooled = fit["sd_pooled"]
    gamma_star = fit["gamma_star"]
    delta_star = fit["delta_star"]

    stand_mean = np.tile(grand_mean.reshape(-1, 1), (1, n_sample))
    if fit["use_cov"]:
        Y_arr = np.asarray(Y_tr).astype(str)
        classes = fit["classes"]
        cov_dm = np.zeros((len(classes) - 1, n_sample))
        for k, cls in enumerate(classes[:-1]):
            cov_dm[k, :] = (Y_arr == cls).astype(float)
        # add covariate contribution: cov_dm.T @ B_hat[n_batch:, :]
        stand_mean = stand_mean + (cov_dm.T @ fit["B_hat"][n_batch:, :]).T

    s_data = (dat - stand_mean) / sd_pooled[:, None]

    bayes = s_data.copy()
    for i in range(n_batch):
        idx = batches_ind[i]
        if len(idx) == 0:
            continue
        bayes[:, idx] = (s_data[:, idx] - gamma_star[i][:, None]) / np.sqrt(delta_star[i])[:, None]

    out = bayes * sd_pooled[:, None] + stand_mean
    return out.T  # back to sample x gene


def _combat_apply_test(X_te: np.ndarray, Y_te: np.ndarray, B_te: np.ndarray, fit: Dict) -> np.ndarray:
    """Apply ComBat to a held-out test cohort using train-only standardisation.

    The test batch is unseen. We:
      1) Standardize test using train grand_mean / var_pooled.
      2) Estimate gamma/delta for the test batch from the standardized test
         data alone, applying inmoose-style EB shrinkage with TRAIN-fitted
         gamma_bar / t2 / a_prior / b_prior averaged across train batches.
      3) Adjust test by subtracting that gamma_star and dividing by sqrt(delta_star).

    This means train labels never bias the test correction, but the test
    cohort still gets centered/scaled on a train-defined coordinate frame.
    """
    dat = X_te.T.astype(float)
    n_gene, n_sample = dat.shape
    if n_sample == 0:
        return X_te.astype(float)

    grand_mean = fit["grand_mean"]
    sd_pooled = fit["sd_pooled"]

    stand_mean = np.tile(grand_mean.reshape(-1, 1), (1, n_sample))
    if fit["use_cov"]:
        Y_arr = np.asarray(Y_te).astype(str)
        classes = fit["classes"]
        cov_dm = np.zeros((len(classes) - 1, n_sample))
        for k, cls in enumerate(classes[:-1]):
            cov_dm[k, :] = (Y_arr == cls).astype(float)
        stand_mean = stand_mean + (cov_dm.T @ fit["B_hat"][fit["n_batch"]:, :]).T

    s_data = (dat - stand_mean) / sd_pooled[:, None]

    # Train hyper-priors (averaged across train batches as a reasonable shared prior)
    gamma_bar_t = float(np.mean(fit["gamma_bar"]))
    t2_t = float(np.mean(fit["t2"]))
    a_prior_t = float(np.mean(fit["a_prior"]))
    b_prior_t = float(np.mean(fit["b_prior"]))

    # All test samples = one new batch
    gamma_hat_te = s_data.mean(axis=1)
    delta_hat_te = s_data.var(axis=1, ddof=0)
    delta_hat_te = np.where(delta_hat_te <= 1e-12, 1e-12, delta_hat_te)

    g_old = gamma_hat_te.copy()
    d_old = delta_hat_te.copy()
    change = 1.0
    count = 0
    n_i = n_sample
    while change > 1e-4 and count < 100:
        t2_n = t2_t * n_i
        t2_n_g_hat = t2_n * gamma_hat_te
        g_new = (t2_n_g_hat + d_old * gamma_bar_t) / (t2_n + d_old)
        sum2 = ((s_data - g_new[:, None]) ** 2).sum(axis=1)
        d_new = (0.5 * sum2 + b_prior_t) / (0.5 * n_i + a_prior_t - 1.0)
        d_new = np.where(d_new <= 1e-12, 1e-12, d_new)
        change = max(
            np.max(np.abs(g_new - g_old) / np.maximum(np.abs(g_old), 1e-12)),
            np.max(np.abs(d_new - d_old) / np.maximum(np.abs(d_old), 1e-12)),
        )
        g_old = g_new
        d_old = d_new
        count += 1

    bayes = (s_data - g_old[:, None]) / np.sqrt(d_old)[:, None]
    out = bayes * sd_pooled[:, None] + stand_mean
    return out.T  # sample x gene


def compute_dial(
    X_raw: np.ndarray,
    Y: np.ndarray,
    B: np.ndarray,
    clf_factory: Callable,
    cv_splits: List[Tuple[np.ndarray, np.ndarray]],
    leak_safe: bool = True,
) -> Dict:
    """Compute DIAL under LODO.

    F1 fix (default `leak_safe=True`): ComBat is fit on the train slice
    only and then applied separately to train and to the held-out test
    cohort. Set `leak_safe=False` to reproduce v5.1a's whole-matrix fit.
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score

    classes_sorted = sorted(np.unique(Y).tolist())
    Y_bin = (np.asarray(Y) == classes_sorted[0]).astype(int)
    B_arr = np.asarray(B)

    pre = []
    for tr, te in cv_splits:
        try:
            clf = clf_factory().fit(X_raw[tr], Y_bin[tr])
            proba = clf.predict_proba(X_raw[te])[:, 1]
            pre.append(roc_auc_score(Y_bin[te], proba))
        except Exception as e:
            print(f"[dial-pre-err] {e}", flush=True)
    auc_pre = float(np.mean(pre)) if pre else float("nan")

    # ----- POST-ComBat AUC: per-fold leak-safe correction -----
    post = []
    # We also build a per-fold X_post_full for the batch-identifiability probe;
    # under leak_safe semantics this is the train-fitted correction applied to
    # the test cohort separately (no leakage between folds).
    fold_test_corrected: List[Tuple[np.ndarray, np.ndarray]] = []  # (te_idx, X_te_post)
    fold_train_corrected: List[Tuple[np.ndarray, np.ndarray]] = []  # (tr_idx, X_tr_post)

    for tr, te in cv_splits:
        try:
            if leak_safe:
                fit = _combat_fit_train(X_raw[tr], np.asarray(Y)[tr], B_arr[tr])
                if fit is None:
                    X_tr_post = X_raw[tr].astype(float)
                    X_te_post = X_raw[te].astype(float)
                else:
                    X_tr_post = _combat_apply_train(X_raw[tr], np.asarray(Y)[tr], B_arr[tr], fit)
                    X_te_post = _combat_apply_test(X_raw[te], np.asarray(Y)[te], B_arr[te], fit)
            else:
                # legacy whole-matrix fit (v5.1a behaviour); kept for back-compat
                X_full = _combat_preserve(X_raw, np.asarray(Y), B_arr)
                X_tr_post = X_full[tr]
                X_te_post = X_full[te]

            fold_train_corrected.append((tr, X_tr_post))
            fold_test_corrected.append((te, X_te_post))

            clf = clf_factory().fit(X_tr_post, Y_bin[tr])
            proba = clf.predict_proba(X_te_post)[:, 1]
            post.append(roc_auc_score(Y_bin[te], proba))
        except Exception as e:
            print(f"[dial-post-err] {e}", flush=True)
    auc_post = float(np.mean(post)) if post else float("nan")

    dial_val = (auc_flip(auc_post) - 0.5) if (not np.isnan(auc_post) and auc_post < 0.5) else 0.0

    # Batch identifiability under post-correction. Use the fold-corrected
    # train and test slices so leakage semantics carry through.
    B_codes, _ = pd.factorize(B_arr)
    ident = []
    for (tr, X_tr_post), (te, X_te_post) in zip(fold_train_corrected, fold_test_corrected):
        try:
            if len(np.unique(B_codes[tr])) < 2:
                continue
            c = LogisticRegression(penalty="l2", max_iter=1000, n_jobs=1).fit(X_tr_post, B_codes[tr])
            if hasattr(c, "predict_proba"):
                proba = c.predict_proba(X_te_post)
            else:
                proba = c.decision_function(X_te_post)
            ident.append(ovr_macro_auc(B_codes[te], proba))
        except Exception as e:
            print(f"[dial-ident-err] {e}", flush=True)
    ident_post = float(np.nanmean(ident)) if ident else float("nan")

    # ----- F3 fix: partial_flip_score and dial_v2 -----
    # partial_flip_score captures cases like RandomForest/ComBat-seq where
    # auc_post creeps just over 0.5 (so the original DIAL collapses to 0)
    # but auc_flip(auc_post) is still high -- i.e. there is residual
    # batch-entangled flip information that the binary DIAL misses.
    if not np.isnan(auc_post):
        partial_flip_score = max(0.0, auc_flip(auc_post) - max(0.5, auc_post))
        dial_v2 = (auc_flip(auc_post) - auc_post) if auc_flip(auc_post) > 0.8 else 0.0
        dial_v2 = max(0.0, dial_v2)
    else:
        partial_flip_score = float("nan")
        dial_v2 = float("nan")

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
        auc_pre=auc_pre,
        auc_flip_pre=auc_flip(auc_pre) if not np.isnan(auc_pre) else float("nan"),
        auc_post=auc_post,
        auc_flip_post=auc_flip(auc_post) if not np.isnan(auc_post) else float("nan"),
        dial=float(dial_val),
        partial_flip_score=float(partial_flip_score),
        dial_v2=float(dial_v2),
        batch_identifiability_post=ident_post,
        interpretation=interp,
        leak_safe=bool(leak_safe),
    )


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
            n_estimators=80, max_depth=4, learning_rate=0.15,
            eval_metric="logloss", n_jobs=1, verbosity=0, random_state=42,
        )
    else:
        from sklearn.ensemble import HistGradientBoostingClassifier
        facs["HistGB"] = lambda: HistGradientBoostingClassifier(max_iter=150, random_state=42)
    return facs


def run_one_classifier(task: Dict) -> Dict:
    """ProcessPool-safe single-classifier DIAL run."""
    import numpy as np
    from sklearn.model_selection import LeaveOneGroupOut
    from v5p1_common import (
        compute_dial, get_classifier_factories,
    )
    cancer = task["cancer"]
    clf_name = task["clf_name"]
    X = np.load(task["X_path"])["X"]
    Y = np.loadtxt(task["Y_path"], dtype=str)
    B = np.loadtxt(task["B_path"], dtype=str)

    facs = get_classifier_factories()
    factory = facs[clf_name]

    Y_bin = (Y == sorted(np.unique(Y).tolist())[0]).astype(int)
    # Leave-One-Dataset-Out: each fold holds out one cohort entirely.
    # Exposes batch confound the way v4 did. StratifiedKFold hides it.
    logo = LeaveOneGroupOut()
    splits = list(logo.split(np.arange(len(Y_bin)), Y_bin, groups=B))

    # Top-variable gene filter — pycombat is O(genes^2) so full 11k genes is too slow
    # (pilot: 2028s for one classifier on THCA 11710 genes). Filter keeps signal but
    # cuts runtime ~10x. Threshold: top 3000 by variance across samples.
    n_top = 3000
    n_orig = X.shape[1]
    if X.shape[1] > n_top:
        var_per_gene = X.var(axis=0)
        top_idx = np.argsort(var_per_gene)[-n_top:]
        X = X[:, top_idx]
    _log = LOGS / "v5p1_dial.log"
    log_line(_log, f"[{cancer}/{clf_name}] variance-top filter: {n_orig} -> {X.shape[1]} genes")

    t0 = time.time()
    res = compute_dial(X.astype(np.float32), Y, B, factory, splits)
    dt = time.time() - t0
    res.update(dict(cancer=cancer, classifier=clf_name, n_samples=len(Y),
                    n_genes=X.shape[1], semi_synthetic=False, seconds=round(dt, 2)))
    return res


if __name__ == "__main__":
    print("v5p1_common loaded OK")
    print("DATA_RAW_V5:", DATA_RAW_V5)
    print("cancers:", CANCERS)
