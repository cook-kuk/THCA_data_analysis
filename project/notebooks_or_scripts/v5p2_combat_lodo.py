#!/usr/bin/env python3
"""v5.2 Task A1 — Proper LODO ComBat with fit/transform separation.

v5.1 applied ComBat to full X (train+test pooled) BEFORE the LODO split,
leaking test-batch information into ComBat parameter estimation.

This module refits a covariate-preserving location-scale ComBat on the
TRAIN cohorts of each LODO fold, then transforms train and test with
the train-derived parameters. For the held-out cohort (unseen batch),
γ and δ are estimated locally from the test cohort alone — the
conservative choice documented in the v5.2 spec.

Model (per gene g):
    X_{gi} = α_g + y_i^T β_g + γ_{g, b(i)} + δ_{g, b(i)} · ε_{gi}

Fit: OLS per gene with design [intercept | Y-onehot | batch-onehot].
     α_g, β_g, γ_{gb} read off coefficients; δ_{gb} = residual SD per batch.

Transform: for each sample i in batch b
    if b ∈ train_batches:  γ̂, δ̂ from fit
    else (held-out):        γ̂ = X̄_b − α_g (no y info — unsupervised);
                            δ̂ = SD over test batch + floor
    X̂_{gi} = α_g + y_i^T β_g + (X_{gi} − α_g − y_i^T β_g − γ̂) · (σ_train / δ̂)

Output: corrected matrix where batch location-scale is removed while the
biological (y) signal is preserved via β_g re-injection on train samples.
Test samples cannot re-inject β_g·y (y is what we predict), so we only
subtract batch shift/scale and rely on β_g·y being encoded in X.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


class LinearComBat:
    def __init__(self, eps: float = 1e-6):
        self.eps = eps
        self.alpha_g: np.ndarray | None = None      # (G,)
        self.beta_g: np.ndarray | None = None       # (G, K_cov)
        self.gamma: dict = {}                        # batch -> (G,)
        self.delta: dict = {}                        # batch -> (G,)
        self.sigma_train: np.ndarray | None = None  # (G,) overall residual SD target
        self.y_levels_: list | None = None
        self.train_batches_: set = set()

    @staticmethod
    def _onehot(v: np.ndarray, levels: list) -> np.ndarray:
        mat = np.zeros((len(v), len(levels)), dtype=np.float64)
        for k, lv in enumerate(levels):
            mat[:, k] = (np.asarray(v) == lv).astype(np.float64)
        return mat

    def fit(self, X: np.ndarray, y: np.ndarray, batch: np.ndarray) -> "LinearComBat":
        X = np.asarray(X, dtype=np.float64)
        y = np.asarray(y)
        batch = np.asarray(batch)
        N, G = X.shape

        y_levels = sorted(np.unique(y).tolist())
        b_levels = sorted(np.unique(batch).tolist())
        self.y_levels_ = y_levels
        self.train_batches_ = set(b_levels)

        # Drop reference level for y and batch to avoid collinearity with intercept.
        # In the single-batch case (LODO 2-cohort fold), B_oh is empty (k_b=0).
        Y_oh = self._onehot(y, y_levels)[:, 1:]           # (N, K_y-1)
        if len(b_levels) >= 2:
            B_oh = self._onehot(batch, b_levels)[:, 1:]   # (N, K_b-1)
        else:
            B_oh = np.zeros((N, 0))
        intercept = np.ones((N, 1), dtype=np.float64)
        D = np.hstack([intercept, Y_oh, B_oh])             # (N, P)
        k_y = Y_oh.shape[1]
        k_b = B_oh.shape[1]

        # Solve OLS: coefs (P, G). lstsq is O(N^2 G) — fine for our sizes.
        coefs, *_ = np.linalg.lstsq(D, X, rcond=None)

        alpha_only = coefs[0:1, :]                         # (1, G) grand-intercept
        beta_cov = coefs[1:1 + k_y, :]                     # (k_y, G) covariate effects (non-ref)
        gamma_bat = coefs[1 + k_y:, :] if k_b > 0 else np.zeros((0, G))   # (k_b, G)

        # Build full per-level arrays (reference level = 0 coefficient).
        full_beta = np.vstack([np.zeros((1, G)), beta_cov])         # (K_y, G)
        full_gamma = np.vstack([np.zeros((1, G)), gamma_bat])       # (K_b, G)

        # Compute residuals and per-batch SD.
        fitted = D @ coefs
        resid = X - fitted
        delta = {}
        for k, b in enumerate(b_levels):
            mask = (batch == b)
            r = resid[mask]
            if r.shape[0] < 2:
                d = np.full(G, np.nan)
            else:
                d = r.std(axis=0, ddof=1)
            d = np.where(np.isfinite(d) & (d > self.eps), d, self.eps)
            delta[b] = d

        self.alpha_g = alpha_only.ravel()                   # (G,)
        # store beta_g as dict: y_level -> (G,)
        self.beta_g = {lv: full_beta[i] for i, lv in enumerate(y_levels)}
        self.gamma = {lv: full_gamma[i] for i, lv in enumerate(b_levels)}
        self.delta = {b: delta[b] for b in b_levels}

        # Target SD: overall residual SD across all train samples (avoid batch-conditional).
        self.sigma_train = np.where(
            np.isfinite(resid.std(axis=0, ddof=1)) & (resid.std(axis=0, ddof=1) > self.eps),
            resid.std(axis=0, ddof=1),
            1.0,
        )
        return self

    def transform(self, X: np.ndarray, batch: np.ndarray,
                  y: np.ndarray | None = None) -> np.ndarray:
        X = np.asarray(X, dtype=np.float64)
        batch = np.asarray(batch)
        if y is not None:
            y = np.asarray(y)
        N, G = X.shape
        out = np.empty_like(X)

        for b in np.unique(batch):
            mask = (batch == b)
            Xb = X[mask]
            if b in self.gamma:
                gamma_b = self.gamma[b]
                delta_b = self.delta[b]
            else:
                # Held-out cohort: no y info, local unsupervised estimate.
                mu_b = Xb.mean(axis=0)
                gamma_b = mu_b - self.alpha_g
                sd_b = Xb.std(axis=0, ddof=1) if Xb.shape[0] > 1 else np.ones(G)
                delta_b = np.where(np.isfinite(sd_b) & (sd_b > self.eps), sd_b, self.eps)

            if y is not None:
                # re-inject covariate effects on train samples (biology preservation)
                beta_part = np.zeros((Xb.shape[0], G))
                for i, yi in enumerate(y[mask]):
                    if yi in self.beta_g:
                        beta_part[i] = self.beta_g[yi]
                bio = self.alpha_g + beta_part
            else:
                bio = np.broadcast_to(self.alpha_g, Xb.shape)

            # z-standardize within batch (against batch's own γ,δ), then rescale to σ_train
            centered = Xb - self.alpha_g - gamma_b
            # Scale fix: divide by delta_b then multiply by sigma_train.
            out[mask] = bio + centered * (self.sigma_train / delta_b)
        return out


def compute_dial_lodo_proper(
    X_raw: np.ndarray,
    Y: np.ndarray,
    B: np.ndarray,
    clf_factory,
    cv_splits: list,
) -> dict:
    """Re-implementation of compute_dial with PROPER LODO ComBat.

    ComBat fit is done separately for each fold on train data only, then
    applied to train+test. Solves the v5.1 information-leak.
    """
    from sklearn.metrics import roc_auc_score
    from sklearn.linear_model import LogisticRegression

    classes_sorted = sorted(np.unique(Y).tolist())
    Y_bin = (np.asarray(Y) == classes_sorted[0]).astype(int)

    # ------- AUC pre (no ComBat) -------
    pre = []
    for tr, te in cv_splits:
        try:
            clf = clf_factory().fit(X_raw[tr], Y_bin[tr])
            proba = clf.predict_proba(X_raw[te])[:, 1]
            pre.append(roc_auc_score(Y_bin[te], proba))
        except Exception as e:
            print(f"[dial-pre-err] {e}", flush=True)
    auc_pre = float(np.mean(pre)) if pre else float("nan")

    # ------- AUC post (proper LODO ComBat) -------
    post = []
    post_per_fold = []   # (held_out_cohort, auc) for transparency
    B_arr = np.asarray(B)
    Y_arr = np.asarray(Y)
    for tr, te in cv_splits:
        try:
            held_out = np.unique(B_arr[te]).tolist()
            # Skip only if train has a single class (classifier cannot fit).
            # Single-batch train is fine: ComBat fit becomes intercept+covariate only,
            # and the held-out cohort is centered locally on transform.
            if len(np.unique(Y_bin[tr])) < 2:
                continue

            combat = LinearComBat().fit(X_raw[tr], Y_arr[tr], B_arr[tr])
            X_tr_adj = combat.transform(X_raw[tr], B_arr[tr], y=Y_arr[tr])
            X_te_adj = combat.transform(X_raw[te], B_arr[te], y=None)

            # Sanitize: replace any nan/inf produced by extreme scaling
            X_tr_adj = np.nan_to_num(X_tr_adj, nan=0.0, posinf=0.0, neginf=0.0)
            X_te_adj = np.nan_to_num(X_te_adj, nan=0.0, posinf=0.0, neginf=0.0)

            clf = clf_factory().fit(X_tr_adj, Y_bin[tr])
            proba = clf.predict_proba(X_te_adj)[:, 1]
            auc = roc_auc_score(Y_bin[te], proba)
            post.append(auc)
            post_per_fold.append((",".join(map(str, held_out)), float(auc)))
        except Exception as e:
            print(f"[dial-post-err] {e}", flush=True)
    auc_post = float(np.mean(post)) if post else float("nan")

    # ------- DIAL flip magnitude -------
    def auc_flip(a: float) -> float:
        return max(a, 1.0 - a)

    dial_val = (auc_flip(auc_post) - 0.5) if (not np.isnan(auc_post) and auc_post < 0.5) else 0.0

    # ------- Batch identifiability post (reuse pattern from v5.1) -------
    b_codes, _ = pd.factorize(B_arr)
    ident = []
    # Apply ComBat per fold for identifiability too (keep methodology consistent)
    for tr, te in cv_splits:
        try:
            if len(np.unique(b_codes[tr])) < 2:
                continue   # batch-identifiability needs >=2 train batches by definition
            if len(np.unique(Y_bin[tr])) < 2:
                continue
            combat = LinearComBat().fit(X_raw[tr], Y_arr[tr], B_arr[tr])
            X_tr_adj = combat.transform(X_raw[tr], B_arr[tr], y=Y_arr[tr])
            X_te_adj = combat.transform(X_raw[te], B_arr[te], y=None)
            X_tr_adj = np.nan_to_num(X_tr_adj)
            X_te_adj = np.nan_to_num(X_te_adj)
            c = LogisticRegression(penalty="l2", max_iter=1000, n_jobs=1).fit(X_tr_adj, b_codes[tr])
            if hasattr(c, "predict_proba"):
                p = c.predict_proba(X_te_adj)
            else:
                p = c.decision_function(X_te_adj)
            # reuse ovr_macro_auc semantics from common
            try:
                from sklearn.metrics import roc_auc_score
                uniq = np.unique(b_codes[te])
                if len(uniq) < 2:
                    continue
                if p.ndim == 1 or p.shape[1] == 1:
                    ident.append(float(roc_auc_score(b_codes[te], p)))
                elif len(uniq) == 2 and p.shape[1] >= 2:
                    ident.append(float(roc_auc_score(b_codes[te], p[:, 1])))
                else:
                    ident.append(float(roc_auc_score(b_codes[te], p, multi_class="ovr", average="macro")))
            except Exception:
                continue
        except Exception as e:
            print(f"[dial-ident-err] {e}", flush=True)
    ident_post = float(np.nanmean(ident)) if ident else float("nan")

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
        batch_identifiability_post=ident_post,
        interpretation=interp,
        per_fold_auc_post=post_per_fold,
    )


if __name__ == "__main__":
    # quick self-test
    rng = np.random.default_rng(0)
    N, G = 120, 50
    batch = np.array(["A"] * 60 + ["B"] * 60)
    y = np.array(["x"] * 30 + ["z"] * 30 + ["x"] * 30 + ["z"] * 30)
    X = rng.normal(0, 1, size=(N, G))
    # Inject batch shift + class signal
    X[batch == "A"] += 3.0
    X[:, 0] += (y == "x").astype(float) * 2.0
    cb = LinearComBat().fit(X, y, batch)
    Xc = cb.transform(X, batch, y=y)
    # held-out batch test
    batch2 = np.array(["C"] * 20)
    y2 = np.array(["x"] * 10 + ["z"] * 10)
    X2 = rng.normal(0, 1, size=(20, G)) + 5.0
    X2c = cb.transform(X2, batch2, y=None)
    print("X train means per batch BEFORE:", X[batch == "A"].mean(0)[:3], X[batch == "B"].mean(0)[:3])
    print("X train means per batch AFTER :", Xc[batch == "A"].mean(0)[:3], Xc[batch == "B"].mean(0)[:3])
    print("X held-out mean BEFORE:", X2.mean(0)[:3])
    print("X held-out mean AFTER :", X2c.mean(0)[:3])
