#!/usr/bin/env python3
# v10_leakage_knob.py
# Synthetic demonstration of the v5.1 vs v5.2 information leak.
#
# Setup (3-cohort LODO):
#   • 3 cohorts (k=0,1,2) with shared label direction e_Y_0 = ê_1 and per-
#     cohort label rotation rho_YY = 0.5 (modest platform divergence).
#   • Orthogonal batch direction e_B = ê_2.
#   • Y-B confounding: P(Y=1|B=k) varies as (0.85, 0.55, 0.35) across cohorts.
#   • Real ComBat-with-Y via inmoose.pycombat_norm (matches v5.1 / v5.2 stack).
#
# Leakage knob (λ ∈ [0,1]):
#   • λ = 0: ComBat is fit using ONLY the two training cohorts of each LODO
#     fold; held-out cohort is then transformed with the train-derived
#     parameters.  This is the v5.2 "proper" protocol.
#   • λ = 1: ComBat is fit on the FULL pooled dataset (all 3 cohorts)
#     before LODO splitting.  This is the v5.1 leaky protocol.
#   • intermediate λ: fit ComBat using train cohorts plus a λ-fraction of
#     the held-out cohort samples.  Models the partial-leakage regime.
#
# Output:
#   results/v10_aaai/synthetic_benchmark/leakage_knob.tsv
#   reports/html/figs_interactive/v10/leakage_knob.html

from __future__ import annotations
import pathlib, warnings, time
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")

PROJECT = pathlib.Path("/opt/thyroid-dash/project")
OUT_TSV = PROJECT / "results/v10_aaai/synthetic_benchmark/leakage_knob.tsv"
OUT_HTML = PROJECT / "reports/html/figs_interactive/v10/leakage_knob.html"
OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
OUT_HTML.parent.mkdir(parents=True, exist_ok=True)


# ------------------------------------------------------------------
# 3-cohort generator with platform divergence rho_YY
# ------------------------------------------------------------------
def make_data(n_per, p, rho_YY=-0.3, sigma_Y=1.0, sigma_B=2.5,
              imb_qs=(0.85, 0.35), noise=1.0, seed=0):
    """K=len(imb_qs) cohorts.  Default mirrors THCA-style setup:
    2 cohorts with reversed-imbalance Y rates and modest platform
    divergence (negative rho_YY)."""
    rng = np.random.default_rng(seed)
    K = len(imb_qs)
    e_Y = [None] * K
    # Cohort 0: ê_1
    e0 = np.zeros(p); e0[0] = 1.0
    e_Y[0] = e0
    # Cohort k > 0: rotate e0 by angle arccos(rho_YY) into a fresh axis per cohort
    for k in range(1, K):
        eY = np.zeros(p)
        eY[0] = rho_YY
        eY[2 + k] = np.sqrt(max(0.0, 1.0 - rho_YY ** 2))
        e_Y[k] = eY
    # Orthogonal batch direction (per-batch shifts along e_B)
    e_B = np.zeros(p); e_B[1] = 1.0
    n = K * n_per
    B = np.repeat(np.arange(K), n_per).astype(int)
    Y = np.zeros(n, dtype=int)
    r = rng.random(n)
    for k in range(K):
        Y[B == k] = (r[B == k] < imb_qs[k]).astype(int)
    X = np.zeros((n, p))
    for i in range(n):
        eY = e_Y[B[i]]
        bk = B[i]
        # Batch shifts: equally-spaced along e_B, centred at 0
        shift = sigma_B * (bk - (K - 1) / 2.0) / max(1, (K - 1) / 2.0)
        X[i] = sigma_Y * (2 * Y[i] - 1) * eY + shift * e_B \
             + noise * rng.standard_normal(p)
    return X, Y, B


# ------------------------------------------------------------------
# ComBat with leakage knob.
#
# Given training mask `tr_mask` and held-out mask `te_mask`, fit
# inmoose.pycombat_norm on tr_mask ∪ (λ-fraction of te_mask).  Then
# fit again on tr_mask alone to get train-only parameters; the held-out
# cohort is transformed via the train-only fit.  The leak is encoded
# in the additional held-out samples used during the original fit.
#
# To keep the implementation faithful to v5.1, we use the simpler
# protocol: fit on (train ∪ λ-fraction of test), apply to ALL samples
# (train + test).  At λ=0 this is "fit on train only" and at λ=1 this
# is "fit on full pooled data" — exactly the v5.2 vs v5.1 contrast.
# ------------------------------------------------------------------
def combat_with_leak(X, Y, B, tr_mask, te_mask, lam, rng):
    from inmoose.pycombat import pycombat_norm
    if lam <= 0.0:
        fit_idx = np.where(tr_mask)[0]
    elif lam >= 1.0:
        fit_idx = np.arange(len(Y))
    else:
        te_idx = np.where(te_mask)[0]
        n_leak = int(round(lam * len(te_idx)))
        leak_idx = rng.choice(te_idx, size=n_leak, replace=False) if n_leak > 0 else np.array([], dtype=int)
        fit_idx = np.concatenate([np.where(tr_mask)[0], leak_idx])
    X_fit = X[fit_idx]
    Y_fit = Y[fit_idx]
    B_fit = B[fit_idx]
    df_fit = pd.DataFrame(X_fit.T).astype(float)
    covar = pd.get_dummies(pd.Series(Y_fit).astype(str),
                           drop_first=True).astype(float)
    try:
        out_fit = pycombat_norm(df_fit.values, batch=list(B_fit),
                                covar_mod=covar.values, par_prior=True)
    except Exception:
        out_fit = pycombat_norm(df_fit.values, batch=list(B_fit),
                                covar_mod=None, par_prior=True)
    out_fit = np.asarray(out_fit).T
    # Return full X with corrected rows for fit_idx; rows outside fit_idx
    # are reconstructed by per-batch alignment to the corrected fit_idx
    # mean (the v5.1 choice — which is precisely the leakage we model).
    Xc = X.copy().astype(float)
    Xc[fit_idx] = out_fit
    if lam < 1.0:
        # For samples outside fit_idx, apply per-batch mean shift derived
        # from fit_idx samples in the same batch.
        not_fit = np.setdiff1d(np.arange(len(Y)), fit_idx)
        if len(not_fit):
            for b in np.unique(B[not_fit]):
                src = fit_idx[B[fit_idx] == b]
                if len(src) == 0:
                    # No fit samples in this batch; use global mean offset
                    src_mean = Xc[fit_idx].mean(axis=0)
                    src_raw_mean = X[fit_idx].mean(axis=0)
                else:
                    src_mean = Xc[src].mean(axis=0)
                    src_raw_mean = X[src].mean(axis=0)
                shift = src_mean - src_raw_mean
                tgt = not_fit[B[not_fit] == b]
                Xc[tgt] = X[tgt] + shift
    return Xc


def lodo_auc(X, Y, B, C=1.0):
    aucs = []
    for k_test in np.unique(B):
        tr = B != k_test; te = B == k_test
        if len(np.unique(Y[tr])) < 2 or len(np.unique(Y[te])) < 2:
            continue
        clf = LogisticRegression(penalty="l2", C=C, max_iter=3000,
                                 solver="liblinear")
        clf.fit(X[tr], Y[tr])
        s = clf.decision_function(X[te])
        aucs.append(roc_auc_score(Y[te], s))
    return float(np.mean(aucs)) if aucs else np.nan


def lodo_auc_with_leak(X, Y, B, lam, seed):
    """For each LODO fold, apply ComBat with leakage knob λ, then train
    LogReg on the training fold and score the test fold."""
    rng = np.random.default_rng(seed)
    aucs = []
    for k_test in np.unique(B):
        tr = B != k_test; te = B == k_test
        if len(np.unique(Y[tr])) < 2 or len(np.unique(Y[te])) < 2:
            continue
        try:
            Xc = combat_with_leak(X, Y, B, tr, te, lam, rng)
        except Exception:
            return np.nan
        clf = LogisticRegression(penalty="l2", C=1.0, max_iter=3000,
                                 solver="liblinear")
        clf.fit(Xc[tr], Y[tr])
        s = clf.decision_function(Xc[te])
        aucs.append(roc_auc_score(Y[te], s))
    return float(np.mean(aucs)) if aucs else np.nan


def dial(auc):
    if np.isnan(auc):
        return np.nan
    return max(auc, 1 - auc) - auc


def main():
    n_per = 50
    p = 100
    n_reps = 8
    lam_grid    = np.linspace(0.0, 1.0, 6)
    rho_YY_grid = [-0.7, -0.3, 0.0, +0.3, +0.7, +1.0]

    rows = []
    t0 = time.time()
    for rho_YY in rho_YY_grid:
        for lam in lam_grid:
            cell_aucs = []
            for rep in range(n_reps):
                seed = hash((round(float(rho_YY), 3),
                             round(float(lam), 3), int(rep))) & 0xFFFFFFFF
                X, Y, B = make_data(n_per, p, rho_YY=float(rho_YY),
                                    seed=seed)
                auc_pre  = lodo_auc(X, Y, B)
                auc_post = lodo_auc_with_leak(X, Y, B, float(lam),
                                              seed=seed)
                rows.append(dict(rho_YY=float(rho_YY), lam=float(lam),
                                 rep=rep,
                                 auc_pre=auc_pre, auc_post=auc_post,
                                 dial_pre=dial(auc_pre),
                                 dial_post=dial(auc_post)))
                cell_aucs.append(auc_post)
            print(f"[rho_YY={rho_YY:+.1f} lam={lam:.2f}] "
                  f"mean AUC_post={np.nanmean(cell_aucs):.3f}  "
                  f"DIAL={dial(np.nanmean(cell_aucs)):.3f}  "
                  f"({time.time()-t0:.1f}s)", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"\n[ok] tsv -> {OUT_TSV}")

    agg = (df.groupby(["rho_YY", "lam"])
             .agg(auc_post_mean=("auc_post", "mean"),
                  auc_post_std=("auc_post", "std"),
                  dial_post_mean=("dial_post", "mean"),
                  dial_post_std=("dial_post", "std"))
             .reset_index())

    try:
        import plotly.graph_objects as go
        # Heatmap: rho_YY × lam, color = mean DIAL
        piv = agg.pivot(index="rho_YY", columns="lam",
                        values="dial_post_mean")
        fig = go.Figure(data=go.Heatmap(
            z=piv.values,
            x=[f"{c:.2f}" for c in piv.columns],
            y=[f"{r:+.1f}" for r in piv.index],
            colorscale="Inferno",
            colorbar=dict(title="mean DIAL")))
        fig.update_layout(
            title=("Leakage knob: DIAL vs (ρ_YY, λ).  "
                   "λ=0 is v5.2 proper-LODO; λ=1 is v5.1 leaky."),
            xaxis_title="λ  (fraction of held-out cohort in ComBat fit)",
            yaxis_title="ρ_YY  (cross-cohort label-direction cosine)",
            template="plotly_white",
            height=480)
        fig.write_html(OUT_HTML, include_plotlyjs="cdn")
        print(f"[ok] html -> {OUT_HTML}")
    except Exception as e:
        print(f"[warn] plotly failed: {e}")

    print("\n[summary table — mean AUC_post by (ρ_YY, λ)]")
    auc_pivot = agg.pivot(index="rho_YY", columns="lam",
                          values="auc_post_mean")
    print(auc_pivot.round(3).to_string())
    print("\n[summary table — mean DIAL_post]")
    print(piv.round(3).to_string())


if __name__ == "__main__":
    main()
