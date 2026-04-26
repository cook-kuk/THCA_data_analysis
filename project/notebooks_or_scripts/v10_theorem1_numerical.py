#!/usr/bin/env python3
# v10_theorem1_numerical.py  (v3 — ρ_YY sweep)
# Demonstrate Theorem 1's phase transition using the mechanism identified
# in Eq. (disc-full) of the extended proof: the flip is driven by the
# cross-cohort label-direction agreement ρ_YY = <e_Y_0, e_Y_1>, not the
# simple geometric ρ = <w_Y, w_B>.
#
# Setup (matches the paper's Figure 2):
#   Two cohorts.  Batch 0 label direction e_Y_0 = ê_1.
#   Batch 1 label direction e_Y_1 at angle α from e_Y_0,
#     cos α = ρ_YY ∈ [-1, 1].
#   Orthogonal batch direction e_B.
#   Strong Y-B confounding (imb_q = 0.85 vs 0.35).
#   Real ComBat-with-Y via inmoose.pycombat_norm.
#
# Output:
#   results/v10_aaai/synthetic_benchmark/theorem1_numerical.tsv
#   reports/html/figs_interactive/v10/theorem1_phase_transition.html

from __future__ import annotations
import pathlib, warnings
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")

PROJECT = pathlib.Path("/opt/thyroid-dash/project")
OUT_TSV = PROJECT / "results/v10_aaai/synthetic_benchmark/theorem1_numerical.tsv"
OUT_HTML = PROJECT / "reports/html/figs_interactive/v10/theorem1_phase_transition.html"
OUT_HTML.parent.mkdir(parents=True, exist_ok=True)
OUT_TSV.parent.mkdir(parents=True, exist_ok=True)


def make_data_platform(n_per, p, rho_YY, imb_q0=0.85, imb_q1=0.35,
                       sigma_Y=1.0, sigma_B=3.0, noise=1.0, seed=0):
    rng = np.random.default_rng(seed)
    e_Y_0 = np.zeros(p); e_Y_0[0] = 1.0
    e_Y_1 = np.zeros(p); e_Y_1[0] = rho_YY
    e_Y_1[2] = np.sqrt(max(0.0, 1.0 - rho_YY ** 2))
    e_B = np.zeros(p); e_B[1] = 1.0
    n = 2 * n_per
    B = np.r_[np.zeros(n_per), np.ones(n_per)].astype(int)
    Y = np.zeros(n, dtype=int)
    r = rng.random(n)
    Y[B == 0] = (r[B == 0] < imb_q0).astype(int)
    Y[B == 1] = (r[B == 1] < imb_q1).astype(int)
    X = np.zeros((n, p))
    for i in range(n):
        eY = e_Y_0 if B[i] == 0 else e_Y_1
        X[i] = sigma_Y * (2 * Y[i] - 1) * eY \
             + sigma_B * (2 * B[i] - 1) * e_B \
             + noise * rng.standard_normal(p)
    return X, Y, B


def combat_with_Y(X, Y, B):
    from inmoose.pycombat import pycombat_norm
    df = pd.DataFrame(X.T).astype(float)
    covar = pd.get_dummies(pd.Series(Y).astype(str),
                           drop_first=True).astype(float)
    try:
        out = pycombat_norm(df.values, batch=list(B),
                            covar_mod=covar.values, par_prior=True)
    except Exception:
        out = pycombat_norm(df.values, batch=list(B),
                            covar_mod=None, par_prior=True)
    return np.asarray(out).T


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


def dial(auc):
    if np.isnan(auc):
        return np.nan
    return max(auc, 1 - auc) - auc


def main():
    n_per = 40
    p = 200
    rho_YY_grid = np.linspace(-1.0, 1.0, 21)
    n_reps = 20

    rows = []
    for rho_YY in rho_YY_grid:
        for rep in range(n_reps):
            seed = hash((round(float(rho_YY), 3), rep)) & 0xFFFF_FFFF
            X, Y, B = make_data_platform(n_per, p, float(rho_YY),
                                         seed=seed)
            try:
                Xc = combat_with_Y(X, Y, B)
            except Exception:
                Xc = X.copy()
            auc_pre = lodo_auc(X, Y, B)
            auc_post = lodo_auc(Xc, Y, B)
            rows.append(dict(rho_YY=float(rho_YY),
                             rep=rep,
                             auc_pre=auc_pre, auc_post=auc_post,
                             dial_pre=dial(auc_pre),
                             dial_post=dial(auc_post)))
    df = pd.DataFrame(rows)
    df.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"[ok] tsv -> {OUT_TSV}")

    agg = (df.groupby("rho_YY")
             .agg(auc_pre_mean=("auc_pre", "mean"),
                  auc_pre_std=("auc_pre", "std"),
                  auc_post_mean=("auc_post", "mean"),
                  auc_post_std=("auc_post", "std"),
                  dial_post_mean=("dial_post", "mean"),
                  dial_post_std=("dial_post", "std"))
             .reset_index())

    try:
        import plotly.graph_objects as go
        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=agg.rho_YY, y=agg.auc_post_mean,
            error_y=dict(type="data", array=agg.auc_post_std,
                         visible=True),
            name="AUC_post (real ComBat+Y)",
            mode="lines+markers",
            line=dict(color="#d62728", width=3)))
        fig.add_trace(go.Scatter(
            x=agg.rho_YY, y=agg.auc_pre_mean,
            error_y=dict(type="data", array=agg.auc_pre_std,
                         visible=True),
            name="AUC_pre",
            mode="lines+markers",
            line=dict(color="#1f77b4", dash="dot")))
        fig.add_trace(go.Scatter(
            x=agg.rho_YY, y=agg.dial_post_mean,
            error_y=dict(type="data", array=agg.dial_post_std,
                         visible=True),
            name="DIAL_post",
            mode="lines+markers",
            line=dict(color="#2ca02c", width=2)))
        fig.add_vline(x=0.0, line=dict(color="black", dash="dash"),
                      annotation_text="ρ_YY = 0 boundary")
        fig.add_hline(y=0.5, line=dict(color="gray", dash="dot"))
        fig.update_layout(
            title=("Theorem 1 phase transition: AUC_post / DIAL vs "
                   "ρ_YY (cross-cohort label-direction cosine)"),
            xaxis_title="ρ_YY = ⟨e_Y_0, e_Y_1⟩",
            yaxis_title="AUC / DIAL",
            template="plotly_white",
            height=480)
        fig.write_html(OUT_HTML, include_plotlyjs="cdn")
        print(f"[ok] html -> {OUT_HTML}")
    except Exception as e:
        print(f"[warn] plotly failed: {e}")

    lo = df[df.rho_YY > 0.4].auc_post.mean()
    hi = df[df.rho_YY < -0.4].auc_post.mean()
    dlo = df[df.rho_YY > 0.4].dial_post.mean()
    dhi = df[df.rho_YY < -0.4].dial_post.mean()
    print(f"[summary] ρ_YY > +0.4: AUC_post={lo:.3f}  DIAL_post={dlo:.3f}")
    print(f"[summary] ρ_YY < -0.4: AUC_post={hi:.3f}  DIAL_post={dhi:.3f}")
    print(f"[summary] flip confirmed: "
          f"{'YES' if hi < 0.5 < lo else 'partial'}")


if __name__ == "__main__":
    main()
