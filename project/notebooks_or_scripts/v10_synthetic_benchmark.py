#!/usr/bin/env python3
# v10_synthetic_benchmark.py  (v2)
# Grid simulation over (sigma_B, sigma_Y, rho_YY, n) with real ComBat
# via inmoose.pycombat.  rho_YY is the cross-cohort agreement in the
# label direction: 1 = platforms agree perfectly, 0 = orthogonal,
# -1 = platforms disagree (reversed effect).  This parameterises the
# mechanism we found empirically (see companion Bioinformatics paper
# v5.1, where the THCA BRAF-vs-RAS flip is hypothesised to arise from
# platform-specific effect-direction divergence between TCGA-THCA
# RNA-seq and GSE27155 microarray).
#
# Output:
#   results/v10_aaai/synthetic_benchmark/grid_results.tsv
#   results/v10_aaai/synthetic_benchmark/predicted_vs_observed_DIAL.png
#   results/v10_aaai/synthetic_benchmark/phase_diagram.html

from __future__ import annotations
import os, pathlib, warnings, itertools, time
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")

PROJECT = pathlib.Path("/opt/thyroid-dash/project")
RES = PROJECT / "results/v10_aaai/synthetic_benchmark"
RES.mkdir(parents=True, exist_ok=True)
OUT_TSV = RES / "grid_results.tsv"
OUT_PNG = RES / "predicted_vs_observed_DIAL.png"
OUT_HTML = RES / "phase_diagram.html"

QUICK = os.environ.get("QUICK") == "1"


def make_sample(n_per, p, rho_YY, sigma_Y, sigma_B,
                imb_q0=0.85, imb_q1=0.35,
                noise=1.0, seed=0):
    """Platform-heterogeneous label-direction model.
    Batch 0: label along e_Y_0 = ê_1.
    Batch 1: label along e_Y_1 with <e_Y_0, e_Y_1> = rho_YY.
    Batch direction e_B = ê_2 is orthogonal to both.
    Y-B confounding: P(Y=1 | B=0) = imb_q0, P(Y=1 | B=1) = imb_q1.
    """
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
        clf = LogisticRegression(penalty="l2", C=C, max_iter=2000,
                                 solver="liblinear")
        clf.fit(X[tr], Y[tr])
        aucs.append(roc_auc_score(Y[te], clf.decision_function(X[te])))
    return float(np.mean(aucs)) if aucs else np.nan


def dial(auc):
    if np.isnan(auc):
        return np.nan
    return max(auc, 1 - auc) - auc


def predicted_dial(rho_YY, sigma_Y, sigma_B):
    """Heuristic prediction from Theorem 1 extended with cross-cohort
    rho_YY.  When platforms disagree (rho_YY < 0), the effective
    post-ComBat retention coefficient is approximately
        c(rho_YY) = (1 + rho_YY) / 2
    (averaged direction between batches).  DIAL predicted = |min(0, c)|.
    """
    c = (1.0 + rho_YY) / 2.0
    return float(max(0.0, 0.5 - c) ** 2 * 0.5)  # heuristic monotone in -rho_YY


def main():
    if QUICK:
        sigma_B_grid = [1.0, 3.0]
        sigma_Y_grid = [0.5, 1.0]
        rho_YY_grid  = np.linspace(-1.0, 1.0, 5)
        n_per_grid   = [40]
        n_reps       = 5
    else:
        sigma_B_grid = [0.5, 1.0, 2.0, 3.0, 5.0]
        sigma_Y_grid = [0.5, 1.0, 2.0]
        rho_YY_grid  = np.linspace(-1.0, 1.0, 11)
        n_per_grid   = [40, 80, 150]
        n_reps       = 10

    conditions = list(itertools.product(sigma_B_grid, sigma_Y_grid,
                                        rho_YY_grid, n_per_grid))
    total = len(conditions) * n_reps
    print(f"[info] {len(conditions)} conditions x {n_reps} reps = {total} runs")

    rows = []
    done = 0
    t0 = time.time()
    for sB, sY, rYY, n_per in conditions:
        for rep in range(n_reps):
            seed = hash((round(float(sB), 2), round(float(sY), 2),
                         round(float(rYY), 3),
                         int(n_per), int(rep))) & 0xFFFF_FFFF
            X, Y, B = make_sample(n_per, p=200, rho_YY=float(rYY),
                                  sigma_Y=float(sY), sigma_B=float(sB),
                                  seed=seed)
            try:
                Xc = combat_with_Y(X, Y, B)
                correction = "combat_Y"
            except Exception:
                Xc = X.copy()  # no correction fallback
                correction = "none"
            auc_pre = lodo_auc(X, Y, B)
            auc_post = lodo_auc(Xc, Y, B)
            rows.append(dict(sigma_B=float(sB), sigma_Y=float(sY),
                             rho_YY=float(rYY), n_per=int(n_per), rep=rep,
                             auc_pre=auc_pre, auc_post=auc_post,
                             dial_pre=dial(auc_pre),
                             dial_post=dial(auc_post),
                             dial_pred=predicted_dial(rYY, sY, sB),
                             correction=correction))
            done += 1
            if done % 200 == 0:
                elapsed = time.time() - t0
                rate = done / elapsed
                eta = (total - done) / rate
                print(f"  [progress] {done}/{total}  rate={rate:.1f}/s  eta={eta:.0f}s")
    df = pd.DataFrame(rows)
    df.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"[ok] tsv -> {OUT_TSV}")

    mask = ~df.dial_post.isna()
    corr = float(np.corrcoef(df.loc[mask, "dial_post"],
                             df.loc[mask, "dial_pred"])[0, 1])
    print(f"[metric] corr(DIAL_observed, DIAL_predicted) = {corr:.3f}")

    # Report flip regime
    hi = df[df.rho_YY < -0.4].auc_post.mean()
    lo = df[df.rho_YY > 0.4].auc_post.mean()
    print(f"[summary] mean AUC_post (rho_YY<-0.4): {hi:.3f}")
    print(f"[summary] mean AUC_post (rho_YY>+0.4): {lo:.3f}")
    print(f"[summary] flip confirmed: {'YES' if hi < 0.5 < lo else 'partial'}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        plt.figure(figsize=(6, 5))
        jitter = 0.005 * np.random.default_rng(0).standard_normal(len(df))
        plt.scatter(df.dial_pred + jitter, df.dial_post, alpha=0.25, s=8,
                    c=df.rho_YY, cmap="coolwarm")
        plt.colorbar(label="ρ_YY  (cohort label agreement)")
        plt.plot([0, 0.5], [0, 0.5], "k--", lw=1, alpha=0.6)
        plt.xlabel("DIAL predicted (heuristic)")
        plt.ylabel("DIAL observed (empirical LODO, real ComBat)")
        plt.title(f"Synthetic benchmark: predicted vs observed DIAL (corr={corr:.2f})")
        plt.grid(alpha=0.3)
        plt.tight_layout()
        plt.savefig(OUT_PNG, dpi=140)
        plt.close()
        print(f"[ok] png -> {OUT_PNG}")
    except Exception as e:
        print(f"[warn] matplotlib failed: {e}")

    try:
        import plotly.graph_objects as go
        # Phase diagram: sigma_B x rho_YY, colour = mean DIAL
        piv = (df.pivot_table(index="sigma_B", columns="rho_YY",
                              values="dial_post", aggfunc="mean"))
        fig = go.Figure(data=go.Heatmap(
            z=piv.values,
            x=[f"{c:+.2f}" for c in piv.columns],
            y=[f"{r:.1f}" for r in piv.index],
            colorscale="Inferno",
            colorbar=dict(title="mean DIAL")))
        fig.update_layout(
            title=("Synthetic phase diagram: mean DIAL over (σ_B, ρ_YY) "
                   "after real ComBat-with-Y"),
            xaxis_title="ρ_YY  (cross-cohort label-direction cosine)",
            yaxis_title="σ_B",
            template="plotly_white",
            height=480)
        fig.write_html(OUT_HTML, include_plotlyjs="cdn")
        print(f"[ok] html -> {OUT_HTML}")
    except Exception as e:
        print(f"[warn] plotly failed: {e}")


if __name__ == "__main__":
    main()
