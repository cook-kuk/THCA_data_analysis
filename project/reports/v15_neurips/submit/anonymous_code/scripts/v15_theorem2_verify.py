#!/usr/bin/env python3
"""v15 Theorem 2 numerical verification.

Simulates shift components (covariate / label / conditional-parallel /
conditional-perp) and regresses observed DIAL against the Theorem 2
decomposition. Output: theorem2_numerical.tsv + interactive plot.

Honest scope: Gaussian class-conditionals, shared covariance, linear
ComBat-style mean-centring correction. Runs ~30s on CPU.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score

PROJECT = Path("/opt/thyroid-dash/project")
RES = PROJECT / "results" / "v15_neurips" / "theory_validation"
FIGS = PROJECT / "reports" / "html" / "figs_interactive" / "v15"
CKPT = PROJECT / "results" / "v15_neurips" / "checkpoints"
RES.mkdir(parents=True, exist_ok=True)
FIGS.mkdir(parents=True, exist_ok=True)
CKPT.mkdir(parents=True, exist_ok=True)

RNG = np.random.default_rng(20260425)


def auc_flip(a: float) -> float:
    return max(a, 1.0 - a)


def compute_dial(X: np.ndarray, Y: np.ndarray, B: np.ndarray | None = None) -> tuple[float, float]:
    """If B given: train on B=0, test on B=1 (source→target). Else 5-fold CV."""
    if B is not None:
        src = B == 0; tgt = B == 1
        if len(np.unique(Y[src])) < 2 or len(np.unique(Y[tgt])) < 2:
            return float("nan"), float("nan")
        clf = LogisticRegression(max_iter=2000).fit(X[src], Y[src])
        proba = clf.predict_proba(X[tgt])[:, 1]
        auc = float(roc_auc_score(Y[tgt], proba))
    else:
        skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=0)
        aucs = []
        for tr, te in skf.split(X, Y):
            clf = LogisticRegression(max_iter=2000).fit(X[tr], Y[tr])
            proba = clf.predict_proba(X[te])[:, 1]
            aucs.append(roc_auc_score(Y[te], proba))
        auc = float(np.mean(aucs))
    dial = (auc_flip(auc) - 0.5) if auc < 0.5 else 0.0
    return dial, auc


def combat_mean_center(X: np.ndarray, B: np.ndarray, Y: np.ndarray | None = None,
                       preserve_subtype: bool = False) -> np.ndarray:
    """Linear ComBat surrogate.

    * Without covariates (preserve_subtype=False, Y=None): subtracts the batch
      mean and adds back the global mean.
    * With subtype preservation: for each (batch, class), the *deviation* of
      the within-batch class mean from the global class mean is removed; the
      global class means are preserved.
    """
    X = X.copy().astype(float)
    gmean = X.mean(axis=0, keepdims=True)

    if not preserve_subtype or Y is None:
        for b in np.unique(B):
            idx = B == b
            X[idx] = X[idx] - X[idx].mean(axis=0, keepdims=True) + gmean
        return X

    for y in np.unique(Y):
        my = idx_y = Y == y
        gmy = X[my].mean(axis=0, keepdims=True)
        for b in np.unique(B):
            m = (B == b) & my
            if m.sum() > 1:
                X[m] = X[m] - X[m].mean(axis=0, keepdims=True) + gmy
    return X


def simulate(
    shift_type: str,
    mag: float,
    p: int = 20,
    n_per_class: int = 150,
    align: float = 0.9,
    rng: np.random.Generator | None = None,
) -> dict:
    """Generate a source+target pair and evaluate DIAL.

    shift_type in {"none","covariate","label","cond_parallel","cond_perp","subspace_aligned"}.
    mag scales the shift.
    align controls the cos(w_Y, w_B) used for the subspace-aligned regime.
    """
    rng = rng or np.random.default_rng()

    # Label direction w_Y = e_0
    wY = np.zeros(p); wY[0] = 1.0
    # Batch direction w_B: rotate wY by angle theta so that cos = align
    wB = align * wY + np.sqrt(max(1 - align**2, 0.0)) * np.eye(p)[1]

    # Class means along wY
    m0, m1 = -wY, +wY

    # Class-ratio skew per batch — this is what makes naive ComBat dangerous.
    # skew in [-0.4, 0.4]: source batch class-1 fraction = 0.5-skew, target = 0.5+skew.
    skew = 0.0
    add0_s = np.zeros(p); add1_s = np.zeros(p)  # source added shift (rare)
    add0_t = np.zeros(p); add1_t = np.zeros(p)

    if shift_type == "none":
        pass
    elif shift_type == "covariate":
        v = mag * np.eye(p)[2]   # orthogonal to wY
        add0_t = add1_t = v
    elif shift_type == "label":
        skew = 0.4 * min(mag, 1.0)   # prior shift magnitude (bounded)
    elif shift_type == "cond_parallel":
        # Class-conditional shift parallel to wY — direct label flip in target.
        add0_t = +mag * wY
        add1_t = -mag * wY
    elif shift_type == "cond_perp":
        ortho = np.eye(p)[3]
        add0_t = -mag * ortho
        add1_t = +mag * ortho
    elif shift_type == "subspace_aligned":
        # The *canonical* THCA-style flip: class imbalance per batch + batch shift
        # along a direction that aligns with wY (align coefficient).
        skew = 0.3 * min(mag, 1.0)
        add0_t = add1_t = mag * wB
    else:
        raise ValueError(shift_type)

    # Sample counts with skew
    ns0 = int(round(n_per_class * (1 + skew)))
    ns1 = int(round(n_per_class * (1 - skew)))
    nt0 = int(round(n_per_class * (1 - skew)))
    nt1 = int(round(n_per_class * (1 + skew)))

    Xs_0 = rng.normal(size=(ns0, p)) + m0 + add0_s
    Xs_1 = rng.normal(size=(ns1, p)) + m1 + add1_s
    Xt_0 = rng.normal(size=(nt0, p)) + m0 + add0_t
    Xt_1 = rng.normal(size=(nt1, p)) + m1 + add1_t

    X_raw = np.vstack([Xs_0, Xs_1, Xt_0, Xt_1])
    Y = np.concatenate([np.zeros(ns0), np.ones(ns1),
                        np.zeros(nt0), np.ones(nt1)]).astype(int)
    B = np.concatenate([np.zeros(ns0 + ns1),
                        np.ones(nt0 + nt1)]).astype(int)
    n0 = nt0; n1 = nt1; prior_shift = skew

    # Naive ComBat (no covariate): the flip-prone variant per Theorem 1.
    X_post = combat_mean_center(X_raw, B, Y=None, preserve_subtype=False)

    # Evaluate source→target transfer AUC.
    dial, auc = compute_dial(X_post, Y, B=B)

    # Theoretical KL components (population-level, closed form for shared-Σ Gaussians)
    Xs_all = np.vstack([Xs_0, Xs_1]); Xt_all = np.vstack([Xt_0, Xt_1])
    delta_cov = 0.5 * float(np.linalg.norm(Xt_all.mean(axis=0) - Xs_all.mean(axis=0)) ** 2)
    delta_lab = 0.0 if prior_shift == 0 else 2 * prior_shift**2  # rough 2nd-order
    diff = (add1_t - add0_t) if shift_type not in ("label", "covariate") else np.zeros(p)
    delta_parallel = 0.5 * (wY @ diff) ** 2
    delta_perp = 0.5 * np.linalg.norm(diff - (wY @ diff) * wY) ** 2

    return dict(
        shift_type=shift_type,
        mag=mag,
        align=align,
        dial=dial,
        auc=auc,
        delta_cov=delta_cov,
        delta_lab=delta_lab,
        delta_parallel=delta_parallel,
        delta_perp=delta_perp,
    )


def main():
    t0 = time.time()
    rows = []
    shift_types = ["none", "covariate", "label",
                   "cond_parallel", "cond_perp", "subspace_aligned"]
    mags = [0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0]
    n_rep = 15

    for st in shift_types:
        for mag in mags:
            for rep in range(n_rep):
                rng = np.random.default_rng(20260425 + hash((st, mag, rep)) % (2**31))
                r = simulate(st, mag, rng=rng)
                r["rep"] = rep
                rows.append(r)

    df = pd.DataFrame(rows)
    out_tsv = RES / "theorem2_numerical.tsv"
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"Wrote {out_tsv} rows={len(df)}")

    # Regression: DIAL ~ delta_parallel (Theorem 2 prediction)
    from sklearn.linear_model import LinearRegression
    Xfeat = df[["delta_cov", "delta_lab", "delta_parallel", "delta_perp"]].values
    y = df["dial"].values
    reg = LinearRegression().fit(Xfeat, y)
    r2 = reg.score(Xfeat, y)
    coef = dict(zip(["delta_cov", "delta_lab", "delta_parallel", "delta_perp"], reg.coef_.tolist()))
    print(f"Linear regression DIAL ~ KL components: R^2={r2:.3f}, coef={coef}")

    # Save checkpoint json
    (CKPT / "task1_theorem2.json").write_text(json.dumps(dict(
        r2=r2, coef=coef, n_rows=len(df), wall=time.time() - t0
    ), indent=2))

    # Interactive figure
    fig = make_subplots(rows=2, cols=3, subplot_titles=shift_types,
                        horizontal_spacing=0.09, vertical_spacing=0.14)
    colour = dict(zip(shift_types,
                      ["#7dd3fc", "#60a5fa", "#a78bfa",
                       "#f59e0b", "#34d399", "#f43f5e"]))
    for idx, st in enumerate(shift_types):
        r, c = idx // 3 + 1, idx % 3 + 1
        d = df[df.shift_type == st]
        g = d.groupby("mag")
        mu, se = g["dial"].mean(), g["dial"].std() / np.sqrt(n_rep)
        fig.add_trace(go.Scatter(
            x=mu.index, y=mu.values,
            error_y=dict(type="data", array=se.values, color="#94a3b8"),
            mode="lines+markers", name=st,
            line=dict(color=colour[st], width=2),
            marker=dict(size=8),
            showlegend=False,
        ), row=r, col=c)
        fig.update_xaxes(title_text="shift magnitude", row=r, col=c,
                         gridcolor="rgba(255,255,255,0.08)")
        fig.update_yaxes(title_text="DIAL", range=[-0.02, 0.5],
                         row=r, col=c, gridcolor="rgba(255,255,255,0.08)")

    fig.update_layout(
        title=dict(
            text=f"Theorem 2 numerical validation — DIAL ~ KL components  (R²={r2:.3f})",
            x=0.01, xanchor="left", font=dict(color="#e2e8f0", size=18),
        ),
        height=600, paper_bgcolor="#0b1220", plot_bgcolor="#0b1220",
        font=dict(color="#e2e8f0"),
        margin=dict(l=60, r=30, t=90, b=50),
    )
    out_html = FIGS / "theorem2_decomposition.html"
    fig.write_html(out_html, include_plotlyjs="cdn",
                   config=dict(displayModeBar=False))
    print(f"Wrote {out_html}")

    # Summary print
    print("\n=== DIAL average by shift_type at max magnitude ===")
    print(df[df.mag == 2.0].groupby("shift_type")["dial"].agg(["mean", "std"]).round(3))


if __name__ == "__main__":
    main()
