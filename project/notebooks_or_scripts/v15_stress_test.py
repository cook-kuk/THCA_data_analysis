#!/usr/bin/env python3
"""v15 Synthetic stress test.

4 shift types × 5 magnitudes × 5 classifier families × 20 reps  (= 2000 cells).
Computes DIAL and target-AUC per cell; outputs phase diagram + ROC for
DIAL as binary detector of the 'cond_parallel' (subspace-aligned) regime.

Runtime budget: ~5 minutes CPU.
"""
from __future__ import annotations
import time, json
from pathlib import Path
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import roc_auc_score, roc_curve, auc as sk_auc

# reuse shared simulator from verify script
import sys
sys.path.insert(0, str(Path(__file__).parent))
from v15_theorem2_verify import (
    simulate as _simulate_for_decomposition,
    combat_mean_center,
)

PROJECT = Path("/opt/thyroid-dash/project")
RES = PROJECT / "results" / "v15_neurips" / "synthetic_stress"
FIGS = PROJECT / "reports" / "html" / "figs_interactive" / "v15"
CKPT = PROJECT / "results" / "v15_neurips" / "checkpoints"
RES.mkdir(parents=True, exist_ok=True)


def auc_flip(a: float) -> float:
    return max(a, 1.0 - a)


def make_classifier(name: str):
    if name == "LogReg":
        return LogisticRegression(max_iter=2000)
    if name == "SVM-RBF":
        return SVC(kernel="rbf", probability=True, gamma="scale")
    if name == "RF":
        return RandomForestClassifier(n_estimators=100, n_jobs=-1, random_state=0)
    if name == "MLP":
        return MLPClassifier(hidden_layer_sizes=(32,), max_iter=400, random_state=0)
    if name == "GB":
        return GradientBoostingClassifier(n_estimators=60, random_state=0)
    raise ValueError(name)


def run_cell(shift_type: str, mag: float, clf_name: str, rng: np.random.Generator, p=20, n=150):
    # Build data the same way the verify simulator does.
    wY = np.zeros(p); wY[0] = 1.0
    wB = 0.9 * wY + np.sqrt(1 - 0.81) * np.eye(p)[1]
    m0, m1 = -wY, +wY

    skew = 0.0
    add0_t = add1_t = np.zeros(p)
    if shift_type == "covariate":
        add0_t = add1_t = mag * np.eye(p)[2]
    elif shift_type == "label":
        skew = 0.4 * min(mag, 1.0)
    elif shift_type == "concept":
        # Target label definition flips (mag=1.0 → full inversion).
        add0_t = +2 * mag * wY; add1_t = -2 * mag * wY
    elif shift_type == "subspace_aligned":
        # Class-conditional residual along span(wY) — our target regime.
        # mag=1 → 1.5x source class separation (full flip).
        add0_t = +1.5 * mag * wY; add1_t = -1.5 * mag * wY
        skew = 0.15

    ns0 = int(n * (1 + skew)); ns1 = int(n * (1 - skew))
    nt0 = int(n * (1 - skew)); nt1 = int(n * (1 + skew))
    Xs_0 = rng.normal(size=(ns0, p)) + m0
    Xs_1 = rng.normal(size=(ns1, p)) + m1
    Xt_0 = rng.normal(size=(nt0, p)) + m0 + add0_t
    Xt_1 = rng.normal(size=(nt1, p)) + m1 + add1_t
    X = np.vstack([Xs_0, Xs_1, Xt_0, Xt_1])
    Y = np.concatenate([np.zeros(ns0), np.ones(ns1), np.zeros(nt0), np.ones(nt1)]).astype(int)
    B = np.concatenate([np.zeros(ns0 + ns1), np.ones(nt0 + nt1)]).astype(int)

    # Naive ComBat (no covariate) — flip-prone.
    Xc = combat_mean_center(X, B, Y=None, preserve_subtype=False)

    src = B == 0; tgt = B == 1
    try:
        clf = make_classifier(clf_name).fit(Xc[src], Y[src])
        proba = clf.predict_proba(Xc[tgt])[:, 1]
        auc_tgt = float(roc_auc_score(Y[tgt], proba))
        # source CV "looks good" AUC, for comparison
        clf_src = make_classifier(clf_name).fit(Xc[src], Y[src])
        proba_src = clf_src.predict_proba(Xc[src])[:, 1]
        auc_src = float(roc_auc_score(Y[src], proba_src))
    except Exception as e:
        return dict(shift_type=shift_type, mag=mag, clf=clf_name,
                    auc_src=float("nan"), auc_tgt=float("nan"),
                    dial=float("nan"), error=str(e)[:80])

    dial = (auc_flip(auc_tgt) - 0.5) if auc_tgt < 0.5 else 0.0
    return dict(shift_type=shift_type, mag=mag, clf=clf_name,
                auc_src=auc_src, auc_tgt=auc_tgt, dial=dial)


def main():
    t0 = time.time()
    shift_types = ["covariate", "label", "concept", "subspace_aligned"]
    mags = [0.0, 0.25, 0.5, 0.75, 1.0]
    clfs = ["LogReg", "SVM-RBF", "RF", "MLP", "GB"]
    n_rep = 20

    rows = []
    total = len(shift_types) * len(mags) * len(clfs) * n_rep
    count = 0
    for st in shift_types:
        for mg in mags:
            for cn in clfs:
                for rep in range(n_rep):
                    count += 1
                    rng = np.random.default_rng(hash((st, mg, cn, rep)) % (2**31))
                    r = run_cell(st, mg, cn, rng)
                    r["rep"] = rep
                    rows.append(r)
                    if count % 200 == 0:
                        print(f"  [{count}/{total}] {st} mg={mg} clf={cn}", flush=True)

    df = pd.DataFrame(rows)
    out_tsv = RES / "stress_test_grid.tsv"
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"Wrote {out_tsv} rows={len(df)}")

    # Summary table
    summ = df.groupby(["shift_type", "mag"]).agg(
        dial_mean=("dial", "mean"), dial_std=("dial", "std"),
        auc_tgt_mean=("auc_tgt", "mean"),
        auc_src_mean=("auc_src", "mean"),
    ).round(3).reset_index()
    summ.to_csv(RES / "stress_test_summary.tsv", sep="\t", index=False)
    print(summ.pivot(index="mag", columns="shift_type", values="dial_mean"))

    # ROC: DIAL as binary detector of subspace_aligned (positive) vs others (negative)
    # Use non-trivial magnitude (>0) rows to avoid baseline confusion.
    non_zero = df[df.mag > 0.0]
    y_true = (non_zero.shift_type == "subspace_aligned").astype(int).values
    scores = non_zero.dial.values
    fpr, tpr, thr = roc_curve(y_true, scores)
    roc_auc = sk_auc(fpr, tpr)

    # DIAL also expected to flag 'concept' (hard to distinguish from ours). Report both.
    y_true_both = non_zero.shift_type.isin(["subspace_aligned", "concept"]).astype(int).values
    fpr2, tpr2, _ = roc_curve(y_true_both, scores)
    roc_auc_both = sk_auc(fpr2, tpr2)

    print(f"ROC: DIAL vs subspace_aligned only AUC = {roc_auc:.3f}")
    print(f"ROC: DIAL vs subspace_aligned+concept AUC = {roc_auc_both:.3f}")

    # Phase diagram (shift_type × magnitude)
    ptab = df.groupby(["shift_type", "mag"])["dial"].mean().unstack("mag")
    fig = make_subplots(rows=1, cols=2,
                        subplot_titles=("Phase diagram: DIAL by shift_type × magnitude",
                                        "ROC — DIAL as shift-type detector"),
                        column_widths=[0.55, 0.45], horizontal_spacing=0.12)
    fig.add_trace(go.Heatmap(
        z=ptab.values, x=[f"{m:g}" for m in ptab.columns], y=ptab.index,
        colorscale=[[0, "#0b1220"], [0.15, "#1e3a8a"], [0.5, "#f59e0b"], [1, "#f43f5e"]],
        colorbar=dict(title="DIAL"), zmin=0, zmax=0.5,
    ), row=1, col=1)

    fig.add_trace(go.Scatter(x=fpr, y=tpr, mode="lines",
                             name=f"subspace_aligned (AUC={roc_auc:.2f})",
                             line=dict(color="#f59e0b", width=3)), row=1, col=2)
    fig.add_trace(go.Scatter(x=fpr2, y=tpr2, mode="lines",
                             name=f"+concept (AUC={roc_auc_both:.2f})",
                             line=dict(color="#34d399", width=3, dash="dash")), row=1, col=2)
    fig.add_trace(go.Scatter(x=[0, 1], y=[0, 1], mode="lines",
                             line=dict(color="rgba(255,255,255,0.3)", dash="dot"),
                             showlegend=False), row=1, col=2)
    fig.update_xaxes(title_text="magnitude", row=1, col=1, gridcolor="rgba(255,255,255,0.1)")
    fig.update_yaxes(title_text="shift_type", row=1, col=1, gridcolor="rgba(255,255,255,0.1)")
    fig.update_xaxes(title_text="FPR", row=1, col=2, gridcolor="rgba(255,255,255,0.1)")
    fig.update_yaxes(title_text="TPR", row=1, col=2, gridcolor="rgba(255,255,255,0.1)")

    fig.update_layout(
        title=dict(text="Synthetic stress test: DIAL discriminates subspace-aligned shift",
                   x=0.01, font=dict(color="#e2e8f0", size=18)),
        height=520, paper_bgcolor="#0b1220", plot_bgcolor="#0b1220",
        font=dict(color="#e2e8f0"),
        legend=dict(x=0.78, y=0.22, bgcolor="rgba(11,18,32,0.6)"),
        margin=dict(l=60, r=30, t=90, b=50),
    )
    out_html = FIGS / "stress_test_phase_diagram.html"
    fig.write_html(out_html, include_plotlyjs="cdn", config=dict(displayModeBar=False))
    print(f"Wrote {out_html}")

    (CKPT / "task2_stress.json").write_text(json.dumps(dict(
        n_rows=len(df), roc_auc_subspace=roc_auc, roc_auc_any=roc_auc_both,
        wall=time.time() - t0,
    ), indent=2))


if __name__ == "__main__":
    main()
