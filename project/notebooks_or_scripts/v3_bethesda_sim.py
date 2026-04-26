"""v3 STEP 6 — Bethesda prevalence sweep simulation.

Synthetic cohort with controlled malignancy prevalence; apply best v3 classifier
and report operating curves + decision curve + surgery-reduction bar chart.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, "/opt/thyroid-dash/project/notebooks_or_scripts")
from v3_common import (DATASET_FILES, FIGS, MAPK_OUTPUT, RESULTS_TABLES,
                        get_tierA67_clean, load_expr, small_n_warning,
                        write_tsv, _expr_sample_id_column)


def load_pools():
    genes = [g for g in get_tierA67_clean() if g not in set(MAPK_OUTPUT)]
    expr = load_expr(DATASET_FILES["TCGA-THCA"])
    sm = pd.read_csv("/opt/thyroid-dash/project/metadata/sample_master.tsv",
                     sep="\t")
    sm_t = sm[sm.dataset == "TCGA-THCA"].copy()

    # malignant pool: tumor w/ BRAF_like OR dedifferentiated
    malig = sm_t[(sm_t.normal_vs_tumor == "tumor") &
                  sm_t.molecular_subtype.isin(["BRAF_like", "dedifferentiated"])]
    # benign pool: tumor RAS_like + normals
    benign = sm_t[((sm_t.normal_vs_tumor == "tumor") &
                    (sm_t.molecular_subtype == "RAS_like")) |
                   (sm_t.normal_vs_tumor == "normal")]

    def make_matrix(frame):
        f = frame.copy()
        f["_expr_id"] = _expr_sample_id_column("TCGA-THCA", f)
        samps = [s for s in f["_expr_id"] if s in expr.columns]
        return expr.reindex(genes).fillna(0.0).T.loc[samps].to_numpy(float)

    M = make_matrix(malig)
    B = make_matrix(benign)
    return genes, M, B


def train_base_classifier(genes, M, B):
    X = np.vstack([M, B])
    y = np.concatenate([np.ones(len(M), int), np.zeros(len(B), int)])
    sc = StandardScaler().fit(X)
    clf = LogisticRegression(penalty="l2", max_iter=1500,
                             class_weight="balanced",
                             solver="liblinear")
    clf.fit(sc.transform(X), y)
    return sc, clf


def simulate_cohort(M, B, prev, n=10000, sigma=0.3, seed=42):
    rng = np.random.default_rng(seed)
    n_m = int(round(prev * n))
    n_b = n - n_m
    # sample with replacement
    idx_m = rng.integers(0, len(M), n_m)
    idx_b = rng.integers(0, len(B), n_b)
    Xm = M[idx_m] + rng.normal(0, sigma, (n_m, M.shape[1]))
    Xb = B[idx_b] + rng.normal(0, sigma, (n_b, B.shape[1]))
    X = np.vstack([Xm, Xb])
    y = np.concatenate([np.ones(n_m, int), np.zeros(n_b, int)])
    perm = rng.permutation(n)
    return X[perm], y[perm]


def sweep_thresholds(y_true, p, prev):
    out = []
    for t in np.arange(0.05, 0.951, 0.01):
        yp = (p >= t).astype(int)
        tp = int(((y_true == 1) & (yp == 1)).sum())
        fp = int(((y_true == 0) & (yp == 1)).sum())
        tn = int(((y_true == 0) & (yp == 0)).sum())
        fn = int(((y_true == 1) & (yp == 0)).sum())
        n = len(y_true)
        se = tp / (tp + fn) if (tp + fn) > 0 else np.nan
        sp = tn / (tn + fp) if (tn + fp) > 0 else np.nan
        npv = tn / (tn + fn) if (tn + fn) > 0 else np.nan
        ppv = tp / (tp + fp) if (tp + fp) > 0 else np.nan
        unnec = fp / n
        missed = fn / n
        nb = tp / n - fp / n * (t / max(1 - t, 1e-6))
        nb_all = prev - (1 - prev) * (t / max(1 - t, 1e-6))
        out.append(dict(threshold=float(t), prev=float(prev),
                         sens=se, spec=sp, npv=npv, ppv=ppv,
                         unnecessary_surgery=unnec,
                         missed_cancer=missed,
                         nb_model=nb, nb_treat_all=nb_all, nb_treat_none=0.0,
                         tp=tp, fp=fp, tn=tn, fn=fn))
    return out


def run():
    genes, M, B = load_pools()
    sc, clf = train_base_classifier(genes, M, B)

    # Isotonic recalibration on 20% holdout of the combined TCGA pool
    rng = np.random.default_rng(0)
    X_all = np.vstack([M, B])
    y_all = np.concatenate([np.ones(len(M), int), np.zeros(len(B), int)])
    perm = rng.permutation(len(y_all))
    cal_idx = perm[:int(0.2 * len(y_all))]
    p_cal = clf.predict_proba(sc.transform(X_all[cal_idx]))[:, 1]
    iso = IsotonicRegression(out_of_bounds="clip")
    iso.fit(p_cal, y_all[cal_idx])

    prevs = [0.10, 0.15, 0.20, 0.25, 0.30]
    all_rows = []
    surgery_rows = []
    for prev in prevs:
        Xs, ys = simulate_cohort(M, B, prev, n=10000, sigma=0.3,
                                  seed=42 + int(prev * 100))
        p_raw = clf.predict_proba(sc.transform(Xs))[:, 1]
        p_iso = iso.transform(p_raw)
        rows = sweep_thresholds(ys, p_iso, prev)
        all_rows.extend(rows)
        # Pick threshold at Se>=0.95 -> report unnecessary surgery
        df_tmp = pd.DataFrame(rows)
        ok = df_tmp[df_tmp["sens"] >= 0.95]
        if len(ok) == 0:
            continue
        chosen = ok.sort_values("threshold", ascending=False).iloc[0]
        baseline_surgery = 1.0  # treat all: 100%
        model_surgery = chosen["tp"] + chosen["fp"]
        total = chosen[["tp", "fp", "tn", "fn"]].sum()
        model_surgery_rate = model_surgery / total
        reduction = baseline_surgery - model_surgery_rate
        surgery_rows.append(dict(prev=prev, chosen_threshold=float(chosen["threshold"]),
                                  sens=float(chosen["sens"]),
                                  spec=float(chosen["spec"]),
                                  npv=float(chosen["npv"]),
                                  ppv=float(chosen["ppv"]),
                                  model_surgery_rate=float(model_surgery_rate),
                                  reduction_vs_treat_all=float(reduction)))

    df = pd.DataFrame(all_rows)
    write_tsv(df, RESULTS_TABLES / "v3_bethesda_sim.tsv")

    surg_df = pd.DataFrame(surgery_rows)
    write_tsv(surg_df, RESULTS_TABLES / "v3_bethesda_surgery_reduction.tsv")

    # Figure 1: operating curves (Se, Sp, NPV vs threshold, one trace per prev)
    fig = make_subplots(rows=2, cols=2,
                         subplot_titles=("Sensitivity", "Specificity",
                                         "NPV", "Unnecessary surgery rate"))
    palette = ["#5eead4", "#38bdf8", "#f472b6", "#fb7185", "#a78bfa"]
    for i, prev in enumerate(prevs):
        sub = df[df.prev == prev].sort_values("threshold")
        color = palette[i]
        fig.add_scatter(x=sub["threshold"], y=sub["sens"], mode="lines",
                         name=f"prev={prev}", line=dict(color=color),
                         legendgroup=f"p{prev}", row=1, col=1)
        fig.add_scatter(x=sub["threshold"], y=sub["spec"], mode="lines",
                         line=dict(color=color),
                         legendgroup=f"p{prev}", showlegend=False, row=1, col=2)
        fig.add_scatter(x=sub["threshold"], y=sub["npv"], mode="lines",
                         line=dict(color=color),
                         legendgroup=f"p{prev}", showlegend=False, row=2, col=1)
        fig.add_scatter(x=sub["threshold"], y=sub["unnecessary_surgery"],
                         mode="lines", line=dict(color=color),
                         legendgroup=f"p{prev}", showlegend=False, row=2, col=2)
    fig.update_layout(title="Bethesda simulation: operating curves across prevalence",
                      template="plotly_dark", height=720)
    fig.write_html(FIGS / "v3_bethesda_operating.html",
                    include_plotlyjs="cdn")

    # Figure 2: decision curve across prevalence (at one representative prev=0.20)
    fig_dc = go.Figure()
    for prev in prevs:
        sub = df[df.prev == prev].sort_values("threshold")
        fig_dc.add_scatter(x=sub["threshold"], y=sub["nb_model"],
                            mode="lines", name=f"model prev={prev}")
    # reference treat-all at prev=0.20
    ref = df[df.prev == 0.20].sort_values("threshold")
    fig_dc.add_scatter(x=ref["threshold"], y=ref["nb_treat_all"],
                        mode="lines", line=dict(color="#94a3b8", dash="dot"),
                        name="treat-all (p=0.20)")
    fig_dc.add_hline(y=0.0, line=dict(color="#94a3b8", dash="dash"),
                      annotation_text="treat-none")
    fig_dc.update_layout(title="Decision curve across prevalences",
                         xaxis_title="Threshold", yaxis_title="Net benefit",
                         template="plotly_dark", height=500)
    fig_dc.write_html(FIGS / "v3_bethesda_decision_curve.html",
                      include_plotlyjs="cdn")

    # Figure 3: surgery reduction bar at Se>=0.95
    fig_sr = go.Figure()
    fig_sr.add_bar(x=[f"{p:.0%}" for p in surg_df["prev"]],
                    y=surg_df["reduction_vs_treat_all"] * 100.0,
                    marker_color="#5eead4",
                    text=[f"{r * 100:.0f}%" for r in surg_df["reduction_vs_treat_all"]],
                    textposition="outside")
    fig_sr.update_layout(
        title="Surgery-reduction at Se≥0.95 vs treat-all (simulated)",
        xaxis_title="Malignancy prevalence",
        yaxis_title="Reduction in surgery rate (%)",
        template="plotly_dark", height=420)
    fig_sr.write_html(FIGS / "v3_bethesda_surgery_reduction.html",
                      include_plotlyjs="cdn")

    return dict(prevs=prevs, best_rows=len(surg_df),
                 out=str(RESULTS_TABLES / "v3_bethesda_sim.tsv"))


def main():
    res = run()
    print(json.dumps(res, indent=2, default=str))


if __name__ == "__main__":
    main()
