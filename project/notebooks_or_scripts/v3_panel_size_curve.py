"""v3 STEP 3 — Panel size curve.

For k in [4,8,16,24,32,40,50,60] and 3 strategies, evaluate:
  - TCGA 5-fold CV AUC + bootstrap CI
  - External AUC on GSE27155 and GSE126698 (isotonic recalibrated)
  - NPV at Se>=0.95

Feature selection happens INSIDE each CV fold.
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
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, "/opt/thyroid-dash/project/notebooks_or_scripts")
from v3_common import (DATASET_FILES, FIGS, MAPK_OUTPUT, RESULTS_TABLES,
                        bootstrap_ci, get_tds16, get_tierA67_clean, load_expr,
                        small_n_warning, write_tsv, _expr_sample_id_column)

KS = [4, 8, 16, 24, 32, 40, 50, 60]
STRATEGIES = ["univariate_d", "tds16_union", "qubo_neal"]


def cohens_d_vec(X, y):
    pos = X[y == 1]
    neg = X[y == 0]
    mu1 = pos.mean(axis=0)
    mu0 = neg.mean(axis=0)
    n1, n0 = len(pos), len(neg)
    s1 = pos.var(axis=0, ddof=1) if n1 > 1 else np.zeros(X.shape[1])
    s0 = neg.var(axis=0, ddof=1) if n0 > 1 else np.zeros(X.shape[1])
    denom = np.sqrt(((n1 - 1) * s1 + (n0 - 1) * s0) / max(n0 + n1 - 2, 1))
    denom[denom == 0] = 1.0
    return (mu1 - mu0) / denom


def select_univariate(X, y, gene_names, k):
    d = np.abs(cohens_d_vec(X, y))
    idx = np.argsort(-d)[:k]
    return idx


def select_tds16_union(X, y, gene_names, k):
    tds = [g for g in get_tds16() if g in gene_names]
    tds_idx = [gene_names.index(g) for g in tds]
    remaining = [i for i in range(len(gene_names)) if i not in set(tds_idx)]
    d = np.abs(cohens_d_vec(X, y))
    rem_sorted = sorted(remaining, key=lambda i: -d[i])
    need = max(k - len(tds_idx), 0)
    chosen = list(tds_idx) + rem_sorted[:need]
    return np.array(chosen[:k])


def _build_mrmr_qubo(X, y, k_target, alpha=0.8):
    """Inlined copy of ml_quantum_full.build_mrmr_qubo to avoid triggering
    ml_quantum_full's global init/logging on import."""
    n = X.shape[1]
    Xs = (X - X.mean(axis=0)) / (X.std(axis=0) + 1e-9)
    ys = (y - y.mean()) / (y.std() + 1e-9)
    rel = np.abs(Xs.T @ ys) / len(ys)
    R = np.abs(np.corrcoef(Xs.T))
    np.fill_diagonal(R, 0.0)
    lam = 3.0 * max(rel.max(), 1e-3)
    linear = {}
    quad = {}
    for i in range(n):
        linear[i] = -alpha * float(rel[i]) + lam * (1.0 - 2.0 * k_target)
    for i in range(n):
        for j in range(i + 1, n):
            quad[(i, j)] = (1.0 - alpha) * float(R[i, j]) + 2.0 * lam
    offset = lam * k_target * k_target
    return linear, quad, offset


def select_qubo_neal(X, y, gene_names, k, seed=42):
    """QUBO mRMR feature selection, solved by neal simulated-annealing.
    Reduced reads/sweeps to keep 120-fold wall-time manageable."""
    try:
        import dimod
        import neal
        linear, quad, offset = _build_mrmr_qubo(X, y, k_target=k, alpha=0.8)
        bqm = dimod.BinaryQuadraticModel.from_qubo(
            {**{(i, i): v for i, v in linear.items()}, **quad}, offset=offset)
        sampler = neal.SimulatedAnnealingSampler()
        ss = sampler.sample(bqm, num_reads=20, num_sweeps=80, seed=seed)
        x_best = ss.first.sample
        xbin = [int(x_best[i]) for i in range(len(linear))]
        if xbin is None:
            return select_univariate(X, y, gene_names, k)
        chosen = np.where(np.array(xbin) > 0)[0]
        if len(chosen) == 0:
            return select_univariate(X, y, gene_names, k)
        if len(chosen) > k:
            # tiebreak by relevance
            d = np.abs(cohens_d_vec(X, y))
            chosen = chosen[np.argsort(-d[chosen])[:k]]
        elif len(chosen) < k:
            # top up
            rest = [i for i in range(X.shape[1]) if i not in set(chosen)]
            d = np.abs(cohens_d_vec(X, y))
            rest_sorted = sorted(rest, key=lambda i: -d[i])
            need = k - len(chosen)
            chosen = np.concatenate([chosen, np.array(rest_sorted[:need])])
        return chosen.astype(int)
    except Exception as e:
        print(f"[qubo_neal fallback] {e}")
        return select_univariate(X, y, gene_names, k)


SELECTORS = {
    "univariate_d": select_univariate,
    "tds16_union": select_tds16_union,
    "qubo_neal": select_qubo_neal,
}


def npv_at_sens(y_true, y_score, sens_target=0.95):
    if len(np.unique(y_true)) < 2:
        return np.nan, np.nan
    fpr, tpr, thr = roc_curve(y_true, y_score)
    ok = tpr >= sens_target
    if not ok.any():
        return np.nan, np.nan
    idx = int(np.argmax(ok))
    t = thr[idx]
    yp = (y_score >= t).astype(int)
    tn = int(((y_true == 0) & (yp == 0)).sum())
    fn = int(((y_true == 1) & (yp == 0)).sum())
    denom = tn + fn
    return (tn / denom if denom > 0 else np.nan), float(t)


def load_all():
    genes = [g for g in get_tierA67_clean() if g not in set(MAPK_OUTPUT)]

    def load(ds):
        expr = load_expr(DATASET_FILES[ds])
        sm = pd.read_csv("/opt/thyroid-dash/project/metadata/sample_master.tsv",
                         sep="\t")
        sm = sm[(sm.dataset == ds) & (sm.normal_vs_tumor == "tumor") &
                 sm.molecular_subtype.isin(["BRAF_like", "RAS_like"])].copy()
        sm["_expr_id"] = _expr_sample_id_column(ds, sm)
        sm = sm[sm["_expr_id"].isin(expr.columns)]
        if sm.empty:
            return None, None
        samples = sm["_expr_id"].tolist()
        y = sm.set_index("_expr_id").loc[samples, "molecular_subtype"].map(
            {"BRAF_like": 1, "RAS_like": 0}).to_numpy(int)
        X = expr.reindex(genes).fillna(0.0).T.loc[samples].to_numpy(float)
        return X, y

    X_tcga, y_tcga = load("TCGA-THCA")
    X_g27, y_g27 = load("GSE27155")
    X_g126, y_g126 = load("GSE126698")
    return genes, (X_tcga, y_tcga), (X_g27, y_g27), (X_g126, y_g126)


def run():
    genes, (X_tcga, y_tcga), (X_g27, y_g27), (X_g126, y_g126) = load_all()
    rows = []

    for strat in STRATEGIES:
        for k in KS:
            if k > len(genes):
                continue
            # 5-fold CV on TCGA
            skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
            fold_aucs = []
            ext_aucs_27 = []
            ext_aucs_126 = []
            ext_npvs_27 = []
            ext_npvs_126 = []
            sel = SELECTORS[strat]
            for tr, te in skf.split(X_tcga, y_tcga):
                idx = sel(X_tcga[tr], y_tcga[tr], genes, k)
                idx = np.asarray(idx, dtype=int)
                Xtr = X_tcga[tr][:, idx]
                Xte = X_tcga[te][:, idx]
                sc = StandardScaler().fit(Xtr)
                clf = LogisticRegression(penalty="l2", max_iter=1000,
                                         class_weight="balanced",
                                         solver="liblinear")
                clf.fit(sc.transform(Xtr), y_tcga[tr])
                p_te = clf.predict_proba(sc.transform(Xte))[:, 1]
                try:
                    fold_aucs.append(roc_auc_score(y_tcga[te], p_te))
                except Exception:
                    pass
                # external
                def ext_eval(X_ext, y_ext):
                    if X_ext is None or y_ext is None or len(np.unique(y_ext)) < 2:
                        return np.nan, np.nan
                    X_e = X_ext[:, idx]
                    p_e_raw = clf.predict_proba(sc.transform(X_e))[:, 1]
                    # isotonic recalibrate on half of external
                    rng = np.random.default_rng(42)
                    perm = rng.permutation(len(y_ext))
                    half = len(y_ext) // 2
                    cal_idx = perm[:half]
                    eval_idx = perm[half:]
                    try:
                        iso = IsotonicRegression(out_of_bounds="clip")
                        iso.fit(p_e_raw[cal_idx], y_ext[cal_idx])
                        p_e_iso = iso.transform(p_e_raw[eval_idx])
                    except Exception:
                        p_e_iso = p_e_raw[eval_idx]
                    try:
                        auc_e = roc_auc_score(y_ext[eval_idx], p_e_iso)
                    except Exception:
                        auc_e = np.nan
                    npv_e, _ = npv_at_sens(y_ext[eval_idx], p_e_iso)
                    return auc_e, npv_e

                a27, n27 = ext_eval(X_g27, y_g27)
                a126, n126 = ext_eval(X_g126, y_g126)
                if not np.isnan(a27):
                    ext_aucs_27.append(a27)
                if not np.isnan(a126):
                    ext_aucs_126.append(a126)
                if not np.isnan(n27):
                    ext_npvs_27.append(n27)
                if not np.isnan(n126):
                    ext_npvs_126.append(n126)

            fold_aucs = np.array(fold_aucs)
            tcga_mean = float(np.mean(fold_aucs)) if fold_aucs.size else np.nan
            lo, hi = bootstrap_ci(fold_aucs) if fold_aucs.size else (np.nan, np.nan)
            rows.append(dict(
                strategy=strat, k=k,
                tcga_mean_auc=tcga_mean,
                tcga_auc_lo=lo, tcga_auc_hi=hi,
                gse27155_mean_auc=float(np.mean(ext_aucs_27)) if ext_aucs_27 else np.nan,
                gse126698_mean_auc=float(np.mean(ext_aucs_126)) if ext_aucs_126 else np.nan,
                gse27155_mean_npv_se95=float(np.mean(ext_npvs_27)) if ext_npvs_27 else np.nan,
                gse126698_mean_npv_se95=float(np.mean(ext_npvs_126)) if ext_npvs_126 else np.nan,
                n_tcga=int(len(y_tcga)),
                n_gse27155=int(len(y_g27)) if y_g27 is not None else 0,
                n_gse126698=int(len(y_g126)) if y_g126 is not None else 0,
            ))
    df = pd.DataFrame(rows)
    out = RESULTS_TABLES / "v3_panel_size_curve.tsv"
    write_tsv(df, out)

    # plateau detection: smallest k where TCGA AUC within 0.005 of max per strat
    plateau = {}
    for strat in STRATEGIES:
        sub = df[df.strategy == strat].sort_values("k")
        if sub.empty:
            continue
        max_auc = sub["tcga_mean_auc"].max()
        plateau[strat] = int(sub[sub.tcga_mean_auc >= max_auc - 0.005].k.min())

    # Identify best (strategy, k)
    best = df.sort_values(by=["gse27155_mean_auc"], ascending=False).head(1).iloc[0]
    best_info = dict(strategy=best["strategy"], k=int(best["k"]),
                     gse27155_auc=float(best["gse27155_mean_auc"]),
                     tcga_auc=float(best["tcga_mean_auc"]))
    (RESULTS_TABLES / "v3_panel_best.json").write_text(
        json.dumps(best_info, indent=2))

    # Figures
    fig = make_subplots(rows=3, cols=1, shared_xaxes=True,
                         subplot_titles=("TCGA 5-fold CV AUC",
                                         "GSE27155 external AUC (isotonic recal)",
                                         "GSE126698 external AUC (isotonic recal)"))
    colors = {"univariate_d": "#5eead4", "tds16_union": "#38bdf8",
              "qubo_neal": "#f472b6"}
    for strat in STRATEGIES:
        sub = df[df.strategy == strat].sort_values("k")
        fig.add_trace(go.Scatter(x=sub["k"], y=sub["tcga_mean_auc"],
                                  mode="lines+markers", name=f"{strat} (TCGA)",
                                  line=dict(color=colors[strat]),
                                  legendgroup=strat), row=1, col=1)
        fig.add_trace(go.Scatter(x=sub["k"], y=sub["gse27155_mean_auc"],
                                  mode="lines+markers", name=f"{strat} (G27)",
                                  line=dict(color=colors[strat], dash="dash"),
                                  legendgroup=strat, showlegend=False),
                       row=2, col=1)
        fig.add_trace(go.Scatter(x=sub["k"], y=sub["gse126698_mean_auc"],
                                  mode="lines+markers", name=f"{strat} (G126)",
                                  line=dict(color=colors[strat], dash="dot"),
                                  legendgroup=strat, showlegend=False),
                       row=3, col=1)
    for strat, k_plat in plateau.items():
        fig.add_vline(x=k_plat, line=dict(dash="dot",
                                           color=colors.get(strat, "#94a3b8")),
                       annotation_text=f"plateau@{strat}", row=1, col=1)
    fig.update_layout(title="Panel size curve (3 strategies, 3 cohorts)",
                      template="plotly_dark", height=820,
                      xaxis3_title="Panel size k")
    fig.write_html(FIGS / "v3_panel_size_curve.html", include_plotlyjs="cdn")

    # NPV heatmap
    piv = df.pivot_table(index="strategy", columns="k",
                          values="gse27155_mean_npv_se95")
    fig2 = go.Figure(data=go.Heatmap(z=piv.values, x=piv.columns, y=piv.index,
                                      colorscale="Teal",
                                      colorbar=dict(title="NPV@Se≥0.95")))
    fig2.update_layout(title="NPV at Se≥0.95 (GSE27155)",
                       template="plotly_dark", height=360)
    fig2.write_html(FIGS / "v3_panel_size_npv_heatmap.html",
                     include_plotlyjs="cdn")

    return dict(rows=len(df), best=best_info, out=str(out), plateau=plateau)


def main():
    res = run()
    print(json.dumps(res, indent=2, default=str))


if __name__ == "__main__":
    main()
