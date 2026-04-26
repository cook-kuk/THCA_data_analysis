"""v3 STEP 2 — Honesty audit.

Tasks A-F, each callable individually, with a main() that runs all in a single
process. The orchestrator parallelises tasks across processes.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy.stats import spearmanr
from sklearn.calibration import CalibratedClassifierCV, calibration_curve
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (average_precision_score, balanced_accuracy_score,
                             brier_score_loss, roc_auc_score, roc_curve)
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, "/opt/thyroid-dash/project/notebooks_or_scripts")
from v3_common import (FIGS, MAPK_OUTPUT, RESULTS_TABLES, bootstrap_ci,
                        build_feature_matrix, get_tds16, get_tierA67_clean,
                        small_n_warning, write_tsv, _expr_sample_id_column)


RNG = np.random.default_rng(42)


# ---------- helpers ----------
def cohens_d(x_pos: np.ndarray, x_neg: np.ndarray) -> float:
    n1, n2 = len(x_pos), len(x_neg)
    if n1 < 2 or n2 < 2:
        return 0.0
    s1, s2 = x_pos.var(ddof=1), x_neg.var(ddof=1)
    denom = np.sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / max(n1 + n2 - 2, 1))
    if denom == 0:
        return 0.0
    return float((x_pos.mean() - x_neg.mean()) / denom)


def cv_auc(X, y, n_splits: int = 5, seed: int = 42):
    if len(np.unique(y)) < 2 or y.sum() < 2 or (len(y) - y.sum()) < 2:
        return float("nan"), np.array([])
    skf = StratifiedKFold(n_splits=min(n_splits, int(min(y.sum(), len(y) - y.sum()))),
                          shuffle=True, random_state=seed)
    proba = np.zeros(len(y), dtype=float)
    fold_aucs = []
    for tr, te in skf.split(X, y):
        sc = StandardScaler()
        X_tr = sc.fit_transform(X[tr])
        X_te = sc.transform(X[te])
        clf = LogisticRegression(penalty="l2", max_iter=1000, class_weight="balanced",
                                 solver="liblinear")
        clf.fit(X_tr, y[tr])
        p = clf.predict_proba(X_te)[:, 1]
        proba[te] = p
        try:
            fold_aucs.append(roc_auc_score(y[te], p))
        except Exception:
            pass
    try:
        overall = roc_auc_score(y, proba)
    except Exception:
        overall = float("nan")
    return overall, np.array(fold_aucs)


def bootstrap_auc_ci(X, y, n_boot=1000, seed=42):
    X = np.asarray(X)
    y = np.asarray(y)
    rng = np.random.default_rng(seed)
    aucs = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(y), len(y))
        yb, Xb = y[idx], X[idx]
        if len(np.unique(yb)) < 2:
            continue
        sc = StandardScaler()
        Xb_s = sc.fit_transform(Xb)
        clf = LogisticRegression(penalty="l2", max_iter=500, solver="liblinear")
        clf.fit(Xb_s, yb)
        try:
            p = clf.predict_proba(Xb_s)[:, 1]
            aucs.append(roc_auc_score(yb, p))
        except Exception:
            pass
    if not aucs:
        return float("nan"), float("nan"), float("nan")
    lo, hi = np.percentile(aucs, [2.5, 97.5])
    return float(np.mean(aucs)), float(lo), float(hi)


# ---------- TASK A: leakage curve ----------
def task_leakage_curve():
    genes = get_tierA67_clean()
    X_df, y = build_feature_matrix("TCGA-THCA", genes)
    X = X_df.to_numpy(float)
    y_arr = y.to_numpy(int)

    # rank genes by |Cohen's d|
    d_abs = np.array([abs(cohens_d(X[y_arr == 1, i], X[y_arr == 0, i]))
                      for i in range(X.shape[1])])
    order = np.argsort(-d_abs)

    rows = []
    for k in [0, 3, 5, 10, 15, 20, 25, 30]:
        keep = np.array([i for i in range(X.shape[1]) if i not in set(order[:k])])
        if keep.size == 0:
            rows.append(dict(k_removed=k, mean_auc=np.nan, auc_lo=np.nan,
                              auc_hi=np.nan, n_features=0))
            continue
        Xk = X[:, keep]
        mean_auc, ci_lo, ci_hi = bootstrap_auc_ci(Xk, y_arr, n_boot=1000, seed=42 + k)
        rows.append(dict(k_removed=k, mean_auc=mean_auc, auc_lo=ci_lo, auc_hi=ci_hi,
                          n_features=int(keep.size),
                          removed_genes=",".join(X_df.columns[order[:k]])))
    df = pd.DataFrame(rows)
    out = RESULTS_TABLES / "v3_leakage_curve.tsv"
    write_tsv(df, out)

    fig = go.Figure()
    fig.add_scatter(x=df["k_removed"], y=df["mean_auc"], mode="lines+markers",
                    name="Mean AUC", line=dict(color="#5eead4"))
    fig.add_scatter(x=df["k_removed"], y=df["auc_hi"], mode="lines", name="95% CI hi",
                    line=dict(dash="dot", color="#5eead4"))
    fig.add_scatter(x=df["k_removed"], y=df["auc_lo"], mode="lines", name="95% CI lo",
                    line=dict(dash="dot", color="#5eead4"),
                    fill="tonexty", fillcolor="rgba(94,234,212,0.15)")
    fig.update_layout(title="Leakage curve: TierA67_clean TCGA AUC vs top-k removed",
                      xaxis_title="Top-k |Cohen's d| genes removed",
                      yaxis_title="Bootstrap AUC (1000 iter)",
                      template="plotly_dark", height=420)
    fig.write_html(FIGS / "v3_leakage_curve.html", include_plotlyjs="cdn")
    return {"task": "A_leakage", "rows": len(df), "out": str(out)}


# ---------- TASK B: MAPK-output ablation ----------
def task_mapk_ablation():
    genes = get_tierA67_clean()
    X_df, y = build_feature_matrix("TCGA-THCA", genes)
    y_arr = y.to_numpy(int)

    rows = []
    for label, keep_cols in [("full_TierA67_clean", list(X_df.columns)),
                              ("minus_MAPK_OUTPUT",
                               [c for c in X_df.columns if c not in set(MAPK_OUTPUT)])]:
        X_sub = X_df[keep_cols].to_numpy(float)
        auc_cv, fold_aucs = cv_auc(X_sub, y_arr)
        lo, hi = bootstrap_ci(fold_aucs) if fold_aucs.size else (np.nan, np.nan)
        rows.append(dict(variant=label, n_features=len(keep_cols),
                          mean_cv_auc=auc_cv, cv_auc_lo=lo, cv_auc_hi=hi))
    df = pd.DataFrame(rows)
    out = RESULTS_TABLES / "v3_mapk_ablation.tsv"
    write_tsv(df, out)

    fig = go.Figure()
    fig.add_bar(x=df["variant"], y=df["mean_cv_auc"],
                error_y=dict(type="data",
                             array=df["cv_auc_hi"] - df["mean_cv_auc"],
                             arrayminus=df["mean_cv_auc"] - df["cv_auc_lo"]),
                marker_color=["#38bdf8", "#f472b6"])
    fig.update_layout(title="MAPK-output ablation (TCGA 5-fold CV AUC)",
                      yaxis_title="AUC (fold-mean, 95% CI)", yaxis_range=[0.5, 1.02],
                      template="plotly_dark", height=420)
    fig.write_html(FIGS / "v3_mapk_ablation.html", include_plotlyjs="cdn")
    return {"task": "B_mapk_ablation", "rows": len(df), "out": str(out)}


# ---------- TASK C: permutation null ----------
def task_permutation_null(n_perm: int = 1000):
    genes = get_tierA67_clean()
    X_df, y = build_feature_matrix("TCGA-THCA", genes)
    X = X_df.to_numpy(float)
    y_arr = y.to_numpy(int)

    real_auc, _ = cv_auc(X, y_arr)

    rng = np.random.default_rng(0)
    null_aucs = []
    for i in range(n_perm):
        yp = y_arr.copy()
        rng.shuffle(yp)
        sc = StandardScaler()
        Xs = sc.fit_transform(X)
        clf = LogisticRegression(penalty="l2", max_iter=400, solver="liblinear")
        clf.fit(Xs, yp)
        try:
            auc = roc_auc_score(yp, clf.predict_proba(Xs)[:, 1])
        except Exception:
            auc = 0.5
        null_aucs.append(auc)
    null_aucs = np.array(null_aucs)
    p_emp = float((np.sum(null_aucs >= real_auc) + 1) / (n_perm + 1))

    df = pd.DataFrame(dict(iteration=np.arange(n_perm), null_auc=null_aucs))
    df_summary = pd.DataFrame([dict(real_auc=real_auc,
                                     null_mean=float(np.mean(null_aucs)),
                                     null_p5=float(np.percentile(null_aucs, 5)),
                                     null_p95=float(np.percentile(null_aucs, 95)),
                                     empirical_p=p_emp,
                                     n_perm=n_perm)])
    out = RESULTS_TABLES / "v3_permutation_null.tsv"
    write_tsv(df_summary, out)
    (RESULTS_TABLES / "v3_permutation_null_samples.tsv").write_text(
        df.to_csv(sep="\t", index=False))

    fig = go.Figure()
    fig.add_histogram(x=null_aucs, nbinsx=40, marker_color="#64748b",
                      name="Shuffled labels")
    fig.add_vline(x=real_auc, line=dict(color="#f472b6", width=3),
                   annotation_text=f"Real AUC = {real_auc:.3f}")
    fig.update_layout(title=f"Permutation null (n={n_perm}, empirical p={p_emp:.4f})",
                      xaxis_title="Training AUC under shuffled labels",
                      yaxis_title="Frequency",
                      template="plotly_dark", height=420)
    fig.write_html(FIGS / "v3_permutation_null.html", include_plotlyjs="cdn")
    return {"task": "C_permutation", "real_auc": real_auc, "p": p_emp,
            "out": str(out)}


# ---------- TASK D: dataset identifiability ----------
def task_dataset_identifiability():
    from v3_common import load_expr, DATASET_FILES
    genes = get_tierA67_clean()
    sm_master = pd.read_csv("/opt/thyroid-dash/project/metadata/sample_master.tsv",
                            sep="\t")

    rows = []
    X_parts = []
    y_parts = []
    for ds in ["TCGA-THCA", "GSE27155", "GSE76039", "GSE126698"]:
        expr = load_expr(DATASET_FILES[ds])
        sm = sm_master[(sm_master.dataset == ds) &
                       (sm_master.normal_vs_tumor == "tumor")].copy()
        sm["_expr_id"] = _expr_sample_id_column(ds, sm)
        samples = [s for s in sm["_expr_id"] if s in expr.columns]
        if not samples:
            continue
        X = expr.reindex(genes).fillna(0.0).T.loc[samples].to_numpy(float)
        X_parts.append(X)
        y_parts.extend([ds] * len(samples))

    X_all = np.vstack(X_parts)
    y_all = np.array(y_parts)
    classes = list(dict.fromkeys(y_all.tolist()))

    # one-vs-rest AUC across 5-fold
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    ovr_auc = {c: [] for c in classes}
    for tr, te in skf.split(X_all, y_all):
        sc = StandardScaler()
        Xtr = sc.fit_transform(X_all[tr])
        Xte = sc.transform(X_all[te])
        for c in classes:
            y_tr_bin = (y_all[tr] == c).astype(int)
            y_te_bin = (y_all[te] == c).astype(int)
            if y_tr_bin.sum() < 2 or (len(y_tr_bin) - y_tr_bin.sum()) < 2:
                continue
            clf = LogisticRegression(penalty="l2", max_iter=500,
                                     solver="liblinear")
            clf.fit(Xtr, y_tr_bin)
            p = clf.predict_proba(Xte)[:, 1]
            try:
                ovr_auc[c].append(roc_auc_score(y_te_bin, p))
            except Exception:
                pass

    for c in classes:
        a = np.array(ovr_auc[c])
        rows.append(dict(dataset=c, n_samples=int((y_all == c).sum()),
                          ovr_auc_mean=float(np.mean(a)) if a.size else np.nan,
                          ovr_auc_std=float(np.std(a)) if a.size else np.nan,
                          ovr_auc_folds=",".join([f"{x:.3f}" for x in a])))

    df = pd.DataFrame(rows)
    out = RESULTS_TABLES / "v3_dataset_identifiability.tsv"
    write_tsv(df, out)

    fig = go.Figure()
    fig.add_bar(x=df["dataset"], y=df["ovr_auc_mean"],
                error_y=dict(type="data", array=df["ovr_auc_std"]),
                marker_color="#fb7185")
    fig.add_hline(y=0.5, line=dict(dash="dash", color="#94a3b8"),
                   annotation_text="chance")
    fig.update_layout(title="Dataset identifiability (OvR LogReg, TierA67_clean)",
                      yaxis_title="AUC (5-fold OvR)", yaxis_range=[0.3, 1.02],
                      template="plotly_dark", height=420)
    fig.write_html(FIGS / "v3_dataset_identifiability.html",
                    include_plotlyjs="cdn")
    return {"task": "D_identifiability", "rows": len(df), "out": str(out)}


# ---------- TASK E: LODO ----------
def lodo_metrics(y_true, y_score):
    if len(np.unique(y_true)) < 2:
        return dict(auc=np.nan, pr_auc=np.nan, bacc=np.nan,
                     npv_at_se95=np.nan, brier=np.nan, youden=np.nan)
    auc = roc_auc_score(y_true, y_score)
    pr_auc = average_precision_score(y_true, y_score)
    fpr, tpr, thr = roc_curve(y_true, y_score)
    j = tpr - fpr
    k = int(np.argmax(j))
    youden = float(thr[k])
    y_pred = (y_score >= youden).astype(int)
    bacc = balanced_accuracy_score(y_true, y_pred)
    # NPV at Se>=0.95
    valid = tpr >= 0.95
    npv_val = np.nan
    if valid.any():
        idx = np.argmax(valid)  # first true
        t = thr[idx]
        yp = (y_score >= t).astype(int)
        tn = int(((y_true == 0) & (yp == 0)).sum())
        fn = int(((y_true == 1) & (yp == 0)).sum())
        if (tn + fn) > 0:
            npv_val = tn / (tn + fn)
    brier = brier_score_loss(y_true, y_score)
    return dict(auc=float(auc), pr_auc=float(pr_auc), bacc=float(bacc),
                 npv_at_se95=float(npv_val) if not np.isnan(npv_val) else np.nan,
                 brier=float(brier), youden=youden)


def task_lodo():
    from v3_common import load_expr, DATASET_FILES
    genes = [g for g in get_tierA67_clean() if g not in set(MAPK_OUTPUT)]

    # Build matrices (tumor-only, BRAF_like/RAS_like) per dataset
    def build_ds(ds):
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

    external = ["GSE27155", "GSE126698", "GSE213647"]
    X_tcga, y_tcga = build_ds("TCGA-THCA")

    rows = []
    ext_data = {}
    for ds in external:
        X, y = build_ds(ds)
        ext_data[ds] = (X, y)

    for c in external:
        X_c, y_c = ext_data[c]
        if X_c is None or y_c is None or len(np.unique(y_c)) < 2:
            rows.append(dict(external=c, mode="train_tcga_only",
                             n_test=0 if y_c is None else int(len(y_c)),
                             **{k: np.nan for k in
                                ["auc", "pr_auc", "bacc",
                                 "npv_at_se95", "brier", "youden"]},
                             note="skipped: <2 classes"))
            rows.append(dict(external=c, mode="train_tcga_plus_others",
                             n_test=0 if y_c is None else int(len(y_c)),
                             **{k: np.nan for k in
                                ["auc", "pr_auc", "bacc",
                                 "npv_at_se95", "brier", "youden"]},
                             note="skipped: <2 classes"))
            continue

        # (1) train TCGA, test C
        sc = StandardScaler().fit(X_tcga)
        clf = LogisticRegression(penalty="l2", max_iter=1000,
                                 class_weight="balanced", solver="liblinear")
        clf.fit(sc.transform(X_tcga), y_tcga)
        p_c = clf.predict_proba(sc.transform(X_c))[:, 1]
        m1 = lodo_metrics(y_c, p_c)
        rows.append(dict(external=c, mode="train_tcga_only",
                          n_test=int(len(y_c)), **m1, note=""))

        # (2) train TCGA + others\C, test C
        X_list = [X_tcga]
        y_list = [y_tcga]
        for other in external:
            if other == c:
                continue
            Xo, yo = ext_data[other]
            if Xo is not None and yo is not None and len(np.unique(yo)) >= 2:
                X_list.append(Xo)
                y_list.append(yo)
        X_tr = np.vstack(X_list)
        y_tr = np.concatenate(y_list)
        sc = StandardScaler().fit(X_tr)
        clf2 = LogisticRegression(penalty="l2", max_iter=1000,
                                  class_weight="balanced", solver="liblinear")
        clf2.fit(sc.transform(X_tr), y_tr)
        p_c2 = clf2.predict_proba(sc.transform(X_c))[:, 1]
        m2 = lodo_metrics(y_c, p_c2)
        rows.append(dict(external=c, mode="train_tcga_plus_others",
                          n_test=int(len(y_c)), **m2,
                          note=small_n_warning(int(len(y_c)))))

    df = pd.DataFrame(rows)
    out = RESULTS_TABLES / "v3_lodo.tsv"
    write_tsv(df, out)

    fig = go.Figure()
    piv = df.pivot_table(index="external", columns="mode", values="auc",
                          aggfunc="first")
    colors = {"train_tcga_only": "#38bdf8",
              "train_tcga_plus_others": "#f472b6"}
    for mode in piv.columns:
        fig.add_bar(name=mode, x=piv.index, y=piv[mode],
                    marker_color=colors.get(mode, "#a3a3a3"))
    fig.add_hline(y=0.5, line=dict(dash="dash", color="#94a3b8"))
    fig.update_layout(barmode="group",
                      title="LODO AUC per external cohort",
                      yaxis_title="AUC", yaxis_range=[0.0, 1.05],
                      template="plotly_dark", height=440)
    fig.write_html(FIGS / "v3_lodo.html", include_plotlyjs="cdn")
    return {"task": "E_lodo", "rows": len(df), "out": str(out)}


# ---------- TASK F: GSE27155 calibration + decision curve ----------
def task_calibration_gse27155():
    from v3_common import load_expr, DATASET_FILES
    genes = [g for g in get_tierA67_clean() if g not in set(MAPK_OUTPUT)]

    # Train on TCGA
    expr_t = load_expr(DATASET_FILES["TCGA-THCA"])
    sm = pd.read_csv("/opt/thyroid-dash/project/metadata/sample_master.tsv",
                     sep="\t")
    sm_t = sm[(sm.dataset == "TCGA-THCA") & (sm.normal_vs_tumor == "tumor")
               & sm.molecular_subtype.isin(["BRAF_like", "RAS_like"])].copy()
    sm_t["_expr_id"] = _expr_sample_id_column("TCGA-THCA", sm_t)
    sm_t = sm_t[sm_t["_expr_id"].isin(expr_t.columns)]
    samps_t = sm_t["_expr_id"].tolist()
    X_tcga = expr_t.reindex(genes).fillna(0.0).T.loc[samps_t].to_numpy(float)
    y_tcga = (sm_t.set_index("_expr_id").loc[samps_t, "molecular_subtype"]
                  .map({"BRAF_like": 1, "RAS_like": 0}).to_numpy(int))

    expr_g = load_expr(DATASET_FILES["GSE27155"])
    sm_g = sm[(sm.dataset == "GSE27155") & (sm.normal_vs_tumor == "tumor")
               & sm.molecular_subtype.isin(["BRAF_like", "RAS_like"])].copy()
    sm_g["_expr_id"] = _expr_sample_id_column("GSE27155", sm_g)
    sm_g = sm_g[sm_g["_expr_id"].isin(expr_g.columns)]
    samps_g = sm_g["_expr_id"].tolist()
    X_g = expr_g.reindex(genes).fillna(0.0).T.loc[samps_g].to_numpy(float)
    y_g = (sm_g.set_index("_expr_id").loc[samps_g, "molecular_subtype"]
               .map({"BRAF_like": 1, "RAS_like": 0}).to_numpy(int))

    sc = StandardScaler().fit(X_tcga)
    base = LogisticRegression(penalty="l2", max_iter=1000,
                              class_weight="balanced", solver="liblinear")
    base.fit(sc.transform(X_tcga), y_tcga)
    p_g_uncal = base.predict_proba(sc.transform(X_g))[:, 1]

    frac_pos_unc, mean_pred_unc = calibration_curve(y_g, p_g_uncal, n_bins=10,
                                                     strategy="quantile")
    brier_unc = brier_score_loss(y_g, p_g_uncal)
    ece_unc = float(np.abs(frac_pos_unc - mean_pred_unc).mean())

    # Isotonic recalibration on GSE27155 5-fold
    from sklearn.isotonic import IsotonicRegression
    from sklearn.model_selection import StratifiedKFold
    p_iso = np.zeros(len(y_g))
    try:
        skf = StratifiedKFold(n_splits=min(5, int(y_g.sum())), shuffle=True,
                               random_state=42)
        for tr, te in skf.split(X_g, y_g):
            iso = IsotonicRegression(out_of_bounds="clip")
            iso.fit(p_g_uncal[tr], y_g[tr])
            p_iso[te] = iso.transform(p_g_uncal[te])
    except Exception:
        p_iso = p_g_uncal.copy()

    from sklearn.linear_model import LogisticRegression as _LR
    p_platt = np.zeros(len(y_g))
    try:
        skf = StratifiedKFold(n_splits=min(5, int(y_g.sum())), shuffle=True,
                               random_state=42)
        for tr, te in skf.split(X_g, y_g):
            platt = _LR(max_iter=500)
            platt.fit(p_g_uncal[tr].reshape(-1, 1), y_g[tr])
            p_platt[te] = platt.predict_proba(
                p_g_uncal[te].reshape(-1, 1))[:, 1]
    except Exception:
        p_platt = p_g_uncal.copy()

    frac_pos_iso, mean_pred_iso = calibration_curve(y_g, p_iso, n_bins=10,
                                                     strategy="quantile")
    brier_iso = brier_score_loss(y_g, p_iso)
    ece_iso = float(np.abs(frac_pos_iso - mean_pred_iso).mean())

    frac_pos_platt, mean_pred_platt = calibration_curve(y_g, p_platt, n_bins=10,
                                                         strategy="quantile")
    brier_platt = brier_score_loss(y_g, p_platt)
    ece_platt = float(np.abs(frac_pos_platt - mean_pred_platt).mean())

    calib_rows = []
    for mode, brier, ece, mean_pred, frac_pos in [
        ("uncal", brier_unc, ece_unc, mean_pred_unc, frac_pos_unc),
        ("isotonic", brier_iso, ece_iso, mean_pred_iso, frac_pos_iso),
        ("platt", brier_platt, ece_platt, mean_pred_platt, frac_pos_platt),
    ]:
        for mp, fp in zip(mean_pred, frac_pos):
            calib_rows.append(dict(mode=mode, bin_mean_pred=float(mp),
                                    bin_frac_pos=float(fp),
                                    brier=float(brier), ece=float(ece),
                                    n_samples=int(len(y_g))))
    df = pd.DataFrame(calib_rows)
    out = RESULTS_TABLES / "v3_calibration_gse27155.tsv"
    write_tsv(df, out)

    fig = go.Figure()
    fig.add_scatter(x=[0, 1], y=[0, 1], mode="lines",
                     line=dict(dash="dash", color="#94a3b8"),
                     name="Perfect calibration")
    for mode, brier, ece, mean_pred, frac_pos, color in [
        ("uncal", brier_unc, ece_unc, mean_pred_unc, frac_pos_unc, "#f472b6"),
        ("isotonic", brier_iso, ece_iso, mean_pred_iso, frac_pos_iso, "#5eead4"),
        ("platt", brier_platt, ece_platt, mean_pred_platt, frac_pos_platt, "#38bdf8"),
    ]:
        fig.add_scatter(x=mean_pred, y=frac_pos, mode="lines+markers",
                         name=f"{mode} (Brier={brier:.3f}, ECE={ece:.3f})",
                         line=dict(color=color))
    fig.update_layout(
        title=f"Reliability: TCGA→GSE27155 (n={len(y_g)}{small_n_warning(len(y_g))})",
        xaxis_title="Predicted probability", yaxis_title="Observed frequency",
        template="plotly_dark", height=460)
    fig.write_html(FIGS / "v3_calibration_gse27155.html", include_plotlyjs="cdn")

    # Decision curve: threshold sweep + net benefit
    dc_rows = []
    for t in np.arange(0.05, 0.501, 0.01):
        yp = (p_iso >= t).astype(int)
        tp = int(((y_g == 1) & (yp == 1)).sum())
        fp = int(((y_g == 0) & (yp == 1)).sum())
        tn = int(((y_g == 0) & (yp == 0)).sum())
        fn = int(((y_g == 1) & (yp == 0)).sum())
        n = len(y_g)
        if n == 0:
            continue
        # Net benefit = TP/N - FP/N * (t/(1-t))
        nb_model = tp / n - fp / n * (t / max(1 - t, 1e-6))
        prev = y_g.mean()
        nb_all = prev - (1 - prev) * (t / max(1 - t, 1e-6))
        nb_none = 0.0
        sens = tp / (tp + fn) if (tp + fn) > 0 else np.nan
        spec = tn / (tn + fp) if (tn + fp) > 0 else np.nan
        dc_rows.append(dict(threshold=float(t), tp=tp, fp=fp, tn=tn, fn=fn,
                             sens=sens, spec=spec,
                             nb_model=float(nb_model), nb_treat_all=float(nb_all),
                             nb_treat_none=float(nb_none)))
    dc_df = pd.DataFrame(dc_rows)
    out2 = RESULTS_TABLES / "v3_decision_curve_gse27155.tsv"
    write_tsv(dc_df, out2)

    fig2 = go.Figure()
    fig2.add_scatter(x=dc_df["threshold"], y=dc_df["nb_model"],
                      mode="lines", name="Isotonic model",
                      line=dict(color="#5eead4", width=3))
    fig2.add_scatter(x=dc_df["threshold"], y=dc_df["nb_treat_all"],
                      mode="lines", name="Treat all",
                      line=dict(color="#f472b6", dash="dot"))
    fig2.add_scatter(x=dc_df["threshold"], y=dc_df["nb_treat_none"],
                      mode="lines", name="Treat none",
                      line=dict(color="#94a3b8", dash="dot"))
    fig2.update_layout(
        title=f"Decision curve (GSE27155, n={len(y_g)}{small_n_warning(len(y_g))})",
        xaxis_title="Threshold probability",
        yaxis_title="Net benefit",
        template="plotly_dark", height=420)
    fig2.write_html(FIGS / "v3_decision_curve_gse27155.html",
                     include_plotlyjs="cdn")
    return {"task": "F_calibration", "out": str(out), "out2": str(out2),
            "n_test": int(len(y_g)), "brier_iso": brier_iso, "ece_iso": ece_iso}


TASKS = {
    "A": task_leakage_curve,
    "B": task_mapk_ablation,
    "C": task_permutation_null,
    "D": task_dataset_identifiability,
    "E": task_lodo,
    "F": task_calibration_gse27155,
}


def run_task(key: str):
    fn = TASKS[key]
    try:
        result = fn()
        return {"key": key, "ok": True, **result}
    except Exception as e:
        import traceback
        tb = traceback.format_exc()
        return {"key": key, "ok": False, "error": str(e), "tb": tb}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--task", default="all", choices=list(TASKS.keys()) + ["all"])
    args = p.parse_args()
    if args.task == "all":
        results = [run_task(k) for k in TASKS]
    else:
        results = [run_task(args.task)]
    print(json.dumps(results, indent=2, default=str))


if __name__ == "__main__":
    main()
