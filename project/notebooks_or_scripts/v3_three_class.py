"""v3 STEP 4 — Three-class (BRAF_like vs RAS_like vs dedifferentiated).

Multinomial LogReg + GradientBoosting, 5-fold CV, per-class OvR AUC,
macro AUC, confusion, per-class P/R/F1. SHAP mean|SHAP| top-15 per class.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (classification_report, confusion_matrix,
                             roc_auc_score)
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, "/opt/thyroid-dash/project/notebooks_or_scripts")
from v3_common import (DATASET_FILES, FIGS, MAPK_OUTPUT, RESULTS_TABLES,
                        get_tierA67_clean, load_expr, small_n_warning,
                        write_tsv, _expr_sample_id_column)

CLASSES = ["BRAF_like", "RAS_like", "dedifferentiated"]


def load_three_class():
    genes = [g for g in get_tierA67_clean() if g not in set(MAPK_OUTPUT)]
    sm = pd.read_csv("/opt/thyroid-dash/project/metadata/sample_master.tsv",
                     sep="\t")

    X_parts, y_parts = [], []
    for ds in ["TCGA-THCA", "GSE126698", "GSE213647", "GSE76039"]:
        if ds not in DATASET_FILES:
            continue
        expr = load_expr(DATASET_FILES[ds])
        sub = sm[(sm.dataset == ds) & (sm.normal_vs_tumor == "tumor") &
                  sm.molecular_subtype.isin(CLASSES)].copy()
        sub["_expr_id"] = _expr_sample_id_column(ds, sub)
        sub = sub[sub["_expr_id"].isin(expr.columns)]
        if sub.empty:
            continue
        samples = sub["_expr_id"].tolist()
        y = sub.set_index("_expr_id").loc[samples, "molecular_subtype"].to_numpy()
        X = expr.reindex(genes).fillna(0.0).T.loc[samples].to_numpy(float)
        X_parts.append(X)
        y_parts.append(y)
    X = np.vstack(X_parts)
    y = np.concatenate(y_parts)
    return X, y, genes


def cv_three_class(X, y, genes, model_name="logreg", n_splits=5, seed=42):
    classes = CLASSES
    class_to_idx = {c: i for i, c in enumerate(classes)}
    y_int = np.array([class_to_idx[v] for v in y])

    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    proba = np.zeros((len(y), len(classes)))
    preds = np.zeros(len(y), dtype=int)

    for tr, te in skf.split(X, y_int):
        sc = StandardScaler().fit(X[tr])
        Xtr = sc.transform(X[tr])
        Xte = sc.transform(X[te])
        if model_name == "logreg":
            clf = LogisticRegression(penalty="l2", max_iter=1500,
                                     class_weight="balanced",
                                     solver="lbfgs")
        else:
            clf = GradientBoostingClassifier(n_estimators=120, max_depth=3,
                                             random_state=seed)
        clf.fit(Xtr, y_int[tr])
        p = clf.predict_proba(Xte)
        # Align columns to class order
        cls = clf.classes_
        for j, c in enumerate(cls):
            proba[te, c] = p[:, j]
        preds[te] = clf.predict(Xte)

    per_class = {}
    for i, c in enumerate(classes):
        y_bin = (y_int == i).astype(int)
        try:
            per_class[c] = float(roc_auc_score(y_bin, proba[:, i]))
        except Exception:
            per_class[c] = float("nan")

    macro_auc = float(np.nanmean(list(per_class.values())))

    conf = confusion_matrix(y_int, preds, labels=list(range(len(classes))))
    rep = classification_report(y_int, preds,
                                 labels=list(range(len(classes))),
                                 target_names=classes,
                                 output_dict=True, zero_division=0)

    class_counts = {c: int((y == c).sum()) for c in classes}

    return dict(per_class_auc=per_class, macro_auc=macro_auc,
                 confusion=conf.tolist(), class_names=classes, report=rep,
                 class_counts=class_counts, proba=proba, y=y_int)


def shap_per_class(X, y, genes, seed=42):
    """Compute mean|SHAP| per class on 200-sample stratified holdout."""
    try:
        import shap
    except Exception:
        return None

    class_to_idx = {c: i for i, c in enumerate(CLASSES)}
    y_int = np.array([class_to_idx[v] for v in y])
    rng = np.random.default_rng(seed)

    # 80/20 split
    perm = rng.permutation(len(y))
    cut = int(len(y) * 0.8)
    tr_idx, te_idx = perm[:cut], perm[cut:]
    # cap test at 200
    if len(te_idx) > 200:
        te_idx = te_idx[:200]

    sc = StandardScaler().fit(X[tr_idx])
    # Use RandomForestClassifier for multiclass SHAP (GradientBoostingClassifier
    # TreeExplainer works only for binary problems).
    from sklearn.ensemble import RandomForestClassifier
    clf = RandomForestClassifier(n_estimators=300, max_depth=6,
                                 random_state=seed, n_jobs=-1)
    clf.fit(sc.transform(X[tr_idx]), y_int[tr_idx])

    try:
        explainer = shap.TreeExplainer(clf)
        sv = explainer.shap_values(sc.transform(X[te_idx]))
    except Exception as e:
        print(f"[shap] fallback: {e}")
        return None

    # sklearn TreeExplainer returns list of n_classes or 3D array
    if isinstance(sv, list):
        mats = sv
    elif isinstance(sv, np.ndarray) and sv.ndim == 3:
        mats = [sv[:, :, i] for i in range(sv.shape[2])]
    else:
        mats = [sv]

    out = {}
    for i, c in enumerate(CLASSES):
        if i < len(mats):
            m = np.abs(mats[i]).mean(axis=0)
            idx = np.argsort(-m)[:15]
            out[c] = [(genes[j], float(m[j])) for j in idx]
        else:
            out[c] = []
    return out


def run():
    X, y, genes = load_three_class()
    counts = {c: int((y == c).sum()) for c in CLASSES}
    print(f"[three_class] class counts: {counts}")

    res_logreg = cv_three_class(X, y, genes, model_name="logreg")
    res_gb = cv_three_class(X, y, genes, model_name="gb")

    # Save summary TSV
    summary_rows = []
    for model_name, res in [("logreg", res_logreg), ("gb", res_gb)]:
        for c in CLASSES:
            rep = res["report"].get(c, {})
            summary_rows.append(dict(
                model=model_name, class_name=c,
                n=counts[c],
                ovr_auc=res["per_class_auc"][c],
                precision=rep.get("precision", np.nan),
                recall=rep.get("recall", np.nan),
                f1=rep.get("f1-score", np.nan),
                small_n_flag=("⚠ small n" if counts[c] < 20 else ""),
            ))
        summary_rows.append(dict(model=model_name, class_name="MACRO", n=len(y),
                                 ovr_auc=res["macro_auc"],
                                 precision=np.nan, recall=np.nan, f1=np.nan,
                                 small_n_flag=""))
    df = pd.DataFrame(summary_rows)
    out = RESULTS_TABLES / "v3_three_class_metrics.tsv"
    write_tsv(df, out)

    # Confusion
    conf_rows = []
    for model_name, res in [("logreg", res_logreg), ("gb", res_gb)]:
        conf = res["confusion"]
        for i, c_true in enumerate(CLASSES):
            for j, c_pred in enumerate(CLASSES):
                conf_rows.append(dict(model=model_name, true=c_true,
                                       pred=c_pred, count=int(conf[i][j])))
    conf_df = pd.DataFrame(conf_rows)
    write_tsv(conf_df, RESULTS_TABLES / "v3_three_class_confusion.tsv")

    # ROC OvR figure (logreg)
    from sklearn.metrics import roc_curve
    fig_roc = go.Figure()
    proba = res_logreg["proba"]
    y_int = res_logreg["y"]
    for i, c in enumerate(CLASSES):
        y_bin = (y_int == i).astype(int)
        if y_bin.sum() == 0 or y_bin.sum() == len(y_bin):
            continue
        fpr, tpr, _ = roc_curve(y_bin, proba[:, i])
        fig_roc.add_scatter(x=fpr, y=tpr, mode="lines",
                             name=f"{c} (AUC={res_logreg['per_class_auc'][c]:.3f}, "
                                  f"n={counts[c]}{small_n_warning(counts[c])})")
    fig_roc.add_scatter(x=[0, 1], y=[0, 1], mode="lines",
                         line=dict(dash="dash", color="#94a3b8"),
                         name="chance")
    fig_roc.update_layout(
        title=f"3-class OvR ROC (LogReg, macro AUC={res_logreg['macro_auc']:.3f})",
        xaxis_title="FPR", yaxis_title="TPR",
        template="plotly_dark", height=480)
    fig_roc.write_html(FIGS / "v3_three_class_roc.html",
                        include_plotlyjs="cdn")

    # Confusion heatmap (logreg)
    conf = np.array(res_logreg["confusion"])
    fig_cm = go.Figure(data=go.Heatmap(z=conf, x=CLASSES, y=CLASSES,
                                        text=conf.astype(str),
                                        texttemplate="%{text}",
                                        colorscale="Viridis"))
    fig_cm.update_layout(title="Confusion matrix (LogReg, 5-fold CV)",
                         xaxis_title="Predicted", yaxis_title="True",
                         template="plotly_dark", height=460)
    fig_cm.write_html(FIGS / "v3_three_class_confusion.html",
                      include_plotlyjs="cdn")

    # SHAP per class
    shap_out = shap_per_class(X, y, genes)
    if shap_out is not None:
        shap_rows = []
        for c, features in shap_out.items():
            for rank, (g, v) in enumerate(features, 1):
                shap_rows.append(dict(class_name=c, rank=rank, gene=g,
                                       mean_abs_shap=v))
        shap_df = pd.DataFrame(shap_rows)
        write_tsv(shap_df, RESULTS_TABLES / "v3_three_class_shap.tsv")

        fig_shap = make_subplots(rows=1, cols=len(CLASSES),
                                  subplot_titles=[f"{c} (n={counts[c]}{small_n_warning(counts[c])})"
                                                  for c in CLASSES])
        for i, c in enumerate(CLASSES):
            sub = shap_df[shap_df.class_name == c].sort_values(
                "mean_abs_shap", ascending=True)
            fig_shap.add_trace(go.Bar(x=sub["mean_abs_shap"], y=sub["gene"],
                                       orientation="h",
                                       marker_color=["#5eead4", "#38bdf8",
                                                     "#f472b6"][i]),
                                row=1, col=i + 1)
        fig_shap.update_layout(title="SHAP top-15 per class (GB, 200-sample holdout)",
                                template="plotly_dark", height=600,
                                showlegend=False)
        fig_shap.write_html(FIGS / "v3_three_class_shap.html",
                             include_plotlyjs="cdn")

    return dict(logreg_macro=res_logreg["macro_auc"],
                 gb_macro=res_gb["macro_auc"],
                 class_counts=counts,
                 out=str(out))


def main():
    res = run()
    print(json.dumps(res, indent=2, default=str))


if __name__ == "__main__":
    main()
