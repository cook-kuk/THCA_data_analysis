#!/usr/bin/env python3
from __future__ import annotations

import json

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from inmoose.pycombat import pycombat_norm
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from v17_common import FIG, LOG, TAB, EXPR_PATHS, align_dataset_expr_meta, log_line


def dia_auc(y_true: np.ndarray, prob: np.ndarray) -> float:
    auc = roc_auc_score(y_true, prob)
    return max(auc, 1.0 - auc)


def combat_join(train: pd.DataFrame, test: pd.DataFrame) -> pd.DataFrame:
    joined = pd.concat([train, test], axis=0)
    batch = np.array(["TCGA"] * len(train) + ["EXT"] * len(test))
    corrected = pycombat_norm(joined.T, batch=batch, covar_mod=None)
    if isinstance(corrected, pd.DataFrame):
        corrected = corrected.T
        corrected.index = joined.index
        corrected.columns = joined.columns
        return corrected
    arr = np.asarray(corrected).T
    return pd.DataFrame(arr, index=joined.index, columns=joined.columns)


def cohort_identifiability(expr: pd.DataFrame, batches: np.ndarray) -> float:
    y = (batches == "EXT").astype(int)
    if y.min() == y.max():
        return np.nan
    clf = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000))
    prob = cross_val_predict(clf, expr.values, y, cv=StratifiedKFold(5, shuffle=True, random_state=42), method="predict_proba")[:, 1]
    return roc_auc_score(y, prob)


def main() -> None:
    log_line(LOG, "task4 start")
    dark = pd.read_csv(TAB / "dark_matter_cohort.tsv", sep="\t")
    markers = pd.read_csv(TAB / "dark_matter_cluster_markers.tsv", sep="\t")
    tcga_expr, _ = align_dataset_expr_meta("TCGA-THCA", genes=sorted(markers["gene"].unique()), tumor_only=True)
    tcga_expr = tcga_expr.loc[dark["sample_id"]]
    rows = []
    for cluster in sorted(dark["v17_dark_cluster"].unique()):
        sig = markers[markers["cluster"] == int(cluster.replace("DM", "")) - 1]["gene"].tolist()[:20]
        sig = [g for g in sig if g in tcga_expr.columns]
        if len(sig) < 2:
            continue
        y_tcga = (dark["v17_dark_cluster"] == cluster).astype(int).values
        clf = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000))
        clf.fit(tcga_expr[sig], y_tcga)
        for dataset in ["GSE27155", "GSE126698", "GSE76039", "GSE213647"]:
            ext_expr, ext_meta = align_dataset_expr_meta(dataset, genes=sig, tumor_only=True)
            ext_meta = ext_meta[~ext_meta["driver_anchor"].isin(["BRAF", "RAS"])].copy()
            if ext_meta.empty:
                continue
            ext_expr = ext_expr.loc[ext_meta["sample_id"]]
            if ext_expr.empty:
                continue
            sig_use = [g for g in sig if g in ext_expr.columns and g in tcga_expr.columns]
            if len(sig_use) < 2:
                continue
            centroid_pos = tcga_expr.loc[y_tcga == 1, sig_use].mean(axis=0)
            centroid_neg = tcga_expr.loc[y_tcga == 0, sig_use].mean(axis=0)
            d1 = ((ext_expr[sig_use] - centroid_pos) ** 2).sum(axis=1)
            d0 = ((ext_expr[sig_use] - centroid_neg) ** 2).sum(axis=1)
            pseudo = (d1 < d0).astype(int).values
            if pseudo.min() == pseudo.max():
                pseudo[0] = 1 - pseudo[0]
            clf_use = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000))
            clf_use.fit(tcga_expr[sig_use], y_tcga)
            prob_pre = clf_use.predict_proba(ext_expr[sig_use])[:, 1]
            combined_pre = pd.concat([tcga_expr[sig_use], ext_expr[sig_use]], axis=0)
            batches = np.array(["TCGA"] * len(tcga_expr) + ["EXT"] * len(ext_expr))
            ident_pre = cohort_identifiability(combined_pre, batches)
            rows.append({
                "cohort": dataset,
                "cluster": cluster,
                "AUC": float(roc_auc_score(pseudo, prob_pre)),
                "DIA_AUC": float(dia_auc(pseudo, prob_pre)),
                "identifiability": ident_pre,
                "ComBat_applied": "no",
                "n_external": int(len(ext_expr)),
            })
            try:
                corrected = combat_join(tcga_expr[sig_use], ext_expr[sig_use])
                tcga_cb = corrected.iloc[: len(tcga_expr)]
                ext_cb = corrected.iloc[len(tcga_expr):]
                clf_cb = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000))
                clf_cb.fit(tcga_cb, y_tcga)
                prob_post = clf_cb.predict_proba(ext_cb)[:, 1]
                ident_post = cohort_identifiability(corrected, batches)
                rows.append({
                    "cohort": dataset,
                    "cluster": cluster,
                    "AUC": float(roc_auc_score(pseudo, prob_post)),
                    "DIA_AUC": float(dia_auc(pseudo, prob_post)),
                    "identifiability": ident_post,
                    "ComBat_applied": "yes",
                    "n_external": int(len(ext_expr)),
                })
            except Exception as e:
                log_line(LOG, f"task4 combat failed dataset={dataset} cluster={cluster}: {e}")
    out = pd.DataFrame(rows)
    out.to_csv(TAB / "dial_audit_v17.tsv", sep="\t", index=False)

    if not out.empty:
        forest = px.scatter(out, x="DIA_AUC", y="cohort", color="ComBat_applied", facet_col="cluster", title="v17 DIAL forest proxy")
        forest.write_html(FIG / "dial_forest_plot.html", include_plotlyjs="cdn")
        pivot = out.pivot_table(index="cohort", columns=["cluster", "ComBat_applied"], values="identifiability")
        heat = px.imshow(pivot, aspect="auto", title="v17 identifiability heatmap")
        heat.write_html(FIG / "dial_identifiability_heatmap.html", include_plotlyjs="cdn")
        comp = px.box(out, x="ComBat_applied", y="identifiability", color="cluster", points="all", title="ComBat comparison of identifiability")
        comp.write_html(FIG / "dial_combat_comparison.html", include_plotlyjs="cdn")
    log_line(LOG, f"task4 done rows={len(out)}")


if __name__ == "__main__":
    main()
