#!/usr/bin/env python3
from __future__ import annotations

import json
import re

import numpy as np
import pandas as pd
import plotly.express as px
import requests
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.neighbors import KNeighborsClassifier
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

from v17_common import align_dataset_expr_meta, load_sample_master
from v17p3_common import FIG, LOG, SEED, TAB, TIERA67, V17TAB, dia_auc, dm_labels, fit_dm_classifier, load_dark_tiera, log_line, preprocess_expr, write_json


def query_mygene(ids: list[str]) -> pd.DataFrame:
    rows = []
    chunk = 800
    for i in range(0, len(ids), chunk):
        sub = ids[i:i + chunk]
        r = requests.post(
            "https://mygene.info/v3/query",
            data={"q": ",".join(sub), "scopes": "ensembl.gene", "fields": "symbol", "species": "human"},
            timeout=120,
        )
        r.raise_for_status()
        rows.extend(r.json())
    out = []
    for rec in rows:
        q = str(rec.get("query", ""))
        sym = rec.get("symbol")
        if q and sym:
            out.append({"ensembl": q, "symbol": str(sym).upper()})
    return pd.DataFrame(out).drop_duplicates()


def recover_gse213647() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    raw = pd.read_csv("/data/thca/data_processed/bulk_rnaseq_v3/GSE213647_v3_log2.tsv", sep="\t")
    first = raw.iloc[:500, 0].astype(str)
    ens_like = first.str.match(r"^ENSG").mean()
    diag = pd.DataFrame([{
        "dataset": "GSE213647",
        "first_column": raw.columns[0],
        "ensembl_like_fraction_top500": float(ens_like),
        "failure_mode": "column named gene_symbol but values are Ensembl IDs" if ens_like > 0.9 else "other",
        "raw_rows": int(raw.shape[0]),
        "raw_samples": int(raw.shape[1] - 1),
    }])
    ids = raw.iloc[:, 0].astype(str).str.replace(r"\\..*$", "", regex=True).tolist()
    mapping = query_mygene(sorted(set(ids)))
    raw["ensembl"] = raw.iloc[:, 0].astype(str).str.replace(r"\\..*$", "", regex=True)
    rec = raw.merge(mapping, on="ensembl", how="left")
    rec = rec.dropna(subset=["symbol"]).drop_duplicates("symbol").set_index("symbol")
    expr = rec.drop(columns=[raw.columns[0], "ensembl"]).T
    expr.index.name = "sample_id"
    meta = load_sample_master()
    meta = meta[(meta["dataset"] == "GSE213647") & (meta["normal_vs_tumor"] != "normal")].copy()
    common = [s for s in meta["sample_id"] if s in expr.index]
    meta = meta.set_index("sample_id").loc[common].reset_index()
    expr = expr.loc[common].apply(pd.to_numeric, errors="coerce")
    return diag, mapping, expr


def main() -> None:
    log_line(LOG, "v17p3 F1 start")
    tcga_expr, dark_meta = load_dark_tiera()
    marker_df = pd.read_csv(V17TAB / "dark_matter_cluster_markers.tsv", sep="\t")
    marker_genes = marker_df.groupby("cluster")["gene"].head(20).astype(str).str.upper().tolist()
    train_genes = sorted(set([g.upper() for g in TIERA67] + marker_genes))
    tcga_expr.columns = tcga_expr.columns.astype(str).str.upper()

    diag, mapping, g213_expr = recover_gse213647()
    mapping.to_csv(TAB / "F1_gene_recovery_mapping.tsv", sep="\t", index=False)
    diag.to_csv(TAB / "F1_gene_id_diagnosis.tsv", sep="\t", index=False)

    rows = []
    pred_rows = []
    datasets = {
        "GSE27155": align_dataset_expr_meta("GSE27155", genes=train_genes, tumor_only=True)[0],
        "GSE76039": align_dataset_expr_meta("GSE76039", genes=train_genes, tumor_only=True)[0],
        "GSE126698": align_dataset_expr_meta("GSE126698", genes=train_genes, tumor_only=True)[0],
        "GSE213647": g213_expr,
    }
    meta_map = {}
    for ds in ["GSE27155", "GSE76039", "GSE126698"]:
        _, meta = align_dataset_expr_meta(ds, genes=train_genes, tumor_only=True)
        meta_map[ds] = meta
    meta_map["GSE213647"] = load_sample_master().query("dataset == 'GSE213647' and normal_vs_tumor != 'normal'").copy()

    classifiers = {
        "LogReg": make_pipeline(StandardScaler(), LogisticRegression(max_iter=4000, random_state=SEED)),
        "RandomForest": RandomForestClassifier(n_estimators=400, random_state=SEED, min_samples_leaf=3),
        "SVM_RBF": make_pipeline(StandardScaler(), SVC(probability=True, kernel="rbf", gamma="scale", random_state=SEED)),
        "kNN": make_pipeline(StandardScaler(), KNeighborsClassifier(n_neighbors=9)),
    }
    y = dm_labels(dark_meta)
    for ds, expr in datasets.items():
        if expr is None or expr.empty:
            continue
        expr = expr.copy()
        expr.columns = expr.columns.astype(str).str.upper()
        use = [g for g in train_genes if g in tcga_expr.columns and g in expr.columns]
        meta = meta_map[ds].copy()
        if "sample_id" in meta.columns:
            common = [s for s in meta["sample_id"] if s in expr.index]
            meta = meta.set_index("sample_id").loc[common].reset_index()
            expr = expr.loc[common]
        if ds in {"GSE27155", "GSE76039", "GSE213647"}:
            meta = meta[~meta["driver_anchor"].isin(["BRAF", "RAS"])].copy()
            expr = expr.loc[meta["sample_id"]]
        if len(use) == 0 or len(expr) == 0:
            continue
        pos = tcga_expr.loc[y == 0, use].mean(axis=0)
        neg = tcga_expr.loc[y == 1, use].mean(axis=0)
        d_dm1 = ((expr[use] - pos) ** 2).sum(axis=1)
        d_dm2 = ((expr[use] - neg) ** 2).sum(axis=1)
        pseudo = (d_dm2 < d_dm1).astype(int).to_numpy()
        if len(np.unique(pseudo)) < 2 and len(pseudo) > 1:
            pseudo[0] = 1 - pseudo[0]
        for name, clf in classifiers.items():
            clf.fit(tcga_expr[use], y)
            prob = clf.predict_proba(expr[use])[:, 1]
            if len(np.unique(pseudo)) < 2:
                auc = np.nan
                dia = np.nan
            else:
                auc = float(roc_auc_score(pseudo, prob))
                dia = float(dia_auc(pseudo, prob))
            rows.append({"cohort": ds, "classifier": name, "n": int(len(expr)), "gene_overlap": int(len(use)), "AUC": auc, "DIA_AUC": dia})
            pred_rows.append(pd.DataFrame({
                "cohort": ds,
                "classifier": name,
                "sample_id": expr.index,
                "prob_dm2": prob,
                "predicted_cluster": np.where(prob >= 0.5, "DM2", "DM1"),
            }))
    out = pd.DataFrame(rows)
    out.to_csv(TAB / "F1_external_5cohort_recovery.tsv", sep="\t", index=False)
    pred = pd.concat(pred_rows, ignore_index=True) if pred_rows else pd.DataFrame()
    pred.to_csv(TAB / "F1_external_predictions.tsv", sep="\t", index=False)

    if not out.empty:
        best = out.sort_values(["cohort", "DIA_AUC"], ascending=[True, False]).drop_duplicates("cohort")
        px.scatter(best, x="DIA_AUC", y="cohort", color="classifier", size="gene_overlap", title="Phase 3 5-cohort external recovery").write_html(
            FIG / "F1_5cohort_forest.html", include_plotlyjs="cdn"
        )
        heat = out.pivot(index="classifier", columns="cohort", values="DIA_AUC")
        px.imshow(heat, text_auto=".2f", aspect="auto", title="Phase 3 DIA-AUC heatmap").write_html(
            FIG / "F1_5cohort_dia_auc_heatmap.html", include_plotlyjs="cdn"
        )
        g213 = pred[(pred["cohort"] == "GSE213647") & (pred["classifier"] == "LogReg")]
        if not g213.empty:
            px.histogram(g213, x="prob_dm2", color="predicted_cluster", nbins=40, title="Recovered GSE213647 DM distribution").write_html(
                FIG / "F1_recovered_gse213647_dm_distribution.html", include_plotlyjs="cdn"
            )
    summary = {
        "gse213647_gene_overlap_tiera": int(sum(g.upper() in g213_expr.columns for g in TIERA67)),
        "robust_cohort_count_dia_auc_gt_0_85": int(out.groupby("cohort")["DIA_AUC"].max().ge(0.85).sum()) if not out.empty else 0,
    }
    write_json(TAB / "F1_summary.json", summary)
    log_line(LOG, f"v17p3 F1 done {summary}")


if __name__ == "__main__":
    main()

