#!/usr/bin/env python3
from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.express as px
import scanpy as sc
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
from sklearn.neighbors import KNeighborsClassifier
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression

from v17p2_common import FIG, LOG, PHASE1_TAB, SCRNA_H5AD, TAB, dm_binary_labels, ensure_overlap, identifiability_auc, load_dark_cohort, log_line, sparse_gene_score, top_marker_sets, write_json
from v17_common import align_dataset_expr_meta


DATASETS = ["GSE27155", "GSE126698", "GSE76039", "GSE213647"]


def build_classifiers():
    return {
        "LogReg": make_pipeline(StandardScaler(), LogisticRegression(max_iter=4000, random_state=42)),
        "RandomForest": RandomForestClassifier(n_estimators=400, random_state=42, min_samples_leaf=3),
        "SVM_RBF": make_pipeline(StandardScaler(), SVC(probability=True, kernel="rbf", gamma="scale", random_state=42)),
        "kNN_marker20": make_pipeline(StandardScaler(), KNeighborsClassifier(n_neighbors=9)),
    }


def main() -> None:
    log_line(LOG, "v17p2 layer4 start")
    tcga_expr, meta = load_dark_cohort()
    y = dm_binary_labels(meta)
    marker_sets = top_marker_sets()
    use_genes = sorted(set(marker_sets["DM1"][:20] + marker_sets["DM2"][:20]))
    tcga_train = tcga_expr[use_genes].copy()

    diag_rows = []
    auc_rows = []
    pred_rows = []
    for ds in DATASETS:
        ext_expr, ext_meta = align_dataset_expr_meta(ds, genes=use_genes, tumor_only=True)
        ext_meta = ext_meta[~ext_meta["driver_anchor"].isin(["BRAF", "RAS"])].copy()
        common = [s for s in ext_meta["sample_id"] if s in ext_expr.index]
        ext_meta = ext_meta.set_index("sample_id").loc[common].reset_index() if common else ext_meta.iloc[0:0].copy()
        ext_expr = ext_expr.loc[common] if common else ext_expr.iloc[0:0].copy()
        overlap = ensure_overlap(tcga_train, ext_expr, use_genes)
        failure = []
        if len(common) == 0:
            failure.append("no_dark_candidates")
        if len(overlap) < 10:
            failure.append("low_gene_overlap")
        diag_rows.append({
            "cohort": ds,
            "n_tumor_candidates": int(len(ext_meta)),
            "n_gene_overlap": int(len(overlap)),
            "failure_mode": ",".join(failure) if failure else "ok",
            "expr_min": float(np.nanmin(ext_expr[overlap].to_numpy(float))) if overlap and len(ext_expr) else np.nan,
            "expr_max": float(np.nanmax(ext_expr[overlap].to_numpy(float))) if overlap and len(ext_expr) else np.nan,
        })
        if len(ext_meta) < 5 or len(overlap) < 10:
            continue
        centroid_dm1 = tcga_train.loc[y == 0, overlap].mean(axis=0)
        centroid_dm2 = tcga_train.loc[y == 1, overlap].mean(axis=0)
        d1 = ((ext_expr[overlap] - centroid_dm1) ** 2).sum(axis=1)
        d2 = ((ext_expr[overlap] - centroid_dm2) ** 2).sum(axis=1)
        pseudo = (d2 < d1).astype(int).to_numpy()
        if len(np.unique(pseudo)) < 2:
            pseudo[0] = 1 - pseudo[0]
        for name, clf in build_classifiers().items():
            clf.fit(tcga_train[overlap], y)
            prob = clf.predict_proba(ext_expr[overlap])[:, 1]
            dia = max(roc_auc_score(pseudo, prob), 1.0 - roc_auc_score(pseudo, prob))
            auc_rows.append({"cohort": ds, "classifier": name, "AUC": float(roc_auc_score(pseudo, prob)), "DIA_AUC": float(dia), "n": int(len(ext_expr))})
            pred = pd.DataFrame({
                "cohort": ds,
                "classifier": name,
                "sample_id": ext_meta["sample_id"].values,
                "pseudo_cluster": np.where(pseudo == 1, "DM2", "DM1"),
                "prob_dm2": prob,
                "predicted_cluster": np.where(prob >= 0.5, "DM2", "DM1"),
                "confidence": np.abs(prob - 0.5) * 2.0,
            })
            pred_rows.append(pred)
    pd.DataFrame(diag_rows).to_csv(TAB / "external_alignment_diagnosis.tsv", sep="\t", index=False)
    auc_df = pd.DataFrame(auc_rows)
    auc_df.to_csv(TAB / "external_classifier_comparison.tsv", sep="\t", index=False)
    pred_df = pd.concat(pred_rows, ignore_index=True) if pred_rows else pd.DataFrame(columns=["cohort", "classifier", "sample_id", "predicted_cluster", "confidence"])
    pred_df.to_csv(TAB / "external_transfer_per_sample.tsv", sep="\t", index=False)

    if not auc_df.empty:
        best = auc_df.sort_values(["cohort", "DIA_AUC"], ascending=[True, False]).drop_duplicates("cohort")
        px.scatter(best, x="DIA_AUC", y="cohort", size="n", color="classifier", title="External 4-cohort transfer").write_html(FIG / "external_4cohort_forest.html", include_plotlyjs="cdn")
        pivot = auc_df.pivot(index="classifier", columns="cohort", values="DIA_AUC")
        px.imshow(pivot, text_auto=".2f", aspect="auto", title="Classifier by cohort DIA-AUC").write_html(FIG / "external_confusion_matrices.html", include_plotlyjs="cdn")

    scrna_rows = []
    if SCRNA_H5AD.exists():
        ad = sc.read_h5ad(SCRNA_H5AD)
        ad.var["gene_symbol"] = ad.var["gene_symbol"].astype(str)
        idx_dm1 = [i for i, g in enumerate(ad.var["gene_symbol"]) if g in marker_sets["DM1"][:20]]
        idx_dm2 = [i for i, g in enumerate(ad.var["gene_symbol"]) if g in marker_sets["DM2"][:20]]
        ad.obs["DM1_score"] = sparse_gene_score(ad.X, idx_dm1)
        ad.obs["DM2_score"] = sparse_gene_score(ad.X, idx_dm2)
        ad.obs["dominant_dm"] = np.where(ad.obs["DM2_score"] > ad.obs["DM1_score"], "DM2", "DM1")
        scrna = ad.obs.reset_index().rename(columns={"index": "cell_id"})
        scrna[["cell_id", "patient_id", "DM1_score", "DM2_score", "dominant_dm"]].to_csv(TAB / "scrna_cell_dm_score.tsv", sep="\t", index=False)
        use = ad[:, sorted(set(idx_dm1 + idx_dm2))].copy()
        sc.pp.scale(use, max_value=10)
        sc.tl.pca(use, n_comps=min(10, use.n_vars - 1, use.n_obs - 1), svd_solver="arpack")
        sc.pp.neighbors(use, n_neighbors=15)
        sc.tl.umap(use, random_state=42)
        um = pd.DataFrame(use.obsm["X_umap"], columns=["UMAP1", "UMAP2"], index=use.obs_names)
        plot_df = scrna.set_index("cell_id").join(um).reset_index()
        px.scatter(plot_df.sample(min(30000, len(plot_df)), random_state=42), x="UMAP1", y="UMAP2", color="dominant_dm", hover_data=["patient_id"], title="scRNA DM signature projection").write_html(
            FIG / "scrna_dm_signature_umap.html", include_plotlyjs="cdn"
        )
        per_pt = plot_df.groupby(["patient_id", "dominant_dm"]).size().reset_index(name="n")
        px.bar(per_pt, x="patient_id", y="n", color="dominant_dm", barmode="stack", title="Per-patient DM heterogeneity").write_html(
            FIG / "scrna_per_patient_heterogeneity.html", include_plotlyjs="cdn"
        )
    else:
        pd.DataFrame(columns=["cell_id", "patient_id", "DM1_score", "DM2_score", "dominant_dm"]).to_csv(TAB / "scrna_cell_dm_score.tsv", sep="\t", index=False)

    summary = {
        "cohorts_tested": int(len(auc_df["cohort"].unique())) if not auc_df.empty else 0,
        "cohorts_dia_auc_gt_0_85": int(auc_df.groupby("cohort")["DIA_AUC"].max().ge(0.85).sum()) if not auc_df.empty else 0,
        "scrna_cells_scored": int(pd.read_csv(TAB / "scrna_cell_dm_score.tsv", sep="\t").shape[0]),
    }
    write_json(TAB / "external_summary.json", summary)
    log_line(LOG, f"v17p2 layer4 done {summary}")


if __name__ == "__main__":
    main()
