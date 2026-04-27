"""
v17p35 AMP-4 — 8-gene RAI decision tool (npj 결정타).

Train LogReg on TCGA-THCA (RAI score binary at median), 5-fold CV, then validate on GSE76039 (PDTC/ATC).
Compare vs BRAF-only baseline (NRI).
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score, roc_curve

PROJECT = Path("/opt/thyroid-dash/project")
TAB = PROJECT / "results/v17p35/tables"
FIG = PROJECT / "results/v17p35/figs"
np.random.seed(42)

# 8-gene RAI panel
GENE_SYMBOLS = ["SLC5A5", "TPO", "TG", "TSHR", "PAX8", "NKX2-1", "FOXE1", "DIO1"]


def load_expr(path):
    df = pd.read_csv(path, sep="\t", index_col=0)
    # try to find genes
    return df


def main():
    # Use A2 full-cohort DM scores (513 TCGA samples, includes BRAF/RAS labelled)
    a2 = pd.read_csv(PROJECT / "results/v17p3/tables/A2_dm_score_full_cohort.tsv", sep="\t")
    a2 = a2[(a2["dm_like"].isin(["DM1_like", "DM2_like"])) &
            (a2["dataset"] == "TCGA-THCA")].copy()
    rai = a2.rename(columns={"sample_id": "sample_id"})[["sample_id", "dm_like", "rai_score_recalc", "driver_anchor"]]
    rai["rai_score"] = rai["rai_score_recalc"]
    print(f"A2 full-cohort DM samples (TCGA): {len(rai)}")
    print(f"  DM1: {(rai['dm_like']=='DM1_like').sum()}, DM2: {(rai['dm_like']=='DM2_like').sum()}")

    # TCGA expression
    tcga = load_expr(PROJECT / "data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_log2.tsv")
    print(f"TCGA expr: {tcga.shape}")
    found = [g for g in GENE_SYMBOLS if g in tcga.index]
    missing = [g for g in GENE_SYMBOLS if g not in tcga.index]
    print(f"Genes found in TCGA: {found}, missing: {missing}")
    # try alternates
    if "NKX2-1" in missing and "NKX2_1" in tcga.index:
        tcga.index = tcga.index.str.replace("NKX2_1", "NKX2-1")
        found.append("NKX2-1")

    # Match by short barcode
    rai["short"] = rai["sample_id"].str[:15]
    tcga_cols_short = {c[:15]: c for c in tcga.columns if c.startswith("TCGA")}
    rai["mapped"] = rai["short"].map(tcga_cols_short)
    rai = rai.dropna(subset=["mapped"])
    print(f"  Mapped to TCGA expr: {len(rai)}")

    common = rai["mapped"].tolist()

    # Build X (n_samples × 8 genes), y = DM1 (1) vs DM2 (0)
    X = tcga.loc[found, common].T.values
    y = (rai["dm_like"].values == "DM1_like").astype(int)  # 1 = DM1-like (MAPK active, dediff)
    median_rai = float(rai["rai_score"].median())
    print(f"  X: {X.shape}, y mean (DM1 fraction): {y.mean():.3f}")

    # 5-fold CV
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    models = {}

    # LogReg
    lr = LogisticRegression(max_iter=2000, random_state=42, C=1.0)
    lr_pred = cross_val_predict(lr, X, y, cv=skf, method="predict_proba")[:, 1]
    auc_lr = roc_auc_score(y, lr_pred)
    models["LogReg"] = {"model": lr, "cv_pred": lr_pred, "auc": auc_lr}
    print(f"  LogReg 5-fold CV AUC = {auc_lr:.3f}")

    # RandomForest
    rf = RandomForestClassifier(n_estimators=300, random_state=42, n_jobs=2)
    rf_pred = cross_val_predict(rf, X, y, cv=skf, method="predict_proba")[:, 1]
    auc_rf = roc_auc_score(y, rf_pred)
    models["RandomForest"] = {"model": rf, "cv_pred": rf_pred, "auc": auc_rf}
    print(f"  RandomForest 5-fold CV AUC = {auc_rf:.3f}")

    # Fit final models on all training data
    lr.fit(X, y)
    rf.fit(X, y)

    # Coefficients
    coef_df = pd.DataFrame({
        "gene": found,
        "logreg_coef": lr.coef_[0],
        "rf_importance": rf.feature_importances_,
    }).sort_values("rf_importance", ascending=False)
    coef_df.to_csv(TAB / "AMP4_8gene_model_coefficients.tsv", sep="\t", index=False)
    print(f"\n{coef_df.to_string(index=False)}")

    # ROC curve points
    fpr_lr, tpr_lr, _ = roc_curve(y, lr_pred)
    fpr_rf, tpr_rf, _ = roc_curve(y, rf_pred)
    roc_df = pd.DataFrame({
        "fpr": np.concatenate([fpr_lr, fpr_rf]),
        "tpr": np.concatenate([tpr_lr, tpr_rf]),
        "model": ["LogReg"] * len(fpr_lr) + ["RandomForest"] * len(fpr_rf),
    })
    roc_df.to_csv(TAB / "AMP4_roc_data.tsv", sep="\t", index=False)

    # ---- BRAF-only baseline (from rai DataFrame which has driver_anchor) ----
    sm_thca = None
    sm_path = PROJECT / "results/v17/tables/sample_master_v17_full.tsv"
    if sm_path.exists():
        sm = pd.read_csv(sm_path, sep="\t")
        sm_thca = sm.set_index("sample_id")

    braf_binary = (rai["driver_anchor"].values == "BRAF").astype(int)
    auc_braf = None
    nri_proxy = None
    if braf_binary.sum() > 5 and braf_binary.sum() < len(braf_binary) - 5:
        braf_lr = LogisticRegression(max_iter=1000, random_state=42)
        braf_pred = cross_val_predict(braf_lr, braf_binary.reshape(-1, 1), y, cv=skf, method="predict_proba")[:, 1]
        auc_braf = roc_auc_score(y, braf_pred)
        nri_proxy = auc_lr - auc_braf
        print(f"\n  BRAF-only baseline AUC = {auc_braf:.3f}")
        print(f"  ΔAUC (8-gene LogReg − BRAF-only) = {nri_proxy:+.3f}")

    # ---- PDTC validation on GSE76039 ----
    pdtc_auc = None
    if sm_thca is not None:
        try:
            gse = load_expr(PROJECT / "data_processed/microarray/GSE76039_microarray_expression_log2.tsv")
            gse_found = [g for g in found if g in gse.index]
            print(f"\n  GSE76039: {gse.shape}, 8-gene present: {len(gse_found)}")

            pdtc_samples = sm_thca[sm_thca["dataset"] == "GSE76039"].copy()
            pdtc_y_full = (pdtc_samples["histology_subtype"] == "PDTC").astype(int)
            common_gse = list(set(pdtc_samples.index) & set(gse.columns))
            if len(common_gse) > 10 and len(gse_found) >= 4:
                X_gse = gse.loc[gse_found, common_gse].T.values
                y_gse = pdtc_y_full.loc[common_gse].values
                X_tcga_sub = tcga.loc[gse_found, common].T.values
                lr2 = LogisticRegression(max_iter=2000, random_state=42)
                lr2.fit(X_tcga_sub, y)
                gse_pred = lr2.predict_proba(X_gse)[:, 1]
                if y_gse.std() > 0:
                    pdtc_auc = float(roc_auc_score(y_gse, gse_pred))
                    print(f"  GSE76039 PDTC vs ATC AUC (8-gene model) = {pdtc_auc:.3f} (n={len(y_gse)})")
        except Exception as e:
            print(f"  PDTC validation fail: {e}")

    # Save summary
    summary = {
        "n_train_samples": int(len(y)),
        "n_genes_used": int(len(found)),
        "genes_used": found,
        "rai_threshold": median_rai,
        "logreg_cv_auc": round(auc_lr, 3),
        "rf_cv_auc": round(auc_rf, 3),
        "braf_only_auc": round(auc_braf, 3) if auc_braf is not None else None,
        "auc_delta_8gene_vs_braf": round(nri_proxy, 3) if nri_proxy is not None else None,
        "pdtc_validation_auc": round(pdtc_auc, 3) if pdtc_auc is not None else None,
    }
    (TAB / "AMP4_summary.json").write_text(json.dumps(summary, indent=2))

    # Also write CV performance table
    cv_perf = pd.DataFrame([
        {"model": "LogReg_8gene", "cv_auc": auc_lr, "n_samples": len(y)},
        {"model": "RandomForest_8gene", "cv_auc": auc_rf, "n_samples": len(y)},
        {"model": "LogReg_BRAF_only", "cv_auc": auc_braf, "n_samples": len(y)},
        {"model": "LogReg_8gene_PDTC_extval", "cv_auc": pdtc_auc, "n_samples": "GSE76039"},
    ])
    cv_perf.to_csv(TAB / "AMP4_cv_performance.tsv", sep="\t", index=False)

    print("\n=== AMP-4 SUMMARY ===")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
