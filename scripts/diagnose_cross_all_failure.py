#!/usr/bin/env python3
"""Diagnose why v0 F_CROSS_all underperforms C_counterfactual."""

from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedShuffleSplit

from cross_neo_v1_common import (
    V0,
    V1,
    build_base_feature_table,
    ensure_v1_dirs,
    fill_by_train_median,
    get_fold_ids,
    load_folds,
    load_master,
    metrics,
)


def group_columns(base: pd.DataFrame) -> dict[str, list[str]]:
    cf = [c for c in base.columns if c.startswith("cf_")]
    qk = [c for c in ["GP_quantum", "VQC", "W7A_QK_only", "W7A_full", "tcr_motif_score"] if c in base.columns]
    geom = [
        c
        for c in base.columns
        if c
        not in {
            "sample_id",
            "label",
            "peptide_mut",
            "hla",
            "hla_supertype",
            "study",
            "near_peptide_cluster",
            "peptide_len",
            *cf,
            *qk,
        }
    ]
    retr = [
        "near_peptide_similarity_train",
        "near_positive_peptide_similarity_train",
        "same_hla_near_peptide_similarity_train",
        "wt_self_similarity",
        "same_hla_train_density",
        "same_hla_positive_rate_train",
        "train_fold_prevalence",
        "exact_peptide_hla_hit_train",
        "exact_peptide_hit_train",
    ]
    return {"C_counterfactual": cf, "QK_fixed": qk, "structure_geometry": geom, "retrieval": retr}


def fit_rf() -> RandomForestClassifier:
    return RandomForestClassifier(
        n_estimators=160,
        max_depth=3,
        min_samples_leaf=5,
        class_weight="balanced",
        random_state=20260509,
        n_jobs=-1,
    )


def main() -> None:
    ensure_v1_dirs()
    master = load_master()
    folds = load_folds()
    base = build_base_feature_table(master)
    groups = group_columns(base)
    retr = pd.read_csv(V0 / "retrieval_features_by_fold.tsv", sep="\t")
    v0m = pd.read_csv(V0 / "metrics_by_split.tsv", sep="\t")
    diagnosis = []
    for name, cols in groups.items():
        present = [c for c in cols if c in base.columns or c in retr.columns]
        miss = base[[c for c in present if c in base.columns]].isna().mean().mean() if present else np.nan
        diagnosis.append({"hypothesis": "feature_group_dimensionality", "feature_group": name, "n_features": len(present), "mean_missing_rate": miss})
    for split in sorted(v0m["split_name"].unique()):
        c = v0m[(v0m["split_name"] == split) & (v0m["feature_group"] == "C_counterfactual") & (v0m["model"] == "rf_secondary")]
        f = v0m[(v0m["split_name"] == split) & (v0m["feature_group"] == "F_CROSS_all") & (v0m["model"] == "rf_secondary")]
        if len(c) and len(f):
            diagnosis.append(
                {
                    "hypothesis": "branch_noise_overwhelms_counterfactual",
                    "split_name": split,
                    "C_AUPRC": float(c.iloc[0]["AUPRC"]),
                    "F_CROSS_all_AUPRC": float(f.iloc[0]["AUPRC"]),
                    "delta_all_minus_C": float(f.iloc[0]["AUPRC"] - c.iloc[0]["AUPRC"]),
                    "C_top10": float(c.iloc[0]["top10_precision"]),
                    "F_top10": float(f.iloc[0]["top10_precision"]),
                }
            )
    # Variance and shortcut diagnostics.
    num_cols = groups["C_counterfactual"] + groups["QK_fixed"] + groups["structure_geometry"]
    for by in ["study", "hla_supertype", "label"]:
        for grp, sub in base.groupby(by, dropna=False):
            if len(sub) < 3:
                continue
            diagnosis.append(
                {
                    "hypothesis": f"feature_variance_by_{by}",
                    "stratum": str(grp),
                    "n": len(sub),
                    "mean_feature_variance": float(sub[num_cols].var(numeric_only=True).mean()),
                    "label_rate": float(sub["label"].mean()),
                }
            )
    pd.DataFrame(diagnosis).to_csv(V1 / "cross_all_failure_diagnosis.tsv", sep="\t", index=False)

    feature_sets = {
        "C only": groups["C_counterfactual"],
        "C + QK": groups["C_counterfactual"] + groups["QK_fixed"],
        "C + structure": groups["C_counterfactual"] + groups["structure_geometry"],
        "C + retrieval": groups["C_counterfactual"] + groups["retrieval"],
        "C + structure + QK": groups["C_counterfactual"] + groups["structure_geometry"] + groups["QK_fixed"],
        "C + retrieval + QK": groups["C_counterfactual"] + groups["retrieval"] + groups["QK_fixed"],
        "all": groups["C_counterfactual"] + groups["structure_geometry"] + groups["QK_fixed"] + groups["retrieval"],
    }
    pred_rows = []
    for (split, fold), fdf in folds.groupby(["split_name", "fold_id"], sort=False):
        train_ids, test_ids = get_fold_ids(master, folds, split, fold)
        rfold = retr[(retr["split_name"] == split) & (retr["fold_id"] == fold)]
        feat = base.merge(rfold[["sample_id", *groups["retrieval"]]], on="sample_id", how="left")
        for c in groups["retrieval"]:
            feat[c] = pd.to_numeric(feat[c], errors="coerce").fillna(0.0)
        train = feat[feat["sample_id"].isin(train_ids)].copy()
        test = feat[feat["sample_id"].isin(test_ids)].copy()
        if len(train) < 10 or len(test) < 2 or train["label"].nunique() < 2 or test["label"].nunique() < 2:
            continue
        for name, cols in feature_sets.items():
            cols = [c for c in cols if c in feat.columns]
            xtr, xte = fill_by_train_median(train, test, cols)
            ytr = train["label"].to_numpy(int)
            yte = test["label"].to_numpy(int)
            clf = fit_rf()
            clf.fit(xtr, ytr)
            pred = clf.predict_proba(xte)[:, 1]
            for sid, y, p in zip(test["sample_id"], yte, pred):
                pred_rows.append({"split_name": split, "fold_id": fold, "feature_set": name, "sample_id": sid, "label": int(y), "score": float(p)})
    pred = pd.DataFrame(pred_rows)
    metric_rows = []
    for (split, fs), sub in pred.groupby(["split_name", "feature_set"]):
        metric_rows.append({"split_name": split, "feature_set": fs, **metrics(sub["label"].to_numpy(), sub["score"].to_numpy())})
    ab = pd.DataFrame(metric_rows).sort_values(["split_name", "AUPRC"], ascending=[True, False])
    ab.to_csv(V1 / "drop_group_ablation.tsv", sep="\t", index=False)

    # Fold-only group permutation on inner validation, never outer test.
    imp_rows = []
    for (split, fold), fdf in list(folds.groupby(["split_name", "fold_id"], sort=False))[:40]:
        train_ids, _ = get_fold_ids(master, folds, split, fold)
        tr0 = base[base["sample_id"].isin(train_ids)].copy()
        if len(tr0) < 25 or tr0["label"].nunique() < 2:
            continue
        sss = StratifiedShuffleSplit(n_splits=1, test_size=0.28, random_state=20260509)
        idx_tr, idx_val = next(sss.split(np.zeros(len(tr0)), tr0["label"].to_numpy(int)))
        tr, val = tr0.iloc[idx_tr].copy(), tr0.iloc[idx_val].copy()
        all_cols = groups["C_counterfactual"] + groups["QK_fixed"] + groups["structure_geometry"]
        xtr, xval = fill_by_train_median(tr, val, all_cols)
        ytr, yval = tr["label"].to_numpy(int), val["label"].to_numpy(int)
        clf = fit_rf()
        clf.fit(xtr, ytr)
        base_score = clf.predict_proba(xval)[:, 1]
        base_ap = metrics(yval, base_score)["AUPRC"]
        val_df = val[all_cols].copy()
        for gname, cols in {k: v for k, v in groups.items() if k != "retrieval"}.items():
            perm = val_df.copy()
            rng = np.random.default_rng(20260509 + len(imp_rows))
            for c in cols:
                if c in perm.columns:
                    perm[c] = rng.permutation(perm[c].to_numpy())
            _, xperm = fill_by_train_median(tr, perm.assign(sample_id=val["sample_id"].values, label=yval), all_cols)
            p = clf.predict_proba(xperm)[:, 1]
            imp_rows.append({"split_name": split, "fold_id": fold, "feature_group": gname, "inner_val_AUPRC": base_ap, "permuted_AUPRC": metrics(yval, p)["AUPRC"], "delta": base_ap - metrics(yval, p)["AUPRC"]})
    imp = pd.DataFrame(imp_rows)
    imp.to_csv(V1 / "feature_group_importance.tsv", sep="\t", index=False)
    lines = [
        "# CROSS-all Failure Diagnosis",
        "",
        "Primary observation: naive concatenation does not dominate counterfactual features, consistent with small-n noisy-branch overfitting.",
        "",
        "## Drop-Group Ablation Top Rows",
        ab.sort_values("AUPRC", ascending=False).head(12).to_markdown(index=False),
        "",
        "## Inner-Fold Permutation Importance",
        imp.groupby("feature_group")["delta"].agg(["count", "mean", "median"]).reset_index().to_markdown(index=False) if len(imp) else "Not available.",
        "",
        "Interpretation: use late/gated fusion rather than raw feature concatenation.",
    ]
    (V1 / "cross_all_failure_report.md").write_text("\n".join(lines) + "\n")
    print(f"[v1-crossall] ablation_rows={len(ab)} importance_rows={len(imp)}")


if __name__ == "__main__":
    main()
