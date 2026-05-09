#!/usr/bin/env python3
"""Fold-safe gated mixture-of-experts for CROSS-Neo v1."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

from cross_neo_v1_common import (
    V1,
    build_base_feature_table,
    ensure_v1_dirs,
    fill_by_train_median,
    fold_context_features,
    get_fold_ids,
    load_folds,
    load_master,
    metrics,
    normalize_weights,
    qkernel,
    weighted_score,
)


EXPERTS = ["expert_counterfactual", "expert_qk_no_anchor", "expert_qk_quantum_only", "expert_structure", "expert_retrieval"]


def group_cols(base: pd.DataFrame) -> dict[str, list[str]]:
    cf = [c for c in base.columns if c.startswith("cf_")]
    q = [c for c in ["GP_quantum", "VQC", "W7A_QK_only", "W7A_full"] if c in base.columns]
    geom = [
        c
        for c in [
            "mean_pLDDT_peptide",
            "min_pLDDT_peptide",
            "mean_pLDDT_HLA",
            "anchor_pLDDT",
            "interface_contacts_8A",
            "interface_contacts_10A",
            "n_buried_residues_8A",
            "radius_of_gyration_peptide",
            "peptide_helicity_proxy",
            "peptide_bulge_proxy",
            "contact_density",
            "structure_missing",
            "structure_low_confidence",
        ]
        if c in base.columns
    ]
    return {"cf": cf, "q": q, "geom": geom}


def rf(max_depth: int = 3) -> RandomForestClassifier:
    return RandomForestClassifier(
        n_estimators=180,
        max_depth=max_depth,
        min_samples_leaf=5,
        class_weight="balanced",
        random_state=20260509,
        n_jobs=-1,
    )


def predict_rf(train: pd.DataFrame, test: pd.DataFrame, cols: list[str]) -> np.ndarray:
    if train["label"].nunique() < 2 or len(cols) == 0:
        return np.full(len(test), train["label"].mean() if len(train) else 0.5)
    xtr, xte = fill_by_train_median(train, test, cols)
    clf = rf()
    clf.fit(xtr, train["label"].to_numpy(int))
    return clf.predict_proba(xte)[:, 1]


def predict_qk(train: pd.DataFrame, test: pd.DataFrame, cols: list[str]) -> np.ndarray:
    if train["label"].nunique() < 2 or len(cols) == 0:
        return np.full(len(test), train["label"].mean() if len(train) else 0.5)
    xtr, xte = fill_by_train_median(train, test, cols)
    scaler = StandardScaler()
    xtr = scaler.fit_transform(xtr)
    xte = scaler.transform(xte)
    clf = SVC(kernel="precomputed", C=1.0, class_weight="balanced")
    clf.fit(qkernel(xtr, xtr), train["label"].to_numpy(int))
    raw = clf.decision_function(qkernel(xte, xtr))
    return 1.0 / (1.0 + np.exp(-np.clip(raw, -30, 30)))


def retrieval_score(ctx: pd.DataFrame) -> np.ndarray:
    s = (
        0.35 * ctx["ctx_same_hla_positive_rate_train"].fillna(ctx["ctx_train_prevalence"]).to_numpy(float)
        + 0.25 * ctx["ctx_near_peptide_similarity"].fillna(0).to_numpy(float)
        + 0.20 * ctx["ctx_exact_or_peptide_hit"].fillna(0).to_numpy(float)
        + 0.10 * ctx["ctx_near_hit"].fillna(0).to_numpy(float)
        + 0.10 * (1 - ctx["ctx_rare_hla"].fillna(1).to_numpy(float))
    )
    return np.clip(s, 0, 1)


def expert_scores(fit_df: pd.DataFrame, pred_df: pd.DataFrame, ref_df: pd.DataFrame, cols: dict[str, list[str]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    ctx = fold_context_features(pred_df, ref_df)
    out = pred_df[["sample_id", "label"]].copy()
    out["expert_counterfactual"] = predict_rf(fit_df, pred_df, cols["cf"])
    out["expert_structure"] = predict_rf(fit_df, pred_df, cols["geom"])
    out["expert_qk_quantum_only"] = predict_qk(fit_df, pred_df, cols["q"])
    out["expert_qk_no_anchor"] = predict_qk(fit_df, pred_df, cols["q"] + cols["geom"][:6])
    out["expert_retrieval"] = retrieval_score(ctx)
    return out, ctx


def inner_oof(train: pd.DataFrame, cols: dict[str, list[str]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    y = train["label"].to_numpy(int)
    n_splits = min(3, np.bincount(y).min() if len(np.unique(y)) == 2 else 1)
    if n_splits < 2:
        exp, ctx = expert_scores(train, train, train, cols)
        return exp, ctx
    rows, ctx_rows = [], []
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=20260509)
    for tr_idx, va_idx in skf.split(np.zeros(len(train)), y):
        fit = train.iloc[tr_idx].copy()
        val = train.iloc[va_idx].copy()
        exp, ctx = expert_scores(fit, val, fit, cols)
        rows.append(exp)
        ctx_rows.append(ctx)
    return pd.concat(rows, ignore_index=True), pd.concat(ctx_rows, ignore_index=True)


def rule_gated_score(exp: pd.DataFrame, ctx: pd.DataFrame) -> tuple[np.ndarray, list[dict[str, float]]]:
    scores = []
    weights = []
    e = exp.set_index("sample_id")
    c = ctx.set_index("sample_id")
    for sid in exp["sample_id"]:
        row = c.loc[sid]
        w = {
            "expert_counterfactual": 0.45,
            "expert_qk_no_anchor": 0.35,
            "expert_qk_quantum_only": 0.10,
            "expert_structure": 0.05,
            "expert_retrieval": 0.05,
        }
        if row["ctx_exact_or_peptide_hit"] or row["ctx_near_hit"]:
            w = {"expert_counterfactual": 0.35, "expert_qk_no_anchor": 0.25, "expert_qk_quantum_only": 0.10, "expert_structure": 0.05, "expert_retrieval": 0.25}
        if row["structure_missing"] == 0 and row["structure_low_confidence"] == 0:
            disagreement = abs(e.loc[sid, "expert_counterfactual"] - e.loc[sid, "expert_qk_no_anchor"])
            if disagreement >= 0.25:
                w["expert_structure"] += 0.10
                w["expert_counterfactual"] -= 0.05
                w["expert_qk_no_anchor"] -= 0.05
        if row["ctx_ood_score"] >= 5:
            w["expert_retrieval"] *= 0.25
            w["expert_structure"] *= 0.5
            w["expert_counterfactual"] += 0.10
            w["expert_qk_no_anchor"] += 0.05
        w = normalize_weights(w)
        weights.append(w)
        scores.append(float(sum(w[k] * e.loc[sid, k] for k in EXPERTS)))
    return np.asarray(scores), weights


def choose_nested_weights(inner_exp: pd.DataFrame) -> tuple[dict[str, float], pd.DataFrame]:
    candidates = [
        {"name": "C_QKno_equal", "expert_counterfactual": 0.5, "expert_qk_no_anchor": 0.5, "expert_qk_quantum_only": 0, "expert_structure": 0, "expert_retrieval": 0},
        {"name": "C_QKquant_equal", "expert_counterfactual": 0.5, "expert_qk_no_anchor": 0, "expert_qk_quantum_only": 0.5, "expert_structure": 0, "expert_retrieval": 0},
        {"name": "C_heavy_QKno", "expert_counterfactual": 0.65, "expert_qk_no_anchor": 0.35, "expert_qk_quantum_only": 0, "expert_structure": 0, "expert_retrieval": 0},
        {"name": "QKno_heavy", "expert_counterfactual": 0.35, "expert_qk_no_anchor": 0.65, "expert_qk_quantum_only": 0, "expert_structure": 0, "expert_retrieval": 0},
        {"name": "C_QKno_structure", "expert_counterfactual": 0.45, "expert_qk_no_anchor": 0.4, "expert_qk_quantum_only": 0, "expert_structure": 0.15, "expert_retrieval": 0},
        {"name": "all_lowdim_equalish", "expert_counterfactual": 0.35, "expert_qk_no_anchor": 0.3, "expert_qk_quantum_only": 0.15, "expert_structure": 0.1, "expert_retrieval": 0.1},
    ]
    rows = []
    for cand in candidates:
        name = cand.pop("name")
        w = normalize_weights(cand)
        s = weighted_score(inner_exp, w)
        mm = metrics(inner_exp["label"].to_numpy(), s)
        rows.append({"candidate": name, "weights": json.dumps(w), **mm})
        cand["name"] = name
    res = pd.DataFrame(rows).sort_values(["AUPRC", "top10_precision"], ascending=False)
    best = json.loads(res.iloc[0]["weights"])
    return best, res


def stacked_meta(inner_exp: pd.DataFrame, inner_ctx: pd.DataFrame, outer_exp: pd.DataFrame, outer_ctx: pd.DataFrame) -> np.ndarray:
    ctx_cols = ["ctx_peptide_len", "ctx_unseen_hla", "ctx_rare_hla", "ctx_near_peptide_similarity", "ctx_clean_no_reference", "structure_missing", "structure_low_confidence", "ctx_ood_score"]
    tr = inner_exp.merge(inner_ctx[["sample_id", *ctx_cols]], on="sample_id", how="left")
    te = outer_exp.merge(outer_ctx[["sample_id", *ctx_cols]], on="sample_id", how="left")
    for df in (tr, te):
        df["expert_margin"] = np.abs(df[EXPERTS].to_numpy(float) - 0.5).mean(axis=1)
        df["expert_range"] = df[EXPERTS].max(axis=1) - df[EXPERTS].min(axis=1)
    cols = EXPERTS + ctx_cols + ["expert_margin", "expert_range"]
    xtr, xte = fill_by_train_median(tr, te, cols)
    if tr["label"].nunique() < 2:
        return np.full(len(te), tr["label"].mean())
    clf = LogisticRegression(max_iter=2000, class_weight="balanced", C=0.5)
    clf.fit(xtr, tr["label"].to_numpy(int))
    return clf.predict_proba(xte)[:, 1]


def main() -> None:
    ensure_v1_dirs()
    master = load_master()
    folds = load_folds()
    base = build_base_feature_table(master)
    cols = group_cols(base)
    pred_rows, weight_rows = [], []
    for (split, fold), _ in folds.groupby(["split_name", "fold_id"], sort=False):
        train_ids, test_ids = get_fold_ids(master, folds, split, fold)
        train = base[base["sample_id"].isin(train_ids)].copy()
        test = base[base["sample_id"].isin(test_ids)].copy()
        if len(train) < 15 or len(test) < 2 or train["label"].nunique() < 2 or test["label"].nunique() < 2:
            continue
        outer_exp, outer_ctx = expert_scores(train, test, train, cols)
        inner_exp, inner_ctx = inner_oof(train, cols)
        nested_w, nested_table = choose_nested_weights(inner_exp)
        for _, row in nested_table.iterrows():
            weight_rows.append({"split_name": split, "fold_id": fold, "variant": "nested_candidate_grid", **row.to_dict()})
        scores = {
            "prespecified_equal_weight_C_QK_no_anchor": weighted_score(outer_exp, {"expert_counterfactual": 0.5, "expert_qk_no_anchor": 0.5}),
            "prespecified_C_0.5_QK_quantum_0.5": weighted_score(outer_exp, {"expert_counterfactual": 0.5, "expert_qk_quantum_only": 0.5}),
            "nested_learned_gate_C_QK_structure": weighted_score(outer_exp, nested_w),
            "stacked_calibrated_meta_model_lowdim_only": stacked_meta(inner_exp, inner_ctx, outer_exp, outer_ctx),
        }
        rule_score, rule_weights = rule_gated_score(outer_exp, outer_ctx)
        scores["rule_gated_C_QK_structure"] = rule_score
        mean_rule = {k: float(np.mean([w[k] for w in rule_weights])) for k in EXPERTS}
        weight_rows.append({"split_name": split, "fold_id": fold, "variant": "rule_gated_mean_weights", "candidate": "rule", "weights": json.dumps(mean_rule)})
        weight_rows.append({"split_name": split, "fold_id": fold, "variant": "nested_selected_weights", "candidate": "selected", "weights": json.dumps(nested_w)})
        for variant, score in scores.items():
            for sid, y, p in zip(test["sample_id"], test["label"], score):
                pred_rows.append({"split_name": split, "fold_id": fold, "sample_id": sid, "label": int(y), "model": variant, "score": float(p)})
    pred = pd.DataFrame(pred_rows)
    pred.to_csv(V1 / "gated_moe_predictions.tsv", sep="\t", index=False)
    metric_rows = []
    for (split, model), sub in pred.groupby(["split_name", "model"]):
        metric_rows.append({"split_name": split, "model": model, **metrics(sub["label"].to_numpy(), sub["score"].to_numpy())})
    met = pd.DataFrame(metric_rows).sort_values(["split_name", "AUPRC"], ascending=[True, False])
    met.to_csv(V1 / "gated_moe_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(weight_rows).to_csv(V1 / "gated_moe_fold_weights.tsv", sep="\t", index=False)
    lines = [
        "# CROSS-Neo v1 Gated MoE Report",
        "",
        "All outer-fold scores are trained from the outer train fold only. Nested learned weights are selected from inner OOF predictions within the outer train fold.",
        "",
        "## Top Gated Models",
        met.sort_values("AUPRC", ascending=False).head(15).to_markdown(index=False),
        "",
        "Prespecified and rule-gated variants are eligible for conservative interpretation; nested/stacked variants are fold-safe but still small-n.",
    ]
    (V1 / "gated_moe_report.md").write_text("\n".join(lines) + "\n")
    print(f"[v1-gated] predictions={len(pred)} metrics={len(met)}")


if __name__ == "__main__":
    main()
