#!/usr/bin/env python3
"""Reviewer-style stress battery for the spatial CNV ecosystem claim.

This script deliberately attacks the claim with controls a reviewer is likely
to ask for: slide-level effects, condition strata, QC/composition residuals,
spatial block permutation, leave-one-slide influence, and QC predictability of
the territory label.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import OneHotEncoder, StandardScaler


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
ARM = ROOT / "project/results/high_impact_topic_pilots_2026_05_08/spatial_cnv_territory_gate_2026_05_08/spatial_cnv_territory_per_spot.tsv.gz"
INFER = ROOT / "project/results/high_impact_topic_pilots_2026_05_08/infercnv_lite_fast_2026_05_08/infercnv_lite_fast_per_spot.tsv.gz"
ST = ROOT / "project/results/high_impact_topic_pilots_2026_05_08/top_topic_validation_2026_05_08/spatial_expression_cnv_per_spot.tsv.gz"
QC = ROOT / "project/results/high_impact_topic_pilots_2026_05_08/infercnv_lite_qc_residual_gate_2026_05_08/spot_qc.tsv.gz"
OUT = ROOT / "project/results/high_impact_topic_pilots_2026_05_08/spatial_cnv_reviewer_stress_battery_2026_05_08"
REPORT = ROOT / "project/reports/2026_05_08_spatial_cnv_reviewer_stress_battery.md"

OUTCOMES = [
    "TACSTD2_z",
    "TROP2_raw",
    "Macrophage_TAM_z",
    "TLS_B_z",
    "Tcell_cytotoxic_z",
    "RAI_thyroid_z",
    "EMT_stress_z",
]
KEY_OUTCOMES = ["TROP2_raw", "TACSTD2_z", "Macrophage_TAM_z", "Tcell_cytotoxic_z", "TLS_B_z"]
QC_AXES = ["log_total_counts", "n_genes_by_counts", "pct_mito"]
COMPOSITION_AXES = ["Epithelial_score", "Proliferation_score", "EMT_score", "CAF_ECM_score", "Hypoxia_score"]
COVARIATE_SETS = ["qc", "qc_composition"]


def fmt(x) -> str:
    if pd.isna(x):
        return "NA"
    if isinstance(x, (float, np.floating)):
        if abs(x) > 0 and abs(x) < 1e-4:
            return f"{x:.2e}"
        return f"{x:.3g}"
    return str(x)


def cohen_d(a, b) -> float:
    a = pd.Series(a).dropna().astype(float)
    b = pd.Series(b).dropna().astype(float)
    if len(a) < 2 or len(b) < 2:
        return np.nan
    pooled = ((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2)
    if pooled <= 0 or np.isnan(pooled):
        return np.nan
    return float((a.mean() - b.mean()) / math.sqrt(pooled))


def residualize(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    if x.shape[1] == 0:
        return y - np.mean(y)
    model = LinearRegression()
    model.fit(x, y)
    return y - model.predict(x)


def beta_after_covariates(y: np.ndarray, high: np.ndarray, cov: np.ndarray) -> tuple[float, float, float]:
    high_res = residualize(high.astype(float), cov)
    y_res = residualize(y.astype(float), cov)
    if np.std(high_res) <= 0 or np.std(y_res) <= 0:
        return np.nan, np.nan, np.nan
    beta = float(np.dot(high_res, y_res) / np.dot(high_res, high_res))
    r, p = stats.pearsonr(high_res, y_res)
    return beta, float(r), float(p)


def sign_test_p(values: pd.Series) -> float:
    v = pd.Series(values).dropna()
    v = v[v != 0]
    if v.empty:
        return np.nan
    k = int((v > 0).sum())
    return float(stats.binomtest(k, n=int(v.shape[0]), p=0.5, alternative="two-sided").pvalue)


def wilcoxon_p(values: pd.Series) -> float:
    v = pd.Series(values).dropna()
    v = v[v != 0]
    if v.shape[0] < 2:
        return np.nan
    try:
        return float(stats.wilcoxon(v).pvalue)
    except ValueError:
        return np.nan


def load_dat() -> pd.DataFrame:
    st_cols = ["sample_id", "spot_id"] + [c for c in COMPOSITION_AXES + OUTCOMES if c]
    st = pd.read_csv(ST, sep="\t", usecols=lambda c: c in st_cols)
    st = st.loc[:, ~st.columns.duplicated()].copy()

    qc = pd.read_csv(QC, sep="\t")
    frames = []

    arm = pd.read_csv(ARM, sep="\t")
    arm = arm.loc[:, ~arm.columns.duplicated()].copy()
    arm["territory_label_set"] = "arm_level_expression_cnv"
    arm["territory"] = arm["cnv_territory"].astype(str)
    arm["territory_high"] = arm["territory"].str.contains("high").astype(int)
    cols = [
        "sample_id",
        "spot_id",
        "condition",
        "array_row",
        "array_col",
        "territory_label_set",
        "territory",
        "territory_high",
    ] + [c for c in OUTCOMES + COMPOSITION_AXES if c in arm.columns]
    frames.append(arm[cols].copy())

    infer = pd.read_csv(INFER, sep="\t")
    infer = infer.loc[:, ~infer.columns.duplicated()].copy()
    infer["territory_label_set"] = "gene_bin_infercnv_lite"
    infer["territory"] = infer["infercnv_t00_territory"].astype(str)
    infer["territory_high"] = infer["infercnv_t00_high"].astype(int)
    cols = [
        "sample_id",
        "spot_id",
        "condition",
        "array_row",
        "array_col",
        "territory_label_set",
        "territory",
        "territory_high",
    ] + [c for c in OUTCOMES + COMPOSITION_AXES if c in infer.columns]
    frames.append(infer[cols].copy())

    out_frames = []
    for tmp in frames:
        missing = [c for c in OUTCOMES + COMPOSITION_AXES if c not in tmp.columns]
        if missing:
            tmp = tmp.merge(st[["sample_id", "spot_id"] + missing], on=["sample_id", "spot_id"], how="left")
        out_frames.append(tmp)

    dat = pd.concat(out_frames, ignore_index=True)
    dat = dat.merge(qc, on=["sample_id", "spot_id", "condition"], how="inner")
    dat = dat[dat["condition"].isin(["PTC", "LPTC", "ATC"])].copy()
    return dat.loc[:, ~dat.columns.duplicated()]


def covariate_matrix(dat: pd.DataFrame, covariate_set: str, include_sample: bool) -> np.ndarray:
    blocks = []
    if include_sample:
        enc = OneHotEncoder(sparse_output=False, handle_unknown="ignore")
        blocks.append(enc.fit_transform(dat[["sample_id"]]))
    if covariate_set in ["qc", "qc_composition"]:
        x = dat[QC_AXES].replace([np.inf, -np.inf], np.nan).fillna(dat[QC_AXES].median(numeric_only=True)).to_numpy(float)
        blocks.append(StandardScaler().fit_transform(x))
    if covariate_set == "qc_composition":
        comp_cols = [c for c in COMPOSITION_AXES if c in dat.columns]
        x = dat[comp_cols].replace([np.inf, -np.inf], np.nan).fillna(dat[comp_cols].median(numeric_only=True)).to_numpy(float)
        blocks.append(StandardScaler().fit_transform(x))
    if not blocks:
        return np.empty((dat.shape[0], 0))
    return np.hstack(blocks)


def slide_high_low_effects(dat: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    for (label_set, sample_id, condition), sub in dat.groupby(["territory_label_set", "sample_id", "condition"]):
        high = sub[sub["territory_high"] == 1]
        low = sub[sub["territory"].str.contains("low")]
        for outcome in [c for c in OUTCOMES if c in sub.columns]:
            if high[outcome].dropna().shape[0] < 5 or low[outcome].dropna().shape[0] < 5:
                continue
            rows.append(
                {
                    "territory_label_set": label_set,
                    "sample_id": sample_id,
                    "condition": condition,
                    "outcome": outcome,
                    "n_high": high[outcome].dropna().shape[0],
                    "n_low": low[outcome].dropna().shape[0],
                    "delta_high_low": high[outcome].mean() - low[outcome].mean(),
                    "cohen_d": cohen_d(high[outcome], low[outcome]),
                    "mw_p": float(stats.mannwhitneyu(high[outcome].dropna(), low[outcome].dropna()).pvalue),
                }
            )
    eff = pd.DataFrame(rows)
    eff.to_csv(OUT / "slide_high_low_effects.tsv", sep="\t", index=False)

    summaries = []
    for keys, sub in eff.groupby(["territory_label_set", "outcome"]):
        label_set, outcome = keys
        for scope, s2 in [("all_cancer_slides", sub)] + [(f"condition_{c}", x) for c, x in sub.groupby("condition")]:
            vals = s2["delta_high_low"].dropna()
            if vals.empty:
                continue
            summaries.append(
                {
                    "territory_label_set": label_set,
                    "outcome": outcome,
                    "scope": scope,
                    "n_slides": vals.shape[0],
                    "median_delta": float(vals.median()),
                    "mean_delta": float(vals.mean()),
                    "q25_delta": float(vals.quantile(0.25)),
                    "q75_delta": float(vals.quantile(0.75)),
                    "n_positive": int((vals > 0).sum()),
                    "sign_test_p": sign_test_p(vals),
                    "wilcoxon_p": wilcoxon_p(vals),
                }
            )
    summary = pd.DataFrame(summaries)
    summary.to_csv(OUT / "slide_high_low_summary.tsv", sep="\t", index=False)
    return eff, summary


def slide_residual_effects(dat: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    for (label_set, sample_id, condition), sub in dat.groupby(["territory_label_set", "sample_id", "condition"]):
        high = sub["territory_high"].astype(float).to_numpy()
        if np.unique(high).size < 2:
            continue
        for cov_set in COVARIATE_SETS:
            cov = covariate_matrix(sub, cov_set, include_sample=False)
            for outcome in [c for c in OUTCOMES if c in sub.columns]:
                mask = sub[outcome].notna().to_numpy()
                if mask.sum() < 30 or np.unique(high[mask]).size < 2:
                    continue
                beta, r, p = beta_after_covariates(sub.loc[mask, outcome].to_numpy(float), high[mask], cov[mask])
                rows.append(
                    {
                        "territory_label_set": label_set,
                        "sample_id": sample_id,
                        "condition": condition,
                        "covariate_set": cov_set,
                        "outcome": outcome,
                        "n": int(mask.sum()),
                        "beta_high": beta,
                        "pearson_r": r,
                        "p": p,
                    }
                )
    eff = pd.DataFrame(rows)
    eff.to_csv(OUT / "slide_residual_effects.tsv", sep="\t", index=False)

    summaries = []
    for keys, sub in eff.groupby(["territory_label_set", "covariate_set", "outcome"]):
        label_set, cov_set, outcome = keys
        for scope, s2 in [("all_cancer_slides", sub)] + [(f"condition_{c}", x) for c, x in sub.groupby("condition")]:
            vals = s2["beta_high"].dropna()
            if vals.empty:
                continue
            summaries.append(
                {
                    "territory_label_set": label_set,
                    "covariate_set": cov_set,
                    "outcome": outcome,
                    "scope": scope,
                    "n_slides": vals.shape[0],
                    "median_beta": float(vals.median()),
                    "mean_beta": float(vals.mean()),
                    "q25_beta": float(vals.quantile(0.25)),
                    "q75_beta": float(vals.quantile(0.75)),
                    "n_positive": int((vals > 0).sum()),
                    "sign_test_p": sign_test_p(vals),
                    "wilcoxon_p": wilcoxon_p(vals),
                }
            )
    summary = pd.DataFrame(summaries)
    summary.to_csv(OUT / "slide_residual_summary.tsv", sep="\t", index=False)
    return eff, summary


def condition_pooled_effects(dat: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (label_set, condition), sub in dat.groupby(["territory_label_set", "condition"]):
        high = sub["territory_high"].astype(float).to_numpy()
        if np.unique(high).size < 2:
            continue
        for cov_set in COVARIATE_SETS:
            cov = covariate_matrix(sub, cov_set, include_sample=True)
            for outcome in [c for c in OUTCOMES if c in sub.columns]:
                mask = sub[outcome].notna().to_numpy()
                if mask.sum() < 50 or np.unique(high[mask]).size < 2:
                    continue
                beta, r, p = beta_after_covariates(sub.loc[mask, outcome].to_numpy(float), high[mask], cov[mask])
                rows.append(
                    {
                        "territory_label_set": label_set,
                        "condition": condition,
                        "covariate_set": cov_set,
                        "outcome": outcome,
                        "n": int(mask.sum()),
                        "beta_high": beta,
                        "pearson_r": r,
                        "p": p,
                    }
                )
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "condition_pooled_residual_effects.tsv", sep="\t", index=False)
    return out


def pooled_beta(dat: pd.DataFrame, label_set: str, outcome: str, cov_set: str, drop_sample: str | None = None) -> tuple[float, float, float, int]:
    sub = dat[dat["territory_label_set"] == label_set].copy()
    if drop_sample is not None:
        sub = sub[sub["sample_id"] != drop_sample].copy()
    mask = sub[outcome].notna().to_numpy()
    high = sub["territory_high"].astype(float).to_numpy()
    if mask.sum() < 50 or np.unique(high[mask]).size < 2:
        return np.nan, np.nan, np.nan, int(mask.sum())
    cov = covariate_matrix(sub, cov_set, include_sample=True)
    beta, r, p = beta_after_covariates(sub.loc[mask, outcome].to_numpy(float), high[mask], cov[mask])
    return beta, r, p, int(mask.sum())


def loso_influence(dat: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for label_set in sorted(dat["territory_label_set"].unique()):
        samples = sorted(dat.loc[dat["territory_label_set"] == label_set, "sample_id"].unique())
        for outcome in [c for c in OUTCOMES if c in dat.columns]:
            full_beta, full_r, full_p, full_n = pooled_beta(dat, label_set, outcome, "qc_composition")
            for sample_id in samples:
                beta, r, p, n = pooled_beta(dat, label_set, outcome, "qc_composition", drop_sample=sample_id)
                rows.append(
                    {
                        "territory_label_set": label_set,
                        "outcome": outcome,
                        "dropped_sample": sample_id,
                        "full_beta": full_beta,
                        "drop_beta": beta,
                        "delta_beta_vs_full": beta - full_beta if not pd.isna(beta) and not pd.isna(full_beta) else np.nan,
                        "drop_r": r,
                        "drop_p": p,
                        "n": n,
                        "sign_changed": bool(np.sign(beta) != np.sign(full_beta)) if not pd.isna(beta) and not pd.isna(full_beta) else False,
                    }
                )
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "loso_influence.tsv", sep="\t", index=False)
    return out


def add_spatial_blocks(dat: pd.DataFrame, q: int = 6) -> pd.DataFrame:
    parts = []
    for (_, sample_id), sub in dat.groupby(["territory_label_set", "sample_id"]):
        sub = sub.copy()
        row_rank = sub["array_row"].rank(method="first")
        col_rank = sub["array_col"].rank(method="first")
        sub["_row_block"] = pd.qcut(row_rank, q=min(q, sub["array_row"].nunique()), labels=False, duplicates="drop")
        sub["_col_block"] = pd.qcut(col_rank, q=min(q, sub["array_col"].nunique()), labels=False, duplicates="drop")
        sub["spatial_block"] = sub["_row_block"].astype(str) + "_" + sub["_col_block"].astype(str)
        parts.append(sub.drop(columns=["_row_block", "_col_block"]))
    return pd.concat(parts, ignore_index=True)


def block_permutation(dat: pd.DataFrame, n_perm: int = 400) -> pd.DataFrame:
    rng = np.random.default_rng(8)
    dat = add_spatial_blocks(dat)
    group_cols = ["territory_label_set", "sample_id", "condition", "spatial_block"]
    agg_cols = {"territory_high": "mean", "log_total_counts": "mean", "n_genes_by_counts": "mean", "pct_mito": "mean"}
    for c in COMPOSITION_AXES + OUTCOMES:
        if c in dat.columns:
            agg_cols[c] = "mean"
    block = dat.groupby(group_cols).agg(agg_cols).reset_index()
    counts = dat.groupby(group_cols).size().reset_index(name="n_spots")
    block = block.merge(counts, on=group_cols, how="left")
    block = block[block["n_spots"] >= 10].copy()

    rows = []
    for label_set, sub0 in block.groupby("territory_label_set"):
        sub0 = sub0.reset_index(drop=True)
        for outcome in [c for c in OUTCOMES if c in sub0.columns]:
            sub = sub0[sub0[outcome].notna()].copy().reset_index(drop=True)
            if sub.shape[0] < 50 or sub["territory_high"].nunique() < 2:
                continue
            cov = covariate_matrix(sub.rename(columns={"territory_high": "territory_high"}), "qc_composition", include_sample=True)
            y = sub[outcome].to_numpy(float)
            high = sub["territory_high"].to_numpy(float)
            obs_beta, obs_r, _ = beta_after_covariates(y, high, cov)
            null = []
            for _ in range(n_perm):
                perm_high = high.copy()
                for sample_id, idx in sub.groupby("sample_id").groups.items():
                    idx_arr = np.asarray(list(idx), dtype=int)
                    perm_high[idx_arr] = rng.permutation(perm_high[idx_arr])
                beta, _, _ = beta_after_covariates(y, perm_high, cov)
                null.append(beta)
            null = np.asarray(null, dtype=float)
            p_two = float((1 + np.sum(np.abs(null) >= abs(obs_beta))) / (n_perm + 1))
            rows.append(
                {
                    "territory_label_set": label_set,
                    "outcome": outcome,
                    "n_blocks": sub.shape[0],
                    "obs_beta": obs_beta,
                    "obs_r": obs_r,
                    "null_mean_beta": float(np.nanmean(null)),
                    "null_sd_beta": float(np.nanstd(null, ddof=1)),
                    "block_perm_p_two_sided": p_two,
                    "n_perm": n_perm,
                }
            )
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "spatial_block_permutation.tsv", sep="\t", index=False)
    return out


def qc_predictability(dat: pd.DataFrame) -> pd.DataFrame:
    rows = []
    feature_cols = QC_AXES + [c for c in COMPOSITION_AXES if c in dat.columns]
    for label_set, d0 in dat.groupby("territory_label_set"):
        for sample_id, sub in d0.groupby("sample_id"):
            y = sub["territory_high"].astype(int).to_numpy()
            if np.unique(y).size < 2 or sub.shape[0] < 50:
                continue
            x = sub[feature_cols].replace([np.inf, -np.inf], np.nan).fillna(sub[feature_cols].median(numeric_only=True)).to_numpy(float)
            x = StandardScaler().fit_transform(x)
            model = LogisticRegression(max_iter=500, class_weight="balanced", C=0.5, solver="liblinear")
            model.fit(x, y)
            prob = model.predict_proba(x)[:, 1]
            rows.append(
                {
                    "territory_label_set": label_set,
                    "mode": "within_slide_apparent",
                    "sample_id": sample_id,
                    "n": sub.shape[0],
                    "AUROC": roc_auc_score(y, prob),
                    "high_rate": y.mean(),
                }
            )

        pred_frames = []
        for held in sorted(d0["sample_id"].unique()):
            train = d0["sample_id"] != held
            test = ~train
            ytr = d0.loc[train, "territory_high"].astype(int).to_numpy()
            yte = d0.loc[test, "territory_high"].astype(int).to_numpy()
            if np.unique(ytr).size < 2 or np.unique(yte).size < 2:
                continue
            x = d0[feature_cols].replace([np.inf, -np.inf], np.nan).fillna(d0[feature_cols].median(numeric_only=True)).to_numpy(float)
            scaler = StandardScaler().fit(x[train])
            xtr = scaler.transform(x[train])
            xte = scaler.transform(x[test])
            model = LogisticRegression(max_iter=500, class_weight="balanced", C=0.5, solver="liblinear")
            model.fit(xtr, ytr)
            prob = model.predict_proba(xte)[:, 1]
            pred_frames.append(pd.DataFrame({"sample_id": held, "observed": yte, "predicted": prob}))
        if pred_frames:
            pred = pd.concat(pred_frames, ignore_index=True)
            rows.append(
                {
                    "territory_label_set": label_set,
                    "mode": "leave_one_slide_out",
                    "sample_id": "pooled",
                    "n": pred.shape[0],
                    "AUROC": roc_auc_score(pred["observed"].astype(int), pred["predicted"]),
                    "high_rate": pred["observed"].mean(),
                }
            )
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "qc_composition_predicts_territory.tsv", sep="\t", index=False)
    return out


def stress_scorecard(
    slide_summary: pd.DataFrame,
    cond: pd.DataFrame,
    block: pd.DataFrame,
    loso: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    ss = slide_summary[
        (slide_summary["covariate_set"] == "qc_composition")
        & (slide_summary["scope"] == "all_cancer_slides")
        & (slide_summary["outcome"].isin(KEY_OUTCOMES))
    ]
    for _, r in ss.iterrows():
        label_set = r["territory_label_set"]
        outcome = r["outcome"]
        csub = cond[
            (cond["territory_label_set"] == label_set)
            & (cond["covariate_set"] == "qc_composition")
            & (cond["outcome"] == outcome)
        ]
        bsub = block[(block["territory_label_set"] == label_set) & (block["outcome"] == outcome)]
        lsub = loso[(loso["territory_label_set"] == label_set) & (loso["outcome"] == outcome)]
        condition_positive = int((csub["beta_high"] > 0).sum()) if not csub.empty else 0
        n_conditions = int(csub["condition"].nunique()) if not csub.empty else 0
        block_p = float(bsub["block_perm_p_two_sided"].iloc[0]) if not bsub.empty else np.nan
        block_beta = float(bsub["obs_beta"].iloc[0]) if not bsub.empty else np.nan
        n_sign_change = int(lsub["sign_changed"].sum()) if not lsub.empty else np.nan
        min_loso_beta = float(lsub["drop_beta"].min()) if not lsub.empty else np.nan
        max_loso_beta = float(lsub["drop_beta"].max()) if not lsub.empty else np.nan

        median_beta = float(r["median_beta"])
        n_positive = int(r["n_positive"])
        n_slides = int(r["n_slides"])
        block_positive = bool(block_beta > 0) if not pd.isna(block_beta) else False
        if median_beta > 0 and not block_positive and not pd.isna(block_beta):
            verdict = "directionally_inconsistent"
        elif median_beta > 0 and n_positive >= 9 and n_conditions == condition_positive and block_p <= 0.05 and n_sign_change == 0:
            verdict = "survives_strong"
        elif median_beta > 0 and n_positive >= 7 and condition_positive >= max(1, n_conditions - 1) and (pd.isna(block_p) or block_p <= 0.15):
            verdict = "survives_moderate"
        elif median_beta > 0:
            verdict = "weak_positive"
        else:
            verdict = "fails_or_reverses"

        rows.append(
            {
                "territory_label_set": label_set,
                "outcome": outcome,
                "median_slide_beta": median_beta,
                "n_positive_slides": n_positive,
                "n_slides": n_slides,
                "slide_sign_test_p": r["sign_test_p"],
                "slide_wilcoxon_p": r["wilcoxon_p"],
                "condition_positive": condition_positive,
                "n_conditions": n_conditions,
                "block_beta": block_beta,
                "block_perm_p": block_p,
                "loso_sign_changes": n_sign_change,
                "min_loso_beta": min_loso_beta,
                "max_loso_beta": max_loso_beta,
                "verdict": verdict,
            }
        )
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "reviewer_stress_scorecard.tsv", sep="\t", index=False)
    return out


def write_report(
    dat: pd.DataFrame,
    raw_summary: pd.DataFrame,
    slide_summary: pd.DataFrame,
    cond: pd.DataFrame,
    block: pd.DataFrame,
    qc_pred: pd.DataFrame,
    scorecard: pd.DataFrame,
) -> None:
    lines = [
        "# Spatial CNV reviewer stress battery - 2026-05-08\n",
        "Purpose: attack the spatial CNV ecosystem claim with slide-level, QC/composition, condition-stratified, spatial-block, and influence analyses.\n",
        "## Verdict\n",
    ]
    for label_set in ["arm_level_expression_cnv", "gene_bin_infercnv_lite"]:
        lines.append(f"### {label_set}")
        q = qc_pred[(qc_pred["territory_label_set"] == label_set) & (qc_pred["mode"] == "within_slide_apparent")]
        loso_q = qc_pred[(qc_pred["territory_label_set"] == label_set) & (qc_pred["mode"] == "leave_one_slide_out")]
        if not q.empty:
            lines.append(f"- QC/composition predicts territory within slide: median apparent AUROC **{fmt(q['AUROC'].median())}**.")
        if not loso_q.empty:
            lines.append(f"- QC/composition leave-one-slide-out AUROC: **{fmt(loso_q['AUROC'].iloc[0])}**.")
        for outcome in KEY_OUTCOMES:
            row = scorecard[(scorecard["territory_label_set"] == label_set) & (scorecard["outcome"] == outcome)]
            if row.empty:
                continue
            r = row.iloc[0]
            lines.append(
                f"- {outcome}: **{r['verdict']}**; median slide beta {fmt(r['median_slide_beta'])}, "
                f"positive slides {int(r['n_positive_slides'])}/{int(r['n_slides'])}, "
                f"condition+ {int(r['condition_positive'])}/{int(r['n_conditions'])}, "
                f"block p {fmt(r['block_perm_p'])}, LOSO sign changes {int(r['loso_sign_changes'])}."
            )
    lines.append("\nInterpretation: treat TROP2/TACSTD2 as QC-sensitive exploratory phenotypes unless they survive slide-level and block-level tests. The most defensible biology is whatever survives as positive across slides, conditions, block permutation, and leave-one-slide influence.\n")

    lines.extend(["## Scorecard\n", "| Label set | Outcome | Verdict | median beta | +slides | +conditions | block beta | block p | LOSO sign changes |", "|---|---|---|---:|---:|---:|---:|---:|---:|"])
    for _, r in scorecard.iterrows():
        lines.append(
            f"| {r['territory_label_set']} | {r['outcome']} | {r['verdict']} | {fmt(r['median_slide_beta'])} | "
            f"{int(r['n_positive_slides'])}/{int(r['n_slides'])} | {int(r['condition_positive'])}/{int(r['n_conditions'])} | "
            f"{fmt(r['block_beta'])} | {fmt(r['block_perm_p'])} | {int(r['loso_sign_changes'])} |"
        )

    lines.extend(["\n## Raw Slide High-Low Summary\n", "| Label set | Outcome | Scope | median delta | +slides | sign p | Wilcoxon p |", "|---|---|---|---:|---:|---:|---:|"])
    show_raw = raw_summary[(raw_summary["outcome"].isin(KEY_OUTCOMES)) & (raw_summary["scope"] == "all_cancer_slides")]
    for _, r in show_raw.iterrows():
        lines.append(
            f"| {r['territory_label_set']} | {r['outcome']} | {r['scope']} | {fmt(r['median_delta'])} | "
            f"{int(r['n_positive'])}/{int(r['n_slides'])} | {fmt(r['sign_test_p'])} | {fmt(r['wilcoxon_p'])} |"
        )

    lines.extend(["\n## QC/Composition Territory Predictability\n", "| Label set | Mode | n | AUROC | high rate |", "|---|---|---:|---:|---:|"])
    qshow = pd.concat(
        [
            qc_pred[qc_pred["mode"] == "leave_one_slide_out"],
            qc_pred[qc_pred["mode"] == "within_slide_apparent"]
            .groupby(["territory_label_set", "mode"])
            .agg(n=("n", "sum"), AUROC=("AUROC", "median"), high_rate=("high_rate", "median"))
            .reset_index()
            .assign(sample_id="median_slide"),
        ],
        ignore_index=True,
    )
    for _, r in qshow.iterrows():
        lines.append(f"| {r['territory_label_set']} | {r['mode']} | {int(r['n'])} | {fmt(r['AUROC'])} | {fmt(r['high_rate'])} |")

    lines.extend(["\n## Output files"])
    for p in sorted(OUT.glob("*")):
        lines.append(f"- `{p}`")
    REPORT.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    dat = load_dat()
    raw_effects, raw_summary = slide_high_low_effects(dat)
    slide_effects, slide_summary = slide_residual_effects(dat)
    cond = condition_pooled_effects(dat)
    loso = loso_influence(dat)
    block = block_permutation(dat, n_perm=400)
    qc_pred = qc_predictability(dat)
    scorecard = stress_scorecard(slide_summary, cond, block, loso)
    write_report(dat, raw_summary, slide_summary, cond, block, qc_pred, scorecard)
    print(f"[write] {OUT}")
    print(f"[write] {REPORT}")


if __name__ == "__main__":
    main()
