#!/usr/bin/env python3
"""Link CLEAN-NeoBench errors to training/data distribution.

This script turns the existing failure-mode tables into a reviewer-facing
diagnostic audit: which methods are vulnerable to source priors, low-prevalence
settings, rare HLA alleles, leakage-risk regions, and external/heldout shift.
It does not train a new predictor.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from common import dataframe_to_markdown, ensure_dir, update_manifest, write_tsv


SLICE_FEATURES = [
    "distribution_partition",
    "source_name",
    "hla_allele_4digit",
    "hla_supertype",
    "peptide_length",
    "leakage_risk_level",
    "split_low_prevalence",
    "exact_peptide_hla_train_overlap",
    "near_peptide_train_overlap",
    "mhc_class",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=".", help="Repository root")
    parser.add_argument("--output-root", required=True, help="CLEAN-NeoBench output directory")
    parser.add_argument("--rare-hla-support", type=int, default=20, help="HLA support threshold for rare allele diagnostics")
    parser.add_argument("--low-hla-support", type=int, default=50, help="HLA support threshold for low-support diagnostics")
    parser.add_argument("--low-prevalence-cutoff", type=float, default=0.10, help="Source prevalence cutoff")
    return parser.parse_args()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def finite_float(value: Any, default: float = math.nan) -> float:
    try:
        x = float(value)
        return x if math.isfinite(x) else default
    except Exception:
        return default


def boolish_series(s: pd.Series) -> pd.Series:
    return s.astype(str).str.lower().isin(["true", "1", "yes", "y"])


def failure_mask(df: pd.DataFrame, token: str) -> pd.Series:
    if "primary_failure_type" not in df.columns:
        return pd.Series(False, index=df.index)
    primary = df["primary_failure_type"].fillna("").astype(str).eq(token)
    tags = df.get("failure_tags", pd.Series("", index=df.index)).fillna("").astype(str).str.contains(token, regex=False)
    return primary | tags


def summarize_slice(g: pd.DataFrame, feature: str, value: str, method: str) -> dict[str, Any]:
    label = pd.to_numeric(g["label"], errors="coerce")
    neg = label.eq(0)
    pos = label.eq(1)
    high_fp = failure_mask(g, "high_ranked_negative_top10pct")
    miss_pos = failure_mask(g, "missed_positive_bottom50pct")
    thr_fp = failure_mask(g, "threshold_false_positive_score_ge_0.5")
    thr_fn = failure_mask(g, "threshold_false_negative_score_lt_0.5")
    n_neg = int(neg.sum())
    n_pos = int(pos.sum())
    fp_rate = float((high_fp & neg).sum() / n_neg) if n_neg else math.nan
    fn_rate = float((miss_pos & pos).sum() / n_pos) if n_pos else math.nan
    threshold_fp_rate = float((thr_fp & neg).sum() / n_neg) if n_neg else math.nan
    threshold_fn_rate = float((thr_fn & pos).sum() / n_pos) if n_pos else math.nan
    parts = [x for x in [fp_rate, fn_rate, threshold_fp_rate, threshold_fn_rate] if math.isfinite(x)]
    error_pressure = float(np.mean(parts)) if parts else math.nan
    return {
        "method_name": method,
        "method_role": g["method_role"].iloc[0] if "method_role" in g else "",
        "method_family": g["method_family"].iloc[0] if "method_family" in g else "",
        "feature": feature,
        "value": value,
        "n_scored": int(len(g)),
        "n_pos": n_pos,
        "n_neg": n_neg,
        "positive_prevalence": float(label.mean()) if label.notna().any() else math.nan,
        "mean_score": float(pd.to_numeric(g["score_for_analysis"], errors="coerce").mean()),
        "median_score": float(pd.to_numeric(g["score_for_analysis"], errors="coerce").median()),
        "mean_rank_percentile": float(pd.to_numeric(g["rank_percentile"], errors="coerce").mean()),
        "high_ranked_negative_top10pct": int((high_fp & neg).sum()),
        "high_ranked_negative_rate_among_neg": fp_rate,
        "missed_positive_bottom50pct": int((miss_pos & pos).sum()),
        "missed_positive_rate_among_pos": fn_rate,
        "threshold_false_positive_score_ge_0.5": int((thr_fp & neg).sum()),
        "threshold_false_positive_rate_among_neg": threshold_fp_rate,
        "threshold_false_negative_score_lt_0.5": int((thr_fn & pos).sum()),
        "threshold_false_negative_rate_among_pos": threshold_fn_rate,
        "error_pressure_score": error_pressure,
        "mean_source_positive_prevalence": float(pd.to_numeric(g.get("source_positive_prevalence", np.nan), errors="coerce").mean()),
        "mean_hla_allele_support_count": float(pd.to_numeric(g.get("hla_allele_support_count", np.nan), errors="coerce").mean()),
    }


def build_error_slice_matrix(scored: pd.DataFrame) -> pd.DataFrame:
    if scored.empty:
        return pd.DataFrame()
    rows = []
    for method, mg in scored.groupby("method_name", dropna=False):
        for feature in SLICE_FEATURES:
            if feature not in mg.columns:
                continue
            values = mg[feature].fillna("NA").astype(str)
            for value, gidx in values.groupby(values).groups.items():
                g = mg.loc[list(gidx)]
                rows.append(summarize_slice(g, feature, value, str(method)))
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.sort_values(["error_pressure_score", "n_scored"], ascending=[False, False])


def top_methods_for_slice(slice_rows: pd.DataFrame, metric: str) -> str:
    if slice_rows.empty or metric not in slice_rows:
        return ""
    view = slice_rows.sort_values(metric, ascending=False)
    pairs = []
    for _, r in view.head(5).iterrows():
        val = finite_float(r.get(metric), math.nan)
        if math.isfinite(val):
            pairs.append(f"{r['method_name']}={val:.3f}")
    return "; ".join(pairs)


def diagnose_shift_error(row: pd.Series) -> str:
    train_n = finite_float(row.get("train_n"), 0)
    ext_n = finite_float(row.get("external_n"), 0)
    abs_delta = finite_float(row.get("absolute_fraction_delta"), 0)
    prev_delta = finite_float(row.get("positive_prevalence_delta"), 0)
    fp = finite_float(row.get("max_high_ranked_negative_rate"), 0)
    fn = finite_float(row.get("max_missed_positive_rate"), 0)
    feature = str(row.get("feature", ""))
    value = str(row.get("value", ""))
    reasons = []
    if train_n == 0 and ext_n > 0:
        reasons.append("external/heldout-only slice; no local train support")
    elif ext_n == 0 and train_n > 0:
        reasons.append("train-only slice; external robustness untested")
    if abs_delta >= 0.20:
        reasons.append("large train-vs-external representation shift")
    if abs(prev_delta) >= 0.30:
        reasons.append("large positive-prevalence shift")
    if fp >= 0.20:
        reasons.append("enriched high-ranked negatives")
    if fn >= 0.40:
        reasons.append("enriched missed positives")
    if feature == "source_name" and "TESLA" in value and fp > 0:
        reasons.append("TESLA-like low-prevalence stress slice")
    if feature == "hla_allele_4digit" and train_n + ext_n < 20:
        reasons.append("rare allele slice")
    return "; ".join(dict.fromkeys(reasons)) if reasons else "monitor"


def build_shift_error_link(shift: pd.DataFrame, matrix: pd.DataFrame) -> pd.DataFrame:
    if shift.empty or matrix.empty:
        return pd.DataFrame()
    agg_rows = []
    for (feature, value), g in matrix.groupby(["feature", "value"], dropna=False):
        agg_rows.append(
            {
                "feature": feature,
                "value": value,
                "n_method_slices": int(g["method_name"].nunique()),
                "mean_error_pressure_score": float(pd.to_numeric(g["error_pressure_score"], errors="coerce").mean()),
                "max_error_pressure_score": float(pd.to_numeric(g["error_pressure_score"], errors="coerce").max()),
                "mean_high_ranked_negative_rate": float(pd.to_numeric(g["high_ranked_negative_rate_among_neg"], errors="coerce").mean()),
                "max_high_ranked_negative_rate": float(pd.to_numeric(g["high_ranked_negative_rate_among_neg"], errors="coerce").max()),
                "mean_missed_positive_rate": float(pd.to_numeric(g["missed_positive_rate_among_pos"], errors="coerce").mean()),
                "max_missed_positive_rate": float(pd.to_numeric(g["missed_positive_rate_among_pos"], errors="coerce").max()),
                "top_fp_methods": top_methods_for_slice(g, "high_ranked_negative_rate_among_neg"),
                "top_fn_methods": top_methods_for_slice(g, "missed_positive_rate_among_pos"),
            }
        )
    agg = pd.DataFrame(agg_rows)
    merged = shift.copy()
    merged["value"] = merged["value"].fillna("NA").astype(str)
    merged = merged.merge(agg, on=["feature", "value"], how="left")
    merged["positive_prevalence_delta"] = (
        pd.to_numeric(merged["external_positive_prevalence"], errors="coerce")
        - pd.to_numeric(merged["train_positive_prevalence"], errors="coerce")
    )
    merged["distribution_shift_score"] = (
        pd.to_numeric(merged["absolute_fraction_delta"], errors="coerce").fillna(0)
        + pd.to_numeric(merged["positive_prevalence_delta"], errors="coerce").abs().fillna(0)
        + np.log1p(pd.to_numeric(merged["external_over_train_fraction_ratio"], errors="coerce").replace([np.inf, -np.inf], 10).fillna(0)).clip(0, 3) / 3
    )
    merged["error_shift_link_score"] = (
        merged["distribution_shift_score"].fillna(0)
        * (1 + pd.to_numeric(merged["max_error_pressure_score"], errors="coerce").fillna(0))
    )
    merged["diagnosis"] = merged.apply(diagnose_shift_error, axis=1)
    return merged.sort_values(["error_shift_link_score", "external_n", "train_n"], ascending=[False, False, False])


def safe_corr(x: pd.Series, y: pd.Series) -> float:
    x = pd.to_numeric(x, errors="coerce")
    y = pd.to_numeric(y, errors="coerce")
    mask = x.notna() & y.notna()
    if mask.sum() < 3 or x[mask].nunique() < 2 or y[mask].nunique() < 2:
        return math.nan
    return float(x[mask].corr(y[mask]))


def rate_for_mask(df: pd.DataFrame, row_mask: pd.Series, failure: str, label_value: int) -> float:
    sub = df[row_mask & pd.to_numeric(df["label"], errors="coerce").eq(label_value)]
    if sub.empty:
        return math.nan
    return float(failure_mask(sub, failure).mean())


def build_method_vulnerability(scored: pd.DataFrame, rare_hla: int, low_hla: int, low_prev: float) -> pd.DataFrame:
    if scored.empty:
        return pd.DataFrame()
    rows = []
    for method, g in scored.groupby("method_name", dropna=False):
        source_summary = (
            g.groupby("source_name", dropna=False)
            .agg(
                n=("candidate_id", "size"),
                mean_score=("score_for_analysis", "mean"),
                source_positive_prevalence=("source_positive_prevalence", "mean"),
            )
            .reset_index()
        )
        score_source_prev_corr = safe_corr(source_summary["mean_score"], source_summary["source_positive_prevalence"])
        source_score_range = float(source_summary["mean_score"].max() - source_summary["mean_score"].min()) if len(source_summary) else math.nan
        source_error_rates = []
        for _, srow in source_summary.iterrows():
            src = srow["source_name"]
            sg = g[g["source_name"].eq(src)]
            neg = pd.to_numeric(sg["label"], errors="coerce").eq(0)
            pos = pd.to_numeric(sg["label"], errors="coerce").eq(1)
            fp_rate = float(failure_mask(sg[neg], "high_ranked_negative_top10pct").mean()) if neg.any() else math.nan
            fn_rate = float(failure_mask(sg[pos], "missed_positive_bottom50pct").mean()) if pos.any() else math.nan
            source_error_rates.extend([x for x in [fp_rate, fn_rate] if math.isfinite(x)])
        source_error_rate_range = float(max(source_error_rates) - min(source_error_rates)) if source_error_rates else math.nan

        source_prev = pd.to_numeric(g.get("source_positive_prevalence", pd.Series(np.nan, index=g.index)), errors="coerce")
        hla_support = pd.to_numeric(g.get("hla_allele_support_count", pd.Series(np.nan, index=g.index)), errors="coerce")
        low_prev_mask = source_prev < low_prev
        rare_hla_mask = hla_support < rare_hla
        low_hla_mask = hla_support < low_hla
        leakage_high_mask = g.get("leakage_risk_level", pd.Series("", index=g.index)).fillna("").astype(str).str.lower().eq("high")
        external_mask = g.get("distribution_partition", pd.Series("", index=g.index)).fillna("").astype(str).ne("local_train_pool")

        low_prev_fp = rate_for_mask(g, low_prev_mask, "high_ranked_negative_top10pct", 0)
        rare_hla_fp = rate_for_mask(g, rare_hla_mask, "high_ranked_negative_top10pct", 0)
        rare_hla_fn = rate_for_mask(g, rare_hla_mask, "missed_positive_bottom50pct", 1)
        low_hla_fn = rate_for_mask(g, low_hla_mask, "missed_positive_bottom50pct", 1)
        high_leakage_fp = rate_for_mask(g, leakage_high_mask, "high_ranked_negative_top10pct", 0)
        external_fp = rate_for_mask(g, external_mask, "high_ranked_negative_top10pct", 0)
        external_fn = rate_for_mask(g, external_mask, "missed_positive_bottom50pct", 1)

        components = [
            abs(score_source_prev_corr) if math.isfinite(score_source_prev_corr) else 0,
            source_score_range if math.isfinite(source_score_range) else 0,
            source_error_rate_range if math.isfinite(source_error_rate_range) else 0,
            low_prev_fp if math.isfinite(low_prev_fp) else 0,
            rare_hla_fp if math.isfinite(rare_hla_fp) else 0,
            rare_hla_fn if math.isfinite(rare_hla_fn) else 0,
            high_leakage_fp if math.isfinite(high_leakage_fp) else 0,
            max(external_fp if math.isfinite(external_fp) else 0, external_fn if math.isfinite(external_fn) else 0),
        ]
        vulnerability_score = float(np.mean(components))
        issues = []
        if math.isfinite(score_source_prev_corr) and abs(score_source_prev_corr) >= 0.45:
            issues.append("source-prior coupling")
        if math.isfinite(low_prev_fp) and low_prev_fp >= 0.15:
            issues.append("low-prevalence false-positive pressure")
        if max(rare_hla_fp if math.isfinite(rare_hla_fp) else 0, rare_hla_fn if math.isfinite(rare_hla_fn) else 0) >= 0.30:
            issues.append("rare-HLA instability")
        if math.isfinite(high_leakage_fp) and high_leakage_fp >= 0.10:
            issues.append("leakage-region false-positive pressure")
        if max(external_fp if math.isfinite(external_fp) else 0, external_fn if math.isfinite(external_fn) else 0) >= 0.30:
            issues.append("external/holdout fragility")
        rows.append(
            {
                "method_name": method,
                "method_role": g["method_role"].iloc[0] if "method_role" in g else "",
                "method_family": g["method_family"].iloc[0] if "method_family" in g else "",
                "n_scored": int(len(g)),
                "source_score_prevalence_corr": score_source_prev_corr,
                "source_score_range": source_score_range,
                "source_error_rate_range": source_error_rate_range,
                "low_prevalence_high_ranked_negative_rate": low_prev_fp,
                "rare_hla_high_ranked_negative_rate": rare_hla_fp,
                "rare_hla_missed_positive_rate": rare_hla_fn,
                "low_hla_missed_positive_rate": low_hla_fn,
                "high_leakage_high_ranked_negative_rate": high_leakage_fp,
                "external_high_ranked_negative_rate": external_fp,
                "external_missed_positive_rate": external_fn,
                "distribution_vulnerability_score": vulnerability_score,
                "primary_distribution_issues": "; ".join(issues) if issues else "no single dominant issue",
            }
        )
    return pd.DataFrame(rows).sort_values("distribution_vulnerability_score", ascending=False)


def build_remediation_plan(vulnerability: pd.DataFrame, link: pd.DataFrame) -> pd.DataFrame:
    rows = [
        {
            "priority": "P0",
            "failure_axis": "public training-corpus overlap unresolved",
            "evidence": "public_clean_comparators_allowed remains zero until corpus files are supplied",
            "recommended_action": "Collect official training corpora and run row-level peptide/HLA overlap audit before using public scores as clean features.",
            "barneo_policy": "Keep public predictors as caveated comparators only.",
        },
        {
            "priority": "P0",
            "failure_axis": "source prevalence shift",
            "evidence": "CEDAR is positive-rich while TESLA/ITSNdb validation slices are low-prevalence.",
            "recommended_action": "Require source-heldout calibration, source-balanced weighting, and low-prevalence top-k stress tests.",
            "barneo_policy": "Penalize low-prevalence false-positive pressure and keep abstention active.",
        },
        {
            "priority": "P1",
            "failure_axis": "HLA allele support imbalance",
            "evidence": "rare/low-support HLA slices show unstable false-positive or missed-positive behavior.",
            "recommended_action": "Add allele support count, Korean-HLA focus calibration, and allele-heldout reliability caps.",
            "barneo_policy": "Cap confidence for rare/underrepresented HLA alleles.",
        },
        {
            "priority": "P1",
            "failure_axis": "high leakage-risk regions",
            "evidence": "exact/near peptide-HLA overlap creates claim risk even when ranking looks strong.",
            "recommended_action": "Audit exact and near peptide-HLA overlaps before any clean benchmark claim.",
            "barneo_policy": "High leakage risk remains abstain/claim-blocked.",
        },
        {
            "priority": "P2",
            "failure_axis": "patient metadata missing",
            "evidence": "PAAD/THCA disease timing, presentation, antigen expression, immune context, and safety fields are mostly absent.",
            "recommended_action": "Add patient-level fields before patient-gated score can be interpreted.",
            "barneo_policy": "Patient-gated rows stay research_triage_only and clinical_use=false.",
        },
    ]
    if not vulnerability.empty:
        top = vulnerability.head(5)
        rows.append(
            {
                "priority": "P1",
                "failure_axis": "method-specific distribution vulnerability",
                "evidence": "; ".join(f"{r.method_name}: {r.primary_distribution_issues}" for r in top.itertuples()),
                "recommended_action": "Use method-specific reliability weights instead of a single global winner.",
                "barneo_policy": "BMA weights plus abstention should remain the primary deployment surface.",
            }
        )
    if not link.empty:
        top_link = link.head(5)
        rows.append(
            {
                "priority": "P1",
                "failure_axis": "shift-linked error slices",
                "evidence": "; ".join(f"{r.feature}={r.value}: {r.diagnosis}" for r in top_link.itertuples()),
                "recommended_action": "Build challenge splits around these slices and track them in every leaderboard.",
                "barneo_policy": "Use these slices as reviewer-safe stress tests, not as hidden post-hoc filters.",
            }
        )
    return pd.DataFrame(rows)


def write_report(output_root: Path, matrix: pd.DataFrame, link: pd.DataFrame, vulnerability: pd.DataFrame, remediation: pd.DataFrame) -> None:
    risky_slices = link[
        link["diagnosis"].fillna("").astype(str).ne("monitor")
    ].head(40) if not link.empty else link
    top_matrix = matrix[
        pd.to_numeric(matrix["n_scored"], errors="coerce").fillna(0) >= 10
    ].head(50) if not matrix.empty else matrix
    text = f"""# CLEAN-NeoBench Distribution Error Audit

## Purpose

This audit asks whether wrong predictions align with training-data distribution. It links method errors to source prevalence, HLA support, peptide length, leakage-risk, low-prevalence stress groups, and external/heldout slices. It does not train a new predictor and does not create a SOTA claim.

## Main Interpretation

The current practical risk is not simply "which method has the highest mean AUPRC." The main risk is that methods can win aggregate metrics while remaining fragile under source shift, low-prevalence settings, HLA underrepresentation, or overlap-risk regions.

## Method Distribution Vulnerability

{dataframe_to_markdown(vulnerability[["method_name", "method_role", "n_scored", "source_score_prevalence_corr", "low_prevalence_high_ranked_negative_rate", "rare_hla_high_ranked_negative_rate", "rare_hla_missed_positive_rate", "high_leakage_high_ranked_negative_rate", "external_high_ranked_negative_rate", "external_missed_positive_rate", "distribution_vulnerability_score", "primary_distribution_issues"]].head(30) if not vulnerability.empty else vulnerability, max_rows=30)}

## Shift-Linked Error Slices

{dataframe_to_markdown(risky_slices[["feature", "value", "train_n", "external_n", "train_positive_prevalence", "external_positive_prevalence", "distribution_shift_score", "max_error_pressure_score", "error_shift_link_score", "top_fp_methods", "top_fn_methods", "diagnosis"]] if not risky_slices.empty else risky_slices, max_rows=40)}

## Error Slice Matrix

{dataframe_to_markdown(top_matrix[["method_name", "feature", "value", "n_scored", "positive_prevalence", "high_ranked_negative_rate_among_neg", "missed_positive_rate_among_pos", "threshold_false_positive_rate_among_neg", "threshold_false_negative_rate_among_pos", "error_pressure_score"]].head(50) if not top_matrix.empty else top_matrix, max_rows=50)}

## Remediation Plan

{dataframe_to_markdown(remediation, max_rows=20)}

## Reviewer-Safe Claim Boundary

Allowed: distribution-aware benchmark audit, source/HLA/low-prevalence failure diagnosis, reliability weighting and abstention guidance.

Forbidden: new SOTA predictor, clinical vaccine selection, external validation proven, quantum advantage, or public pretrained tools as clean baselines without row-level training overlap audit.
"""
    (output_root / "CLEAN_NEOBENCH_DISTRIBUTION_ERROR_AUDIT.md").write_text(text.strip() + "\n")

    kr = f"""# CLEAN-NeoBench Distribution Error Audit KR

## 한 줄 결론

지금 싸움은 단순 AUPRC 1등 싸움이 아니다. 실제로는 **source prevalence shift**, **rare/low-support HLA**, **high leakage-risk region**, **low-prevalence TESLA-like setting**에서 누가 망가지는지 보는 싸움이다.

## method별 취약성

{dataframe_to_markdown(vulnerability[["method_name", "distribution_vulnerability_score", "primary_distribution_issues", "low_prevalence_high_ranked_negative_rate", "rare_hla_missed_positive_rate", "external_missed_positive_rate"]].head(25) if not vulnerability.empty else vulnerability, max_rows=25)}

## 분포 shift와 error가 같이 보이는 slice

{dataframe_to_markdown(risky_slices[["feature", "value", "train_n", "external_n", "distribution_shift_score", "max_error_pressure_score", "diagnosis"]].head(30) if not risky_slices.empty else risky_slices, max_rows=30)}

## 실용 판단

- aggregate leaderboard에서 이긴 method라도 source/HLA slice에서 무너지면 clean claim 금지.
- BAR-Neo-BMA는 이 취약성을 weight와 abstention으로 흡수하는 controller 위치가 맞다.
- QK는 feature/selector/fallback로만 사용하고 quantum advantage 표현은 금지.
- public pretrained tool은 training-corpus overlap audit 전까지 caveated comparator다.
- PAAD/THCA patient gate는 patient metadata 들어오기 전까지 demo/triage-only다.

## 다음 행동

1. public training corpus를 넣고 row-level overlap audit.
2. source-heldout/HLA-heldout challenge split을 고정.
3. rare-HLA/low-prevalence slice의 manual review queue를 우선 검토.
4. patient metadata가 있는 PAAD/THCA demo 1-2개를 붙여 실제 gate를 검증.
"""
    (output_root / "CLEAN_NEOBENCH_DISTRIBUTION_ERROR_AUDIT_KR.md").write_text(kr.strip() + "\n")


def main() -> None:
    args = parse_args()
    output_root = Path(args.output_root)
    ensure_dir(output_root)

    scored = read_tsv(output_root / "clean_neobench_failure_scored_candidates.tsv")
    shift = read_tsv(output_root / "clean_neobench_train_external_distribution_shift.tsv")
    if scored.empty:
        raise SystemExit("clean_neobench_failure_scored_candidates.tsv is required; run analyze_failure_modes.py first")

    matrix = build_error_slice_matrix(scored)
    link = build_shift_error_link(shift, matrix)
    vulnerability = build_method_vulnerability(scored, args.rare_hla_support, args.low_hla_support, args.low_prevalence_cutoff)
    remediation = build_remediation_plan(vulnerability, link)

    write_tsv(matrix, output_root / "clean_neobench_error_slice_matrix.tsv")
    write_tsv(link, output_root / "clean_neobench_train_distribution_error_link.tsv")
    write_tsv(vulnerability, output_root / "clean_neobench_method_distribution_vulnerability.tsv")
    write_tsv(remediation, output_root / "clean_neobench_distribution_remediation_plan.tsv")
    write_report(output_root, matrix, link, vulnerability, remediation)

    outputs = [
        "clean_neobench_error_slice_matrix.tsv",
        "clean_neobench_train_distribution_error_link.tsv",
        "clean_neobench_method_distribution_vulnerability.tsv",
        "clean_neobench_distribution_remediation_plan.tsv",
        "CLEAN_NEOBENCH_DISTRIBUTION_ERROR_AUDIT.md",
        "CLEAN_NEOBENCH_DISTRIBUTION_ERROR_AUDIT_KR.md",
    ]
    manifest_path = output_root / "run_manifest.json"
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
    else:
        manifest = {}
    manifest.setdefault("output_files", [])
    for name in outputs:
        if name not in manifest["output_files"]:
            manifest["output_files"].append(name)
    manifest.setdefault("summary", {})
    manifest["summary"].update(
        {
            "n_distribution_error_slice_rows": int(len(matrix)),
            "n_distribution_error_link_rows": int(len(link)),
            "n_method_distribution_vulnerability_rows": int(len(vulnerability)),
        }
    )
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    update_manifest(
        output_root,
        "distribution_error_audit",
        {
            "outputs": outputs,
            "n_error_slice_rows": int(len(matrix)),
            "n_shift_error_link_rows": int(len(link)),
            "n_method_vulnerability_rows": int(len(vulnerability)),
            "warnings": [
                "Distribution-error audit is diagnostic only; it does not establish a new SOTA predictor."
            ],
        },
    )
    print(
        "[clean-neobench-distribution-error] "
        f"slices={len(matrix)} shift_links={len(link)} methods={len(vulnerability)}"
    )


if __name__ == "__main__":
    main()
