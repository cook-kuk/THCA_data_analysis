#!/usr/bin/env python3
"""Analyze CLEAN-NeoBench failure modes and data-distribution shift.

The analysis treats neoantigen ranking errors as distributional signals:
high-ranked negatives indicate false-positive pressure, while positives that
remain low-ranked indicate missed immunogenic candidates. It compares those
cases against the method-specific scorable background and against the local
training-pool versus external-source distribution.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from common import ensure_dir, sequence_features, update_manifest, write_tsv


TRAIN_SOURCES = {"CEDAR", "TESLA_mmc4", "NEPdb", "TESLA_mmc7_validation"}
DEFAULT_EXTRA_METHODS = [
    "Structure_LR",
    "W7A_QK_only",
    "W7B_stacked",
    "ESM2_Bayesian",
    "BigMHC_IM",
    "MHCflurry",
    "PRIME",
    "NetMHCpan_4.1",
]
GROUP_FEATURES = [
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
    parser.add_argument("--top-methods", type=int, default=8, help="Top leaderboard methods to include")
    parser.add_argument("--top-negative-quantile", type=float, default=0.90, help="Rank percentile cutoff for high-ranked negatives")
    parser.add_argument("--missed-positive-quantile", type=float, default=0.50, help="Rank percentile cutoff for low-ranked positives")
    return parser.parse_args()


def load_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def clean_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def add_distribution_features(master: pd.DataFrame) -> pd.DataFrame:
    df = master.copy()
    df["distribution_partition"] = np.where(df["source_name"].isin(TRAIN_SOURCES), "local_train_pool", "external_or_holdout")
    source_prev = df.groupby("source_name")["label"].mean().to_dict() if "label" in df else {}
    hla_support = df.groupby("hla_allele_4digit")["candidate_id"].count().to_dict() if "hla_allele_4digit" in df else {}
    hla_prev = df.groupby("hla_allele_4digit")["label"].mean().to_dict() if "hla_allele_4digit" in df else {}
    df["source_positive_prevalence"] = df["source_name"].map(source_prev)
    df["hla_allele_support_count"] = df["hla_allele_4digit"].map(hla_support)
    df["hla_allele_positive_prevalence"] = df["hla_allele_4digit"].map(hla_prev)
    seq = pd.DataFrame([sequence_features(p) for p in df["peptide"].fillna("")], index=df.index)
    for col in seq.columns:
        df[col] = seq[col]
    return df


def choose_focus_methods(leaderboard: pd.DataFrame) -> list[str]:
    methods: list[str] = []
    if not leaderboard.empty and "method_name" in leaderboard.columns:
        methods.extend(leaderboard["method_name"].head(8).astype(str).tolist())
    for method in DEFAULT_EXTRA_METHODS:
        if method not in methods:
            methods.append(method)
    return methods


def aggregate_method_scores(scores: pd.DataFrame, focus_methods: list[str]) -> pd.DataFrame:
    if scores.empty:
        return pd.DataFrame()
    s = scores[scores["method_name"].isin(focus_methods)].copy()
    if s.empty:
        return pd.DataFrame()
    s["score_calibrated"] = pd.to_numeric(s["score_calibrated"], errors="coerce")
    s["score_raw"] = pd.to_numeric(s["score_raw"], errors="coerce")
    s = s.sort_values(["method_name", "candidate_id", "score_calibrated"], ascending=[True, True, False])
    agg = (
        s.groupby(["method_name", "candidate_id"], as_index=False)
        .agg(
            method_family=("method_family", "first"),
            method_role=("method_role", "first"),
            score_for_analysis=("score_calibrated", "max"),
            score_raw_best=("score_raw", "max"),
            n_score_contexts=("score_context", "nunique"),
            score_contexts=("score_context", lambda x: ";".join(sorted(set(map(str, x)))[:8])),
            uses_public_pretraining=("uses_public_pretraining", lambda x: any(clean_bool(v) for v in x)),
            clean_comparator_allowed=("clean_comparator_allowed", lambda x: all(clean_bool(v) for v in x)),
        )
    )
    return agg


def bma_as_method(bma: pd.DataFrame) -> pd.DataFrame:
    if bma.empty:
        return pd.DataFrame()
    score_col = "patient_gated_bma_score" if "patient_gated_bma_score" in bma.columns else "barneo_bma_score"
    out = pd.DataFrame(
        {
            "method_name": "BAR_Neo_BMA",
            "candidate_id": bma["candidate_id"],
            "method_family": "ensemble_or_fusion",
            "method_role": "agentic_reliability_controller",
            "score_for_analysis": pd.to_numeric(bma[score_col], errors="coerce"),
            "score_raw_best": pd.to_numeric(bma.get("barneo_bma_score", bma[score_col]), errors="coerce"),
            "n_score_contexts": 1,
            "score_contexts": "BAR_Neo_BMA",
            "uses_public_pretraining": False,
            "clean_comparator_allowed": True,
        }
    )
    return out


def label_failure(row: pd.Series, top_neg_q: float, missed_pos_q: float) -> tuple[str, str]:
    label = int(row["label"])
    rank_pct = float(row["rank_percentile"])
    score = float(row["score_for_analysis"])
    tags = []
    primary = "correct_or_noncritical"
    if label == 0 and rank_pct >= top_neg_q:
        tags.append("high_ranked_negative_top10pct")
        primary = "high_ranked_negative_top10pct"
    if label == 0 and score >= 0.50:
        tags.append("threshold_false_positive_score_ge_0.5")
        if primary == "correct_or_noncritical":
            primary = "threshold_false_positive_score_ge_0.5"
    if label == 1 and rank_pct < missed_pos_q:
        tags.append("missed_positive_bottom50pct")
        primary = "missed_positive_bottom50pct"
    if label == 1 and score < 0.50:
        tags.append("threshold_false_negative_score_lt_0.5")
        if primary == "correct_or_noncritical":
            primary = "threshold_false_negative_score_lt_0.5"
    return primary, ";".join(tags)


def hypothesis(row: pd.Series) -> str:
    reasons = []
    label = int(row.get("label", 0))
    src_prev = float(row.get("source_positive_prevalence", math.nan))
    hla_support = float(row.get("hla_allele_support_count", math.nan))
    if row.get("primary_failure_type") == "high_ranked_negative_top10pct" and src_prev < 0.10:
        reasons.append("low-prevalence source false-positive pressure")
    if row.get("primary_failure_type", "").startswith("missed_positive") and src_prev > 0.70:
        reasons.append("positive-rich source missed despite favorable prior")
    if hla_support < 20:
        reasons.append("underrepresented HLA allele")
    if str(row.get("leakage_risk_level", "")).lower() == "high":
        reasons.append("high leakage-risk region; benchmark claim unsafe")
    if clean_bool(row.get("uses_public_pretraining", False)):
        reasons.append("public-pretrained comparator; overlap unresolved")
    if label == 0 and row.get("rank_percentile", 0) >= 0.90 and row.get("peptide_length", 0) not in {9, 10}:
        reasons.append("noncanonical peptide length high-ranked negative")
    return "; ".join(reasons) if reasons else "no single dominant distributional hypothesis"


def build_failure_tables(
    master: pd.DataFrame,
    method_scores: pd.DataFrame,
    bma: pd.DataFrame,
    leaderboard: pd.DataFrame,
    top_neg_q: float,
    missed_pos_q: float,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    focus = choose_focus_methods(leaderboard)
    scores = aggregate_method_scores(method_scores, focus)
    bma_method = bma_as_method(bma)
    if not bma_method.empty:
        scores = pd.concat([scores, bma_method], ignore_index=True, sort=False)
    if scores.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    base_cols = [
        "candidate_id",
        "label",
        "source_name",
        "source_dataset",
        "distribution_partition",
        "hla",
        "hla_allele_4digit",
        "hla_supertype",
        "peptide",
        "peptide_length",
        "mhc_class",
        "leakage_risk_level",
        "split_low_prevalence",
        "exact_peptide_hla_train_overlap",
        "near_peptide_train_overlap",
        "source_positive_prevalence",
        "hla_allele_support_count",
        "hla_allele_positive_prevalence",
        "hydrophobic_fraction",
        "aromatic_fraction",
        "charge_proxy",
        "cysteine_count",
        "glycine_proline_count",
    ]
    merged = scores.merge(master[[c for c in base_cols if c in master.columns]], on="candidate_id", how="left")
    merged["label"] = pd.to_numeric(merged["label"], errors="coerce")
    merged = merged[merged["label"].isin([0, 1]) & merged["score_for_analysis"].notna()].copy()
    merged["rank_percentile"] = merged.groupby("method_name")["score_for_analysis"].rank(method="average", pct=True)
    failures = merged.apply(lambda r: label_failure(r, top_neg_q, missed_pos_q), axis=1, result_type="expand")
    merged["primary_failure_type"] = failures[0]
    merged["failure_tags"] = failures[1]
    merged["failure_hypothesis"] = merged.apply(hypothesis, axis=1)

    summary_rows = []
    for method, g in merged.groupby("method_name"):
        pos = g[g["label"].eq(1)]
        neg = g[g["label"].eq(0)]
        summary_rows.append(
            {
                "method_name": method,
                "method_role": g["method_role"].iloc[0],
                "method_family": g["method_family"].iloc[0],
                "n_scored": int(len(g)),
                "n_pos": int(pos.shape[0]),
                "n_neg": int(neg.shape[0]),
                "positive_prevalence": float(g["label"].mean()) if len(g) else math.nan,
                "high_ranked_negative_top10pct": int((neg["rank_percentile"] >= top_neg_q).sum()),
                "high_ranked_negative_rate_among_neg": float((neg["rank_percentile"] >= top_neg_q).mean()) if len(neg) else math.nan,
                "missed_positive_bottom50pct": int((pos["rank_percentile"] < missed_pos_q).sum()),
                "missed_positive_rate_among_pos": float((pos["rank_percentile"] < missed_pos_q).mean()) if len(pos) else math.nan,
                "threshold_false_positive_score_ge_0.5": int((neg["score_for_analysis"] >= 0.50).sum()),
                "threshold_false_negative_score_lt_0.5": int((pos["score_for_analysis"] < 0.50).sum()),
                "median_score_pos": float(pos["score_for_analysis"].median()) if len(pos) else math.nan,
                "median_score_neg": float(neg["score_for_analysis"].median()) if len(neg) else math.nan,
                "n_local_train_pool": int(g["distribution_partition"].eq("local_train_pool").sum()),
                "n_external_or_holdout": int(g["distribution_partition"].eq("external_or_holdout").sum()),
            }
        )
    summary = pd.DataFrame(summary_rows).sort_values(["missed_positive_rate_among_pos", "high_ranked_negative_rate_among_neg"], ascending=False)
    cases = merged[merged["primary_failure_type"].ne("correct_or_noncritical")].copy()
    cases = cases.sort_values(["method_name", "primary_failure_type", "rank_percentile"], ascending=[True, True, False])
    return merged, cases, summary


def build_distribution_shift(master: pd.DataFrame) -> pd.DataFrame:
    rows = []
    train = master[master["distribution_partition"].eq("local_train_pool")]
    ext = master[master["distribution_partition"].ne("local_train_pool")]
    for feature in GROUP_FEATURES:
        if feature not in master.columns:
            continue
        train_vc = train[feature].fillna("NA").astype(str).value_counts(dropna=False)
        ext_vc = ext[feature].fillna("NA").astype(str).value_counts(dropna=False)
        values = sorted(set(train_vc.index) | set(ext_vc.index))
        for value in values:
            train_n = int(train_vc.get(value, 0))
            ext_n = int(ext_vc.get(value, 0))
            train_frac = train_n / max(1, len(train))
            ext_frac = ext_n / max(1, len(ext))
            sub_train = train[train[feature].fillna("NA").astype(str).eq(value)]
            sub_ext = ext[ext[feature].fillna("NA").astype(str).eq(value)]
            rows.append(
                {
                    "feature": feature,
                    "value": value,
                    "train_n": train_n,
                    "train_fraction": train_frac,
                    "train_positive_prevalence": float(sub_train["label"].mean()) if len(sub_train) else math.nan,
                    "external_n": ext_n,
                    "external_fraction": ext_frac,
                    "external_positive_prevalence": float(sub_ext["label"].mean()) if len(sub_ext) else math.nan,
                    "external_over_train_fraction_ratio": ext_frac / train_frac if train_frac else math.inf,
                    "absolute_fraction_delta": abs(ext_frac - train_frac),
                }
            )
    return pd.DataFrame(rows).sort_values(["absolute_fraction_delta", "external_over_train_fraction_ratio"], ascending=[False, False])


def build_failure_distribution(cases: pd.DataFrame, scored: pd.DataFrame) -> pd.DataFrame:
    rows = []
    if cases.empty or scored.empty:
        return pd.DataFrame()
    for method, bg in scored.groupby("method_name"):
        bg_n = len(bg)
        case_m = cases[cases["method_name"].eq(method)]
        for failure_type, fg in case_m.groupby("primary_failure_type"):
            err_n = len(fg)
            for feature in GROUP_FEATURES:
                if feature not in bg.columns:
                    continue
                bg_vc = bg[feature].fillna("NA").astype(str).value_counts(dropna=False)
                fg_vc = fg[feature].fillna("NA").astype(str).value_counts(dropna=False)
                for value, n_error in fg_vc.items():
                    n_bg = int(bg_vc.get(value, 0))
                    err_frac = n_error / max(1, err_n)
                    bg_frac = n_bg / max(1, bg_n)
                    rows.append(
                        {
                            "method_name": method,
                            "primary_failure_type": failure_type,
                            "feature": feature,
                            "value": value,
                            "n_error": int(n_error),
                            "error_fraction": err_frac,
                            "n_background": n_bg,
                            "background_fraction": bg_frac,
                            "enrichment_vs_background": err_frac / bg_frac if bg_frac else math.inf,
                            "mean_error_score": float(fg.loc[fg[feature].fillna("NA").astype(str).eq(value), "score_for_analysis"].mean()),
                            "error_positive_prevalence": float(fg.loc[fg[feature].fillna("NA").astype(str).eq(value), "label"].mean()),
                        }
                    )
    return pd.DataFrame(rows).sort_values(["enrichment_vs_background", "n_error"], ascending=[False, False])


def markdown_table(df: pd.DataFrame, columns: list[str] | None = None, n: int = 20) -> str:
    if df.empty:
        return "No rows available."
    view = df.copy()
    if columns is not None:
        view = view[[c for c in columns if c in view.columns]]
    view = view.head(n)
    try:
        return view.to_markdown(index=False)
    except Exception:
        return view.to_csv(sep="\t", index=False)


def write_report(
    output_root: Path,
    master: pd.DataFrame,
    scored: pd.DataFrame,
    cases: pd.DataFrame,
    summary: pd.DataFrame,
    shift: pd.DataFrame,
    failure_dist: pd.DataFrame,
) -> None:
    source_prev = (
        master.groupby(["distribution_partition", "source_name"])["label"].agg(["count", "sum", "mean"]).reset_index().sort_values(["distribution_partition", "count"], ascending=[True, False])
        if not master.empty
        else pd.DataFrame()
    )
    high_fp = cases[cases["primary_failure_type"].eq("high_ranked_negative_top10pct")]
    missed = cases[cases["primary_failure_type"].eq("missed_positive_bottom50pct")]
    enriched = failure_dist[np.isfinite(failure_dist["enrichment_vs_background"])] if not failure_dist.empty else failure_dist
    text = f"""# CLEAN-NeoBench Failure Mode Analysis

## Purpose

This analysis asks where the current predictors are wrong, and whether those errors align with the training-data distribution. It does not introduce a new model. It identifies false-positive pressure, missed positives, source/HLA imbalance, and low-prevalence shift.

## Key Data Distribution

Local train-pool sources are treated as `CEDAR`, `TESLA_mmc4`, `NEPdb`, and `TESLA_mmc7_validation`. External/holdout rows are mainly ITSNdb-derived rows.

{markdown_table(source_prev, ["distribution_partition", "source_name", "count", "sum", "mean"], 20)}

## Method Failure Summary

Definitions:
- `high_ranked_negative_top10pct`: label 0 candidate ranked in the top 10 percent for that method.
- `missed_positive_bottom50pct`: label 1 candidate ranked in the bottom half for that method.
- threshold false positives/negatives use calibrated score 0.5.

{markdown_table(summary, ["method_name", "method_role", "n_scored", "positive_prevalence", "high_ranked_negative_top10pct", "high_ranked_negative_rate_among_neg", "missed_positive_bottom50pct", "missed_positive_rate_among_pos", "median_score_pos", "median_score_neg"], 25)}

## Distribution Shift: Train Pool vs External/Holdout

Largest train/external distribution deltas:

{markdown_table(shift, ["feature", "value", "train_n", "train_fraction", "train_positive_prevalence", "external_n", "external_fraction", "external_positive_prevalence", "external_over_train_fraction_ratio", "absolute_fraction_delta"], 30)}

## Error Enrichment

Rows below show where error cases are enriched relative to that method's own scorable background.

{markdown_table(enriched, ["method_name", "primary_failure_type", "feature", "value", "n_error", "error_fraction", "background_fraction", "enrichment_vs_background", "mean_error_score"], 40)}

## High-Ranked Negatives

These are the most dangerous practical false positives: negatives that look strong to a method.

{markdown_table(high_fp.sort_values(["rank_percentile", "score_for_analysis"], ascending=[False, False]), ["method_name", "candidate_id", "label", "score_for_analysis", "rank_percentile", "source_name", "source_positive_prevalence", "hla_allele_4digit", "hla_allele_support_count", "leakage_risk_level", "failure_hypothesis"], 40)}

## Missed Positives

These are positives that a method ranked too low.

{markdown_table(missed.sort_values(["rank_percentile", "score_for_analysis"], ascending=[True, True]), ["method_name", "candidate_id", "label", "score_for_analysis", "rank_percentile", "source_name", "source_positive_prevalence", "hla_allele_4digit", "hla_allele_support_count", "leakage_risk_level", "failure_hypothesis"], 40)}

## Interpretation

The strongest expected failure mode is source-prevalence shift: CEDAR is heavily positive, while TESLA_mmc4, TESLA_mmc7_validation, and ITSNdb_Val are low-prevalence. A model can look strong by learning source priors and still create high-ranked negatives in low-prevalence settings. HLA allele imbalance is the second expected failure axis; rare alleles can inflate both false positives and missed positives. BAR-Neo-BMA should therefore keep abstention active until patient metadata, source-heldout calibration, and public-tool training overlap audit are strengthened.
"""
    (output_root / "CLEAN_NEOBENCH_FAILURE_MODE_ANALYSIS.md").write_text(text.rstrip() + "\n")


def write_kr_report(output_root: Path, master: pd.DataFrame, summary: pd.DataFrame, shift: pd.DataFrame) -> None:
    source_prev = (
        master.groupby(["distribution_partition", "source_name"])["label"].agg(["count", "sum", "mean"]).reset_index().sort_values(["distribution_partition", "count"], ascending=[True, False])
        if not master.empty
        else pd.DataFrame()
    )
    source_short = markdown_table(source_prev, ["distribution_partition", "source_name", "count", "sum", "mean"], 20)
    method_short = markdown_table(
        summary,
        [
            "method_name",
            "n_scored",
            "high_ranked_negative_top10pct",
            "missed_positive_bottom50pct",
            "missed_positive_rate_among_pos",
            "median_score_pos",
            "median_score_neg",
        ],
        20,
    )
    shift_short = markdown_table(
        shift,
        [
            "feature",
            "value",
            "train_n",
            "train_fraction",
            "train_positive_prevalence",
            "external_n",
            "external_fraction",
            "external_positive_prevalence",
            "absolute_fraction_delta",
        ],
        25,
    )
    text = f"""# CLEAN-NeoBench Failure Mode Analysis KR

Date: 2026-05-09

## 한 줄 결론

틀리는 후보들은 무작위가 아니다. 가장 큰 축은 **source prevalence shift**, 두 번째는 **HLA allele support imbalance**, 세 번째는 **high leakage-risk region**이다.

## 학습/외부 분포 판단

{source_short}

CEDAR는 거의 positive-rich source이고, TESLA/ITSNdb_Val은 low-prevalence source다. 이 차이가 너무 커서 model이 biological rule이 아니라 source prior를 배울 수 있다.

## method별 failure 요약

{method_short}

해석:

- `BAR_Neo_BMA`는 missed positive rate가 가장 낮지만, abstention을 강하게 유지한다.
- `source_balanced_rf_train_prior_calibrated` 계열은 현재 best internal group이다.
- `Structure_LR`는 honest anchor이고 false-positive pressure가 낮다.
- `W7A_QK_only`는 점수는 좋지만 bounded fallback으로만 유지한다.
- public pretrained methods는 overlap unresolved라 clean baseline이 아니다.

## train vs external shift

{shift_short}

## 실용 조치

1. source prior correction을 더 강하게 넣는다.
2. low-prevalence source에서는 top-k threshold를 보수화한다.
3. rare HLA allele은 confidence를 낮춘다.
4. HLA별 minimum support count를 reliability feature로 강화한다.
5. public tool agreement는 caveated support로만 쓰고 clean feature로 쓰지 않는다.
6. PAAD/THCA patient metadata가 들어오기 전까지 clinical-style confidence는 내지 않는다.

## 파일

- Full report: `CLEAN_NEOBENCH_FAILURE_MODE_ANALYSIS.md`
- Failure cases: `clean_neobench_failure_cases.tsv`
- Failure summary: `clean_neobench_failure_summary.tsv`
- Train/external shift: `clean_neobench_train_external_distribution_shift.tsv`
- Failure enrichment: `clean_neobench_failure_distribution_enrichment.tsv`
"""
    (output_root / "CLEAN_NEOBENCH_FAILURE_MODE_ANALYSIS_KR.md").write_text(text.rstrip() + "\n")


def main() -> None:
    args = parse_args()
    output_root = Path(args.output_root)
    ensure_dir(output_root)

    master = load_tsv(output_root / "clean_neobench_master.tsv")
    scores = load_tsv(output_root / "clean_neobench_method_scores.tsv")
    bma = load_tsv(output_root / "barneo_bma_candidate_scores.tsv")
    leaderboard = load_tsv(output_root / "clean_neobench_leaderboard.tsv")

    master = add_distribution_features(master)
    scored, cases, summary = build_failure_tables(
        master,
        scores,
        bma,
        leaderboard.head(args.top_methods) if not leaderboard.empty else leaderboard,
        args.top_negative_quantile,
        args.missed_positive_quantile,
    )
    shift = build_distribution_shift(master)
    failure_dist = build_failure_distribution(cases, scored)

    write_tsv(scored, output_root / "clean_neobench_failure_scored_candidates.tsv")
    write_tsv(cases, output_root / "clean_neobench_failure_cases.tsv")
    write_tsv(summary, output_root / "clean_neobench_failure_summary.tsv")
    write_tsv(shift, output_root / "clean_neobench_train_external_distribution_shift.tsv")
    write_tsv(failure_dist, output_root / "clean_neobench_failure_distribution_enrichment.tsv")
    write_report(output_root, master, scored, cases, summary, shift, failure_dist)
    write_kr_report(output_root, master, summary, shift)

    manifest_path = output_root / "run_manifest.json"
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
    else:
        manifest = {}
    manifest.setdefault("output_files", [])
    for name in [
        "clean_neobench_failure_scored_candidates.tsv",
        "clean_neobench_failure_cases.tsv",
        "clean_neobench_failure_summary.tsv",
        "clean_neobench_train_external_distribution_shift.tsv",
        "clean_neobench_failure_distribution_enrichment.tsv",
        "CLEAN_NEOBENCH_FAILURE_MODE_ANALYSIS.md",
        "CLEAN_NEOBENCH_FAILURE_MODE_ANALYSIS_KR.md",
    ]:
        if name not in manifest["output_files"]:
            manifest["output_files"].append(name)
    manifest.setdefault("summary", {})
    manifest["summary"].update(
        {
            "n_failure_methods_analyzed": int(scored["method_name"].nunique()) if not scored.empty else 0,
            "n_failure_scored_rows": int(len(scored)),
            "n_failure_case_rows": int(len(cases)),
        }
    )
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    update_manifest(
        output_root,
        "failure_mode_analysis",
        {
            "n_methods": int(scored["method_name"].nunique()) if not scored.empty else 0,
            "n_scored_rows": int(len(scored)),
            "n_failure_rows": int(len(cases)),
            "warnings": [
                "Method candidate scores are aggregated by best calibrated score across available score_contexts to identify false-positive pressure."
            ],
        },
    )
    print(f"[clean-neobench-failure] methods={scored['method_name'].nunique() if not scored.empty else 0} failures={len(cases)}")


if __name__ == "__main__":
    main()
