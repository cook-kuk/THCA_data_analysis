#!/usr/bin/env python3
"""Build BAR-Neo-X interpretability and claim-safe reranking outputs.

This layer does not replace BAR-Neo. It decomposes each candidate score into
auditable components and adds a claim-safe review score that deliberately
downranks leakage-heavy rows.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from common import average_precision_score, ensure_dir, roc_auc_score, update_manifest, write_tsv


TOP_KS = [10, 20, 50, 100]
SCORE_COLUMNS = [
    "barneo_score",
    "barneo_bma_score",
    "patient_gated_score",
    "patient_gated_bma_score",
    "barneo_x_discovery_score",
    "barneo_x_claim_safe_score",
]
OUTPUT_FILES = [
    "barneo_x_candidate_scores.tsv",
    "barneo_x_candidate_explanations.tsv",
    "barneo_x_component_attributions.tsv",
    "barneo_x_ablation_audit.tsv",
    "barneo_x_metric_audit.tsv",
    "barneo_x_topk_safety_audit.tsv",
    "BAR_NEO_X_INTERPRETABILITY_BOOST_REPORT.md",
    "BAR_NEO_X_INTERPRETABILITY_BOOST_SUMMARY.json",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=".", help="Repository root")
    parser.add_argument("--output-root", required=True, help="Output directory")
    return parser.parse_args()


def load(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def num(value: Any, default: float = 0.0) -> float:
    try:
        v = float(value)
        return v if np.isfinite(v) else default
    except Exception:
        return default


def truthy(value: Any) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def empty(value: Any) -> bool:
    s = str(value).strip().lower()
    return s in {"", "na", "nan", "none", "null", "n/a"}


def best_internal_scores(scores: pd.DataFrame) -> pd.DataFrame:
    if scores.empty:
        return pd.DataFrame(columns=["candidate_id", "best_clean_internal_score", "best_anchor_score", "best_public_score"])
    internal = scores[scores["method_role"].isin(["anchor", "internal_candidate", "bounded_fallback"])]
    anchor = scores[scores["method_role"].eq("anchor")]
    public = scores[scores["method_role"].eq("caveated_public_comparator")]
    out = pd.DataFrame({"candidate_id": scores["candidate_id"].drop_duplicates()})
    for name, df in [
        ("best_clean_internal_score", internal),
        ("best_anchor_score", anchor),
        ("best_public_score", public),
    ]:
        vals = df.groupby("candidate_id")["score_calibrated"].max().rename(name)
        out = out.merge(vals, on="candidate_id", how="left")
    return out


def candidate_table(output_root: Path) -> pd.DataFrame:
    master = load(output_root / "clean_neobench_master.tsv")
    scores = load(output_root / "clean_neobench_method_scores.tsv")
    barneo = load(output_root / "barneo_candidate_scores.tsv")
    bma = load(output_root / "barneo_bma_candidate_scores.tsv")
    internal = best_internal_scores(scores)

    cols = [
        "candidate_id",
        "source_name",
        "study_id",
        "patient_id",
        "peptide",
        "hla_allele_4digit",
        "label",
        "mhc_class",
        "cancer_type",
        "disease_context",
        "leakage_risk_level",
        "exact_peptide_train_overlap",
        "exact_peptide_hla_train_overlap",
        "near_peptide_train_overlap",
        "study_train_overlap",
        "patient_train_overlap",
        "public_tool_training_overlap_any",
        "split_low_prevalence",
        "split_korean_hla_focus",
    ]
    df = master[[c for c in cols if c in master.columns]].copy()
    keep_bar = [
        "candidate_id",
        "barneo_score",
        "patient_gated_score",
        "confidence_score",
        "confidence_bin",
        "abstain",
        "abstention_reason_primary",
        "disease_gate",
        "presentation_gate",
        "antigen_gate",
        "immune_context_gate",
        "model_confidence_gate",
    ]
    keep_bma = [
        "candidate_id",
        "barneo_bma_score",
        "patient_gated_bma_score",
        "posterior_uncertainty",
        "expert_disagreement_range",
        "bma_confidence_score",
        "bma_confidence_bin",
        "bma_abstain",
        "bma_abstention_reason_primary",
        "selected_method_count",
        "selected_methods",
    ]
    df = df.merge(barneo[[c for c in keep_bar if c in barneo.columns]], on="candidate_id", how="left")
    df = df.merge(bma[[c for c in keep_bma if c in bma.columns]], on="candidate_id", how="left")
    df = df.merge(internal, on="candidate_id", how="left")
    return df


def contribution_dict(row: pd.Series) -> dict[str, float]:
    barneo = num(row.get("barneo_score"), 0.5)
    bma = num(row.get("barneo_bma_score"), barneo)
    confidence = num(row.get("confidence_score"), 0.5)
    internal = num(row.get("best_clean_internal_score"), barneo)
    anchor = num(row.get("best_anchor_score"), 0.5)
    leak = str(row.get("leakage_risk_level", "")).lower()
    method_count = num(row.get("selected_method_count"), 0)
    disagreement = num(row.get("expert_disagreement_range"), 0)
    uncertainty = num(row.get("posterior_uncertainty"), 0)
    gates = [
        num(row.get("disease_gate"), 0.8),
        num(row.get("presentation_gate"), 0.8),
        num(row.get("antigen_gate"), 0.8),
        num(row.get("immune_context_gate"), 0.8),
        num(row.get("model_confidence_gate"), 0.8),
    ]
    patient_gate = float(np.prod(np.clip(gates, 0, 1)))

    pos = {
        "BAR-Neo model evidence": 0.42 * barneo,
        "BMA expert consensus": 0.24 * bma,
        "clean internal predictor support": 0.16 * internal,
        "anchor support": 0.06 * anchor,
        "confidence support": 0.12 * confidence,
    }
    neg = {
        "posterior disagreement penalty": -0.16 * min(1.0, disagreement),
        "posterior uncertainty penalty": -0.12 * min(1.0, uncertainty),
        "patient gate incompleteness penalty": -0.15 * (1.0 - patient_gate),
    }
    if leak == "high":
        neg["high leakage risk penalty"] = -0.24
    elif leak == "medium":
        neg["medium leakage risk penalty"] = -0.10
    if truthy(row.get("exact_peptide_hla_train_overlap")):
        neg["exact peptide-HLA overlap penalty"] = -0.18
    elif truthy(row.get("exact_peptide_train_overlap")):
        neg["exact peptide overlap penalty"] = -0.12
    if truthy(row.get("near_peptide_train_overlap")):
        neg["near peptide overlap penalty"] = -0.06
    if truthy(row.get("study_train_overlap")):
        neg["study overlap penalty"] = -0.05
    if truthy(row.get("patient_train_overlap")):
        neg["patient overlap penalty"] = -0.10
    if truthy(row.get("public_tool_training_overlap_any")):
        neg["public pretrained overlap caveat penalty"] = -0.08
    if str(row.get("split_low_prevalence", "")) == "low_prevalence":
        neg["low-prevalence source penalty"] = -0.05
    if method_count < 3:
        neg["sparse expert support penalty"] = -0.06
    return {**pos, **neg}


def add_scores(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in df.iterrows():
        contrib = contribution_dict(row)
        positive = sum(v for v in contrib.values() if v > 0)
        penalty = -sum(v for v in contrib.values() if v < 0)
        discovery = np.clip(positive - 0.55 * penalty, 0, 1)
        leak = str(row.get("leakage_risk_level", "")).lower()
        claim_gate = 0.50 if leak == "high" else 0.75 if leak == "medium" else 1.0
        if truthy(row.get("exact_peptide_hla_train_overlap")):
            claim_gate *= 0.80
        if truthy(row.get("public_tool_training_overlap_any")):
            claim_gate *= 0.85
        metadata_cap = 0.75 if empty(row.get("cancer_type")) or empty(row.get("disease_context")) else 1.0
        claim_safe = np.clip(discovery * claim_gate * metadata_cap, 0, 1)
        sorted_pos = sorted(((k, v) for k, v in contrib.items() if v > 0), key=lambda x: x[1], reverse=True)
        sorted_neg = sorted(((k, v) for k, v in contrib.items() if v < 0), key=lambda x: x[1])
        out = row.to_dict()
        out.update(
            {
                "barneo_x_positive_component_sum": positive,
                "barneo_x_penalty_sum": penalty,
                "barneo_x_discovery_score": float(discovery),
                "barneo_x_claim_safe_score": float(claim_safe),
                "barneo_x_claim_gate": float(claim_gate),
                "barneo_x_metadata_cap": float(metadata_cap),
                "barneo_x_top_positive_factors": "; ".join(f"{k}={v:.3f}" for k, v in sorted_pos[:4]),
                "barneo_x_top_negative_factors": "; ".join(f"{k}={v:.3f}" for k, v in sorted_neg[:5]),
                "barneo_x_primary_action": primary_action(row, claim_safe),
            }
        )
        rows.append(out)
    out_df = pd.DataFrame(rows)
    out_df["barneo_x_discovery_rank"] = out_df["barneo_x_discovery_score"].rank(ascending=False, method="first").astype(int)
    out_df["barneo_x_claim_safe_rank"] = out_df["barneo_x_claim_safe_score"].rank(ascending=False, method="first").astype(int)
    return out_df.sort_values("barneo_x_claim_safe_rank")


def primary_action(row: pd.Series, claim_safe: float) -> str:
    if str(row.get("leakage_risk_level", "")).lower() == "high":
        return "manual_review_only_high_leakage"
    if claim_safe >= 0.50:
        return "priority_review_candidate"
    if claim_safe >= 0.30:
        return "secondary_review_candidate"
    return "abstain_or_metadata_required"


def metric_rows(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    y = pd.to_numeric(df.get("label"), errors="coerce")
    valid = y.isin([0, 1])
    for score_col in SCORE_COLUMNS:
        if score_col not in df.columns:
            continue
        s = pd.to_numeric(df[score_col], errors="coerce")
        mask = valid & s.notna()
        if mask.sum() < 10 or y[mask].nunique() < 2:
            continue
        yy = y[mask].astype(int).to_numpy()
        ss = s[mask].astype(float).to_numpy()
        order = np.argsort(-ss)
        row = {
            "score": score_col,
            "n": int(mask.sum()),
            "apparent_auprc": float(average_precision_score(yy, ss)) if average_precision_score else np.nan,
            "apparent_auroc": float(roc_auc_score(yy, ss)) if roc_auc_score else np.nan,
            "top10_precision": float(yy[order[:10]].mean()) if len(order) >= 10 else np.nan,
            "top20_precision": float(yy[order[:20]].mean()) if len(order) >= 20 else np.nan,
            "top10_high_leakage_fraction": float(
                df.loc[mask].iloc[order[:10]]["leakage_risk_level"].astype(str).str.lower().eq("high").mean()
            )
            if len(order) >= 10 and "leakage_risk_level" in df.columns
            else np.nan,
        }
        rows.append(row)
    return pd.DataFrame(rows)


def score_metric_row(df: pd.DataFrame, score_col: str) -> dict[str, Any] | None:
    y = pd.to_numeric(df.get("label"), errors="coerce")
    s = pd.to_numeric(df.get(score_col), errors="coerce")
    valid = y.isin([0, 1]) & s.notna()
    if valid.sum() < 10 or y[valid].nunique() < 2:
        return None
    yy = y[valid].astype(int).to_numpy()
    ss = s[valid].astype(float).to_numpy()
    ranked = df.loc[valid].iloc[np.argsort(-ss)]
    top10 = ranked.head(10)
    top20 = ranked.head(20)
    return {
        "score": score_col,
        "n": int(valid.sum()),
        "apparent_auprc": float(average_precision_score(yy, ss)) if average_precision_score else np.nan,
        "apparent_auroc": float(roc_auc_score(yy, ss)) if roc_auc_score else np.nan,
        "top10_precision": float(pd.to_numeric(top10.get("label"), errors="coerce").mean()) if len(top10) else np.nan,
        "top20_precision": float(pd.to_numeric(top20.get("label"), errors="coerce").mean()) if len(top20) else np.nan,
        "top10_high_leakage_fraction": float(top10["leakage_risk_level"].astype(str).str.lower().eq("high").mean()) if len(top10) else np.nan,
        "top20_high_leakage_fraction": float(top20["leakage_risk_level"].astype(str).str.lower().eq("high").mean()) if len(top20) else np.nan,
    }


def component_attribution_rows(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in df.iterrows():
        contrib = contribution_dict(row)
        total_positive = sum(v for v in contrib.values() if v > 0)
        total_penalty = -sum(v for v in contrib.values() if v < 0)
        for name, value in contrib.items():
            if abs(value) < 1e-12:
                continue
            sign = "penalty" if "penalty" in name or value < 0 else "positive"
            denominator = total_positive if sign == "positive" else total_penalty
            rows.append(
                {
                    "candidate_id": row.get("candidate_id"),
                    "peptide": row.get("peptide"),
                    "hla_allele_4digit": row.get("hla_allele_4digit"),
                    "label": row.get("label"),
                    "leakage_risk_level": row.get("leakage_risk_level"),
                    "barneo_x_claim_safe_rank": row.get("barneo_x_claim_safe_rank"),
                    "component_name": name,
                    "component_sign": sign,
                    "component_value": float(value),
                    "component_abs_value": float(abs(value)),
                    "component_fraction_within_sign": float(abs(value) / denominator) if denominator else 0.0,
                }
            )
    return pd.DataFrame(rows).sort_values(["barneo_x_claim_safe_rank", "component_sign", "component_abs_value"], ascending=[True, True, False])


def ablation_rows(df: pd.DataFrame) -> pd.DataFrame:
    audit = df.copy()
    audit["barneo_x_no_penalty_score"] = np.clip(
        pd.to_numeric(audit["barneo_x_positive_component_sum"], errors="coerce")
        * pd.to_numeric(audit["barneo_x_claim_gate"], errors="coerce")
        * pd.to_numeric(audit["barneo_x_metadata_cap"], errors="coerce"),
        0,
        1,
    )
    audit["barneo_x_no_leakage_gate_score"] = np.clip(
        pd.to_numeric(audit["barneo_x_discovery_score"], errors="coerce")
        * pd.to_numeric(audit["barneo_x_metadata_cap"], errors="coerce"),
        0,
        1,
    )
    audit["barneo_x_no_metadata_cap_score"] = np.clip(
        pd.to_numeric(audit["barneo_x_discovery_score"], errors="coerce")
        * pd.to_numeric(audit["barneo_x_claim_gate"], errors="coerce"),
        0,
        1,
    )
    audit["barneo_x_positive_only_score"] = np.clip(pd.to_numeric(audit["barneo_x_positive_component_sum"], errors="coerce"), 0, 1)
    rows = []
    for col, interpretation in [
        ("barneo_score", "raw apparent BAR-Neo score; high label metric but unsafe for clean top-k claims"),
        ("barneo_x_positive_only_score", "positive evidence only; exposes why penalties are needed"),
        ("barneo_x_no_penalty_score", "claim gate and metadata cap retained, but component penalties removed"),
        ("barneo_x_no_leakage_gate_score", "component penalties retained, but high/medium leakage gate removed"),
        ("barneo_x_no_metadata_cap_score", "component penalties and leakage gate retained, but missing metadata cap removed"),
        ("barneo_x_claim_safe_score", "full BAR-Neo-X reviewer-facing score"),
    ]:
        row = score_metric_row(audit, col)
        if row is None:
            continue
        row["interpretation"] = interpretation
        rows.append(row)
    return pd.DataFrame(rows)


def topk_rows(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for score_col in ["barneo_x_discovery_score", "barneo_x_claim_safe_score"]:
        ranked = df.sort_values(score_col, ascending=False)
        for k in TOP_KS:
            sub = ranked.head(k)
            label = pd.to_numeric(sub.get("label"), errors="coerce")
            rows.append(
                {
                    "score": score_col,
                    "top_k": k,
                    "precision": float(label.mean()) if label.notna().any() else np.nan,
                    "high_leakage_fraction": float(sub["leakage_risk_level"].astype(str).str.lower().eq("high").mean()),
                    "median_claim_safe_score": float(pd.to_numeric(sub["barneo_x_claim_safe_score"], errors="coerce").median()),
                    "priority_review_rows": int(sub["barneo_x_primary_action"].eq("priority_review_candidate").sum()),
                }
            )
    return pd.DataFrame(rows)


def write_report(output_root: Path, scores: pd.DataFrame, metrics: pd.DataFrame, topk: pd.DataFrame, ablation: pd.DataFrame) -> None:
    top_claim = scores.sort_values("barneo_x_claim_safe_score", ascending=False).head(20)
    cols = [
        "candidate_id",
        "peptide",
        "hla_allele_4digit",
        "label",
        "leakage_risk_level",
        "barneo_x_claim_safe_score",
        "barneo_x_primary_action",
        "barneo_x_top_positive_factors",
        "barneo_x_top_negative_factors",
    ]
    lines = [
        "# BAR-Neo-X Interpretability Boost",
        "",
        "## Position",
        "",
        "BAR-Neo-X is a post-hoc explanation and claim-safe reranking layer. It keeps the high-performing BAR-Neo/BMA evidence, decomposes each candidate into positive and negative score components, and downranks rows that are not safe for clean benchmark claims.",
        "",
        "It is not a new public SOTA claim and not a quantum-advantage claim.",
        "",
        "## Apparent Candidate-Level Metrics",
        "",
        metrics.round(4).to_markdown(index=False) if not metrics.empty else "No metrics available.",
        "",
        "## Top-K Safety Audit",
        "",
        topk.round(4).to_markdown(index=False) if not topk.empty else "No top-k audit available.",
        "",
        "## Ablation Audit",
        "",
        ablation.round(4).to_markdown(index=False) if not ablation.empty else "No ablation audit available.",
        "",
        "## Top Claim-Safe Review Rows",
        "",
        top_claim[[c for c in cols if c in top_claim.columns]].round(4).to_markdown(index=False),
        "",
        "## Interpretation",
        "",
        "Use `barneo_x_discovery_score` for internal hypothesis discovery. Use `barneo_x_claim_safe_score` for reviewer-facing candidate review, because it explicitly penalizes leakage, overlap, sparse expert support, missing patient context, and expert disagreement.",
    ]
    (output_root / "BAR_NEO_X_INTERPRETABILITY_BOOST_REPORT.md").write_text("\n".join(lines) + "\n")


def update_output_manifest(output_root: Path, summary: dict[str, Any]) -> None:
    path = output_root / "run_manifest.json"
    if path.exists():
        try:
            manifest = json.loads(path.read_text())
        except Exception:
            manifest = {}
    else:
        manifest = {}
    manifest.setdefault("output_files", [])
    for name in OUTPUT_FILES:
        if name not in manifest["output_files"]:
            manifest["output_files"].append(name)
    manifest.setdefault("summary", {})
    manifest["summary"].update(
        {
            "n_barneo_x_candidates": int(summary["n_candidates"]),
            "n_barneo_x_priority_review_candidates": int(summary["n_priority_review_candidates"]),
            "barneo_x_top_claim_safe_score": float(summary["top_claim_safe_score"]),
            "n_barneo_x_component_attribution_rows": int(summary["n_component_attribution_rows"]),
        }
    )
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def main() -> None:
    args = parse_args()
    output_root = Path(args.output_root)
    ensure_dir(output_root)
    scored = add_scores(candidate_table(output_root))
    explanation_cols = [
        "candidate_id",
        "peptide",
        "hla_allele_4digit",
        "label",
        "leakage_risk_level",
        "barneo_score",
        "barneo_bma_score",
        "barneo_x_discovery_score",
        "barneo_x_claim_safe_score",
        "barneo_x_discovery_rank",
        "barneo_x_claim_safe_rank",
        "barneo_x_positive_component_sum",
        "barneo_x_penalty_sum",
        "barneo_x_claim_gate",
        "barneo_x_metadata_cap",
        "barneo_x_primary_action",
        "barneo_x_top_positive_factors",
        "barneo_x_top_negative_factors",
        "selected_methods",
        "abstention_reason_primary",
        "bma_abstention_reason_primary",
    ]
    write_tsv(scored[[c for c in explanation_cols if c in scored.columns]], output_root / "barneo_x_candidate_explanations.tsv")
    component_attributions = component_attribution_rows(scored)
    write_tsv(component_attributions, output_root / "barneo_x_component_attributions.tsv")
    write_tsv(
        scored[
            [
                "candidate_id",
                "barneo_x_discovery_score",
                "barneo_x_claim_safe_score",
                "barneo_x_discovery_rank",
                "barneo_x_claim_safe_rank",
                "barneo_x_primary_action",
            ]
        ].sort_values("barneo_x_claim_safe_rank"),
        output_root / "barneo_x_candidate_scores.tsv",
    )
    metrics = metric_rows(scored)
    topk = topk_rows(scored)
    ablation = ablation_rows(scored)
    write_tsv(metrics, output_root / "barneo_x_metric_audit.tsv")
    write_tsv(topk, output_root / "barneo_x_topk_safety_audit.tsv")
    write_tsv(ablation, output_root / "barneo_x_ablation_audit.tsv")
    write_report(output_root, scored, metrics, topk, ablation)
    summary = {
        "n_candidates": int(len(scored)),
        "n_priority_review_candidates": int(scored["barneo_x_primary_action"].eq("priority_review_candidate").sum()),
        "top_claim_safe_score": float(scored["barneo_x_claim_safe_score"].max()),
        "n_component_attribution_rows": int(len(component_attributions)),
    }
    (output_root / "BAR_NEO_X_INTERPRETABILITY_BOOST_SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n")
    update_output_manifest(output_root, summary)
    update_manifest(
        output_root,
        "barneo_x_interpretability_boost",
        {
            "outputs": OUTPUT_FILES,
            "n_candidates": int(len(scored)),
            "n_priority_review_candidates": int(scored["barneo_x_primary_action"].eq("priority_review_candidate").sum()),
            "n_component_attribution_rows": int(len(component_attributions)),
            "warnings": [
                "BAR-Neo-X metrics are candidate-level apparent audits, not an external SOTA claim.",
                "Claim-safe score intentionally downranks high-leakage rows even when labels are positive.",
            ],
        },
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
