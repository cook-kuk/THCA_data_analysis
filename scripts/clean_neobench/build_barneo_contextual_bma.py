#!/usr/bin/env python3
"""Build context-aware BAR-Neo-BMA scores.

Global BMA weights are useful, but real neoantigen review needs different
expert reliability in low-prevalence, rare-HLA, external/holdout, and
leakage-risk settings. This script applies benchmark-derived contextual
penalties and bonuses without training a new black-box model.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from common import KOREAN_HLA_ALLELES, dataframe_to_markdown, ensure_dir, update_manifest, write_tsv


TRAIN_SOURCES = {"CEDAR", "TESLA_mmc4", "NEPdb", "TESLA_mmc7_validation"}
PUBLIC_ROLE = "caveated_public_comparator"
UNCERTAINTY_ROLE = "uncertainty_only"
FALLBACK_ROLE = "bounded_fallback"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=".", help="Repository root")
    parser.add_argument("--output-root", required=True, help="CLEAN-NeoBench output directory")
    parser.add_argument("--max-methods", type=int, default=7, help="Maximum contextual experts per candidate")
    return parser.parse_args()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def num(value: Any, default: float = 0.0) -> float:
    try:
        x = float(value)
        return x if math.isfinite(x) else default
    except Exception:
        return default


def truthy(value: Any) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def candidate_context_features(master: pd.DataFrame) -> pd.DataFrame:
    df = master.copy()
    df["source_positive_prevalence"] = df["source_name"].map(df.groupby("source_name")["label"].mean())
    df["source_n"] = df["source_name"].map(df.groupby("source_name")["candidate_id"].count())
    df["hla_allele_support_count"] = df["hla_allele_4digit"].map(df.groupby("hla_allele_4digit")["candidate_id"].count())
    df["hla_allele_positive_prevalence"] = df["hla_allele_4digit"].map(df.groupby("hla_allele_4digit")["label"].mean())
    df["distribution_partition"] = np.where(df["source_name"].isin(TRAIN_SOURCES), "local_train_pool", "external_or_holdout")
    df["is_low_prevalence_context"] = (
        df["split_low_prevalence"].astype(str).eq("low_prevalence")
        | pd.to_numeric(df["source_positive_prevalence"], errors="coerce").lt(0.10)
    )
    df["is_rare_hla_context"] = pd.to_numeric(df["hla_allele_support_count"], errors="coerce").lt(20)
    df["is_low_support_hla_context"] = pd.to_numeric(df["hla_allele_support_count"], errors="coerce").lt(50)
    df["is_external_context"] = df["distribution_partition"].ne("local_train_pool")
    df["is_high_leakage_context"] = df["leakage_risk_level"].astype(str).str.lower().eq("high")
    df["is_korean_hla_context"] = df["hla_allele_4digit"].isin(KOREAN_HLA_ALLELES)
    return df


def context_labels(row: pd.Series) -> list[str]:
    labels = ["global"]
    labels.append(f"source:{row.get('source_name', 'NA')}")
    labels.append(f"hla:{row.get('hla_allele_4digit', 'NA')}")
    if truthy(row.get("is_low_prevalence_context")):
        labels.append("low_prevalence")
    if truthy(row.get("is_rare_hla_context")):
        labels.append("rare_hla")
    elif truthy(row.get("is_low_support_hla_context")):
        labels.append("low_support_hla")
    if truthy(row.get("is_external_context")):
        labels.append("external_or_holdout")
    if truthy(row.get("is_high_leakage_context")):
        labels.append("high_leakage_review_only")
    if truthy(row.get("is_korean_hla_context")):
        labels.append("korean_hla_focus")
    return [str(x) for x in labels if str(x) not in {"", "nan", "None"}]


def build_perf_maps(split_metrics: pd.DataFrame) -> dict[tuple[str, str, str], float]:
    """Return a simple 0-1 performance utility by context/method."""

    perf: dict[tuple[str, str, str], float] = {}
    if split_metrics.empty:
        return perf
    sm = split_metrics.copy()
    sm["AUPRC"] = pd.to_numeric(sm["AUPRC"], errors="coerce")
    sm["top10_precision"] = pd.to_numeric(sm["top10_precision"], errors="coerce")
    sm["calibration_ece"] = pd.to_numeric(sm["calibration_ece"], errors="coerce")
    sm["context_utility"] = (
        0.55 * sm["AUPRC"].fillna(0)
        + 0.35 * sm["top10_precision"].fillna(0)
        + 0.10 * (1.0 - sm["calibration_ece"].fillna(0.25).clip(0, 1))
    ).clip(0, 1)
    for _, r in sm.iterrows():
        contract = str(r.get("split_contract", ""))
        group = str(r.get("split_group", ""))
        method = str(r.get("method_name", ""))
        util = num(r.get("context_utility"), math.nan)
        if not math.isfinite(util):
            continue
        if contract == "source_heldout":
            perf[("source", group, method)] = util
        if contract in {"hla_heldout", "korean_hla_focus"}:
            perf[("hla", group, method)] = max(util, perf.get(("hla", group, method), 0.0))
        if contract == "low_prevalence_heldout":
            perf[("low_prevalence", group, method)] = util
    return perf


def vulnerability_map(vulnerability: pd.DataFrame) -> dict[str, dict[str, Any]]:
    if vulnerability.empty:
        return {}
    return {str(r["method_name"]): r.to_dict() for _, r in vulnerability.iterrows()}


def method_penalty_for_context(method: pd.Series, candidate: pd.Series, vuln: dict[str, dict[str, Any]], perf: dict[tuple[str, str, str], float]) -> tuple[float, list[str]]:
    name = str(method["method_name"])
    role = str(method.get("method_role", ""))
    family = str(method.get("method_family", ""))
    v = vuln.get(name, {})
    multiplier = 1.0
    reasons: list[str] = []

    source = str(candidate.get("source_name", ""))
    hla = str(candidate.get("hla_allele_4digit", ""))
    source_perf = perf.get(("source", source, name), math.nan)
    hla_perf = perf.get(("hla", hla, name), math.nan)
    if math.isfinite(source_perf):
        multiplier *= 0.65 + 0.70 * source_perf
        reasons.append(f"source-heldout utility={source_perf:.3f}")
    if math.isfinite(hla_perf):
        multiplier *= 0.70 + 0.60 * hla_perf
        reasons.append(f"HLA-heldout utility={hla_perf:.3f}")

    source_corr = abs(num(v.get("source_score_prevalence_corr"), 0.0))
    if source_corr >= 0.45:
        multiplier *= math.exp(-0.20 * source_corr)
        reasons.append("source-prior coupling penalty")

    if truthy(candidate.get("is_low_prevalence_context")):
        fp = num(v.get("low_prevalence_high_ranked_negative_rate"), 0.0)
        multiplier *= math.exp(-1.10 * fp)
        if fp > 0:
            reasons.append(f"low-prevalence FP penalty={fp:.3f}")
        if family == "internal_baseline":
            multiplier *= 1.08
            reasons.append("internal-baseline low-prevalence bonus")

    if truthy(candidate.get("is_rare_hla_context")) or truthy(candidate.get("is_low_support_hla_context")):
        rare_fp = num(v.get("rare_hla_high_ranked_negative_rate"), 0.0)
        rare_fn = num(v.get("rare_hla_missed_positive_rate"), 0.0)
        multiplier *= math.exp(-0.90 * max(rare_fp, rare_fn))
        if max(rare_fp, rare_fn) > 0:
            reasons.append(f"rare-HLA instability penalty={max(rare_fp, rare_fn):.3f}")
        if role == UNCERTAINTY_ROLE:
            multiplier *= 1.12
            reasons.append("uncertainty branch sparse-HLA bonus")

    if truthy(candidate.get("is_external_context")):
        ext_fp = num(v.get("external_high_ranked_negative_rate"), 0.0)
        ext_fn = num(v.get("external_missed_positive_rate"), 0.0)
        multiplier *= math.exp(-0.85 * max(ext_fp, ext_fn))
        if max(ext_fp, ext_fn) > 0:
            reasons.append(f"external fragility penalty={max(ext_fp, ext_fn):.3f}")

    if truthy(candidate.get("is_high_leakage_context")):
        leak_fp = num(v.get("high_leakage_high_ranked_negative_rate"), 0.0)
        multiplier *= math.exp(-0.80 * leak_fp)
        reasons.append("high-leakage review-only context")

    if truthy(candidate.get("is_korean_hla_context")):
        if math.isfinite(hla_perf):
            multiplier *= 1.05
            reasons.append("Korean-HLA observed split utility bonus")
        elif role == PUBLIC_ROLE:
            multiplier *= 0.92
            reasons.append("Korean-HLA public-overlap unresolved penalty")

    if role == PUBLIC_ROLE:
        multiplier *= 0.55
        reasons.append("public pretrained caveated support only")
    if role == FALLBACK_ROLE:
        multiplier *= 0.82
        reasons.append("bounded fallback contribution cap")
    if role == UNCERTAINTY_ROLE:
        multiplier *= 0.70
        reasons.append("uncertainty-only not ranking headline")

    return max(0.02, min(3.0, multiplier)), reasons


def context_weight_preview(
    weights: pd.DataFrame,
    master: pd.DataFrame,
    vuln: dict[str, dict[str, Any]],
    perf: dict[tuple[str, str, str], float],
    max_methods: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    representatives = []
    contexts = [
        "global",
        "low_prevalence",
        "rare_hla",
        "low_support_hla",
        "external_or_holdout",
        "high_leakage_review_only",
        "korean_hla_focus",
    ]
    for source in sorted(master["source_name"].dropna().astype(str).unique()):
        contexts.append(f"source:{source}")
    for allele in sorted(set(master["hla_allele_4digit"].dropna().astype(str)) & set(KOREAN_HLA_ALLELES)):
        contexts.append(f"hla:{allele}")

    base = pd.Series(
        {
            "source_name": "GLOBAL",
            "hla_allele_4digit": "GLOBAL",
            "is_low_prevalence_context": False,
            "is_rare_hla_context": False,
            "is_low_support_hla_context": False,
            "is_external_context": False,
            "is_high_leakage_context": False,
            "is_korean_hla_context": False,
        }
    )
    for ctx in contexts:
        rep = base.copy()
        if ctx.startswith("source:"):
            rep["source_name"] = ctx.split(":", 1)[1]
        if ctx.startswith("hla:"):
            rep["hla_allele_4digit"] = ctx.split(":", 1)[1]
        if ctx == "global":
            rep["source_name"] = "GLOBAL"
            rep["hla_allele_4digit"] = "GLOBAL"
        rep["is_low_prevalence_context"] = ctx == "low_prevalence"
        rep["is_rare_hla_context"] = ctx == "rare_hla"
        rep["is_low_support_hla_context"] = ctx in {"rare_hla", "low_support_hla"}
        rep["is_external_context"] = ctx == "external_or_holdout" or (
            ctx.startswith("source:") and str(rep.get("source_name", "")) not in TRAIN_SOURCES
        )
        rep["is_high_leakage_context"] = ctx == "high_leakage_review_only"
        rep["is_korean_hla_context"] = ctx == "korean_hla_focus" or str(rep.get("hla_allele_4digit", "")) in KOREAN_HLA_ALLELES
        alpha_rows = []
        for _, m in weights.iterrows():
            mult, reasons = method_penalty_for_context(m, rep, vuln, perf)
            alpha = num(m.get("posterior_alpha"), 1.0) * mult
            alpha_rows.append({**m.to_dict(), "context_label": ctx, "context_multiplier": mult, "contextual_alpha": alpha, "context_adjustment_reasons": "; ".join(dict.fromkeys(reasons))})
        df = pd.DataFrame(alpha_rows)
        total = df["contextual_alpha"].sum()
        df["contextual_weight"] = df["contextual_alpha"] / total if total else 0.0
        df = df.sort_values("contextual_weight", ascending=False)
        rows.append(df)
        selected = df.head(max_methods).copy()
        selected["selected_rank"] = range(1, len(selected) + 1)
        representatives.append(selected)
    return pd.concat(rows, ignore_index=True), pd.concat(representatives, ignore_index=True)


def select_candidate_methods(
    candidate: pd.Series,
    method_rows: pd.DataFrame,
    available_scores: pd.Series,
    vuln: dict[str, dict[str, Any]],
    perf: dict[tuple[str, str, str], float],
    max_methods: int,
    clean_only: bool,
) -> pd.DataFrame:
    rows = []
    for _, m in method_rows.iterrows():
        name = str(m["method_name"])
        if name not in available_scores.index or pd.isna(available_scores.get(name)):
            continue
        role = str(m.get("method_role", ""))
        if clean_only and role == PUBLIC_ROLE:
            continue
        mult, reasons = method_penalty_for_context(m, candidate, vuln, perf)
        alpha = num(m.get("posterior_alpha"), 1.0) * mult
        if role == FALLBACK_ROLE:
            alpha = min(alpha, 0.12 * method_rows["posterior_alpha"].sum())
        if role == UNCERTAINTY_ROLE:
            alpha = min(alpha, 0.08 * method_rows["posterior_alpha"].sum())
        rows.append(
            {
                "method_name": name,
                "method_role": role,
                "method_family": m.get("method_family", ""),
                "score": float(available_scores.get(name)),
                "contextual_alpha": alpha,
                "context_multiplier": mult,
                "uses_public_pretraining": bool(m.get("uses_public_pretraining", False)),
                "context_adjustment_reasons": "; ".join(dict.fromkeys(reasons)),
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out = out.sort_values("contextual_alpha", ascending=False).head(max_methods).copy()
    total = out["contextual_alpha"].sum()
    out["contextual_weight"] = out["contextual_alpha"] / total if total else 1.0 / len(out)
    return out


def weighted_score(selection: pd.DataFrame) -> tuple[float, float, float]:
    if selection.empty:
        return math.nan, math.nan, math.nan
    scores = selection["score"].astype(float).to_numpy()
    weights = selection["contextual_weight"].astype(float).to_numpy()
    weights = weights / weights.sum() if weights.sum() else np.ones(len(weights)) / len(weights)
    score = float(np.dot(weights, scores))
    uncertainty = float(np.sqrt(max(0.0, np.dot(weights, (scores - score) ** 2))))
    disagreement = float(scores.max() - scores.min()) if len(scores) else 0.0
    return score, uncertainty, disagreement


def candidate_scores(master: pd.DataFrame, scores: pd.DataFrame, weights: pd.DataFrame, vuln: dict[str, dict[str, Any]], perf: dict[tuple[str, str, str], float], max_methods: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    pivot = scores.pivot_table(index="candidate_id", columns="method_name", values="score_calibrated", aggfunc="max")
    rows = []
    audit_rows = []
    method_rows = weights.copy()
    for _, cand in master.iterrows():
        cid = cand["candidate_id"]
        available = pivot.loc[cid].dropna() if cid in pivot.index else pd.Series(dtype=float)
        sel_all = select_candidate_methods(cand, method_rows, available, vuln, perf, max_methods, clean_only=False)
        sel_clean = select_candidate_methods(cand, method_rows, available, vuln, perf, max_methods, clean_only=True)
        score, unc, disagree = weighted_score(sel_all)
        clean_score, clean_unc, clean_disagree = weighted_score(sel_clean)
        if not math.isfinite(score):
            score = 0.0
            unc = 1.0
            disagree = 1.0
        if not math.isfinite(clean_score):
            clean_score = math.nan
            clean_unc = math.nan
            clean_disagree = math.nan

        reasons = []
        if truthy(cand.get("is_high_leakage_context")):
            reasons.append("High leakage risk blocks clean contextual claim")
        if truthy(cand.get("is_low_prevalence_context")):
            reasons.append("Low-prevalence source; false-positive pressure audit required")
        if truthy(cand.get("is_rare_hla_context")):
            reasons.append("Rare HLA allele underrepresented in benchmark")
        elif truthy(cand.get("is_low_support_hla_context")):
            reasons.append("Low HLA allele support in benchmark")
        if truthy(cand.get("is_external_context")):
            reasons.append("External/holdout source context")
        public_used = bool((sel_all["method_role"].eq(PUBLIC_ROLE)).any()) if not sel_all.empty else False
        if public_used:
            reasons.append("Public pretrained expert used only as caveated support")
        if len(sel_all) < 3:
            reasons.append("Sparse contextual expert support")
        if disagree > 0.45:
            reasons.append("High contextual expert disagreement")
        if unc > 0.20:
            reasons.append("High contextual posterior uncertainty")

        confidence = max(0.0, min(1.0, 1.0 - 0.75 * unc - 0.35 * disagree))
        if truthy(cand.get("is_high_leakage_context")):
            confidence = min(confidence, 0.38)
        if truthy(cand.get("is_low_prevalence_context")):
            confidence = min(confidence, 0.55)
        if truthy(cand.get("is_rare_hla_context")):
            confidence = min(confidence, 0.40)
        if public_used:
            confidence = min(confidence, 0.65)
        if len(sel_all) < 3:
            confidence = min(confidence, 0.45)

        clean_claim_allowed = (
            len(reasons) == 0
            and not public_used
            and math.isfinite(clean_score)
            and confidence >= 0.55
        )
        contexts = context_labels(cand)
        rows.append(
            {
                "candidate_id": cid,
                "context_labels": ";".join(contexts),
                "contextual_bma_score": score,
                "clean_contextual_bma_score": clean_score,
                "contextual_confidence_score": confidence,
                "contextual_confidence_bin": "high" if confidence >= 0.70 else "medium" if confidence >= 0.45 else "low",
                "contextual_posterior_uncertainty": unc,
                "contextual_expert_disagreement": disagree,
                "clean_contextual_posterior_uncertainty": clean_unc,
                "clean_contextual_expert_disagreement": clean_disagree,
                "contextual_abstain": bool(reasons) or confidence < 0.45,
                "contextual_abstention_reason_primary": reasons[0] if reasons else "",
                "contextual_abstention_reason_all": "; ".join(dict.fromkeys(reasons)),
                "contextual_clean_claim_allowed": clean_claim_allowed,
                "selected_contextual_methods": ",".join(sel_all["method_name"].astype(str)) if not sel_all.empty else "",
                "selected_clean_contextual_methods": ",".join(sel_clean["method_name"].astype(str)) if not sel_clean.empty else "",
                "selected_contextual_method_count": int(len(sel_all)),
                "selected_clean_contextual_method_count": int(len(sel_clean)),
                "caveated_public_expert_used": public_used,
                "source_name": cand.get("source_name", ""),
                "hla_allele_4digit": cand.get("hla_allele_4digit", ""),
                "peptide": cand.get("peptide", ""),
                "label": cand.get("label", np.nan),
                "source_positive_prevalence": cand.get("source_positive_prevalence", np.nan),
                "hla_allele_support_count": cand.get("hla_allele_support_count", np.nan),
                "leakage_risk_level": cand.get("leakage_risk_level", ""),
                "split_low_prevalence": cand.get("split_low_prevalence", ""),
            }
        )
        if not sel_all.empty:
            a = sel_all.copy()
            a.insert(0, "candidate_id", cid)
            a["clean_only_selection"] = False
            audit_rows.append(a)
        if not sel_clean.empty:
            c = sel_clean.copy()
            c.insert(0, "candidate_id", cid)
            c["clean_only_selection"] = True
            audit_rows.append(c)

    out = pd.DataFrame(rows)
    out["contextual_bma_rank_global"] = out["contextual_bma_score"].rank(ascending=False, method="first").astype(int)
    out["clean_contextual_bma_rank_global"] = out["clean_contextual_bma_score"].rank(ascending=False, method="first", na_option="bottom").astype(int)
    out = out.sort_values("contextual_bma_rank_global")
    audit = pd.concat(audit_rows, ignore_index=True, sort=False) if audit_rows else pd.DataFrame()
    return out, audit


def build_context_summary(candidates: pd.DataFrame) -> pd.DataFrame:
    if candidates.empty:
        return pd.DataFrame()
    rows = []
    context_set = sorted({ctx for labels in candidates["context_labels"].fillna("") for ctx in str(labels).split(";") if ctx})
    for ctx in context_set:
        g = candidates[candidates["context_labels"].fillna("").str.contains(ctx, regex=False)]
        label = pd.to_numeric(g["label"], errors="coerce")
        rows.append(
            {
                "context_label": ctx,
                "n_candidates": int(len(g)),
                "n_pos": int(label.fillna(0).sum()),
                "positive_prevalence": float(label.mean()) if label.notna().any() else math.nan,
                "mean_contextual_bma_score": float(g["contextual_bma_score"].mean()),
                "mean_clean_contextual_bma_score": float(pd.to_numeric(g["clean_contextual_bma_score"], errors="coerce").mean()),
                "mean_confidence": float(g["contextual_confidence_score"].mean()),
                "abstention_rate": float(g["contextual_abstain"].mean()),
                "clean_claim_allowed": int(g["contextual_clean_claim_allowed"].sum()),
                "caveated_public_used_rate": float(g["caveated_public_expert_used"].mean()),
            }
        )
    return pd.DataFrame(rows).sort_values(["abstention_rate", "n_candidates"], ascending=[False, False])


def write_report(output_root: Path, context_weights: pd.DataFrame, selected: pd.DataFrame, scores: pd.DataFrame, summary: pd.DataFrame) -> None:
    top_scores = scores.head(30)
    clean_prompts = scores[scores["contextual_clean_claim_allowed"].astype(bool)].sort_values("clean_contextual_bma_score", ascending=False)
    text = f"""# BAR-Neo Contextual BMA Report

## Purpose

Contextual BAR-Neo-BMA adjusts method posterior weights by candidate context: source, HLA support, Korean-HLA focus, low-prevalence stress, external/holdout status, and leakage risk. It is a reliability controller, not a new deep predictor.

## Practical Algorithm

1. Start from global BAR-Neo-BMA method posterior weights.
2. Apply benchmark-derived penalties from distribution-error audit.
3. Add source-heldout and HLA-heldout utility bonuses when available.
4. Penalize public pretrained methods as caveated support until row-level training-corpus overlap audit is complete.
5. Compute both all-expert contextual score and clean-contextual score that excludes public pretrained comparators.
6. Abstain when context risk, sparse support, public caveats, leakage, or disagreement make a clean claim unsafe.

## Contextual Method Weight Preview

{dataframe_to_markdown(context_weights[["context_label", "method_name", "method_role", "method_family", "contextual_weight", "context_multiplier", "context_adjustment_reasons"]].head(50) if not context_weights.empty else context_weights, max_rows=50)}

## Context Selector Preview

{dataframe_to_markdown(selected[["context_label", "selected_rank", "method_name", "method_role", "contextual_weight", "context_adjustment_reasons"]].head(80) if not selected.empty else selected, max_rows=80)}

## Candidate Scores

{dataframe_to_markdown(top_scores[["candidate_id", "source_name", "hla_allele_4digit", "contextual_bma_score", "clean_contextual_bma_score", "contextual_confidence_score", "contextual_abstain", "contextual_abstention_reason_primary", "selected_contextual_methods"]], max_rows=30)}

## Clean Contextual Claim Prompts

{dataframe_to_markdown(clean_prompts[["candidate_id", "source_name", "hla_allele_4digit", "clean_contextual_bma_score", "contextual_confidence_score", "selected_clean_contextual_methods"]].head(30) if not clean_prompts.empty else clean_prompts, max_rows=30)}

## Context Summary

{dataframe_to_markdown(summary, max_rows=60)}

## Claim Boundary

Allowed: context-aware reliability weighting, agentic ensemble triage, caveated all-expert support, clean internal score view, abstention.

Forbidden: clinical vaccine selection, new SOTA predictor, quantum advantage, public pretrained tools as clean baselines without row-level overlap audit.
"""
    (output_root / "BAR_NEO_CONTEXTUAL_BMA_REPORT.md").write_text(text.strip() + "\n")

    kr = f"""# BAR-Neo Contextual BMA KR

## 한 줄 결론

이제 global ensemble이 아니라 **상황별 ensemble**이다. CEDAR/TESLA/외부 source, rare-HLA, Korean-HLA, leakage-risk에 따라 method weight가 바뀐다.

## 제일 중요한 점

- public pretrained method는 점수 계산에는 caveated support로 들어갈 수 있지만 clean claim에는 쓰지 않는다.
- `clean_contextual_bma_score`는 public comparator를 제외한 내부/anchor/fallback view다.
- high leakage, rare HLA, low prevalence, external source는 confidence cap과 abstention reason을 만든다.
- QK 계열은 bounded fallback contribution으로만 유지한다.

## 상위 후보

{dataframe_to_markdown(top_scores[["candidate_id", "source_name", "hla_allele_4digit", "contextual_bma_score", "clean_contextual_bma_score", "contextual_confidence_score", "contextual_abstention_reason_primary"]], max_rows=25)}

## context 요약

{dataframe_to_markdown(summary.head(40), max_rows=40)}
"""
    (output_root / "BAR_NEO_CONTEXTUAL_BMA_KR.md").write_text(kr.strip() + "\n")


def main() -> None:
    args = parse_args()
    output_root = Path(args.output_root)
    ensure_dir(output_root)

    master = read_tsv(output_root / "clean_neobench_master.tsv")
    method_scores = read_tsv(output_root / "clean_neobench_method_scores.tsv")
    weights = read_tsv(output_root / "barneo_bma_method_weights.tsv")
    split_metrics = read_tsv(output_root / "clean_neobench_split_metrics.tsv")
    vulnerability = read_tsv(output_root / "clean_neobench_method_distribution_vulnerability.tsv")
    if master.empty or method_scores.empty or weights.empty:
        raise SystemExit("master, method scores, and BMA method weights are required")

    master_ctx = candidate_context_features(master)
    perf = build_perf_maps(split_metrics)
    vuln = vulnerability_map(vulnerability)
    context_weights, selected_preview = context_weight_preview(weights, master_ctx, vuln, perf, args.max_methods)
    candidates, selector_audit = candidate_scores(master_ctx, method_scores, weights, vuln, perf, args.max_methods)
    summary = build_context_summary(candidates)

    write_tsv(context_weights, output_root / "barneo_contextual_bma_method_weights.tsv")
    write_tsv(selected_preview, output_root / "barneo_contextual_bma_selector_preview.tsv")
    write_tsv(candidates, output_root / "barneo_contextual_bma_candidate_scores.tsv")
    write_tsv(selector_audit, output_root / "barneo_contextual_bma_candidate_selector_audit.tsv")
    write_tsv(summary, output_root / "barneo_contextual_bma_context_summary.tsv")
    write_report(output_root, context_weights, selected_preview, candidates, summary)

    outputs = [
        "barneo_contextual_bma_method_weights.tsv",
        "barneo_contextual_bma_selector_preview.tsv",
        "barneo_contextual_bma_candidate_scores.tsv",
        "barneo_contextual_bma_candidate_selector_audit.tsv",
        "barneo_contextual_bma_context_summary.tsv",
        "BAR_NEO_CONTEXTUAL_BMA_REPORT.md",
        "BAR_NEO_CONTEXTUAL_BMA_KR.md",
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
            "n_contextual_bma_candidates": int(len(candidates)),
            "n_contextual_bma_abstain": int(candidates["contextual_abstain"].sum()) if not candidates.empty else 0,
            "n_contextual_bma_clean_claim_allowed": int(candidates["contextual_clean_claim_allowed"].sum()) if not candidates.empty else 0,
            "n_contextual_bma_context_weight_rows": int(len(context_weights)),
        }
    )
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    update_manifest(
        output_root,
        "barneo_contextual_bma",
        {
            "outputs": outputs,
            "n_candidates": int(len(candidates)),
            "n_abstain": int(candidates["contextual_abstain"].sum()) if not candidates.empty else 0,
            "n_clean_claim_allowed": int(candidates["contextual_clean_claim_allowed"].sum()) if not candidates.empty else 0,
            "warnings": [
                "Contextual BMA is a reliability weighting layer, not a new SOTA predictor.",
                "Public pretrained experts are caveated support and excluded from clean_contextual_bma_score.",
            ],
        },
    )
    print(
        "[barneo-contextual-bma] "
        f"candidates={len(candidates)} abstain={int(candidates['contextual_abstain'].sum()) if not candidates.empty else 0} "
        f"clean_claim_allowed={int(candidates['contextual_clean_claim_allowed'].sum()) if not candidates.empty else 0}"
    )


if __name__ == "__main__":
    main()
