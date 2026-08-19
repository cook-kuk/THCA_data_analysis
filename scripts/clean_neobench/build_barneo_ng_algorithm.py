#!/usr/bin/env python3
"""Apply BAR-Neo-NG: claim-gated benchmark-adaptive neoantigen triage."""

from __future__ import annotations

import argparse
import html
import json
from pathlib import Path

import numpy as np
import pandas as pd

from common import dataframe_to_markdown, ensure_dir, update_manifest, write_tsv


HTML_NAME = "barneo_ng_algorithm_results_2026_05_10.html"
OUTPUTS = [
    "barneo_ng_algorithm_scores.tsv",
    "barneo_ng_t1_candidates.tsv",
    "barneo_ng_decision_summary.tsv",
    "BAR_NEO_NG_ALGORITHM_CARD.md",
    "BAR_NEO_NG_RESULTS_KR.md",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo-root", default=".", help="Repository root")
    p.add_argument("--output-root", required=True, help="CLEAN-NeoBench output root")
    p.add_argument("--hub-root", default="project/papers_hub_2026_05_04", help="HTML hub root")
    return p.parse_args()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def num(df: pd.DataFrame, col: str, default: float = 0.0) -> pd.Series:
    if col not in df.columns:
        return pd.Series(default, index=df.index, dtype=float)
    return pd.to_numeric(df[col], errors="coerce").fillna(default)


def bool_series(df: pd.DataFrame, col: str, default: bool = False) -> pd.Series:
    if col not in df.columns:
        return pd.Series(default, index=df.index, dtype=bool)
    return df[col].astype(str).str.lower().isin(["true", "1", "yes"])


def truth_value(value: object) -> bool:
    """Parse booleans robustly after TSV round trips."""
    return str(value).strip().lower() in {"true", "1", "yes"}


def clip01(x: pd.Series) -> pd.Series:
    return pd.to_numeric(x, errors="coerce").fillna(0.0).clip(0.0, 1.0)


def build_ng_scores(stress: pd.DataFrame, kill: pd.DataFrame, t1: pd.DataFrame) -> pd.DataFrame:
    if stress.empty:
        return pd.DataFrame()
    df = stress.copy()
    kill_cols = [
        "candidate_id",
        "reviewer_kill_disposition",
        "reviewer_kill_reason",
        "exact_peptide_train_overlap",
        "exact_peptide_hla_train_overlap",
        "near_peptide_train_overlap",
        "nearest_neighbor_similarity",
        "source_negatives_at_or_above",
        "source_hla_negatives_at_or_above",
        "top_method",
    ]
    if not kill.empty:
        df = df.merge(kill[[c for c in kill_cols if c in kill.columns]], on="candidate_id", how="left")
    else:
        df["reviewer_kill_disposition"] = ""
    if not t1.empty:
        tcols = [
            "candidate_id",
            "peptide_hla_assay_design_ready",
            "assay_design_readiness_score",
            "translational_metadata_score",
            "impact_readiness_score",
            "impact_readiness_tier",
            "missing_high_impact_metadata",
        ]
        df = df.merge(t1[[c for c in tcols if c in t1.columns]], on="candidate_id", how="left")

    discovery = clip01(num(df, "stress_guarded_discovery_score"))
    clean = clip01(num(df, "stress_guarded_claim_safe_score"))
    confidence = clip01(num(df, "stress_guarded_confidence"))
    disagreement = clip01(num(df, "stress_guarded_expert_disagreement"))
    agreement = (1.0 - disagreement).clip(0.0, 1.0)
    clean_method_support = (num(df, "stress_guarded_clean_method_count") / 25.0).clip(0.0, 1.0)
    hla_support = np.log1p(num(df, "hla_allele_support_count")) / np.log1p(max(float(num(df, "hla_allele_support_count").max()), 1.0))
    hla_support = hla_support.fillna(0.0).clip(0.0, 1.0)
    public_penalty = (1.0 - 2.0 * num(df, "public_weight_fraction")).clip(0.25, 1.0)
    fallback_penalty = (1.0 - 1.5 * num(df, "fallback_weight_fraction")).clip(0.35, 1.0)

    leakage = df.get("leakage_risk_level", pd.Series("", index=df.index)).astype(str).str.lower()
    leakage_penalty = pd.Series(1.0, index=df.index)
    leakage_penalty[leakage.eq("medium")] = 0.72
    leakage_penalty[leakage.eq("high")] = 0.25
    exact_overlap = bool_series(df, "exact_peptide_train_overlap") | bool_series(df, "exact_peptide_hla_train_overlap")
    near_overlap = bool_series(df, "near_peptide_train_overlap") | num(df, "nearest_neighbor_similarity").ge(0.80)
    leakage_penalty[exact_overlap] = np.minimum(leakage_penalty[exact_overlap], 0.20)
    leakage_penalty[near_overlap & ~exact_overlap] = np.minimum(leakage_penalty[near_overlap & ~exact_overlap], 0.70)
    mhc_penalty = pd.Series(1.0, index=df.index)
    mhc_penalty[df.get("mhc_class", pd.Series("", index=df.index)).astype(str).str.upper().eq("II")] = 0.20

    assay_ready = bool_series(df, "peptide_hla_assay_design_ready")
    assay_score = clip01(num(df, "assay_design_readiness_score"))
    translational_meta = clip01(num(df, "translational_metadata_score"))
    t1_bonus = pd.Series(0.0, index=df.index)
    t1_pass = df.get("reviewer_kill_disposition", pd.Series("", index=df.index)).eq("passes_current_reviewer_kill_audit")
    t1_bonus[t1_pass] = 0.08
    t1_bonus[t1_pass & assay_ready] = 0.12

    base = (
        0.30 * clean
        + 0.18 * discovery
        + 0.16 * confidence
        + 0.12 * agreement
        + 0.10 * clean_method_support
        + 0.07 * hla_support
        + 0.05 * assay_score
        + 0.02 * translational_meta
        + t1_bonus
    )
    gated = base * leakage_penalty * public_penalty * fallback_penalty * mhc_penalty

    # Reviewer-safe caps: high evidence can rank for discovery, but cannot become a clean claim through a blocked gate.
    cap = pd.Series(0.88, index=df.index)
    cap[leakage.eq("high") | exact_overlap] = 0.24
    cap[near_overlap & ~exact_overlap] = np.minimum(cap[near_overlap & ~exact_overlap], 0.56)
    cap[df.get("reviewer_kill_disposition", pd.Series("", index=df.index)).eq("manual_audit_required")] = 0.56
    cap[df.get("reviewer_kill_disposition", pd.Series("", index=df.index)).eq("blocked_from_clean_claim")] = 0.22
    cap[df.get("mhc_class", pd.Series("", index=df.index)).astype(str).str.upper().eq("II")] = 0.18
    cap[t1_pass & assay_ready] = 0.78
    cap[t1_pass & assay_ready & translational_meta.lt(0.2)] = 0.72
    ng_score = np.minimum(gated, cap).clip(0.0, 1.0)

    decisions = []
    reasons = []
    for i, r in df.iterrows():
        cid = r.get("candidate_id", "")
        disposition = str(r.get("reviewer_kill_disposition", ""))
        reason = []
        if str(r.get("mhc_class", "")).upper() == "II":
            decisions.append("abstain_mhc_ii_separate_benchmark")
            reasons.append("MHC-II candidate excluded from class-I BAR-Neo-NG decision")
            continue
        if disposition == "blocked_from_clean_claim" or bool(exact_overlap.loc[i]) or str(r.get("leakage_risk_level", "")).lower() == "high":
            decisions.append("blocked_from_clean_claim")
            reasons.append(str(r.get("reviewer_kill_reason", "high leakage or exact overlap")))
            continue
        if disposition == "manual_audit_required" or bool(near_overlap.loc[i]):
            decisions.append("T2_resolve_before_claim")
            reasons.append(str(r.get("reviewer_kill_reason", "near-neighbor or unresolved audit caveat")))
            continue
        if disposition == "passes_current_reviewer_kill_audit" and bool(assay_ready.loc[i]):
            decisions.append("T1_phla_assay_design_candidate")
            if float(translational_meta.loc[i]) < 0.2:
                reason.append("pHLA assay-ready but translational metadata blocked")
            reason.append("passed reviewer kill-audit")
            reason.append("low public/fallback dependence")
            reasons.append("; ".join(reason))
            continue
        if truth_value(r.get("stress_guarded_abstain", False)):
            decisions.append("abstain")
            reasons.append(str(r.get("stress_guarded_reason_primary", "stress-guarded abstention")))
            continue
        if float(ng_score.loc[i]) >= 0.45:
            decisions.append("watchlist_manual_review")
            reasons.append("moderate BAR-Neo-NG score; not a clean T1 claim")
        else:
            decisions.append("low_priority_or_abstain")
            reasons.append("insufficient claim-safe evidence under BAR-Neo-NG gates")

    out_cols = [
        "candidate_id",
        "source_name",
        "hla_allele_4digit",
        "peptide",
        "label",
        "mhc_class",
        "leakage_risk_level",
        "stress_guarded_rank_global",
        "stress_guarded_action",
        "stress_guarded_reason_primary",
        "reviewer_kill_disposition",
        "reviewer_kill_reason",
        "impact_readiness_tier",
        "missing_high_impact_metadata",
    ]
    out = df[[c for c in out_cols if c in df.columns]].copy()
    out["barneo_ng_score"] = ng_score
    out["barneo_ng_discovery_component"] = discovery
    out["barneo_ng_clean_component"] = clean
    out["barneo_ng_confidence_component"] = confidence
    out["barneo_ng_agreement_component"] = agreement
    out["barneo_ng_method_support_component"] = clean_method_support
    out["barneo_ng_hla_support_component"] = hla_support
    out["barneo_ng_assay_component"] = assay_score
    out["barneo_ng_translational_metadata_component"] = translational_meta
    out["barneo_ng_leakage_penalty"] = leakage_penalty
    out["barneo_ng_public_penalty"] = public_penalty
    out["barneo_ng_fallback_penalty"] = fallback_penalty
    out["barneo_ng_cap"] = cap
    out["barneo_ng_decision"] = decisions
    out["barneo_ng_reason"] = reasons
    out = out.sort_values(["barneo_ng_score", "stress_guarded_rank_global"], ascending=[False, True])
    out["barneo_ng_rank_global"] = range(1, len(out) + 1)
    if "hla_allele_4digit" in out.columns:
        out["barneo_ng_rank_within_hla"] = out.groupby("hla_allele_4digit")["barneo_ng_score"].rank(ascending=False, method="first")
    return out


def write_reports(output_root: Path, hub_root: Path, scores: pd.DataFrame) -> None:
    if scores.empty:
        summary = pd.DataFrame()
        t1 = pd.DataFrame()
    else:
        summary = scores["barneo_ng_decision"].value_counts().rename_axis("barneo_ng_decision").reset_index(name="n")
        t1 = scores[scores["barneo_ng_decision"].eq("T1_phla_assay_design_candidate")].copy()
    top_cols = [
        "candidate_id",
        "barneo_ng_rank_global",
        "barneo_ng_score",
        "barneo_ng_decision",
        "source_name",
        "hla_allele_4digit",
        "peptide",
        "label",
        "leakage_risk_level",
        "reviewer_kill_disposition",
        "impact_readiness_tier",
        "barneo_ng_reason",
    ]
    card = """# BAR-Neo-NG Algorithm Card

## What It Is

BAR-Neo-NG is a claim-gated BAR-Neo layer. It combines stress-guarded benchmark reliability, reviewer kill-audit status, public/fallback dependency penalties, HLA support, assay readiness, and translational metadata gates.

## What It Does Not Claim

It is not a new SOTA neoantigen predictor, not clinical vaccine selection, not external validation, and not a QK headline claim.

## Decision Logic

- High leakage or exact train overlap caps/blocks clean claims.
- Near-neighbor caveats become `T2_resolve_before_claim`.
- Kill-audit pass plus pHLA assay readiness becomes `T1_phla_assay_design_candidate`.
- Missing translational metadata caps the T1 score and blocks patient/translational claims.
- Public pretrained support is penalized unless row-level overlap audit is clean.
"""
    (output_root / "BAR_NEO_NG_ALGORITHM_CARD.md").write_text(card.strip() + "\n")

    kr = f"""# BAR-Neo-NG Results KR

## 한 줄 결론

우리 새 알고리즘 `BAR-Neo-NG`를 실제 후보 전체에 적용했다. 핵심은 점수만 올리는 게 아니라, leakage/public/fallback/near-neighbor/metadata gate를 걸어 **T1 pHLA assay-design 후보만 안전하게 올리고 나머지는 수동검증/차단/abstain**으로 보내는 것이다.

## Decision summary

{dataframe_to_markdown(summary, max_rows=20)}

## Top BAR-Neo-NG rows

{dataframe_to_markdown(scores[[c for c in top_cols if c in scores.columns]].head(25) if not scores.empty else scores, max_rows=25)}

## T1 BAR-Neo-NG candidates

{dataframe_to_markdown(t1[[c for c in top_cols if c in t1.columns]] if not t1.empty else t1, max_rows=20)}

## Claim boundary

BAR-Neo-NG는 research triage / pHLA assay-design ranking이다. clinical vaccine selection, new SOTA predictor, external validation proven, QK headline claim은 금지한다.
"""
    (output_root / "BAR_NEO_NG_RESULTS_KR.md").write_text(kr.strip() + "\n")

    table_cols = [c for c in top_cols if c in scores.columns]
    rows = []
    for _, r in scores.head(30).iterrows():
        rows.append("<tr>" + "".join(f"<td>{html.escape(str(r[c]))}</td>" for c in table_cols) + "</tr>")
    head = "".join(f"<th>{html.escape(c)}</th>" for c in table_cols)
    stat_html = "".join(f"<div><b>{int(row.n):,}</b><span>{html.escape(str(row.barneo_ng_decision))}</span></div>" for _, row in summary.iterrows())
    page = f"""<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>BAR-Neo-NG Algorithm Results</title>
<style>
body{{margin:0;background:#0d1117;color:#e6edf3;font-family:JetBrains Mono,ui-monospace,Menlo,monospace;line-height:1.55}}
header{{padding:44px 34px 24px;border-bottom:1px solid #2a3441}}h1{{font-family:Georgia,serif;font-size:46px;margin:0 0 10px}}.lead{{color:#c8d1dc;max-width:980px}}
.stats{{display:grid;grid-template-columns:repeat(4,minmax(140px,1fr));gap:10px;margin-top:20px}}.stats div{{border:1px solid #2a3441;background:#101820;padding:12px}}.stats b{{display:block;color:#5eead4;font-size:24px}}.stats span{{color:#9aa7b4;font-size:12px}}
main{{max-width:1320px;margin:auto;padding:24px}}section{{border:1px solid #2a3441;background:#101820;padding:18px;overflow:auto}}table{{width:100%;border-collapse:collapse;font-size:12px}}th,td{{border-bottom:1px solid #2a3441;padding:8px;text-align:left;vertical-align:top;white-space:nowrap}}th{{color:#e3b341}}.path{{padding:0 24px 36px;color:#9aa7b4;font-size:12px}}
</style></head><body><header><h1>BAR-Neo-NG Algorithm Results</h1><p class="lead">Claim-gated BAR-Neo algorithm layer: benchmark reliability plus reviewer kill-audit, assay readiness, public/fallback penalties, and metadata gates.</p><div class="stats">{stat_html}</div></header><main><section><table><thead><tr>{head}</tr></thead><tbody>{''.join(rows)}</tbody></table></section></main><p class="path">Source: {html.escape(str(output_root / 'barneo_ng_algorithm_scores.tsv'))}</p></body></html>"""
    hub_root.mkdir(parents=True, exist_ok=True)
    (hub_root / HTML_NAME).write_text(page)


def main() -> None:
    args = parse_args()
    repo_root = Path(args.repo_root).resolve()
    output_root = Path(args.output_root)
    if not output_root.is_absolute():
        output_root = repo_root / output_root
    hub_root = Path(args.hub_root)
    if not hub_root.is_absolute():
        hub_root = repo_root / hub_root
    ensure_dir(output_root)
    ensure_dir(hub_root)

    stress = read_tsv(output_root / "barneo_stress_guarded_candidate_scores.tsv")
    kill = read_tsv(output_root / "barneo_high_impact_reviewer_kill_audit.tsv")
    t1 = read_tsv(output_root / "barneo_t1_translational_readiness.tsv")
    scores = build_ng_scores(stress, kill, t1)
    t1_scores = scores[scores["barneo_ng_decision"].eq("T1_phla_assay_design_candidate")].copy() if not scores.empty else pd.DataFrame()
    summary = scores["barneo_ng_decision"].value_counts().rename_axis("barneo_ng_decision").reset_index(name="n") if not scores.empty else pd.DataFrame()

    write_tsv(scores, output_root / "barneo_ng_algorithm_scores.tsv")
    write_tsv(t1_scores, output_root / "barneo_ng_t1_candidates.tsv")
    write_tsv(summary, output_root / "barneo_ng_decision_summary.tsv")
    write_reports(output_root, hub_root, scores)

    outputs = OUTPUTS + [str(hub_root / HTML_NAME)]
    manifest_path = output_root / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    manifest.setdefault("output_files", [])
    for out in outputs:
        if str(out) not in manifest["output_files"]:
            manifest["output_files"].append(str(out))
    manifest.setdefault("summary", {})
    manifest["summary"].update(
        {
            "n_barneo_ng_scores": int(len(scores)),
            "n_barneo_ng_t1_candidates": int(len(t1_scores)),
            "n_barneo_ng_decision_rows": int(len(summary)),
            "barneo_ng_top_score": float(scores["barneo_ng_score"].max()) if not scores.empty else 0.0,
            "barneo_ng_html": str(hub_root / HTML_NAME),
        }
    )
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    update_manifest(
        output_root,
        "barneo_ng_algorithm",
        {
            "outputs": outputs,
            "n_scores": int(len(scores)),
            "warnings": ["BAR-Neo-NG is claim-gated research triage, not SOTA or clinical selection."],
        },
    )
    print(f"[barneo-ng] scores={len(scores)} t1={len(t1_scores)} top={manifest['summary']['barneo_ng_top_score']:.4f}")


if __name__ == "__main__":
    main()
