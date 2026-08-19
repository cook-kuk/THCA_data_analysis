#!/usr/bin/env python3
"""Build high-impact evidence cards for stress-guarded BAR-Neo leads."""

from __future__ import annotations

import argparse
import html
import json
from pathlib import Path

import numpy as np
import pandas as pd

from common import dataframe_to_markdown, ensure_dir, update_manifest, write_tsv


PAGE_NAME = "barneo_high_impact_candidate_pack_2026_05_10.html"
OUTPUTS = [
    "barneo_high_impact_lead_candidates.tsv",
    "barneo_high_impact_topk_audit.tsv",
    "barneo_high_impact_ablation_audit.tsv",
    "barneo_high_impact_method_support_matrix.tsv",
    "BAR_NEO_HIGH_IMPACT_CANDIDATE_PACK.md",
    "BAR_NEO_HIGH_IMPACT_CANDIDATE_PACK_KR.md",
    PAGE_NAME,
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo-root", default=".", help="Repository root")
    p.add_argument("--output-root", required=True, help="CLEAN-NeoBench output root")
    p.add_argument("--hub-root", default="project/papers_hub_2026_05_04", help="HTML hub root")
    p.add_argument("--lead-n", type=int, default=12, help="Number of leads to profile")
    return p.parse_args()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def num(s: pd.Series, default: float = 0.0) -> pd.Series:
    return pd.to_numeric(s, errors="coerce").fillna(default)


def fmt(x: object) -> str:
    try:
        if pd.isna(x):
            return "NA"
        f = float(x)
        if abs(f) >= 10 and f.is_integer():
            return f"{int(f):,}"
        return f"{f:.3f}"
    except Exception:
        return str(x)


def build_topk(scores: pd.DataFrame) -> pd.DataFrame:
    if scores.empty:
        return pd.DataFrame()
    s = scores.sort_values("stress_guarded_rank_global").copy()
    s["label_num"] = pd.to_numeric(s.get("label"), errors="coerce")
    total_pos = int(s["label_num"].eq(1).sum())
    rows = []
    for k in [1, 3, 4, 5, 10, 20, 50, 100]:
        top = s.head(k)
        n = len(top)
        pos = int(top["label_num"].eq(1).sum())
        rows.append(
            {
                "slice": f"top_{k}_global",
                "n": n,
                "n_pos": pos,
                "precision": pos / n if n else np.nan,
                "recall_of_all_positives": pos / total_pos if total_pos else np.nan,
                "mean_score": float(num(top["stress_guarded_final_review_score"]).mean()) if n else np.nan,
                "mean_confidence": float(num(top["stress_guarded_confidence"]).mean()) if n else np.nan,
                "n_high_leakage": int(top.get("is_high_leakage", pd.Series(False, index=top.index)).astype(bool).sum()),
                "mean_public_weight_fraction": float(num(top.get("public_weight_fraction", pd.Series(index=top.index))).mean()) if n else np.nan,
                "mean_fallback_weight_fraction": float(num(top.get("fallback_weight_fraction", pd.Series(index=top.index))).mean()) if n else np.nan,
                "claim_boundary": "benchmark-label audit, not external validation",
            }
        )
    for action in ["claim_safe_candidate_after_manual_audit", "priority_review_candidate", "support_review_candidate"]:
        top = s[s["stress_guarded_action"].eq(action)]
        n = len(top)
        pos = int(top["label_num"].eq(1).sum()) if n else 0
        rows.append(
            {
                "slice": action,
                "n": n,
                "n_pos": pos,
                "precision": pos / n if n else np.nan,
                "recall_of_all_positives": pos / total_pos if total_pos else np.nan,
                "mean_score": float(num(top["stress_guarded_final_review_score"]).mean()) if n else np.nan,
                "mean_confidence": float(num(top["stress_guarded_confidence"]).mean()) if n else np.nan,
                "n_high_leakage": int(top.get("is_high_leakage", pd.Series(False, index=top.index)).astype(bool).sum()) if n else 0,
                "mean_public_weight_fraction": float(num(top.get("public_weight_fraction", pd.Series(index=top.index))).mean()) if n else np.nan,
                "mean_fallback_weight_fraction": float(num(top.get("fallback_weight_fraction", pd.Series(index=top.index))).mean()) if n else np.nan,
                "claim_boundary": "benchmark-label audit, not external validation",
            }
        )
    return pd.DataFrame(rows)


def build_leads(scores: pd.DataFrame, lead_n: int) -> pd.DataFrame:
    if scores.empty:
        return pd.DataFrame()
    s = scores.sort_values("stress_guarded_rank_global").copy()
    ordered = pd.concat(
        [
            s[s["stress_guarded_action"].eq("claim_safe_candidate_after_manual_audit")],
            s[s["stress_guarded_action"].eq("priority_review_candidate")],
            s[s["stress_guarded_action"].eq("support_review_candidate")],
            s,
        ],
        ignore_index=True,
    ).drop_duplicates("candidate_id")
    leads = ordered.head(lead_n).copy()
    leads["benchmark_label_status"] = np.where(
        pd.to_numeric(leads.get("label"), errors="coerce").eq(1),
        "known_positive_in_current_benchmark",
        "not_positive_or_unlabeled_in_current_benchmark",
    )
    leads["public_qk_independence_note"] = np.where(
        (num(leads.get("public_weight_fraction", pd.Series(index=leads.index))) <= 0.05)
        & (num(leads.get("fallback_weight_fraction", pd.Series(index=leads.index))) <= 0.08),
        "not driven by caveated public or fallback support",
        "inspect public/fallback support before claim",
    )
    leads["manual_audit_priority"] = np.select(
        [
            leads["stress_guarded_action"].eq("claim_safe_candidate_after_manual_audit"),
            leads["stress_guarded_action"].eq("priority_review_candidate"),
        ],
        ["A_claim_safe_manual_audit", "B_priority_review"],
        default="C_support_or_watchlist",
    )
    leads["claim_boundary"] = "research triage only; manual biological and overlap audit required"
    keep = [
        "candidate_id",
        "manual_audit_priority",
        "stress_guarded_rank_global",
        "stress_guarded_final_review_score",
        "stress_guarded_claim_safe_score",
        "stress_guarded_discovery_score",
        "stress_guarded_confidence",
        "stress_guarded_action",
        "benchmark_label_status",
        "source_name",
        "hla_allele_4digit",
        "peptide",
        "label",
        "stress_guarded_clean_method_count",
        "public_weight_fraction",
        "fallback_weight_fraction",
        "stress_guarded_expert_disagreement",
        "stress_guarded_reason_primary",
        "public_qk_independence_note",
        "claim_boundary",
    ]
    return leads[[c for c in keep if c in leads.columns]]


def build_support(leads: pd.DataFrame, contrib: pd.DataFrame) -> pd.DataFrame:
    if leads.empty or contrib.empty:
        return pd.DataFrame()
    c = contrib[contrib["candidate_id"].isin(leads["candidate_id"])].copy()
    if c.empty:
        return pd.DataFrame()
    c["candidate_method_contribution"] = num(c["candidate_method_contribution"])
    c["candidate_method_weight"] = num(c["candidate_method_weight"])
    c["score_for_fusion"] = num(c["score_for_fusion"])
    rows = []
    for cid, g in c.groupby("candidate_id", sort=False):
        g = g.sort_values("candidate_method_contribution", ascending=False)
        internal = g[g["method_role"].isin(["anchor", "internal_candidate"])]
        public = g[g["method_role"].eq("caveated_public_comparator")]
        fallback = g[g["method_role"].eq("bounded_fallback")]
        rows.append(
            {
                "candidate_id": cid,
                "top_internal_methods": ", ".join(internal["method_name"].astype(str).head(8)),
                "top_public_methods": ", ".join(public["method_name"].astype(str).head(5)),
                "top_fallback_methods": ", ".join(fallback["method_name"].astype(str).head(5)),
                "n_internal_methods_in_top_contrib": int(internal["method_name"].nunique()),
                "n_public_methods_in_top_contrib": int(public["method_name"].nunique()),
                "n_fallback_methods_in_top_contrib": int(fallback["method_name"].nunique()),
                "top_method": g["method_name"].iloc[0],
                "top_method_role": g["method_role"].iloc[0],
                "top_method_score": float(g["score_for_fusion"].iloc[0]),
                "top_method_weight": float(g["candidate_method_weight"].iloc[0]),
                "top_method_contribution": float(g["candidate_method_contribution"].iloc[0]),
            }
        )
    return pd.DataFrame(rows)


def build_ablation(leads: pd.DataFrame, support: pd.DataFrame) -> pd.DataFrame:
    if leads.empty:
        return pd.DataFrame()
    out = leads.merge(support, on="candidate_id", how="left") if not support.empty else leads.copy()
    out["clean_only_proxy_score"] = num(out["stress_guarded_claim_safe_score"]) * num(out["stress_guarded_confidence"], 1.0)
    out["public_removed_expected_delta"] = -num(out.get("public_weight_fraction", pd.Series(index=out.index)))
    out["fallback_removed_expected_delta"] = -num(out.get("fallback_weight_fraction", pd.Series(index=out.index)))
    out["ablation_interpretation"] = np.where(
        (num(out.get("public_weight_fraction", pd.Series(index=out.index))) <= 0.05)
        & (num(out.get("fallback_weight_fraction", pd.Series(index=out.index))) <= 0.08)
        & (num(out.get("stress_guarded_clean_method_count", pd.Series(index=out.index))) >= 10),
        "robust to public/fallback removal by design; clean internal support dominates",
        "requires manual ablation review before strong claim",
    )
    keep = [
        "candidate_id",
        "manual_audit_priority",
        "stress_guarded_final_review_score",
        "stress_guarded_claim_safe_score",
        "clean_only_proxy_score",
        "stress_guarded_clean_method_count",
        "public_weight_fraction",
        "fallback_weight_fraction",
        "top_internal_methods",
        "top_method",
        "top_method_role",
        "ablation_interpretation",
    ]
    return out[[c for c in keep if c in out.columns]]


def html_table(df: pd.DataFrame, n: int = 30) -> str:
    if df.empty:
        return "<p class='muted'>No rows available.</p>"
    view = df.head(n)
    head = "".join(f"<th>{html.escape(c)}</th>" for c in view.columns)
    rows = []
    for _, r in view.iterrows():
        rows.append("<tr>" + "".join(f"<td>{html.escape(fmt(r[c]))}</td>" for c in view.columns) + "</tr>")
    return f"<div class='table-wrap'><table><thead><tr>{head}</tr></thead><tbody>{''.join(rows)}</tbody></table></div>"


def write_reports(output_root: Path, hub_root: Path, leads: pd.DataFrame, topk: pd.DataFrame, ablation: pd.DataFrame, support: pd.DataFrame) -> None:
    md = f"""# BAR-Neo High-Impact Candidate Pack

## Purpose

Reviewer-facing evidence cards for stress-guarded BAR-Neo leads: top-k behavior, public/fallback ablation proxy, and clean internal support.

## Top-K Audit

{dataframe_to_markdown(topk, max_rows=20)}

## Lead Candidates

{dataframe_to_markdown(leads, max_rows=30)}

## Public/Fallback Ablation Proxy

{dataframe_to_markdown(ablation, max_rows=30)}

## Method Support Matrix

{dataframe_to_markdown(support, max_rows=30)}

## Claim Boundary

Allowed: high-impact research triage and manual-audit evidence. Forbidden: clinical vaccine selection, new SOTA predictor, external validation proven, public clean baseline without row-level audit, or quantum advantage.
"""
    (output_root / "BAR_NEO_HIGH_IMPACT_CANDIDATE_PACK.md").write_text(md.strip() + "\n")
    kr = f"""# BAR-Neo High-Impact Candidate Pack KR

## 한 줄 결론

stress-guarded 후보 중 manual-audit claim-safe lead를 evidence card로 고정했다. 이 패키지는 public/QK 의존이 낮은지, clean internal support가 충분한지, top-k가 benchmark positive를 잡는지 보여준다.

## Top-k / action audit

{dataframe_to_markdown(topk, max_rows=20)}

## Lead candidates

{dataframe_to_markdown(leads, max_rows=30)}

## Public/QK 제거 관점 ablation proxy

{dataframe_to_markdown(ablation, max_rows=30)}

## Internal support matrix

{dataframe_to_markdown(support, max_rows=30)}

## 해석

- claim-safe manual-audit 후보는 현재 benchmark label 기준 positive인지 확인한다.
- public/fallback weight fraction이 낮으면 public pretrained/QK가 후보를 끌어올린 것이 아니라는 방어 논리가 생긴다.
- 그래도 clinical/SOTA claim은 금지다. 다음 단계는 row-level public overlap audit, peptide/HLA/source manual audit, PAAD/THCA metadata 연결이다.
"""
    (output_root / "BAR_NEO_HIGH_IMPACT_CANDIDATE_PACK_KR.md").write_text(kr.strip() + "\n")

    n_claim = int(leads["manual_audit_priority"].eq("A_claim_safe_manual_audit").sum()) if not leads.empty else 0
    n_claim_pos = int(
        leads[leads["manual_audit_priority"].eq("A_claim_safe_manual_audit")]["benchmark_label_status"].eq("known_positive_in_current_benchmark").sum()
    ) if not leads.empty else 0
    css = "body{margin:0;background:#0d1117;color:#e6edf3;font-family:JetBrains Mono,monospace;line-height:1.55}header{padding:42px 34px;background:#101820;border-bottom:1px solid #2a3441}h1{font-family:Georgia,serif;font-size:46px;margin:0 0 12px}.lead{color:#c8d1dc;max-width:980px}.stats{display:grid;grid-template-columns:repeat(4,minmax(140px,1fr));gap:12px}.stat{border:1px solid #2a3441;background:#151b23;padding:14px}.stat b{display:block;color:#5eead4;font-size:24px}.stat span{color:#9aa7b4;font-size:12px}main{max-width:1360px;margin:auto;padding:28px}section{border-bottom:1px solid #2a3441;padding:22px 0}.table-wrap{overflow:auto;border:1px solid #2a3441;background:#101820}table{width:100%;border-collapse:collapse;font-size:12px}th,td{border-bottom:1px solid #2a3441;padding:8px 9px;text-align:left;white-space:nowrap}th{color:#e3b341;background:#111821}.muted{color:#9aa7b4}"
    page = f"""<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>BAR-Neo High-Impact Candidate Pack</title><style>{css}</style></head><body>
<header><h1>BAR-Neo High-Impact Candidate Pack</h1><p class="lead">Stress-guarded lead-candidate evidence cards. Research triage only: no clinical selection, no SOTA claim, no public-clean-baseline claim.</p>
<div class="stats"><div class="stat"><b>{len(leads):,}</b><span>Lead rows</span></div><div class="stat"><b>{n_claim:,}</b><span>Claim-safe manual-audit</span></div><div class="stat"><b>{n_claim_pos:,}</b><span>Benchmark positives among claim-safe</span></div><div class="stat"><b>0</b><span>Public clean baseline claims</span></div></div></header>
<main><section><h2>Top-k audit</h2>{html_table(topk, 20)}</section><section><h2>Lead candidates</h2>{html_table(leads, 30)}</section><section><h2>Ablation proxy</h2>{html_table(ablation, 30)}</section><section><h2>Method support matrix</h2>{html_table(support, 30)}</section><p class="muted">Output root: {html.escape(str(output_root))}</p></main></body></html>"""
    (hub_root / PAGE_NAME).write_text(page)


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
    scores = read_tsv(output_root / "barneo_stress_guarded_candidate_scores.tsv")
    contrib = read_tsv(output_root / "barneo_stress_guarded_expert_contributions.tsv")
    topk = build_topk(scores)
    leads = build_leads(scores, args.lead_n)
    support = build_support(leads, contrib)
    ablation = build_ablation(leads, support)

    write_tsv(leads, output_root / "barneo_high_impact_lead_candidates.tsv")
    write_tsv(topk, output_root / "barneo_high_impact_topk_audit.tsv")
    write_tsv(ablation, output_root / "barneo_high_impact_ablation_audit.tsv")
    write_tsv(support, output_root / "barneo_high_impact_method_support_matrix.tsv")
    write_reports(output_root, hub_root, leads, topk, ablation, support)

    claim = leads[leads["manual_audit_priority"].eq("A_claim_safe_manual_audit")] if not leads.empty else pd.DataFrame()
    manifest_path = output_root / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    manifest.setdefault("output_files", [])
    for name in OUTPUTS:
        if name not in manifest["output_files"]:
            manifest["output_files"].append(name)
    manifest.setdefault("summary", {})
    manifest["summary"].update(
        {
            "n_high_impact_leads": int(len(leads)),
            "n_high_impact_claim_safe_manual_audit": int(len(claim)),
            "n_high_impact_claim_safe_benchmark_positives": int(claim["benchmark_label_status"].eq("known_positive_in_current_benchmark").sum()) if not claim.empty else 0,
            "high_impact_html": str(hub_root / PAGE_NAME),
        }
    )
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    update_manifest(
        output_root,
        "barneo_high_impact_candidate_pack",
        {"outputs": OUTPUTS, "n_leads": int(len(leads)), "warnings": ["High-impact pack is research triage and manual-audit support only."]},
    )
    print(
        "[barneo-high-impact] "
        f"leads={len(leads)} claim_safe={len(claim)} "
        f"claim_safe_benchmark_pos={int(claim['benchmark_label_status'].eq('known_positive_in_current_benchmark').sum()) if not claim.empty else 0}"
    )


if __name__ == "__main__":
    main()
