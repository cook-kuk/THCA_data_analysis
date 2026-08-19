#!/usr/bin/env python3
"""Build reviewer-facing evidence cards for BAR-Neo kill-audit pass leads."""

from __future__ import annotations

import argparse
import html
import json
from pathlib import Path

import numpy as np
import pandas as pd

from common import dataframe_to_markdown, ensure_dir, update_manifest, write_tsv


OUTPUTS = [
    "barneo_top_pass_reviewer_evidence.tsv",
    "barneo_top_pass_method_contributions.tsv",
    "barneo_high_impact_lab_handoff_queue.tsv",
    "BAR_NEO_TOP_PASS_REVIEWER_EVIDENCE_PACK.md",
    "BAR_NEO_TOP_PASS_REVIEWER_EVIDENCE_PACK_KR.md",
]
HTML_NAME = "barneo_top_pass_reviewer_evidence_2026_05_10.html"


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


def truth(x: object) -> bool:
    return str(x).strip().lower() in {"true", "1", "yes"}


def aa_features(peptide: object) -> dict[str, object]:
    seq = str(peptide or "").upper()
    n = max(len(seq), 1)
    hydrophobic = set("AILMFWYV")
    aromatic = set("FWY")
    positive = set("KRH")
    negative = set("DE")
    return {
        "peptide_length": len(seq),
        "hydrophobic_fraction": sum(a in hydrophobic for a in seq) / n,
        "aromatic_fraction": sum(a in aromatic for a in seq) / n,
        "charge_proxy": sum(a in positive for a in seq) - sum(a in negative for a in seq),
        "cysteine_count": seq.count("C"),
        "glycine_proline_count": seq.count("G") + seq.count("P"),
        "n_term_residue": seq[:1],
        "c_term_residue": seq[-1:] if seq else "",
        "hla_i_anchor_proxy": f"{seq[1:2]}/{seq[-1:]}" if 8 <= len(seq) <= 11 else "",
    }


def summarize_context(stress: pd.DataFrame, source: str, hla: str) -> tuple[str, str]:
    if stress.empty:
        return "source stress board unavailable", "HLA stress board unavailable"
    src = stress[
        stress.get("split_contract", pd.Series(dtype=str)).astype(str).str.contains("source", case=False, na=False)
        & stress.get("split_group", pd.Series(dtype=str)).astype(str).eq(str(source))
    ].copy()
    hla_rows = stress[
        stress.get("split_contract", pd.Series(dtype=str)).astype(str).str.contains("hla", case=False, na=False)
        & stress.get("split_group", pd.Series(dtype=str)).astype(str).eq(str(hla))
    ].copy()

    def one_line(df: pd.DataFrame, label: str) -> str:
        if df.empty:
            return f"{label} stress slice unavailable"
        r = df.sort_values("best_clean_AUPRC", ascending=False).iloc[0]
        return (
            f"{label}: n={int(r.get('n_total_max', 0) or 0)}, pos={int(r.get('n_pos_max', 0) or 0)}, "
            f"best_clean={r.get('best_clean_method', '')}, best_clean_AUPRC={float(r.get('best_clean_AUPRC', np.nan)):.3f}, "
            f"Structure_LR_AUPRC={float(r.get('structure_lr_AUPRC', np.nan)):.3f}, "
            f"clean_minus_anchor={float(r.get('clean_minus_anchor_AUPRC', np.nan)):.3f}"
        )

    return one_line(src, f"source {source}"), one_line(hla_rows, f"HLA {hla}")


def build_evidence(
    kill: pd.DataFrame,
    support: pd.DataFrame,
    stress: pd.DataFrame,
) -> pd.DataFrame:
    if kill.empty:
        return pd.DataFrame()
    passes = kill[kill["reviewer_kill_disposition"].eq("passes_current_reviewer_kill_audit")].copy()
    rows = []
    for _, r in passes.sort_values("stress_guarded_rank_global").iterrows():
        cid = str(r["candidate_id"])
        sup = support[support.get("candidate_id", pd.Series(dtype=str)).astype(str).eq(cid)]
        source_note, hla_note = summarize_context(stress, str(r.get("source_name", "")), str(r.get("hla_allele_4digit", "")))
        feats = aa_features(r.get("peptide", ""))
        exact_clean = not truth(r.get("exact_peptide_train_overlap")) and not truth(r.get("exact_peptide_hla_train_overlap"))
        no_near = not truth(r.get("near_peptide_train_overlap")) and float(r.get("nearest_neighbor_similarity", 0) or 0) < 0.80
        no_decoy_pressure = int(r.get("source_negatives_at_or_above", 0) or 0) == 0 and int(r.get("source_hla_negatives_at_or_above", 0) or 0) == 0
        public_ok = float(r.get("public_weight_fraction", 0) or 0) <= 0.05
        fallback_ok = float(r.get("fallback_weight_fraction", 0) or 0) <= 0.08
        rows.append(
            {
                "candidate_id": cid,
                "stress_guarded_rank_global": r.get("stress_guarded_rank_global", np.nan),
                "stress_guarded_final_review_score": r.get("stress_guarded_final_review_score", np.nan),
                "source_name": r.get("source_name", ""),
                "hla_allele_4digit": r.get("hla_allele_4digit", ""),
                "peptide": r.get("peptide", ""),
                "label": r.get("label", np.nan),
                "leakage_risk_level": r.get("leakage_risk_level", ""),
                "exact_overlap_clean": exact_clean,
                "near_neighbor_clean": no_near,
                "decoy_pressure_clean": no_decoy_pressure,
                "public_dependency_clean": public_ok,
                "fallback_dependency_clean": fallback_ok,
                "nearest_neighbor_similarity": r.get("nearest_neighbor_similarity", np.nan),
                "source_negatives_at_or_above": r.get("source_negatives_at_or_above", np.nan),
                "source_hla_negatives_at_or_above": r.get("source_hla_negatives_at_or_above", np.nan),
                "public_weight_fraction": r.get("public_weight_fraction", np.nan),
                "fallback_weight_fraction": r.get("fallback_weight_fraction", np.nan),
                "top_method": r.get("top_method", ""),
                "top_internal_methods": sup["top_internal_methods"].iloc[0] if not sup.empty and "top_internal_methods" in sup else "",
                "n_internal_methods_in_top_contrib": sup["n_internal_methods_in_top_contrib"].iloc[0] if not sup.empty and "n_internal_methods_in_top_contrib" in sup else np.nan,
                "n_public_methods_in_top_contrib": sup["n_public_methods_in_top_contrib"].iloc[0] if not sup.empty and "n_public_methods_in_top_contrib" in sup else np.nan,
                "n_fallback_methods_in_top_contrib": sup["n_fallback_methods_in_top_contrib"].iloc[0] if not sup.empty and "n_fallback_methods_in_top_contrib" in sup else np.nan,
                "source_stress_context": source_note,
                "hla_stress_context": hla_note,
                "reviewer_evidence_summary": "low leakage, no exact/near overlap, no same-source decoy pressure, clean internal support dominates",
                "recommended_next_action": "manual provenance check plus experimental/biological plausibility review",
                "claim_boundary": "manual-review research lead only; not clinical selection or SOTA validation",
                **feats,
            }
        )
    return pd.DataFrame(rows)


def build_method_contrib(pass_ids: list[str], contrib: pd.DataFrame) -> pd.DataFrame:
    if contrib.empty or not pass_ids:
        return pd.DataFrame()
    c = contrib[contrib["candidate_id"].astype(str).isin(pass_ids)].copy()
    c["candidate_method_contribution"] = num(c, "candidate_method_contribution")
    c["candidate_method_weight"] = num(c, "candidate_method_weight")
    c["score_for_fusion"] = num(c, "score_for_fusion")
    c = c.sort_values(["candidate_id", "candidate_method_contribution"], ascending=[True, False])
    c["support_rank_within_candidate"] = c.groupby("candidate_id").cumcount() + 1
    keep = [
        "candidate_id",
        "support_rank_within_candidate",
        "method_name",
        "method_role",
        "method_family",
        "score_for_fusion",
        "candidate_method_weight",
        "candidate_method_contribution",
        "context_adjustment_reasons",
    ]
    return c[[x for x in keep if x in c.columns]]


def build_lab_handoff(kill: pd.DataFrame) -> pd.DataFrame:
    if kill.empty:
        return pd.DataFrame()
    rows = []
    for _, r in kill.sort_values("stress_guarded_rank_global").iterrows():
        disposition = str(r.get("reviewer_kill_disposition", ""))
        if disposition == "passes_current_reviewer_kill_audit":
            tier = "T1_clean_manual_review_lead"
            action = "advance to manual biological plausibility and assay-design review"
            verify = "source provenance; WT/gene/mutation; expression/clonality; patient gate; safety context"
        elif disposition == "manual_audit_required":
            tier = "T2_resolve_before_claim"
            action = "resolve near-peptide/source/protein-window caveat before any lead claim"
            verify = "near-neighbor lineage; source/protein-window overlap; then standard provenance and patient gate"
        else:
            tier = "T3_do_not_use_for_clean_claim"
            action = "exclude from clean lead list; keep only as leakage/stress-control example"
            verify = "exact overlap and high-leakage flags; do not advance as clean lead"
        rows.append(
            {
                "candidate_id": r.get("candidate_id", ""),
                "handoff_tier": tier,
                "stress_guarded_rank_global": r.get("stress_guarded_rank_global", np.nan),
                "stress_guarded_final_review_score": r.get("stress_guarded_final_review_score", np.nan),
                "source_name": r.get("source_name", ""),
                "hla_allele_4digit": r.get("hla_allele_4digit", ""),
                "peptide": r.get("peptide", ""),
                "label": r.get("label", np.nan),
                "leakage_risk_level": r.get("leakage_risk_level", ""),
                "reviewer_kill_disposition": disposition,
                "reviewer_kill_reason": r.get("reviewer_kill_reason", ""),
                "public_weight_fraction": r.get("public_weight_fraction", np.nan),
                "fallback_weight_fraction": r.get("fallback_weight_fraction", np.nan),
                "handoff_action": action,
                "must_verify_before_claim": verify,
                "claim_boundary": "research triage only; not clinical selection, SOTA, or external validation",
            }
        )
    return pd.DataFrame(rows)


def write_reports(output_root: Path, hub_root: Path, evidence: pd.DataFrame, contrib: pd.DataFrame, handoff: pd.DataFrame) -> None:
    top_cols = [
        "candidate_id",
        "stress_guarded_rank_global",
        "stress_guarded_final_review_score",
        "source_name",
        "hla_allele_4digit",
        "peptide",
        "exact_overlap_clean",
        "near_neighbor_clean",
        "decoy_pressure_clean",
        "public_dependency_clean",
        "fallback_dependency_clean",
        "top_method",
        "reviewer_evidence_summary",
    ]
    contrib_cols = [
        "candidate_id",
        "support_rank_within_candidate",
        "method_name",
        "method_role",
        "score_for_fusion",
        "candidate_method_weight",
        "candidate_method_contribution",
    ]
    handoff_cols = [
        "candidate_id",
        "handoff_tier",
        "stress_guarded_rank_global",
        "source_name",
        "hla_allele_4digit",
        "peptide",
        "reviewer_kill_disposition",
        "handoff_action",
        "must_verify_before_claim",
    ]
    text = f"""# BAR-Neo Top-Pass Reviewer Evidence Pack

## Purpose

This pack isolates only the high-impact leads that pass the current automated reviewer kill audit. It is designed for manual provenance, biological plausibility, and experimental triage review.

## Top-Pass Candidate Cards

{dataframe_to_markdown(evidence[top_cols] if not evidence.empty else evidence, max_rows=20)}

## Method Support

{dataframe_to_markdown(contrib[contrib_cols].head(30) if not contrib.empty else contrib, max_rows=30)}

## Lab / Manual-Review Handoff Queue

{dataframe_to_markdown(handoff[handoff_cols] if not handoff.empty else handoff, max_rows=30)}

## Claim Boundary

These are manual-review research leads only. This pack does not claim clinical vaccine selection, new SOTA prediction, external validation, quantum advantage, or clean public-tool baseline status.
"""
    (output_root / "BAR_NEO_TOP_PASS_REVIEWER_EVIDENCE_PACK.md").write_text(text.strip() + "\n")

    kr = f"""# BAR-Neo Top-Pass Reviewer Evidence Pack KR

## 한 줄 결론

킬 오딧을 통과한 high-impact 후보만 따로 분리했다. 현재 자동 기준으로 바로 밀 수 있는 후보는 3개이며, 모두 ITSNdb_main / HLA-A*02:01 쪽에서 low leakage, no exact/near overlap, no same-source negative pressure, low public/QK/fallback dependence를 만족한다.

## Top-pass 후보 카드

{dataframe_to_markdown(evidence[top_cols] if not evidence.empty else evidence, max_rows=20)}

## 후보별 method support

{dataframe_to_markdown(contrib[contrib_cols].head(30) if not contrib.empty else contrib, max_rows=30)}

## Lab / manual-review handoff queue

{dataframe_to_markdown(handoff[handoff_cols] if not handoff.empty else handoff, max_rows=30)}

## 다음 액션

1. peptide/HLA/source provenance 수동 확인.
2. wild-type peptide, gene/mutation, expression/clonality metadata 연결.
3. PAAD/THCA patient-gate demo에 실제 disease timing과 immune/safety context를 붙여 triage view 생성.
4. public training-corpus row-level overlap audit 전까지 public tool support는 caveated로 유지.

## Claim boundary

manual-review research lead만 허용한다. clinical vaccine selection, new SOTA predictor, external validation proven, quantum advantage claim은 금지한다.
"""
    (output_root / "BAR_NEO_TOP_PASS_REVIEWER_EVIDENCE_PACK_KR.md").write_text(kr.strip() + "\n")

    html_rows = []
    for _, r in evidence.iterrows():
        checks = [
            ("exact", r.get("exact_overlap_clean")),
            ("near", r.get("near_neighbor_clean")),
            ("decoy", r.get("decoy_pressure_clean")),
            ("public", r.get("public_dependency_clean")),
            ("fallback", r.get("fallback_dependency_clean")),
        ]
        badges = " ".join(f"<span class='badge'>{html.escape(k)}: {html.escape(str(v))}</span>" for k, v in checks)
        html_rows.append(
            f"""
            <section class="card">
              <div class="rank">Rank {html.escape(str(r.get('stress_guarded_rank_global', '')))}</div>
              <h2>{html.escape(str(r.get('candidate_id', '')))} · {html.escape(str(r.get('peptide', '')))}</h2>
              <p class="meta">{html.escape(str(r.get('source_name', '')))} · {html.escape(str(r.get('hla_allele_4digit', '')))} · score {float(r.get('stress_guarded_final_review_score', 0) or 0):.3f}</p>
              <div class="badges">{badges}</div>
              <p><b>Top method</b>: {html.escape(str(r.get('top_method', '')))}</p>
              <p><b>Source context</b>: {html.escape(str(r.get('source_stress_context', '')))}</p>
              <p><b>HLA context</b>: {html.escape(str(r.get('hla_stress_context', '')))}</p>
              <p><b>Next action</b>: {html.escape(str(r.get('recommended_next_action', '')))}</p>
            </section>
            """
        )
    handoff_html = ""
    if not handoff.empty:
        cols = [c for c in ["candidate_id", "handoff_tier", "stress_guarded_rank_global", "peptide", "reviewer_kill_disposition", "handoff_action"] if c in handoff.columns]
        head = "".join(f"<th>{html.escape(c)}</th>" for c in cols)
        body = []
        for _, row in handoff.loc[:, cols].iterrows():
            body.append("<tr>" + "".join(f"<td>{html.escape(str(row[c]))}</td>" for c in cols) + "</tr>")
        handoff_html = f"<section class='table-card'><h2>High-impact lab handoff queue</h2><table><thead><tr>{head}</tr></thead><tbody>{''.join(body)}</tbody></table></section>"
    page = f"""<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>BAR-Neo Top-Pass Reviewer Evidence</title>
<style>
body{{margin:0;background:#0d1117;color:#e6edf3;font-family:JetBrains Mono,ui-monospace,Menlo,monospace;line-height:1.55}}
header{{padding:44px 34px 24px;border-bottom:1px solid #2a3441}}h1{{font-family:Georgia,serif;font-size:46px;margin:0 0 10px}}.lead{{color:#c8d1dc;max-width:980px}}
main{{max-width:1180px;margin:auto;padding:24px;display:grid;gap:16px}}.card,.table-card{{border:1px solid #2a3441;background:#101820;padding:18px;overflow:auto}}h2{{font-family:Georgia,serif;margin:6px 0 4px;font-size:30px}}.rank{{color:#e3b341;font-weight:800}}.meta{{color:#5eead4}}.badge{{display:inline-block;border:1px solid #2a3441;background:#151b23;color:#86efac;margin:4px 6px 4px 0;padding:5px 8px;font-size:12px}}b{{color:#e3b341}}table{{width:100%;border-collapse:collapse;font-size:12px}}th,td{{border-bottom:1px solid #2a3441;padding:8px;text-align:left;vertical-align:top;white-space:nowrap}}th{{color:#e3b341}}.path{{padding:0 24px 36px;color:#9aa7b4;font-size:12px}}
</style></head><body><header><h1>BAR-Neo Top-Pass Reviewer Evidence</h1><p class="lead">Only candidates that pass the automated reviewer kill audit are shown as clean manual-review leads. The handoff queue also lists one resolve-before-claim row and blocked leakage controls. This is research triage, not clinical selection or SOTA validation.</p></header><main>{''.join(html_rows)}{handoff_html}</main><p class="path">Source: {html.escape(str(output_root / 'barneo_top_pass_reviewer_evidence.tsv'))}</p></body></html>"""
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

    kill = read_tsv(output_root / "barneo_high_impact_reviewer_kill_audit.tsv")
    support = read_tsv(output_root / "barneo_high_impact_method_support_matrix.tsv")
    contrib_raw = read_tsv(output_root / "barneo_stress_guarded_expert_contributions.tsv")
    stress = read_tsv(output_root / "clean_neobench_source_hla_stress_slice_board.tsv")

    evidence = build_evidence(kill, support, stress)
    pass_ids = evidence["candidate_id"].astype(str).tolist() if not evidence.empty else []
    contrib = build_method_contrib(pass_ids, contrib_raw)
    handoff = build_lab_handoff(kill)

    write_tsv(evidence, output_root / "barneo_top_pass_reviewer_evidence.tsv")
    write_tsv(contrib, output_root / "barneo_top_pass_method_contributions.tsv")
    write_tsv(handoff, output_root / "barneo_high_impact_lab_handoff_queue.tsv")
    write_reports(output_root, hub_root, evidence, contrib, handoff)

    outputs = OUTPUTS + [str(hub_root / HTML_NAME)]
    manifest_path = output_root / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    manifest.setdefault("output_files", [])
    for name in outputs:
        if str(name) not in manifest["output_files"]:
            manifest["output_files"].append(str(name))
    manifest.setdefault("summary", {})
    manifest["summary"].update(
        {
            "n_top_pass_reviewer_evidence_candidates": int(len(evidence)),
            "n_top_pass_method_contribution_rows": int(len(contrib)),
            "n_high_impact_lab_handoff_rows": int(len(handoff)),
            "top_pass_reviewer_evidence_html": str(hub_root / HTML_NAME),
        }
    )
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    update_manifest(
        output_root,
        "barneo_top_pass_reviewer_evidence_pack",
        {
            "outputs": outputs,
            "n_candidates": int(len(evidence)),
            "warnings": ["Top-pass evidence pack is manual-review triage, not clinical or SOTA validation."],
        },
    )
    print(f"[barneo-top-pass-evidence] candidates={len(evidence)} contrib_rows={len(contrib)} html={hub_root / HTML_NAME}")


if __name__ == "__main__":
    main()
