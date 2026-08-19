#!/usr/bin/env python3
"""Build translational-readiness evidence pack for BAR-Neo T1 leads."""

from __future__ import annotations

import argparse
import html
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from common import dataframe_to_markdown, ensure_dir, normalize_empty, update_manifest, write_tsv


HTML_NAME = "barneo_t1_translational_readiness_2026_05_10.html"
OUTPUTS = [
    "barneo_t1_translational_readiness.tsv",
    "barneo_t1_metadata_gap_matrix.tsv",
    "barneo_t1_assay_handoff.tsv",
    "BAR_NEO_T1_TRANSLATIONAL_READINESS_PACK.md",
    "BAR_NEO_T1_TRANSLATIONAL_READINESS_PACK_KR.md",
]

REQUIRED_FIELDS = [
    ("gene", "antigen_identity", "needed to connect peptide to tumor biology and mutation provenance"),
    ("mutation_id", "antigen_identity", "needed to verify mutant origin and avoid peptide-only claims"),
    ("wt_peptide", "antigen_identity", "needed for mutant-vs-wildtype specificity review"),
    ("source_protein_window", "provenance", "needed for source/protein-window overlap and manual provenance audit"),
    ("patient_id", "patient_gate", "needed for patient-level split/gate and leakage control"),
    ("cancer_type", "patient_gate", "needed for PAAD/THCA disease-context triage"),
    ("disease_context", "patient_gate", "needed to decide MRD/resected/high-risk vs low-priority context"),
    ("expression_tpm", "antigen_evidence", "needed for antigen expression support"),
    ("mutant_expression", "antigen_evidence", "needed for mutant allele expression support"),
    ("vaf", "antigen_evidence", "needed for clonality/subclonality review"),
    ("clonality", "antigen_evidence", "needed for clonal antigen prioritization"),
    ("hla_loh", "presentation_safety", "needed for presentation intactness review"),
    ("b2m_status", "presentation_safety", "needed for antigen presentation intactness review"),
    ("antigen_processing_status", "presentation_safety", "needed for processing/presentation review"),
    ("immune_context_score", "immune_context", "needed for vaccine-response plausibility"),
    ("ifng_score", "immune_context", "needed for inflamed/excluded immune-context review"),
    ("cytolytic_score", "immune_context", "needed for immune effector context review"),
    ("treatment_context", "patient_gate", "needed for research-triage timing and combination context"),
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo-root", default=".", help="Repository root")
    p.add_argument("--output-root", required=True, help="CLEAN-NeoBench output root")
    p.add_argument("--hub-root", default="project/papers_hub_2026_05_04", help="HTML hub root")
    return p.parse_args()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def truth(value: Any) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def present(value: Any) -> bool:
    return normalize_empty(value) != ""


def value_of(row: pd.Series, col: str) -> Any:
    return row[col] if col in row.index else np.nan


def build_gap_matrix(t1: pd.DataFrame, master: pd.DataFrame) -> pd.DataFrame:
    rows = []
    if t1.empty:
        return pd.DataFrame()
    merged = t1.merge(master, on="candidate_id", how="left", suffixes=("", "_master"))
    for _, r in merged.iterrows():
        for field, group, why in REQUIRED_FIELDS:
            v = value_of(r, field)
            rows.append(
                {
                    "candidate_id": r.get("candidate_id", ""),
                    "peptide": r.get("peptide", r.get("peptide_master", "")),
                    "hla_allele_4digit": r.get("hla_allele_4digit", r.get("hla_allele_4digit_master", "")),
                    "field": field,
                    "gate_group": group,
                    "metadata_present": present(v),
                    "current_value": normalize_empty(v),
                    "why_needed": why,
                    "claim_if_missing": "benchmark/manual-review only; no patient/translational lead claim",
                }
            )
    return pd.DataFrame(rows)


def build_readiness(t1: pd.DataFrame, master: pd.DataFrame, gap: pd.DataFrame) -> pd.DataFrame:
    if t1.empty:
        return pd.DataFrame()
    merged = t1.merge(master, on="candidate_id", how="left", suffixes=("", "_master"))
    rows = []
    for _, r in merged.sort_values("stress_guarded_rank_global").iterrows():
        cid = str(r.get("candidate_id", ""))
        peptide = normalize_empty(r.get("peptide", r.get("peptide_master", "")))
        hla = normalize_empty(r.get("hla_allele_4digit", r.get("hla_allele_4digit_master", "")))
        mhc_class = normalize_empty(r.get("mhc_class", r.get("mhc_class_master", "")))
        peptide_len = int(float(r.get("peptide_length_master", r.get("peptide_length", len(peptide))) or len(peptide) or 0))
        benchmark_checks = [
            truth(r.get("exact_overlap_clean")),
            truth(r.get("near_neighbor_clean")),
            truth(r.get("decoy_pressure_clean")),
            truth(r.get("public_dependency_clean")),
            truth(r.get("fallback_dependency_clean")),
        ]
        benchmark_evidence_score = sum(benchmark_checks) / max(len(benchmark_checks), 1)
        peptide_hla_assay_ready = bool(peptide and hla and mhc_class == "I" and 8 <= peptide_len <= 11 and benchmark_evidence_score >= 0.95)
        assay_design_score = 1.0 if peptide_hla_assay_ready else 0.35 if peptide and hla else 0.0

        g = gap[gap["candidate_id"].astype(str).eq(cid)].copy()
        by_group = g.groupby("gate_group")["metadata_present"].mean().to_dict() if not g.empty else {}
        antigen_identity_score = float(by_group.get("antigen_identity", 0.0))
        antigen_evidence_score = float(by_group.get("antigen_evidence", 0.0))
        patient_gate_score = float(by_group.get("patient_gate", 0.0))
        presentation_safety_score = float(by_group.get("presentation_safety", 0.0))
        immune_context_score = float(by_group.get("immune_context", 0.0))
        translational_metadata_score = (
            0.28 * antigen_identity_score
            + 0.24 * antigen_evidence_score
            + 0.18 * patient_gate_score
            + 0.18 * presentation_safety_score
            + 0.12 * immune_context_score
        )
        overall = 0.35 * benchmark_evidence_score + 0.25 * assay_design_score + 0.40 * translational_metadata_score
        missing = g[~g["metadata_present"]]["field"].astype(str).tolist() if not g.empty else []
        if peptide_hla_assay_ready and translational_metadata_score < 0.35:
            tier = "assay_design_ready_metadata_blocked"
            action = "design pHLA-level review/assay plan, but block patient/translational claim until metadata are linked"
        elif peptide_hla_assay_ready:
            tier = "translational_review_ready"
            action = "advance to integrated biological and patient-gated review"
        else:
            tier = "not_assay_ready"
            action = "resolve peptide/HLA/class/benchmark evidence before assay planning"
        rows.append(
            {
                "candidate_id": cid,
                "stress_guarded_rank_global": r.get("stress_guarded_rank_global", np.nan),
                "stress_guarded_final_review_score": r.get("stress_guarded_final_review_score", np.nan),
                "source_name": r.get("source_name", ""),
                "hla_allele_4digit": hla,
                "peptide": peptide,
                "mhc_class": mhc_class,
                "peptide_length": peptide_len,
                "benchmark_evidence_score": benchmark_evidence_score,
                "peptide_hla_assay_design_ready": peptide_hla_assay_ready,
                "assay_design_readiness_score": assay_design_score,
                "antigen_identity_score": antigen_identity_score,
                "antigen_evidence_score": antigen_evidence_score,
                "patient_gate_score": patient_gate_score,
                "presentation_safety_score": presentation_safety_score,
                "immune_context_score": immune_context_score,
                "translational_metadata_score": translational_metadata_score,
                "impact_readiness_score": overall,
                "impact_readiness_tier": tier,
                "missing_high_impact_metadata": ", ".join(missing[:18]),
                "recommended_next_action": action,
                "claim_boundary": "manual-review research triage; no clinical/SOTA/external-validation claim",
            }
        )
    return pd.DataFrame(rows)


def build_assay_handoff(readiness: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, r in readiness.sort_values("stress_guarded_rank_global").iterrows():
        anchor_proxy = ""
        seq = str(r.get("peptide", ""))
        if 8 <= len(seq) <= 11:
            anchor_proxy = f"P2={seq[1:2]}; C-term={seq[-1:]}"
        rows.append(
            {
                "candidate_id": r.get("candidate_id", ""),
                "assay_handoff_tier": "pHLA_assay_design_candidate" if truth(r.get("peptide_hla_assay_design_ready")) else "hold",
                "peptide": r.get("peptide", ""),
                "hla_allele_4digit": r.get("hla_allele_4digit", ""),
                "peptide_length": r.get("peptide_length", np.nan),
                "anchor_proxy": anchor_proxy,
                "suggested_manual_checks": "synthesize peptide feasibility; verify HLA binding/presentation; add WT peptide comparator; link gene/mutation/expression before patient claim",
                "minimum_next_data_packet": "WT peptide, gene, mutation_id, expression/mutant_expression, VAF/clonality, patient disease context, HLA-LOH/B2M",
                "claim_boundary": "pHLA assay-design candidate only until patient/translational metadata are linked",
            }
        )
    return pd.DataFrame(rows)


def write_reports(output_root: Path, hub_root: Path, readiness: pd.DataFrame, gap: pd.DataFrame, assay: pd.DataFrame) -> None:
    readiness_cols = [
        "candidate_id",
        "stress_guarded_rank_global",
        "stress_guarded_final_review_score",
        "peptide",
        "hla_allele_4digit",
        "benchmark_evidence_score",
        "peptide_hla_assay_design_ready",
        "translational_metadata_score",
        "impact_readiness_tier",
        "recommended_next_action",
    ]
    assay_cols = [
        "candidate_id",
        "assay_handoff_tier",
        "peptide",
        "hla_allele_4digit",
        "anchor_proxy",
        "minimum_next_data_packet",
    ]
    gap_summary = (
        gap.groupby(["gate_group", "field"])["metadata_present"].mean().reset_index().rename(columns={"metadata_present": "present_fraction"})
        if not gap.empty
        else pd.DataFrame()
    )
    md = f"""# BAR-Neo T1 Translational Readiness Pack

## Purpose

This pack separates pHLA assay-design readiness from patient/translational readiness for the BAR-Neo T1 leads that passed the reviewer kill audit.

## T1 Readiness

{dataframe_to_markdown(readiness[readiness_cols] if not readiness.empty else readiness, max_rows=20)}

## Assay Handoff

{dataframe_to_markdown(assay[assay_cols] if not assay.empty else assay, max_rows=20)}

## Metadata Gap Summary

{dataframe_to_markdown(gap_summary, max_rows=80)}

## Claim Boundary

The T1 rows are pHLA assay-design/manual-review candidates. They are not clinical vaccine selections, not new SOTA validation, and not patient-level translational claims until WT/gene/mutation/expression/clonality/patient/safety metadata are linked.
"""
    (output_root / "BAR_NEO_T1_TRANSLATIONAL_READINESS_PACK.md").write_text(md.strip() + "\n")

    kr = f"""# BAR-Neo T1 Translational Readiness Pack KR

## 한 줄 결론

T1 3개는 benchmark/reviewer-kill 기준에서는 강하다. 다만 현재 metadata로는 **pHLA assay-design 후보**까지가 안전하고, patient/translational lead claim은 WT peptide, gene/mutation, expression/clonality, patient context, presentation/safety metadata가 붙기 전까지 막아야 한다.

## T1 readiness table

{dataframe_to_markdown(readiness[readiness_cols] if not readiness.empty else readiness, max_rows=20)}

## Assay handoff

{dataframe_to_markdown(assay[assay_cols] if not assay.empty else assay, max_rows=20)}

## Metadata gap summary

{dataframe_to_markdown(gap_summary, max_rows=80)}

## 다음 액션

1. T1 3개에 WT peptide, gene, mutation_id, source protein/window를 붙인다.
2. expression/mutant_expression, VAF/clonality를 붙여 antigen evidence gate를 연다.
3. PAAD/THCA patient disease context, HLA-LOH/B2M, immune/safety context를 붙여 patient gate를 연다.
4. 그 전까지는 pHLA assay-design/manual-review 후보라고만 표현한다.

## Claim boundary

clinical vaccine selection, new SOTA predictor, external validation proven, quantum advantage claim은 금지한다.
"""
    (output_root / "BAR_NEO_T1_TRANSLATIONAL_READINESS_PACK_KR.md").write_text(kr.strip() + "\n")

    cards = []
    for _, r in readiness.iterrows():
        cards.append(
            f"""
            <section class="card">
              <div class="rank">T1 rank {html.escape(str(r.get('stress_guarded_rank_global', '')))}</div>
              <h2>{html.escape(str(r.get('candidate_id', '')))} · {html.escape(str(r.get('peptide', '')))}</h2>
              <p class="meta">{html.escape(str(r.get('hla_allele_4digit', '')))} · benchmark {float(r.get('benchmark_evidence_score', 0) or 0):.2f} · assay-ready {html.escape(str(r.get('peptide_hla_assay_design_ready', '')))}</p>
              <p><b>Tier:</b> {html.escape(str(r.get('impact_readiness_tier', '')))}</p>
              <p><b>Missing metadata:</b> {html.escape(str(r.get('missing_high_impact_metadata', '')))}</p>
              <p><b>Action:</b> {html.escape(str(r.get('recommended_next_action', '')))}</p>
            </section>
            """
        )
    html_text = f"""<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>BAR-Neo T1 Translational Readiness</title>
<style>
body{{margin:0;background:#0d1117;color:#e6edf3;font-family:JetBrains Mono,ui-monospace,Menlo,monospace;line-height:1.55}}
header{{padding:44px 34px 24px;border-bottom:1px solid #2a3441}}h1{{font-family:Georgia,serif;font-size:46px;margin:0 0 10px}}.lead{{color:#c8d1dc;max-width:980px}}
main{{max-width:1180px;margin:auto;padding:24px;display:grid;gap:16px}}.card{{border:1px solid #2a3441;background:#101820;padding:18px}}h2{{font-family:Georgia,serif;margin:6px 0 4px;font-size:30px}}.rank{{color:#e3b341;font-weight:800}}.meta{{color:#5eead4}}b{{color:#e3b341}}.path{{padding:0 24px 36px;color:#9aa7b4;font-size:12px}}
</style></head><body><header><h1>BAR-Neo T1 Translational Readiness</h1><p class="lead">T1 candidates are benchmark-clean manual-review leads and pHLA assay-design candidates. Patient/translational claims remain blocked until metadata gaps are filled.</p></header><main>{''.join(cards)}</main><p class="path">Source: {html.escape(str(output_root / 'barneo_t1_translational_readiness.tsv'))}</p></body></html>"""
    hub_root.mkdir(parents=True, exist_ok=True)
    (hub_root / HTML_NAME).write_text(html_text)


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

    top_pass = read_tsv(output_root / "barneo_top_pass_reviewer_evidence.tsv")
    master = read_tsv(output_root / "clean_neobench_master.tsv")
    if top_pass.empty:
        readiness = pd.DataFrame()
        gap = pd.DataFrame()
        assay = pd.DataFrame()
    else:
        gap = build_gap_matrix(top_pass, master)
        readiness = build_readiness(top_pass, master, gap)
        assay = build_assay_handoff(readiness)

    write_tsv(readiness, output_root / "barneo_t1_translational_readiness.tsv")
    write_tsv(gap, output_root / "barneo_t1_metadata_gap_matrix.tsv")
    write_tsv(assay, output_root / "barneo_t1_assay_handoff.tsv")
    write_reports(output_root, hub_root, readiness, gap, assay)

    manifest_path = output_root / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    outputs = OUTPUTS + [str(hub_root / HTML_NAME)]
    manifest.setdefault("output_files", [])
    for out in outputs:
        if str(out) not in manifest["output_files"]:
            manifest["output_files"].append(str(out))
    manifest.setdefault("summary", {})
    manifest["summary"].update(
        {
            "n_t1_translational_readiness_candidates": int(len(readiness)),
            "n_t1_phla_assay_design_ready": int(readiness["peptide_hla_assay_design_ready"].astype(bool).sum()) if not readiness.empty else 0,
            "n_t1_translational_review_ready": int(readiness["impact_readiness_tier"].eq("translational_review_ready").sum()) if not readiness.empty else 0,
            "n_t1_metadata_gap_rows": int(len(gap)),
            "t1_translational_readiness_html": str(hub_root / HTML_NAME),
        }
    )
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    update_manifest(
        output_root,
        "barneo_t1_translational_readiness_pack",
        {
            "outputs": outputs,
            "n_candidates": int(len(readiness)),
            "warnings": ["T1 translational pack distinguishes pHLA assay readiness from patient/translational claims."],
        },
    )
    print(
        "[barneo-t1-translational] "
        f"candidates={len(readiness)} assay_ready={manifest['summary']['n_t1_phla_assay_design_ready']} "
        f"translational_ready={manifest['summary']['n_t1_translational_review_ready']}"
    )


if __name__ == "__main__":
    main()
