#!/usr/bin/env python3
"""Build execution-ready upgrade pack for BAR-Neo T1 assay candidates."""

from __future__ import annotations

import argparse
import html
import json
from pathlib import Path

import pandas as pd

from common import dataframe_to_markdown, ensure_dir, update_manifest, write_tsv


HTML_NAME = "barneo_t1_impact_upgrade_plan_2026_05_10.html"
OUTPUTS = [
    "barneo_t1_metadata_intake_template.tsv",
    "barneo_t1_assay_design_matrix.tsv",
    "barneo_t1_go_nogo_criteria.tsv",
    "BAR_NEO_T1_IMPACT_UPGRADE_PLAN.md",
    "BAR_NEO_T1_IMPACT_UPGRADE_PLAN_KR.md",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo-root", default=".", help="Repository root")
    p.add_argument("--output-root", required=True, help="CLEAN-NeoBench output root")
    p.add_argument("--hub-root", default="project/papers_hub_2026_05_04", help="HTML hub root")
    return p.parse_args()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def build_metadata_intake(t1: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, r in t1.sort_values("stress_guarded_rank_global").iterrows():
        rows.append(
            {
                "candidate_id": r.get("candidate_id", ""),
                "peptide": r.get("peptide", ""),
                "hla_allele_4digit": r.get("hla_allele_4digit", ""),
                "current_tier": r.get("impact_readiness_tier", ""),
                "wt_peptide": "",
                "gene": "",
                "mutation_id": "",
                "protein_id": "",
                "source_protein_window": "",
                "expression_tpm": "",
                "mutant_expression": "",
                "vaf": "",
                "clonality": "",
                "patient_id": "",
                "sample_id": "",
                "cancer_type": "",
                "disease_context": "",
                "treatment_context": "",
                "hla_loh": "",
                "b2m_status": "",
                "antigen_processing_status": "",
                "immune_context_score": "",
                "ifng_score": "",
                "cytolytic_score": "",
                "evidence_file_path": "",
                "data_owner": "",
                "manual_review_note": "",
                "validation_rule": "fill source-backed values only; leave blank rather than infer",
                "ready_when": "gene/mutation/WT peptide/expression/clonality/patient/presentation fields are source-backed",
            }
        )
    return pd.DataFrame(rows)


def build_assay_design(t1: pd.DataFrame) -> pd.DataFrame:
    assays = [
        (
            "peptide_synthesis_qc",
            "mutant peptide synthesis and QC",
            "peptide identity, purity, solubility",
            "peptide feasible and QC-passing",
            "pHLA assay planning",
        ),
        (
            "hla_binding_stability",
            "HLA-A*02:01 binding/stability assay or source-backed binding evidence",
            "binding/stability evidence for mutant peptide",
            "binding/presentation support without relying on unaudited public predictors",
            "presentation gate",
        ),
        (
            "wt_specificity_comparator",
            "WT peptide comparator review",
            "WT peptide sequence and mutant-vs-WT delta",
            "mutant-supporting specificity; no strong WT safety red flag",
            "safety gate",
        ),
        (
            "tumor_antigen_evidence",
            "tumor antigen expression/clonality check",
            "gene, mutation, expression, mutant expression, VAF/clonality",
            "source-backed antigen evidence present",
            "translational gate",
        ),
        (
            "presentation_intactness",
            "presentation intactness check",
            "HLA-LOH, B2M, antigen processing status if available",
            "no known presentation hard fail",
            "patient gate",
        ),
        (
            "immunogenicity_screen",
            "research immunogenicity screen",
            "T-cell activation/multimer/functional readout in research setting",
            "signal above negative controls with WT comparator interpreted",
            "manual research triage",
        ),
    ]
    rows = []
    for _, r in t1.sort_values("stress_guarded_rank_global").iterrows():
        for rank, (assay_id, assay, required_input, pass_rule, gate) in enumerate(assays, 1):
            rows.append(
                {
                    "candidate_id": r.get("candidate_id", ""),
                    "peptide": r.get("peptide", ""),
                    "hla_allele_4digit": r.get("hla_allele_4digit", ""),
                    "assay_priority": rank,
                    "assay_id": assay_id,
                    "assay_or_review": assay,
                    "required_input": required_input,
                    "pass_rule": pass_rule,
                    "gate_unlocked": gate,
                    "current_status": "ready_to_plan" if assay_id in {"peptide_synthesis_qc", "hla_binding_stability"} else "blocked_until_metadata_linked",
                    "claim_boundary": "research assay planning only; not clinical selection",
                }
            )
    return pd.DataFrame(rows)


def build_go_nogo(t1: pd.DataFrame) -> pd.DataFrame:
    gates = [
        ("G0_benchmark_cleanliness", "BAR-Neo kill-audit pass, low leakage, no near/exact overlap", "pass", "keep as T1 manual-review candidate"),
        ("G1_provenance_identity", "gene, mutation_id, WT peptide, source protein/window source-backed", "blocked", "no translational claim until filled"),
        ("G2_antigen_evidence", "expression/mutant_expression/VAF/clonality source-backed", "blocked", "no patient prioritization until filled"),
        ("G3_presentation_safety", "HLA-LOH, B2M, antigen processing and WT comparator reviewed", "blocked", "no vaccine-first framing until reviewed"),
        ("G4_patient_context", "disease context, treatment timing, patient/sample ID linked", "blocked", "research triage only until linked"),
        ("G5_research_assay_signal", "binding/presentation/immunogenicity readout supports mutant peptide", "not_started", "upgrade to experimental lead only after signal"),
    ]
    rows = []
    for _, r in t1.sort_values("stress_guarded_rank_global").iterrows():
        for gate_id, criterion, current, action in gates:
            rows.append(
                {
                    "candidate_id": r.get("candidate_id", ""),
                    "peptide": r.get("peptide", ""),
                    "hla_allele_4digit": r.get("hla_allele_4digit", ""),
                    "gate_id": gate_id,
                    "criterion": criterion,
                    "current_status": current,
                    "go_action_if_pass": action,
                    "no_go_action_if_fail": "deprioritize or keep as benchmark-only example",
                    "claim_boundary": "stage-gated research triage",
                }
            )
    return pd.DataFrame(rows)


def write_reports(output_root: Path, hub_root: Path, intake: pd.DataFrame, assay: pd.DataFrame, gonogo: pd.DataFrame) -> None:
    intake_cols = [
        "candidate_id",
        "peptide",
        "hla_allele_4digit",
        "current_tier",
        "wt_peptide",
        "gene",
        "mutation_id",
        "expression_tpm",
        "vaf",
        "patient_id",
        "validation_rule",
    ]
    assay_cols = [
        "candidate_id",
        "assay_priority",
        "assay_id",
        "assay_or_review",
        "required_input",
        "current_status",
    ]
    gate_summary = (
        gonogo.groupby(["gate_id", "current_status"]).size().reset_index(name="n")
        if not gonogo.empty
        else pd.DataFrame()
    )
    md = f"""# BAR-Neo T1 Impact Upgrade Plan

## Purpose

This plan turns the three BAR-Neo T1 pHLA assay-design candidates into a concrete metadata intake, assay-design matrix, and go/no-go gate table.

## Metadata Intake Template

{dataframe_to_markdown(intake[intake_cols] if not intake.empty else intake, max_rows=20)}

## Assay Design Matrix

{dataframe_to_markdown(assay[assay_cols] if not assay.empty else assay, max_rows=40)}

## Gate Summary

{dataframe_to_markdown(gate_summary, max_rows=40)}

## Claim Boundary

This is an impact upgrade for research execution. It does not claim clinical vaccine selection, new SOTA prediction, external validation, or clean public-baseline status.
"""
    (output_root / "BAR_NEO_T1_IMPACT_UPGRADE_PLAN.md").write_text(md.strip() + "\n")

    kr = f"""# BAR-Neo T1 Impact Upgrade Plan KR

## 한 줄 결론

T1 3개를 “좋아 보이는 후보”에서 “바로 데이터/실험팀에 넘길 수 있는 실행 패킷”으로 바꿨다. 지금은 pHLA assay-design 후보이며, translational claim으로 올리려면 metadata intake와 go/no-go gate를 통과해야 한다.

## Metadata intake template

{dataframe_to_markdown(intake[intake_cols] if not intake.empty else intake, max_rows=20)}

## Assay design matrix

{dataframe_to_markdown(assay[assay_cols] if not assay.empty else assay, max_rows=40)}

## Go/no-go gate summary

{dataframe_to_markdown(gate_summary, max_rows=40)}

## 제일 임팩트 큰 다음 액션

1. T1 3개에 WT peptide/gene/mutation/source window를 먼저 채운다.
2. expression/VAF/clonality를 붙이면 “benchmark-only”에서 “antigen evidence-supported”로 올라간다.
3. HLA-LOH/B2M/presentation status와 WT comparator를 붙이면 safety/presentation gate를 연다.
4. 이후 research assay readout이 생기면 experimental lead로 격상 가능하다.

## Claim boundary

clinical vaccine selection, new SOTA predictor, external validation proven, quantum advantage claim은 금지한다.
"""
    (output_root / "BAR_NEO_T1_IMPACT_UPGRADE_PLAN_KR.md").write_text(kr.strip() + "\n")

    cards = []
    for cid, group in assay.groupby("candidate_id", sort=False):
        first = group.iloc[0]
        cards.append(
            f"""
            <section class="card">
              <h2>{html.escape(str(cid))} · {html.escape(str(first.get('peptide', '')))}</h2>
              <p class="meta">{html.escape(str(first.get('hla_allele_4digit', '')))} · pHLA assay-design candidate</p>
              <ol>{''.join(f"<li><b>{html.escape(str(row.assay_id))}</b>: {html.escape(str(row.current_status))}</li>" for _, row in group.iterrows())}</ol>
            </section>
            """
        )
    page = f"""<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>BAR-Neo T1 Impact Upgrade Plan</title>
<style>
body{{margin:0;background:#0d1117;color:#e6edf3;font-family:JetBrains Mono,ui-monospace,Menlo,monospace;line-height:1.55}}
header{{padding:44px 34px 24px;border-bottom:1px solid #2a3441}}h1{{font-family:Georgia,serif;font-size:46px;margin:0 0 10px}}.lead{{color:#c8d1dc;max-width:980px}}
main{{max-width:1180px;margin:auto;padding:24px;display:grid;gap:16px}}.card{{border:1px solid #2a3441;background:#101820;padding:18px}}h2{{font-family:Georgia,serif;margin:6px 0 4px;font-size:30px}}.meta{{color:#5eead4}}b{{color:#e3b341}}li{{margin:8px 0}}.path{{padding:0 24px 36px;color:#9aa7b4;font-size:12px}}
</style></head><body><header><h1>BAR-Neo T1 Impact Upgrade Plan</h1><p class="lead">Execution packet for the three T1 pHLA assay-design candidates. Metadata and assay gates explicitly block patient/translational claims until filled.</p></header><main>{''.join(cards)}</main><p class="path">Source: {html.escape(str(output_root / 'barneo_t1_assay_design_matrix.tsv'))}</p></body></html>"""
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

    t1 = read_tsv(output_root / "barneo_t1_translational_readiness.tsv")
    intake = build_metadata_intake(t1)
    assay = build_assay_design(t1)
    gonogo = build_go_nogo(t1)

    write_tsv(intake, output_root / "barneo_t1_metadata_intake_template.tsv")
    write_tsv(assay, output_root / "barneo_t1_assay_design_matrix.tsv")
    write_tsv(gonogo, output_root / "barneo_t1_go_nogo_criteria.tsv")
    write_reports(output_root, hub_root, intake, assay, gonogo)

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
            "n_t1_metadata_intake_rows": int(len(intake)),
            "n_t1_assay_design_rows": int(len(assay)),
            "n_t1_go_nogo_rows": int(len(gonogo)),
            "t1_impact_upgrade_html": str(hub_root / HTML_NAME),
        }
    )
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    update_manifest(
        output_root,
        "barneo_t1_impact_upgrade_pack",
        {
            "outputs": outputs,
            "n_intake_rows": int(len(intake)),
            "warnings": ["Impact upgrade pack is research execution planning, not clinical selection."],
        },
    )
    print(f"[barneo-t1-impact-upgrade] intake={len(intake)} assay={len(assay)} gates={len(gonogo)} html={hub_root / HTML_NAME}")


if __name__ == "__main__":
    main()
