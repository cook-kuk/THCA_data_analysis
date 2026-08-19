#!/usr/bin/env python3
"""Build a reviewer-safe Nature-grade evolution pack for CLEAN-NeoBench/BAR-Neo."""

from __future__ import annotations

import argparse
import html
import json
from pathlib import Path

import pandas as pd

from common import dataframe_to_markdown, ensure_dir, update_manifest, write_tsv


HTML_NAME = "barneo_nature_grade_evolution_2026_05_10.html"
OUTPUTS = [
    "clean_neobench_nature_claim_ladder.tsv",
    "clean_neobench_nature_validation_gap_table.tsv",
    "clean_neobench_nature_figure_blueprint.tsv",
    "clean_neobench_nature_reviewer_attack_matrix.tsv",
    "clean_neobench_nature_30day_sprint.tsv",
    "BAR_NEO_NATURE_GRADE_EVOLUTION_PLAN.md",
    "BAR_NEO_NATURE_GRADE_EVOLUTION_PLAN_KR.md",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo-root", default=".", help="Repository root")
    p.add_argument("--output-root", required=True, help="CLEAN-NeoBench output root")
    p.add_argument("--hub-root", default="project/papers_hub_2026_05_04", help="HTML hub root")
    return p.parse_args()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def manifest_summary(output_root: Path) -> dict[str, object]:
    path = output_root / "run_manifest.json"
    if not path.exists():
        return {}
    try:
        manifest = json.loads(path.read_text())
    except Exception:
        return {}
    summary = manifest.get("summary", {})
    return summary if isinstance(summary, dict) else {}


def val(summary: dict[str, object], key: str, default: object = 0) -> object:
    return summary.get(key, default)


def build_claim_ladder(summary: dict[str, object]) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "claim_level": "L1_current_safe_core",
                "claim": "Leakage-aware AI neoantigen predictor benchmarking framework",
                "current_status": "supported",
                "evidence_now": f"{val(summary, 'n_candidates')} candidates; {val(summary, 'n_methods')} methods; {val(summary, 'n_metric_rows')} split metrics",
                "required_to_upgrade": "freeze input manifests and add public training-corpus row-level overlap files",
                "forbidden_overclaim": "new SOTA predictor",
            },
            {
                "claim_level": "L2_current_safe_algorithm",
                "claim": "Benchmark-adaptive reliability ranking with abstention",
                "current_status": "supported",
                "evidence_now": f"BAR-Neo/BMA/stress-guarded layers; reviewer kill-audit pass={val(summary, 'n_reviewer_kill_pass')}",
                "required_to_upgrade": "prospective lockbox or independently held-out source validation",
                "forbidden_overclaim": "external validation proven",
            },
            {
                "claim_level": "L3_current_execution_ready",
                "claim": "pHLA assay-design candidates for manual research review",
                "current_status": "supported_with_boundary",
                "evidence_now": f"T1 assay-ready={val(summary, 'n_t1_phla_assay_design_ready')}; translational-ready={val(summary, 'n_t1_translational_review_ready')}",
                "required_to_upgrade": "WT peptide, gene, mutation, expression, clonality, patient context, HLA-LOH/B2M",
                "forbidden_overclaim": "clinical vaccine selection",
            },
            {
                "claim_level": "L4_near_term_upgrade",
                "claim": "Source/HLA robust reliability framework",
                "current_status": "partially_supported",
                "evidence_now": "source/HLA stress boards and win/loss reports exist",
                "required_to_upgrade": "pre-register source-heldout/HLA-heldout success thresholds and rerun after public overlap audit",
                "forbidden_overclaim": "Class I and Class II unified predictor",
            },
            {
                "claim_level": "L5_future_translational",
                "claim": "Translational neoantigen research triage system",
                "current_status": "blocked",
                "evidence_now": "patient-gated demo exists but patient metadata are sparse",
                "required_to_upgrade": "real PAAD/THCA patient metadata and at least one source-backed antigen/presentation/immune context layer",
                "forbidden_overclaim": "clinical utility",
            },
            {
                "claim_level": "L6_future_experimental",
                "claim": "Experimental support for selected T1 candidates",
                "current_status": "not_started",
                "evidence_now": "assay design matrix created",
                "required_to_upgrade": "binding/presentation/immunogenicity readout with WT comparator and controls",
                "forbidden_overclaim": "validated vaccine target without assay and safety evidence",
            },
        ]
    )


def build_validation_gap(summary: dict[str, object], public_audit: pd.DataFrame, t1: pd.DataFrame) -> pd.DataFrame:
    public_unresolved = int(val(summary, "n_public_methods_overlap_unresolved", 0) or 0)
    t1_meta_ready = int(val(summary, "n_t1_translational_review_ready", 0) or 0)
    return pd.DataFrame(
        [
            {
                "gap_id": "GAP01_public_overlap_audit",
                "blocks_claim": "public tools as clean comparators",
                "current_status": f"{public_unresolved} public methods unresolved; clean public comparators={val(summary, 'n_public_clean_comparators_allowed')}",
                "minimum_fix": "drop public training corpora into project/data/public_training_corpora and rerun row-level audit",
                "success_criterion": "clean_comparator_allowed_after_row_audit=true or method remains caveated",
                "impact_if_fixed": "turns comparator section from caveated to auditable",
            },
            {
                "gap_id": "GAP02_t1_antigen_identity",
                "blocks_claim": "T1 antigen identity",
                "current_status": f"T1 translational-ready={t1_meta_ready}; current T1 rows={len(t1)}",
                "minimum_fix": "fill WT peptide, gene, mutation_id, protein/source window for each T1 row",
                "success_criterion": "all T1 rows have source-backed antigen identity fields",
                "impact_if_fixed": "upgrades T1 from pHLA assay-design candidate to antigen identity-supported lead",
            },
            {
                "gap_id": "GAP03_t1_antigen_expression",
                "blocks_claim": "tumor antigen evidence",
                "current_status": "expression/VAF/clonality absent in current T1 master rows",
                "minimum_fix": "link expression_tpm, mutant_expression, VAF and clonality or mark not available",
                "success_criterion": "source-backed expression/clonality evidence or explicit no-data caveat",
                "impact_if_fixed": "enables translational triage score rather than benchmark-only score",
            },
            {
                "gap_id": "GAP04_presentation_safety",
                "blocks_claim": "presentation and safety gate",
                "current_status": "HLA-LOH/B2M/processing/WT comparator absent",
                "minimum_fix": "link HLA-LOH, B2M, antigen processing status and WT comparator review",
                "success_criterion": "no known presentation hard fail; WT comparator caveat resolved",
                "impact_if_fixed": "enables stronger pHLA-to-patient triage framing",
            },
            {
                "gap_id": "GAP05_independent_lockbox",
                "blocks_claim": "external/source-heldout robustness",
                "current_status": "source/HLA stress exists but not prospective lockbox",
                "minimum_fix": "freeze candidate/method manifests and evaluate on held-out source or new audit drop",
                "success_criterion": "predefined AUPRC/top-k/calibration thresholds met without modifying method weights",
                "impact_if_fixed": "moves from framework demonstration toward stronger validation paper",
            },
        ]
    )


def build_figure_blueprint() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "figure": "Fig1",
                "title": "CLEAN-NeoBench benchmark contract",
                "panel_plan": "candidate schema; leakage flags; split contracts; public-tool caveat logic",
                "source_outputs": "clean_neobench_master.tsv; clean_neobench_overlap_flags.tsv; CLEAN_NEOBENCH_METHOD_CARD.md",
                "nature_grade_message": "the contribution is a leakage-aware benchmark contract, not another black-box predictor",
            },
            {
                "figure": "Fig2",
                "title": "Method zoo under strict split contracts",
                "panel_plan": "overall leaderboard; clean internal board; source/HLA heldout; calibration",
                "source_outputs": "clean_neobench_leaderboard.tsv; clean_neobench_split_metrics.tsv; clean_neobench_winloss_method_summary.tsv",
                "nature_grade_message": "performance is shown by split and failure mode, not pooled AUROC",
            },
            {
                "figure": "Fig3",
                "title": "BAR-Neo reliability and abstention",
                "panel_plan": "BMA weights; contextual weights; abstention funnel; failure-aware downweighting",
                "source_outputs": "barneo_bma_method_weights.tsv; barneo_contextual_bma_method_weights.tsv; barneo_failure_aware_candidate_scores.tsv",
                "nature_grade_message": "benchmark behavior becomes calibrated reliability and abstention",
            },
            {
                "figure": "Fig4",
                "title": "Reviewer kill-audit and top-pass evidence",
                "panel_plan": "12 high-impact leads; pass/manual/blocked disposition; T1 method support",
                "source_outputs": "barneo_high_impact_reviewer_kill_audit.tsv; barneo_top_pass_reviewer_evidence.tsv",
                "nature_grade_message": "the system kills attractive but leaky candidates instead of hiding them",
            },
            {
                "figure": "Fig5",
                "title": "T1 translational readiness and execution plan",
                "panel_plan": "assay-ready vs metadata-blocked; intake template; assay matrix; go/no-go gates",
                "source_outputs": "barneo_t1_translational_readiness.tsv; barneo_t1_assay_design_matrix.tsv; barneo_t1_go_nogo_criteria.tsv",
                "nature_grade_message": "T1 candidates are converted into executable research tasks with claim boundaries",
            },
            {
                "figure": "ExtendedData",
                "title": "Distribution error and public overlap audit",
                "panel_plan": "source prevalence shift; low-prevalence errors; public corpus audit intake",
                "source_outputs": "clean_neobench_method_distribution_vulnerability.tsv; clean_neobench_public_overlap_audit_intake.tsv",
                "nature_grade_message": "failure analysis and comparator caveats are first-class evidence",
            },
        ]
    )


def build_reviewer_attack_matrix() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "reviewer_attack": "Your public baselines are contaminated",
                "defense_now": "public pretrained methods are caveated; clean allowed count remains 0",
                "remaining_gap": "row-level public training corpora are not yet dropped",
                "next_artifact": "clean_neobench_public_training_row_overlap_summary.tsv",
            },
            {
                "reviewer_attack": "Your model is just learning source labels",
                "defense_now": "source-heldout/study-heldout metrics and source/HLA stress boards are explicit",
                "remaining_gap": "prospective lockbox or new source drop still needed",
                "next_artifact": "frozen source-heldout validation report",
            },
            {
                "reviewer_attack": "Top candidates are leakage artifacts",
                "defense_now": "reviewer kill-audit blocks 8/12 high-impact rows and keeps only 3 T1",
                "remaining_gap": "manual provenance confirmation for T1 source identity",
                "next_artifact": "T1 metadata intake template completion",
            },
            {
                "reviewer_attack": "High score does not mean usable antigen",
                "defense_now": "T1 translational readiness separates pHLA assay-design from patient/translational claims",
                "remaining_gap": "WT/gene/mutation/expression/clonality missing",
                "next_artifact": "T1 antigen identity and expression evidence table",
            },
            {
                "reviewer_attack": "This is clinical vaccine selection",
                "defense_now": "all reports explicitly state research triage and assay-design boundary",
                "remaining_gap": "avoid clinical language in manuscript and demos",
                "next_artifact": "claim boundary checklist",
            },
        ]
    )


def build_sprint() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "priority": 1,
                "workstream": "T1 antigen identity",
                "task": "Fill WT peptide, gene, mutation_id, source_protein_window for CNV0_02407/CNV0_02504/CNV0_02410",
                "expected_output": "barneo_t1_metadata_intake_template.tsv completed",
                "upgrade_unlocked": "T1 antigen identity-supported leads",
            },
            {
                "priority": 2,
                "workstream": "Public overlap audit",
                "task": "Collect/drop public training corpora for NetMHCpan/MHCflurry/BigMHC/PRIME/MixMHCpred/etc.",
                "expected_output": "public training row-overlap audit summary",
                "upgrade_unlocked": "auditable comparator framing",
            },
            {
                "priority": 3,
                "workstream": "T1 assay planning",
                "task": "Prepare peptide synthesis/HLA-A*02:01 binding-stability planning table",
                "expected_output": "assay design matrix with owner/status fields",
                "upgrade_unlocked": "research assay execution packet",
            },
            {
                "priority": 4,
                "workstream": "Lockbox validation",
                "task": "Freeze runner outputs and evaluate on new held-out source or locked patient/study split",
                "expected_output": "source/HLA lockbox validation report",
                "upgrade_unlocked": "stronger external robustness claim",
            },
            {
                "priority": 5,
                "workstream": "Manuscript architecture",
                "task": "Build figure panels from current TSVs without writing voice-protected prose",
                "expected_output": "figure blueprint and factual result blocks",
                "upgrade_unlocked": "Nature Methods/NBME-style resource paper scaffold",
            },
        ]
    )


def write_reports(output_root: Path, hub_root: Path, tables: dict[str, pd.DataFrame], summary: dict[str, object]) -> None:
    claim = tables["claim_ladder"]
    gap = tables["validation_gap"]
    figs = tables["figure_blueprint"]
    attack = tables["reviewer_attack"]
    sprint = tables["sprint"]

    md = f"""# BAR-Neo Nature-Grade Evolution Plan

## Position

CLEAN-NeoBench + BAR-Neo is currently a leakage-aware benchmarking and reliability-triage framework. The Nature-grade path is to strengthen validation, provenance, and execution gates without overclaiming SOTA or clinical utility.

## Claim Ladder

{dataframe_to_markdown(claim, max_rows=20)}

## Validation Gaps

{dataframe_to_markdown(gap, max_rows=20)}

## Figure Blueprint

{dataframe_to_markdown(figs, max_rows=20)}

## Reviewer Attack Matrix

{dataframe_to_markdown(attack, max_rows=20)}

## 30-Day Sprint

{dataframe_to_markdown(sprint, max_rows=20)}

## Claim Boundary

Allowed: leakage-aware benchmark contract, benchmark-adaptive reliability ranking, pHLA assay-design research triage. Forbidden: clinical vaccine selection, new SOTA predictor, external validation proven, quantum advantage, clean public baselines without row-level audit.
"""
    (output_root / "BAR_NEO_NATURE_GRADE_EVOLUTION_PLAN.md").write_text(md.strip() + "\n")

    kr = f"""# BAR-Neo Nature-Grade Evolution Plan KR

## 한 줄 결론

지금 Nature급으로 진화시키는 길은 더 큰 모델이 아니라 **claim ladder + validation gate + T1 execution packet**이다. 현재는 Nature Methods/Nature Biomedical Engineering 스타일의 benchmark/resource core까지 왔고, T1 3개는 pHLA assay-design 후보로 올라갔다. Translational/Nature Cancer급 claim은 metadata와 assay/lockbox가 붙어야 한다.

## 현재 핵심 숫자

| 항목 | 값 |
|---|---:|
| candidates | {summary.get('n_candidates', 'NA')} |
| methods | {summary.get('n_methods', 'NA')} |
| split metrics | {summary.get('n_metric_rows', 'NA')} |
| reviewer kill pass | {summary.get('n_reviewer_kill_pass', 'NA')} |
| T1 pHLA assay ready | {summary.get('n_t1_phla_assay_design_ready', 'NA')} |
| T1 translational ready | {summary.get('n_t1_translational_review_ready', 'NA')} |
| public clean comparators allowed | {summary.get('n_public_clean_comparators_allowed', 'NA')} |

## Claim ladder

{dataframe_to_markdown(claim, max_rows=20)}

## Validation gaps

{dataframe_to_markdown(gap, max_rows=20)}

## Figure blueprint

{dataframe_to_markdown(figs, max_rows=20)}

## Reviewer attack matrix

{dataframe_to_markdown(attack, max_rows=20)}

## 30-day sprint

{dataframe_to_markdown(sprint, max_rows=20)}

## 결론

지금은 `Nature-grade methods/resource core`다. 다음 업그레이드는 T1 3개에 WT/gene/mutation/expression/clonality/patient/presentation metadata를 붙이고, public overlap audit과 source/HLA lockbox를 완료하는 것이다.

## Claim boundary

clinical vaccine selection, new SOTA predictor, external validation proven, quantum advantage claim은 금지한다.
"""
    (output_root / "BAR_NEO_NATURE_GRADE_EVOLUTION_PLAN_KR.md").write_text(kr.strip() + "\n")

    sections = [
        ("Claim Ladder", claim, ["claim_level", "claim", "current_status", "required_to_upgrade"]),
        ("Validation Gaps", gap, ["gap_id", "blocks_claim", "current_status", "minimum_fix"]),
        ("Figure Blueprint", figs, ["figure", "title", "nature_grade_message"]),
        ("Reviewer Attack Matrix", attack, ["reviewer_attack", "defense_now", "remaining_gap"]),
        ("30-Day Sprint", sprint, ["priority", "workstream", "task", "upgrade_unlocked"]),
    ]
    body = []
    for title, df, cols in sections:
        use = [c for c in cols if c in df.columns]
        head = "".join(f"<th>{html.escape(c)}</th>" for c in use)
        rows = []
        for _, row in df.loc[:, use].iterrows():
            rows.append("<tr>" + "".join(f"<td>{html.escape(str(row[c]))}</td>" for c in use) + "</tr>")
        body.append(f"<section><h2>{html.escape(title)}</h2><table><thead><tr>{head}</tr></thead><tbody>{''.join(rows)}</tbody></table></section>")

    page = f"""<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>BAR-Neo Nature-Grade Evolution</title>
<style>
body{{margin:0;background:#0d1117;color:#e6edf3;font-family:JetBrains Mono,ui-monospace,Menlo,monospace;line-height:1.55}}
header{{padding:44px 34px 24px;border-bottom:1px solid #2a3441}}h1{{font-family:Georgia,serif;font-size:46px;margin:0 0 10px}}.lead{{color:#c8d1dc;max-width:980px}}
.stats{{display:grid;grid-template-columns:repeat(5,minmax(120px,1fr));gap:10px;margin-top:20px}}.stats div{{border:1px solid #2a3441;background:#101820;padding:12px}}.stats b{{display:block;color:#5eead4;font-size:24px}}.stats span{{color:#9aa7b4;font-size:12px}}
main{{max-width:1320px;margin:auto;padding:24px;display:grid;gap:18px}}section{{border:1px solid #2a3441;background:#101820;padding:18px;overflow:auto}}h2{{font-family:Georgia,serif;font-size:30px;margin:0 0 10px}}table{{width:100%;border-collapse:collapse;font-size:12px}}th,td{{border-bottom:1px solid #2a3441;padding:8px;text-align:left;vertical-align:top}}th{{color:#e3b341}}.path{{padding:0 24px 36px;color:#9aa7b4;font-size:12px}}
</style></head><body><header><h1>BAR-Neo Nature-Grade Evolution</h1><p class="lead">Claim ladder, validation gaps, figure blueprint, reviewer attack matrix, and 30-day sprint for evolving CLEAN-NeoBench/BAR-Neo without SOTA or clinical overclaiming.</p><div class="stats">
<div><b>{summary.get('n_reviewer_kill_pass', 'NA')}</b><span>Kill-audit pass</span></div>
<div><b>{summary.get('n_t1_phla_assay_design_ready', 'NA')}</b><span>T1 assay-ready</span></div>
<div><b>{summary.get('n_t1_translational_review_ready', 'NA')}</b><span>T1 translational-ready</span></div>
<div><b>{summary.get('n_methods', 'NA')}</b><span>Methods</span></div>
<div><b>{summary.get('n_public_clean_comparators_allowed', 'NA')}</b><span>Clean public comparators</span></div>
</div></header><main>{''.join(body)}</main><p class="path">Source: {html.escape(str(output_root / 'BAR_NEO_NATURE_GRADE_EVOLUTION_PLAN_KR.md'))}</p></body></html>"""
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

    summary = manifest_summary(output_root)
    public_audit = read_tsv(output_root / "clean_neobench_public_tool_overlap_audit.tsv")
    t1 = read_tsv(output_root / "barneo_t1_translational_readiness.tsv")
    tables = {
        "claim_ladder": build_claim_ladder(summary),
        "validation_gap": build_validation_gap(summary, public_audit, t1),
        "figure_blueprint": build_figure_blueprint(),
        "reviewer_attack": build_reviewer_attack_matrix(),
        "sprint": build_sprint(),
    }
    write_tsv(tables["claim_ladder"], output_root / "clean_neobench_nature_claim_ladder.tsv")
    write_tsv(tables["validation_gap"], output_root / "clean_neobench_nature_validation_gap_table.tsv")
    write_tsv(tables["figure_blueprint"], output_root / "clean_neobench_nature_figure_blueprint.tsv")
    write_tsv(tables["reviewer_attack"], output_root / "clean_neobench_nature_reviewer_attack_matrix.tsv")
    write_tsv(tables["sprint"], output_root / "clean_neobench_nature_30day_sprint.tsv")
    write_reports(output_root, hub_root, tables, summary)

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
            "n_nature_claim_ladder_rows": int(len(tables["claim_ladder"])),
            "n_nature_validation_gap_rows": int(len(tables["validation_gap"])),
            "n_nature_figure_blueprint_rows": int(len(tables["figure_blueprint"])),
            "n_nature_reviewer_attack_rows": int(len(tables["reviewer_attack"])),
            "n_nature_30day_sprint_rows": int(len(tables["sprint"])),
            "nature_grade_evolution_html": str(hub_root / HTML_NAME),
        }
    )
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    update_manifest(
        output_root,
        "barneo_nature_grade_evolution_pack",
        {
            "outputs": outputs,
            "n_claim_ladder_rows": int(len(tables["claim_ladder"])),
            "warnings": ["Nature-grade evolution plan is a claim architecture and validation roadmap, not a SOTA or clinical claim."],
        },
    )
    print(
        "[barneo-nature-evolution] "
        f"claims={len(tables['claim_ladder'])} gaps={len(tables['validation_gap'])} "
        f"figures={len(tables['figure_blueprint'])} sprint={len(tables['sprint'])} html={hub_root / HTML_NAME}"
    )


if __name__ == "__main__":
    main()
