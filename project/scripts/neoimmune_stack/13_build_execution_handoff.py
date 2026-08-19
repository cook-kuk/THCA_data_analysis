#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from neoimmune_common import ensure_run_dir, safe_read_table, write_md, write_tsv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()
    outdir = ensure_run_dir(args.outdir)
    canon = safe_read_table(outdir / "data" / "canonical_candidates.tsv.gz")
    clean = safe_read_table(outdir / "metrics" / "leaderboard_clean_track_strict_no_overlap.tsv")
    prod = safe_read_table(outdir / "metrics" / "leaderboard_production_track_strict_no_overlap.tsv")
    top20 = safe_read_table(outdir / "predictions" / "patient_top20_candidates.tsv")
    failure = safe_read_table(outdir / "predictions" / "failure_cases.tsv")

    req = pd.DataFrame(
        [
            ["patient_id", "required", "patient-level top-N endpoint", "de-identified stable ID", "blocks primary endpoint if missing"],
            ["sample_id", "required", "link tumor, normal, blood, assay", "stable sample ID", "blocks reproducibility if missing"],
            ["hla_typing", "required", "HLA restriction and external predictor runs", "4-digit class I minimum; class II optional", "blocks patient-specific ranking"],
            ["mutation_table", "required", "candidate generation", "gene, mutation_id, protein change, DNA change", "blocks candidate table"],
            ["mutant_peptide", "required", "model input", "8-11mer class I; longer/class II optional", "blocks all predictors"],
            ["wildtype_peptide", "strongly_required", "self-similarity and delta features", "same window as mutant where possible", "weakens immunogenicity branch"],
            ["expression_tpm", "required", "filter silent targets", "RNA-seq TPM or normalized expression", "high false-positive risk if missing"],
            ["vaf", "required", "clonality proxy", "tumor VAF", "subclonal target risk if missing"],
            ["clonality", "recommended", "patient-prioritized targeting", "clonal/subclonal estimate", "ranking uncertainty if missing"],
            ["hla_loh_status", "required", "presentation failure risk", "LOH yes/no/unknown", "false-positive presentation risk"],
            ["b2m_status", "required", "APM failure risk", "WT/mutated/loss/unknown", "false-positive presentation risk"],
            ["tap1_tap2_expr", "recommended", "APM activity", "TPM or z-score", "context uncertainty"],
            ["ms_immunopeptidomics", "optional_validation", "presentation validation", "present/absent if available", "needed for presentation claim"],
            ["t_cell_assay_result", "optional_validation", "immunogenicity validation", "positive/negative assay", "needed for immunogenicity claim"],
            ["clinical_context", "recommended", "interpret patient actionability", "tumor type, stage, treatment, biopsy time", "weak clinical handoff"],
        ],
        columns=["field", "priority", "why_needed", "acceptable_format", "risk_if_missing"],
    )
    write_tsv(req, outdir / "reports" / "hospital_data_request_sheet.tsv")

    review = pd.DataFrame(
        [
            ["candidate_id", "auto", "canonical candidate ID"],
            ["patient_id", "hospital", "de-identified patient ID"],
            ["rank_within_patient", "auto", "top-N position"],
            ["mutant_peptide", "hospital/auto", "peptide to synthesize/test"],
            ["hla_allele", "hospital", "restriction allele"],
            ["local_clean_score", "auto", "clean local model score"],
            ["production_stack_score", "auto", "practical prioritization score"],
            ["presentation_support", "auto/vendor", "NetMHC/MHCflurry/BigMHC/PRIME support summary"],
            ["tcr_self_similarity_support", "auto", "Wave8/TCR/self-similarity branch evidence"],
            ["structure_quantum_support", "auto", "Structure_LR and quantum branch evidence"],
            ["expression_vaf_gate", "hospital/auto", "expression and clonality suitability"],
            ["apm_hla_loh_gate", "hospital/auto", "presentation failure risk"],
            ["model_disagreement", "auto", "whether candidate is robust or controversial"],
            ["wetlab_decision", "human", "test / hold / reject / rescue-control"],
            ["claim_boundary", "auto", "what can and cannot be claimed"],
        ],
        columns=["field", "owner", "description"],
    )
    write_tsv(review, outdir / "reports" / "wetlab_candidate_review_form.tsv")

    matrix_rows = [
        ["Structure_LR", "clean", "structure proxy", "yes", "yes", "baseline clean method"],
        ["Wave8_TCR_SelfSim_full", "clean", "TCR-visible immunogenicity/self-similarity", "yes", "yes", "best biological story; strict no-overlap signal should be emphasized cautiously"],
        ["ESM2_Bayesian", "clean", "sequence representation", "yes", "yes", "strong existing artifact; needs strict/patient-heldout confirmation"],
        ["W7A_QK_only", "clean", "quantum kernel", "yes", "yes", "must survive source-heldout before big claim"],
        ["W7B_stacked", "clean", "local stack", "yes_if_local_only", "yes", "strong existing artifact; audit training scope"],
        ["MHCflurry", "external", "presentation", "no", "yes", "frozen comparator only"],
        ["BigMHC_IM", "external", "immunogenicity public predictor", "no", "yes", "strong comparator; public training contamination risk"],
        ["PRIME", "external", "presentation/immunogenicity", "no", "yes", "license-sensitive frozen comparator"],
        ["NetMHCpan", "external", "binding/presentation", "no", "yes", "binding is not immunogenicity"],
        ["Latest LLM", "analysis infra", "rationale and audit copilot", "no", "no_score_feature", "never predicts label; summarizes structured evidence only"],
    ]
    write_tsv(pd.DataFrame(matrix_rows, columns=["component", "source_type", "role", "clean_track_feature_allowed", "production_track_feature_allowed", "claim_boundary"]), outdir / "reports" / "algorithm_integration_matrix.tsv")

    attack = """# Reviewer Attack Defense Table KR

| Expected attack | Answer | Figure/table support | Claim boundary |
|---|---|---|---|
| 이건 binding predictor 조합 아닌가? | clean track과 production track을 분리했다. clean track은 public predictor score를 feature로 쓰지 않는다. | `model_registry.tsv`, `leaderboard_clean_track_strict_no_overlap.tsv` | production stack은 실용 레이어, clean novelty는 local-only track에서만 주장 |
| leakage가 심한 public benchmark 아닌가? | 맞다. 그래서 exact peptide, pHLA, study, patient, HLA, public-tool training risk를 별도 audit하고 strict no-existing-overlap view를 만들었다. | `leakage_audit.md`, `leakage_summary.tsv` | random split 숫자는 smoke test로만 사용 |
| patient-level ranking이라면서 patient ID가 없는데? | 현재 public scaffold는 patient ID가 부족하므로 patient-level endpoint를 claim하지 않는다. 병원 데이터 수집이 다음 gate다. | `patient_topN_report.md`, `hospital_data_request_sheet.tsv` | primary endpoint는 병원 per-patient table 확보 후 |
| clinical vaccine efficacy를 예측하나? | 아니다. 후보 우선순위화 시스템이다. efficacy claim은 금지한다. | `claim_ladder.tsv`, `claim_ladder.png` | no clinical efficacy claim |
| presentation을 증명했나? | 아니다. NetMHC/MHCflurry/BigMHC는 predictor다. MS 없으면 presentation confirmed라고 말하지 않는다. | `algorithm_integration_matrix.tsv` | MS 필요 |
| immunogenicity를 증명했나? | T-cell assay label이 있는 경우에만 benchmark label로 쓴다. 새로운 후보 immunogenicity는 wet-lab 전까지 hypothesis다. | `failure_case_audit.md` | T-cell assay 필요 |
| quantum kernel은 진짜 살아남나? | 현 단계에서는 promising branch다. source-heldout/no-overlap에서 독립적으로 살아남아야 한다. | `leaderboard_clean_track_strict_no_overlap.tsv` | overclaim 금지 |
| Wave8/TCR이 왜 중요한가? | vaccine bottleneck은 binding만이 아니라 TCR-visible immunogenicity다. strict no-overlap clean view에서 TCR/self-sim branch가 상위권이다. | strict clean leaderboard | 다른 candidate set comparator와 직접 비교 금지 |
| LLM 쓰면 환각/black box 아닌가? | LLM은 점수 feature가 아니다. structured evidence를 읽어 contradiction/rationale/report만 생성한다. | `latest_llm_copilot_blueprint.md` | no LLM-as-predictor claim |
"""
    write_md(attack, outdir / "reports" / "reviewer_attack_defense_table_kr.md")

    action = """# 7-Day Execution Plan KR

## Day 1
- 병원 데이터 필드 확정: `hospital_data_request_sheet.tsv` 그대로 전달.
- 내부 알고리즘 owner 확정: Structure/Wave8/ESM2/Quantum/BAR-Neo.

## Day 2
- patient-level candidate table 샘플 5명으로 스모크 테스트.
- HLA typing, expression, VAF, HLA LOH, B2M/APM 필수 컬럼 QC.

## Day 3
- external predictor frozen run: NetMHCpan, MHCflurry, BigMHC, PRIME.
- clean track에는 외부 predictor score가 절대 들어가지 않는지 audit.

## Day 4
- local algorithms 재실행: Structure_LR, Wave8_TCR_SelfSim_full, ESM2_Bayesian, quantum_kernel_no_anchor_gamma1.0.
- source/patient/no-overlap split 점검.

## Day 5
- patient top-20 후보 생성.
- high-score/low-expression, HLA-LOH, high-disagreement 후보 제거 또는 control로 이동.

## Day 6
- wet-lab 후보 회의: top candidates + rescue candidates + negative/disagreement controls.
- MS/T-cell assay 가능성 판단.

## Day 7
- paper/사업 pitch 업데이트.
- claim ladder 기준으로 초록/노랑/빨강 claim 분리.
"""
    write_md(action, outdir / "reports" / "next_7_day_execution_plan_kr.md")

    # Create a wet-lab queue scaffold from current top20, with hard warnings.
    q = top20.copy()
    for c in ["wetlab_priority", "test_type", "reason_to_test", "reason_to_hold", "claim_boundary"]:
        q[c] = ""
    if not q.empty:
        q["wetlab_priority"] = np.where(pd.to_numeric(q.get("model_disagreement_score"), errors="coerce").fillna(1) < 0.25, "review", "hold_or_disagreement_control")
        q["test_type"] = "format_scaffold_only_until_real_patient_data"
        q["reason_to_test"] = "high ranked candidate in current scaffold; requires real patient context before ordering assay"
        q["reason_to_hold"] = "public scaffold lacks real patient ID/expression/VAF/APM context"
        q["claim_boundary"] = "not a validated vaccine candidate"
    write_tsv(q, outdir / "predictions" / "wetlab_candidate_queue_scaffold.tsv")

    prompts_dir = outdir / "reports" / "llm_prompts"
    prompts_dir.mkdir(parents=True, exist_ok=True)
    (prompts_dir / "candidate_rationale_prompt.md").write_text(
        """You are an evidence auditor for cancer vaccine candidate ranking.
Input is a structured row with scores, leakage flags, expression, VAF, HLA LOH/APM, and labels.
Return:
1. one-sentence rationale,
2. top three supporting facts,
3. top three reasons to reject or hold,
4. exact claim boundary.
Do not invent missing wet-lab validation.
Do not call the candidate clinically validated.
""",
        encoding="utf-8",
    )
    (prompts_dir / "reviewer_defense_prompt.md").write_text(
        """Given the NeoImmune-Stack report tables, produce reviewer attacks and concise answers.
Rules:
- separate clean science from production stack,
- call out leakage,
- never equate HLA binding with immunogenicity,
- do not use LLM output as a model score.
""",
        encoding="utf-8",
    )

    fig_table = pd.DataFrame(
        [
            ["Figure 1", "System architecture", "neoimmune_stack_system_map.png", "How CLEAN-Neo++ separates clean local algorithms from production frozen public tools", "architecture only; no clinical validation"],
            ["Figure 2", "Strict clean leaderboard", "strict_clean_leaderboard.png + leaderboard_clean_track_strict_no_overlap.tsv", "Local algorithm behavior after excluding existing overlap flags", "small/public subset; needs patient/source-heldout confirmation"],
            ["Figure 3", "Strict production leaderboard", "strict_production_leaderboard.png + leaderboard_production_track_strict_no_overlap.tsv", "Practical stack and public comparator behavior", "production integration, not clean novelty"],
            ["Figure 4", "Claim ladder", "claim_ladder.png", "What is strong now vs what needs patient/wet-lab data", "explicitly forbids efficacy/presentation/immunogenicity overclaim"],
            ["Figure 5", "Patient data gap", "patient_data_gap.png", "Why hospital data is the next critical gate", "patient endpoint not claimable from current public scaffold"],
            ["Table 1", "Model registry", "registry/model_registry.tsv", "Public and local model roles and allowed tracks", "runnability/license status can change"],
            ["Table 2", "Canonical candidate schema", "data/canonical_candidates.tsv.gz", "Standardized candidate universe", "labels are heterogeneous public labels"],
            ["Table 3", "Leakage summary", "metrics/leakage_summary.tsv", "Leakage controls and flagged counts", "available flags only; does not prove no hidden contamination"],
            ["Table 4", "Algorithm integration matrix", "reports/algorithm_integration_matrix.tsv", "How each algorithm enters clean/production tracks", "does not imply every branch is validated"],
            ["Table 5", "Hospital data request", "reports/hospital_data_request_sheet.tsv", "Minimum fields to unlock patient-level endpoint", "required before top-N claim"],
        ],
        columns=["item", "title", "file", "role_in_story", "claim_boundary"],
    )
    write_tsv(fig_table, outdir / "reports" / "paper_figure_table_manifest_kr.tsv")

    results_scaffold = """# Manuscript Results Scaffold KR

## Result 1. Candidate universe and track separation
We standardized candidate records from existing CLEAN-Neobench/cross-neo/neoantigen hub artifacts into a canonical schema. The key design is a hard split between a clean-science track and a production-stack track.

Claim boundary: this is data/model integration, not clinical validation.

## Result 2. Public predictor scores are useful but cannot define clean novelty
MHCflurry, BigMHC, PRIME, and NetMHCpan-style artifacts can be collected as frozen comparators/features. They are excluded from clean-science training.

Claim boundary: binding/presentation predictors are not immunogenicity proof.

## Result 3. Local algorithms provide the clean scientific core
Structure_LR, Wave8/TCR-self-similarity, ESM2_Bayesian, W7A/W7B, and quantum-kernel families are registered as local branches. In strict no-existing-overlap views, TCR/self-similarity branches become a key biological story.

Claim boundary: strict public subset signal is promising but must be confirmed under true source-heldout/patient-heldout splits.

## Result 4. Leakage audit changes the interpretation of performance
The pipeline reports exact peptide, peptide-HLA, mutant-WT, study, patient, HLA, source-window, and public predictor training-contamination risk where available.

Claim boundary: random split and existing-artifact performance are smoke-test evidence only.

## Result 5. Production stack is strong as an operating layer
The production stack combines local branches with frozen external predictors and patient-context fields when available.

Claim boundary: production stacking can support prioritization, not a clean algorithm novelty claim.

## Result 6. Patient-level top-N is the right endpoint but needs hospital data
The code writes patient top-20 outputs and patient-level metrics, but the current public integrated table lacks real patient IDs.

Claim boundary: patient-level Recall@20/hit rate is not claimable until real patient-level candidate sets are obtained.

## Result 7. Wet-lab handoff is now explicit
The package defines what to test, what to hold, and what metadata blocks a candidate. This converts algorithm output into a collaborator-facing assay queue.

Claim boundary: candidates are hypotheses, not validated vaccine products.
"""
    write_md(results_scaffold, outdir / "reports" / "manuscript_results_scaffold_kr.md")

    business = """# NeoImmune-Stack Business Brief KR

## Product sentence
NeoImmune-Stack is an AI operating layer that turns patient-specific mutation candidates into leakage-audited, evidence-labeled, wet-lab-ready top-N vaccine candidate queues.

## Why this is commercially sharper than another binding predictor
- Most public tools score binding/presentation; the bottleneck is immunogenicity and patient-level prioritization.
- Hospitals need a decision layer that explains why a candidate should be synthesized/tested.
- The system can absorb existing public tools without being locked to one vendor predictor.
- Local algorithms remain the differentiating science layer.

## MVP
Input: patient mutation table, HLA typing, expression, VAF/clonality, HLA LOH/APM.
Output: top-20 candidate queue, evidence rationale, rejection reasons, wet-lab test plan, claim boundary.

## First paid/partner pilot
5-10 patients with matched sequencing and HLA typing. Run all frozen public tools plus local CLEAN-Neo++ branches. Select top 20 plus controls. Validate a subset by MS/T-cell assay if available.

## Do not sell as
- validated vaccine efficacy prediction
- standalone clinical diagnostic
- LLM-generated vaccine design
"""
    write_md(business, outdir / "reports" / "business_product_brief_kr.md")

    print(outdir / "reports" / "reviewer_attack_defense_table_kr.md")


if __name__ == "__main__":
    main()
