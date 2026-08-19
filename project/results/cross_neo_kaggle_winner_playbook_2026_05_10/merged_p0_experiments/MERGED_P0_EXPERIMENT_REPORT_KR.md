# CROSS-Neo Kaggle-merged P0 experiment report

Generated: 2026-05-10T12:00:53

## One-line result

Kaggle winner playbook에서 바로 옮길 수 있는 부분은 새 거대 모델이 아니라 **decision controller**였다. BAR-Neo 점수는 discovery와 claim-safe를 분리했고, TCR/MD/wetlab plate는 올리되 patient metadata와 overlap/leakage는 hard cap으로 잠갔다.

## Activated tricks

- OpenVaccine/Ribonanza style: heterogeneous ensemble은 유지하되 validation stress axis에서 무너지는 expert는 cap.
- BELKA/cheminformatics style: high-throughput 후보는 plate-first queue로 만들고, claim은 gate 통과 후보만 허용.
- Open Problems perturbation style: source/HLA/patient distribution shift를 별도 contract로 잠금.
- MoA/CAFA style: noisy label과 method-family disagreement를 별도 audit queue로 뽑음.
- PANDA pathology style: top-line score보다 reproducible operating point와 audit trail을 우선.

## Output counts

- candidate rows: 2715
- method rows: 100
- wetlab plate rows: 13
- TCR/structure disagreement rows: 649
- label-noise audit rows: 2715

## Candidate action counts

```json
{
  "audit_blocked_overlap_or_leakage": 2305,
  "watchlist_or_abstain": 390,
  "tcr_structure_disagreement_review": 17,
  "wetlab_priority_with_claim_boundary": 3
}
```

## Plate tier counts

```json
{
  "TIER_3_TCR_DIAGNOSTIC_ONLY": 11,
  "TIER_1_TCR_MD_CONTROL_READY": 2
}
```

## Generalization contract counts

```json
{
  "internal_support_only_until_stress_axes_pass": 77,
  "benchmark_claim_with_low_prevalence_caveat": 14,
  "diagnostic_or_ablation_only": 9
}
```

## Patient gate counts

```json
{
  "research_triage_only_missing_patient_context": 1750
}
```

## Top claim-safe/controller candidates

| candidate_id   | peptide         | hla_allele_4digit   |   bma_v2_claim_safe_score |   bma_v2_discovery_score | bma_v2_action                    | primary_claim_blocker   |
|:---------------|:----------------|:--------------------|--------------------------:|-------------------------:|:---------------------------------|:------------------------|
| CNV0_00169     | AALEDTLAETEAR   | HLA-B*37:01         |                     0.420 |                    0.701 | audit_blocked_overlap_or_leakage | patient context gate    |
| CNV0_01132     | GADGVGKSAL      | HLA-C*08:02         |                     0.420 |                    0.699 | audit_blocked_overlap_or_leakage | patient context gate    |
| CNV0_00804     | SRSYTSGPGSRISSS | HLA-B*27:09         |                     0.420 |                    0.699 | audit_blocked_overlap_or_leakage | patient context gate    |
| CNV0_00688     | YYYGIKDLATVFF   | HLA-A*24:02         |                     0.420 |                    0.698 | audit_blocked_overlap_or_leakage | patient context gate    |
| CNV0_00837     | EQFLDGDGWTSR    | HLA-A*34:02         |                     0.420 |                    0.696 | audit_blocked_overlap_or_leakage | patient context gate    |
| CNV0_00874     | FPSEYLSSHLEA    | HLA-B*54:01         |                     0.420 |                    0.692 | audit_blocked_overlap_or_leakage | patient context gate    |
| CNV0_00694     | AEIEGLKGQRASLE  | HLA-B*35:02         |                     0.420 |                    0.692 | audit_blocked_overlap_or_leakage | patient context gate    |
| CNV0_00494     | QEMASVEAAVSL    | HLA-A*02:01         |                     0.420 |                    0.690 | audit_blocked_overlap_or_leakage | patient context gate    |

## P0 wetlab plate v2

| row_id     | peptide    | hla_4digit   |   plate_v2_order_score | plate_v2_tier               | plate_v2_assay_bundle                                          |
|:-----------|:-----------|:-------------|-----------------------:|:----------------------------|:---------------------------------------------------------------|
| CNV0_01132 | GADGVGKSAL | HLA-C*08:02  |                  0.674 | TIER_1_TCR_MD_CONTROL_READY | mutant_vs_wt_decoy_pMHC_binding+TCR_tetramer+T_cell_activation |
| CNV0_01151 | HMTEVVRHC  | HLA-A*02:01  |                  0.603 | TIER_1_TCR_MD_CONTROL_READY | mutant_vs_wt_decoy_pMHC_binding+TCR_tetramer+T_cell_activation |
| CNV0_02448 | ILDKVLVHL  | HLA-A*02:01  |                  0.501 | TIER_3_TCR_DIAGNOSTIC_ONLY  | low-cost pMHC binding screen before TCR/MD expansion           |
| CNV0_02398 | KLILWRGLK  | HLA-A*03:01  |                  0.498 | TIER_3_TCR_DIAGNOSTIC_ONLY  | low-cost pMHC binding screen before TCR/MD expansion           |
| CNV0_01439 | SYLDSGIHF  | HLA-A*24:02  |                  0.441 | TIER_3_TCR_DIAGNOSTIC_ONLY  | low-cost pMHC binding screen before TCR/MD expansion           |
| CNV0_02478 | YVDFREYEYY | HLA-A*01:01  |                  0.435 | TIER_3_TCR_DIAGNOSTIC_ONLY  | low-cost pMHC binding screen before TCR/MD expansion           |
| CNV0_02479 | ILDTAGREEY | HLA-A*01:01  |                  0.411 | TIER_3_TCR_DIAGNOSTIC_ONLY  | low-cost pMHC binding screen before TCR/MD expansion           |
| CNV0_01118 | GADGVGKSA  | HLA-C*08:02  |                  0.360 | TIER_3_TCR_DIAGNOSTIC_ONLY  | low-cost pMHC binding screen before TCR/MD expansion           |

## Highest label-noise / source-conflict audits

| candidate_id   | peptide     | hla_allele_4digit   |   label_noise_priority_score | label_noise_bucket                      | label_noise_next_action                                                 |
|:---------------|:------------|:--------------------|-----------------------------:|:----------------------------------------|:------------------------------------------------------------------------|
| CNV0_00048     | TSGTPSCSW   | HLA-B*58:01         |                        0.741 | possible_false_negative_or_assay_miss   | manual assay/source review; consider rescue in wetlab if controls exist |
| CNV0_00711     | RREPPHLARNF | HLA-C*07:02         |                        0.739 | possible_false_negative_or_assay_miss   | manual assay/source review; consider rescue in wetlab if controls exist |
| CNV0_00682     | NYRREPPHL   | HLA-C*07:02         |                        0.736 | possible_false_negative_or_assay_miss   | manual assay/source review; consider rescue in wetlab if controls exist |
| CNV0_00674     | FKHSLSVSL   | HLA-C*07:02         |                        0.734 | possible_false_negative_or_assay_miss   | manual assay/source review; consider rescue in wetlab if controls exist |
| CNV0_00006     | VKEDPKWEF   | HLA-C*07:02         |                        0.732 | possible_false_negative_or_assay_miss   | manual assay/source review; consider rescue in wetlab if controls exist |
| CNV0_00027     | GNNDVKEDP   | HLA-C*07:02         |                        0.731 | possible_false_negative_or_assay_miss   | manual assay/source review; consider rescue in wetlab if controls exist |
| CNV0_01132     | GADGVGKSAL  | HLA-C*08:02         |                        0.683 | positive_label_overlap_or_leakage_audit | remove from clean claim set until overlap provenance is resolved        |
| CNV0_00058     | LLTQDLVQEK  | HLA-A*11:01         |                        0.678 | possible_false_negative_or_assay_miss   | manual assay/source review; consider rescue in wetlab if controls exist |

## Decision

1. `claim_safe_decision_controller`는 작동한다. 현재 데이터에서는 patient metadata와 overlap cap이 가장 강한 제한자라서, 높은 discovery 후보도 clinical/patient claim으로 자동 승격하지 않는다.
2. Wetlab은 `plate_v2_tier` 기준으로 진행한다. Tier 1은 TCR/MD/control-readiness가 같이 올라온 후보, Tier 2는 pMHC/DL 또는 benchmark claim 후보, Tier 3-4는 audit reserve다.
3. 다음 paper-blocking 실험은 label-noise/source-conflict 상위 후보와 TCR-structure disagreement 상위 후보의 수동 provenance audit이다. 새 모델 학습은 이 audit 후에만 의미 있다.

## Files

- `barneo_bma_v2_diversity_validity_scores.tsv`
- `barneo_bma_v2_method_weights.tsv`
- `cross_neo_wetlab_plate_v2.tsv`
- `cross_neo_generalization_contract_v2.tsv`
- `tcr_structure_disagreement_queue.tsv`
- `patient_context_gate_v2.tsv`
- `cross_neo_label_noise_audit_queue.tsv`
- `merged_p0_experiment_summary.json`
