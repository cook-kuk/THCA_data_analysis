# BAR-Neo T1 Impact Upgrade Plan KR

## 한 줄 결론

T1 3개를 “좋아 보이는 후보”에서 “바로 데이터/실험팀에 넘길 수 있는 실행 패킷”으로 바꿨다. 지금은 pHLA assay-design 후보이며, translational claim으로 올리려면 metadata intake와 go/no-go gate를 통과해야 한다.

## Metadata intake template

| candidate_id   | peptide   | hla_allele_4digit   | current_tier                        | wt_peptide   | gene   | mutation_id   | expression_tpm   | vaf   | patient_id   | validation_rule                                               |
|:---------------|:----------|:--------------------|:------------------------------------|:-------------|:-------|:--------------|:-----------------|:------|:-------------|:--------------------------------------------------------------|
| CNV0_02407     | KLMNIQQKL | HLA-A*02:01         | assay_design_ready_metadata_blocked |              |        |               |                  |       |              | fill source-backed values only; leave blank rather than infer |
| CNV0_02504     | LLVDLAEEL | HLA-A*02:01         | assay_design_ready_metadata_blocked |              |        |               |                  |       |              | fill source-backed values only; leave blank rather than infer |
| CNV0_02410     | MLGEQLFPL | HLA-A*02:01         | assay_design_ready_metadata_blocked |              |        |               |                  |       |              | fill source-backed values only; leave blank rather than infer |

## Assay design matrix

| candidate_id   |   assay_priority | assay_id                  | assay_or_review                                                       | required_input                                                    | current_status                |
|:---------------|-----------------:|:--------------------------|:----------------------------------------------------------------------|:------------------------------------------------------------------|:------------------------------|
| CNV0_02407     |                1 | peptide_synthesis_qc      | mutant peptide synthesis and QC                                       | peptide identity, purity, solubility                              | ready_to_plan                 |
| CNV0_02407     |                2 | hla_binding_stability     | HLA-A*02:01 binding/stability assay or source-backed binding evidence | binding/stability evidence for mutant peptide                     | ready_to_plan                 |
| CNV0_02407     |                3 | wt_specificity_comparator | WT peptide comparator review                                          | WT peptide sequence and mutant-vs-WT delta                        | blocked_until_metadata_linked |
| CNV0_02407     |                4 | tumor_antigen_evidence    | tumor antigen expression/clonality check                              | gene, mutation, expression, mutant expression, VAF/clonality      | blocked_until_metadata_linked |
| CNV0_02407     |                5 | presentation_intactness   | presentation intactness check                                         | HLA-LOH, B2M, antigen processing status if available              | blocked_until_metadata_linked |
| CNV0_02407     |                6 | immunogenicity_screen     | research immunogenicity screen                                        | T-cell activation/multimer/functional readout in research setting | blocked_until_metadata_linked |
| CNV0_02504     |                1 | peptide_synthesis_qc      | mutant peptide synthesis and QC                                       | peptide identity, purity, solubility                              | ready_to_plan                 |
| CNV0_02504     |                2 | hla_binding_stability     | HLA-A*02:01 binding/stability assay or source-backed binding evidence | binding/stability evidence for mutant peptide                     | ready_to_plan                 |
| CNV0_02504     |                3 | wt_specificity_comparator | WT peptide comparator review                                          | WT peptide sequence and mutant-vs-WT delta                        | blocked_until_metadata_linked |
| CNV0_02504     |                4 | tumor_antigen_evidence    | tumor antigen expression/clonality check                              | gene, mutation, expression, mutant expression, VAF/clonality      | blocked_until_metadata_linked |
| CNV0_02504     |                5 | presentation_intactness   | presentation intactness check                                         | HLA-LOH, B2M, antigen processing status if available              | blocked_until_metadata_linked |
| CNV0_02504     |                6 | immunogenicity_screen     | research immunogenicity screen                                        | T-cell activation/multimer/functional readout in research setting | blocked_until_metadata_linked |
| CNV0_02410     |                1 | peptide_synthesis_qc      | mutant peptide synthesis and QC                                       | peptide identity, purity, solubility                              | ready_to_plan                 |
| CNV0_02410     |                2 | hla_binding_stability     | HLA-A*02:01 binding/stability assay or source-backed binding evidence | binding/stability evidence for mutant peptide                     | ready_to_plan                 |
| CNV0_02410     |                3 | wt_specificity_comparator | WT peptide comparator review                                          | WT peptide sequence and mutant-vs-WT delta                        | blocked_until_metadata_linked |
| CNV0_02410     |                4 | tumor_antigen_evidence    | tumor antigen expression/clonality check                              | gene, mutation, expression, mutant expression, VAF/clonality      | blocked_until_metadata_linked |
| CNV0_02410     |                5 | presentation_intactness   | presentation intactness check                                         | HLA-LOH, B2M, antigen processing status if available              | blocked_until_metadata_linked |
| CNV0_02410     |                6 | immunogenicity_screen     | research immunogenicity screen                                        | T-cell activation/multimer/functional readout in research setting | blocked_until_metadata_linked |

## Go/no-go gate summary

| gate_id                  | current_status   |   n |
|:-------------------------|:-----------------|----:|
| G0_benchmark_cleanliness | pass             |   3 |
| G1_provenance_identity   | blocked          |   3 |
| G2_antigen_evidence      | blocked          |   3 |
| G3_presentation_safety   | blocked          |   3 |
| G4_patient_context       | blocked          |   3 |
| G5_research_assay_signal | not_started      |   3 |

## 제일 임팩트 큰 다음 액션

1. T1 3개에 WT peptide/gene/mutation/source window를 먼저 채운다.
2. expression/VAF/clonality를 붙이면 “benchmark-only”에서 “antigen evidence-supported”로 올라간다.
3. HLA-LOH/B2M/presentation status와 WT comparator를 붙이면 safety/presentation gate를 연다.
4. 이후 research assay readout이 생기면 experimental lead로 격상 가능하다.

## Claim boundary

clinical vaccine selection, new SOTA predictor, external validation proven, quantum advantage claim은 금지한다.
