# BAR-Neo T1 Translational Readiness Pack

## Purpose

This pack separates pHLA assay-design readiness from patient/translational readiness for the BAR-Neo T1 leads that passed the reviewer kill audit.

## T1 Readiness

| candidate_id   |   stress_guarded_rank_global |   stress_guarded_final_review_score | peptide   | hla_allele_4digit   |   benchmark_evidence_score | peptide_hla_assay_design_ready   |   translational_metadata_score | impact_readiness_tier               | recommended_next_action                                                                              |
|:---------------|-----------------------------:|------------------------------------:|:----------|:--------------------|---------------------------:|:---------------------------------|-------------------------------:|:------------------------------------|:-----------------------------------------------------------------------------------------------------|
| CNV0_02407     |                            1 |                            0.551702 | KLMNIQQKL | HLA-A*02:01         |                          1 | True                             |                              0 | assay_design_ready_metadata_blocked | design pHLA-level review/assay plan, but block patient/translational claim until metadata are linked |
| CNV0_02504     |                            2 |                            0.541321 | LLVDLAEEL | HLA-A*02:01         |                          1 | True                             |                              0 | assay_design_ready_metadata_blocked | design pHLA-level review/assay plan, but block patient/translational claim until metadata are linked |
| CNV0_02410     |                            3 |                            0.511978 | MLGEQLFPL | HLA-A*02:01         |                          1 | True                             |                              0 | assay_design_ready_metadata_blocked | design pHLA-level review/assay plan, but block patient/translational claim until metadata are linked |

## Assay Handoff

| candidate_id   | assay_handoff_tier          | peptide   | hla_allele_4digit   | anchor_proxy   | minimum_next_data_packet                                                                                         |
|:---------------|:----------------------------|:----------|:--------------------|:---------------|:-----------------------------------------------------------------------------------------------------------------|
| CNV0_02407     | pHLA_assay_design_candidate | KLMNIQQKL | HLA-A*02:01         | P2=L; C-term=L | WT peptide, gene, mutation_id, expression/mutant_expression, VAF/clonality, patient disease context, HLA-LOH/B2M |
| CNV0_02504     | pHLA_assay_design_candidate | LLVDLAEEL | HLA-A*02:01         | P2=L; C-term=L | WT peptide, gene, mutation_id, expression/mutant_expression, VAF/clonality, patient disease context, HLA-LOH/B2M |
| CNV0_02410     | pHLA_assay_design_candidate | MLGEQLFPL | HLA-A*02:01         | P2=L; C-term=L | WT peptide, gene, mutation_id, expression/mutant_expression, VAF/clonality, patient disease context, HLA-LOH/B2M |

## Metadata Gap Summary

| gate_group          | field                     |   present_fraction |
|:--------------------|:--------------------------|-------------------:|
| antigen_evidence    | clonality                 |                  0 |
| antigen_evidence    | expression_tpm            |                  0 |
| antigen_evidence    | mutant_expression         |                  0 |
| antigen_evidence    | vaf                       |                  0 |
| antigen_identity    | gene                      |                  0 |
| antigen_identity    | mutation_id               |                  0 |
| antigen_identity    | wt_peptide                |                  0 |
| immune_context      | cytolytic_score           |                  0 |
| immune_context      | ifng_score                |                  0 |
| immune_context      | immune_context_score      |                  0 |
| patient_gate        | cancer_type               |                  0 |
| patient_gate        | disease_context           |                  0 |
| patient_gate        | patient_id                |                  0 |
| patient_gate        | treatment_context         |                  0 |
| presentation_safety | antigen_processing_status |                  0 |
| presentation_safety | b2m_status                |                  0 |
| presentation_safety | hla_loh                   |                  0 |
| provenance          | source_protein_window     |                  0 |

## Claim Boundary

The T1 rows are pHLA assay-design/manual-review candidates. They are not clinical vaccine selections, not new SOTA validation, and not patient-level translational claims until WT/gene/mutation/expression/clonality/patient/safety metadata are linked.
