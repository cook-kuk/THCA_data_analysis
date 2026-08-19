# BAR-Neo T1 Impact Upgrade Plan

## Purpose

This plan turns the three BAR-Neo T1 pHLA assay-design candidates into a concrete metadata intake, assay-design matrix, and go/no-go gate table.

## Metadata Intake Template

| candidate_id   | peptide   | hla_allele_4digit   | current_tier                        | wt_peptide   | gene   | mutation_id   | expression_tpm   | vaf   | patient_id   | validation_rule                                               |
|:---------------|:----------|:--------------------|:------------------------------------|:-------------|:-------|:--------------|:-----------------|:------|:-------------|:--------------------------------------------------------------|
| CNV0_02407     | KLMNIQQKL | HLA-A*02:01         | assay_design_ready_metadata_blocked |              |        |               |                  |       |              | fill source-backed values only; leave blank rather than infer |
| CNV0_02504     | LLVDLAEEL | HLA-A*02:01         | assay_design_ready_metadata_blocked |              |        |               |                  |       |              | fill source-backed values only; leave blank rather than infer |
| CNV0_02410     | MLGEQLFPL | HLA-A*02:01         | assay_design_ready_metadata_blocked |              |        |               |                  |       |              | fill source-backed values only; leave blank rather than infer |

## Assay Design Matrix

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

## Gate Summary

| gate_id                  | current_status   |   n |
|:-------------------------|:-----------------|----:|
| G0_benchmark_cleanliness | pass             |   3 |
| G1_provenance_identity   | blocked          |   3 |
| G2_antigen_evidence      | blocked          |   3 |
| G3_presentation_safety   | blocked          |   3 |
| G4_patient_context       | blocked          |   3 |
| G5_research_assay_signal | not_started      |   3 |

## Claim Boundary

This is an impact upgrade for research execution. It does not claim clinical vaccine selection, new SOTA prediction, external validation, or clean public-baseline status.
