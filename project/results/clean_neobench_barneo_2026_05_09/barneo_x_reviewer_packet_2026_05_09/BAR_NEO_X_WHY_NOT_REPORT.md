# BAR-Neo-X Why-Not Audit

## Position

This audit explains why a candidate is priority, secondary, blocked, or metadata-required. It is a reviewer-facing triage explanation layer, not a clinical decision layer.

## Headline Counts

- Hard claim-blocked by leakage or identity overlap: 2304
- Rescuable by patient/disease metadata completion: 410
- Priority review now: 1
- Would cross priority threshold if metadata were complete: 23

## Primary Blocker Summary

| primary_blocker                              |   n_candidates |   n_positive_label |   median_claim_safe_score |   n_would_be_priority_if_metadata_complete | top_rescue_lane                                                | example_required_next_data                                                                         |
|:---------------------------------------------|---------------:|-------------------:|--------------------------:|-------------------------------------------:|:---------------------------------------------------------------|:---------------------------------------------------------------------------------------------------|
| claim_blocked_by_leakage_or_identity_overlap |           2304 |               1143 |                 0.0102881 |                                          0 | clean_claim_not_rescuable_without_new_holdout_or_overlap_audit | row-level training-corpus audit, de-duplicated split, independent external cohort                  |
| patient_metadata_incomplete                  |            410 |                 35 |                 0.094695  |                                         22 | rescuable_by_patient_metadata_completion                       | cancer type, disease context, tumor stage, expression, VAF/clonality, HLA-LOH, B2M, immune context |
| none_priority_candidate                      |              1 |                  1 |                 0.522455  |                                          1 | priority_review_now                                            | complete patient metadata and independent validation before translational claim                    |

## Reason Summary

| reason                                  |   n_candidates |   n_positive_label |   median_claim_safe_score | top_rescue_lane                                                |
|:----------------------------------------|---------------:|-------------------:|--------------------------:|:---------------------------------------------------------------|
| missing_patient_disease_metadata        |           2715 |               1179 |                0.0863286  | clean_claim_not_rescuable_without_new_holdout_or_overlap_audit |
| model_abstention                        |           2707 |               1179 |                0.086226   | clean_claim_not_rescuable_without_new_holdout_or_overlap_audit |
| low_claim_safe_score                    |           2676 |               1152 |                0.0852815  | clean_claim_not_rescuable_without_new_holdout_or_overlap_audit |
| high_leakage_claim_blocker              |           2304 |               1143 |                0.0102881  | clean_claim_not_rescuable_without_new_holdout_or_overlap_audit |
| near_peptide_overlap_caution            |           2300 |               1122 |                0.0100759  | clean_claim_not_rescuable_without_new_holdout_or_overlap_audit |
| exact_peptide_hla_overlap_claim_blocker |           2274 |               1114 |                0.00762889 | clean_claim_not_rescuable_without_new_holdout_or_overlap_audit |
| study_overlap_caution                   |           2086 |               1039 |                0.0117117  | clean_claim_not_rescuable_without_new_holdout_or_overlap_audit |
| benchmark_negative_label                |           1536 |                  0 |                0          | clean_claim_not_rescuable_without_new_holdout_or_overlap_audit |
| low_prevalence_source                   |           1035 |                 48 |                0          | clean_claim_not_rescuable_without_new_holdout_or_overlap_audit |
| expert_disagreement                     |            493 |                106 |                0.0937277  | rescuable_by_patient_metadata_completion                       |
| posterior_uncertainty                   |            317 |                 72 |                0.0939802  | rescuable_by_patient_metadata_completion                       |
| public_training_overlap_caveat          |            213 |                103 |                0          | clean_claim_not_rescuable_without_new_holdout_or_overlap_audit |
| sparse_expert_support                   |            200 |                 87 |                0          | clean_claim_not_rescuable_without_new_holdout_or_overlap_audit |
| medium_leakage_caution                  |             23 |                  5 |                0.0423284  | rescuable_by_patient_metadata_completion                       |
| underrepresented_hla_or_sparse_allele   |              3 |                  1 |                0.423746   | rescuable_by_patient_metadata_completion                       |
| exact_peptide_overlap_claim_caveat      |              1 |                  1 |                0.249595   | rescuable_by_patient_metadata_completion                       |

## Rescue Lane Summary

| rescue_lane                                                    |   n_candidates |   n_positive_label |   median_claim_safe_score | example_required_next_data                                                                         | claim_boundary                                         |
|:---------------------------------------------------------------|---------------:|-------------------:|--------------------------:|:---------------------------------------------------------------------------------------------------|:-------------------------------------------------------|
| clean_claim_not_rescuable_without_new_holdout_or_overlap_audit |           2304 |               1143 |                 0.0102881 | row-level training-corpus audit, de-duplicated split, independent external cohort                  | do not use as clean top-k benchmark evidence           |
| rescuable_by_patient_metadata_completion                       |            410 |                 35 |                 0.094695  | cancer type, disease context, tumor stage, expression, VAF/clonality, HLA-LOH, B2M, immune context | can be review-prioritized after metadata cap is lifted |
| priority_review_now                                            |              1 |                  1 |                 0.522455  | complete patient metadata and independent validation before translational claim                    | reviewer-facing research triage only                   |

## Top Claim-Safe Rows With Why-Not Context

| candidate_id   | peptide    | hla_allele_4digit   |   label | leakage_risk_level   |   barneo_x_claim_safe_score | primary_blocker             | rescue_lane                              | required_next_data                                                                                 |
|:---------------|:-----------|:--------------------|--------:|:---------------------|----------------------------:|:----------------------------|:-----------------------------------------|:---------------------------------------------------------------------------------------------------|
| CNV0_02407     | KLMNIQQKL  | HLA-A*02:01         |       1 | low                  |                      0.5225 | none_priority_candidate     | priority_review_now                      | complete patient metadata and independent validation before translational claim                    |
| CNV0_02503     | RLSDFSEQL  | HLA-A*02:01         |       1 | low                  |                      0.5    | patient_metadata_incomplete | rescuable_by_patient_metadata_completion | cancer type, disease context, tumor stage, expression, VAF/clonality, HLA-LOH, B2M, immune context |
| CNV0_02410     | MLGEQLFPL  | HLA-A*02:01         |       1 | low                  |                      0.4972 | patient_metadata_incomplete | rescuable_by_patient_metadata_completion | cancer type, disease context, tumor stage, expression, VAF/clonality, HLA-LOH, B2M, immune context |
| CNV0_02509     | SLLRSLENV  | HLA-A*02:01         |       1 | low                  |                      0.4957 | patient_metadata_incomplete | rescuable_by_patient_metadata_completion | cancer type, disease context, tumor stage, expression, VAF/clonality, HLA-LOH, B2M, immune context |
| CNV0_02409     | IILVAVPHV  | HLA-A*02:01         |       1 | low                  |                      0.4946 | patient_metadata_incomplete | rescuable_by_patient_metadata_completion | cancer type, disease context, tumor stage, expression, VAF/clonality, HLA-LOH, B2M, immune context |
| CNV0_02504     | LLVDLAEEL  | HLA-A*02:01         |       1 | low                  |                      0.4907 | patient_metadata_incomplete | rescuable_by_patient_metadata_completion | cancer type, disease context, tumor stage, expression, VAF/clonality, HLA-LOH, B2M, immune context |
| CNV0_02420     | KLANPLPYT  | HLA-A*02:01         |       0 | low                  |                      0.479  | patient_metadata_incomplete | rescuable_by_patient_metadata_completion | cancer type, disease context, tumor stage, expression, VAF/clonality, HLA-LOH, B2M, immune context |
| CNV0_02450     | KELEGILLL  | HLA-B*44:03         |       1 | low                  |                      0.4684 | patient_metadata_incomplete | rescuable_by_patient_metadata_completion | cancer type, disease context, tumor stage, expression, VAF/clonality, HLA-LOH, B2M, immune context |
| CNV0_02424     | LADEAEVYL  | HLA-A*02:01         |       1 | low                  |                      0.4663 | patient_metadata_incomplete | rescuable_by_patient_metadata_completion | cancer type, disease context, tumor stage, expression, VAF/clonality, HLA-LOH, B2M, immune context |
| CNV0_02452     | LLDGFLATV  | HLA-A*02:01         |       1 | low                  |                      0.4579 | patient_metadata_incomplete | rescuable_by_patient_metadata_completion | cancer type, disease context, tumor stage, expression, VAF/clonality, HLA-LOH, B2M, immune context |
| CNV0_02408     | FLYNLLTRV  | HLA-A*02:01         |       1 | low                  |                      0.4517 | patient_metadata_incomplete | rescuable_by_patient_metadata_completion | cancer type, disease context, tumor stage, expression, VAF/clonality, HLA-LOH, B2M, immune context |
| CNV0_02478     | YVDFREYEYY | HLA-A*01:01         |       1 | low                  |                      0.4462 | patient_metadata_incomplete | rescuable_by_patient_metadata_completion | cancer type, disease context, tumor stage, expression, VAF/clonality, HLA-LOH, B2M, immune context |
| CNV0_02448     | ILDKVLVHL  | HLA-A*02:01         |       1 | low                  |                      0.446  | patient_metadata_incomplete | rescuable_by_patient_metadata_completion | cancer type, disease context, tumor stage, expression, VAF/clonality, HLA-LOH, B2M, immune context |
| CNV0_02508     | YILKYSVFL  | HLA-A*02:01         |       1 | low                  |                      0.4438 | patient_metadata_incomplete | rescuable_by_patient_metadata_completion | cancer type, disease context, tumor stage, expression, VAF/clonality, HLA-LOH, B2M, immune context |
| CNV0_02507     | WVLALFDEV  | HLA-A*02:01         |       1 | low                  |                      0.4426 | patient_metadata_incomplete | rescuable_by_patient_metadata_completion | cancer type, disease context, tumor stage, expression, VAF/clonality, HLA-LOH, B2M, immune context |
| CNV0_02397     | LPIQYEPVL  | HLA-B*35:03         |       0 | low                  |                      0.4324 | patient_metadata_incomplete | rescuable_by_patient_metadata_completion | cancer type, disease context, tumor stage, expression, VAF/clonality, HLA-LOH, B2M, immune context |
| CNV0_02400     | ETSKQVTRW  | HLA-A*25:01         |       0 | low                  |                      0.4237 | patient_metadata_incomplete | rescuable_by_patient_metadata_completion | cancer type, disease context, tumor stage, expression, VAF/clonality, HLA-LOH, B2M, immune context |
| CNV0_02425     | GIVEGLITT  | HLA-A*02:01         |       1 | low                  |                      0.4124 | patient_metadata_incomplete | rescuable_by_patient_metadata_completion | cancer type, disease context, tumor stage, expression, VAF/clonality, HLA-LOH, B2M, immune context |
| CNV0_02460     | GIVEGLITTV | HLA-A*02:01         |       1 | low                  |                      0.4078 | patient_metadata_incomplete | rescuable_by_patient_metadata_completion | cancer type, disease context, tumor stage, expression, VAF/clonality, HLA-LOH, B2M, immune context |
| CNV0_02480     | ILDTAGKEEY | HLA-A*01:01         |       1 | medium               |                      0.3905 | patient_metadata_incomplete | rescuable_by_patient_metadata_completion | cancer type, disease context, tumor stage, expression, VAF/clonality, HLA-LOH, B2M, immune context |
| CNV0_02709     | ALPVALPSL  | HLA-A*02:01         |       1 | low                  |                      0.3902 | patient_metadata_incomplete | rescuable_by_patient_metadata_completion | cancer type, disease context, tumor stage, expression, VAF/clonality, HLA-LOH, B2M, immune context |
| CNV0_02401     | YIDERFERY  | HLA-A*01:01         |       0 | low                  |                      0.3851 | patient_metadata_incomplete | rescuable_by_patient_metadata_completion | cancer type, disease context, tumor stage, expression, VAF/clonality, HLA-LOH, B2M, immune context |
| CNV0_02458     | IIGAGPAEV  | HLA-A*02:01         |       0 | low                  |                      0.3798 | patient_metadata_incomplete | rescuable_by_patient_metadata_completion | cancer type, disease context, tumor stage, expression, VAF/clonality, HLA-LOH, B2M, immune context |
| CNV0_02710     | SLLSGLLRA  | HLA-A*02:01         |       1 | low                  |                      0.3732 | patient_metadata_incomplete | rescuable_by_patient_metadata_completion | cancer type, disease context, tumor stage, expression, VAF/clonality, HLA-LOH, B2M, immune context |
| CNV0_02477     | VVVGAVGVG  | HLA-B*35:01         |       1 | low                  |                      0.3639 | patient_metadata_incomplete | rescuable_by_patient_metadata_completion | cancer type, disease context, tumor stage, expression, VAF/clonality, HLA-LOH, B2M, immune context |
| CNV0_02572     | TTYSPIGEK  | HLA-A*03:01         |       0 | low                  |                      0.3587 | patient_metadata_incomplete | rescuable_by_patient_metadata_completion | cancer type, disease context, tumor stage, expression, VAF/clonality, HLA-LOH, B2M, immune context |
| CNV0_02405     | KMIGNHLWV  | HLA-A*02:01         |       1 | low                  |                      0.3568 | patient_metadata_incomplete | rescuable_by_patient_metadata_completion | cancer type, disease context, tumor stage, expression, VAF/clonality, HLA-LOH, B2M, immune context |
| CNV0_02712     | DKESEEEVS  | HLA-C*12:03         |       1 | low                  |                      0.3553 | patient_metadata_incomplete | rescuable_by_patient_metadata_completion | cancer type, disease context, tumor stage, expression, VAF/clonality, HLA-LOH, B2M, immune context |
| CNV0_02713     | AVCPWTWLR  | HLA-A*11:01         |       1 | low                  |                      0.3477 | patient_metadata_incomplete | rescuable_by_patient_metadata_completion | cancer type, disease context, tumor stage, expression, VAF/clonality, HLA-LOH, B2M, immune context |
| CNV0_02417     | LLSIIFFPA  | HLA-A*02:01         |       0 | low                  |                      0.3402 | patient_metadata_incomplete | rescuable_by_patient_metadata_completion | cancer type, disease context, tumor stage, expression, VAF/clonality, HLA-LOH, B2M, immune context |

## Claim Boundary

High leakage, exact peptide-HLA overlap, or patient overlap rows are not clean benchmark claims without a new holdout or row-level overlap audit. Missing metadata rows can be rescued only by completing patient/disease/presentation context.
