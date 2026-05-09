# MD Positive Control Selection

## Boundary

Positive controls are sanity controls for MD runtime, contact analysis, and visualization. They are not cancer specificity or neoantigen validation evidence.

## Summary

- candidate control rows: 20
- selected manifest rows updated: 14

## Selected Rows

| batch_id                                                              | peptide         | hla_4digit   | complex_kind   | sequence   | sequence_status                                                       | readiness                           | blocked_by                                                     |
|:----------------------------------------------------------------------|:----------------|:-------------|:---------------|:-----------|:----------------------------------------------------------------------|:------------------------------------|:---------------------------------------------------------------|
| P01_HMTEVVRHC_HLA_A_02_01_same_hla_positive_control_pMHC_3x10ns       | HMTEVVRHC       | HLA-A*02:01  | pMHC           | ELAGIGILTV | selected_positive_control_exact_hla; tcr_row=TCRREG_0537478; pdb=5nqk | ready_existing_positive_control_pdb |                                                                |
| P01_HMTEVVRHC_HLA_A_02_01_same_hla_positive_control_TCR_pMHC_3x10ns   | HMTEVVRHC       | HLA-A*02:01  | TCR-pMHC       | ELAGIGILTV | selected_positive_control_exact_hla; tcr_row=TCRREG_0537478; pdb=5nqk | ready_existing_positive_control_pdb |                                                                |
| P02_GADGVGKSAL_HLA_C_08_02_same_hla_positive_control_pMHC_3x10ns      | GADGVGKSAL      | HLA-C*08:02  | pMHC           | GADGVGKSA  | selected_positive_control_exact_hla; tcr_row=TCRREG_0351706; pdb=6ULR | ready_existing_positive_control_pdb |                                                                |
| P02_GADGVGKSAL_HLA_C_08_02_same_hla_positive_control_TCR_pMHC_3x10ns  | GADGVGKSAL      | HLA-C*08:02  | TCR-pMHC       | FGDVGSTLF  | selected_positive_control_exact_hla; tcr_row=TCRREG_0053410; pdb=nan  | structure_generation_required       | prepare_positive_control_structure                             |
| P03_MAWSLGVLVALPFPL_HLA_B_40_01_same_hla_positive_control_pMHC_3x10ns | MAWSLGVLVALPFPL | HLA-B*40:01  | pMHC           | nan        | not_selected                                                          | blocked                             | select_nonidentical_known_positive_same_or_similar_hla_control |
| P04_MAWSLGVLVALPFPL_HLA_A_02_01_same_hla_positive_control_pMHC_3x10ns | MAWSLGVLVALPFPL | HLA-A*02:01  | pMHC           | nan        | not_selected                                                          | blocked                             | select_nonidentical_known_positive_same_or_similar_hla_control |
| P05_FSLVFLVYSVFKNNV_HLA_A_11_01_same_hla_positive_control_pMHC_3x10ns | FSLVFLVYSVFKNNV | HLA-A*11:01  | pMHC           | nan        | not_selected                                                          | blocked                             | select_nonidentical_known_positive_same_or_similar_hla_control |
| P06_MAWSLGVLVALPFPL_HLA_B_18_01_same_hla_positive_control_pMHC_3x10ns | MAWSLGVLVALPFPL | HLA-B*18:01  | pMHC           | nan        | not_selected                                                          | blocked                             | select_nonidentical_known_positive_same_or_similar_hla_control |
| P07_MAWSLGVLVALPFPL_HLA_A_25_01_same_hla_positive_control_pMHC_3x10ns | MAWSLGVLVALPFPL | HLA-A*25:01  | pMHC           | nan        | not_selected                                                          | blocked                             | select_nonidentical_known_positive_same_or_similar_hla_control |
| P08_MSSLAATTFHWKKCR_HLA_B_07_02_same_hla_positive_control_pMHC_3x10ns | MSSLAATTFHWKKCR | HLA-B*07:02  | pMHC           | nan        | not_selected                                                          | blocked                             | select_nonidentical_known_positive_same_or_similar_hla_control |
| P09_MSSLAATTFHWKKCR_HLA_A_68_01_same_hla_positive_control_pMHC_3x10ns | MSSLAATTFHWKKCR | HLA-A*68:01  | pMHC           | nan        | not_selected                                                          | blocked                             | select_nonidentical_known_positive_same_or_similar_hla_control |
| P10_MSSLAATTFHWKKCR_HLA_B_35_03_same_hla_positive_control_pMHC_3x10ns | MSSLAATTFHWKKCR | HLA-B*35:03  | pMHC           | nan        | not_selected                                                          | blocked                             | select_nonidentical_known_positive_same_or_similar_hla_control |
| P11_FSLVFLVYSVFKNNV_HLA_B_38_01_same_hla_positive_control_pMHC_3x10ns | FSLVFLVYSVFKNNV | HLA-B*38:01  | pMHC           | nan        | not_selected                                                          | blocked                             | select_nonidentical_known_positive_same_or_similar_hla_control |
| P12_FSLVFLVYSVFKNNV_HLA_B_52_01_same_hla_positive_control_pMHC_3x10ns | FSLVFLVYSVFKNNV | HLA-B*52:01  | pMHC           | nan        | not_selected                                                          | blocked                             | select_nonidentical_known_positive_same_or_similar_hla_control |
