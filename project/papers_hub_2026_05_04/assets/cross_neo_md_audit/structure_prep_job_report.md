# Counterfactual Structure Preparation Jobs

## Boundary

This file prepares future MD inputs. It does not launch new simulations and does not establish immunogenicity.

## Summary

- total prep jobs: 56
- OpenMM-ready existing/extracted template jobs: 7
- blocked or structure-generation-required jobs: 49

## Ready Jobs

| prep_id     | batch_id                                                            | peptide    | hla_4digit   | control_type                         | complex_kind   | sequence   | structure_route                                    | input_template_pdb                                                                                                                                                                           |
|:------------|:--------------------------------------------------------------------|:-----------|:-------------|:-------------------------------------|:---------------|:-----------|:---------------------------------------------------|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| MDPREP_0000 | P01_HMTEVVRHC_HLA_A_02_01_mutant_pmhc_3x10ns                        | HMTEVVRHC  | HLA-A*02:01  | mutant                               | pMHC           | HMTEVVRHC  | openmm_prepare_existing_candidate_template         | /home/seungho/personal/THCA_data_analysis/project/results/cross_neo_v2_sota_sprint_2026_05_09/tcr_extension/md_escalation/p0_structures/pilot_complexes/6VRN_HMTEVVRHC_HLA-A0201_chainP.pdb  |
| MDPREP_0001 | P01_HMTEVVRHC_HLA_A_02_01_mutant_tcr_pmhc_3x10ns                    | HMTEVVRHC  | HLA-A*02:01  | mutant                               | TCR-pMHC       | HMTEVVRHC  | openmm_prepare_existing_candidate_template         | /home/seungho/personal/THCA_data_analysis/project/results/cross_neo_v2_sota_sprint_2026_05_09/tcr_extension/md_escalation/p0_structures/pilot_complexes/6VRN_HMTEVVRHC_HLA-A0201_chainP.pdb  |
| MDPREP_0006 | P01_HMTEVVRHC_HLA_A_02_01_same_hla_positive_control_pMHC_3x10ns     | HMTEVVRHC  | HLA-A*02:01  | same_or_similar_hla_positive_control | pMHC           | ELAGIGILTV | openmm_prepare_existing_extracted_positive_control | /home/seungho/personal/THCA_data_analysis/project/results/cross_neo_md_audit_2026_05_10/counterfactual_design/positive_control_pdbs/extracted/5NQK_ELAGIGILTV_HLA-A0201_chainP.pdb           |
| MDPREP_0007 | P01_HMTEVVRHC_HLA_A_02_01_same_hla_positive_control_TCR_pMHC_3x10ns | HMTEVVRHC  | HLA-A*02:01  | same_or_similar_hla_positive_control | TCR-pMHC       | ELAGIGILTV | openmm_prepare_existing_extracted_positive_control | /home/seungho/personal/THCA_data_analysis/project/results/cross_neo_md_audit_2026_05_10/counterfactual_design/positive_control_pdbs/extracted/5NQK_ELAGIGILTV_HLA-A0201_chainP.pdb           |
| MDPREP_0008 | P02_GADGVGKSAL_HLA_C_08_02_mutant_pmhc_3x10ns                       | GADGVGKSAL | HLA-C*08:02  | mutant                               | pMHC           | GADGVGKSAL | openmm_prepare_existing_candidate_template         | /home/seungho/personal/THCA_data_analysis/project/results/cross_neo_v2_sota_sprint_2026_05_09/tcr_extension/md_escalation/p0_structures/pilot_complexes/6UON_GADGVGKSAL_HLA-C0802_chainF.pdb |
| MDPREP_0009 | P02_GADGVGKSAL_HLA_C_08_02_mutant_tcr_pmhc_3x10ns                   | GADGVGKSAL | HLA-C*08:02  | mutant                               | TCR-pMHC       | GADGVGKSAL | openmm_prepare_existing_candidate_template         | /home/seungho/personal/THCA_data_analysis/project/results/cross_neo_v2_sota_sprint_2026_05_09/tcr_extension/md_escalation/p0_structures/pilot_complexes/6UON_GADGVGKSAL_HLA-C0802_chainF.pdb |
| MDPREP_0014 | P02_GADGVGKSAL_HLA_C_08_02_same_hla_positive_control_pMHC_3x10ns    | GADGVGKSAL | HLA-C*08:02  | same_or_similar_hla_positive_control | pMHC           | GADGVGKSA  | openmm_prepare_existing_template                   | /home/seungho/personal/THCA_data_analysis/project/results/cross_neo_md_audit_2026_05_10/counterfactual_design/positive_control_pdbs/extracted/6ULR_GADGVGKSA_HLA-C0802_chainC.pdb            |

## Blocked Job Counts

| structure_route                                  | blocked_by                                                     |   n |
|:-------------------------------------------------|:---------------------------------------------------------------|----:|
| blocked_missing_template_or_structure            | requires_pmhc_structure_generation                             |  10 |
| blocked_missing_template_or_structure            | select_nonidentical_known_positive_same_or_similar_hla_control |  10 |
| blocked_missing_wt_sequence                      | missing_wt_sequence                                            |  12 |
| confirm_wt_then_thread_or_predict_structure      | confirm_wt_sequence_before_simulation                          |   2 |
| predict_decoy_pmhc_structure                     | prepare_decoy_structure_from_template                          |  10 |
| predict_or_fetch_positive_control_structure      | prepare_positive_control_structure                             |   1 |
| thread_decoy_on_candidate_template_then_pdbfixer | prepare_decoy_structure_from_template                          |   4 |
