# MD Control Readiness Report

## Current Live Gate

- HMTEVVRHC / HLA-A*02:01 primary 10 ns status: 6600.0 ps, 300.40 K, complete=False.
- Remote watcher says running=True, replicate_watcher_running=True.
- Do not launch duplicate HMTEVVRHC primary jobs while the current primary is running.

## Ready Structure Jobs

- OpenMM-ready structure/template jobs: 7

| prep_id     | batch_id                                                            | peptide    | hla_4digit   | control_type                         | complex_kind   | sequence   | structure_route                                    |
|:------------|:--------------------------------------------------------------------|:-----------|:-------------|:-------------------------------------|:---------------|:-----------|:---------------------------------------------------|
| MDPREP_0000 | P01_HMTEVVRHC_HLA_A_02_01_mutant_pmhc_3x10ns                        | HMTEVVRHC  | HLA-A*02:01  | mutant                               | pMHC           | HMTEVVRHC  | openmm_prepare_existing_candidate_template         |
| MDPREP_0001 | P01_HMTEVVRHC_HLA_A_02_01_mutant_tcr_pmhc_3x10ns                    | HMTEVVRHC  | HLA-A*02:01  | mutant                               | TCR-pMHC       | HMTEVVRHC  | openmm_prepare_existing_candidate_template         |
| MDPREP_0006 | P01_HMTEVVRHC_HLA_A_02_01_same_hla_positive_control_pMHC_3x10ns     | HMTEVVRHC  | HLA-A*02:01  | same_or_similar_hla_positive_control | pMHC           | ELAGIGILTV | openmm_prepare_existing_extracted_positive_control |
| MDPREP_0007 | P01_HMTEVVRHC_HLA_A_02_01_same_hla_positive_control_TCR_pMHC_3x10ns | HMTEVVRHC  | HLA-A*02:01  | same_or_similar_hla_positive_control | TCR-pMHC       | ELAGIGILTV | openmm_prepare_existing_extracted_positive_control |
| MDPREP_0008 | P02_GADGVGKSAL_HLA_C_08_02_mutant_pmhc_3x10ns                       | GADGVGKSAL | HLA-C*08:02  | mutant                               | pMHC           | GADGVGKSAL | openmm_prepare_existing_candidate_template         |
| MDPREP_0009 | P02_GADGVGKSAL_HLA_C_08_02_mutant_tcr_pmhc_3x10ns                   | GADGVGKSAL | HLA-C*08:02  | mutant                               | TCR-pMHC       | GADGVGKSAL | openmm_prepare_existing_candidate_template         |
| MDPREP_0014 | P02_GADGVGKSAL_HLA_C_08_02_same_hla_positive_control_pMHC_3x10ns    | GADGVGKSAL | HLA-C*08:02  | same_or_similar_hla_positive_control | pMHC           | GADGVGKSA  | openmm_prepare_existing_template                   |

## Positive Control Extraction

| target_peptide   | control_peptide   | control_hla_4digit   | pdb_id   | status                          | selected_chains   |   contacts_tcr_peptide_initial |
|:-----------------|:------------------|:---------------------|:---------|:--------------------------------|:------------------|-------------------------------:|
| HMTEVVRHC        | ELAGIGILTV        | HLA-A*02:01          | 5NQK     | ready_TCR_pMHC_positive_control | A|B|H|L|P         |                             40 |
| HMTEVVRHC        | AAGIGILTV         | HLA-A*02:01          | 3QDJ     | ready_TCR_pMHC_positive_control | A|B|C|D|E         |                             31 |
| HMTEVVRHC        | ALGIGILTV         | HLA-A*02:01          | 4EUP     | ready_TCR_pMHC_positive_control | A|B|C|I|J         |                             62 |
| GADGVGKSAL       | GADGVGKSA         | HLA-C*08:02          | 6ULR     | incomplete_positive_control     | A|B|C             |                              0 |

## Blocked Work

- Blocked or structure-generation-required jobs: 49

| structure_route                                  | blocked_by                                                     |   n |
|:-------------------------------------------------|:---------------------------------------------------------------|----:|
| blocked_missing_template_or_structure            | requires_pmhc_structure_generation                             |  10 |
| blocked_missing_template_or_structure            | select_nonidentical_known_positive_same_or_similar_hla_control |  10 |
| blocked_missing_wt_sequence                      | missing_wt_sequence                                            |  12 |
| confirm_wt_then_thread_or_predict_structure      | confirm_wt_sequence_before_simulation                          |   2 |
| predict_decoy_pmhc_structure                     | prepare_decoy_structure_from_template                          |  10 |
| predict_or_fetch_positive_control_structure      | prepare_positive_control_structure                             |   1 |
| thread_decoy_on_candidate_template_then_pdbfixer | prepare_decoy_structure_from_template                          |   4 |

## Candidate Control Sequences

| row_id     | peptide         | hla_4digit   | candidate_sequence   | tentative_wt_sequence   | wt_status                                            | anchor_preserved_decoy   | template_pdb_id   |
|:-----------|:----------------|:-------------|:---------------------|:------------------------|:-----------------------------------------------------|:-------------------------|:------------------|
| CNV0_01151 | HMTEVVRHC       | HLA-A*02:01  | HMTEVVRHC            | nan                     | missing_wt_sequence_blocked                          | HMTHRVVEC                | 6VRN              |
| CNV0_01132 | GADGVGKSAL      | HLA-C*08:02  | GADGVGKSAL           | GAGGVGKSAL              | inferred_kras_g12d_like_requires_manual_confirmation | GASGVKGADL               | 6UON              |
| CNV0_01388 | MAWSLGVLVALPFPL | HLA-B*40:01  | MAWSLGVLVALPFPL      | nan                     | missing_wt_sequence_blocked                          | VAWPLGAVLMSPLFL          | nan               |
| CNV0_01384 | MAWSLGVLVALPFPL | HLA-A*02:01  | MAWSLGVLVALPFPL      | nan                     | missing_wt_sequence_blocked                          | VAWPLGAVLMSPLFL          | nan               |
| CNV0_01146 | FSLVFLVYSVFKNNV | HLA-A*11:01  | FSLVFLVYSVFKNNV      | nan                     | missing_wt_sequence_blocked                          | NSFYVNKLFSVLFVV          | nan               |
| CNV0_01163 | MAWSLGVLVALPFPL | HLA-B*18:01  | MAWSLGVLVALPFPL      | nan                     | missing_wt_sequence_blocked                          | VAWPLGAVLMSPLFL          | nan               |
| CNV0_01164 | MAWSLGVLVALPFPL | HLA-A*25:01  | MAWSLGVLVALPFPL      | nan                     | missing_wt_sequence_blocked                          | VAWPLGAVLMSPLFL          | nan               |
| CNV0_01276 | MSSLAATTFHWKKCR | HLA-B*07:02  | MSSLAATTFHWKKCR      | nan                     | missing_wt_sequence_blocked                          | TSKTHWLCASAKMFR          | nan               |
| CNV0_00930 | MSSLAATTFHWKKCR | HLA-A*68:01  | MSSLAATTFHWKKCR      | nan                     | missing_wt_sequence_blocked                          | TSKTHWLCASAKMFR          | nan               |
| CNV0_00931 | MSSLAATTFHWKKCR | HLA-B*35:03  | MSSLAATTFHWKKCR      | nan                     | missing_wt_sequence_blocked                          | TSKTHWLCASAKMFR          | nan               |
| CNV0_01391 | FSLVFLVYSVFKNNV | HLA-B*38:01  | FSLVFLVYSVFKNNV      | nan                     | missing_wt_sequence_blocked                          | NSFYVNKLFSVLFVV          | nan               |
| CNV0_01144 | FSLVFLVYSVFKNNV | HLA-B*52:01  | FSLVFLVYSVFKNNV      | nan                     | missing_wt_sequence_blocked                          | NSFYVNKLFSVLFVV          | nan               |

## Decision

- Best completed MD-supported candidate remains `GADGVGKSAL / HLA-C*08:02`.
- Best pending high-value candidate remains `HMTEVVRHC / HLA-A*02:01`; wait for the active 10 ns DCD before trajectory-derived claims.
- The next useful compute batch is controls, not more duplicate candidate-only runs: WT confirmation, decoy structures, and positive controls.
- Mutant-specific recognition and WT cross-reactivity claims remain blocked until WT/decoy/control trajectories are completed and analyzed.
