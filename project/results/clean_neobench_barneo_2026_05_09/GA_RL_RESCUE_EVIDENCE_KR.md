# GA/RL Rescue Evidence KR

## 한 줄 결론

회복 큐 3개는 exact/near/public overlap 없이 남았지만, 아직 patient/translational claim에는 필요한 metadata가 비어 있다. 이 패키지는 그 missingness를 바로 보여준다.

## Claim boundary

Allowed:

- rescue evidence intake
- metadata completion planning
- review-ready uplift queue

Forbidden:

- clinical vaccine selection
- new SOTA predictor
- public clean comparator claim without audit

## Summary

| evidence_grade   |   n_candidates |   mean_missing |   mean_priority |   mean_source_prev |   mean_hla_prev |
|:-----------------|---------------:|---------------:|----------------:|-------------------:|----------------:|
| C                |              3 |             14 |        0.563031 |           0.451605 |        0.348348 |

## Evidence queue

| candidate_id   | peptide   | hla_allele_4digit   | source_name   |   label |   ga_rl_barneo_priority_score |   top_nn_similarity |   source_positive_prevalence |   hla_positive_prevalence |   present_metadata_count |   missing_metadata_count | present_metadata   | missing_metadata                                                                                                                                                                     | evidence_grade   | evidence_action                             | evidence_boundary                            | rescue_note                                                    |
|:---------------|:----------|:--------------------|:--------------|--------:|------------------------------:|--------------------:|-----------------------------:|--------------------------:|-------------------------:|-------------------------:|:-------------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:-----------------|:--------------------------------------------|:---------------------------------------------|:---------------------------------------------------------------|
| CNV0_02448     | ILDKVLVHL | HLA-A*02:01         | ITSNdb_main   |       1 |                      0.62     |            0.555556 |                    0.648241  |                  0.348348 |                        0 |                       14 | none               | wt_peptide;gene;mutation_id;source_protein_window;expression_tpm;vaf;clonality;patient_id;disease_context;hla_loh;b2m_status;antigen_processing_status;tumor_stage;treatment_context | C                | intake missing metadata before T1 promotion | research triage only; not clinical selection | High-GA watchlist candidate without exact/near/public overlap. |
| CNV0_02708     | ALDPLLLRI | HLA-A*02:01         | ITSNdb_Val    |       1 |                      0.605917 |            0.555556 |                    0.0583333 |                  0.348348 |                        0 |                       14 | none               | wt_peptide;gene;mutation_id;source_protein_window;expression_tpm;vaf;clonality;patient_id;disease_context;hla_loh;b2m_status;antigen_processing_status;tumor_stage;treatment_context | C                | intake missing metadata before T1 promotion | research triage only; not clinical selection | High-GA watchlist candidate without exact/near/public overlap. |
| CNV0_02405     | KMIGNHLWV | HLA-A*02:01         | ITSNdb_main   |       1 |                      0.463175 |            0.4      |                    0.648241  |                  0.348348 |                        0 |                       14 | none               | wt_peptide;gene;mutation_id;source_protein_window;expression_tpm;vaf;clonality;patient_id;disease_context;hla_loh;b2m_status;antigen_processing_status;tumor_stage;treatment_context | C                | intake missing metadata before T1 promotion | research triage only; not clinical selection | High-GA watchlist candidate without exact/near/public overlap. |

## Missingness note

The rescue watchlist remains promising because it has no exact/near/public overlap, but the missing patient-level and biology fields keep it in research triage, not T1 translational claim.

## Output files

- `/home/seungho/personal/THCA_data_analysis/project/results/clean_neobench_barneo_2026_05_09/ga_rl_rescue_evidence.tsv`
- `/home/seungho/personal/THCA_data_analysis/project/results/clean_neobench_barneo_2026_05_09/ga_rl_rescue_evidence_summary.tsv`
- `/home/seungho/personal/THCA_data_analysis/project/results/clean_neobench_barneo_2026_05_09/figures/fig_ga_rl_rescue_evidence_missingness.png`
