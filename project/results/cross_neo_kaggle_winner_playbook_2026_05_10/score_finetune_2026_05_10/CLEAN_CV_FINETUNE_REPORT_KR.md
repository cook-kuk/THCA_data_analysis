# CROSS-Neo clean-CV scoring + finetune report

Generated: 2026-05-10T12:08:31

## Result

Clean low-leakage subset에서 source/HLA held-out OOF로 lightweight stacker를 튜닝했다. Best config는 `classical_stack_mix` / `C=0.3`이고, 선택 기준은 mean AUPRC, worst-axis AUPRC, top-10 precision, calibration penalty를 같이 본 것이다.

## Best fine-tuned CV score

| feature_set         |   c_value | protocol           |   n_features |   n_rows |   n_pos |   auprc |   auroc |   brier |   ece |   top10_precision |   top20_precision |
|:--------------------|----------:|:-------------------|-------------:|---------:|--------:|--------:|--------:|--------:|------:|------------------:|------------------:|
| classical_stack_mix |     0.300 | source_heldout_oof |           23 |      410 |      35 |   0.261 |   0.454 |   0.604 | 0.727 |             0.600 |             0.350 |
| classical_stack_mix |     0.300 | hla_heldout_oof    |           23 |      410 |      35 |   0.224 |   0.823 |   0.152 | 0.198 |             0.200 |             0.200 |

## Reference scores on same clean subset

These are not refit OOF models; they are context/reference scores on the same clean rows.

| feature_set                                 | protocol                         |   auprc |   auroc |   brier |   ece |   top10_precision |   top20_precision |
|:--------------------------------------------|:---------------------------------|--------:|--------:|--------:|------:|------------------:|------------------:|
| reference_stress_guarded_final_review_score | clean_subset_reference_not_refit |   0.658 |   0.937 |   0.078 | 0.169 |             0.900 |             0.700 |
| reference_bma_v2_discovery_score            | clean_subset_reference_not_refit |   0.620 |   0.892 |   0.158 | 0.328 |             0.800 |             0.650 |
| reference_stress_guarded_discovery_score    | clean_subset_reference_not_refit |   0.591 |   0.914 |   0.110 | 0.234 |             0.700 |             0.800 |
| reference_stress_guarded_clean_score        | clean_subset_reference_not_refit |   0.547 |   0.899 |   0.109 | 0.235 |             0.700 |             0.800 |
| reference_validity_dag_cap                  | clean_subset_reference_not_refit |   0.535 |   0.898 |   0.121 | 0.260 |             0.700 |             0.700 |
| reference_bma_v2_claim_safe_score           | clean_subset_reference_not_refit |   0.532 |   0.898 |   0.121 | 0.259 |             0.700 |             0.700 |

## Top rescored candidates

| candidate_id   | peptide     | hla_allele_4digit   |   label | source_name   | leakage_risk_level   |   finetuned_clean_stack_score |   bma_v2_discovery_score |   bma_v2_claim_safe_score | finetune_action                            | primary_claim_blocker   |
|:---------------|:------------|:--------------------|--------:|:--------------|:---------------------|------------------------------:|-------------------------:|--------------------------:|:-------------------------------------------|:------------------------|
| CNV0_01622     | ILDTAGHEEY  | HLA-A*01:01         |       1 | TESLA_mmc4    | high                 |                         0.979 |                    0.365 |                     0.272 | high_finetune_score_but_overlap_blocked    | antigen model gate      |
| CNV0_01859     | FLTKLVHLV   | HLA-A*02:01         |       0 | TESLA_mmc4    | high                 |                         0.979 |                    0.371 |                     0.277 | high_finetune_score_but_overlap_blocked    | antigen model gate      |
| CNV0_01732     | FMATYEINV   | HLA-A*02:01         |       0 | TESLA_mmc4    | high                 |                         0.977 |                    0.358 |                     0.265 | high_finetune_score_but_overlap_blocked    | antigen model gate      |
| CNV0_01960     | MARAIRVRTF  | HLA-B*08:01         |       0 | TESLA_mmc4    | high                 |                         0.976 |                    0.376 |                     0.281 | high_finetune_score_but_overlap_blocked    | antigen model gate      |
| CNV0_02079     | FLCEILRSMSI | HLA-A*02:01         |       0 | TESLA_mmc4    | high                 |                         0.973 |                    0.363 |                     0.269 | high_finetune_score_but_overlap_blocked    | antigen model gate      |
| CNV0_02452     | LLDGFLATV   | HLA-A*02:01         |       1 | ITSNdb_main   | low                  |                         0.972 |                    0.631 |                     0.420 | high_finetune_score_manual_candidate_audit | patient context gate    |
| CNV0_01642     | LMKKRDNL    | HLA-B*08:01         |       0 | TESLA_mmc4    | high                 |                         0.968 |                    0.374 |                     0.278 | high_finetune_score_but_overlap_blocked    | antigen model gate      |
| CNV0_01546     | PSAEVEMTFY  | HLA-A*01:01         |       0 | TESLA_mmc4    | high                 |                         0.967 |                    0.370 |                     0.275 | high_finetune_score_but_overlap_blocked    | antigen model gate      |
| CNV0_02398     | KLILWRGLK   | HLA-A*03:01         |       1 | ITSNdb_main   | medium               |                         0.962 |                    0.575 |                     0.420 | high_finetune_score_manual_candidate_audit | patient context gate    |
| CNV0_02077     | LTEQYNEKY   | HLA-A*01:01         |       0 | TESLA_mmc4    | high                 |                         0.960 |                    0.389 |                     0.291 | high_finetune_score_but_overlap_blocked    | antigen model gate      |
| CNV0_01220     | SSDSQEENY   | HLA-A*01:01         |       0 | NEPdb         | high                 |                         0.959 |                    0.268 |                     0.239 | high_finetune_score_but_overlap_blocked    | antigen model gate      |
| CNV0_02477     | VVVGAVGVG   | HLA-B*35:01         |       1 | ITSNdb_main   | low                  |                         0.959 |                    0.337 |                     0.322 | high_finetune_score_manual_candidate_audit | antigen model gate      |
| CNV0_01908     | VLFDRLSKLA  | HLA-A*02:01         |       0 | TESLA_mmc4    | high                 |                         0.959 |                    0.358 |                     0.266 | high_finetune_score_but_overlap_blocked    | antigen model gate      |
| CNV0_01994     | VLFDRLSKL   | HLA-A*02:01         |       0 | TESLA_mmc4    | high                 |                         0.958 |                    0.361 |                     0.268 | high_finetune_score_but_overlap_blocked    | antigen model gate      |
| CNV0_02007     | TMMCVSRNEL  | HLA-A*02:01         |       0 | TESLA_mmc4    | high                 |                         0.957 |                    0.364 |                     0.271 | high_finetune_score_but_overlap_blocked    | antigen model gate      |

## Top positive coefficients

| method_name                                   |   coefficient |   abs_coefficient |
|:----------------------------------------------|--------------:|------------------:|
| sourceheld_counterfactual_rf                  |        -1.040 |             1.040 |
| source_balanced_rf                            |         0.863 |             0.863 |
| source_balanced_plus_pu_rf                    |         0.777 |             0.777 |
| F_topk_diversification_hla_cluster            |        -0.736 |             0.736 |
| B_rank_normalization_hla_supertype_train_only |        -0.520 |             0.520 |
| C_source_balanced_training                    |        -0.460 |             0.460 |
| E_pu_style_conservative_ranker                |         0.447 |             0.447 |
| bagging_pu_rf                                 |         0.426 |             0.426 |
| decoy_focal_positive_pattern                  |         0.393 |             0.393 |
| positive_only_centroid                        |         0.359 |             0.359 |
| pu_logistic_weighted                          |        -0.349 |             0.349 |
| pu_weighted_rf                                |        -0.327 |             0.327 |
| pairwise_positive_over_unlabeled              |         0.327 |             0.327 |
| source_balanced_rf_train_prior_calibrated     |        -0.299 |             0.299 |
| hard_decoy_sequence_only_source_stress        |         0.280 |             0.280 |

## Underperforming components to cap/drop/recalibrate

| method_name                | finetune_component_action                 | underperformance_severity   |   clean_train_coverage |   source_min_AUPRC |   hla_min_AUPRC |   mean_ECE |
|:---------------------------|:------------------------------------------|:----------------------------|-----------------------:|-------------------:|----------------:|-----------:|
| Structure_LR               | drop_from_finetune_missing_clean_coverage | high                        |                  0.249 |              0.515 |           0.333 |      0.166 |
| W7B_stacked                | drop_from_finetune_missing_clean_coverage | high                        |                  0.256 |              0.464 |           0.298 |      0.219 |
| anchor_lr                  | drop_from_finetune_missing_clean_coverage | high                        |                  0.215 |              0.477 |           0.097 |      0.170 |
| Stack_mean_E1              | drop_from_finetune_missing_clean_coverage | high                        |                  0.256 |              0.516 |           0.256 |      0.232 |
| Stack_LR_inmaster_E3a      | drop_from_finetune_missing_clean_coverage | high                        |                  0.256 |              0.458 |           0.203 |      0.230 |
| Stack_median_E2            | drop_from_finetune_missing_clean_coverage | high                        |                  0.256 |              0.480 |           0.209 |      0.227 |
| Wave8_TCR_SelfSim_no_exact | drop_from_finetune_missing_clean_coverage | high                        |                  0.256 |              0.466 |           0.193 |      0.171 |
| RF_biophys                 | drop_from_finetune_missing_clean_coverage | high                        |                  0.256 |              0.417 |           0.174 |      0.248 |
| Wave8_TCR_motif_only       | drop_from_finetune_missing_clean_coverage | high                        |                  0.256 |              0.532 |           0.147 |      0.165 |
| Wave8_SelfSim_no_exact     | drop_from_finetune_missing_clean_coverage | high                        |                  0.256 |              0.420 |           0.150 |      0.180 |
| Wave8_TCR_SelfSim_full     | drop_from_finetune_missing_clean_coverage | high                        |                  0.256 |              0.484 |           0.143 |      0.182 |
| Wave8_TCR_only             | drop_from_finetune_missing_clean_coverage | high                        |                  0.256 |              0.556 |           0.158 |      0.203 |
| W7A_full                   | drop_from_finetune_missing_clean_coverage | high                        |                  0.256 |              0.508 |           0.265 |      0.253 |
| Wave8_SelfSim_full         | drop_from_finetune_missing_clean_coverage | high                        |                  0.256 |              0.410 |           0.133 |      0.193 |
| kNN                        | drop_from_finetune_missing_clean_coverage | high                        |                  0.256 |              0.498 |           0.272 |      0.217 |

## Decision

1. 점수는 `finetuned_clean_stack_score`로 냈다. 이 점수는 clean-CV tuned score이고, overlap/patient gate를 통과하지 못한 후보는 여전히 claim-safe가 아니다.
2. 너무 안 나오는 축은 새 모델을 바로 키우기보다 `underperforming_component_finetune_actions.tsv`의 cap/drop/recalibrate 지시대로 줄인다.
3. P0 wetlab 후보는 fine-tune score가 높아도 `overlap_clean_cap`과 `patient_context_gate`가 잠긴 상태라 실험 우선순위와 논문 claim을 분리한다.

## Files

- `clean_cv_finetune_metrics.tsv`
- `fine_tuned_candidate_scores.tsv`
- `fine_tuned_stack_coefficients.tsv`
- `underperforming_component_finetune_actions.tsv`
- `clean_cv_finetune_summary.json`
