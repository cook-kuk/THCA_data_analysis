# CLEAN-NeoBench Win/Loss Report KR

## 한 줄 결론

이제 “누가 이기고 지는지”가 split별로 고정됐다. 단, 이김은 benchmark slice 기준이고 SOTA/clinical claim이 아니다.

## Structure_LR 대비 요약

| method_name                                       | method_role                |   n_wins_vs_anchor |   n_losses_vs_anchor |   median_delta_AUPRC_vs_anchor |   source_heldout_median_delta_AUPRC |   hla_heldout_median_delta_AUPRC |
|:--------------------------------------------------|:---------------------------|-------------------:|---------------------:|-------------------------------:|------------------------------------:|---------------------------------:|
| W7B_stacked                                       | internal_candidate         |                 16 |                   13 |                     0.0946738  |                          0.172817   |                       0.0946738  |
| W7A_QK_only                                       | bounded_fallback           |                 21 |                    3 |                     0.088293   |                          0.088293   |                       0.0904311  |
| source_balanced_plus_pu_rf_train_prior_calibrated | internal_candidate         |                 12 |                   11 |                     0.0689085  |                         -0.325587   |                       0.0697866  |
| W7A_full                                          | internal_candidate         |                 17 |                    7 |                     0.0642506  |                          0.13852    |                       0.0642506  |
| pu_weighted_rf_train_prior_calibrated             | internal_candidate         |                 12 |                   11 |                     0.0570796  |                         -0.311408   |                       0.0629524  |
| source_balanced_rf_train_prior_calibrated         | internal_candidate         |                 12 |                   11 |                     0.056627   |                         -0.304438   |                       0.0798837  |
| MHCflurry                                         | caveated_public_comparator |                 15 |                   12 |                     0.0146306  |                          0.0362674  |                      -0.00700622 |
| Stack_mean_E1                                     | internal_candidate         |                 14 |                    9 |                     0.00137235 |                          0.0801305  |                       0          |
| BigMHC_IM                                         | caveated_public_comparator |                 14 |                   14 |                     0          |                          0.107489   |                       0.0335532  |
| Stack_LR_inmaster_E3a                             | internal_candidate         |                 14 |                   12 |                     0          |                          0.138726   |                       0          |
| VQC                                               | bounded_fallback           |                 10 |                   14 |                    -0.00720437 |                          0.0302225  |                      -0.101626   |
| Stack_median_E2                                   | internal_candidate         |                 10 |                   15 |                    -0.0159156  |                          0.00745875 |                       0.00650619 |
| kNN                                               | internal_candidate         |                 13 |                   16 |                    -0.0167897  |                          0.172817   |                      -0.0549693  |
| DeepImmuno                                        | caveated_public_comparator |                 11 |                   15 |                    -0.0219207  |                          0.0766955  |                      -0.151984   |
| MultiTask_A_only                                  | internal_candidate         |                 12 |                   15 |                    -0.023461   |                          0.0962295  |                      -0.134162   |
| rule_gated_C_QK_structure                         | bounded_fallback           |                  7 |                   12 |                    -0.0379278  |                         -0.089575   |                      -0.0408529  |
| MultiTask_AC                                      | internal_candidate         |                 14 |                   15 |                    -0.048964   |                          0.107341   |                      -0.122081   |
| prespecified_equal_weight_C_QK_no_anchor          | bounded_fallback           |                  6 |                   12 |                    -0.0509     |                         -0.101556   |                      -0.0438519  |
| prespecified_lr_qk_quantum_only_w0.5              | bounded_fallback           |                  2 |                   14 |                    -0.052332   |                         -0.0271753  |                      -0.012797   |
| MultiTask_AB                                      | internal_candidate         |                 12 |                   17 |                    -0.0557468  |                          0.107341   |                      -0.15717    |
| MultiTask_AD                                      | internal_candidate         |                 10 |                   17 |                    -0.0574946  |                          0.0874993  |                      -0.15717    |
| Wave8_TCR_SelfSim_full                            | internal_candidate         |                  5 |                   24 |                    -0.0579351  |                         -0.0539842  |                      -0.304945   |
| Wave8_TCR_SelfSim_no_exact                        | internal_candidate         |                 12 |                   18 |                    -0.0606859  |                         -0.048724   |                      -0.195884   |
| GP_quantum                                        | bounded_fallback           |                  6 |                   20 |                    -0.0619544  |                         -0.0619544  |                       0          |
| sourceheld_counterfactual_rf                      | internal_candidate         |                  8 |                   13 |                    -0.068192   |                         -0.30402    |                      -0.00393129 |
| Wave8_TCR_motif_only                              | internal_candidate         |                 10 |                   19 |                    -0.0726853  |                          0.0247566  |                      -0.244207   |
| TSCAPE_TITANiAN                                   | caveated_public_comparator |                 10 |                   19 |                    -0.0735363  |                         -0.0358582  |                      -0.0755385  |
| MHCnuggets_2                                      | caveated_public_comparator |                  9 |                   18 |                    -0.0740465  |                          0.00436087 |                      -0.147975   |
| MultiTask_full                                    | internal_candidate         |                 12 |                   15 |                    -0.0761678  |                          0.107341   |                      -0.152681   |
| TENT                                              | internal_candidate         |                 10 |                   19 |                    -0.0816529  |                          0.0962295  |                      -0.249464   |

## 실용 포지션

| method_name                                       | method_role                |   median_delta_AUPRC_vs_anchor | recommended_position                                                       |
|:--------------------------------------------------|:---------------------------|-------------------------------:|:---------------------------------------------------------------------------|
| W7B_stacked                                       | internal_candidate         |                     0.0946738  | beats anchor in matched split median; still needs source/HLA stress review |
| W7A_QK_only                                       | bounded_fallback           |                     0.088293   | bounded fallback/fusion component                                          |
| source_balanced_plus_pu_rf_train_prior_calibrated | internal_candidate         |                     0.0689085  | beats anchor in matched split median; still needs source/HLA stress review |
| W7A_full                                          | internal_candidate         |                     0.0642506  | beats anchor in matched split median; still needs source/HLA stress review |
| pu_weighted_rf_train_prior_calibrated             | internal_candidate         |                     0.0570796  | beats anchor in matched split median; still needs source/HLA stress review |
| source_balanced_rf_train_prior_calibrated         | internal_candidate         |                     0.056627   | beats anchor in matched split median; still needs source/HLA stress review |
| MHCflurry                                         | caveated_public_comparator |                     0.0146306  | field comparator only; public overlap audit required                       |
| Stack_mean_E1                                     | internal_candidate         |                     0.00137235 | beats anchor in matched split median; still needs source/HLA stress review |
| BigMHC_IM                                         | caveated_public_comparator |                     0          | field comparator only; public overlap audit required                       |
| Stack_LR_inmaster_E3a                             | internal_candidate         |                     0          | does not beat anchor median; retain as comparator or support branch        |
| VQC                                               | bounded_fallback           |                    -0.00720437 | bounded fallback/fusion component                                          |
| Stack_median_E2                                   | internal_candidate         |                    -0.0159156  | does not beat anchor median; retain as comparator or support branch        |
| kNN                                               | internal_candidate         |                    -0.0167897  | does not beat anchor median; retain as comparator or support branch        |
| DeepImmuno                                        | caveated_public_comparator |                    -0.0219207  | field comparator only; public overlap audit required                       |
| MultiTask_A_only                                  | internal_candidate         |                    -0.023461   | does not beat anchor median; retain as comparator or support branch        |
| rule_gated_C_QK_structure                         | bounded_fallback           |                    -0.0379278  | bounded fallback/fusion component                                          |
| MultiTask_AC                                      | internal_candidate         |                    -0.048964   | does not beat anchor median; retain as comparator or support branch        |
| prespecified_equal_weight_C_QK_no_anchor          | bounded_fallback           |                    -0.0509     | bounded fallback/fusion component                                          |
| prespecified_lr_qk_quantum_only_w0.5              | bounded_fallback           |                    -0.052332   | bounded fallback/fusion component                                          |
| MultiTask_AB                                      | internal_candidate         |                    -0.0557468  | does not beat anchor median; retain as comparator or support branch        |
| MultiTask_AD                                      | internal_candidate         |                    -0.0574946  | does not beat anchor median; retain as comparator or support branch        |
| Wave8_TCR_SelfSim_full                            | internal_candidate         |                    -0.0579351  | does not beat anchor median; retain as comparator or support branch        |
| Wave8_TCR_SelfSim_no_exact                        | internal_candidate         |                    -0.0606859  | does not beat anchor median; retain as comparator or support branch        |
| GP_quantum                                        | bounded_fallback           |                    -0.0619544  | bounded fallback/fusion component                                          |
| sourceheld_counterfactual_rf                      | internal_candidate         |                    -0.068192   | does not beat anchor median; retain as comparator or support branch        |

## 해석

- `Structure_LR`는 honest local anchor다.
- 내부 adaptive/PU/source-balanced 계열이 여러 split에서 anchor를 이길 수 있다.
- public pretrained methods는 이겨도 caveated comparator다.
- QK 계열은 bounded fallback/fusion component로만 둔다.
- contextual/failure-aware/claim-safe layer가 clean claim을 강하게 막는 현재 상태가 reviewer-safe하다.
