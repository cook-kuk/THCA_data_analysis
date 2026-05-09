# CLEAN-NeoBench / BAR-Neo Current Status KR

## 한 줄 결론

현재 위치는 `new SOTA predictor`가 아니라 **leakage-aware AI neoantigen benchmark + benchmark-adaptive reliability/abstention controller**다. 실용적으로는 agent가 여러 expert/ensemble을 불러오고, BAR-Neo-BMA가 benchmark 성능과 실패 패턴으로 가중치를 조절한 뒤, claim-safe layer가 위험한 후보를 abstain/manual review로 보낸다.

## 현재 run 숫자

| 항목 | 값 |
|---|---:|
| candidates | 2,715 |
| methods | 100 |
| BMA weighted methods | 100 |
| split metric rows | 7,144 |
| BAR-Neo-BMA candidates | 2,715 |
| BAR-Neo-BMA abstain | 2,708 |
| BAR-Neo-BMA non-abstain | 7 |
| contextual clean claims allowed | 0 |
| BAR-Neo-X priority review candidates | 2 |
| BAR-Neo-X top claim-safe score | 0.529 |
| challenge rows | 885 |
| visual dashboard figures | 12 |
| public clean comparators allowed | 0 |

## 우리 알고리즘 위치

- `Structure_LR`: honest local anchor.
- `W7B_stacked`, `W7A_full`, source-balanced/PU RF: 현재 clean internal 후보군.
- `BAR-Neo`: 후보 intrinsic + method rank + reliability + uncertainty + leakage/source/HLA context로 candidate reliability score와 abstention reason을 낸다.
- `BAR-Neo-BMA`: 실용 agent ensemble controller. expert utility, top-k, calibration, source collapse, role prior, public caveat penalty, QK bounded fallback penalty를 가중치에 반영한다.
- `BAR-Neo-X`: leakage/claim-safe reranking과 설명 layer. 지금 clean claim은 막고 priority review만 만든다.
- QK 계열: 성능 좋은 slice가 있어도 **bounded fallback/fusion component**다. headline 또는 quantum advantage claim이 아니다.

## Structure_LR 대비 승패 핵심표

| method_name                                       | method_role                | method_disposition              |   n_matched_anchor_splits |   n_wins_vs_anchor |   n_losses_vs_anchor |   n_ties_vs_anchor |   win_rate_vs_anchor |   median_delta_AUPRC_vs_anchor |   source_heldout_median_delta_AUPRC |   hla_heldout_median_delta_AUPRC |
|:--------------------------------------------------|:---------------------------|:--------------------------------|--------------------------:|-------------------:|---------------------:|-------------------:|---------------------:|-------------------------------:|------------------------------------:|---------------------------------:|
| W7B_stacked                                       | internal_candidate         | clean_internal_candidate        |                        29 |                 16 |                   13 |                  0 |            0.551724  |                     0.0946738  |                          0.172817   |                       0.0946738  |
| W7A_QK_only                                       | bounded_fallback           | bounded_fallback_not_headline   |                        29 |                 21 |                    3 |                  5 |            0.724138  |                     0.088293   |                          0.088293   |                       0.0904311  |
| source_balanced_plus_pu_rf_train_prior_calibrated | internal_candidate         | clean_internal_candidate        |                        23 |                 12 |                   11 |                  0 |            0.521739  |                     0.0689085  |                         -0.325587   |                       0.0697866  |
| W7A_full                                          | internal_candidate         | clean_internal_candidate        |                        29 |                 17 |                    7 |                  5 |            0.586207  |                     0.0642506  |                          0.13852    |                       0.0642506  |
| pu_weighted_rf_train_prior_calibrated             | internal_candidate         | clean_internal_candidate        |                        23 |                 12 |                   11 |                  0 |            0.521739  |                     0.0570796  |                         -0.311408   |                       0.0629524  |
| source_balanced_rf_train_prior_calibrated         | internal_candidate         | clean_internal_candidate        |                        23 |                 12 |                   11 |                  0 |            0.521739  |                     0.056627   |                         -0.304438   |                       0.0798837  |
| MHCflurry                                         | caveated_public_comparator | caveated_public_not_clean_claim |                        30 |                 15 |                   12 |                  3 |            0.5       |                     0.0146306  |                          0.0362674  |                      -0.00700622 |
| Stack_mean_E1                                     | internal_candidate         | clean_internal_candidate        |                        29 |                 14 |                    9 |                  6 |            0.482759  |                     0.00137235 |                          0.0801305  |                       0          |
| BigMHC_IM                                         | caveated_public_comparator | caveated_public_not_clean_claim |                        30 |                 14 |                   14 |                  2 |            0.466667  |                     0          |                          0.107489   |                       0.0335532  |
| Stack_LR_inmaster_E3a                             | internal_candidate         | clean_internal_candidate        |                        29 |                 14 |                   12 |                  3 |            0.482759  |                     0          |                          0.138726   |                       0          |
| VQC                                               | bounded_fallback           | bounded_fallback_not_headline   |                        29 |                 10 |                   14 |                  5 |            0.344828  |                    -0.00720437 |                          0.0302225  |                      -0.101626   |
| Stack_median_E2                                   | internal_candidate         | clean_internal_candidate        |                        29 |                 10 |                   15 |                  4 |            0.344828  |                    -0.0159156  |                          0.00745875 |                       0.00650619 |
| kNN                                               | internal_candidate         | clean_internal_candidate        |                        29 |                 13 |                   16 |                  0 |            0.448276  |                    -0.0167897  |                          0.172817   |                      -0.0549693  |
| DeepImmuno                                        | caveated_public_comparator | caveated_public_not_clean_claim |                        29 |                 11 |                   15 |                  3 |            0.37931   |                    -0.0219207  |                          0.0766955  |                      -0.151984   |
| MultiTask_A_only                                  | internal_candidate         | clean_internal_candidate        |                        29 |                 12 |                   15 |                  2 |            0.413793  |                    -0.023461   |                          0.0962295  |                      -0.134162   |
| rule_gated_C_QK_structure                         | bounded_fallback           | bounded_fallback_not_headline   |                        21 |                  7 |                   12 |                  2 |            0.333333  |                    -0.0379278  |                         -0.089575   |                      -0.0408529  |
| MultiTask_AC                                      | internal_candidate         | clean_internal_candidate        |                        29 |                 14 |                   15 |                  0 |            0.482759  |                    -0.048964   |                          0.107341   |                      -0.122081   |
| prespecified_equal_weight_C_QK_no_anchor          | bounded_fallback           | bounded_fallback_not_headline   |                        21 |                  6 |                   12 |                  3 |            0.285714  |                    -0.0509     |                         -0.101556   |                      -0.0438519  |
| prespecified_lr_qk_quantum_only_w0.5              | bounded_fallback           | bounded_fallback_not_headline   |                        21 |                  2 |                   14 |                  5 |            0.0952381 |                    -0.052332   |                         -0.0271753  |                      -0.012797   |
| MultiTask_AB                                      | internal_candidate         | clean_internal_candidate        |                        29 |                 12 |                   17 |                  0 |            0.413793  |                    -0.0557468  |                          0.107341   |                      -0.15717    |

## Clean internal 후보만 보면

| method_name                                       |   n_wins_vs_anchor |   n_losses_vs_anchor |   win_rate_vs_anchor |   median_delta_AUPRC_vs_anchor |   source_heldout_median_delta_AUPRC |   hla_heldout_median_delta_AUPRC |
|:--------------------------------------------------|-------------------:|---------------------:|---------------------:|-------------------------------:|------------------------------------:|---------------------------------:|
| W7B_stacked                                       |                 16 |                   13 |             0.551724 |                     0.0946738  |                          0.172817   |                       0.0946738  |
| source_balanced_plus_pu_rf_train_prior_calibrated |                 12 |                   11 |             0.521739 |                     0.0689085  |                         -0.325587   |                       0.0697866  |
| W7A_full                                          |                 17 |                    7 |             0.586207 |                     0.0642506  |                          0.13852    |                       0.0642506  |
| pu_weighted_rf_train_prior_calibrated             |                 12 |                   11 |             0.521739 |                     0.0570796  |                         -0.311408   |                       0.0629524  |
| source_balanced_rf_train_prior_calibrated         |                 12 |                   11 |             0.521739 |                     0.056627   |                         -0.304438   |                       0.0798837  |
| Stack_mean_E1                                     |                 14 |                    9 |             0.482759 |                     0.00137235 |                          0.0801305  |                       0          |
| Stack_LR_inmaster_E3a                             |                 14 |                   12 |             0.482759 |                     0          |                          0.138726   |                       0          |
| Stack_median_E2                                   |                 10 |                   15 |             0.344828 |                    -0.0159156  |                          0.00745875 |                       0.00650619 |
| kNN                                               |                 13 |                   16 |             0.448276 |                    -0.0167897  |                          0.172817   |                      -0.0549693  |
| MultiTask_A_only                                  |                 12 |                   15 |             0.413793 |                    -0.023461   |                          0.0962295  |                      -0.134162   |
| MultiTask_AC                                      |                 14 |                   15 |             0.482759 |                    -0.048964   |                          0.107341   |                      -0.122081   |
| MultiTask_AB                                      |                 12 |                   17 |             0.413793 |                    -0.0557468  |                          0.107341   |                      -0.15717    |
| MultiTask_AD                                      |                 10 |                   17 |             0.344828 |                    -0.0574946  |                          0.0874993  |                      -0.15717    |
| Wave8_TCR_SelfSim_full                            |                  5 |                   24 |             0.166667 |                    -0.0579351  |                         -0.0539842  |                      -0.304945   |
| Wave8_TCR_SelfSim_no_exact                        |                 12 |                   18 |             0.4      |                    -0.0606859  |                         -0.048724   |                      -0.195884   |

## QK / fallback 계열 판단

| method_name                              |   n_wins_vs_anchor |   n_losses_vs_anchor |   win_rate_vs_anchor |   median_delta_AUPRC_vs_anchor |   source_heldout_median_delta_AUPRC |   hla_heldout_median_delta_AUPRC |
|:-----------------------------------------|-------------------:|---------------------:|---------------------:|-------------------------------:|------------------------------------:|---------------------------------:|
| W7A_QK_only                              |                 21 |                    3 |            0.724138  |                     0.088293   |                           0.088293  |                        0.0904311 |
| VQC                                      |                 10 |                   14 |            0.344828  |                    -0.00720437 |                           0.0302225 |                       -0.101626  |
| rule_gated_C_QK_structure                |                  7 |                   12 |            0.333333  |                    -0.0379278  |                          -0.089575  |                       -0.0408529 |
| prespecified_equal_weight_C_QK_no_anchor |                  6 |                   12 |            0.285714  |                    -0.0509     |                          -0.101556  |                       -0.0438519 |
| prespecified_lr_qk_quantum_only_w0.5     |                  2 |                   14 |            0.0952381 |                    -0.052332   |                          -0.0271753 |                       -0.012797  |
| GP_quantum                               |                  6 |                   20 |            0.206897  |                    -0.0619544  |                          -0.0619544 |                        0         |
| v0_fixed_late_fusion_C_QK_quantum        |                  0 |                   17 |            0         |                    -0.0931558  |                          -0.0960211 |                       -0.0791871 |
| prespecified_C_0.5_QK_quantum_0.5        |                  0 |                   17 |            0         |                    -0.0948908  |                          -0.0966102 |                       -0.0844871 |
| prespecified_rf_qk_quantum_only_w0.5     |                  0 |                   19 |            0         |                    -0.104584   |                          -0.10715   |                       -0.0845496 |
| source_prespecified_rf_qk_compact_w0.5   |                  5 |                   17 |            0.217391  |                    -0.113573   |                          -0.291362  |                       -0.0502046 |
| H_anchor_rf_plus_qk_compact_w0.5         |                  2 |                   21 |            0.0869565 |                    -0.115749   |                          -0.30821   |                       -0.113906  |
| nested_learned_gate_C_QK_structure       |                  2 |                   17 |            0.0952381 |                    -0.118569   |                          -0.158793  |                       -0.116568  |

판단: `W7A_QK_only`는 여러 matched split에서 Structure_LR을 이기지만, 역할은 fallback/fusion이다. 보고서에서는 이 결과를 ranking headline이나 quantum claim으로 쓰지 않는다.

## Public pretrained tool 판단

| method_name     |   n_wins_vs_anchor |   n_losses_vs_anchor |   win_rate_vs_anchor |   median_delta_AUPRC_vs_anchor | method_disposition              |
|:----------------|-------------------:|---------------------:|---------------------:|-------------------------------:|:--------------------------------|
| MHCflurry       |                 15 |                   12 |             0.5      |                      0.0146306 | caveated_public_not_clean_claim |
| BigMHC_IM       |                 14 |                   14 |             0.466667 |                      0         | caveated_public_not_clean_claim |
| DeepImmuno      |                 11 |                   15 |             0.37931  |                     -0.0219207 | caveated_public_not_clean_claim |
| TSCAPE_TITANiAN |                 10 |                   19 |             0.344828 |                     -0.0735363 | caveated_public_not_clean_claim |
| MHCnuggets_2    |                  9 |                   18 |             0.310345 |                     -0.0740465 | caveated_public_not_clean_claim |
| TransPHLA       |                  5 |                   22 |             0.172414 |                     -0.104131  | caveated_public_not_clean_claim |
| NetMHCstabpan   |                 10 |                   19 |             0.344828 |                     -0.116494  | caveated_public_not_clean_claim |
| NetMHCpan_4.1   |                  7 |                   23 |             0.233333 |                     -0.133333  | caveated_public_not_clean_claim |
| PRIME           |                  5 |                   23 |             0.166667 |                     -0.140627  | caveated_public_not_clean_claim |

| method_name     | training_overlap_audit_status           | clean_comparator_allowed_after_audit   | reviewer_disposition            |
|:----------------|:----------------------------------------|:---------------------------------------|:--------------------------------|
| BigMHC_IM       | unresolved_no_row_level_training_corpus | False                                  | caveated_public_comparator_only |
| DeepImmuno      | unresolved_no_row_level_training_corpus | False                                  | caveated_public_comparator_only |
| MHCflurry       | unresolved_no_row_level_training_corpus | False                                  | caveated_public_comparator_only |
| MHCnuggets_2    | unresolved_no_row_level_training_corpus | False                                  | caveated_public_comparator_only |
| NetMHCpan_4.1   | unresolved_no_row_level_training_corpus | False                                  | caveated_public_comparator_only |
| NetMHCstabpan   | unresolved_no_row_level_training_corpus | False                                  | caveated_public_comparator_only |
| PRIME           | unresolved_no_row_level_training_corpus | False                                  | caveated_public_comparator_only |
| TSCAPE_TITANiAN | unresolved_no_row_level_training_corpus | False                                  | caveated_public_comparator_only |
| TransPHLA       | unresolved_no_row_level_training_corpus | False                                  | caveated_public_comparator_only |

판단: public tool이 일부 split에서 좋아도 row-level training-corpus overlap audit 전까지는 caveated comparator다.

## 어디서 틀리는지

- source shift: 일부 RF/PU/source-balanced 계열은 overall AUPRC가 높지만 source-heldout median delta가 크게 음수다.
- low prevalence: TESLA-like 저 prevalence slice에서 false-positive pressure가 커진다.
- rare / underrepresented HLA: HLA-heldout 또는 Korean-HLA focus에서 method별 collapse가 갈린다.
- public/internal disagreement: public pretrained predictor가 높고 clean internal support가 약하면 manual review 또는 abstain.
- patient gate: PAAD/THCA patient-level fields가 부족해서 현재 patient-gated demo는 clinical-use가 아니라 research triage only다.

## Challenge / distribution audit

| challenge_axis                       |   n_unique_candidates |   positive_prevalence | recommended_split_contract                          | success_metric                                                                     |
|:-------------------------------------|----------------------:|----------------------:|:----------------------------------------------------|:-----------------------------------------------------------------------------------|
| missed_positive_rescue_watchlist     |                   120 |              1        | HLA_heldout + source_heldout                        | recall@20, mean positive rank, missed-positive rate                                |
| korean_hla_focus_stress              |                   120 |              0.95     | korean_hla_focus                                    | allele-level AUPRC, top10 precision, positive rank                                 |
| external_holdout_fragility           |                   120 |              0.791667 | source_heldout_external                             | source-heldout AUPRC, top-k survival, confidence-risk curve                        |
| patient_gate_metadata_blocker        |                   120 |              0.966667 | patient_gated_PAAD_THCA_demo                        | metadata completion rate, gate pass/fail concordance, patient-gated rank stability |
| rare_hla_support_gap                 |                   120 |              0.991667 | HLA_heldout + Korean_HLA_focus_if_applicable        | allele-specific AUPRC, top10 precision, calibration ECE                            |
| high_score_claim_blocked             |                   120 |              0.858333 | exact_phla_holdout + near_peptide_cluster_holdout   | drop in top10 precision after removing overlap-risk rows                           |
| low_prevalence_false_positive_stress |                   120 |              0        | low_prevalence_heldout                              | top10 precision, top20 precision, false-positive rate among negatives              |
| public_internal_disagreement         |                    43 |              0.767442 | public_overlap_audit + clean_internal_only_ablation | delta AUPRC/top10 between public-inclusive and clean-internal-only views           |
| claim_safe_priority_review           |                     2 |              1        | exact_phla + source_heldout + HLA_heldout           | AUPRC, top10 precision, calibration ECE, zero unresolved overlap                   |

## Patient-gated blocker

| field                     | gate_group      | required_for   | metadata_status    | reviewer_safe_default_if_missing         |
|:--------------------------|:----------------|:---------------|:-------------------|:-----------------------------------------|
| patient_id                | identity        | both           | present_but_empty  | unknown patient context                  |
| cancer_type               | disease_context | both           | present_but_empty  | no disease-specific confidence           |
| disease_context           | disease_context | both           | present_but_empty  | disease gate unknown                     |
| tumor_burden_context      | disease_context | PAAD           | absent_from_master | PAAD vaccine-first confidence capped     |
| disease_status            | disease_context | both           | absent_from_master | disease gate unknown                     |
| tumor_stage               | disease_context | both           | present_but_empty  | stage-based gate unknown                 |
| histology                 | disease_context | THCA           | absent_from_master | routine low-risk PTC cannot be separated |
| driver_mutation           | antigen         | both           | absent_from_master | driver-aware strategy unavailable        |
| kras_mutation             | antigen         | PAAD           | absent_from_master | KRAS-public branch unavailable           |
| braf_v600e                | disease_context | THCA           | absent_from_master | BRAF-specific THCA strategy unavailable  |
| ret_ntrk_fusion           | disease_context | THCA           | absent_from_master | fusion-targeted backbone unknown         |
| tert_status               | disease_context | THCA           | absent_from_master | risk context incomplete                  |
| hla_loh                   | presentation    | both           | present_but_empty  | presentation escape unknown              |
| b2m_status                | presentation    | both           | present_but_empty  | presentation escape unknown              |
| hla_expression            | presentation    | both           | absent_from_master | presentation gate unknown                |
| antigen_processing_status | presentation    | both           | present_but_empty  | processing gate unknown                  |
| expression_tpm            | antigen         | both           | present_but_empty  | antigen evidence weak                    |
| mutant_expression         | antigen         | both           | present_but_empty  | antigen evidence weak                    |
| vaf                       | antigen         | both           | present_but_empty  | clonality unknown                        |
| clonality                 | antigen         | both           | present_but_empty  | antigen durability unknown               |
| immune_context_score      | immune_context  | both           | present_but_empty  | immune gate unknown                      |
| tls_score                 | immune_context  | THCA           | present_but_empty  | THCA immune-active branch weak           |
| ifng_score                | immune_context  | both           | present_but_empty  | immune gate unknown                      |
| cytolytic_score           | immune_context  | both           | present_but_empty  | immune gate unknown                      |
| autoimmune_history        | safety          | both           | absent_from_master | safety gate unknown                      |

## 바로 볼 파일

- HTML dossier: `project/papers_hub_2026_05_04/clean_neobench_barneo_dossier_2026_05_09.html`
- Visual dashboard: `project/papers_hub_2026_05_04/clean_neobench_visual_dashboard_2026_05_10.html`
- Win/loss report: `project/results/clean_neobench_barneo_2026_05_09/CLEAN_NEOBENCH_WINLOSS_REPORT_KR.md`
- Win/loss table: `project/results/clean_neobench_barneo_2026_05_09/clean_neobench_winloss_method_summary.tsv`
- Split winner board: `project/results/clean_neobench_barneo_2026_05_09/clean_neobench_split_winner_board.tsv`
- Distribution error audit: `project/results/clean_neobench_barneo_2026_05_09/CLEAN_NEOBENCH_DISTRIBUTION_ERROR_AUDIT_KR.md`
- Challenge pack: `project/results/clean_neobench_barneo_2026_05_09/CLEAN_NEOBENCH_CHALLENGE_PACK_KR.md`

## 다음 우선순위

1. Public tool training-corpus row-level overlap audit.
2. Source-heldout collapse가 큰 RF/PU/source-balanced 모델의 source reweighting 또는 abstention 강화.
3. Rare/Korean-HLA focus split에서 W7B/W7A/source-balanced RF를 다시 stress-test.
4. PAAD/THCA patient-gated demo에 실제 disease timing, expression/clonality, immune context, safety metadata 연결.
5. MHC-II는 별도 benchmark로 분리.

## Claim boundary

Allowed: leakage-aware AI neoantigen predictor benchmarking framework, benchmark-adaptive reliability ranking, calibrated candidate prioritization with abstention, reviewer-safe comparison, research triage framework.

Forbidden: clinical vaccine selection, new SOTA predictor, external validation proven, quantum advantage, public tools as clean baselines without overlap audit, Class I and Class II unified predictor.
