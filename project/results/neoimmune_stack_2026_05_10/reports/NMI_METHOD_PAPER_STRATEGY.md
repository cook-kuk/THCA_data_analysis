# NMI Method Paper Strategy

**NMI = Neoantigen Multimodal Immunogenicity.**

## Decision
Use NMI as the method-paper core, not the public-predictor-containing GA/RL production controller.

## Why this is cleaner
- No NetMHCpan, MHCflurry, BigMHC, PRIME, HLApollo, BAR-Neo, KG-GA, or product-impact scores are used as NMI training features.
- NMI uses peptide/HLA/self/TCR/tumor context plus local branches only: Structure, ESM2, Wave8/TCR-self-similarity, quantum/QK, and W7 stacks.
- Public predictors are reported only as frozen comparators.

## Best current NMI row
- `NMI all-branch mean` / `nmi_eligible_locked_formula`: AUPRC=0.801, AUROC=0.852, patient_hit@34=1.000.

## Leaderboard
| track                  | model                                                  | split                           | features                                                       | clean_track_allowed   |     n |   positives |     AUROC |    AUPRC |   Precision@10 |   Precision@20 |   Precision@34 |   Recall@20 |   Recall@34 |   patient_hit_rate@20 |   patient_recall@20 |   patient_hit_rate@34 |   patient_recall@34 |   patients_evaluated |   calibration_brier |   Recall@10 |   Precision@50 |   Recall@50 | score_col                  |
|:-----------------------|:-------------------------------------------------------|:--------------------------------|:---------------------------------------------------------------|:----------------------|------:|------------:|----------:|---------:|---------------:|---------------:|---------------:|------------:|------------:|----------------------:|--------------------:|----------------------:|--------------------:|---------------------:|--------------------:|------------:|---------------:|------------:|:---------------------------|
| nmi_clean_method_track | NMI all-branch mean                                    | nmi_eligible_locked_formula     | locked clean local branch formula on NMI-eligible rows         | True                  |   319 |         136 | 0.852178  | 0.801046 |            0.8 |           0.9  |      0.911765  | 0.132353    | 0.227941    |              1        |          0.132353   |              1        |           0.227941  |                    1 |            0.18316  | 0.0588235   |           0.92 | 0.338235    | NMI_all_branch_mean        |
| nmi_clean_method_track | NMI TCR-QK-Structure                                   | nmi_eligible_locked_formula     | locked TCR/quantum/structure/ESM2 formula on NMI-eligible rows | True                  |   311 |         130 | 0.846026  | 0.800592 |            1   |           0.95 |      0.911765  | 0.146154    | 0.238462    |              1        |          0.146154   |              1        |           0.238462  |                    1 |            0.181379 | 0.0769231   |           0.88 | 0.338462    | NMI_tcr_qk_structure       |
| nmi_clean_method_track | NMI equal top5                                         | nmi_eligible_locked_formula     | W7B/ESM2/QK/full/Structure equal-weight formula                | True                  |   319 |         136 | 0.808824  | 0.765184 |            1   |           0.95 |      0.911765  | 0.139706    | 0.227941    |              1        |          0.139706   |              1        |           0.227941  |                    1 |            0.178526 | 0.0735294   |           0.84 | 0.308824    | NMI_equal_top5             |
| nmi_clean_method_track | NMI W7B-ESM2-QK                                        | nmi_eligible_locked_formula     | W7B/ESM2/QK weighted formula                                   | True                  |   311 |         130 | 0.818997  | 0.765107 |            1   |           0.95 |      0.882353  | 0.146154    | 0.230769    |              1        |          0.146154   |              1        |           0.230769  |                    1 |            0.172651 | 0.0769231   |           0.86 | 0.330769    | NMI_w7b_esm2_qk            |
| nmi_clean_method_track | same-set public comparator: BigMHC_IM                  | nmi_eligible_same_rows          | frozen public predictor on NMI-eligible rows                   | False                 |   319 |         136 | 0.747308  | 0.626858 |            0.5 |           0.6  |      0.676471  | 0.0882353   | 0.169118    |              1        |          0.0882353  |              1        |           0.169118  |                    1 |            0.219336 | 0.0367647   |           0.74 | 0.272059    | BigMHC_IM                  |
| nmi_clean_method_track | same-set public comparator: MHCflurry_2.0_presentation | nmi_eligible_same_rows          | frozen public predictor on NMI-eligible rows                   | False                 |   319 |         136 | 0.641916  | 0.592932 |            0.8 |           0.8  |      0.794118  | 0.117647    | 0.198529    |              1        |          0.117647   |              1        |           0.198529  |                    1 |            0.401156 | 0.0588235   |           0.64 | 0.235294    | MHCflurry_2.0_presentation |
| nmi_clean_method_track | same-set public comparator: NetMHCpan_4.1_EL           | nmi_eligible_same_rows          | frozen public predictor on NMI-eligible rows                   | False                 |   319 |         136 | 0.56688   | 0.493205 |            0.5 |           0.55 |      0.470588  | 0.0808824   | 0.117647    |              1        |          0.0808824  |              1        |           0.117647  |                    1 |            0.573532 | 0.0367647   |           0.52 | 0.191176    | NetMHCpan_4.1_EL           |
| nmi_clean_method_track | same-set public comparator: PRIME                      | nmi_eligible_same_rows          | frozen public predictor on NMI-eligible rows                   | False                 |   319 |         136 | 0.519949  | 0.439568 |            0.4 |           0.45 |      0.411765  | 0.0661765   | 0.102941    |              1        |          0.0661765  |              1        |           0.102941  |                    1 |            0.452571 | 0.0294118   |           0.44 | 0.161765    | PRIME                      |
| nmi_clean_method_track | public comparator: BigMHC_IM                           | retrospective_public_comparator | frozen public predictor score                                  | False                 |  2902 |        1550 | 0.679993  | 0.708425 |            0.7 |           0.8  |      0.823529  | 0.0103226   | 0.0180645   |              1        |          0.0103226  |              1        |           0.0180645 |                    1 |            0.306495 | 0.00451613  |           0.8  | 0.0258065   | BigMHC_IM                  |
| nmi_clean_method_track | public comparator: MHCflurry_2.0_presentation          | retrospective_public_comparator | frozen public predictor score                                  | False                 |  2755 |        1493 | 0.656119  | 0.700977 |            0.8 |           0.85 |      0.764706  | 0.0113865   | 0.0174146   |              1        |          0.0113865  |              1        |           0.0174146 |                    1 |            0.323978 | 0.00535834  |           0.76 | 0.0254521   | MHCflurry_2.0_presentation |
| nmi_clean_method_track | public comparator: PRIME                               | retrospective_public_comparator | frozen public predictor score                                  | False                 |  2899 |        1550 | 0.589324  | 0.592204 |            0.4 |           0.5  |      0.470588  | 0.00645161  | 0.0103226   |              1        |          0.00645161 |              1        |           0.0103226 |                    1 |            0.265079 | 0.00258065  |           0.42 | 0.0135484   | PRIME                      |
| nmi_clean_method_track | public comparator: NetMHCpan_4.1_EL                    | retrospective_public_comparator | frozen public predictor score                                  | False                 |   319 |         136 | 0.56688   | 0.493205 |            0.5 |           0.55 |      0.470588  | 0.0808824   | 0.117647    |              1        |          0.0808824  |              1        |           0.117647  |                    1 |            0.573532 | 0.0367647   |           0.52 | 0.191176    | NetMHCpan_4.1_EL           |
| nmi_clean_method_track | NMI-branch LR                                          | source_grouped_oof_fallback     | transparent own branch fusion, source-grouped OOF              | True                  | 22090 |        3988 | 0.346265  | 0.133619 |            0.3 |           0.4  |      0.323529  | 0.00200602  | 0.00275827  |              0.712644 |          0.734113   |              0.758621 |           0.756254  |                   87 |            0.774782 | 0.000752257 |           0.24 | 0.00300903  | NMI_branch_LR_source_oof   |
| nmi_clean_method_track | NMI-transparent formula                                | source_grouped_oof_fallback     | fixed clean branch-weight formula                              | True                  | 22090 |        3988 | 0.22977   | 0.131393 |            0.8 |           0.85 |      0.852941  | 0.00426279  | 0.00727182  |              0.701149 |          0.712278   |              0.735632 |           0.757045  |                   87 |            0.251599 | 0.00200602  |           0.8  | 0.0100301   | NMI_transparent_score      |
| nmi_clean_method_track | NMI-core LR                                            | source_grouped_oof_fallback     | peptide/HLA/self/TCR/tumor context only                        | True                  | 22090 |        3988 | 0.288904  | 0.120352 |            0.1 |           0.15 |      0.235294  | 0.000752257 | 0.00200602  |              0.712644 |          0.73159    |              0.770115 |           0.777599  |                   87 |            0.765906 | 0.000250752 |           0.18 | 0.00225677  | NMI_core_LR_source_oof     |
| nmi_clean_method_track | NMI-clean consensus                                    | source_grouped_oof_fallback     | peptide/HLA/self/TCR/tumor context + own branch fusion         | True                  | 22090 |        3988 | 0.0976974 | 0.10542  |            0.8 |           0.5  |      0.323529  | 0.00250752  | 0.00275827  |              0.724138 |          0.735685   |              0.758621 |           0.759463  |                   87 |            0.521807 | 0.00200602  |           0.24 | 0.00300903  | NMI_clean_consensus        |
| nmi_clean_method_track | NMI-branch ExtraTrees                                  | source_grouped_oof_fallback     | own branch fusion, source-grouped OOF                          | True                  | 22090 |        3988 | 0.0964827 | 0.104758 |            0.7 |           0.55 |      0.382353  | 0.00275827  | 0.00325978  |              0.747126 |          0.760939   |              0.804598 |           0.79989   |                   87 |            0.575714 | 0.00175527  |           0.26 | 0.00325978  | NMI_branch_ET_source_oof   |
| nmi_clean_method_track | NMI-branch RF                                          | source_grouped_oof_fallback     | own branch fusion, source-grouped OOF                          | True                  | 22090 |        3988 | 0.0797281 | 0.101116 |            0   |           0.05 |      0.0294118 | 0.000250752 | 0.000250752 |              0.735632 |          0.743305   |              0.747126 |           0.753721  |                   87 |            0.61532  | 0           |           0.06 | 0.000752257 | NMI_branch_RF_source_oof   |

## Branch feature manifest
| branch_feature                           |   n_available | allowed_in_nmi   | reason                                                  |
|:-----------------------------------------|--------------:|:-----------------|:--------------------------------------------------------|
| ESM2_Bayesian                            |           311 | True             | local/own branch score; public predictor tokens blocked |
| GP_quantum                               |           319 | True             | local/own branch score; public predictor tokens blocked |
| Stack_LR_OOF_E3b                         |           106 | True             | local/own branch score; public predictor tokens blocked |
| Stack_LR_inmaster_E3a                    |           319 | True             | local/own branch score; public predictor tokens blocked |
| Stack_mean_E1                            |           319 | True             | local/own branch score; public predictor tokens blocked |
| Stack_median_E2                          |           319 | True             | local/own branch score; public predictor tokens blocked |
| Structure_LR                             |           311 | True             | local/own branch score; public predictor tokens blocked |
| VQC                                      |           319 | True             | local/own branch score; public predictor tokens blocked |
| W7A_QK_only                              |           319 | True             | local/own branch score; public predictor tokens blocked |
| W7A_full                                 |           319 | True             | local/own branch score; public predictor tokens blocked |
| W7B_stacked                              |           319 | True             | local/own branch score; public predictor tokens blocked |
| Wave8_SelfSim_full                       |           319 | True             | local/own branch score; public predictor tokens blocked |
| Wave8_SelfSim_no_exact                   |           319 | True             | local/own branch score; public predictor tokens blocked |
| Wave8_TCR_SelfSim_full                   |           319 | True             | local/own branch score; public predictor tokens blocked |
| Wave8_TCR_SelfSim_no_exact               |           319 | True             | local/own branch score; public predictor tokens blocked |
| Wave8_TCR_motif_only                     |           319 | True             | local/own branch score; public predictor tokens blocked |
| Wave8_TCR_only                           |           319 | True             | local/own branch score; public predictor tokens blocked |
| hard_decoy_rule_aux_C_QK_no_anchor       |            89 | True             | local/own branch score; public predictor tokens blocked |
| nested_learned_gate_C_QK_structure       |            89 | True             | local/own branch score; public predictor tokens blocked |
| nested_lr_qk_no_anchor_train_selected    |            89 | True             | local/own branch score; public predictor tokens blocked |
| nested_lr_qk_quantum_only_train_selected |            89 | True             | local/own branch score; public predictor tokens blocked |
| nested_rf_qk_no_anchor_train_selected    |            89 | True             | local/own branch score; public predictor tokens blocked |
| nested_rf_qk_quantum_only_train_selected |            89 | True             | local/own branch score; public predictor tokens blocked |
| prespecified_C_0.5_QK_quantum_0.5        |            89 | True             | local/own branch score; public predictor tokens blocked |
| prespecified_equal_weight_C_QK_no_anchor |            89 | True             | local/own branch score; public predictor tokens blocked |
| prespecified_lr_qk_no_anchor_w0.5        |            89 | True             | local/own branch score; public predictor tokens blocked |
| prespecified_lr_qk_quantum_only_w0.5     |            89 | True             | local/own branch score; public predictor tokens blocked |
| prespecified_rf_qk_no_anchor_w0.5        |            89 | True             | local/own branch score; public predictor tokens blocked |
| prespecified_rf_qk_quantum_only_w0.5     |            89 | True             | local/own branch score; public predictor tokens blocked |
| qk_clean_no_tcr_gamma1                   |            89 | True             | local/own branch score; public predictor tokens blocked |
| qk_no_anchor                             |            89 | True             | local/own branch score; public predictor tokens blocked |
| qk_no_anchor_gamma1                      |            89 | True             | local/own branch score; public predictor tokens blocked |
| qk_quantum_only                          |            89 | True             | local/own branch score; public predictor tokens blocked |
| qk_quantum_only_gamma1                   |            89 | True             | local/own branch score; public predictor tokens blocked |
| v0_QK_no_anchor_fallback                 |            89 | True             | local/own branch score; public predictor tokens blocked |
| v0_QK_quantum_only_fallback              |            89 | True             | local/own branch score; public predictor tokens blocked |
| v0_fixed_late_fusion_C_QK_no_anchor      |            89 | True             | local/own branch score; public predictor tokens blocked |
| v0_fixed_late_fusion_C_QK_quantum        |            89 | True             | local/own branch score; public predictor tokens blocked |

## Method-paper claim boundary
- Allowed: NMI is a leakage-aware clean multimodal immunogenicity ranker evaluated retrospectively on public labels.
- Allowed: NMI separates clean-method evidence from production stacking and public predictor comparators.
- Forbidden: clinical vaccine efficacy, verified antigen presentation without MS, immunogenicity without T-cell assay labels, or universal superiority across all cohorts.

## Paper title candidates
1. NMI: a leakage-aware multimodal immunogenicity ranker for patient-level neoantigen prioritization
2. CLEAN-NMI: separating TCR-visible immunogenicity from HLA presentation in neoantigen ranking
3. Neoantigen Multimodal Immunogenicity modeling under leakage-controlled patient-level evaluation

## Next method-paper gate
Run NMI on a locked external patient-level set with pre-registered splits: leave-study-out, leave-HLA-supertype-out, and hospital-heldout when available.
