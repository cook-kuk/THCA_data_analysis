# CROSS-Neo v2.1 Selective Ensemble/Fallback Pipeline

Status: **policy hypothesis to freeze before true external testing**. This is not an external-validation result.

## Routing Logic

| split_or_regime                      | conservative                            | topk                                    | reason                                                                                           | abstain   |
|:-------------------------------------|:----------------------------------------|:----------------------------------------|:-------------------------------------------------------------------------------------------------|:----------|
| exact_peptide_hla_holdout            | fast_nested_esm2_qk_gate                | fast_nested_esm2_qk_gate                | exact/reference-clean setting where ESM2+QK gate improved locked AUPRC                           | False     |
| near_peptide_cluster_holdout         | rule_gate_rf_qk_fallback_train_selected | rule_gate_rf_qk_fallback_train_selected | near-neighbor OOD setting where ESM2 fusion harmed AUPRC                                         | False     |
| hla_stratified_group_5fold           | v2_cf_plm_lr                            | fast_nested_esm2_qk_gate                | clean counterfactual model wins AUPRC; ESM2 gate is only top-k auxiliary                         | False     |
| hla_supertype_heldout                | v2_multimodal_lr                        | v2_multimodal_lr                        | multimodal LR wins HLA-supertype AUPRC; ESM2 does not improve enough                             | False     |
| source_heldout_CEDAR                 | v2_multimodal_lr                        | v2_multimodal_lr                        | CEDAR is high-prevalence descriptive source; avoid diagnostic-only PLM pilot as claim model      | False     |
| source_heldout_NEPdb                 | v2_groupdro_proxy_cf_lr                 | v2_groupdro_proxy_cf_lr                 | source-robust proxy recovers NEPdb top-k better than clean anchor                                | False     |
| source_heldout_TESLA_mmc4            | v2_multimodal_esm2_650m_lr_fast         | v2_multimodal_esm2_650m_lr_fast         | severe low-prevalence source shift; ESM2-650M gives best top10/top20 recovery                    | False     |
| source_heldout_TESLA_mmc7_validation | v2_groupdro_proxy_cf_esm2_35m_lr_fast   | v2_groupdro_proxy_cf_esm2_35m_lr_fast   | extreme low-positive source shift; only weak top20 recovery, require abstention caveat           | True      |
| low_prevalence_stress_split          | v2_multimodal_lr                        | v2_multimodal_lr                        | low-prevalence aggregate stress; conservative multimodal model avoids overfitting to ESM2 rescue | True      |

## Internal Retrospective Comparison

| split_name                           | pre_policy_best                         |   pre_policy_AUPRC |   pre_policy_top10 | conservative_expert                     |   conservative_AUPRC |   conservative_top10 | topk_expert                             |   topk_AUPRC |   topk_top10 | abstention_recommended   | routing_reason                                                                                   |
|:-------------------------------------|:----------------------------------------|-------------------:|-------------------:|:----------------------------------------|---------------------:|---------------------:|:----------------------------------------|-------------:|-------------:|:-------------------------|:-------------------------------------------------------------------------------------------------|
| exact_peptide_hla_holdout            | fast_nested_esm2_qk_gate                |          0.636885  |                0.7 | fast_nested_esm2_qk_gate                |            0.636885  |                  0.7 | fast_nested_esm2_qk_gate                |    0.636885  |          0.7 | False                    | exact/reference-clean setting where ESM2+QK gate improved locked AUPRC                           |
| near_peptide_cluster_holdout         | rule_gate_rf_qk_fallback_train_selected |          0.518669  |                0.7 | rule_gate_rf_qk_fallback_train_selected |            0.518669  |                  0.7 | rule_gate_rf_qk_fallback_train_selected |    0.518669  |          0.7 | False                    | near-neighbor OOD setting where ESM2 fusion harmed AUPRC                                         |
| hla_stratified_group_5fold           | v2_cf_plm_lr                            |          0.563197  |                0.5 | v2_cf_plm_lr                            |            0.563197  |                  0.5 | fast_nested_esm2_qk_gate                |    0.533327  |          0.7 | False                    | clean counterfactual model wins AUPRC; ESM2 gate is only top-k auxiliary                         |
| hla_supertype_heldout                | v2_multimodal_lr                        |          0.58542   |                0.6 | v2_multimodal_lr                        |            0.58542   |                  0.6 | v2_multimodal_lr                        |    0.58542   |          0.6 | False                    | multimodal LR wins HLA-supertype AUPRC; ESM2 does not improve enough                             |
| source_heldout_CEDAR                 | v2_plm_lr_frozen_hash_pilot             |          0.957745  |                1   | v2_multimodal_lr                        |            0.946896  |                  1   | v2_multimodal_lr                        |    0.946896  |          1   | False                    | CEDAR is high-prevalence descriptive source; avoid diagnostic-only PLM pilot as claim model      |
| source_heldout_NEPdb                 | source_qk_compact_gamma1                |          0.406426  |                0.4 | v2_groupdro_proxy_cf_lr                 |            0.281349  |                  0.5 | v2_groupdro_proxy_cf_lr                 |    0.281349  |          0.5 | False                    | source-robust proxy recovers NEPdb top-k better than clean anchor                                |
| source_heldout_TESLA_mmc4            | v2_plm_lr_frozen_hash_pilot             |          0.0923091 |                0   | v2_multimodal_esm2_650m_lr_fast         |            0.0834861 |                  0.2 | v2_multimodal_esm2_650m_lr_fast         |    0.0834861 |          0.2 | False                    | severe low-prevalence source shift; ESM2-650M gives best top10/top20 recovery                    |
| source_heldout_TESLA_mmc7_validation | v2_groupdro_proxy_cf_esm2_35m_lr_fast   |          0.0387279 |                0   | v2_groupdro_proxy_cf_esm2_35m_lr_fast   |            0.0387279 |                  0   | v2_groupdro_proxy_cf_esm2_35m_lr_fast   |    0.0387279 |          0   | True                     | extreme low-positive source shift; only weak top20 recovery, require abstention caveat           |
| low_prevalence_stress_split          | v2_esm2_650m_lr_fast                    |          0.0814706 |                0.1 | v2_multimodal_lr                        |            0.0615867 |                  0.1 | v2_multimodal_lr                        |    0.0615867 |          0.1 | True                     | low-prevalence aggregate stress; conservative multimodal model avoids overfitting to ESM2 rescue |

## What This Means

- Use ESM2+QK gate only in exact/reference-clean regimes.
- Do not use ESM2 gate in near-peptide cluster holdout; it harms AUPRC.
- Use source-robust fallback for NEPdb-like source shift.
- Use ESM2-650M fallback for TESLA_mmc4-like severe low-prevalence source shift.
- For TESLA_mmc7-like extreme low-positive settings, report abstention/top20 only unless new evidence appears.

## Claim Boundary

- Allowed after freezing and re-testing: selective ensemble pipeline improves internal exact/top-k rescue and defines abstention for source-shift.
- Forbidden now: external validation, SOTA predictor, clinical vaccine selection, quantum advantage.

Workbook: `/home/seungho/personal/THCA_data_analysis/project/results/cross_neo_v2_sota_sprint_2026_05_09/CROSS_Neo_v2_1_selective_ensemble_pipeline.xlsx`
