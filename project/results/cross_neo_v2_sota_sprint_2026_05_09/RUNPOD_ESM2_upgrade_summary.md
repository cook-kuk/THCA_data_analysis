# RunPod ESM2 Upgrade Summary

RunPod A6000 was used to generate frozen ESM2-35M, ESM2-150M, and ESM2-650M projected features and fast train-fold-safe logistic evaluations.

## Primary Locked Splits

| split_name                   | base_best_model                         |   base_best_AUPRC |   base_best_top10 | runpod_esm2_best_model          |   runpod_esm2_best_AUPRC |   runpod_esm2_best_top10 | combined_best_model                     |   combined_best_AUPRC |   combined_best_top10 |
|:-----------------------------|:----------------------------------------|------------------:|------------------:|:--------------------------------|-------------------------:|-------------------------:|:----------------------------------------|----------------------:|----------------------:|
| exact_peptide_hla_holdout    | prespecified_rf_qk_no_anchor_w0.5       |          0.586965 |               0.7 | v2_cf_esm2_150m_lr_fast         |                 0.540566 |                      0.8 | prespecified_rf_qk_no_anchor_w0.5       |              0.586965 |                   0.7 |
| near_peptide_cluster_holdout | rule_gate_rf_qk_fallback_train_selected |          0.518669 |               0.7 | v2_cf_esm2_650m_lr_fast         |                 0.420486 |                      0.7 | rule_gate_rf_qk_fallback_train_selected |              0.518669 |                   0.7 |
| hla_stratified_group_5fold   | v2_cf_plm_lr                            |          0.563197 |               0.5 | v2_multimodal_esm2_650m_lr_fast |                 0.476056 |                      0.6 | v2_cf_plm_lr                            |              0.563197 |                   0.5 |
| hla_supertype_heldout        | v2_multimodal_lr                        |          0.58542  |               0.6 | v2_multimodal_esm2_650m_lr_fast |                 0.514148 |                      0.6 | v2_multimodal_lr                        |              0.58542  |                   0.6 |

## Source-Heldout Best Rows

| split_name                           | model_name                              |   n |   n_pos |   prevalence |     AUPRC |   top10_precision |   top20_precision | claim_status                  |
|:-------------------------------------|:----------------------------------------|----:|--------:|-------------:|----------:|------------------:|------------------:|:------------------------------|
| source_heldout_CEDAR                 | v2_plm_lr_frozen_hash_pilot             | 909 |     851 |    0.936194  | 0.957745  |               1   |              1    | diagnostic                    |
| source_heldout_CEDAR                 | v2_multimodal_lr                        | 909 |     851 |    0.936194  | 0.946896  |               1   |              1    | reviewer_safe_internal_locked |
| source_heldout_CEDAR                 | v2_cf_plm_lr                            | 909 |     851 |    0.936194  | 0.946877  |               1   |              1    | reviewer_safe_internal_locked |
| source_heldout_CEDAR                 | v2_rule_moe_cf_plm_structure            | 909 |     851 |    0.936194  | 0.946425  |               1   |              1    | reviewer_safe_internal_locked |
| source_heldout_CEDAR                 | v2_cf_plm_rf_hla_ranknorm               | 909 |     851 |    0.936194  | 0.944727  |               1   |              1    | reviewer_safe_internal_locked |
| source_heldout_NEPdb                 | v2_groupdro_proxy_cf_lr                 | 572 |     151 |    0.263986  | 0.281349  |               0.5 |              0.45 | reviewer_safe_internal_locked |
| source_heldout_NEPdb                 | source_qk_compact_gamma1                | 886 |     354 |    0.399549  | 0.406426  |               0.4 |              0.4  | diagnostic                    |
| source_heldout_NEPdb                 | v2_multimodal_lr                        | 572 |     151 |    0.263986  | 0.30528   |               0.3 |              0.4  | reviewer_safe_internal_locked |
| source_heldout_NEPdb                 | v2_cf_plm_lr                            | 572 |     151 |    0.263986  | 0.305266  |               0.3 |              0.4  | reviewer_safe_internal_locked |
| source_heldout_NEPdb                 | v2_plm_lr_frozen_hash_pilot             | 572 |     151 |    0.263986  | 0.277885  |               0.3 |              0.25 | diagnostic                    |
| source_heldout_TESLA_mmc4            | v2_multimodal_esm2_650m_lr_fast         | 605 |      37 |    0.061157  | 0.0834861 |               0.2 |              0.15 | reviewer_safe_internal_locked |
| source_heldout_TESLA_mmc4            | v2_cf_esm2_650m_lr_fast                 | 605 |      37 |    0.061157  | 0.0834368 |               0.2 |              0.15 | reviewer_safe_internal_locked |
| source_heldout_TESLA_mmc4            | v2_esm2_650m_lr_fast                    | 605 |      37 |    0.061157  | 0.0894906 |               0.2 |              0.1  | reviewer_safe_internal_locked |
| source_heldout_TESLA_mmc4            | v2_groupdro_proxy_cf_esm2_650m_lr_fast  | 605 |      37 |    0.061157  | 0.081167  |               0.2 |              0.1  | reviewer_safe_internal_locked |
| source_heldout_TESLA_mmc4            | v2_source_balanced_cf_esm2_650m_lr_fast | 605 |      37 |    0.061157  | 0.081167  |               0.2 |              0.1  | reviewer_safe_internal_locked |
| source_heldout_TESLA_mmc7_validation | v2_groupdro_proxy_cf_esm2_35m_lr_fast   | 310 |       4 |    0.0129032 | 0.0387279 |               0   |              0.05 | reviewer_safe_internal_locked |
| source_heldout_TESLA_mmc7_validation | v2_source_balanced_cf_esm2_35m_lr_fast  | 310 |       4 |    0.0129032 | 0.0387279 |               0   |              0.05 | reviewer_safe_internal_locked |
| source_heldout_TESLA_mmc7_validation | v2_cf_esm2_35m_lr_fast                  | 310 |       4 |    0.0129032 | 0.0371271 |               0   |              0.05 | reviewer_safe_internal_locked |
| source_heldout_TESLA_mmc7_validation | v2_multimodal_esm2_35m_lr_fast          | 310 |       4 |    0.0129032 | 0.0370173 |               0   |              0.05 | reviewer_safe_internal_locked |
| source_heldout_TESLA_mmc7_validation | v2_esm2_150m_lr_fast                    | 310 |       4 |    0.0129032 | 0.0385312 |               0   |              0    | reviewer_safe_internal_locked |

## Decision Impact

- ESM2-150M improves exact peptide-HLA top10 to 0.8 but does not beat the v1 fusion AUPRC leader.
- ESM2-650M improves near-peptide top10 to 0.7 but still does not clearly beat the v1 rule-gated fusion AUPRC.
- ESM2-650M gives useful HLA-supertype support but does not overturn the benchmark-only verdict.
- The public training corpus overlap blocker remains unresolved.

Workbook: `/home/seungho/personal/THCA_data_analysis/project/results/cross_neo_v2_sota_sprint_2026_05_09/CROSS_Neo_v2_all_results_summary_plus_runpod_esm2.xlsx`
