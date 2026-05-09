# Fast ESM2+QK Gate Report

| split_name                   | pre_gate_best                           |   pre_gate_AUPRC |   pre_gate_top10 |   gate_AUPRC |   gate_top10 | combined_best                           |   combined_AUPRC |   combined_top10 |
|:-----------------------------|:----------------------------------------|-----------------:|-----------------:|-------------:|-------------:|:----------------------------------------|-----------------:|-----------------:|
| exact_peptide_hla_holdout    | prespecified_rf_qk_no_anchor_w0.5       |         0.586965 |              0.7 |     0.636885 |          0.7 | fast_nested_esm2_qk_gate                |         0.636885 |              0.7 |
| near_peptide_cluster_holdout | rule_gate_rf_qk_fallback_train_selected |         0.518669 |              0.7 |     0.411047 |          0.4 | rule_gate_rf_qk_fallback_train_selected |         0.518669 |              0.7 |
| hla_stratified_group_5fold   | v2_cf_plm_lr                            |         0.563197 |              0.5 |     0.533327 |          0.7 | v2_cf_plm_lr                            |         0.563197 |              0.5 |
| hla_supertype_heldout        | v2_multimodal_lr                        |         0.58542  |              0.6 |     0.516795 |          0.5 | v2_multimodal_lr                        |         0.58542  |              0.6 |

Selection used same-split OOF predictions on outer-train rows only.

## Selected Experts

| split_name                   | fold_id   | base                                    | esm                             |    w |   train_AUPRC |   train_top10 |
|:-----------------------------|:----------|:----------------------------------------|:--------------------------------|-----:|--------------:|--------------:|
| exact_peptide_hla_holdout    | pmhc_00   | rule_gate_rf_qk_fallback_train_selected | v2_cf_esm2_650m_lr_fast         | 0.5  |      0.670501 |           0.7 |
| exact_peptide_hla_holdout    | pmhc_01   | rule_gate_rf_qk_fallback_train_selected | v2_multimodal_esm2_650m_lr_fast | 0.5  |      0.75032  |           0.8 |
| exact_peptide_hla_holdout    | pmhc_02   | prespecified_rf_qk_no_anchor_w0.5       | v2_multimodal_esm2_650m_lr_fast | 0.5  |      0.654771 |           0.7 |
| exact_peptide_hla_holdout    | pmhc_03   | rule_gate_rf_qk_fallback_train_selected | v2_multimodal_esm2_650m_lr_fast | 0.5  |      0.663275 |           0.7 |
| exact_peptide_hla_holdout    | pmhc_04   | rule_gate_rf_qk_fallback_train_selected | v2_cf_esm2_150m_lr_fast         | 0.5  |      0.675124 |           0.9 |
| near_peptide_cluster_holdout | near_00   | rule_gate_rf_qk_fallback_train_selected | v2_cf_esm2_650m_lr_fast         | 0.75 |      0.523257 |           0.5 |
| near_peptide_cluster_holdout | near_01   | rule_gate_rf_qk_fallback_train_selected | v2_cf_esm2_150m_lr_fast         | 1    |      0.558806 |           0.6 |
| near_peptide_cluster_holdout | near_02   | prespecified_rf_qk_no_anchor_w0.5       | v2_cf_esm2_650m_lr_fast         | 0.75 |      0.58948  |           0.6 |
| near_peptide_cluster_holdout | near_03   | rule_gate_rf_qk_fallback_train_selected | v2_cf_esm2_650m_lr_fast         | 0.75 |      0.594934 |           0.6 |
| near_peptide_cluster_holdout | near_04   | rule_gate_rf_qk_fallback_train_selected | v2_multimodal_esm2_650m_lr_fast | 0.5  |      0.559134 |           0.6 |
| hla_stratified_group_5fold   | hla_00    | prespecified_rf_qk_no_anchor_w0.5       | v2_multimodal_esm2_650m_lr_fast | 0.5  |      0.65674  |           0.7 |
| hla_stratified_group_5fold   | hla_01    | prespecified_rf_qk_no_anchor_w0.5       | v2_cf_esm2_150m_lr_fast         | 0.75 |      0.553    |           0.6 |
| hla_stratified_group_5fold   | hla_02    | prespecified_rf_qk_no_anchor_w0.5       | v2_multimodal_esm2_650m_lr_fast | 0.75 |      0.627424 |           0.6 |
| hla_stratified_group_5fold   | hla_03    | prespecified_rf_qk_no_anchor_w0.5       | v2_multimodal_esm2_650m_lr_fast | 0.75 |      0.656889 |           0.7 |
| hla_stratified_group_5fold   | hla_04    | prespecified_rf_qk_no_anchor_w0.5       | v2_multimodal_esm2_650m_lr_fast | 0.5  |      0.605035 |           0.7 |
| hla_supertype_heldout        | super_00  | rule_gate_rf_qk_fallback_train_selected | v2_multimodal_esm2_650m_lr_fast | 0.5  |      0.63631  |           0.7 |
| hla_supertype_heldout        | super_01  | v2_multimodal_lr                        | v2_multimodal_esm2_650m_lr_fast | 0.5  |      0.567778 |           0.3 |
| hla_supertype_heldout        | super_02  | v2_multimodal_lr                        | v2_cf_esm2_150m_lr_fast         | 1    |      0.65329  |           0.6 |
| hla_supertype_heldout        | super_09  | rule_gate_rf_qk_fallback_train_selected | v2_multimodal_esm2_650m_lr_fast | 0.5  |      0.678194 |           0.7 |

Workbook: `/home/seungho/personal/THCA_data_analysis/project/results/cross_neo_v2_sota_sprint_2026_05_09/CROSS_Neo_v2_all_results_summary_plus_runpod_esm2_fast_gate.xlsx`
