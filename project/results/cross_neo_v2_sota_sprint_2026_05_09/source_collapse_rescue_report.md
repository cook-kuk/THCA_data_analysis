# Source Collapse Rescue Report

## Best Source-Heldout Rows

| source                | model_name                        |   n |   n_pos |   prevalence |     AUPRC |   top10_precision |   top20_precision | claim_status                  |
|:----------------------|:----------------------------------|----:|--------:|-------------:|----------:|------------------:|------------------:|:------------------------------|
| CEDAR                 | v2_plm_lr_frozen_hash_pilot       | 909 |     851 |    0.936194  | 0.957745  |               1   |              1    | diagnostic                    |
| CEDAR                 | v2_multimodal_lr                  | 909 |     851 |    0.936194  | 0.946896  |               1   |              1    | reviewer_safe_internal_locked |
| CEDAR                 | v2_cf_plm_lr                      | 909 |     851 |    0.936194  | 0.946877  |               1   |              1    | reviewer_safe_internal_locked |
| NEPdb                 | v2_groupdro_proxy_cf_lr           | 572 |     151 |    0.263986  | 0.281349  |               0.5 |              0.45 | reviewer_safe_internal_locked |
| NEPdb                 | source_qk_compact_gamma1          | 886 |     354 |    0.399549  | 0.406426  |               0.4 |              0.4  | diagnostic                    |
| NEPdb                 | v2_multimodal_lr                  | 572 |     151 |    0.263986  | 0.30528   |               0.3 |              0.4  | reviewer_safe_internal_locked |
| TESLA_mmc4            | v2_groupdro_proxy_cf_lr           | 605 |      37 |    0.061157  | 0.077088  |               0.1 |              0.15 | reviewer_safe_internal_locked |
| TESLA_mmc4            | v2_class_balanced_cf_plm_hgb      | 605 |      37 |    0.061157  | 0.0841894 |               0.1 |              0.05 | reviewer_safe_internal_locked |
| TESLA_mmc4            | v2_cf_plm_rf_hla_ranknorm         | 605 |      37 |    0.061157  | 0.0792775 |               0.1 |              0.05 | reviewer_safe_internal_locked |
| TESLA_mmc7_validation | anchor_rf                         | 319 |       6 |    0.0188088 | 0.0378096 |               0   |              0    | reviewer_safe_internal_locked |
| TESLA_mmc7_validation | v2_counterfactual_lr_hla_ranknorm | 310 |       4 |    0.0129032 | 0.0332845 |               0   |              0    | reviewer_safe_internal_locked |
| TESLA_mmc7_validation | v2_cf_plm_rf_hla_ranknorm         | 310 |       4 |    0.0129032 | 0.0277166 |               0   |              0    | reviewer_safe_internal_locked |

- NEPdb nonzero top10 models: 8
- TESLA nonzero top20 models: 10

## Interpretation

- This remains an internal/source-heldout stress test, not external validation.
- Imported v1 diagnostic rows can have different source-heldout cardinalities from v2 split files; compare source rows descriptively unless the n/test definition matches.
- Nonzero top-k recovery is necessary but not sufficient for SOTA; public overlap remains a claim blocker unless resolved.
- If TESLA positives remain in the low-score tail across model families, the result should be framed as source-shift benchmark evidence rather than hidden.
