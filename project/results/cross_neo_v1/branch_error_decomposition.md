# CROSS-Neo v0 Branch Error Decomposition

Branch ranks are computed inside each v0 split/fold. Late-fusion rows are reconstructed from fold-safe v0 predictions.

QK rescue cases: 28
QK harm cases: 73
Structure harm cases: 67

## Rank Correlation Snapshot
| split_name                       | expert_a                | expert_b                                                    |   n |   spearman_rank_corr |
|:---------------------------------|:------------------------|:------------------------------------------------------------|----:|---------------------:|
| exact_peptide_hla_holdout        | C_counterfactual_rf     | D_structure_geometry_rf                                     |  89 |           -0.161224  |
| hla_supertype_heldout            | D_structure_geometry_rf | qk_quantum_only_gamma1                                      |  72 |           -0.115971  |
| hla_supertype_heldout            | D_structure_geometry_rf | late_prespecified_equal_weight_Cw0.5_qk_quantum_only_gamma1 |  72 |           -0.0883973 |
| exact_peptide_hla_holdout        | D_structure_geometry_rf | late_prespecified_equal_weight_Cw0.5_qk_no_anchor_gamma1    |  89 |           -0.0774576 |
| near_peptide_cluster_holdout     | D_structure_geometry_rf | late_prespecified_equal_weight_Cw0.5_qk_quantum_only_gamma1 |  89 |           -0.0729304 |
| exact_peptide_hla_holdout        | D_structure_geometry_rf | late_prespecified_equal_weight_Cw0.5_qk_quantum_only_gamma1 |  89 |           -0.0689317 |
| near_peptide_cluster_holdout     | D_structure_geometry_rf | E_quantum_fixed_rf                                          |  89 |           -0.0597137 |
| near_peptide_cluster_holdout     | C_counterfactual_rf     | D_structure_geometry_rf                                     |  89 |           -0.032369  |
| near_peptide_cluster_holdout     | D_structure_geometry_rf | qk_quantum_only_gamma1                                      |  89 |           -0.0248282 |
| repeated_stratified_5x5_internal | D_structure_geometry_rf | late_prespecified_equal_weight_Cw0.5_qk_quantum_only_gamma1 | 445 |           -0.022987  |
| repeated_stratified_5x5_internal | D_structure_geometry_rf | qk_quantum_only_gamma1                                      | 445 |           -0.0212737 |
| hla_supertype_heldout            | C_counterfactual_rf     | qk_quantum_only_gamma1                                      |  72 |           -0.0168784 |

## Oracle Upper Bound Snapshot
| split_name                       | strategy                                                    |   label_informed_oracle |   n |   n_pos |   prevalence |    AUPRC |    AUROC |     Brier |   top5_precision |   recall_at_5 |   enrichment_at_5 |   top10_precision |   recall_at_10 |   enrichment_at_10 |
|:---------------------------------|:------------------------------------------------------------|------------------------:|----:|--------:|-------------:|---------:|---------:|----------:|-----------------:|--------------:|------------------:|------------------:|---------------:|-------------------:|
| near_peptide_cluster_holdout     | oracle_label_informed_best_branch                           |                       1 |  89 |      21 |     0.235955 | 1        | 1        | 0.0788941 |              1   |      0.238095 |            4.2381 |          1        |      0.47619   |             4.2381 |
| exact_peptide_hla_holdout        | oracle_label_informed_best_branch                           |                       1 |  89 |      21 |     0.235955 | 0.997835 | 0.9993   | 0.0722659 |              1   |      0.238095 |            4.2381 |          1        |      0.47619   |             4.2381 |
| hla_supertype_heldout            | oracle_label_informed_best_branch                           |                       1 |  72 |      15 |     0.208333 | 0.995833 | 0.99883  | 0.0724658 |              1   |      0.333333 |            4.8    |          1        |      0.666667  |             4.8    |
| hla_stratified_group_5fold       | oracle_label_informed_best_branch                           |                       1 |  89 |      21 |     0.235955 | 0.993592 | 0.997899 | 0.076083  |              1   |      0.238095 |            4.2381 |          1        |      0.47619   |             4.2381 |
| repeated_stratified_5x5_internal | oracle_label_informed_best_branch                           |                       1 | 445 |     105 |     0.235955 | 0.988053 | 0.995742 | 0.0774897 |              1   |      0.047619 |            4.2381 |          1        |      0.0952381 |             4.2381 |
| study_heldout                    | C_counterfactual_rf                                         |                       0 |   9 |       5 |     0.555556 | 0.711111 | 0.45     | 0.311723  |              0.4 |      0.4      |            0.72   |          0.555556 |      1         |             1      |
| study_heldout                    | F_CROSS_all_rf                                              |                       0 |   9 |       5 |     0.555556 | 0.619286 | 0.4      | 0.331862  |              0.4 |      0.4      |            0.72   |          0.555556 |      1         |             1      |
| study_heldout                    | late_prespecified_equal_weight_Cw0.5_qk_quantum_only_gamma1 |                       0 |   9 |       5 |     0.555556 | 0.596825 | 0.3      | 0.31676   |              0.4 |      0.4      |            0.72   |          0.555556 |      1         |             1      |

Interpretation: oracle rows are label-informed upper bounds for complementarity diagnosis only, not deployable models.
