# CROSS-Neo v1 Gated MoE Report

All outer-fold scores are trained from the outer train fold only. Nested learned weights are selected from inner OOF predictions within the outer train fold.

## Top Gated Models
| split_name                       | model                                     |   n |   n_pos |   prevalence |    AUPRC |    AUROC |    Brier |   top5_precision |   recall_at_5 |   enrichment_at_5 |   top10_precision |   recall_at_10 |   enrichment_at_10 |
|:---------------------------------|:------------------------------------------|----:|--------:|-------------:|---------:|---------:|---------:|-----------------:|--------------:|------------------:|------------------:|---------------:|-------------------:|
| exact_peptide_hla_holdout        | rule_gated_C_QK_structure                 |  89 |      21 |     0.235955 | 0.647073 | 0.778711 | 0.182511 |              1   |      0.238095 |           4.2381  |          0.7      |      0.333333  |            2.96667 |
| study_heldout                    | stacked_calibrated_meta_model_lowdim_only |   9 |       5 |     0.555556 | 0.625397 | 0.4      | 0.283042 |              0.4 |      0.4      |           0.72    |          0.555556 |      1         |            1       |
| hla_supertype_heldout            | nested_learned_gate_C_QK_structure        |  72 |      15 |     0.208333 | 0.611634 | 0.816374 | 0.173832 |              0.8 |      0.266667 |           3.84    |          0.5      |      0.333333  |            2.4     |
| exact_peptide_hla_holdout        | nested_learned_gate_C_QK_structure        |  89 |      21 |     0.235955 | 0.610034 | 0.77451  | 0.182866 |              0.8 |      0.190476 |           3.39048 |          0.8      |      0.380952  |            3.39048 |
| exact_peptide_hla_holdout        | prespecified_equal_weight_C_QK_no_anchor  |  89 |      21 |     0.235955 | 0.594992 | 0.766106 | 0.187094 |              0.8 |      0.190476 |           3.39048 |          0.6      |      0.285714  |            2.54286 |
| study_heldout                    | prespecified_C_0.5_QK_quantum_0.5         |   9 |       5 |     0.555556 | 0.576825 | 0.25     | 0.31791  |              0.4 |      0.4      |           0.72    |          0.555556 |      1         |            1       |
| study_heldout                    | nested_learned_gate_C_QK_structure        |   9 |       5 |     0.555556 | 0.563492 | 0.2      | 0.283    |              0.2 |      0.2      |           0.36    |          0.555556 |      1         |            1       |
| study_heldout                    | prespecified_equal_weight_C_QK_no_anchor  |   9 |       5 |     0.555556 | 0.563492 | 0.2      | 0.283    |              0.2 |      0.2      |           0.36    |          0.555556 |      1         |            1       |
| study_heldout                    | rule_gated_C_QK_structure                 |   9 |       5 |     0.555556 | 0.563492 | 0.2      | 0.297193 |              0.2 |      0.2      |           0.36    |          0.555556 |      1         |            1       |
| repeated_stratified_5x5_internal | rule_gated_C_QK_structure                 | 445 |     105 |     0.235955 | 0.561522 | 0.748964 | 0.184974 |              1   |      0.047619 |           4.2381  |          1        |      0.0952381 |            4.2381  |
| repeated_stratified_5x5_internal | prespecified_equal_weight_C_QK_no_anchor  | 445 |     105 |     0.235955 | 0.546597 | 0.740448 | 0.188814 |              1   |      0.047619 |           4.2381  |          1        |      0.0952381 |            4.2381  |
| hla_supertype_heldout            | rule_gated_C_QK_structure                 |  72 |      15 |     0.208333 | 0.536984 | 0.746199 | 0.185919 |              0.8 |      0.266667 |           3.84    |          0.5      |      0.333333  |            2.4     |
| hla_stratified_group_5fold       | prespecified_equal_weight_C_QK_no_anchor  |  89 |      21 |     0.235955 | 0.533398 | 0.752101 | 0.191442 |              0.8 |      0.190476 |           3.39048 |          0.5      |      0.238095  |            2.11905 |
| hla_stratified_group_5fold       | rule_gated_C_QK_structure                 |  89 |      21 |     0.235955 | 0.530409 | 0.762605 | 0.188401 |              0.6 |      0.142857 |           2.54286 |          0.5      |      0.238095  |            2.11905 |
| near_peptide_cluster_holdout     | rule_gated_C_QK_structure                 |  89 |      21 |     0.235955 | 0.51374  | 0.717087 | 0.186081 |              0.8 |      0.190476 |           3.39048 |          0.6      |      0.285714  |            2.54286 |

Prespecified and rule-gated variants are eligible for conservative interpretation; nested/stacked variants are fold-safe but still small-n.
