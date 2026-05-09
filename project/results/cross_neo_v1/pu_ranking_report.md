# CROSS-Neo v1 PU Ranking Report

Negatives are treated as unlabeled/ambiguous for weighted and bagged PU variants. No public predictor scores are used.

## Top PU Rows
| split_name                       | model                            |   n |   n_pos |   prevalence |    AUPRC |    AUROC |    Brier |   top5_precision |   recall_at_5 |   enrichment_at_5 |   top10_precision |   recall_at_10 |   enrichment_at_10 |
|:---------------------------------|:---------------------------------|----:|--------:|-------------:|---------:|---------:|---------:|-----------------:|--------------:|------------------:|------------------:|---------------:|-------------------:|
| study_heldout                    | positive_only_centroid           |   9 |       5 |     0.555556 | 0.82619  | 0.75     | 0.458469 |              0.6 |     0.6       |           1.08    |          0.555556 |      1         |            1       |
| study_heldout                    | pu_logistic_weighted             |   9 |       5 |     0.555556 | 0.761111 | 0.55     | 0.335175 |              0.6 |     0.6       |           1.08    |          0.555556 |      1         |            1       |
| study_heldout                    | pairwise_positive_over_unlabeled |   9 |       5 |     0.555556 | 0.694444 | 0.5      | 0.331432 |              0.6 |     0.6       |           1.08    |          0.555556 |      1         |            1       |
| study_heldout                    | bagging_pu_rf                    |   9 |       5 |     0.555556 | 0.605397 | 0.35     | 0.283942 |              0.4 |     0.4       |           0.72    |          0.555556 |      1         |            1       |
| exact_peptide_hla_holdout        | pu_logistic_weighted             |  89 |      21 |     0.235955 | 0.545008 | 0.698179 | 0.195375 |              0.8 |     0.190476  |           3.39048 |          0.7      |      0.333333  |            2.96667 |
| hla_supertype_heldout            | bagging_pu_rf                    |  72 |      15 |     0.208333 | 0.516225 | 0.785965 | 0.182186 |              0.8 |     0.266667  |           3.84    |          0.5      |      0.333333  |            2.4     |
| hla_stratified_group_5fold       | bagging_pu_rf                    |  89 |      21 |     0.235955 | 0.515924 | 0.741597 | 0.189625 |              0.8 |     0.190476  |           3.39048 |          0.6      |      0.285714  |            2.54286 |
| hla_stratified_group_5fold       | pu_logistic_weighted             |  89 |      21 |     0.235955 | 0.50624  | 0.69888  | 0.191647 |              0.6 |     0.142857  |           2.54286 |          0.7      |      0.333333  |            2.96667 |
| exact_peptide_hla_holdout        | bagging_pu_rf                    |  89 |      21 |     0.235955 | 0.503928 | 0.718487 | 0.189546 |              0.8 |     0.190476  |           3.39048 |          0.6      |      0.285714  |            2.54286 |
| repeated_stratified_5x5_internal | pu_logistic_weighted             | 445 |     105 |     0.235955 | 0.501967 | 0.674062 | 0.218355 |              1   |     0.047619  |           4.2381  |          1        |      0.0952381 |            4.2381  |
| near_peptide_cluster_holdout     | pu_logistic_weighted             |  89 |      21 |     0.235955 | 0.501837 | 0.679272 | 0.216103 |              0.8 |     0.190476  |           3.39048 |          0.6      |      0.285714  |            2.54286 |
| hla_supertype_heldout            | pu_logistic_weighted             |  72 |      15 |     0.208333 | 0.500442 | 0.74269  | 0.182448 |              0.8 |     0.266667  |           3.84    |          0.5      |      0.333333  |            2.4     |
| exact_peptide_hla_holdout        | pairwise_positive_over_unlabeled |  89 |      21 |     0.235955 | 0.497385 | 0.686275 | 0.289911 |              0.8 |     0.190476  |           3.39048 |          0.6      |      0.285714  |            2.54286 |
| hla_stratified_group_5fold       | pairwise_positive_over_unlabeled |  89 |      21 |     0.235955 | 0.494193 | 0.689076 | 0.2711   |              0.8 |     0.190476  |           3.39048 |          0.5      |      0.238095  |            2.11905 |
| near_peptide_cluster_holdout     | pairwise_positive_over_unlabeled |  89 |      21 |     0.235955 | 0.491819 | 0.70098  | 0.29687  |              0.8 |     0.190476  |           3.39048 |          0.5      |      0.238095  |            2.11905 |
| hla_supertype_heldout            | pairwise_positive_over_unlabeled |  72 |      15 |     0.208333 | 0.482859 | 0.762573 | 0.242131 |              0.6 |     0.2       |           2.88    |          0.5      |      0.333333  |            2.4     |
| repeated_stratified_5x5_internal | pairwise_positive_over_unlabeled | 445 |     105 |     0.235955 | 0.443567 | 0.652213 | 0.310883 |              0.8 |     0.0380952 |           3.39048 |          0.9      |      0.0857143 |            3.81429 |
| near_peptide_cluster_holdout     | bagging_pu_rf                    |  89 |      21 |     0.235955 | 0.441012 | 0.642157 | 0.19795  |              0.6 |     0.142857  |           2.54286 |          0.4      |      0.190476  |            1.69524 |

PU models are promotion candidates only if they improve top-k precision, not merely AUROC.
