# CROSS-all Failure Diagnosis

Primary observation: naive concatenation does not dominate counterfactual features, consistent with small-n noisy-branch overfitting.

## Drop-Group Ablation Top Rows
| split_name                 | feature_set        |   n |   n_pos |   prevalence |    AUPRC |    AUROC |    Brier |   top5_precision |   recall_at_5 |   enrichment_at_5 |   top10_precision |   recall_at_10 |   enrichment_at_10 |
|:---------------------------|:-------------------|----:|--------:|-------------:|---------:|---------:|---------:|-----------------:|--------------:|------------------:|------------------:|---------------:|-------------------:|
| study_heldout              | C + QK             |   9 |       5 |     0.555556 | 0.811111 | 0.6      | 0.293852 |              0.6 |      0.6      |           1.08    |          0.555556 |       1        |            1       |
| study_heldout              | C + retrieval      |   9 |       5 |     0.555556 | 0.711111 | 0.45     | 0.312755 |              0.4 |      0.4      |           0.72    |          0.555556 |       1        |            1       |
| study_heldout              | C + retrieval + QK |   9 |       5 |     0.555556 | 0.645397 | 0.45     | 0.292341 |              0.6 |      0.6      |           1.08    |          0.555556 |       1        |            1       |
| study_heldout              | C only             |   9 |       5 |     0.555556 | 0.611111 | 0.35     | 0.315256 |              0.4 |      0.4      |           0.72    |          0.555556 |       1        |            1       |
| study_heldout              | C + structure + QK |   9 |       5 |     0.555556 | 0.605397 | 0.35     | 0.328645 |              0.4 |      0.4      |           0.72    |          0.555556 |       1        |            1       |
| study_heldout              | all                |   9 |       5 |     0.555556 | 0.605397 | 0.35     | 0.329384 |              0.4 |      0.4      |           0.72    |          0.555556 |       1        |            1       |
| study_heldout              | C + structure      |   9 |       5 |     0.555556 | 0.605397 | 0.35     | 0.342509 |              0.4 |      0.4      |           0.72    |          0.555556 |       1        |            1       |
| exact_peptide_hla_holdout  | C + retrieval      |  89 |      21 |     0.235955 | 0.563002 | 0.757003 | 0.173105 |              0.8 |      0.190476 |           3.39048 |          0.6      |       0.285714 |            2.54286 |
| exact_peptide_hla_holdout  | C + QK             |  89 |      21 |     0.235955 | 0.529267 | 0.718487 | 0.173208 |              0.8 |      0.190476 |           3.39048 |          0.7      |       0.333333 |            2.96667 |
| hla_stratified_group_5fold | C only             |  89 |      21 |     0.235955 | 0.519461 | 0.702381 | 0.17952  |              0.6 |      0.142857 |           2.54286 |          0.6      |       0.285714 |            2.54286 |
| hla_supertype_heldout      | C + retrieval + QK |  72 |      15 |     0.208333 | 0.506766 | 0.769591 | 0.167758 |              0.6 |      0.2      |           2.88    |          0.5      |       0.333333 |            2.4     |
| exact_peptide_hla_holdout  | C + structure      |  89 |      21 |     0.235955 | 0.506174 | 0.692577 | 0.17925  |              0.8 |      0.190476 |           3.39048 |          0.6      |       0.285714 |            2.54286 |

## Inner-Fold Permutation Importance
| feature_group      |   count |       mean |      median |
|:-------------------|--------:|-----------:|------------:|
| C_counterfactual   |      40 | 0.0993132  |  0.115006   |
| QK_fixed           |      40 | 0.00952539 |  0.00416667 |
| structure_geometry |      40 | 0.00174225 | -0.00580808 |

Interpretation: use late/gated fusion rather than raw feature concatenation.
