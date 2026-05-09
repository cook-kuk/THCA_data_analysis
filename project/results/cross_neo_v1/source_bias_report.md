# CROSS-Neo v1 Source Bias Report

This is a source-heldout stress analysis, not external validation.

## Source Diagnostics
| source                |   n |   n_pos |   prevalence |   median_peptide_len |   hla_nunique |   hla_top_fraction |   wt_available_rate |   source_protein_available_rate |   mean_counterfactual_feature |   mean_qk_feature |
|:----------------------|----:|--------:|-------------:|---------------------:|--------------:|-------------------:|--------------------:|--------------------------------:|------------------------------:|------------------:|
| CEDAR                 | 909 |     851 |    0.936194  |                    9 |            95 |           0.152915 |                   0 |                               0 |                      0.277109 |        nan        |
| NEPdb                 | 572 |     151 |    0.263986  |                   10 |            53 |           0.256993 |                   0 |                               0 |                      0.296739 |          0.423916 |
| TESLA_mmc4            | 605 |      37 |    0.061157  |                    9 |            13 |           0.393388 |                   0 |                               0 |                      0.284486 |        nan        |
| TESLA_mmc7_validation | 310 |       4 |    0.0129032 |                   10 |             7 |           0.277419 |                   0 |                               0 |                      0.287777 |          0.684156 |
| ITSNdb strict         |  89 |      21 |    0.235955  |                    9 |            13 |           0.370787 |                   0 |                               0 |                    nan        |        nan        |

## Source Predictability
| task                                 |    n |   accuracy |   macro_f1 |
|:-------------------------------------|-----:|-----------:|-----------:|
| predict_source_from_allowed_features | 2396 |   0.521285 |   0.359127 |

## Corrected Source-Heldout Metrics
| heldout_study         | model                                             |   n |   n_pos |   prevalence |     AUPRC |    AUROC |    Brier |   top5_precision |   recall_at_5 |   enrichment_at_5 |   top10_precision |   recall_at_10 |   enrichment_at_10 |
|:----------------------|:--------------------------------------------------|----:|--------:|-------------:|----------:|---------:|---------:|-----------------:|--------------:|------------------:|------------------:|---------------:|-------------------:|
| CEDAR                 | source_balanced_plus_pu_rf                        | 909 |     851 |    0.936194  | 0.937603  | 0.483792 | 0.107227 |              1   |    0.00587544 |           1.06816 |               1   |      0.0117509 |            1.06816 |
| CEDAR                 | source_balanced_plus_pu_rf_train_prior_calibrated | 909 |     851 |    0.936194  | 0.937603  | 0.483792 | 0.383033 |              1   |    0.00587544 |           1.06816 |               1   |      0.0117509 |            1.06816 |
| CEDAR                 | pu_weighted_rf                                    | 909 |     851 |    0.936194  | 0.934897  | 0.468111 | 0.100235 |              1   |    0.00587544 |           1.06816 |               1   |      0.0117509 |            1.06816 |
| CEDAR                 | pu_weighted_rf_train_prior_calibrated             | 909 |     851 |    0.936194  | 0.934897  | 0.468111 | 0.34965  |              1   |    0.00587544 |           1.06816 |               1   |      0.0117509 |            1.06816 |
| CEDAR                 | source_balanced_rf                                | 909 |     851 |    0.936194  | 0.931176  | 0.433263 | 0.377184 |              1   |    0.00587544 |           1.06816 |               1   |      0.0117509 |            1.06816 |
| CEDAR                 | source_balanced_rf_train_prior_calibrated         | 909 |     851 |    0.936194  | 0.931176  | 0.433263 | 0.751465 |              1   |    0.00587544 |           1.06816 |               1   |      0.0117509 |            1.06816 |
| NEPdb                 | source_balanced_plus_pu_rf                        | 572 |     151 |    0.263986  | 0.225328  | 0.448349 | 0.412326 |              0   |    0          |           0       |               0   |      0         |            0       |
| NEPdb                 | source_balanced_plus_pu_rf_train_prior_calibrated | 572 |     151 |    0.263986  | 0.225328  | 0.448349 | 0.448148 |              0   |    0          |           0       |               0   |      0         |            0       |
| NEPdb                 | pu_weighted_rf                                    | 572 |     151 |    0.263986  | 0.220262  | 0.433625 | 0.497423 |              0   |    0          |           0       |               0   |      0         |            0       |
| NEPdb                 | pu_weighted_rf_train_prior_calibrated             | 572 |     151 |    0.263986  | 0.220262  | 0.433625 | 0.529746 |              0   |    0          |           0       |               0   |      0         |            0       |
| NEPdb                 | source_balanced_rf                                | 572 |     151 |    0.263986  | 0.207305  | 0.391295 | 0.247724 |              0   |    0          |           0       |               0   |      0         |            0       |
| NEPdb                 | source_balanced_rf_train_prior_calibrated         | 572 |     151 |    0.263986  | 0.207305  | 0.391295 | 0.262984 |              0   |    0          |           0       |               0   |      0         |            0       |
| TESLA_mmc7_validation | source_balanced_rf                                | 310 |       4 |    0.0129032 | 0.142665  | 0.693627 | 0.173151 |              0.2 |    0.25       |          15.5     |               0.1 |      0.25      |            7.75    |
| TESLA_mmc7_validation | source_balanced_rf_train_prior_calibrated         | 310 |       4 |    0.0129032 | 0.142665  | 0.693627 | 0.226145 |              0.2 |    0.25       |          15.5     |               0.1 |      0.25      |            7.75    |
| TESLA_mmc4            | pu_weighted_rf                                    | 605 |      37 |    0.061157  | 0.0775619 | 0.609678 | 0.643127 |              0   |    0          |           0       |               0   |      0         |            0       |
| TESLA_mmc4            | pu_weighted_rf_train_prior_calibrated             | 605 |      37 |    0.061157  | 0.0775619 | 0.609678 | 0.739354 |              0   |    0          |           0       |               0   |      0         |            0       |
| TESLA_mmc4            | source_balanced_plus_pu_rf_train_prior_calibrated | 605 |      37 |    0.061157  | 0.0730341 | 0.58641  | 0.637203 |              0   |    0          |           0       |               0   |      0         |            0       |
| TESLA_mmc4            | source_balanced_plus_pu_rf                        | 605 |      37 |    0.061157  | 0.0730341 | 0.58641  | 0.515843 |              0   |    0          |           0       |               0   |      0         |            0       |

## v0 Source-Heldout Reference
| heldout_study         | method                                   |   n |   n_pos |   prevalence |     AUPRC |    AUROC |    Brier |   top5_precision |   recall_at_5 |   enrichment_at_5 |   top10_precision |   recall_at_10 |   enrichment_at_10 |
|:----------------------|:-----------------------------------------|----:|--------:|-------------:|----------:|---------:|---------:|-----------------:|--------------:|------------------:|------------------:|---------------:|-------------------:|
| CEDAR                 | sourceheld_counterfactual_rf             | 913 |     851 |    0.932092  | 0.918082  | 0.416341 | 0.3709   |              0.8 |    0.00470035 |          0.858284 |               0.9 |      0.0105758 |            0.96557 |
| NEPdb                 | sourceheld_counterfactual_rf             | 886 |     354 |    0.399549  | 0.40534   | 0.49521  | 0.261428 |              0   |    0          |          0        |               0   |      0         |            0       |
| TESLA_mmc4            | sourceheld_counterfactual_rf             | 610 |      37 |    0.0606557 | 0.0650109 | 0.549267 | 0.24215  |              0   |    0          |          0        |               0   |      0         |            0       |
| TESLA_mmc7_validation | sourceheld_counterfactual_rf             | 319 |       6 |    0.0188088 | 0.0378096 | 0.655485 | 0.212548 |              0   |    0          |          0        |               0   |      0         |            0       |
| CEDAR                 | sourceheld_prespecified_late_fusion_w0.5 | 913 |     851 |    0.932092  | 0.921816  | 0.439028 | 0.413353 |              1   |    0.00587544 |          1.07286  |               1   |      0.0117509 |            1.07286 |
| NEPdb                 | sourceheld_prespecified_late_fusion_w0.5 | 886 |     354 |    0.399549  | 0.395509  | 0.479297 | 0.293726 |              0   |    0          |          0        |               0   |      0         |            0       |
| TESLA_mmc4            | sourceheld_prespecified_late_fusion_w0.5 | 610 |      37 |    0.0606557 | 0.0551107 | 0.464931 | 0.31171  |              0   |    0          |          0        |               0   |      0         |            0       |
| TESLA_mmc7_validation | sourceheld_prespecified_late_fusion_w0.5 | 319 |       6 |    0.0188088 | 0.0168969 | 0.294462 | 0.241671 |              0   |    0          |          0        |               0   |      0         |            0       |
| CEDAR                 | sourceheld_qk_compact_gamma1             | 913 |     851 |    0.932092  | 0.926744  | 0.461317 | 0.527638 |              0.8 |    0.00470035 |          0.858284 |               0.9 |      0.0105758 |            0.96557 |
| NEPdb                 | sourceheld_qk_compact_gamma1             | 886 |     354 |    0.399549  | 0.406426  | 0.482637 | 0.401408 |              0.4 |    0.00564972 |          1.00113  |               0.4 |      0.0112994 |            1.00113 |
| TESLA_mmc4            | sourceheld_qk_compact_gamma1             | 610 |      37 |    0.0606557 | 0.0516625 | 0.437904 | 0.476714 |              0   |    0          |          0        |               0   |      0         |            0       |
| TESLA_mmc7_validation | sourceheld_qk_compact_gamma1             | 319 |       6 |    0.0188088 | 0.0145965 | 0.180511 | 0.350284 |              0   |    0          |          0        |               0   |      0         |            0       |

Interpretation: source labels are highly imbalanced and source identity is partly predictable from allowed features. Any source-heldout top-k collapse should be treated as assay/source shift unless corrected models show stable top-k rescue.
