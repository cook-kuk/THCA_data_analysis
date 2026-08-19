# BioDarwin public-router v3

## What changed
- v2 used a single GA/RL router on the full public union.
- v3 does leave-one-bundle-out GA/RL, then bags the fold champions, then blends bagged vs full-data router.
- This is the more generalizable public router scaffold.

## Summary
| generated_at        |   elapsed_s |   rows | binary_bundles                                                   |   holdout_folds |   best_alpha_bagged |   alpha_fitness |   v2_mean_soft_ap |   v3_mean_full_ap |   v3_mean_bagged_ap |   v3_mean_blend_ap |   v3_mean_blend_auc | outputs                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |
|:--------------------|------------:|-------:|:-----------------------------------------------------------------|----------------:|--------------------:|----------------:|------------------:|------------------:|--------------------:|-------------------:|--------------------:|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| 2026-05-11T07:39:04 |      190.94 |   3804 | ['cedar_partial', 'itsndb', 'nepdb', 'tesla_mmc4', 'tesla_mmc7'] |               5 |                 0.2 |        0.905301 |          0.918611 |          0.918611 |            0.914185 |           0.920498 |            0.956964 | {'rows': '/home/seungho/personal/THCA_data_analysis/project/results/biodarwin_external_mega_comparison_2026_05_11/public_router_v3_loso_ga_rl/biodarwin_public_router_v3_rows.tsv', 'bundle_summary': '/home/seungho/personal/THCA_data_analysis/project/results/biodarwin_external_mega_comparison_2026_05_11/public_router_v3_loso_ga_rl/biodarwin_public_router_v3_bundle_summary.tsv', 'fold_summary': '/home/seungho/personal/THCA_data_analysis/project/results/biodarwin_external_mega_comparison_2026_05_11/public_router_v3_loso_ga_rl/biodarwin_public_router_v3_fold_summary.tsv', 'alpha_grid': '/home/seungho/personal/THCA_data_analysis/project/results/biodarwin_external_mega_comparison_2026_05_11/public_router_v3_loso_ga_rl/biodarwin_public_router_v3_alpha_grid.tsv', 'v1_eval': '/home/seungho/personal/THCA_data_analysis/project/results/biodarwin_external_mega_comparison_2026_05_11/public_router_v3_loso_ga_rl/biodarwin_public_router_v3_v1_eval.tsv', 'trace': '/home/seungho/personal/THCA_data_analysis/project/results/biodarwin_external_mega_comparison_2026_05_11/public_router_v3_loso_ga_rl/biodarwin_public_router_v3_fold_trace.tsv', 'report': '/home/seungho/personal/THCA_data_analysis/project/results/biodarwin_external_mega_comparison_2026_05_11/public_router_v3_loso_ga_rl/PUBLIC_ROUTER_REPORT_V3.md'} |

## Fold champions
| holdout_bundle   | train_bundles                                                         |   train_fitness |   train_mean_ap |   train_mean_auc |   train_min_ap |   train_std_ap |   holdout_n |   holdout_n_pos |   holdout_n_neg |   holdout_AUPRC |   holdout_AUROC |   holdout_top10_precision |   train_AUPRC |   train_AUROC |
|:-----------------|:----------------------------------------------------------------------|----------------:|----------------:|-----------------:|---------------:|---------------:|------------:|----------------:|----------------:|----------------:|----------------:|--------------------------:|--------------:|--------------:|
| cedar_partial    | dbpepneo2_mhci,itsndb,mcpas,neodb,nepdb,tesla_mmc4,tesla_mmc7         |        0.880411 |        0.90314  |         0.957011 |       0.794866 |      0.0805864 |         909 |             851 |              58 |        0.996778 |        0.955367 |                       1   |      0.858038 |      0.845027 |
| itsndb           | cedar_partial,dbpepneo2_mhci,mcpas,neodb,nepdb,tesla_mmc4,tesla_mmc7  |        0.9277   |        0.957548 |         0.981259 |       0.883374 |      0.0471945 |         319 |             136 |             183 |        0.725229 |        0.807015 |                       0.8 |      0.908144 |      0.854956 |
| nepdb            | cedar_partial,dbpepneo2_mhci,itsndb,mcpas,neodb,tesla_mmc4,tesla_mmc7 |        0.88629  |        0.91503  |         0.950525 |       0.789638 |      0.0885577 |         560 |             140 |             420 |        0.954472 |        0.983486 |                       1   |      0.857698 |      0.740961 |
| tesla_mmc4       | cedar_partial,dbpepneo2_mhci,itsndb,mcpas,neodb,nepdb,tesla_mmc7      |        0.904679 |        0.93867  |         0.951755 |       0.798328 |      0.0825988 |         605 |              37 |             568 |        0.618808 |        0.960839 |                       0.8 |      0.974847 |      0.947587 |
| tesla_mmc7       | cedar_partial,dbpepneo2_mhci,itsndb,mcpas,neodb,nepdb,tesla_mmc4      |        0.915354 |        0.923285 |         0.948897 |       0.793511 |      0.0774623 |         310 |               4 |             306 |        0.486607 |        0.98366  |                       0.3 |      0.931464 |      0.872837 |

## Binary bundle summary
| bundle         |   n |   n_pos |   n_neg |   full_AUPRC |   full_AUROC |   bagged_AUPRC |   bagged_AUROC |   blend_AUPRC |   blend_AUROC |   blend_top10_precision | dominant_bagged_hard_expert     |   dominant_bagged_hard_share |
|:---------------|----:|--------:|--------:|-------------:|-------------:|---------------:|---------------:|--------------:|--------------:|------------------------:|:--------------------------------|-----------------------------:|
| cedar_partial  | 909 |     851 |      58 |     0.996942 |     0.958001 |       0.996964 |       0.958122 |      0.996947 |      0.958041 |                     1   | RF_biophys                      |                     1        |
| dbpepneo2_mhci | 467 |     467 |       0 |   nan        |   nan        |     nan        |     nan        |    nan        |    nan        |                     1   | RF_biophys                      |                     1        |
| itsndb         | 319 |     136 |     183 |     0.78318  |     0.853182 |       0.791652 |       0.857642 |      0.78494  |      0.853785 |                     1   | BioDarwin_public_anchor_noEL_v4 |                     0.739812 |
| mcpas          |  82 |      82 |       0 |   nan        |   nan        |     nan        |     nan        |    nan        |    nan        |                     1   | RF_biophys                      |                     0.585366 |
| neodb          | 552 |     552 |       0 |   nan        |   nan        |     nan        |     nan        |    nan        |    nan        |                     1   | RF_biophys                      |                     0.922101 |
| nepdb          | 560 |     140 |     420 |     0.957747 |     0.984745 |       0.95849  |       0.984932 |      0.95835  |      0.984983 |                     1   | RF_biophys                      |                     1        |
| tesla_mmc4     | 605 |      37 |     568 |     0.855186 |     0.9872   |       0.873817 |       0.989627 |      0.862251 |      0.988009 |                     1   | BioDarwin_public_anchor_noEL_v4 |                     0.527273 |
| tesla_mmc7     | 310 |       4 |     306 |     1        |     1        |       0.95     |       0.999183 |      1        |      1        |                     0.4 | RF_biophys                      |                     0.667742 |

## v2 comparison
| bundle        |   soft_AUPRC |   soft_AUROC |   hard_AUPRC |   hard_AUROC |
|:--------------|-------------:|-------------:|-------------:|-------------:|
| cedar_partial |     0.996942 |     0.958001 |     0.996746 |     0.955529 |
| itsndb        |     0.78318  |     0.853182 |     0.769488 |     0.817944 |
| nepdb         |     0.957747 |     0.984745 |     0.918404 |     0.970952 |
| tesla_mmc4    |     0.855186 |     0.9872   |     0.236755 |     0.868553 |
| tesla_mmc7    |     1        |     1        |     0.101361 |     0.930556 |

## v1 comparison
| bundle        | chosen_algorithm   | router_class    |   n |   n_pos |   n_neg |   positive_rate | binary_metric_possible   | low_power   |    AUPRC |    AUROC |   top5_precision |   top5_hits |   top10_precision |   top10_hits |   top20_precision |   top20_hits |   top50_precision |   top50_hits |
|:--------------|:-------------------|:----------------|----:|--------:|--------:|----------------:|:-------------------------|:------------|---------:|---------:|-----------------:|------------:|------------------:|-------------:|------------------:|-------------:|------------------:|-------------:|
| cedar_partial | RF_biophys         | performance-max | 909 |     851 |      58 |       0.936194  | True                     | False       | 0.996746 | 0.955529 |              1   |           5 |               1   |           10 |              1    |           20 |              1    |           50 |
| nepdb         | RF_biophys         | performance-max | 560 |     140 |     420 |       0.25      | True                     | False       | 0.918404 | 0.970952 |              1   |           5 |               1   |           10 |              1    |           20 |              1    |           50 |
| itsndb        | ESM2_Bayesian      | performance-max | 311 |     130 |     181 |       0.418006  | True                     | False       | 0.691231 | 0.784063 |              1   |           5 |               0.7 |            7 |              0.75 |           15 |              0.74 |           37 |
| tesla_mmc4    | RF_biophys         | performance-max | 605 |      37 |     568 |       0.061157  | True                     | False       | 0.364116 | 0.883946 |              0.4 |           2 |               0.5 |            5 |              0.45 |            9 |              0.38 |           19 |
| tesla_mmc7    | BigMHC_IM          | performance-max | 310 |       4 |     306 |       0.0129032 | True                     | True        | 0.570076 | 0.99183  |              0.6 |           3 |               0.3 |            3 |              0.2  |            4 |              0.08 |            4 |
| cedar_partial | MHCflurry          | reviewer-safe   | 845 |     787 |      58 |       0.931361  | True                     | False       | 0.972759 | 0.737173 |              1   |           5 |               1   |           10 |              1    |           20 |              1    |           50 |
| nepdb         | PRIME              | reviewer-safe   | 483 |     140 |     343 |       0.289855  | True                     | False       | 0.511098 | 0.71893  |              0.8 |           4 |               0.8 |            8 |              0.7  |           14 |              0.54 |           27 |
| itsndb        | DeepImmuno         | reviewer-safe   | 319 |     136 |     183 |       0.426332  | True                     | False       | 0.610914 | 0.701905 |              0.8 |           4 |               0.8 |            8 |              0.7  |           14 |              0.58 |           29 |
| tesla_mmc4    | PRIME              | reviewer-safe   | 605 |      37 |     568 |       0.061157  | True                     | False       | 0.193735 | 0.805434 |              0.2 |           1 |               0.2 |            2 |              0.2  |            4 |              0.16 |            8 |

## Alpha grid
|   alpha_bagged |   fitness |
|---------------:|----------:|
|            0   |  0.90393  |
|            0.1 |  0.90475  |
|            0.2 |  0.905301 |
|            0.3 |  0.900461 |
|            0.4 |  0.900913 |
|            0.5 |  0.901611 |
|            0.6 |  0.901922 |
|            0.7 |  0.902527 |
|            0.8 |  0.902949 |
|            0.9 |  0.903116 |
|            1   |  0.903227 |

## Notes
- The reported champion is the bagged/full blend with the best public-bundle fitness.
- The hard router is only for inspection; the soft blend is the deployable score.