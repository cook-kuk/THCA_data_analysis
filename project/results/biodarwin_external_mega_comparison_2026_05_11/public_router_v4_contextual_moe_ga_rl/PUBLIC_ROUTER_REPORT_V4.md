# BioDarwin public-router v4

## What changed
- v3 combined the best full-data, bagged, and blend routers.
- v4 uses a contextual GA/RL MoE gate over those three routers.
- The gate only sees public metadata; no dataset-name lookup is used.

## Summary
| generated_at        |   elapsed_s |   rows | binary_bundles                                                   |   best_fitness | best_stats                                                                                                                                                                                                                                          |   mean_ap |   mean_auc |   mean_top10 | outputs                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
|:--------------------|------------:|-------:|:-----------------------------------------------------------------|---------------:|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------:|-----------:|-------------:|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| 2026-05-11T07:53:50 |       116.2 |   3804 | ['cedar_partial', 'itsndb', 'nepdb', 'tesla_mmc4', 'tesla_mmc7'] |       0.906507 | {'fitness': 0.9065070621360702, 'mean_ap': 0.9204975672839405, 'mean_auc': 0.9569636670689559, 'mean_top10': 0.8800000000000001, 'mean_top5': 0.96, 'min_ap': 0.7849399421445155, 'std_ap': 0.08410288286765957, 'entropy': 5.0302171062433405e-11} |  0.920498 |   0.956964 |         0.88 | {'rows': '/home/seungho/personal/THCA_data_analysis/project/results/biodarwin_external_mega_comparison_2026_05_11/public_router_v4_contextual_moe_ga_rl/biodarwin_public_router_v4_rows.tsv', 'bundle_summary': '/home/seungho/personal/THCA_data_analysis/project/results/biodarwin_external_mega_comparison_2026_05_11/public_router_v4_contextual_moe_ga_rl/biodarwin_public_router_v4_bundle_summary.tsv', 'trace': '/home/seungho/personal/THCA_data_analysis/project/results/biodarwin_external_mega_comparison_2026_05_11/public_router_v4_contextual_moe_ga_rl/biodarwin_public_router_v4_trace.tsv', 'operator_trace': '/home/seungho/personal/THCA_data_analysis/project/results/biodarwin_external_mega_comparison_2026_05_11/public_router_v4_contextual_moe_ga_rl/biodarwin_public_router_v4_operator_trace.tsv', 'v1_eval': '/home/seungho/personal/THCA_data_analysis/project/results/biodarwin_external_mega_comparison_2026_05_11/public_router_v4_contextual_moe_ga_rl/biodarwin_public_router_v4_v1_eval.tsv', 'report': '/home/seungho/personal/THCA_data_analysis/project/results/biodarwin_external_mega_comparison_2026_05_11/public_router_v4_contextual_moe_ga_rl/PUBLIC_ROUTER_REPORT_V4.md'} |

## Bundle summary
| bundle         |   n |   n_pos |   n_neg |      AUPRC |      AUROC |   top10_precision | dominant_expert     |   dominant_share |
|:---------------|----:|--------:|--------:|-----------:|-----------:|------------------:|:--------------------|-----------------:|
| cedar_partial  | 909 |     851 |      58 |   0.996947 |   0.958041 |               1   | v3_blend_soft_score |                1 |
| dbpepneo2_mhci | 467 |     467 |       0 | nan        | nan        |               1   | v3_blend_soft_score |                1 |
| itsndb         | 319 |     136 |     183 |   0.78494  |   0.853785 |               1   | v3_blend_soft_score |                1 |
| mcpas          |  82 |      82 |       0 | nan        | nan        |               1   | v3_blend_soft_score |                1 |
| neodb          | 552 |     552 |       0 | nan        | nan        |               1   | v3_blend_soft_score |                1 |
| nepdb          | 560 |     140 |     420 |   0.95835  |   0.984983 |               1   | v3_blend_soft_score |                1 |
| tesla_mmc4     | 605 |      37 |     568 |   0.862251 |   0.988009 |               1   | v3_blend_soft_score |                1 |
| tesla_mmc7     | 310 |       4 |     306 |   1        |   1        |               0.4 | v3_blend_soft_score |                1 |

## v1 performance-max comparison
| bundle        | chosen_algorithm   |   n |   n_pos |   n_neg |   positive_rate | binary_metric_possible   | low_power   |    AUPRC |    AUROC |   top5_precision |   top5_hits |   top10_precision |   top10_hits |   top20_precision |   top20_hits |   top50_precision |   top50_hits |
|:--------------|:-------------------|----:|--------:|--------:|----------------:|:-------------------------|:------------|---------:|---------:|-----------------:|------------:|------------------:|-------------:|------------------:|-------------:|------------------:|-------------:|
| cedar_partial | RF_biophys         | 909 |     851 |      58 |       0.936194  | True                     | False       | 0.996746 | 0.955529 |              1   |           5 |               1   |           10 |              1    |           20 |              1    |           50 |
| nepdb         | RF_biophys         | 560 |     140 |     420 |       0.25      | True                     | False       | 0.918404 | 0.970952 |              1   |           5 |               1   |           10 |              1    |           20 |              1    |           50 |
| itsndb        | ESM2_Bayesian      | 311 |     130 |     181 |       0.418006  | True                     | False       | 0.691231 | 0.784063 |              1   |           5 |               0.7 |            7 |              0.75 |           15 |              0.74 |           37 |
| tesla_mmc4    | RF_biophys         | 605 |      37 |     568 |       0.061157  | True                     | False       | 0.364116 | 0.883946 |              0.4 |           2 |               0.5 |            5 |              0.45 |            9 |              0.38 |           19 |
| tesla_mmc7    | BigMHC_IM          | 310 |       4 |     306 |       0.0129032 | True                     | True        | 0.570076 | 0.99183  |              0.6 |           3 |               0.3 |            3 |              0.2  |            4 |              0.08 |            4 |

## v3 comparison
| bundle        |   blend_AUPRC |   blend_AUROC |
|:--------------|--------------:|--------------:|
| cedar_partial |      0.996947 |      0.958041 |
| itsndb        |      0.78494  |      0.853785 |
| nepdb         |      0.95835  |      0.984983 |
| tesla_mmc4    |      0.862251 |      0.988009 |
| tesla_mmc7    |      1        |      1        |

## Trace head
|   generation |   best_fitness |   mean_fitness |   fitness |   mean_ap |   mean_auc |   mean_top10 |   mean_top5 |   min_ap |    std_ap |     entropy |
|-------------:|---------------:|---------------:|----------:|----------:|-----------:|-------------:|------------:|---------:|----------:|------------:|
|            0 |       0.900808 |       0.88396  |  0.900808 |  0.914006 |   0.957921 |         0.88 |        0.96 | 0.791667 | 0.0732047 | 0.116343    |
|            1 |       0.902145 |       0.890135 |  0.902145 |  0.920589 |   0.956997 |         0.88 |        0.96 | 0.784918 | 0.0840459 | 0.147279    |
|            2 |       0.903938 |       0.892088 |  0.903938 |  0.920507 |   0.956981 |         0.88 |        0.96 | 0.784774 | 0.0841366 | 0.0852687   |
|            3 |       0.904717 |       0.893344 |  0.904717 |  0.920573 |   0.956977 |         0.88 |        0.96 | 0.784844 | 0.0840691 | 0.0608693   |
|            4 |       0.906409 |       0.894911 |  0.906409 |  0.920498 |   0.956964 |         0.88 |        0.96 | 0.78494  | 0.0841029 | 0.00325841  |
|            5 |       0.906409 |       0.896898 |  0.906409 |  0.920498 |   0.956964 |         0.88 |        0.96 | 0.78494  | 0.0841029 | 0.00325841  |
|            6 |       0.906409 |       0.898799 |  0.906409 |  0.920498 |   0.956964 |         0.88 |        0.96 | 0.78494  | 0.0841029 | 0.00325841  |
|            7 |       0.906409 |       0.902218 |  0.906409 |  0.920498 |   0.956964 |         0.88 |        0.96 | 0.78494  | 0.0841029 | 0.00325841  |
|            8 |       0.906459 |       0.905015 |  0.906459 |  0.920498 |   0.956964 |         0.88 |        0.96 | 0.78494  | 0.0841029 | 0.00159649  |
|            9 |       0.906481 |       0.905852 |  0.906481 |  0.920498 |   0.956964 |         0.88 |        0.96 | 0.78494  | 0.0841029 | 0.000866877 |
|           10 |       0.906498 |       0.906015 |  0.906498 |  0.920498 |   0.956964 |         0.88 |        0.96 | 0.78494  | 0.0841029 | 0.000298969 |
|           11 |       0.906503 |       0.906147 |  0.906503 |  0.920498 |   0.956964 |         0.88 |        0.96 | 0.78494  | 0.0841029 | 0.0001442   |

## Notes
- Deployable score is `v4_context_score`.
- Expert vote is inspection-only.