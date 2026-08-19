# BioDarwin 96-well GA/RL search v1

Generated: 2026-05-11T01:51:32

## Locked-control comparison
| algorithm                    |   n |   positives |   positive_rate |    AUPRC |    AUROC |   top5_hits |   top5_precision |   top10_hits |   top10_precision |   top24_hits |   top24_precision |   rank_sum |
|:-----------------------------|----:|------------:|----------------:|---------:|---------:|------------:|-----------------:|-------------:|------------------:|-------------:|------------------:|-----------:|
| BioDarwin_wetlab_search_v1   |  25 |           9 |            0.36 | 0.669359 | 0.701389 |           4 |              0.8 |            5 |               0.5 |            9 |             0.375 |        139 |
| BioDarwin_public_anchor_v3   |  25 |           9 |            0.36 | 0.643014 | 0.680556 |           4 |              0.8 |            5 |               0.5 |            9 |             0.375 |        139 |
| BioDarwin_mode_bank_v4       |  25 |           9 |            0.36 | 0.643014 | 0.680556 |           4 |              0.8 |            5 |               0.5 |            9 |             0.375 |        139 |
| BioDarwin_regime_router_v1   |  25 |           9 |            0.36 | 0.630193 | 0.659722 |           4 |              0.8 |            4 |               0.4 |            9 |             0.375 |        139 |
| BigMHC_IM                    |  25 |           9 |            0.36 | 0.591997 | 0.645833 |           4 |              0.8 |            4 |               0.4 |            9 |             0.375 |        139 |
| GA_public_feature_adapter_v0 |  25 |           9 |            0.36 | 0.524623 | 0.604167 |           3 |              0.6 |            4 |               0.4 |            9 |             0.375 |        139 |

## Best genome
| feature                          |      weight | transform   |
|:---------------------------------|------------:|:------------|
| score_disagreement               | 0.202587    | sigmoid     |
| topline_vote_count               | 0.200592    | square      |
| RF_biophys_score                 | 0.11801     | sqrt        |
| BioDarwin_public_anchor_v3       | 0.0991227   | square      |
| BioDarwin_RF_anchor_heuristic_v1 | 0.0948069   | sqrt        |
| BigMHC_IM                        | 0.0831203   | sqrt        |
| BioDarwin_regime_router_v1       | 0.0715625   | sqrt        |
| score_consensus                  | 0.0617016   | sqrt        |
| BioDarwin_mode_bank_v4           | 0.0615057   | sqrt        |
| BigMHC_EL                        | 0.00611523  | sqrt        |
| GA_public_feature_adapter_v0     | 0.000876466 | square      |

## Summary
| generated_at        |   elapsed_s |   plate_rows |   locked_metric_rows |   best_fitness | best_metrics                                                                                                                                                                                                                                                                                                            | outputs                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
|:--------------------|------------:|-------------:|---------------------:|---------------:|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| 2026-05-11T01:51:32 |       18.22 |           96 |                   25 |       0.664546 | {'n': 25, 'positives': 9, 'positive_rate': 0.36, 'AUPRC': 0.6693585816392835, 'AUROC': 0.7013888888888888, 'top5_hits': 4, 'top5_precision': 0.8, 'top10_hits': 5, 'top10_precision': 0.5, 'top24_hits': 9, 'top24_precision': 0.375, 'rank_sum': 139.0, 'sparsity': 0.8525117509095712, 'fitness': 0.6645462307910095} | {'plate': '/home/seungho/personal/THCA_data_analysis/project/results/biodarwin_96well_ga_rl_search_2026_05_11/biodarwin_96well_ga_rl_search_v1.tsv', 'metrics': '/home/seungho/personal/THCA_data_analysis/project/results/biodarwin_96well_ga_rl_search_2026_05_11/biodarwin_96well_ga_rl_search_v1_metrics.tsv', 'genome': '/home/seungho/personal/THCA_data_analysis/project/results/biodarwin_96well_ga_rl_search_2026_05_11/biodarwin_96well_ga_rl_search_v1_genome.tsv', 'trace': '/home/seungho/personal/THCA_data_analysis/project/results/biodarwin_96well_ga_rl_search_2026_05_11/biodarwin_96well_ga_rl_search_v1_trace.tsv', 'html': '/home/seungho/personal/THCA_data_analysis/project/papers_hub_2026_05_04/biodarwin_96well_ga_rl_search_v1.html'} |