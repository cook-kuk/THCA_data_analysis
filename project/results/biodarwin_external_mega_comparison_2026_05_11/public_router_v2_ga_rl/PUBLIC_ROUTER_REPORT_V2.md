# BioDarwin public-router v2

## What changed
- v1 was a dataset-family lookup table.
- v2 is a GA/RL mixture-of-experts router that uses only public metadata for gating.
- The experts are local predictor outputs; the gate sees peptide length, overlap, peptide priors, HLA family, source family, and validation type.

## Champion summary
| generated_at        |   elapsed_s |   rows | bundles                                                                                              | binary_bundles                                                   | experts                                                                                                                                                                                          |   best_fitness | best_stats                                                                                                                                                                                                                                                            |   temperature | outputs                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |
|:--------------------|------------:|-------:|:-----------------------------------------------------------------------------------------------------|:-----------------------------------------------------------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------:|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|--------------:|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| 2026-05-11T07:24:10 |      113.25 |   3804 | ['cedar_partial', 'dbpepneo2_mhci', 'itsndb', 'mcpas', 'neodb', 'nepdb', 'tesla_mmc4', 'tesla_mmc7'] | ['cedar_partial', 'itsndb', 'nepdb', 'tesla_mmc4', 'tesla_mmc7'] | ['BigMHC_IM', 'RF_biophys', 'LR_biophys', 'MHCflurry', 'PRIME', 'ESM2_Bayesian', 'DeepImmuno', 'TransPHLA', 'NetMHCpan', 'GP_quantum', 'VQC', 'Structure_LR', 'BioDarwin_public_anchor_noEL_v4'] |       0.893948 | {'fitness': 0.8939482030481778, 'mean_ap': 0.9186111749635243, 'mean_auc': 0.9566256224461608, 'mean_top10': 0.8800000000000001, 'mean_top5': 0.96, 'min_ap': 0.7831796219787388, 'std_ap': 0.0856223991853251, 'entropy': 0.3992904128245127, 'bundle_eval_rows': 5} |      0.675398 | {'rows': '/home/seungho/personal/THCA_data_analysis/project/results/biodarwin_external_mega_comparison_2026_05_11/public_router_v2_ga_rl/biodarwin_public_router_v2_rows.tsv', 'bundle_summary': '/home/seungho/personal/THCA_data_analysis/project/results/biodarwin_external_mega_comparison_2026_05_11/public_router_v2_ga_rl/biodarwin_public_router_v2_bundle_summary.tsv', 'single_metrics': '/home/seungho/personal/THCA_data_analysis/project/results/biodarwin_external_mega_comparison_2026_05_11/public_router_v2_ga_rl/biodarwin_public_router_v2_single_expert_metrics.tsv', 'best_by_bundle': '/home/seungho/personal/THCA_data_analysis/project/results/biodarwin_external_mega_comparison_2026_05_11/public_router_v2_ga_rl/biodarwin_public_router_v2_best_single_by_bundle.tsv', 'genome': '/home/seungho/personal/THCA_data_analysis/project/results/biodarwin_external_mega_comparison_2026_05_11/public_router_v2_ga_rl/biodarwin_public_router_v2_genome.tsv', 'bias': '/home/seungho/personal/THCA_data_analysis/project/results/biodarwin_external_mega_comparison_2026_05_11/public_router_v2_ga_rl/biodarwin_public_router_v2_bias.tsv', 'trace': '/home/seungho/personal/THCA_data_analysis/project/results/biodarwin_external_mega_comparison_2026_05_11/public_router_v2_ga_rl/biodarwin_public_router_v2_trace.tsv', 'operator_trace': '/home/seungho/personal/THCA_data_analysis/project/results/biodarwin_external_mega_comparison_2026_05_11/public_router_v2_ga_rl/biodarwin_public_router_v2_operator_trace.tsv', 'report': '/home/seungho/personal/THCA_data_analysis/project/results/biodarwin_external_mega_comparison_2026_05_11/public_router_v2_ga_rl/PUBLIC_ROUTER_REPORT_V2.md'} |

## Binary public bundle comparison
| bundle        |   n |   n_pos |   n_neg | binary_metric_possible   |   soft_AUPRC |   soft_AUROC |   soft_top10_precision |   hard_AUPRC |   hard_AUROC |   hard_top10_precision | dominant_hard_expert            |   dominant_hard_share | algorithm     |    AUPRC |    AUROC |   top10_precision | chosen_algorithm   |   AUPRC_v1 |   AUROC_v1 |   top10_precision_v1 |
|:--------------|----:|--------:|--------:|:-------------------------|-------------:|-------------:|-----------------------:|-------------:|-------------:|-----------------------:|:--------------------------------|----------------------:|:--------------|---------:|---------:|------------------:|:-------------------|-----------:|-----------:|---------------------:|
| cedar_partial | 909 |     851 |      58 | True                     |     0.996942 |     0.958001 |                    1   |     0.996746 |     0.955529 |                    1   | RF_biophys                      |              1        | RF_biophys    | 0.996746 | 0.955529 |               1   | RF_biophys         |   0.996746 |   0.955529 |                  1   |
| itsndb        | 319 |     136 |     183 | True                     |     0.78318  |     0.853182 |                    1   |     0.769488 |     0.817944 |                    1   | BioDarwin_public_anchor_noEL_v4 |              0.721003 | ESM2_Bayesian | 0.691231 | 0.784063 |               0.7 | ESM2_Bayesian      |   0.691231 |   0.784063 |                  0.7 |
| nepdb         | 560 |     140 |     420 | True                     |     0.957747 |     0.984745 |                    1   |     0.918404 |     0.970952 |                    1   | RF_biophys                      |              0.998214 | RF_biophys    | 0.918404 | 0.970952 |               1   | RF_biophys         |   0.918404 |   0.970952 |                  1   |
| tesla_mmc4    | 605 |      37 |     568 | True                     |     0.855186 |     0.9872   |                    1   |     0.236755 |     0.868553 |                    0.4 | RF_biophys                      |              0.487603 | RF_biophys    | 0.364116 | 0.883946 |               0.5 | RF_biophys         |   0.364116 |   0.883946 |                  0.5 |
| tesla_mmc7    | 310 |       4 |     306 | True                     |     1        |     1        |                    0.4 |     0.101361 |     0.930556 |                    0   | RF_biophys                      |              0.683871 | BigMHC_IM     | 0.570076 | 0.99183  |               0.3 | BigMHC_IM          |   0.570076 |   0.99183  |                  0.3 |

## v1 performance-max router
| bundle        | chosen_algorithm   |   n |   n_pos |   n_neg |   positive_rate | binary_metric_possible   | low_power   |    AUPRC |    AUROC |   top5_precision |   top5_hits |   top10_precision |   top10_hits |   top20_precision |   top20_hits |   top50_precision |   top50_hits |
|:--------------|:-------------------|----:|--------:|--------:|----------------:|:-------------------------|:------------|---------:|---------:|-----------------:|------------:|------------------:|-------------:|------------------:|-------------:|------------------:|-------------:|
| cedar_partial | RF_biophys         | 909 |     851 |      58 |       0.936194  | True                     | False       | 0.996746 | 0.955529 |              1   |           5 |               1   |           10 |              1    |           20 |              1    |           50 |
| nepdb         | RF_biophys         | 560 |     140 |     420 |       0.25      | True                     | False       | 0.918404 | 0.970952 |              1   |           5 |               1   |           10 |              1    |           20 |              1    |           50 |
| itsndb        | ESM2_Bayesian      | 311 |     130 |     181 |       0.418006  | True                     | False       | 0.691231 | 0.784063 |              1   |           5 |               0.7 |            7 |              0.75 |           15 |              0.74 |           37 |
| tesla_mmc4    | RF_biophys         | 605 |      37 |     568 |       0.061157  | True                     | False       | 0.364116 | 0.883946 |              0.4 |           2 |               0.5 |            5 |              0.45 |            9 |              0.38 |           19 |
| tesla_mmc7    | BigMHC_IM          | 310 |       4 |     306 |       0.0129032 | True                     | True        | 0.570076 | 0.99183  |              0.6 |           3 |               0.3 |            3 |              0.2  |            4 |              0.08 |            4 |

## v1 reviewer-safe router
| bundle        | chosen_algorithm   |   n |   n_pos |   n_neg |   positive_rate | binary_metric_possible   | low_power   |    AUPRC |    AUROC |   top5_precision |   top5_hits |   top10_precision |   top10_hits |   top20_precision |   top20_hits |   top50_precision |   top50_hits |
|:--------------|:-------------------|----:|--------:|--------:|----------------:|:-------------------------|:------------|---------:|---------:|-----------------:|------------:|------------------:|-------------:|------------------:|-------------:|------------------:|-------------:|
| cedar_partial | MHCflurry          | 845 |     787 |      58 |        0.931361 | True                     | False       | 0.972759 | 0.737173 |              1   |           5 |               1   |           10 |               1   |           20 |              1    |           50 |
| nepdb         | PRIME              | 483 |     140 |     343 |        0.289855 | True                     | False       | 0.511098 | 0.71893  |              0.8 |           4 |               0.8 |            8 |               0.7 |           14 |              0.54 |           27 |
| itsndb        | DeepImmuno         | 319 |     136 |     183 |        0.426332 | True                     | False       | 0.610914 | 0.701905 |              0.8 |           4 |               0.8 |            8 |               0.7 |           14 |              0.58 |           29 |
| tesla_mmc4    | PRIME              | 605 |      37 |     568 |        0.061157 | True                     | False       | 0.193735 | 0.805434 |              0.2 |           1 |               0.2 |            2 |               0.2 |            4 |              0.16 |            8 |

## Best single expert per bundle
| bundle        | algorithm     |   n |   n_pos |   n_neg |   positive_rate | binary_metric_possible   | low_power   |    AUPRC |    AUROC |   top5_precision |   top5_hits |   top10_precision |   top10_hits |   top20_precision |   top20_hits |   top50_precision |   top50_hits |
|:--------------|:--------------|----:|--------:|--------:|----------------:|:-------------------------|:------------|---------:|---------:|-----------------:|------------:|------------------:|-------------:|------------------:|-------------:|------------------:|-------------:|
| cedar_partial | RF_biophys    | 909 |     851 |      58 |       0.936194  | True                     | False       | 0.996746 | 0.955529 |              1   |           5 |               1   |           10 |              1    |           20 |              1    |           50 |
| itsndb        | ESM2_Bayesian | 311 |     130 |     181 |       0.418006  | True                     | False       | 0.691231 | 0.784063 |              1   |           5 |               0.7 |            7 |              0.75 |           15 |              0.74 |           37 |
| nepdb         | RF_biophys    | 560 |     140 |     420 |       0.25      | True                     | False       | 0.918404 | 0.970952 |              1   |           5 |               1   |           10 |              1    |           20 |              1    |           50 |
| tesla_mmc4    | RF_biophys    | 605 |      37 |     568 |       0.061157  | True                     | False       | 0.364116 | 0.883946 |              0.4 |           2 |               0.5 |            5 |              0.45 |            9 |              0.38 |           19 |
| tesla_mmc7    | BigMHC_IM     | 310 |       4 |     306 |       0.0129032 | True                     | True        | 0.570076 | 0.99183  |              0.6 |           3 |               0.3 |            3 |              0.2  |            4 |              0.08 |            4 |

## Single-expert baseline table
| bundle        | algorithm                       |   n |   n_pos |   n_neg |   positive_rate | binary_metric_possible   | low_power   |       AUPRC |      AUROC |   top5_precision |   top5_hits |   top10_precision |   top10_hits |   top20_precision |   top20_hits |   top50_precision |   top50_hits |
|:--------------|:--------------------------------|----:|--------:|--------:|----------------:|:-------------------------|:------------|------------:|-----------:|-----------------:|------------:|------------------:|-------------:|------------------:|-------------:|------------------:|-------------:|
| cedar_partial | RF_biophys                      | 909 |     851 |      58 |       0.936194  | True                     | False       |   0.996746  |   0.955529 |              1   |           5 |               1   |           10 |              1    |           20 |              1    |           50 |
| cedar_partial | LR_biophys                      | 909 |     851 |      58 |       0.936194  | True                     | False       |   0.988603  |   0.848495 |              1   |           5 |               1   |           10 |              1    |           20 |              1    |           50 |
| cedar_partial | MHCflurry                       | 845 |     787 |      58 |       0.931361  | True                     | False       |   0.972759  |   0.737173 |              1   |           5 |               1   |           10 |              1    |           20 |              1    |           50 |
| cedar_partial | PRIME                           | 898 |     840 |      58 |       0.935412  | True                     | False       |   0.967722  |   0.683929 |              1   |           5 |               1   |           10 |              1    |           20 |              0.98 |           49 |
| cedar_partial | TransPHLA                       | 736 |     678 |      58 |       0.921196  | True                     | False       |   0.964784  |   0.728105 |              1   |           5 |               1   |           10 |              1    |           20 |              1    |           50 |
| cedar_partial | BigMHC_IM                       | 898 |     840 |      58 |       0.935412  | True                     | False       |   0.954816  |   0.542262 |              1   |           5 |               1   |           10 |              1    |           20 |              1    |           50 |
| cedar_partial | BioDarwin_public_anchor_noEL_v4 | 909 |     851 |      58 |       0.936194  | True                     | False       |   0.954569  |   0.534787 |              1   |           5 |               1   |           10 |              1    |           20 |              1    |           50 |
| cedar_partial | ESM2_Bayesian                   |   0 |       0 |       0 |     nan         | False                    | True        | nan         | nan        |            nan   |           0 |             nan   |            0 |            nan    |            0 |            nan    |            0 |
| cedar_partial | DeepImmuno                      |   0 |       0 |       0 |     nan         | False                    | True        | nan         | nan        |            nan   |           0 |             nan   |            0 |            nan    |            0 |            nan    |            0 |
| cedar_partial | NetMHCpan                       |   0 |       0 |       0 |     nan         | False                    | True        | nan         | nan        |            nan   |           0 |             nan   |            0 |            nan    |            0 |            nan    |            0 |
| cedar_partial | GP_quantum                      |   0 |       0 |       0 |     nan         | False                    | True        | nan         | nan        |            nan   |           0 |             nan   |            0 |            nan    |            0 |            nan    |            0 |
| cedar_partial | VQC                             |   0 |       0 |       0 |     nan         | False                    | True        | nan         | nan        |            nan   |           0 |             nan   |            0 |            nan    |            0 |            nan    |            0 |
| cedar_partial | Structure_LR                    |   0 |       0 |       0 |     nan         | False                    | True        | nan         | nan        |            nan   |           0 |             nan   |            0 |            nan    |            0 |            nan    |            0 |
| itsndb        | ESM2_Bayesian                   | 311 |     130 |     181 |       0.418006  | True                     | False       |   0.691231  |   0.784063 |              1   |           5 |               0.7 |            7 |              0.75 |           15 |              0.74 |           37 |
| itsndb        | RF_biophys                      | 319 |     136 |     183 |       0.426332  | True                     | False       |   0.650115  |   0.753134 |              0.8 |           4 |               0.7 |            7 |              0.75 |           15 |              0.7  |           35 |
| itsndb        | LR_biophys                      | 319 |     136 |     183 |       0.426332  | True                     | False       |   0.645996  |   0.692462 |              1   |           5 |               1   |           10 |              1    |           20 |              0.68 |           34 |
| itsndb        | BioDarwin_public_anchor_noEL_v4 | 319 |     136 |     183 |       0.426332  | True                     | False       |   0.631129  |   0.745018 |              0.4 |           2 |               0.5 |            5 |              0.6  |           12 |              0.72 |           36 |
| itsndb        | BigMHC_IM                       | 319 |     136 |     183 |       0.426332  | True                     | False       |   0.628389  |   0.74767  |              0.4 |           2 |               0.5 |            5 |              0.6  |           12 |              0.74 |           37 |
| itsndb        | Structure_LR                    | 311 |     130 |     181 |       0.418006  | True                     | False       |   0.618189  |   0.689375 |              1   |           5 |               0.9 |            9 |              0.85 |           17 |              0.62 |           31 |
| itsndb        | DeepImmuno                      | 319 |     136 |     183 |       0.426332  | True                     | False       |   0.610914  |   0.701905 |              0.8 |           4 |               0.8 |            8 |              0.7  |           14 |              0.58 |           29 |
| itsndb        | MHCflurry                       | 319 |     136 |     183 |       0.426332  | True                     | False       |   0.58833   |   0.636733 |              0.8 |           4 |               0.8 |            8 |              0.8  |           16 |              0.7  |           35 |
| itsndb        | GP_quantum                      | 319 |     136 |     183 |       0.426332  | True                     | False       |   0.544901  |   0.58699  |              1   |           5 |               0.9 |            9 |              0.8  |           16 |              0.54 |           27 |
| itsndb        | VQC                             | 319 |     136 |     183 |       0.426332  | True                     | False       |   0.507296  |   0.598521 |              0.2 |           1 |               0.5 |            5 |              0.7  |           14 |              0.56 |           28 |
| itsndb        | TransPHLA                       | 318 |     135 |     183 |       0.424528  | True                     | False       |   0.49779   |   0.5797   |              0.4 |           2 |               0.4 |            4 |              0.45 |            9 |              0.54 |           27 |
| itsndb        | NetMHCpan                       | 319 |     136 |     183 |       0.426332  | True                     | False       |   0.493205  |   0.56688  |              0.6 |           3 |               0.5 |            5 |              0.55 |           11 |              0.52 |           26 |
| itsndb        | PRIME                           | 319 |     136 |     183 |       0.426332  | True                     | False       |   0.480376  |   0.567302 |              0.8 |           4 |               0.4 |            4 |              0.55 |           11 |              0.48 |           24 |
| nepdb         | RF_biophys                      | 560 |     140 |     420 |       0.25      | True                     | False       |   0.918404  |   0.970952 |              1   |           5 |               1   |           10 |              1    |           20 |              1    |           50 |
| nepdb         | BigMHC_IM                       | 486 |     140 |     346 |       0.288066  | True                     | False       |   0.762466  |   0.86602  |              1   |           5 |               0.9 |            9 |              0.95 |           19 |              0.9  |           45 |
| nepdb         | BioDarwin_public_anchor_noEL_v4 | 560 |     140 |     420 |       0.25      | True                     | False       |   0.761612  |   0.884082 |              1   |           5 |               0.9 |            9 |              0.95 |           19 |              0.9  |           45 |
| nepdb         | LR_biophys                      | 560 |     140 |     420 |       0.25      | True                     | False       |   0.709765  |   0.855629 |              1   |           5 |               1   |           10 |              1    |           20 |              0.88 |           44 |
| nepdb         | MHCflurry                       | 402 |     140 |     262 |       0.348259  | True                     | False       |   0.552479  |   0.693511 |              0.6 |           3 |               0.7 |            7 |              0.6  |           12 |              0.66 |           33 |
| nepdb         | PRIME                           | 483 |     140 |     343 |       0.289855  | True                     | False       |   0.511098  |   0.71893  |              0.8 |           4 |               0.8 |            8 |              0.7  |           14 |              0.54 |           27 |
| nepdb         | TransPHLA                       | 550 |     139 |     411 |       0.252727  | True                     | False       |   0.411373  |   0.718786 |              0.6 |           3 |               0.5 |            5 |              0.55 |           11 |              0.46 |           23 |
| nepdb         | ESM2_Bayesian                   |   0 |       0 |       0 |     nan         | False                    | True        | nan         | nan        |            nan   |           0 |             nan   |            0 |            nan    |            0 |            nan    |            0 |
| nepdb         | DeepImmuno                      |   0 |       0 |       0 |     nan         | False                    | True        | nan         | nan        |            nan   |           0 |             nan   |            0 |            nan    |            0 |            nan    |            0 |
| nepdb         | NetMHCpan                       |   0 |       0 |       0 |     nan         | False                    | True        | nan         | nan        |            nan   |           0 |             nan   |            0 |            nan    |            0 |            nan    |            0 |
| nepdb         | GP_quantum                      |   0 |       0 |       0 |     nan         | False                    | True        | nan         | nan        |            nan   |           0 |             nan   |            0 |            nan    |            0 |            nan    |            0 |
| nepdb         | VQC                             |   0 |       0 |       0 |     nan         | False                    | True        | nan         | nan        |            nan   |           0 |             nan   |            0 |            nan    |            0 |            nan    |            0 |
| nepdb         | Structure_LR                    |   0 |       0 |       0 |     nan         | False                    | True        | nan         | nan        |            nan   |           0 |             nan   |            0 |            nan    |            0 |            nan    |            0 |
| tesla_mmc4    | RF_biophys                      | 605 |      37 |     568 |       0.061157  | True                     | False       |   0.364116  |   0.883946 |              0.4 |           2 |               0.5 |            5 |              0.45 |            9 |              0.38 |           19 |
| tesla_mmc4    | BioDarwin_public_anchor_noEL_v4 | 605 |      37 |     568 |       0.061157  | True                     | False       |   0.297665  |   0.866102 |              0.4 |           2 |               0.4 |            4 |              0.4  |            8 |              0.32 |           16 |
| tesla_mmc4    | BigMHC_IM                       | 605 |      37 |     568 |       0.061157  | True                     | False       |   0.291388  |   0.864151 |              0.4 |           2 |               0.4 |            4 |              0.35 |            7 |              0.32 |           16 |
| tesla_mmc4    | PRIME                           | 605 |      37 |     568 |       0.061157  | True                     | False       |   0.193735  |   0.805434 |              0.2 |           1 |               0.2 |            2 |              0.2  |            4 |              0.16 |            8 |
| tesla_mmc4    | MHCflurry                       | 599 |      37 |     562 |       0.0617696 | True                     | False       |   0.188558  |   0.75642  |              0.2 |           1 |               0.1 |            1 |              0.2  |            4 |              0.26 |           13 |
| tesla_mmc4    | TransPHLA                       | 605 |      37 |     568 |       0.061157  | True                     | False       |   0.1132    |   0.722069 |              0   |           0 |               0.2 |            2 |              0.2  |            4 |              0.14 |            7 |
| tesla_mmc4    | LR_biophys                      | 605 |      37 |     568 |       0.061157  | True                     | False       |   0.103764  |   0.669014 |              0   |           0 |               0   |            0 |              0    |            0 |              0.1  |            5 |
| tesla_mmc4    | ESM2_Bayesian                   |   0 |       0 |       0 |     nan         | False                    | True        | nan         | nan        |            nan   |           0 |             nan   |            0 |            nan    |            0 |            nan    |            0 |
| tesla_mmc4    | DeepImmuno                      |   0 |       0 |       0 |     nan         | False                    | True        | nan         | nan        |            nan   |           0 |             nan   |            0 |            nan    |            0 |            nan    |            0 |
| tesla_mmc4    | NetMHCpan                       |   0 |       0 |       0 |     nan         | False                    | True        | nan         | nan        |            nan   |           0 |             nan   |            0 |            nan    |            0 |            nan    |            0 |
| tesla_mmc4    | GP_quantum                      |   0 |       0 |       0 |     nan         | False                    | True        | nan         | nan        |            nan   |           0 |             nan   |            0 |            nan    |            0 |            nan    |            0 |
| tesla_mmc4    | VQC                             |   0 |       0 |       0 |     nan         | False                    | True        | nan         | nan        |            nan   |           0 |             nan   |            0 |            nan    |            0 |            nan    |            0 |
| tesla_mmc4    | Structure_LR                    |   0 |       0 |       0 |     nan         | False                    | True        | nan         | nan        |            nan   |           0 |             nan   |            0 |            nan    |            0 |            nan    |            0 |
| tesla_mmc7    | BigMHC_IM                       | 310 |       4 |     306 |       0.0129032 | True                     | True        |   0.570076  |   0.99183  |              0.6 |           3 |               0.3 |            3 |              0.2  |            4 |              0.08 |            4 |
| tesla_mmc7    | BioDarwin_public_anchor_noEL_v4 | 310 |       4 |     306 |       0.0129032 | True                     | True        |   0.570076  |   0.99183  |              0.6 |           3 |               0.3 |            3 |              0.2  |            4 |              0.08 |            4 |
| tesla_mmc7    | RF_biophys                      | 310 |       4 |     306 |       0.0129032 | True                     | True        |   0.313725  |   0.947712 |              0.4 |           2 |               0.2 |            2 |              0.15 |            3 |              0.06 |            3 |
| tesla_mmc7    | PRIME                           | 310 |       4 |     306 |       0.0129032 | True                     | True        |   0.222222  |   0.901144 |              0.2 |           1 |               0.2 |            2 |              0.1  |            2 |              0.06 |            3 |
| tesla_mmc7    | MHCflurry                       | 310 |       4 |     306 |       0.0129032 | True                     | True        |   0.214057  |   0.966503 |              0   |           0 |               0.2 |            2 |              0.15 |            3 |              0.08 |            4 |
| tesla_mmc7    | LR_biophys                      | 310 |       4 |     306 |       0.0129032 | True                     | True        |   0.0417557 |   0.701797 |              0   |           0 |               0   |            0 |              0.05 |            1 |              0.02 |            1 |
| tesla_mmc7    | TransPHLA                       | 310 |       4 |     306 |       0.0129032 | True                     | True        |   0.036211  |   0.814951 |              0   |           0 |               0   |            0 |              0    |            0 |              0.02 |            1 |
| tesla_mmc7    | ESM2_Bayesian                   |   0 |       0 |       0 |     nan         | False                    | True        | nan         | nan        |            nan   |           0 |             nan   |            0 |            nan    |            0 |            nan    |            0 |
| tesla_mmc7    | DeepImmuno                      |   0 |       0 |       0 |     nan         | False                    | True        | nan         | nan        |            nan   |           0 |             nan   |            0 |            nan    |            0 |            nan    |            0 |
| tesla_mmc7    | NetMHCpan                       |   0 |       0 |       0 |     nan         | False                    | True        | nan         | nan        |            nan   |           0 |             nan   |            0 |            nan    |            0 |            nan    |            0 |
| tesla_mmc7    | GP_quantum                      |   0 |       0 |       0 |     nan         | False                    | True        | nan         | nan        |            nan   |           0 |             nan   |            0 |            nan    |            0 |            nan    |            0 |
| tesla_mmc7    | VQC                             |   0 |       0 |       0 |     nan         | False                    | True        | nan         | nan        |            nan   |           0 |             nan   |            0 |            nan    |            0 |            nan    |            0 |
| tesla_mmc7    | Structure_LR                    |   0 |       0 |       0 |     nan         | False                    | True        | nan         | nan        |            nan   |           0 |             nan   |            0 |            nan    |            0 |            nan    |            0 |

## Genome
| expert                          | feature                    |      weight |
|:--------------------------------|:---------------------------|------------:|
| BigMHC_IM                       | hla_A                      |  0.305551   |
| BigMHC_IM                       | validation_unknown         |  0.2151     |
| BigMHC_IM                       | length_sq_norm             |  0.210436   |
| BigMHC_IM                       | hla_other                  |  0.191163   |
| BigMHC_IM                       | seq_helper_prior           |  0.17444    |
| BigMHC_IM                       | source_improve             |  0.0913481  |
| BigMHC_IM                       | source_other               |  0.0730744  |
| BigMHC_IM                       | charged_balance            |  0.0704837  |
| BigMHC_IM                       | source_cedar               |  0.0541773  |
| BigMHC_IM                       | bias                       |  0.0320392  |
| BigMHC_IM                       | source_mcpas               |  0.0152317  |
| BigMHC_IM                       | hydrophobic                |  0          |
| BigMHC_IM                       | validation_held_out        | -0.00337924 |
| BigMHC_IM                       | length_norm                | -0.0060065  |
| BigMHC_IM                       | helper_promiscuity_proxy   | -0.0104036  |
| BigMHC_IM                       | source_dbpepneo2           | -0.0180487  |
| BigMHC_IM                       | hla_C                      | -0.0241983  |
| BigMHC_IM                       | source_tesla               | -0.0942526  |
| BigMHC_IM                       | slp_processability_proxy   | -0.196353   |
| BigMHC_IM                       | source_itsndb              | -0.244378   |
| BigMHC_IM                       | validation_partial_overlap | -0.283006   |
| BigMHC_IM                       | aromatic                   | -0.332963   |
| BigMHC_IM                       | seq_diversity              | -0.416658   |
| BigMHC_IM                       | seq_cytotoxic_prior        | -0.483407   |
| BigMHC_IM                       | validation_external_test   | -0.543745   |
| BigMHC_IM                       | hla_B                      | -0.652335   |
| BigMHC_IM                       | source_neodb               | -0.683605   |
| BigMHC_IM                       | source_nepdb               | -1.09424    |
| BigMHC_IM                       | overlap_norm               | -1.29503    |
| BioDarwin_public_anchor_noEL_v4 | helper_promiscuity_proxy   |  0.867697   |
| BioDarwin_public_anchor_noEL_v4 | seq_cytotoxic_prior        |  0.552871   |
| BioDarwin_public_anchor_noEL_v4 | hla_A                      |  0.427939   |
| BioDarwin_public_anchor_noEL_v4 | bias                       |  0.336027   |
| BioDarwin_public_anchor_noEL_v4 | source_cedar               |  0.311818   |
| BioDarwin_public_anchor_noEL_v4 | seq_helper_prior           |  0.270427   |
| BioDarwin_public_anchor_noEL_v4 | charged_balance            |  0.259005   |
| BioDarwin_public_anchor_noEL_v4 | source_improve             |  0.230193   |
| BioDarwin_public_anchor_noEL_v4 | hydrophobic                |  0.145024   |
| BioDarwin_public_anchor_noEL_v4 | source_dbpepneo2           |  0.0974488  |
| BioDarwin_public_anchor_noEL_v4 | source_itsndb              |  0.0734806  |
| BioDarwin_public_anchor_noEL_v4 | source_nepdb               |  0.058087   |
| BioDarwin_public_anchor_noEL_v4 | length_norm                |  0.0559569  |
| BioDarwin_public_anchor_noEL_v4 | seq_diversity              |  0.0455193  |
| BioDarwin_public_anchor_noEL_v4 | source_other               |  0.0444912  |
| BioDarwin_public_anchor_noEL_v4 | validation_held_out        |  0.0419857  |
| BioDarwin_public_anchor_noEL_v4 | source_mcpas               |  0.0321946  |
| BioDarwin_public_anchor_noEL_v4 | hla_C                      | -0.00579486 |
| BioDarwin_public_anchor_noEL_v4 | validation_unknown         | -0.0344043  |
| BioDarwin_public_anchor_noEL_v4 | length_sq_norm             | -0.0713919  |
| BioDarwin_public_anchor_noEL_v4 | source_tesla               | -0.0908309  |
| BioDarwin_public_anchor_noEL_v4 | validation_external_test   | -0.0996735  |
| BioDarwin_public_anchor_noEL_v4 | hla_B                      | -0.105949   |
| BioDarwin_public_anchor_noEL_v4 | hla_other                  | -0.194667   |
| BioDarwin_public_anchor_noEL_v4 | slp_processability_proxy   | -0.285323   |
| BioDarwin_public_anchor_noEL_v4 | aromatic                   | -0.337515   |
| BioDarwin_public_anchor_noEL_v4 | validation_partial_overlap | -0.454729   |
| BioDarwin_public_anchor_noEL_v4 | source_neodb               | -0.576517   |
| BioDarwin_public_anchor_noEL_v4 | overlap_norm               | -1.40837    |
| DeepImmuno                      | overlap_norm               |  0.705037   |
| DeepImmuno                      | seq_diversity              |  0.559901   |
| DeepImmuno                      | length_sq_norm             |  0.429108   |
| DeepImmuno                      | slp_processability_proxy   |  0.418395   |
| DeepImmuno                      | source_itsndb              |  0.36877    |
| DeepImmuno                      | source_nepdb               |  0.334321   |
| DeepImmuno                      | hla_B                      |  0.198696   |
| DeepImmuno                      | source_mcpas               |  0.192751   |
| DeepImmuno                      | source_tesla               |  0.143699   |
| DeepImmuno                      | source_neodb               |  0.0855308  |
| DeepImmuno                      | source_cedar               |  0.0535421  |
| DeepImmuno                      | source_dbpepneo2           |  0          |
| DeepImmuno                      | length_norm                | -0.0423632  |
| DeepImmuno                      | hla_A                      | -0.0461144  |
| DeepImmuno                      | source_improve             | -0.0543229  |
| DeepImmuno                      | seq_cytotoxic_prior        | -0.0628317  |
| DeepImmuno                      | seq_helper_prior           | -0.0776421  |
| DeepImmuno                      | aromatic                   | -0.0831382  |
| DeepImmuno                      | bias                       | -0.107213   |
| DeepImmuno                      | hla_C                      | -0.130638   |
| DeepImmuno                      | validation_partial_overlap | -0.141028   |
| DeepImmuno                      | validation_external_test   | -0.175421   |
| DeepImmuno                      | validation_held_out        | -0.1873     |
| DeepImmuno                      | charged_balance            | -0.199278   |
| DeepImmuno                      | helper_promiscuity_proxy   | -0.211412   |
| DeepImmuno                      | validation_unknown         | -0.242737   |
| DeepImmuno                      | hydrophobic                | -0.305112   |
| DeepImmuno                      | source_other               | -0.380991   |
| DeepImmuno                      | hla_other                  | -0.8267     |
| ESM2_Bayesian                   | seq_cytotoxic_prior        |  0.889763   |
| ESM2_Bayesian                   | slp_processability_proxy   |  0.746599   |
| ESM2_Bayesian                   | source_neodb               |  0.601132   |
| ESM2_Bayesian                   | hla_C                      |  0.56832    |
| ESM2_Bayesian                   | source_other               |  0.438484   |
| ESM2_Bayesian                   | validation_unknown         |  0.327939   |
| ESM2_Bayesian                   | hla_B                      |  0.28465    |
| ESM2_Bayesian                   | aromatic                   |  0.225451   |
| ESM2_Bayesian                   | seq_diversity              |  0.196266   |
| ESM2_Bayesian                   | length_sq_norm             |  0.186556   |
| ESM2_Bayesian                   | source_itsndb              |  0.173179   |
| ESM2_Bayesian                   | source_nepdb               |  0.112077   |
| ESM2_Bayesian                   | charged_balance            |  0.0986507  |
| ESM2_Bayesian                   | source_dbpepneo2           |  0.0657908  |
| ESM2_Bayesian                   | length_norm                |  0.00859244 |
| ESM2_Bayesian                   | source_improve             |  0.00527186 |
| ESM2_Bayesian                   | source_mcpas               |  0          |
| ESM2_Bayesian                   | source_cedar               | -0.0206672  |
| ESM2_Bayesian                   | hla_A                      | -0.0247702  |
| ESM2_Bayesian                   | validation_held_out        | -0.0399594  |
| ESM2_Bayesian                   | hla_other                  | -0.0442061  |
| ESM2_Bayesian                   | helper_promiscuity_proxy   | -0.07394    |
| ESM2_Bayesian                   | validation_partial_overlap | -0.222818   |
| ESM2_Bayesian                   | hydrophobic                | -0.239623   |
| ESM2_Bayesian                   | seq_helper_prior           | -0.287485   |
| ESM2_Bayesian                   | source_tesla               | -0.35559    |
| ESM2_Bayesian                   | overlap_norm               | -0.417983   |
| ESM2_Bayesian                   | validation_external_test   | -0.444901   |
| ESM2_Bayesian                   | bias                       | -0.554016   |
| GP_quantum                      | source_tesla               |  0.579305   |
| GP_quantum                      | length_sq_norm             |  0.368149   |
| GP_quantum                      | seq_cytotoxic_prior        |  0.281277   |
| GP_quantum                      | seq_helper_prior           |  0.26904    |

## Bias
| expert                          |         bias |
|:--------------------------------|-------------:|
| BigMHC_IM                       | -0.0143273   |
| RF_biophys                      | -0.0139286   |
| LR_biophys                      |  0.0802489   |
| MHCflurry                       | -0.000311444 |
| PRIME                           | -0.0585311   |
| ESM2_Bayesian                   | -0.0424101   |
| DeepImmuno                      |  0.0514782   |
| TransPHLA                       | -0.00727983  |
| NetMHCpan                       |  0.00280821  |
| GP_quantum                      | -0.0454367   |
| VQC                             | -0.078222    |
| Structure_LR                    | -0.078222    |
| BioDarwin_public_anchor_noEL_v4 |  0.0332879   |

## Operator trace head
|   generation |   best_fitness |   mean_fitness |   fitness |   mean_ap |   mean_auc |   mean_top10 |   mean_top5 |   min_ap |   std_ap |   entropy |   bundle_eval_rows |
|-------------:|---------------:|---------------:|----------:|----------:|-----------:|-------------:|------------:|---------:|---------:|----------:|-------------------:|
|            0 |       0.723108 |       0.611286 |  0.723108 |  0.722223 |   0.891215 |         0.78 |        0.8  | 0.590278 | 0.148251 |  0.553564 |                  5 |
|            1 |       0.763913 |       0.664361 |  0.763913 |  0.793984 |   0.882441 |         0.8  |        0.84 | 0.545415 | 0.145368 |  0.531148 |                  5 |
|            2 |       0.807687 |       0.68739  |  0.807687 |  0.838061 |   0.882293 |         0.82 |        0.92 | 0.668321 | 0.141074 |  0.620641 |                  5 |
|            3 |       0.807687 |       0.713365 |  0.807687 |  0.838061 |   0.882293 |         0.82 |        0.92 | 0.668321 | 0.141074 |  0.620641 |                  5 |
|            4 |       0.813702 |       0.734398 |  0.813702 |  0.833906 |   0.88331  |         0.86 |        0.96 | 0.662086 | 0.131358 |  0.611053 |                  5 |
|            5 |       0.813702 |       0.737026 |  0.813702 |  0.833906 |   0.88331  |         0.86 |        0.96 | 0.662086 | 0.131358 |  0.611053 |                  5 |
|            6 |       0.826991 |       0.747496 |  0.826991 |  0.858077 |   0.928764 |         0.84 |        0.92 | 0.66354  | 0.141826 |  0.632141 |                  5 |
|            7 |       0.826991 |       0.75932  |  0.826991 |  0.858077 |   0.928764 |         0.84 |        0.92 | 0.66354  | 0.141826 |  0.632141 |                  5 |
|            8 |       0.827967 |       0.765906 |  0.827967 |  0.85202  |   0.888027 |         0.84 |        0.96 | 0.688871 | 0.124243 |  0.500829 |                  5 |
|            9 |       0.843602 |       0.786799 |  0.843602 |  0.861613 |   0.9244   |         0.86 |        0.96 | 0.722411 | 0.110388 |  0.557426 |                  5 |
|           10 |       0.85673  |       0.792347 |  0.85673  |  0.880042 |   0.941952 |         0.86 |        0.96 | 0.70554  | 0.126139 |  0.438111 |                  5 |
|           11 |       0.85673  |       0.79649  |  0.85673  |  0.880042 |   0.941952 |         0.86 |        0.96 | 0.70554  | 0.126139 |  0.438111 |                  5 |

## Notes
- Helper bundles with a single class are retained in the row table but excluded from binary AUPRC/AUROC fitness.
- The hard router is the discrete expert choice implied by the learned gate; the soft score is the primary champion metric.
- If the next external batch shifts again, freeze the current weights and rerun the GA on the new holdout.