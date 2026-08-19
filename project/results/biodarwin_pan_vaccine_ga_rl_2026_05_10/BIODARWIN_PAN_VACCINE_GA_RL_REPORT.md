# BioDarwin Pan-Vaccine GA/RL algorithm

Generated: 2026-05-11T00:51:57

## Claim boundary

- This is the algorithm artifact, not the public-source scout.
- GA/RL fitness used academic training sources only; validation-like and industrial public locked rows are reported after freezing.
- Class I is scored directly; Class II is represented as a helper axis and scout lane until curated Class-II labels are added.
- Industrial/public patent rows remain do-not-train and manual-QA required.

## Outputs

- Metrics: `/home/seungho/personal/THCA_data_analysis/project/results/biodarwin_pan_vaccine_ga_rl_2026_05_10/biodarwin_pan_vaccine_metrics.tsv`
- Scores: `/home/seungho/personal/THCA_data_analysis/project/results/biodarwin_pan_vaccine_ga_rl_2026_05_10/biodarwin_academic_candidate_scores.tsv`
- Industrial scores: `/home/seungho/personal/THCA_data_analysis/project/results/biodarwin_pan_vaccine_ga_rl_2026_05_10/biodarwin_industrial_locked_scores.tsv`
- Class-II scout: `/home/seungho/personal/THCA_data_analysis/project/results/biodarwin_pan_vaccine_ga_rl_2026_05_10/biodarwin_classII_scout_scores.tsv`

## Frozen metrics

| split                    | algorithm                 |    n |   positives |   positive_rate |    AUPRC |    AUROC |   top5_precision |   top5_hits |   top10_precision |   top10_hits |   top24_precision |   top24_hits |   top96_precision |   top96_hits |
|:-------------------------|:--------------------------|-----:|------------:|----------------:|---------:|---------:|-----------------:|------------:|------------------:|-------------:|------------------:|-------------:|------------------:|-------------:|
| academic_train_sources   | BigMHC_IM                 | 2285 |        1168 |       0.51116   | 0.698129 | 0.68254  |              0.8 |           4 |               0.9 |            9 |          0.791667 |           19 |         0.875     |           84 |
| academic_train_sources   | BigMHC_EL                 | 2285 |        1168 |       0.51116   | 0.700143 | 0.713472 |              0.6 |           3 |               0.7 |            7 |          0.666667 |           16 |         0.760417  |           73 |
| academic_train_sources   | CROSS_core                | 2285 |        1168 |       0.51116   | 0.882405 | 0.862176 |              1   |           5 |               1   |           10 |          0.916667 |           22 |         0.947917  |           91 |
| academic_train_sources   | CROSS_stress              | 2285 |        1168 |       0.51116   | 0.882405 | 0.862176 |              1   |           5 |               1   |           10 |          0.916667 |           22 |         0.947917  |           91 |
| academic_train_sources   | CROSS_BMA                 | 2285 |        1168 |       0.51116   | 0.872902 | 0.848838 |              1   |           5 |               1   |           10 |          0.916667 |           22 |         0.927083  |           89 |
| academic_train_sources   | CROSS_finetuned           | 2285 |        1168 |       0.51116   | 0.855428 | 0.853874 |              1   |           5 |               1   |           10 |          0.916667 |           22 |         0.875     |           84 |
| academic_train_sources   | CROSS_integrated          | 2285 |        1168 |       0.51116   | 0.884649 | 0.87363  |              1   |           5 |               1   |           10 |          1        |           24 |         1         |           96 |
| academic_train_sources   | CROSS_claimsafe           | 2285 |        1168 |       0.51116   | 0.867512 | 0.838301 |              1   |           5 |               1   |           10 |          1        |           24 |         1         |           96 |
| academic_train_sources   | BioDarwin_PublicFallback  | 2285 |        1168 |       0.51116   | 0.734333 | 0.704165 |              1   |           5 |               1   |           10 |          0.958333 |           23 |         0.90625   |           87 |
| academic_train_sources   | BioDarwin_PanVax          | 2285 |        1168 |       0.51116   | 0.917214 | 0.898229 |              1   |           5 |               1   |           10 |          1        |           24 |         1         |           96 |
| academic_validation_like | BigMHC_IM                 |  430 |          11 |       0.0255814 | 0.333705 | 0.827728 |              0.6 |           3 |               0.4 |            4 |          0.25     |            6 |         0.09375   |            9 |
| academic_validation_like | BigMHC_EL                 |  430 |          11 |       0.0255814 | 0.30759  | 0.714689 |              0.4 |           2 |               0.3 |            3 |          0.208333 |            5 |         0.0833333 |            8 |
| academic_validation_like | CROSS_core                |  430 |          11 |       0.0255814 | 0.650789 | 0.87698  |              1   |           5 |               0.7 |            7 |          0.291667 |            7 |         0.0729167 |            7 |
| academic_validation_like | CROSS_stress              |  430 |          11 |       0.0255814 | 0.650789 | 0.87698  |              1   |           5 |               0.7 |            7 |          0.291667 |            7 |         0.0729167 |            7 |
| academic_validation_like | CROSS_BMA                 |  430 |          11 |       0.0255814 | 0.659761 | 0.867433 |              1   |           5 |               0.7 |            7 |          0.291667 |            7 |         0.0729167 |            7 |
| academic_validation_like | CROSS_finetuned           |  430 |          11 |       0.0255814 | 0.605625 | 0.879149 |              1   |           5 |               0.6 |            6 |          0.291667 |            7 |         0.09375   |            9 |
| academic_validation_like | CROSS_integrated          |  430 |          11 |       0.0255814 | 0.713841 | 0.946192 |              1   |           5 |               0.7 |            7 |          0.333333 |            8 |         0.09375   |            9 |
| academic_validation_like | CROSS_claimsafe           |  430 |          11 |       0.0255814 | 0.84648  | 0.991972 |              1   |           5 |               0.8 |            8 |          0.416667 |           10 |         0.114583  |           11 |
| academic_validation_like | BioDarwin_PublicFallback  |  430 |          11 |       0.0255814 | 0.539681 | 0.879366 |              0.8 |           4 |               0.6 |            6 |          0.291667 |            7 |         0.0833333 |            8 |
| academic_validation_like | BioDarwin_PanVax          |  430 |          11 |       0.0255814 | 0.886307 | 0.972445 |              1   |           5 |               0.9 |            9 |          0.416667 |           10 |         0.104167  |           10 |
| industrial_locked_v0     | BigMHC_IM                 |   25 |           9 |       0.36      | 0.578121 | 0.611111 |              0.8 |           4 |               0.4 |            4 |          0.375    |            9 |         0.36      |            9 |
| industrial_locked_v0     | BigMHC_EL                 |   25 |           9 |       0.36      | 0.416128 | 0.451389 |              0.2 |           1 |               0.2 |            2 |          0.375    |            9 |         0.36      |            9 |
| industrial_locked_v0     | GA_public_feature_adapter |   25 |           9 |       0.36      | 0.515396 | 0.590278 |              0.6 |           3 |               0.4 |            4 |          0.375    |            9 |         0.36      |            9 |
| industrial_locked_v0     | BioDarwin_PublicFallback  |   25 |           9 |       0.36      | 0.494322 | 0.541667 |              0.6 |           3 |               0.4 |            4 |          0.375    |            9 |         0.36      |            9 |
| industrial_locked_v0     | BioDarwin_PanVax          |   25 |           9 |       0.36      | 0.494322 | 0.541667 |              0.6 |           3 |               0.4 |            4 |          0.375    |            9 |         0.36      |            9 |
