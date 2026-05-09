# CLEAN-NeoBench Failure Mode Analysis KR

Date: 2026-05-09

## 한 줄 결론

틀리는 후보들은 무작위가 아니다. 가장 큰 축은 **source prevalence shift**, 두 번째는 **HLA allele support imbalance**, 세 번째는 **high leakage-risk region**이다.

## 학습/외부 분포 판단

| distribution_partition   | source_name           |   count |   sum |      mean |
|:-------------------------|:----------------------|--------:|------:|----------:|
| external_or_holdout      | ITSNdb_main           |     199 |   129 | 0.648241  |
| external_or_holdout      | ITSNdb_Val            |     120 |     7 | 0.0583333 |
| local_train_pool         | CEDAR                 |     909 |   851 | 0.936194  |
| local_train_pool         | TESLA_mmc4            |     605 |    37 | 0.061157  |
| local_train_pool         | NEPdb                 |     572 |   151 | 0.263986  |
| local_train_pool         | TESLA_mmc7_validation |     310 |     4 | 0.0129032 |

CEDAR는 거의 positive-rich source이고, TESLA/ITSNdb_Val은 low-prevalence source다. 이 차이가 너무 커서 model이 biological rule이 아니라 source prior를 배울 수 있다.

## method별 failure 요약

| method_name                                       |   n_scored |   high_ranked_negative_top10pct |   missed_positive_bottom50pct |   missed_positive_rate_among_pos |   median_score_pos |   median_score_neg |
|:--------------------------------------------------|-----------:|--------------------------------:|------------------------------:|---------------------------------:|-------------------:|-------------------:|
| NetMHCpan_4.1                                     |        319 |                              23 |                            66 |                        0.485294  |           0.367347 |          0.333333  |
| PRIME                                             |        319 |                              16 |                            56 |                        0.411765  |           0.430387 |          0.351852  |
| Structure_LR                                      |        311 |                               5 |                            50 |                        0.384615  |           0.445418 |          0.390244  |
| MHCflurry                                         |        319 |                              10 |                            52 |                        0.382353  |           0.39403  |          0.267974  |
| BigMHC_IM                                         |        319 |                              15 |                            40 |                        0.294118  |           0.473831 |          0.301587  |
| W7B_stacked                                       |        319 |                               5 |                            40 |                        0.294118  |           0.6      |          0.270833  |
| W7A_QK_only                                       |        319 |                               5 |                            37 |                        0.272059  |           0.5      |          0.32      |
| W7A_full                                          |        319 |                               6 |                            35 |                        0.257353  |           0.545455 |          0.285714  |
| ESM2_Bayesian                                     |        311 |                               9 |                            33 |                        0.253846  |           0.512571 |          0.288462  |
| sourceheld_counterfactual_rf                      |       2396 |                              47 |                           238 |                        0.228188  |           0.607143 |          0.319008  |
| pu_weighted_rf_train_prior_calibrated             |       2396 |                              18 |                           129 |                        0.123682  |           0.93186  |          0.126506  |
| source_balanced_plus_pu_rf_train_prior_calibrated |       2396 |                              17 |                           127 |                        0.121764  |           0.93     |          0.107822  |
| source_balanced_rf_train_prior_calibrated         |       2396 |                              15 |                           122 |                        0.11697   |           0.918605 |          0.112082  |
| BAR_Neo_BMA                                       |       2715 |                              19 |                           107 |                        0.0907549 |           0.379384 |          0.0594981 |

해석:

- `BAR_Neo_BMA`는 missed positive rate가 가장 낮지만, abstention을 강하게 유지한다.
- `source_balanced_rf_train_prior_calibrated` 계열은 현재 best internal group이다.
- `Structure_LR`는 honest anchor이고 false-positive pressure가 낮다.
- `W7A_QK_only`는 점수는 좋지만 bounded fallback으로만 유지한다.
- public pretrained methods는 overlap unresolved라 clean baseline이 아니다.

## train vs external shift

| feature                         | value                 |   train_n |   train_fraction |   train_positive_prevalence |   external_n |   external_fraction |   external_positive_prevalence |   absolute_fraction_delta |
|:--------------------------------|:----------------------|----------:|-----------------:|----------------------------:|-------------:|--------------------:|-------------------------------:|--------------------------:|
| distribution_partition          | external_or_holdout   |         0 |        0         |                nan          |          319 |           1         |                      0.426332  |                 1         |
| distribution_partition          | local_train_pool      |      2396 |        1         |                  0.435309   |            0 |           0         |                    nan         |                 1         |
| source_name                     | ITSNdb_main           |         0 |        0         |                nan          |          199 |           0.623824  |                      0.648241  |                 0.623824  |
| source_name                     | CEDAR                 |       909 |        0.379382  |                  0.936194   |            0 |           0         |                    nan         |                 0.379382  |
| source_name                     | ITSNdb_Val            |         0 |        0         |                nan          |          120 |           0.376176  |                      0.0583333 |                 0.376176  |
| peptide_length                  | 9                     |      1174 |        0.489983  |                  0.491482   |          270 |           0.846395  |                      0.348148  |                 0.356412  |
| exact_peptide_hla_train_overlap | False                 |       305 |        0.127295  |                  0.00983607 |          136 |           0.426332  |                      0.455882  |                 0.299037  |
| exact_peptide_hla_train_overlap | True                  |      2091 |        0.872705  |                  0.49737    |          183 |           0.573668  |                      0.404372  |                 0.299037  |
| near_peptide_train_overlap      | True                  |      2102 |        0.877295  |                  0.494767   |          198 |           0.62069   |                      0.414141  |                 0.256606  |
| near_peptide_train_overlap      | False                 |       294 |        0.122705  |                  0.0102041  |          121 |           0.37931   |                      0.446281  |                 0.256606  |
| source_name                     | TESLA_mmc4            |       605 |        0.252504  |                  0.061157   |            0 |           0         |                    nan         |                 0.252504  |
| source_name                     | NEPdb                 |       572 |        0.238731  |                  0.263986   |            0 |           0         |                    nan         |                 0.238731  |
| leakage_risk_level              | high                  |      2091 |        0.872705  |                  0.49737    |          213 |           0.667712  |                      0.483568  |                 0.204993  |
| leakage_risk_level              | low                   |       294 |        0.122705  |                  0.0102041  |           94 |           0.294671  |                      0.297872  |                 0.171966  |
| source_name                     | TESLA_mmc7_validation |       310 |        0.129382  |                  0.0129032  |            0 |           0         |                    nan         |                 0.129382  |
| peptide_length                  | 10                    |       674 |        0.281302  |                  0.339763   |           49 |           0.153605  |                      0.857143  |                 0.127697  |
| hla_allele_4digit               | HLA-A*02:01           |       555 |        0.231636  |                  0.313514   |          111 |           0.347962  |                      0.522523  |                 0.116326  |
| hla_supertype                   | A02                   |       588 |        0.245409  |                  0.352041   |          113 |           0.354232  |                      0.530973  |                 0.108823  |
| peptide_length                  | 11                    |       222 |        0.0926544 |                  0.490991   |            0 |           0         |                    nan         |                 0.0926544 |
| hla_allele_4digit               | HLA-A*26:01           |        72 |        0.0300501 |                  0.0138889  |           30 |           0.0940439 |                      0         |                 0.0639938 |
| hla_supertype                   | A26                   |        72 |        0.0300501 |                  0.0138889  |           30 |           0.0940439 |                      0         |                 0.0639938 |
| hla_supertype                   | A33                   |       109 |        0.0454925 |                  1          |            0 |           0         |                    nan         |                 0.0454925 |
| hla_allele_4digit               | HLA-A*23:01           |       104 |        0.0434057 |                  0.144231   |            0 |           0         |                    nan         |                 0.0434057 |
| hla_supertype                   | A23                   |       104 |        0.0434057 |                  0.144231   |            0 |           0         |                    nan         |                 0.0434057 |
| hla_allele_4digit               | HLA-A*01:01           |       218 |        0.090985  |                  0.133028   |           42 |           0.131661  |                      0.142857  |                 0.0406765 |

## 실용 조치

1. source prior correction을 더 강하게 넣는다.
2. low-prevalence source에서는 top-k threshold를 보수화한다.
3. rare HLA allele은 confidence를 낮춘다.
4. HLA별 minimum support count를 reliability feature로 강화한다.
5. public tool agreement는 caveated support로만 쓰고 clean feature로 쓰지 않는다.
6. PAAD/THCA patient metadata가 들어오기 전까지 clinical-style confidence는 내지 않는다.

## 파일

- Full report: `CLEAN_NEOBENCH_FAILURE_MODE_ANALYSIS.md`
- Failure cases: `clean_neobench_failure_cases.tsv`
- Failure summary: `clean_neobench_failure_summary.tsv`
- Train/external shift: `clean_neobench_train_external_distribution_shift.tsv`
- Failure enrichment: `clean_neobench_failure_distribution_enrichment.tsv`
