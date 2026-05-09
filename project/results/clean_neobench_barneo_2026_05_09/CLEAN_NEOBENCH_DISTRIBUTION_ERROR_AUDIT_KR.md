# CLEAN-NeoBench Distribution Error Audit KR

## 한 줄 결론

지금 싸움은 단순 AUPRC 1등 싸움이 아니다. 실제로는 **source prevalence shift**, **rare/low-support HLA**, **high leakage-risk region**, **low-prevalence TESLA-like setting**에서 누가 망가지는지 보는 싸움이다.

## method별 취약성

| method_name                                       |   distribution_vulnerability_score | primary_distribution_issues                                                                                     |   low_prevalence_high_ranked_negative_rate |   rare_hla_missed_positive_rate |   external_missed_positive_rate |
|:--------------------------------------------------|-----------------------------------:|:----------------------------------------------------------------------------------------------------------------|-------------------------------------------:|--------------------------------:|--------------------------------:|
| NetMHCpan_4.1                                     |                           0.426447 | source-prior coupling; rare-HLA instability; leakage-region false-positive pressure; external/holdout fragility |                                   0        |                       0.789474  |                       0.419355  |
| PRIME                                             |                           0.376537 | source-prior coupling; rare-HLA instability; external/holdout fragility                                         |                                   0        |                       0.526316  |                       0.435484  |
| pu_weighted_rf_train_prior_calibrated             |                           0.363811 | source-prior coupling                                                                                           |                                   0        |                       0.0522388 |                     nan         |
| source_balanced_plus_pu_rf_train_prior_calibrated |                           0.363488 | source-prior coupling                                                                                           |                                   0        |                       0.0472637 |                     nan         |
| source_balanced_rf_train_prior_calibrated         |                           0.36002  | source-prior coupling                                                                                           |                                   0        |                       0.0472637 |                     nan         |
| MHCflurry                                         |                           0.313558 | source-prior coupling; rare-HLA instability; external/holdout fragility                                         |                                   0        |                       0.684211  |                       0.33871   |
| BigMHC_IM                                         |                           0.312275 | source-prior coupling; rare-HLA instability; external/holdout fragility                                         |                                   0        |                       0.526316  |                       0.322581  |
| sourceheld_counterfactual_rf                      |                           0.308387 | source-prior coupling                                                                                           |                                   0.020595 |                       0.20398   |                     nan         |
| W7A_QK_only                                       |                           0.305821 | source-prior coupling; rare-HLA instability; external/holdout fragility                                         |                                   0        |                       0.473684  |                       0.33871   |
| BAR_Neo_BMA                                       |                           0.302075 | source-prior coupling                                                                                           |                                   0        |                       0.0332542 |                       0.0514706 |
| Structure_LR                                      |                           0.277101 | source-prior coupling; external/holdout fragility                                                               |                                   0        |                       0.230769  |                       0.440678  |
| W7A_full                                          |                           0.261161 | source-prior coupling; external/holdout fragility                                                               |                                   0        |                       0.210526  |                       0.33871   |
| ESM2_Bayesian                                     |                           0.258438 | source-prior coupling; external/holdout fragility                                                               |                                   0        |                       0         |                       0.440678  |
| W7B_stacked                                       |                           0.248239 | source-prior coupling; external/holdout fragility                                                               |                                   0        |                       0.105263  |                       0.483871  |

## 분포 shift와 error가 같이 보이는 slice

| feature                         | value               |   train_n |   external_n |   distribution_shift_score |   max_error_pressure_score | diagnosis                                                                                                                    |
|:--------------------------------|:--------------------|----------:|-------------:|---------------------------:|---------------------------:|:-----------------------------------------------------------------------------------------------------------------------------|
| distribution_partition          | external_or_holdout |         0 |          319 |                   1.7993   |                   0.405485 | external/heldout-only slice; no local train support; large train-vs-external representation shift; enriched missed positives |
| hla_allele_4digit               | HLA-C*07:01         |         5 |            1 |                   1.10677  |                   1        | large positive-prevalence shift; enriched missed positives; rare allele slice                                                |
| hla_allele_4digit               | HLA-A*25:01         |         2 |            1 |                   1.02207  |                   1        | large positive-prevalence shift; enriched high-ranked negatives; enriched missed positives; rare allele slice                |
| hla_allele_4digit               | HLA-B*35:03         |         2 |            1 |                   1.02207  |                   1        | large positive-prevalence shift; enriched high-ranked negatives; enriched missed positives; rare allele slice                |
| hla_allele_4digit               | HLA-B*56:01         |         2 |            1 |                   1.02207  |                   1        | large positive-prevalence shift; enriched missed positives; rare allele slice                                                |
| hla_supertype                   | A25                 |         2 |            1 |                   1.02207  |                   1        | large positive-prevalence shift; enriched high-ranked negatives; enriched missed positives                                   |
| hla_supertype                   | B56                 |         2 |            1 |                   1.02207  |                   1        | large positive-prevalence shift; enriched missed positives                                                                   |
| source_name                     | ITSNdb_main         |         0 |          199 |                   1.42312  |                   0.423835 | external/heldout-only slice; no local train support; large train-vs-external representation shift; enriched missed positives |
| hla_allele_4digit               | HLA-A*68:01         |        37 |            4 |                   0.98483  |                   1        | large positive-prevalence shift; enriched missed positives                                                                   |
| exact_peptide_hla_train_overlap | False               |       305 |          136 |                   1.23508  |                   0.581954 | large train-vs-external representation shift; large positive-prevalence shift; enriched missed positives                     |
| hla_supertype                   | C07                 |        25 |            1 |                   0.974867 |                   1        | large positive-prevalence shift; enriched missed positives                                                                   |
| hla_allele_4digit               | HLA-B*15:01         |        61 |            3 |                   0.973303 |                   1        | large positive-prevalence shift; enriched missed positives                                                                   |
| near_peptide_train_overlap      | False               |       294 |          121 |                   1.1623   |                   0.584192 | large train-vs-external representation shift; large positive-prevalence shift; enriched missed positives                     |
| leakage_risk_level              | medium              |        11 |           12 |                   1.1892   |                   0.514286 | large positive-prevalence shift; enriched high-ranked negatives; enriched missed positives                                   |
| source_name                     | ITSNdb_Val          |         0 |          120 |                   1.17547  |                   0.392857 | external/heldout-only slice; no local train support; large train-vs-external representation shift; enriched missed positives |
| hla_allele_4digit               | HLA-A*02:17         |         0 |            1 |                   0.802433 |                   1        | external/heldout-only slice; no local train support; enriched missed positives; rare allele slice                            |
| hla_allele_4digit               | HLA-B*39:06         |         0 |            1 |                   0.802433 |                   1        | external/heldout-only slice; no local train support; enriched missed positives; rare allele slice                            |
| hla_allele_4digit               | HLA-B*41:02         |         0 |            1 |                   0.802433 |                   1        | external/heldout-only slice; no local train support; enriched missed positives; rare allele slice                            |
| hla_supertype                   | B39                 |         1 |            1 |                   0.716503 |                   1        | enriched missed positives                                                                                                    |
| hla_supertype                   | B15                 |        86 |            3 |                   0.708709 |                   1        | large positive-prevalence shift; enriched missed positives                                                                   |
| hla_allele_4digit               | HLA-B*27:05         |        63 |           21 |                   0.933663 |                   0.492647 | large positive-prevalence shift; enriched high-ranked negatives; enriched missed positives                                   |
| distribution_partition          | local_train_pool    |      2396 |            0 |                   1        |                   0.389877 | train-only slice; external robustness untested; large train-vs-external representation shift; enriched missed positives      |
| hla_supertype                   | A68                 |        54 |            4 |                   0.694487 |                   1        | large positive-prevalence shift; enriched missed positives                                                                   |
| hla_allele_4digit               | HLA-C*06:02         |        16 |            1 |                   0.694336 |                   1        | large positive-prevalence shift; enriched missed positives; rare allele slice                                                |
| hla_supertype                   | C06                 |        16 |            1 |                   0.694336 |                   1        | large positive-prevalence shift; enriched missed positives                                                                   |
| leakage_risk_level              | low                 |       294 |           94 |                   0.867704 |                   0.584192 | enriched missed positives                                                                                                    |
| hla_supertype                   | B27                 |        77 |           21 |                   0.890086 |                   0.492647 | large positive-prevalence shift; enriched high-ranked negatives; enriched missed positives                                   |
| hla_allele_4digit               | HLA-C*12:03         |        16 |            1 |                   0.631836 |                   1        | large positive-prevalence shift; enriched missed positives; rare allele slice                                                |
| hla_allele_4digit               | HLA-C*05:01         |        21 |            1 |                   0.631362 |                   1        | large positive-prevalence shift; enriched missed positives                                                                   |
| hla_supertype                   | C05                 |        21 |            1 |                   0.631362 |                   1        | large positive-prevalence shift; enriched missed positives                                                                   |

## 실용 판단

- aggregate leaderboard에서 이긴 method라도 source/HLA slice에서 무너지면 clean claim 금지.
- BAR-Neo-BMA는 이 취약성을 weight와 abstention으로 흡수하는 controller 위치가 맞다.
- QK는 feature/selector/fallback로만 사용하고 quantum advantage 표현은 금지.
- public pretrained tool은 training-corpus overlap audit 전까지 caveated comparator다.
- PAAD/THCA patient gate는 patient metadata 들어오기 전까지 demo/triage-only다.

## 다음 행동

1. public training corpus를 넣고 row-level overlap audit.
2. source-heldout/HLA-heldout challenge split을 고정.
3. rare-HLA/low-prevalence slice의 manual review queue를 우선 검토.
4. patient metadata가 있는 PAAD/THCA demo 1-2개를 붙여 실제 gate를 검증.
