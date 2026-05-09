# BAR-Neo Contextual BMA KR

## 한 줄 결론

이제 global ensemble이 아니라 **상황별 ensemble**이다. CEDAR/TESLA/외부 source, rare-HLA, Korean-HLA, leakage-risk에 따라 method weight가 바뀐다.

## 제일 중요한 점

- public pretrained method는 점수 계산에는 caveated support로 들어갈 수 있지만 clean claim에는 쓰지 않는다.
- `clean_contextual_bma_score`는 public comparator를 제외한 내부/anchor/fallback view다.
- high leakage, rare HLA, low prevalence, external source는 confidence cap과 abstention reason을 만든다.
- QK 계열은 bounded fallback contribution으로만 유지한다.

## 상위 후보

| candidate_id   | source_name   | hla_allele_4digit   |   contextual_bma_score |   clean_contextual_bma_score |   contextual_confidence_score | contextual_abstention_reason_primary            |
|:---------------|:--------------|:--------------------|-----------------------:|-----------------------------:|------------------------------:|:------------------------------------------------|
| CNV0_00733     | CEDAR         | HLA-B*07:02         |               0.938137 |                     0.938137 |                          0.38 | High leakage risk blocks clean contextual claim |
| CNV0_00490     | CEDAR         | HLA-B*07:02         |               0.937412 |                     0.937412 |                          0.38 | High leakage risk blocks clean contextual claim |
| CNV0_00900     | CEDAR         | HLA-B*27:05         |               0.848266 |                     0.848266 |                          0.38 | High leakage risk blocks clean contextual claim |
| CNV0_00346     | CEDAR         | HLA-B*58:01         |               0.83399  |                     0.83399  |                          0.38 | High leakage risk blocks clean contextual claim |
| CNV0_00727     | CEDAR         | HLA-B*27:05         |               0.831236 |                     0.831236 |                          0.38 | High leakage risk blocks clean contextual claim |
| CNV0_00771     | CEDAR         | HLA-A*24:02         |               0.829741 |                     0.829741 |                          0.38 | High leakage risk blocks clean contextual claim |
| CNV0_00328     | CEDAR         | HLA-A*24:02         |               0.828399 |                     0.828399 |                          0.38 | High leakage risk blocks clean contextual claim |
| CNV0_00852     | CEDAR         | HLA-B*27:05         |               0.82526  |                     0.82526  |                          0.38 | High leakage risk blocks clean contextual claim |
| CNV0_00688     | CEDAR         | HLA-A*24:02         |               0.822144 |                     0.822144 |                          0.38 | High leakage risk blocks clean contextual claim |
| CNV0_00753     | CEDAR         | HLA-A*03:01         |               0.818736 |                     0.818736 |                          0.38 | High leakage risk blocks clean contextual claim |
| CNV0_00902     | CEDAR         | HLA-B*27:05         |               0.814439 |                     0.814439 |                          0.38 | High leakage risk blocks clean contextual claim |
| CNV0_00123     | CEDAR         | HLA-A*23:01         |               0.814074 |                     0.814074 |                          0.38 | High leakage risk blocks clean contextual claim |
| CNV0_00071     | CEDAR         | HLA-A*33:03         |               0.812983 |                     0.812983 |                          0.38 | High leakage risk blocks clean contextual claim |
| CNV0_00074     | CEDAR         | HLA-A*33:03         |               0.812584 |                     0.812584 |                          0.38 | High leakage risk blocks clean contextual claim |
| CNV0_00611     | CEDAR         | HLA-A*33:03         |               0.812523 |                     0.812523 |                          0.38 | High leakage risk blocks clean contextual claim |
| CNV0_00604     | CEDAR         | HLA-A*33:03         |               0.811824 |                     0.811824 |                          0.38 | High leakage risk blocks clean contextual claim |
| CNV0_00626     | CEDAR         | HLA-A*33:03         |               0.81157  |                     0.81157  |                          0.38 | High leakage risk blocks clean contextual claim |
| CNV0_00619     | CEDAR         | HLA-A*33:03         |               0.811271 |                     0.811271 |                          0.38 | High leakage risk blocks clean contextual claim |
| CNV0_00728     | CEDAR         | HLA-B*27:05         |               0.811088 |                     0.811088 |                          0.38 | High leakage risk blocks clean contextual claim |
| CNV0_00649     | CEDAR         | HLA-A*33:03         |               0.811034 |                     0.811034 |                          0.38 | High leakage risk blocks clean contextual claim |
| CNV0_00216     | CEDAR         | HLA-A*33:03         |               0.810679 |                     0.810679 |                          0.38 | High leakage risk blocks clean contextual claim |
| CNV0_00633     | CEDAR         | HLA-A*33:03         |               0.81052  |                     0.81052  |                          0.38 | High leakage risk blocks clean contextual claim |
| CNV0_00812     | CEDAR         | HLA-B*44:03         |               0.809775 |                     0.809775 |                          0.38 | High leakage risk blocks clean contextual claim |
| CNV0_00630     | CEDAR         | HLA-A*33:03         |               0.80969  |                     0.80969  |                          0.38 | High leakage risk blocks clean contextual claim |
| CNV0_00310     | CEDAR         | HLA-B*40:01         |               0.809236 |                     0.809236 |                          0.38 | High leakage risk blocks clean contextual claim |

## context 요약

| context_label                |   n_candidates |   n_pos |   positive_prevalence |   mean_contextual_bma_score |   mean_clean_contextual_bma_score |   mean_confidence |   abstention_rate |   clean_claim_allowed |   caveated_public_used_rate |
|:-----------------------------|---------------:|--------:|----------------------:|----------------------------:|----------------------------------:|------------------:|------------------:|----------------------:|----------------------------:|
| global                       |           2715 |    1179 |            0.434254   |                   0.397749  |                          0.426496 |         0.385487  |                 1 |                     0 |                           0 |
| high_leakage_review_only     |           2304 |    1143 |            0.496094   |                   0.423404  |                          0.459935 |         0.349818  |                 1 |                     0 |                           0 |
| low_prevalence               |           1035 |      48 |            0.0463768  |                   0.20607   |                          0.230327 |         0.391311  |                 1 |                     0 |                           0 |
| korean_hla_focus             |           1031 |     377 |            0.365664   |                   0.376787  |                          0.413263 |         0.374453  |                 1 |                     0 |                           0 |
| source:CEDAR                 |            909 |     851 |            0.936194   |                   0.724647  |                          0.724647 |         0.38      |                 1 |                     0 |                           0 |
| hla:HLA-A*02:01              |            666 |     232 |            0.348348   |                   0.372288  |                          0.408474 |         0.375843  |                 1 |                     0 |                           0 |
| source:TESLA_mmc4            |            605 |      37 |            0.061157   |                   0.242706  |                          0.242706 |         0.38      |                 1 |                     0 |                           0 |
| source:NEPdb                 |            572 |     151 |            0.263986   |                   0.270652  |                          0.270652 |         0.38      |                 1 |                     0 |                           0 |
| rare_hla                     |            522 |     421 |            0.806513   |                   0.564311  |                          0.578724 |         0.37069   |                 1 |                     0 |                           0 |
| external_or_holdout          |            319 |     136 |            0.426332   |                   0.180188  |                          0.422646 |         0.26416   |                 1 |                     0 |                           0 |
| source:TESLA_mmc7_validation |            310 |       4 |            0.0129032  |                   0.200171  |                          0.200171 |         0.547258  |                 1 |                     0 |                           0 |
| hla:HLA-A*01:01              |            260 |      35 |            0.134615   |                   0.242765  |                          0.280528 |         0.366815  |                 1 |                     0 |                           0 |
| source:ITSNdb_main           |            199 |     129 |            0.648241   |                   0.266771  |                          0.4247   |         0.39603   |                 1 |                     0 |                           0 |
| hla:HLA-A*03:01              |            148 |      52 |            0.351351   |                   0.375727  |                          0.386163 |         0.405137  |                 1 |                     0 |                           0 |
| hla:HLA-A*11:01              |            144 |      40 |            0.277778   |                   0.322937  |                          0.390781 |         0.317522  |                 1 |                     0 |                           0 |
| low_support_hla              |            132 |      83 |            0.628788   |                   0.415308  |                          0.438565 |         0.364275  |                 1 |                     0 |                           0 |
| source:ITSNdb_Val            |            120 |       7 |            0.0583333  |                   0.0366036 |                          0.399312 |         0.0454739 |                 1 |                     0 |                           0 |
| hla:HLA-B*08:01              |            113 |       8 |            0.0707965  |                   0.220285  |                          0.222252 |         0.460081  |                 1 |                     0 |                           0 |
| hla:HLA-A*23:01              |            104 |      15 |            0.144231   |                   0.303332  |                          0.303332 |         0.520577  |                 1 |                     0 |                           0 |
| hla:HLA-A*26:01              |            102 |       1 |            0.00980392 |                   0.155802  |                          0.214754 |         0.315686  |                 1 |                     0 |                           0 |
| hla:HLA-A*33:03              |             97 |      97 |            1          |                   0.722658  |                          0.722658 |         0.38      |                 1 |                     0 |                           0 |
| hla:HLA-B*44:02              |             95 |      25 |            0.263158   |                   0.290608  |                          0.290608 |         0.426526  |                 1 |                     0 |                           0 |
| hla:HLA-B*27:05              |             84 |      46 |            0.547619   |                   0.508909  |                          0.521322 |         0.442412  |                 1 |                     0 |                           0 |
| hla:HLA-A*24:02              |             74 |      36 |            0.486486   |                   0.441255  |                          0.4599   |         0.378874  |                 1 |                     0 |                           0 |
| hla:HLA-B*15:01              |             64 |      12 |            0.1875     |                   0.245877  |                          0.253808 |         0.487656  |                 1 |                     0 |                           0 |
| hla:HLA-B*07:02              |             60 |      29 |            0.483333   |                   0.37904   |                          0.39899  |         0.381787  |                 1 |                     0 |                           0 |
| hla:HLA-B*40:01              |             50 |      47 |            0.94       |                   0.668138  |                          0.681774 |         0.3724    |                 1 |                     0 |                           0 |
| hla:HLA-A*68:01              |             41 |      12 |            0.292683   |                   0.284142  |                          0.291246 |         0.372878  |                 1 |                     0 |                           0 |
| hla:HLA-B*38:01              |             26 |      20 |            0.769231   |                   0.472897  |                          0.472897 |         0.38      |                 1 |                     0 |                           0 |
| hla:HLA-B*44:03              |             22 |      19 |            0.863636   |                   0.545209  |                          0.571172 |         0.379265  |                 1 |                     0 |                           0 |
| hla:HLA-C*05:01              |             22 |      11 |            0.5        |                   0.398707  |                          0.417693 |         0.362727  |                 1 |                     0 |                           0 |
| hla:HLA-B*35:01              |             21 |      21 |            1          |                   0.481395  |                          0.594664 |         0.313927  |                 1 |                     0 |                           0 |
| hla:HLA-C*07:02              |             19 |       1 |            0.0526316  |                   0.327214  |                          0.327214 |         0.38      |                 1 |                     0 |                           0 |
| hla:HLA-B*49:01              |             18 |      18 |            1          |                   0.676204  |                          0.676204 |         0.38      |                 1 |                     0 |                           0 |
| hla:HLA-A*68:02              |             17 |      17 |            1          |                   0.684055  |                          0.684055 |         0.38      |                 1 |                     0 |                           0 |
| hla:HLA-C*06:02              |             17 |       8 |            0.470588   |                   0.3766    |                          0.400138 |         0.357647  |                 1 |                     0 |                           0 |
| hla:HLA-C*12:03              |             17 |       9 |            0.529412   |                   0.43016   |                          0.43016  |         0.381176  |                 1 |                     0 |                           0 |
| hla:HLA-B*45:01              |             16 |      16 |            1          |                   0.705478  |                          0.705478 |         0.38      |                 1 |                     0 |                           0 |
| hla:HLA-B*57:01              |             15 |       8 |            0.533333   |                   0.499175  |                          0.499175 |         0.38      |                 1 |                     0 |                           0 |
| hla:HLA-B*51:01              |             14 |      14 |            1          |                   0.636285  |                          0.68523  |         0.352857  |                 1 |                     0 |                           0 |
