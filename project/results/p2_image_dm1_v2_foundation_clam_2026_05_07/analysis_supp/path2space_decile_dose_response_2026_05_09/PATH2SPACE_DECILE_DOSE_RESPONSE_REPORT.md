# Path2Space predicted-decile dose-response controls

## Verdict

- GSE250521 raw top-minus-bottom observed DM1/RAI: **1.479**; decile Spearman **0.407**.
- GSE250521 coord+QC residual top-minus-bottom: **0.128**; decile Spearman **0.206**.
- GSE230424 raw top-minus-bottom observed DM1/low-RAI: **0.699**; decile Spearman **0.607**.
- GSE230424 coord+QC residual top-minus-bottom: **0.159**; decile Spearman **0.222**.

## Interpretation

Predicted-decile dose response checks whether higher predicted DM1/RAI bins carry progressively higher observed DM1/RAI within each slide/sample. This complements top-decile hotspot overlap by avoiding a single threshold. Positive top-minus-bottom separation supports ranking calibration; weaker residual slopes preserve the QC/smoothness caveat.

## Overall Summary

| evidence                                    |     n |   decile_obs_spearman |   decile_obs_p |   top_minus_bottom |   monotonic_increases |   top_decile_mean_obs |   bottom_decile_mean_obs |   group |
|:--------------------------------------------|------:|----------------------:|---------------:|-------------------:|----------------------:|----------------------:|-------------------------:|--------:|
| GSE250521 UNI raw DM1/RAI                   |  3200 |                0.4071 |              0 |             1.4793 |                     9 |                0.913  |                  -0.5663 |     nan |
| GSE250521 UNI coord+QC residual DM1/RAI     |  3200 |                0.2063 |              0 |             0.1281 |                     8 |                0.073  |                  -0.0551 |     nan |
| GSE230424 H&E raw DM1/low-RAI               | 15489 |                0.6069 |              0 |             0.6986 |                     9 |                0.3877 |                  -0.3109 |     nan |
| GSE230424 H&E coord+QC residual DM1/low-RAI | 15489 |                0.2224 |              0 |             0.1586 |                     8 |                0.0724 |                  -0.0861 |     nan |

## Group-Level Summary

| evidence                                    |    n |   decile_obs_spearman |   decile_obs_p |   top_minus_bottom |   monotonic_increases |   top_decile_mean_obs |   bottom_decile_mean_obs | group             |
|:--------------------------------------------|-----:|----------------------:|---------------:|-------------------:|----------------------:|----------------------:|-------------------------:|:------------------|
| GSE250521 UNI raw DM1/RAI                   |  200 |                0.55   |         0      |             2.4099 |                     8 |                1.5057 |                  -0.9042 | GSM7980860_N-1    |
| GSE250521 UNI raw DM1/RAI                   |  200 |                0.4011 |         0      |             1.8142 |                     7 |                1.2733 |                  -0.5408 | GSM7980861_N-2    |
| GSE250521 UNI raw DM1/RAI                   |  200 |                0.2626 |         0.0002 |             0.7212 |                     5 |                0.4069 |                  -0.3143 | GSM7980862_N-3    |
| GSE250521 UNI raw DM1/RAI                   |  200 |                0.313  |         0      |             1.7422 |                     6 |                1.2335 |                  -0.5087 | GSM7980863_N-4    |
| GSE250521 UNI raw DM1/RAI                   |  200 |                0.4654 |         0      |             2.5285 |                     6 |                0.577  |                  -1.9515 | GSM7980864_PTC-1  |
| GSE250521 UNI raw DM1/RAI                   |  200 |                0.6089 |         0      |             2.236  |                     7 |                1.659  |                  -0.577  | GSM7980865_PTC-2  |
| GSE250521 UNI raw DM1/RAI                   |  200 |                0.1844 |         0.0089 |             0.6273 |                     3 |                0.4967 |                  -0.1306 | GSM7980866_PTC-3  |
| GSE250521 UNI raw DM1/RAI                   |  200 |                0.3492 |         0      |             1.3577 |                     5 |                1.1755 |                  -0.1821 | GSM7980867_PTC-4  |
| GSE250521 UNI raw DM1/RAI                   |  200 |                0.4854 |         0      |             1.9098 |                     7 |                1.4595 |                  -0.4503 | GSM7980868_LPTC-1 |
| GSE250521 UNI raw DM1/RAI                   |  200 |                0.655  |         0      |             2.0703 |                     6 |                0.5899 |                  -1.4803 | GSM7980869_LPTC-2 |
| GSE250521 UNI raw DM1/RAI                   |  200 |                0.5261 |         0      |             2.1489 |                     5 |                0.9531 |                  -1.1959 | GSM7980870_LPTC-3 |
| GSE250521 UNI raw DM1/RAI                   |  200 |                0.3175 |         0      |             1.591  |                     5 |                1.4185 |                  -0.1725 | GSM7980871_LPTC-4 |
| GSE250521 UNI raw DM1/RAI                   |  200 |                0.0302 |         0.6716 |             0.0252 |                     5 |               -0.0145 |                  -0.0396 | GSM7980872_ATC-1  |
| GSE250521 UNI raw DM1/RAI                   |  200 |               -0.0433 |         0.5424 |            -0.0577 |                     4 |               -0.0534 |                   0.0042 | GSM7980873_ATC-2  |
| GSE250521 UNI raw DM1/RAI                   |  200 |                0.0806 |         0.2565 |             0.2827 |                     4 |                0.0267 |                  -0.256  | GSM7980874_ATC-3  |
| GSE250521 UNI raw DM1/RAI                   |  200 |                0.5279 |         0      |             2.2617 |                     5 |                1.9004 |                  -0.3613 | GSM7980875_ATC-4  |
| GSE250521 UNI coord+QC residual DM1/RAI     |  200 |                0.4408 |         0      |             0.313  |                     7 |                0.2053 |                  -0.1078 | GSM7980860_N-1    |
| GSE250521 UNI coord+QC residual DM1/RAI     |  200 |                0.1944 |         0.0058 |             0.0816 |                     5 |                0.0286 |                  -0.0531 | GSM7980861_N-2    |
| GSE250521 UNI coord+QC residual DM1/RAI     |  200 |                0.0063 |         0.9298 |            -0.0163 |                     3 |               -0.0011 |                   0.0152 | GSM7980862_N-3    |
| GSE250521 UNI coord+QC residual DM1/RAI     |  200 |                0.1383 |         0.0508 |             0.1856 |                     6 |                0.0926 |                  -0.093  | GSM7980863_N-4    |
| GSE250521 UNI coord+QC residual DM1/RAI     |  200 |                0.4111 |         0      |             0.3099 |                     5 |                0.0906 |                  -0.2193 | GSM7980864_PTC-1  |
| GSE250521 UNI coord+QC residual DM1/RAI     |  200 |                0.4792 |         0      |             0.3121 |                     8 |                0.2052 |                  -0.1069 | GSM7980865_PTC-2  |
| GSE250521 UNI coord+QC residual DM1/RAI     |  200 |                0.1789 |         0.0113 |             0.0775 |                     4 |                0.0555 |                  -0.022  | GSM7980866_PTC-3  |
| GSE250521 UNI coord+QC residual DM1/RAI     |  200 |                0.1147 |         0.1059 |             0.1173 |                     5 |                0.1023 |                  -0.0151 | GSM7980867_PTC-4  |
| GSE250521 UNI coord+QC residual DM1/RAI     |  200 |                0.1897 |         0.0071 |             0.1636 |                     6 |                0.1325 |                  -0.0311 | GSM7980868_LPTC-1 |
| GSE250521 UNI coord+QC residual DM1/RAI     |  200 |                0.2011 |         0.0043 |             0.1193 |                     5 |                0.018  |                  -0.1013 | GSM7980869_LPTC-2 |
| GSE250521 UNI coord+QC residual DM1/RAI     |  200 |                0.2887 |         0      |             0.108  |                     6 |                0.0549 |                  -0.0531 | GSM7980870_LPTC-3 |
| GSE250521 UNI coord+QC residual DM1/RAI     |  200 |                0.1901 |         0.007  |             0.0997 |                     5 |                0.0273 |                  -0.0724 | GSM7980871_LPTC-4 |
| GSE250521 UNI coord+QC residual DM1/RAI     |  200 |               -0.1007 |         0.1558 |            -0.0413 |                     3 |               -0.0235 |                   0.0178 | GSM7980872_ATC-1  |
| GSE250521 UNI coord+QC residual DM1/RAI     |  200 |                0.0462 |         0.5157 |             0.0094 |                     5 |                0.0097 |                   0.0002 | GSM7980873_ATC-2  |
| GSE250521 UNI coord+QC residual DM1/RAI     |  200 |               -0.0807 |         0.2562 |            -0.0138 |                     5 |               -0.0192 |                  -0.0054 | GSM7980874_ATC-3  |
| GSE250521 UNI coord+QC residual DM1/RAI     |  200 |                0.2492 |         0.0004 |             0.2242 |                     5 |                0.189  |                  -0.0351 | GSM7980875_ATC-4  |
| GSE230424 H&E raw DM1/low-RAI               | 3647 |                0.5599 |         0      |             0.7903 |                     9 |                0.4572 |                  -0.333  | P1                |
| GSE230424 H&E raw DM1/low-RAI               | 3154 |                0.6048 |         0      |             0.752  |                     9 |                0.3873 |                  -0.3647 | P2                |
| GSE230424 H&E raw DM1/low-RAI               | 4062 |                0.7848 |         0      |             0.872  |                     9 |                0.4366 |                  -0.4353 | P3                |
| GSE230424 H&E raw DM1/low-RAI               | 4626 |                0.4413 |         0      |             0.4374 |                     9 |                0.2901 |                  -0.1472 | P4                |
| GSE230424 H&E coord+QC residual DM1/low-RAI | 3647 |                0.1634 |         0      |             0.135  |                     7 |                0.0966 |                  -0.0384 | P1                |
| GSE230424 H&E coord+QC residual DM1/low-RAI | 3154 |                0.3606 |         0      |             0.2261 |                     8 |                0.0803 |                  -0.1458 | P2                |
| GSE230424 H&E coord+QC residual DM1/low-RAI | 4062 |                0.1557 |         0      |             0.1194 |                     5 |                0.0096 |                  -0.1098 | P3                |
| GSE230424 H&E coord+QC residual DM1/low-RAI | 4626 |                0.2202 |         0      |             0.1655 |                     8 |                0.1032 |                  -0.0623 | P4                |

## Decile Means

| evidence                                    |   pred_decile |    n |   mean_obs_group_centered |   median_obs_group_centered |   mean_pred |
|:--------------------------------------------|--------------:|-----:|--------------------------:|----------------------------:|------------:|
| GSE250521 UNI raw DM1/RAI                   |             1 |  320 |                   -0.5663 |                     -0.4627 |     -0.757  |
| GSE250521 UNI raw DM1/RAI                   |             2 |  320 |                   -0.4499 |                     -0.4018 |     -0.5003 |
| GSE250521 UNI raw DM1/RAI                   |             3 |  320 |                   -0.3495 |                     -0.2675 |     -0.3542 |
| GSE250521 UNI raw DM1/RAI                   |             4 |  320 |                   -0.1696 |                     -0.1141 |     -0.2213 |
| GSE250521 UNI raw DM1/RAI                   |             5 |  320 |                   -0.0804 |                     -0.0979 |     -0.1028 |
| GSE250521 UNI raw DM1/RAI                   |             6 |  320 |                   -0.0386 |                     -0.0034 |      0.0082 |
| GSE250521 UNI raw DM1/RAI                   |             7 |  320 |                    0.0037 |                      0.1366 |      0.1264 |
| GSE250521 UNI raw DM1/RAI                   |             8 |  320 |                    0.1928 |                      0.1473 |      0.2777 |
| GSE250521 UNI raw DM1/RAI                   |             9 |  320 |                    0.5448 |                      0.3333 |      0.4664 |
| GSE250521 UNI raw DM1/RAI                   |            10 |  320 |                    0.913  |                      0.7076 |      0.8244 |
| GSE250521 UNI coord+QC residual DM1/RAI     |             1 |  320 |                   -0.0551 |                     -0.0377 |     -0.5862 |
| GSE250521 UNI coord+QC residual DM1/RAI     |             2 |  320 |                   -0.0423 |                     -0.0307 |     -0.3511 |
| GSE250521 UNI coord+QC residual DM1/RAI     |             3 |  320 |                   -0.0345 |                     -0.0223 |     -0.2176 |
| GSE250521 UNI coord+QC residual DM1/RAI     |             4 |  320 |                   -0.0194 |                     -0.0211 |     -0.1225 |
| GSE250521 UNI coord+QC residual DM1/RAI     |             5 |  320 |                   -0.0195 |                     -0.0136 |     -0.0444 |
| GSE250521 UNI coord+QC residual DM1/RAI     |             6 |  320 |                    0.0064 |                     -0.0007 |      0.0365 |
| GSE250521 UNI coord+QC residual DM1/RAI     |             7 |  320 |                    0.0128 |                      0.011  |      0.1318 |
| GSE250521 UNI coord+QC residual DM1/RAI     |             8 |  320 |                    0.0207 |                      0.0115 |      0.2326 |
| GSE250521 UNI coord+QC residual DM1/RAI     |             9 |  320 |                    0.058  |                      0.0319 |      0.3726 |
| GSE250521 UNI coord+QC residual DM1/RAI     |            10 |  320 |                    0.073  |                      0.0474 |      0.6416 |
| GSE230424 H&E raw DM1/low-RAI               |             1 | 1547 |                   -0.3109 |                     -0.3047 |     -0.9317 |
| GSE230424 H&E raw DM1/low-RAI               |             2 | 1549 |                   -0.2692 |                     -0.2753 |     -0.6869 |
| GSE230424 H&E raw DM1/low-RAI               |             3 | 1549 |                   -0.2255 |                     -0.234  |     -0.525  |
| GSE230424 H&E raw DM1/low-RAI               |             4 | 1548 |                   -0.1419 |                     -0.1631 |     -0.3692 |
| GSE230424 H&E raw DM1/low-RAI               |             5 | 1551 |                   -0.0689 |                     -0.072  |     -0.1973 |
| GSE230424 H&E raw DM1/low-RAI               |             6 | 1548 |                    0.0252 |                      0.0092 |     -0.0212 |
| GSE230424 H&E raw DM1/low-RAI               |             7 | 1548 |                    0.1236 |                      0.135  |      0.1648 |
| GSE230424 H&E raw DM1/low-RAI               |             8 | 1549 |                    0.1971 |                      0.2432 |      0.366  |
| GSE230424 H&E raw DM1/low-RAI               |             9 | 1549 |                    0.2818 |                      0.3391 |      0.6037 |
| GSE230424 H&E raw DM1/low-RAI               |            10 | 1551 |                    0.3877 |                      0.4281 |      1.0272 |
| GSE230424 H&E coord+QC residual DM1/low-RAI |             1 | 1547 |                   -0.0861 |                     -0.0852 |     -0.5473 |
| GSE230424 H&E coord+QC residual DM1/low-RAI |             2 | 1549 |                   -0.0572 |                     -0.0584 |     -0.332  |
| GSE230424 H&E coord+QC residual DM1/low-RAI |             3 | 1549 |                   -0.0584 |                     -0.057  |     -0.2193 |
| GSE230424 H&E coord+QC residual DM1/low-RAI |             4 | 1548 |                   -0.0261 |                     -0.03   |     -0.1218 |
| GSE230424 H&E coord+QC residual DM1/low-RAI |             5 | 1551 |                   -0.0054 |                     -0.013  |     -0.032  |
| GSE230424 H&E coord+QC residual DM1/low-RAI |             6 | 1548 |                    0.0281 |                      0.0144 |      0.0553 |
| GSE230424 H&E coord+QC residual DM1/low-RAI |             7 | 1548 |                    0.0341 |                      0.0276 |      0.138  |
| GSE230424 H&E coord+QC residual DM1/low-RAI |             8 | 1549 |                    0.0395 |                      0.0316 |      0.2247 |
| GSE230424 H&E coord+QC residual DM1/low-RAI |             9 | 1549 |                    0.0591 |                      0.0427 |      0.3321 |
| GSE230424 H&E coord+QC residual DM1/low-RAI |            10 | 1551 |                    0.0724 |                      0.051  |      0.5525 |
