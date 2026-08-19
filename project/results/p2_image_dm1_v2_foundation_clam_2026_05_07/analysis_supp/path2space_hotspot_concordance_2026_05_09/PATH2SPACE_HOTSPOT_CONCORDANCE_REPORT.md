# Path2Space-inspired hotspot concordance

## Verdict

- GSE250521 top-decile hotspot precision: **0.350**; lift **3.50x**; Fisher OR **6.92**; permutation p = **0.0010**.
- GSE230424 top-decile hotspot precision: **0.333**; lift **3.32x**; Fisher OR **6.22**; permutation p = **0.0010**.
- GSE230424 coord+QC residual top-decile precision: **0.163**; lift **1.63x**; Fisher OR **1.90**; permutation p = **0.0010**.

## Interpretation

The hotspot test asks whether predicted high-DM1/RAI regions recover observed high-DM1/RAI regions within each slide/sample, not merely whether continuous scores are correlated. Top-decile enrichment supports spatial localization. The residual-target version is expectedly weaker but is the strictest localization test after coordinate+QC removal.

## Summary

| evidence                                    |   quantile |   top_percent |     n |   tp |   fp |   fn |    tn |   obs_hotspot_rate |   pred_hotspot_rate |   precision |   recall |   precision_lift |   jaccard |   expected_tp_independent |   tp_over_expected |   fisher_odds_ratio |   fisher_p_greater |   permutation_p_greater_equal_tp |
|:--------------------------------------------|-----------:|--------------:|------:|-----:|-----:|-----:|------:|-------------------:|--------------------:|------------:|---------:|-----------------:|----------:|--------------------------:|-------------------:|--------------------:|-------------------:|---------------------------------:|
| GSE250521 UNI DM1/RAI smoothed              |        0.8 |            20 |  3200 |  282 |  358 |  358 |  2202 |             0.2    |              0.2    |      0.4406 |   0.4406 |           2.2031 |    0.2826 |                    128    |             2.2031 |              4.8451 |                  0 |                            0.001 |
| GSE250521 UNI DM1/RAI smoothed              |        0.9 |            10 |  3200 |  112 |  208 |  208 |  2672 |             0.1    |              0.1    |      0.35   |   0.35   |           3.5    |    0.2121 |                     32    |             3.5    |              6.9172 |                  0 |                            0.001 |
| GSE230424 H&E DM1/low-RAI smoothed          |        0.8 |            20 | 15489 | 1550 | 1550 | 1550 | 10839 |             0.2001 |              0.2001 |      0.5    |   0.5    |           2.4982 |    0.3333 |                    620.44 |             2.4982 |              6.9929 |                  0 |                            0.001 |
| GSE230424 H&E DM1/low-RAI smoothed          |        0.9 |            10 | 15489 |  516 | 1035 | 1035 | 12903 |             0.1001 |              0.1001 |      0.3327 |   0.3327 |           3.3224 |    0.1995 |                    155.31 |             3.3224 |              6.2153 |                  0 |                            0.001 |
| GSE230424 H&E DM1/low-RAI coord+QC residual |        0.8 |            20 | 15489 |  882 | 2218 | 2218 | 10171 |             0.2001 |              0.2001 |      0.2845 |   0.2845 |           1.4216 |    0.1659 |                    620.44 |             1.4216 |              1.8235 |                  0 |                            0.001 |
| GSE230424 H&E DM1/low-RAI coord+QC residual |        0.9 |            10 | 15489 |  253 | 1298 | 1298 | 12640 |             0.1001 |              0.1001 |      0.1631 |   0.1631 |           1.629  |    0.0888 |                    155.31 |             1.629  |              1.8981 |                  0 |                            0.001 |

## Per-Slide/Sample Detail

| evidence                                    |   quantile | group             |    n |   tp |   fp |   fn |   tn |   precision |   recall |   jaccard |   fisher_odds_ratio |   fisher_p_greater |
|:--------------------------------------------|-----------:|:------------------|-----:|-----:|-----:|-----:|-----:|------------:|---------:|----------:|--------------------:|-------------------:|
| GSE250521 UNI DM1/RAI smoothed              |        0.8 | GSM7980860_N-1    |  200 |   21 |   19 |   19 |  141 |      0.525  |   0.525  |    0.3559 |              8.2022 |             0      |
| GSE250521 UNI DM1/RAI smoothed              |        0.8 | GSM7980861_N-2    |  200 |   22 |   18 |   18 |  142 |      0.55   |   0.55   |    0.3793 |              9.642  |             0      |
| GSE250521 UNI DM1/RAI smoothed              |        0.8 | GSM7980862_N-3    |  200 |    9 |   31 |   31 |  129 |      0.225  |   0.225  |    0.1268 |              1.2081 |             0.4028 |
| GSE250521 UNI DM1/RAI smoothed              |        0.8 | GSM7980863_N-4    |  200 |   20 |   20 |   20 |  140 |      0.5    |   0.5    |    0.3333 |              7      |             0      |
| GSE250521 UNI DM1/RAI smoothed              |        0.8 | GSM7980864_PTC-1  |  200 |   11 |   29 |   29 |  131 |      0.275  |   0.275  |    0.1594 |              1.7134 |             0.1354 |
| GSE250521 UNI DM1/RAI smoothed              |        0.8 | GSM7980865_PTC-2  |  200 |   28 |   12 |   12 |  148 |      0.7    |   0.7    |    0.5385 |             28.7778 |             0      |
| GSE250521 UNI DM1/RAI smoothed              |        0.8 | GSM7980866_PTC-3  |  200 |   20 |   20 |   20 |  140 |      0.5    |   0.5    |    0.3333 |              7      |             0      |
| GSE250521 UNI DM1/RAI smoothed              |        0.8 | GSM7980867_PTC-4  |  200 |   19 |   21 |   21 |  139 |      0.475  |   0.475  |    0.3115 |              5.9887 |             0      |
| GSE250521 UNI DM1/RAI smoothed              |        0.8 | GSM7980868_LPTC-1 |  200 |   26 |   14 |   14 |  146 |      0.65   |   0.65   |    0.4815 |             19.3673 |             0      |
| GSE250521 UNI DM1/RAI smoothed              |        0.8 | GSM7980869_LPTC-2 |  200 |   11 |   29 |   29 |  131 |      0.275  |   0.275  |    0.1594 |              1.7134 |             0.1354 |
| GSE250521 UNI DM1/RAI smoothed              |        0.8 | GSM7980870_LPTC-3 |  200 |   16 |   24 |   24 |  136 |      0.4    |   0.4    |    0.25   |              3.7778 |             0.0008 |
| GSE250521 UNI DM1/RAI smoothed              |        0.8 | GSM7980871_LPTC-4 |  200 |   24 |   16 |   16 |  144 |      0.6    |   0.6    |    0.4286 |             13.5    |             0      |
| GSE250521 UNI DM1/RAI smoothed              |        0.8 | GSM7980872_ATC-1  |  200 |   10 |   30 |   30 |  130 |      0.25   |   0.25   |    0.1429 |              1.4444 |             0.2491 |
| GSE250521 UNI DM1/RAI smoothed              |        0.8 | GSM7980873_ATC-2  |  200 |    9 |   31 |   31 |  129 |      0.225  |   0.225  |    0.1268 |              1.2081 |             0.4028 |
| GSE250521 UNI DM1/RAI smoothed              |        0.8 | GSM7980874_ATC-3  |  200 |    8 |   32 |   32 |  128 |      0.2    |   0.2    |    0.1111 |              1      |             0.5771 |
| GSE250521 UNI DM1/RAI smoothed              |        0.8 | GSM7980875_ATC-4  |  200 |   28 |   12 |   12 |  148 |      0.7    |   0.7    |    0.5385 |             28.7778 |             0      |
| GSE250521 UNI DM1/RAI smoothed              |        0.9 | GSM7980860_N-1    |  200 |    9 |   11 |   11 |  169 |      0.45   |   0.45   |    0.2903 |             12.5702 |             0      |
| GSE250521 UNI DM1/RAI smoothed              |        0.9 | GSM7980861_N-2    |  200 |    8 |   12 |   12 |  168 |      0.4    |   0.4    |    0.25   |              9.3333 |             0.0001 |
| GSE250521 UNI DM1/RAI smoothed              |        0.9 | GSM7980862_N-3    |  200 |    4 |   16 |   16 |  164 |      0.2    |   0.2    |    0.1111 |              2.5625 |             0.1222 |
| GSE250521 UNI DM1/RAI smoothed              |        0.9 | GSM7980863_N-4    |  200 |    7 |   13 |   13 |  167 |      0.35   |   0.35   |    0.2121 |              6.9172 |             0.0012 |
| GSE250521 UNI DM1/RAI smoothed              |        0.9 | GSM7980864_PTC-1  |  200 |    3 |   17 |   17 |  163 |      0.15   |   0.15   |    0.0811 |              1.692  |             0.3213 |
| GSE250521 UNI DM1/RAI smoothed              |        0.9 | GSM7980865_PTC-2  |  200 |   11 |    9 |    9 |  171 |      0.55   |   0.55   |    0.3793 |             23.2222 |             0      |
| GSE250521 UNI DM1/RAI smoothed              |        0.9 | GSM7980866_PTC-3  |  200 |   12 |    8 |    8 |  172 |      0.6    |   0.6    |    0.4286 |             32.25   |             0      |
| GSE250521 UNI DM1/RAI smoothed              |        0.9 | GSM7980867_PTC-4  |  200 |   10 |   10 |   10 |  170 |      0.5    |   0.5    |    0.3333 |             17      |             0      |
| GSE250521 UNI DM1/RAI smoothed              |        0.9 | GSM7980868_LPTC-1 |  200 |   12 |    8 |    8 |  172 |      0.6    |   0.6    |    0.4286 |             32.25   |             0      |
| GSE250521 UNI DM1/RAI smoothed              |        0.9 | GSM7980869_LPTC-2 |  200 |    0 |   20 |   20 |  160 |      0      |   0      |    0      |              0      |             1      |
| GSE250521 UNI DM1/RAI smoothed              |        0.9 | GSM7980870_LPTC-3 |  200 |    2 |   18 |   18 |  162 |      0.1    |   0.1    |    0.0526 |              1      |             0.6218 |
| GSE250521 UNI DM1/RAI smoothed              |        0.9 | GSM7980871_LPTC-4 |  200 |    7 |   13 |   13 |  167 |      0.35   |   0.35   |    0.2121 |              6.9172 |             0.0012 |
| GSE250521 UNI DM1/RAI smoothed              |        0.9 | GSM7980872_ATC-1  |  200 |    6 |   14 |   14 |  166 |      0.3    |   0.3    |    0.1765 |              5.0816 |             0.0073 |
| GSE250521 UNI DM1/RAI smoothed              |        0.9 | GSM7980873_ATC-2  |  200 |    1 |   19 |   19 |  161 |      0.05   |   0.05   |    0.0256 |              0.446  |             0.8915 |
| GSE250521 UNI DM1/RAI smoothed              |        0.9 | GSM7980874_ATC-3  |  200 |    3 |   17 |   17 |  163 |      0.15   |   0.15   |    0.0811 |              1.692  |             0.3213 |
| GSE250521 UNI DM1/RAI smoothed              |        0.9 | GSM7980875_ATC-4  |  200 |   17 |    3 |    3 |  177 |      0.85   |   0.85   |    0.7391 |            334.333  |             0      |
| GSE230424 H&E DM1/low-RAI smoothed          |        0.8 | P1                | 3647 |  365 |  365 |  365 | 2552 |      0.5    |   0.5    |    0.3333 |              6.9918 |             0      |
| GSE230424 H&E DM1/low-RAI smoothed          |        0.8 | P2                | 3154 |  245 |  386 |  386 | 2137 |      0.3883 |   0.3883 |    0.2409 |              3.514  |             0      |
| GSE230424 H&E DM1/low-RAI smoothed          |        0.8 | P3                | 4062 |  471 |  342 |  342 | 2907 |      0.5793 |   0.5793 |    0.4078 |             11.7061 |             0      |
| GSE230424 H&E DM1/low-RAI smoothed          |        0.8 | P4                | 4626 |  469 |  457 |  457 | 3243 |      0.5065 |   0.5065 |    0.3391 |              7.2826 |             0      |
| GSE230424 H&E DM1/low-RAI smoothed          |        0.9 | P1                | 3647 |  108 |  257 |  257 | 3025 |      0.2959 |   0.2959 |    0.1736 |              4.9463 |             0      |
| GSE230424 H&E DM1/low-RAI smoothed          |        0.9 | P2                | 3154 |   99 |  217 |  217 | 2621 |      0.3133 |   0.3133 |    0.1857 |              5.5104 |             0      |
| GSE230424 H&E DM1/low-RAI smoothed          |        0.9 | P3                | 4062 |  130 |  277 |  277 | 3378 |      0.3194 |   0.3194 |    0.1901 |              5.7233 |             0      |
| GSE230424 H&E DM1/low-RAI smoothed          |        0.9 | P4                | 4626 |  179 |  284 |  284 | 3879 |      0.3866 |   0.3866 |    0.2396 |              8.6087 |             0      |
| GSE230424 H&E DM1/low-RAI coord+QC residual |        0.8 | P1                | 3647 |  213 |  517 |  517 | 2400 |      0.2918 |   0.2918 |    0.1708 |              1.9125 |             0      |
| GSE230424 H&E DM1/low-RAI coord+QC residual |        0.8 | P2                | 3154 |  193 |  438 |  438 | 2085 |      0.3059 |   0.3059 |    0.1805 |              2.0976 |             0      |
| GSE230424 H&E DM1/low-RAI coord+QC residual |        0.8 | P3                | 4062 |  131 |  682 |  682 | 2567 |      0.1611 |   0.1611 |    0.0876 |              0.723  |             0.9994 |
| GSE230424 H&E DM1/low-RAI coord+QC residual |        0.8 | P4                | 4626 |  345 |  581 |  581 | 3119 |      0.3726 |   0.3726 |    0.2289 |              3.1877 |             0      |
| GSE230424 H&E DM1/low-RAI coord+QC residual |        0.9 | P1                | 3647 |   73 |  292 |  292 | 2990 |      0.2    |   0.2    |    0.1111 |              2.5599 |             0      |
| GSE230424 H&E DM1/low-RAI coord+QC residual |        0.9 | P2                | 3154 |   36 |  280 |  280 | 2558 |      0.1139 |   0.1139 |    0.0604 |              1.1746 |             0.2214 |
| GSE230424 H&E DM1/low-RAI coord+QC residual |        0.9 | P3                | 4062 |   18 |  389 |  389 | 3266 |      0.0442 |   0.0442 |    0.0226 |              0.3885 |             1      |
| GSE230424 H&E DM1/low-RAI coord+QC residual |        0.9 | P4                | 4626 |  126 |  337 |  337 | 3826 |      0.2721 |   0.2721 |    0.1575 |              4.2448 |             0      |
