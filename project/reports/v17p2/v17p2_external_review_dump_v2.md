# v17p2 External Review Dump V2

## 1. Layer 1
- mean off-diagonal ARI: 0.662
- bootstrap persistence median: 0.903
- permutation p-value: 0.0033
- 50% subsample ARI: 0.757

## 2. Layer 2
- DM1/DM2 size: 109 / 69
- age Welch p: 6.19e-07
- stage high Fisher p: 0.092
- OS events: 8
- Cox table:
| covariate   |      coef |   exp(coef) |   se(coef) |   coef lower 95% |   coef upper 95% |   exp(coef) lower 95% |   exp(coef) upper 95% |   cmp to |         z |           p |   -log2(p) |
|:------------|----------:|------------:|-----------:|-----------------:|-----------------:|----------------------:|----------------------:|---------:|----------:|------------:|-----------:|
| cluster_bin | -0.297116 |    0.742958 |  0.812884  |       -1.89034   |         1.29611  |              0.151021 |               3.65504 |        0 | -0.365508 | 0.714732    |   0.484526 |
| age         |  0.161835 |    1.17567  |  0.0484088 |        0.0669556 |         0.256715 |              1.06925  |               1.29268 |        0 |  3.34309  | 0.000828502 |  10.2372   |
| sex_male    |  0.125102 |    1.13326  |  0.827432  |       -1.49664   |         1.74684  |              0.223882 |               5.73645 |        0 |  0.151193 | 0.879823    |   0.184715 |
| stage_high  |  1.20589  |    3.33974  |  1.16641   |       -1.08023   |         3.49202  |              0.339518 |              32.8522  |        0 |  1.03385  | 0.301206    |   1.73118  |

## 3. Layer 3
- hallmark/pathway proxy FDR<0.05 count: 10
- MAPK activity DM2-DM1: -1.080
- top pathways:
| Term                               |       NES |   FDR q-val | gene_set_source   |
|:-----------------------------------|----------:|------------:|:------------------|
| Hallmark_Apoptosis                 | -1.39145  | 6.08834e-25 | proxy             |
| Hallmark_Inflammatory_Response     | -1.24453  | 7.27583e-22 | proxy             |
| Hallmark_EMT                       | -1.15523  | 1.15713e-17 | proxy             |
| Hallmark_TGF_beta                  | -1.20013  | 3.99996e-17 | proxy             |
| Hallmark_Interferon_Gamma          | -1.09781  | 1.69677e-16 | proxy             |
| Hallmark_Oxidative_Phosphorylation |  0.849306 | 1.37826e-07 | proxy             |
| Hallmark_E2F_Targets               | -0.610389 | 0.000101877 | proxy             |
| Hallmark_G2M_Checkpoint            | -0.417601 | 0.00543268  | proxy             |

## 4. Layer 4
- bulk external tested cohorts: 2
- robust cohorts (DIA-AUC>0.85): 2
- scRNA cells scored: 66015
- best transfer:
| cohort   | classifier   |      AUC |   DIA_AUC |   n |
|:---------|:-------------|---------:|----------:|----:|
| GSE76039 | LogReg       | 0.990196 |  0.990196 |  37 |
| GSE27155 | LogReg       | 0.96875  |  0.96875  |  54 |
| GSE76039 | RandomForest | 0.901961 |  0.901961 |  37 |
| GSE27155 | RandomForest | 0.876736 |  0.876736 |  54 |

## 5. Layer 5
- ComBat threshold (identifiability < 0.5): 1.0
- pancancer tested: 3
- max pan-cancer DIA-AUC: 0.9911807292352245

## 6. Self-assessment
| Item | Score | Note |
|---|---:|---|
| Statistical robustness | 4/5 | persistence strong, ARI moderate |
| Clinical stratification | 2/5 | age strong, survival weak |
| Biological interpretability | 4/5 | pathway proxy clear |
| External validation | 3/5 | 2 bulk cohorts + scRNA |
| Method audit | 3/5 | THCA audit strong, pan-cancer partial |
| Wet validation | 0/5 | 없음 |

## 7. Venue
- npj Precision Oncology: 55%
- Genome Medicine: 25%
- Bioinformatics: 75%
- Nature Communications: 8%

## 8. 다음 결정 3개
1. GSE213647 gene harmonization을 더 깊게 밀어 bulk robust cohort를 3개 이상으로 회복할지.
2. TERT/raw fusion recovery를 별도 1주 sprint로 분리할지.
3. 지금 바로 manuscript draft로 들어갈지, 아니면 OS/clinical metadata를 더 보강할지.
