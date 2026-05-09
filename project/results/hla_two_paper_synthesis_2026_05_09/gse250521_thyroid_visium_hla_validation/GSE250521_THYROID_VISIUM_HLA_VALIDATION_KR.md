# GSE250521 thyroid Visium HLA/AP spatial 외부검증

## 결론

GSE250521 Visium raw matrix를 직접 다운로드해 16개 slide를 재분석했다. QC 후 구성은 N slides=4, spots=14,740, PTC slides=4, spots=16,223, LPTC slides=4, spots=13,634, ATC slides=4, spots=11,276 이다.

이 데이터는 HT-specific이 아니다. Paper 2에서는 **thyroid cancer progression/spatial ecology generalization**으로만 사용하고, HT-overlap의 핵심 증거는 GSE138198/GSE163203에 둔다. HLA allele/genotype/risk claim은 금지한다.

## Stage trend across N/PTC/LPTC/ATC

| module | Spearman rho | p | FDR | linear slope |
|---|---:|---:|---:|---:|
| CD74_MIF_axis | 0.76 | 0.000571 | 0.00274 | 0.42 |
| HLA_I | 0.75 | 0.000782 | 0.00274 | 0.52 |
| Myeloid_DC | 0.70 | 0.00237 | 0.00552 | 0.47 |
| T_IFNG | 0.62 | 0.0107 | 0.0179 | 0.52 |
| HLA_II_AP | 0.61 | 0.0128 | 0.0179 | 0.39 |
| AP_TLS_composite | 0.52 | 0.0383 | 0.0447 | 0.22 |
| B_TLS | 0.06 | 0.823 | 0.823 | 0.05 |

## Stage-specific contrasts vs normal thyroid

| contrast | module | Cohen's d | delta score | exact p | AUC |
|---|---|---:|---:|---:|---:|
| ATC_vs_N | Myeloid_DC | 3.07 | 1.56 | 0.042 | 1.00 |
| ATC_vs_N | HLA_I | 2.94 | 1.77 | 0.042 | 1.00 |
| ATC_vs_N | HLA_II_AP | 2.24 | 1.22 | 0.070 | 0.94 |
| ATC_vs_N | T_IFNG | 2.01 | 1.74 | 0.042 | 1.00 |
| ATC_vs_N | AP_TLS_composite | 1.81 | 0.79 | 0.099 | 0.88 |
| ATC_vs_N | B_TLS | 0.98 | 0.37 | 0.296 | 0.69 |
| LPTC_vs_N | HLA_II_AP | 2.99 | 1.28 | 0.042 | 1.00 |
| LPTC_vs_N | Myeloid_DC | 2.90 | 1.27 | 0.042 | 1.00 |
| LPTC_vs_N | AP_TLS_composite | 2.43 | 0.72 | 0.042 | 1.00 |
| LPTC_vs_N | HLA_I | 1.82 | 0.65 | 0.070 | 0.94 |
| LPTC_vs_N | T_IFNG | 0.61 | 0.18 | 0.437 | 0.62 |
| LPTC_vs_N | B_TLS | 0.38 | 0.16 | 0.803 | 0.44 |
| PTC_vs_N | Myeloid_DC | 1.94 | 1.23 | 0.042 | 1.00 |
| PTC_vs_N | HLA_II_AP | 1.61 | 0.98 | 0.099 | 0.88 |
| PTC_vs_N | HLA_I | 1.39 | 0.72 | 0.127 | 0.81 |
| PTC_vs_N | AP_TLS_composite | 1.29 | 0.89 | 0.099 | 0.88 |
| PTC_vs_N | B_TLS | 1.00 | 0.81 | 0.127 | 0.81 |
| PTC_vs_N | T_IFNG | 0.42 | 0.20 | 0.662 | 0.56 |

## Spatial burden above N p90

| stage | module | fraction above N p90 | mean raw spot score |
|---|---|---:|---:|
| N | AP_TLS_composite | 0.101 | 0.24 |
| N | B_TLS | 0.103 | 0.18 |
| N | HLA_I | 0.102 | 0.81 |
| N | HLA_II_AP | 0.099 | 0.29 |
| N | T_IFNG | 0.101 | 0.03 |
| PTC | AP_TLS_composite | 0.472 | 0.46 |
| PTC | B_TLS | 0.324 | 0.30 |
| PTC | HLA_I | 0.450 | 1.04 |
| PTC | HLA_II_AP | 0.546 | 0.61 |
| PTC | T_IFNG | 0.193 | 0.06 |
| LPTC | AP_TLS_composite | 0.533 | 0.44 |
| LPTC | B_TLS | 0.156 | 0.17 |
| LPTC | HLA_I | 0.472 | 1.04 |
| LPTC | HLA_II_AP | 0.670 | 0.70 |
| LPTC | T_IFNG | 0.192 | 0.06 |
| ATC | AP_TLS_composite | 0.604 | 0.46 |
| ATC | B_TLS | 0.120 | 0.19 |
| ATC | HLA_I | 0.879 | 1.46 |
| ATC | HLA_II_AP | 0.738 | 0.73 |
| ATC | T_IFNG | 0.668 | 0.31 |

## 논문 반영 포인트

1. Positive 방향이면 Paper 2 supplement external validation panel에 `GSE250521 Visium progression spatial generalization`으로 넣는다.
2. HT-specific 결론에는 직접 쓰지 않는다.
3. `HLA allele`, `genotype`, `risk allele` 문구를 쓰지 않는다.

## 산출물

- Tables: `project/results/hla_two_paper_synthesis_2026_05_09/gse250521_thyroid_visium_hla_validation/tables/`
- Figures: `project/results/hla_two_paper_synthesis_2026_05_09/gse250521_thyroid_visium_hla_validation/figures/`
- Source data: `project/data/external/GSE250521/GSE250521_RAW.tar`
