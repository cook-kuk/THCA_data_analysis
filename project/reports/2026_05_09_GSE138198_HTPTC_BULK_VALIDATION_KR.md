# GSE138198 HT/PTC bulk 외부검증

## 결론

GSE138198은 HT, PTC with HT background, PTC without HT background, mPTC, normal thyroid를 모두 포함하는 bulk expression 자료다. 이번 재분석에서는 TN n=3, HT n=13, PTCwithoutHT n=6, PTCwithHT n=8, mPTC n=6 구조로, 154개 probe가 86개 target gene에 매핑되었다.

Paper 2에는 이 결과를 **GSE286332와 GSE163203 사이를 잇는 bulk validation**으로 추가할 수 있다. 핵심 주장은 암 코호트 HLA allele association이 아니라, `PTCwithHT`에서 HLA-II/AP 및 TLS/B-cell expression module이 `PTCwithoutHT`보다 높은지 검증하는 것이다.

## Primary contrast: PTCwithHT vs PTCwithoutHT

| module | Cohen's d | delta score | exact p | all-test FDR | AUC |
|---|---:|---:|---:|---:|---:|
| T_IFNG | 1.23 | 0.79 | 0.0373 | 0.0895 | 0.90 |
| HLA_I | 1.19 | 0.86 | 0.0509 | 0.0965 | 0.81 |
| HLA_II_AP | 1.05 | 0.61 | 0.0819 | 0.1404 | 0.77 |
| AP_TLS_composite | 1.00 | 0.56 | 0.0902 | 0.1473 | 0.75 |
| B_TLS | 0.82 | 0.51 | 0.1548 | 0.2064 | 0.71 |
| Myeloid_DC | 0.59 | 0.35 | 0.3006 | 0.3732 | 0.73 |

## Disease-context contrasts

| contrast | module | Cohen's d | exact p | AUC |
|---|---|---:|---:|---:|
| HT_vs_TN | AP_TLS_composite | 2.90 | 0.0036 | 1.00 |
| HT_vs_TN | B_TLS | 2.26 | 0.0053 | 0.97 |
| HT_vs_TN | HLA_II_AP | 3.36 | 0.0036 | 1.00 |
| PTCwithHT_vs_TN | AP_TLS_composite | 2.87 | 0.0120 | 1.00 |
| PTCwithHT_vs_TN | B_TLS | 1.17 | 0.1084 | 0.88 |
| PTCwithHT_vs_TN | HLA_II_AP | 4.38 | 0.0120 | 1.00 |
| PTCwithoutHT_vs_TN | AP_TLS_composite | 1.85 | 0.0353 | 1.00 |
| PTCwithoutHT_vs_TN | B_TLS | 0.29 | 0.8000 | 0.56 |
| PTCwithoutHT_vs_TN | HLA_II_AP | 3.20 | 0.0235 | 1.00 |

## 논문 반영 포인트

1. Paper 2 Result에 `independent HT/PTC bulk validation in GSE138198`를 추가한다.
2. Figure ladder는 `GSE286332 bulk -> GSE138198 HT/PTC bulk -> GSE163203 scRNA -> TCGA-THCA network -> spatial TLS` 순서가 가장 강하다.
3. 문구는 반드시 `HLA-II antigen-presentation expression module`로 제한한다. 이 자료는 germline HLA allele을 직접 검정하지 않는다.

## 산출물

- Tables: `project/results/hla_two_paper_synthesis_2026_05_09/gse138198_htptc_bulk_validation/tables/`
- Figures: `project/results/hla_two_paper_synthesis_2026_05_09/gse138198_htptc_bulk_validation/figures/`
- Source data: `project/data/external/GSE138198/`
