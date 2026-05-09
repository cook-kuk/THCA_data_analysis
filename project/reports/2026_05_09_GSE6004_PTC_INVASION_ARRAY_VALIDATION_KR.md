# GSE6004 PTC center/invasion array 외부검증

## 결론

GSE6004는 PTC 조직에서 normal thyroid, tumor center, invasive area를 비교한 GPL570 microarray 자료다. 이번 분석은 GEO series matrix와 GPL570 annotation을 직접 파싱해 HLA/AP/TLS target gene module을 재계산했다. 표본 구성은 Normal n=4, Center n=7, Invasion n=7 이다.

Paper 2에는 이 결과를 positive evidence로 쓰지 않는다. HT-specific cohort가 아니고, invasion-vs-normal 방향도 HLA/AP/TLS 상승이 아니므로 **negative/specificity stress-test**로만 보관한다. 핵심 외부검증은 GSE138198, GSE163203, GSE213647, GSE248205에 둔다.

## Primary contrast: invasion vs normal

| module | Cohen's d | delta score | exact permutation p | FDR | AUC |
|---|---:|---:|---:|---:|---:|
| HLA_I | -0.24 | -0.27 | 0.704 | 0.829 | 0.32 |
| HLA_II_AP | -0.42 | -0.42 | 0.514 | 0.829 | 0.39 |
| T_IFNG | -0.50 | -0.54 | 0.453 | 0.829 | 0.43 |
| AP_TLS_composite | -0.59 | -0.57 | 0.384 | 0.829 | 0.43 |
| Myeloid_DC | -0.63 | -0.60 | 0.356 | 0.829 | 0.43 |
| B_TLS | -0.74 | -0.72 | 0.326 | 0.829 | 0.39 |
| CD74_MIF_axis | -0.75 | -0.40 | 0.272 | 0.829 | 0.32 |

## Paired complete triplets

T2/T7/T8/T18 네 환자는 normal-center-invasion triplet이 모두 있어 paired delta를 보조로 계산했다.

| module | n pairs | mean invasion-normal delta | Wilcoxon p |
|---|---:|---:|---:|
| HLA_I | 4 | 0.05 | 0.875 |
| HLA_II_AP | 4 | -0.09 | 0.875 |
| AP_TLS_composite | 4 | -0.29 | 0.375 |
| CD74_MIF_axis | 4 | -0.31 | 0.625 |
| Myeloid_DC | 4 | -0.42 | 0.625 |
| T_IFNG | 4 | -0.45 | 0.250 |
| B_TLS | 4 | -0.48 | 0.375 |

## 논문 반영 포인트

1. Main text에는 넣지 않는다. Reviewer가 "모든 PTC geography에서 항상 HLA/AP/TLS가 올라가냐"라고 물을 때 보조 stress-test로만 사용한다.
2. Supplement에도 넣는다면 "not universal across tumor geography datasets"라는 한계/특이성 패널로 둔다.
3. `HLA allele`, `genotype`, `risk allele`, `HT-specific` 표현은 쓰지 않는다.

## 산출물

- Mapped target genes: 86
- Tables: `project/results/hla_two_paper_synthesis_2026_05_09/gse6004_ptc_invasion_array_validation/tables/`
- Figures: `project/results/hla_two_paper_synthesis_2026_05_09/gse6004_ptc_invasion_array_validation/figures/`
- Source data: `project/data/external/GSE6004/`
