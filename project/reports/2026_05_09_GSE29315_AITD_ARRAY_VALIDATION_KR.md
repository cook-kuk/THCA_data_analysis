# GSE29315 AITD array 외부검증

## 결론

GSE29315는 Affymetrix U95 array 기반 thyroid neoplasia/thyroiditis cohort이며, Hashimoto n=6, hyperplasia n=8가 포함된다. GPL8300에서 target HLA/AP/TLS gene 63개가 매핑되었다.

이 결과는 Paper 4의 allele genetics를 직접 반복검증하는 자료가 아니라, HT 조직에서 HLA/AP/TLS expression axis가 올라가는지 확인하는 보조 기전검증이다.

## Primary contrast: Hashimoto vs Hyperplasia

| module | Cohen's d | delta score | exact p | AUC | all-test FDR |
|---|---:|---:|---:|---:|---:|
| HLA_II_AP | 5.45 | 2.29 | 0.0007 | 1.00 | 0.0023 |
| HLA_I | 4.40 | 1.97 | 0.0007 | 1.00 | 0.0023 |
| AP_TLS_composite | 3.52 | 2.37 | 0.0007 | 1.00 | 0.0023 |
| T_IFNG | 3.11 | 2.57 | 0.0010 | 0.98 | 0.0023 |
| Myeloid_DC | 2.98 | 1.76 | 0.0013 | 0.96 | 0.0023 |
| CD74_MIF_axis | 2.72 | 0.99 | 0.0013 | 0.96 | 0.0023 |
| B_TLS | 2.22 | 2.44 | 0.0010 | 0.98 | 0.0023 |

## 논문 반영 포인트

1. GSE248205 spatial 결과의 독립 array corroboration으로 supplement에 배치한다.
2. 오래된 platform이라 gene coverage가 제한적이므로, main claim은 GSE248205와 genetics 자료에 둔다.
3. genotype/allele replication이라고 쓰지 않는다.

## 산출물

- Tables: `project/results/hla_two_paper_synthesis_2026_05_09/gse29315_aitd_array_validation/tables/`
- Figures: `project/results/hla_two_paper_synthesis_2026_05_09/gse29315_aitd_array_validation/figures/`
- Source data: `project/data/external/GSE29315/`
