# GSE163203 HT-PTC single-cell 외부검증

## 결론

GSE163203 원자료를 내려받아 기술분할 샘플을 생물학적 샘플 단위로 병합했다. 분석 단위는 종양 `PTCwithHT` 3명, `PTCwithoutHT` 5명, 그리고 HT 인접조직 2명이다. 총 110,000 cells를 pseudobulk와 cell-level module detection으로 다시 점수화했다.

이 결과는 Paper 2의 HLA 주장을 **HLA allele/genotype이 아니라 HLA/AP expression module**로 유지하면서 강화한다. 특히 HT 동반 PTC에서 HLA-II/AP와 B/TLS 축이 같은 방향으로 움직이는지 독립 single-cell 자료에서 검정하는 역할이다.

## 핵심 효과

| module | Cohen's d | delta score | exact p | FDR | AUC |
|---|---:|---:|---:|---:|---:|
| Myeloid_DC | 2.84 | 1.28 | 0.035 | 0.070 | 1.00 |
| AP_TLS_composite | 2.78 | 1.29 | 0.035 | 0.070 | 1.00 |
| B_TLS | 2.59 | 1.27 | 0.035 | 0.070 | 1.00 |
| HLA_II_AP | 2.57 | 1.31 | 0.053 | 0.079 | 0.93 |

## 논문 반영 포인트

1. Paper 2 Result에 `independent single-cell validation in GSE163203` 소절을 추가한다.
2. GSE286332 bulk/HT-overlap, TCGA-THCA network, scRNA cell-type concordance, spatial TLS validation 뒤에 GSE163203을 붙이면 “bulk -> TCGA -> single-cell -> spatial” 검증 ladder가 완성된다.
3. 리뷰어 방어문은 명확해야 한다: 이 분석은 HLA 유전자 발현과 antigen-presentation/TLS 생태계 검증이며, 암 코호트에서 HLA allele association을 주장하지 않는다.

## 산출물

- Tables: `project/results/hla_two_paper_synthesis_2026_05_09/gse163203_htptc_scrna_validation/tables/`
- Figures: `project/results/hla_two_paper_synthesis_2026_05_09/gse163203_htptc_scrna_validation/figures/`
- Source data: `project/data/external/GSE163203/`
