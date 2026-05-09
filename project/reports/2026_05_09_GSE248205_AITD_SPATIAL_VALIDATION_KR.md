# GSE248205 AITD spatial HLA/AP 검증

## 결론

GSE248205는 2024 Nature Communications AITD spatial transcriptomics 자료로, 이번 재분석에서는 Control n=2, HT n=3, GD n=3, 총 16,985 in-tissue spots를 사용했다.

Paper 4에는 이 데이터를 **HLA allele genetics의 조직 기전 보강**으로 붙일 수 있다. 즉, 유전좌위 주장은 Paper 4의 HLA/GWAS 자료가 담당하고, GSE248205는 HT/GD 조직에서 HLA-II antigen-presentation, CD74/MIF, TLS/B-cell 축이 공간적으로 확장되는지를 보여준다.

Paper 2에는 암 allele 주장이 아니라 `HT background antigen-presentation ecosystem`을 설명하는 외부 tissue context로만 제한해서 넣어야 한다.

## Sample-level contrasts

| contrast | module | Cohen's d | delta score | exact p | AUC |
|---|---|---:|---:|---:|---:|
| AITD_vs_Control | AP_TLS_composite | 2.26 | 1.38 | 0.069 | 1.00 |
| AITD_vs_Control | B_TLS | 2.57 | 1.61 | 0.069 | 1.00 |
| AITD_vs_Control | CD74_MIF_axis | 3.62 | 1.83 | 0.069 | 1.00 |
| AITD_vs_Control | HLA_II_AP | 1.77 | 1.15 | 0.138 | 0.83 |
| GD_vs_Control | AP_TLS_composite | 1.89 | 0.88 | 0.273 | 1.00 |
| GD_vs_Control | B_TLS | 3.38 | 1.05 | 0.182 | 1.00 |
| GD_vs_Control | CD74_MIF_axis | 2.68 | 1.53 | 0.273 | 1.00 |
| GD_vs_Control | HLA_II_AP | 1.06 | 0.71 | 0.455 | 0.67 |
| HT_vs_Control | AP_TLS_composite | 12.24 | 1.88 | 0.182 | 1.00 |
| HT_vs_Control | B_TLS | 9.19 | 2.17 | 0.182 | 1.00 |
| HT_vs_Control | CD74_MIF_axis | 8.95 | 2.13 | 0.182 | 1.00 |
| HT_vs_Control | HLA_II_AP | 21.54 | 1.59 | 0.182 | 1.00 |

## Spatial burden

`fraction_above_control_p90`는 control spot 분포의 90th percentile보다 높은 spot 비율이다.

| group | module | fraction above control p90 | mean spot score |
|---|---|---:|---:|
| Control | AP_TLS_composite | 0.104 | 1.66 |
| GD | AP_TLS_composite | 0.629 | 3.75 |
| HT | AP_TLS_composite | 0.993 | 5.76 |
| Control | B_TLS | 0.108 | 1.15 |
| GD | B_TLS | 0.627 | 3.99 |
| HT | B_TLS | 0.997 | 6.62 |
| Control | CD74_MIF_axis | 0.093 | 1.45 |
| GD | CD74_MIF_axis | 0.635 | 3.46 |
| HT | CD74_MIF_axis | 0.946 | 4.22 |
| Control | HLA_II_AP | 0.096 | 2.16 |
| GD | HLA_II_AP | 0.437 | 3.50 |
| HT | HLA_II_AP | 0.964 | 4.91 |

## 논문 반영 포인트

1. Paper 4 Discussion/Mechanism figure에 `spatial AITD tissue validation` 패널을 추가한다.
2. GD와 HT를 분리해 보여주면 class-II/AP 공통축과 disease-specific tissue architecture를 동시에 주장할 수 있다.
3. 이 자료는 genotype 자료가 아니므로, allele association replication이라고 쓰면 안 된다.

## 산출물

- Tables: `project/results/hla_two_paper_synthesis_2026_05_09/gse248205_aitd_spatial_validation/tables/`
- Figures: `project/results/hla_two_paper_synthesis_2026_05_09/gse248205_aitd_spatial_validation/figures/`
- Source data: `project/data/external/GSE248205/`
