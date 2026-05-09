# GSE213647 Korean thyroid bulk RNA-seq 외부검증

## 결론

GSE213647은 한국 연구진이 생성한 thyroid bulk RNA-seq 632 samples 자료다. 이번 재분석은 STAR `ReadsPerGene.out.tab` 원자료에서 HLA/AP/TLS target gene만 직접 추출했다. 분석 sample 구성은 Normal/ATC n=11; Normal/Follicular neoplasm n=11; Normal/Nodular hyperplasia n=7; Normal/PD n=3; Normal/PTC n=230; Tumor/ATC n=16; Tumor/PD n=5; Tumor/PTC n=349 이다.

Paper 2에는 이 결과를 **Korean thyroid expression generalization**으로 넣을 수 있다. 단, HT-specific 자료가 아니며 HLA allele/genotype 주장이 아니다.

## Primary contrast: PTC tumor vs PTC normal

| module | Cohen's d | delta score | Mann-Whitney p | FDR | AUC |
|---|---:|---:|---:|---:|---:|
| HLA_II_AP | 1.39 | 0.95 | 3.83e-41 | 4.03e-40 | 0.83 |
| Myeloid_DC | 1.23 | 0.81 | 7.95e-37 | 3.34e-36 | 0.81 |
| HLA_I | 1.22 | 0.87 | 1.45e-36 | 5.08e-36 | 0.81 |
| AP_TLS_composite | 0.89 | 0.58 | 1.47e-21 | 3.85e-21 | 0.73 |
| CD74_MIF_axis | 0.70 | 0.58 | 7.85e-15 | 1.37e-14 | 0.69 |
| T_IFNG | 0.68 | 0.49 | 7.27e-16 | 1.39e-15 | 0.70 |
| B_TLS | 0.28 | 0.21 | 3.95e-05 | 5.92e-05 | 0.60 |

## 추가 구조

Across all 632 samples, AP/TLS composite와 thyrocyte differentiation module의 Spearman rho=-0.28, p=8.33e-13. 이 값은 Korean cohort에서 immune/AP axis와 thyroid differentiation state가 어떻게 함께 움직이는지 보여주는 보조 결과다.

## 논문 반영 포인트

1. Paper 2 supplement 또는 external validation panel에 Korean bulk RNA-seq generalization으로 배치한다.
2. HT-overlap claim의 핵심은 GSE138198/GSE163203에 두고, GSE213647은 ancestry/geography-independent thyroid expression robustness로 쓴다.
3. `HLA allele`, `genotype`, `risk allele` 문구를 쓰지 않는다.

## 산출물

- Tables: `project/results/hla_two_paper_synthesis_2026_05_09/gse213647_korean_bulk_validation/tables/`
- Figures: `project/results/hla_two_paper_synthesis_2026_05_09/gse213647_korean_bulk_validation/figures/`
- Source data: `project/data/external/GSE213647/`
