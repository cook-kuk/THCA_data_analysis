# GSE184362/GSE191288 scRNA HLA/AP generalization

## 결론

두 dataset 모두 GEO sample metadata와 raw filename에서 HT/WHT/WOHT label을 확인할 수 없었다. 따라서 HT-specific evidence로 쓰면 안 된다.

대신 Paper 2 supplement에서 **low-tier scRNA expression generalization / boundary stress-test**로만 보관한다. GSE184362 구성은 Paratumor n=6, Tumor n=7, Metastasis n=10이고, GSE191288은 Tumor n=6, NonTumor n=1이다.

## GSE184362 sample-level contrasts

| contrast | module | Cohen's d | delta score | exact p | AUC |
|---|---|---:|---:|---:|---:|
| Metastasis_vs_Tumor | HLA_II_AP | -0.73 | -0.23 | 0.154 | 0.39 |
| Metastasis_vs_Tumor | AP_TLS_composite | -0.76 | -0.18 | 0.154 | 0.39 |
| Metastasis_vs_Tumor | HLA_I | nan | -0.70 | 0.500 | nan |
| Metastasis_vs_Tumor | B_TLS | nan | -0.21 | 1.000 | nan |
| Metastasis_vs_Tumor | T_IFNG | nan | -0.64 | 0.750 | nan |
| Metastasis_vs_Tumor | Myeloid_DC | nan | -1.61 | 0.750 | nan |
| Metastasis_vs_Tumor | CD74_MIF_axis | nan | -0.65 | 0.500 | nan |
| Tumor_vs_Paratumor | T_IFNG | 0.65 | 0.29 | 0.714 | nan |
| Tumor_vs_Paratumor | Myeloid_DC | 0.58 | 0.71 | 1.000 | nan |
| Tumor_vs_Paratumor | CD74_MIF_axis | 0.30 | 0.08 | 1.000 | nan |
| Tumor_vs_Paratumor | HLA_I | 0.21 | 0.03 | 1.000 | nan |
| Tumor_vs_Paratumor | HLA_II_AP | 0.17 | 0.08 | 0.731 | 0.50 |
| Tumor_vs_Paratumor | AP_TLS_composite | 0.03 | 0.01 | 1.000 | 0.48 |
| Tumor_vs_Paratumor | B_TLS | -0.37 | -0.10 | 1.000 | nan |

## GSE184362 paired tumor-paratumor

| module | n pairs | mean tumor-paratumor delta | Wilcoxon p |
|---|---:|---:|---:|
| HLA_II_AP | 6 | 0.12 | 0.500 |
| AP_TLS_composite | 6 | 0.04 | 1.000 |
| HLA_I | 6 | nan | nan |
| B_TLS | 6 | nan | nan |
| T_IFNG | 6 | nan | nan |
| Myeloid_DC | 6 | nan | nan |
| CD74_MIF_axis | 6 | nan | nan |

## GSE191288 qualitative only

GSE191288은 non-tumor가 1개뿐이어서 p-value를 해석하지 않는다.

| module | n tumor | n non-tumor | tumor - non-tumor delta |
|---|---:|---:|---:|
| Myeloid_DC | 6 | 1 | 1.36 |
| T_IFNG | 6 | 1 | 1.09 |
| HLA_I | 6 | 1 | 1.02 |
| HLA_II_AP | 6 | 1 | 1.01 |
| AP_TLS_composite | 6 | 1 | 0.92 |
| B_TLS | 6 | 1 | 0.83 |
| CD74_MIF_axis | 6 | 1 | -0.41 |

## 논문 반영 포인트

1. Main claim에는 쓰지 않는다.
2. Supplement registry에는 `HT label unavailable; expression-only scRNA generalization`으로 둔다.
3. `HLA allele`, `genotype`, `risk allele`, `HT-specific` 표현을 쓰지 않는다.

## 산출물

- Tables: `project/results/hla_two_paper_synthesis_2026_05_09/gse184362_gse191288_scrna_generalization/tables/`
- Figures: `project/results/hla_two_paper_synthesis_2026_05_09/gse184362_gse191288_scrna_generalization/figures/`
- Source data: `project/data/external/GSE184362/`, `project/data/external/GSE191288/`
