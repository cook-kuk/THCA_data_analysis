# TCGA Preliminary Result Text (KR)

## 계산 내용

TCGA-THCA primary tumor 505개에서 RAI differentiation, HLA/APM immune visibility, CD8 exclusion, drug-delivery failure proxy, aggressive dedifferentiation 등 5개 치료취약성 축과 세부 module score를 계산했다.

## 핵심 결과

환자 수준 vulnerability label은 다음과 같이 분리되었다: {'RAI-readable differentiated': 127, 'Mixed/Other': 125, 'Drug-delivery barrier-high': 105, 'HLA-visible inflamed': 44, 'CD8-excluded myeloid/CAF-high': 42, 'HLA-low immune-invisible': 32, 'RAI-low dedifferentiated': 30}. Mixed/Other는 125/505명으로 남아 있으나 전체의 24.8% 수준이며, 제안서에서는 별도 세분화보다 non-dominant/intermediate group으로 두는 것이 안전하다.

주요 축 간 Spearman 상관은 다음과 같아 축들이 완전 중복되지 않음을 보인다.


|                                    |   rai_differentiation_score |   immune_visibility_score |   cd8_exclusion_proxy |   drug_delivery_failure_proxy |   aggressive_dedifferentiation_score |
|:-----------------------------------|----------------------------:|--------------------------:|----------------------:|------------------------------:|-------------------------------------:|
| rai_differentiation_score          |                        1    |                     -0.45 |                 -0.43 |                         -0.66 |                                -0.78 |
| immune_visibility_score            |                       -0.45 |                      1    |                  0.31 |                          0.52 |                                 0.49 |
| cd8_exclusion_proxy                |                       -0.43 |                      0.31 |                  1    |                          0.67 |                                 0.55 |
| drug_delivery_failure_proxy        |                       -0.66 |                      0.52 |                  0.67 |                          1    |                                 0.72 |
| aggressive_dedifferentiation_score |                       -0.78 |                      0.49 |                  0.55 |                          0.72 |                                 1    |

## 임상/driver 연관

탐색적 clinical association이 계산되었으며, FDR < 0.10인 상위 결과는 아래와 같다. 이는 예후모델 주장이 아니라 제안서용 hypothesis-generation 근거이다.


| axis                      | variable                    | test     |    effect |           p |         fdr |   n |
|:--------------------------|:----------------------------|:---------|----------:|------------:|------------:|----:|
| rai_differentiation_score | os_days                     | spearman | -0.105637 | 0.0176774   | 0.0820738   | 504 |
| rai_differentiation_score | OS.time                     | spearman | -0.102917 | 0.0207126   | 0.0897547   | 505 |
| rai_differentiation_score | ajcc_pathologic_tumor_stage | kruskal  | 27.8373   | 1.34567e-05 | 0.000145781 | 501 |
| rai_differentiation_score | stage                       | kruskal  | 29.6299   | 1.65091e-06 | 3.57697e-05 | 455 |
| immune_visibility_score   | ajcc_pathologic_tumor_stage | kruskal  | 28.602    | 9.41599e-06 | 0.000145781 | 501 |
| immune_visibility_score   | stage                       | kruskal  | 25.6427   | 1.13297e-05 | 0.000145781 | 455 |
| cd8_exclusion_proxy       | ajcc_pathologic_tumor_stage | kruskal  | 26.177    | 2.91481e-05 | 0.000236829 | 501 |
| cd8_exclusion_proxy       | stage                       | kruskal  | 21.1665   | 9.72171e-05 | 0.000631911 | 455 |

## 조심스러운 해석

TCGA bulk는 공간 heterogeneity를 증명하지 못하며, immune/stromal score는 세포조성 및 activation state가 섞인 bulk inference이다.

## 삼성 제안서 문장

TCGA-THCA 505명 공개 bulk transcriptome 분석에서 갑상선암은 단일 driver mutation 또는 평균 발현값으로 설명되지 않고, RAI 분화상태, HLA/APM 면역가시성, CD8 exclusion, stromal/drug-delivery barrier, proliferative dedifferentiation이 분리된 치료취약성 상태로 나뉘었다.

## Figure caption

TCGA-THCA public bulk pilot defines patient-level therapeutic vulnerability states across RAI differentiation, immune visibility, stromal exclusion, delivery-failure proxy, and aggressive dedifferentiation axes. Results are hypothesis-generating and require spatial and functional validation.
