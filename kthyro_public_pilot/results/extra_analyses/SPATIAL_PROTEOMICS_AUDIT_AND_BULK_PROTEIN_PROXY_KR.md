# K-Thyro 추가 분석: 공개 spatial proteomics 감사 + bulk protein proxy 검증

작성일: 2026-05-09  
작업 디렉토리: `/home/seungho/personal/THCA_data_analysis/kthyro_public_pilot`

## 한 줄 결론

현재 즉시 재분석 가능한 **thyroid cancer public raw spatial proteomics matrix**는 확인되지 않았다. 대신 로컬에 이미 정리된 461개 thyroid proteomics 샘플을 이용해 K-Thyro 치료취약성 축을 protein level에서 검산했으며, **RAI differentiation protein axis는 ATC에서 강하게 감소하고, myeloid/TGFB barrier protein axis는 ATC에서 강하게 증가**했다.

따라서 제안서 표현은 다음이 가장 안전하다.

> Public spatial transcriptomics에서 치료취약성 niche의 공간 구조가 관찰되었고, 독립 bulk proteomics에서 핵심 축의 단백질 방향성이 재현되었다. 다만 spatial proteomics는 본 과제에서 FFPE-mIHC/GeoMx protein panel로 직접 구축해야 할 validation layer이다.

## 분석 입력

### 1. Public spatial proteomics audit

출력 파일:

- `results/extra_analyses/tables/spatial_proteomics_public_resource_audit.tsv`
- `results/extra_analyses/tables/literature_iHC_spatial_marker_support_candidates.tsv`

감사 결과:

| 리소스 | 판정 | K-Thyro에서의 사용 |
|---|---|---|
| GSE301163 | GeoMx DSP spatial RNA matrix. Spatial은 맞지만 protein matrix는 아님. | spatial ROI transcriptomics 보조 근거 |
| Donati et al. Virchows Archiv 2025 | thyroid anaplastic carcinoma에서 DSP feasibility를 보인 문헌. Public count matrix는 확인 안 됨. | GeoMx/DSP feasibility 근거 |
| Haq/Bychkov/Mete/Jeon/Jung Endocrine Pathology 2025 | spatial transcriptomics + IHC validation 문헌. Raw spatial protein matrix는 아님. | extended marker panel 후보 |
| Local proteogenomic_v1 protein table | bulk proteomics, thyroid-specific, 461 samples. Spatial은 아님. | orthogonal protein-level validation |

### 2. Bulk proteomics proxy

입력 파일:

- `/home/seungho/personal/THCA_data_analysis/project/results/proteogenomic_v1/paper3_mun2025_dediff_layer/protein_target_genes_S1C.tsv`

출력 파일:

- `results/extra_analyses/tables/bulk_proteomics_kthyro_module_scores.tsv`
- `results/extra_analyses/tables/bulk_proteomics_kthyro_module_coverage.tsv`
- `results/extra_analyses/tables/bulk_proteomics_kthyro_module_contrasts.tsv`
- `results/extra_analyses/tables/bulk_proteomics_kthyro_dediff_trends.tsv`
- `results/extra_analyses/figures/bulk_proteomics_kthyro_axis_boxplots.png`
- `results/extra_analyses/figures/bulk_proteomics_kthyro_contrast_heatmap.png`

샘플 수: 461개.

사용 가능한 protein module coverage:

| Protein module | Found proteins |
|---|---|
| RAI differentiation protein | SLC5A5, TPO, TG, TSHR, PAX8, FOXE1, DIO1 |
| HLA-I/APM protein | B2M, TAP1, TAP2, PSMB8, PSMB9 |
| HLA-II/APC protein | CD74 |
| CD8 cytotoxic protein | CD8A, PRF1 |
| IFN/APM activation protein | STAT1, IRF1, PSMB8, PSMB9 |
| Checkpoint/suppression protein | TOX, LAG3, CD274, PDCD1LG2 |
| Myeloid suppressive protein | ITGAM, MRC1, CD163, TGFB1 |
| TLS/B-cell protein | CXCL13, CCL21, JCHAIN, MS4A1, CD79A, MZB1 |
| TGFB barrier single protein | TGFB1 |

## 핵심 결과

### ATC vs PTC protein contrast

출처: `bulk_proteomics_kthyro_module_contrasts.tsv`

| Axis | ATC mean | PTC mean | Cohen's d | FDR q | 해석 |
|---|---:|---:|---:|---:|---|
| RAI differentiation protein | -0.697 | 0.034 | -1.989 | 4.23e-35 | ATC에서 RAI/thyroid differentiation protein 축 강한 저하 |
| Myeloid suppressive protein | 0.797 | -0.288 | 1.679 | 1.99e-23 | ATC에서 myeloid/TGFB 관련 protein 축 강한 상승 |
| TGFB barrier single protein | 0.889 | -0.188 | 1.172 | 1.81e-19 | ATC에서 TGFB barrier protein 상승 |
| HLA-I/APM protein | 0.473 | 0.004 | 0.643 | 1.91e-05 | ATC에서 HLA/APM protein 축 상승 방향 |
| IFN/APM activation protein | 0.440 | 0.068 | 0.585 | 5.44e-05 | ATC에서 IFN/APM protein activation 방향 |
| Checkpoint/suppression protein | 0.180 | -0.080 | 0.495 | 1.86e-02 | ATC에서 checkpoint/suppression protein 축 상승 방향 |
| CD8 cytotoxic protein | 0.111 | -0.046 | 0.212 | 0.280 | 약함 |
| HLA-II/APC protein | 0.326 | 0.064 | 0.303 | 0.0698 | 경계선 |
| TLS/B-cell protein | 0.018 | -0.022 | 0.079 | 0.902 | 지지 약함 |

### Dedifferentiation trend across PTC → PDTC → ATC

출처: `bulk_proteomics_kthyro_dediff_trends.tsv`

| Axis | Spearman rho | FDR q | 방향 |
|---|---:|---:|---|
| RAI differentiation protein | -0.420 | 4.01e-20 | dedifferentiation과 함께 감소 |
| Myeloid suppressive protein | 0.408 | 3.36e-19 | dedifferentiation과 함께 증가 |
| TGFB barrier single protein | 0.342 | 1.33e-13 | dedifferentiation과 함께 증가 |
| Protein aggressive visibility score | 0.323 | 2.79e-12 | dedifferentiation과 함께 증가 |
| HLA-I/APM protein | 0.093 | 0.0909 | 약한 증가, FDR 기준 경계 |

## 제안서에서 바로 쓸 메시지

### Strong message

> 공개 spatial transcriptomics에서 치료취약성 niche가 공간적으로 응집되고, 독립 thyroid proteomics 461개 샘플에서 RAI differentiation 단백질 축 저하와 myeloid/TGFB barrier 단백질 축 상승이 재현되었다. 이는 K-Thyro의 핵심 가설이 RNA 발현 score에만 의존하지 않고 protein-level validation으로 확장 가능함을 보여준다.

### 더 조심해야 하는 문장

> 공개 spatial proteomics로 검증되었다.

이 문장은 아직 쓰면 안 된다. 현재 확인된 것은 `bulk proteomics proxy`이지 `spatial protein niche map`이 아니다.

## Allowed claim / Forbidden claim

### Allowed

- Public thyroid spatial proteomics raw matrix는 아직 확보되지 않았다.
- Public spatial transcriptomics는 spatial niche coherence를 지지한다.
- Local thyroid bulk proteomics proxy는 RAI protein loss와 myeloid/TGFB barrier protein gain을 지지한다.
- Spatial proteomics/mIHC는 본 과제에서 새로 구축해야 할 핵심 validation layer다.

### Forbidden

- Public spatial proteomics에서 K-Thyro niche가 검증되었다.
- Protein-level spatial localization이 이미 입증되었다.
- HLA/APM protein score가 peptide presentation 또는 immunotherapy response를 직접 증명한다.
- TGFB/myeloid protein increase가 drug-delivery failure를 직접 증명한다.

## 삼성 제안서에 넣을 위치

추천 위치:

1. Public-data pilot GO decision 뒤의 보강 박스
2. Spatial niche coherence 장표의 하단 claim-boundary 박스
3. Experimental validation plan 장표에서 “왜 mIHC/GeoMx protein layer가 필요한가” 근거
4. 30억/3년 execution logic에서 “public gap → 우리가 생성할 데이터” 논리

추천 장표 문장:

> 공개 데이터는 spatial RNA niche와 bulk protein 방향성을 지지하지만, 치료결정에 필요한 protein-level spatial localization은 아직 비어 있다. K-Thyro는 FFPE-mIHC/GeoMx protein panel과 fresh tissue functional assay로 이 gap을 직접 메운다.

## Figure 설명

### `bulk_proteomics_kthyro_axis_boxplots.png`

내용:

- PTC/PDTC/ATC 그룹별 protein module score boxplot.
- RAI differentiation protein은 PTC에서 높고 ATC에서 낮다.
- Myeloid suppressive/TGFB barrier protein은 ATC에서 높다.

Caption:

> Independent thyroid bulk proteomics supports protein-level directionality of key K-Thyro axes, but does not provide spatial localization.

### `bulk_proteomics_kthyro_contrast_heatmap.png`

내용:

- ATC vs PTC, PDTC vs PTC, ATC vs PDTC 등 contrast별 Cohen's d heatmap.
- 가장 강한 방향성은 ATC vs PTC에서 RAI protein loss와 myeloid/TGFB barrier protein gain.

Caption:

> Protein contrasts provide orthogonal support for dedifferentiation and stromal-suppressive axes; spatial protein mapping remains a required validation aim.

## 참고한 공개 소스

- GSE301163 / GEO metadata: thyroid spatial transcriptomics dataset, not a protein matrix.
- Donati et al., Virchows Archiv 2025, “Digital spatial profiling for pathologists”: GeoMx DSP feasibility including thyroid anaplastic carcinoma.
- Haq et al., Endocrine Pathology 2025, “Identification of Specific Biomarkers for Anaplastic Thyroid Carcinoma Through Spatial Transcriptomic and Immunohistochemical Profiling”: spatial transcriptomics plus IHC marker validation in ATC.

