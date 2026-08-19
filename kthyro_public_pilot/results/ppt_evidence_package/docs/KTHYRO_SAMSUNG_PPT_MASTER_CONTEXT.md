# K-Thyro Samsung PPT Master Context

작업 위치: `/home/seungho/personal/THCA_data_analysis/kthyro_public_pilot`

목적: 삼성미래기술육성사업 제안서/PPT에 넣을 **public-data pilot evidence package**의 전체 상황 정리.  
기준: 로컬에서 실제 생성된 `results/` 산출물 기반. 새 임상 주장이나 drug overclaim 금지.

---

## 0. 한 줄 결론

**GO.** 공개데이터 pilot 결과, 갑상선암은 평균 driver mutation만으로 치료취약성을 설명하기 어렵고, **RAI 분화상태, HLA/APM 면역가시성, CD8 exclusion, myeloid/CAF barrier, drug-delivery failure proxy**가 환자 수준과 조직 공간 수준에서 분리된 치료취약성 축으로 나타난다.

PPT 앞 메시지는 drug가 아니다.

> **TCGA 505명 + GSE250521 spatial 16 slide에서 치료취약성 niche가 분리되고, 그 niche가 공간적으로 응집된다.**

Drug 파트 문장:

> **Public drug resources nominate perturbation classes for validation, not final therapies.**

Clinical boundary:

> **Public-data pilot supports hypothesis generation and validation design, not clinical deployment.**

---

## 1. 핵심 숫자

### TCGA-THCA bulk

- Public TCGA-THCA primary tumor: **505명**
- Vulnerability label 분포:
  - RAI-readable differentiated: **127**
  - Mixed/Other: **125**
  - Drug-delivery barrier-high: **105**
  - HLA-visible inflamed: **44**
  - CD8-excluded myeloid/CAF-high: **42**
  - HLA-low immune-invisible: **32**
  - RAI-low dedifferentiated: **30**
- Driver group:
  - BRAF: **274**
  - RAS: **54**
  - No BRAF/RAS call: **177**

### Spatial transcriptomics

- Dataset: **GSE250521**
- Slides: **16**
- Spots: **57,144**
- Conditions:
  - normal: **4**
  - PTC: **4**
  - locally advanced PTC: **4**
  - ATC: **4**
- KNN same-niche coherence:
  - **16/16 slides z-score > 2**
  - z-score min/median/max: **6.09 / 22.00 / 43.16**

### External validation expansion

- Public validation matrix v2: **301 tests**
  - strong: **144**
  - moderate: **30**
  - weak: **73**
  - weak_or_none: **52**
  - opposite: **2**
- Exact K-Thyro external GPL570 rescoring:
  - Cohorts: GSE29265, GSE33630, GSE53157, GSE65144
  - Samples: **206**
  - Contrast tests: **168**
  - strong: **65**
  - moderate: **12**
  - prespecified direction: **42 match / 2 opposite**

### Axis separability

- TCGA patient-level axis-pair correlations: median abs rho **0.609**
- External GPL570 axis-pair correlations: median abs rho **0.637**
- Spatial slide-level spot correlations: mostly much lower; no slide had axis-pair abs rho >= 0.75.
- Interpretation: axes are biologically related, but not collapsed into one signal.

### Mutation-only insufficiency

Driver mutation is related to vulnerability labels, but not enough.

- BRAF group: **274 patients**, split into **7 vulnerability labels**
  - largest label: Drug-delivery barrier-high **91/274 = 33.2%**
  - mutation-only ambiguity fraction: **66.8%**
- No BRAF/RAS call group: **177 patients**, split into **7 labels**
  - largest label: RAI-readable differentiated **76/177 = 42.9%**
  - ambiguity fraction: **57.1%**
- RAS group: **54 patients**, split into **3 labels**
  - largest label: RAI-readable differentiated **42/54 = 77.8%**
  - ambiguity fraction: **22.2%**
- Driver-label Cramer’s V: **0.464**, chi-square p **8.05e-40**
  - Interpretation: driver association exists, but does not determine therapeutic vulnerability state.

Driver group variance explained:

| Axis | eta² by driver/no-call group | Unexplained fraction |
|---|---:|---:|
| Drug-delivery failure proxy | 0.332 | 0.668 |
| RAI differentiation | 0.309 | 0.691 |
| Aggressive dedifferentiation | 0.273 | 0.727 |
| Myeloid/CAF barrier | 0.182 | 0.818 |
| CD8 exclusion | 0.162 | 0.838 |
| HLA-I/APM | 0.156 | 0.844 |
| Immune visibility | 0.125 | 0.875 |
| Proliferation | 0.107 | 0.893 |

---

## 2. 제안서 핵심 문장

### 한국어 핵심 문장

> 공개데이터 pilot은 갑상선암 치료반응이 단일 driver mutation이나 평균 bulk expression으로 설명되지 않으며, RAI 분화상태, HLA/APM 면역가시성, CD8 exclusion, myeloid/CAF barrier, drug-delivery failure proxy가 환자 및 조직 공간 수준에서 분리된 치료취약성 축으로 나타남을 보여준다. 특히 GSE250521 공간전사체 16개 slide에서 치료취약성 niche가 비무작위적으로 응집되어, 본 과제의 핵심 가설인 ‘갑상선암 공간 치료취약성 생태계’가 공개데이터 수준에서 지지되었다.

### 삼성육성과제용 메시지

> 본 과제는 갑상선암을 평균적인 driver mutation 질환이 아니라, 수술 전후 perturbation과 공간오믹스·병리·기능실험으로 읽는 치료취약성 생태계로 재정의한다.

### 영어 technical caption

> Public-data pilot analyses support separable therapeutic vulnerability axes in thyroid cancer across TCGA bulk expression, public spatial transcriptomics, external thyroid cohorts, single-cell references, and public perturbation resources. These results support hypothesis generation and experimental validation design, not clinical deployment.

---

## 3. PPT 10장 구성

### Slide 1. Public-data pilot GO decision

**Title:** 공개데이터 pilot 결과: K-Thyro는 GO

**One-line message:**  
TCGA 505명과 GSE250521 spatial 16 slide에서 치료취약성 축이 분리되고, 공간 niche가 비무작위적으로 응집된다.

**Bullets:**
- TCGA-THCA 505명에서 7개 vulnerability label 도출
- GSE250521 16 slide, 57,144 spots에서 치료취약성 niche mapping
- 16/16 slide에서 KNN same-niche coherence z-score > 2
- External public validation matrix v2: 301 tests, strong/moderate 174 tests
- 최종 판단: **GO**, 단 public-data pilot은 hypothesis generation 단계

**Figure to use:**
- `kthyro_public_pilot/results/figures/proposal/figure8_pilot_conclusion_slide.png`
- 보강용: `kthyro_public_pilot/results/public_validation_expansion/figures/public_validation_expansion_v2_support_counts.png`

**Caption:**  
공개데이터 pilot은 치료취약성 축의 존재와 공간 응집성을 지지하지만, 임상 예측 모델이나 치료효과 검증을 의미하지 않는다.

---

### Slide 2. Data sources and analysis flow

**Title:** 공개데이터 기반 K-Thyro pilot 분석 흐름

**One-line message:**  
Bulk, spatial, scRNA, GeoMx, external cohort, drug resource를 하나의 치료취약성 scoring framework로 통합했다.

**Bullets:**
- TCGA-THCA: patient-level vulnerability scoring
- GSE250521 spatial: spot-level niche mapping and slide-level coherence
- External GPL570 cohorts: exact K-Thyro gene set reproducibility check
- scRNA/GeoMx: cell-type and ROI-level sanity check
- DepMap/PRISM: perturbation class hypothesis nomination

**Figure to use:**
- `kthyro_public_pilot/results/figures/proposal/figure1_public_pilot_data_flow_dark.png`
- 또는 `kthyro_public_pilot/results/figures/proposal/figure1_public_pilot_data_flow.png`

**Caption:**  
데이터 종류별 resolution과 claim boundary가 다르므로, 각 layer는 독립적 검증이 아니라 상보적 근거로 해석한다.

---

### Slide 3. TCGA 505-patient vulnerability landscape

**Title:** TCGA 505명에서 치료취약성 축이 환자별로 분리된다

**One-line message:**  
평균 bulk mutation만으로는 RAI, immune visibility, barrier, proliferation 상태를 모두 설명할 수 없다.

**Bullets:**
- TCGA primary tumor 505명에서 K-Thyro gene set score 계산
- RAI differentiation, HLA/APM, CD8/cytotoxic, myeloid/CAF, drug-delivery failure proxy, aggressive dedifferentiation 축 생성
- 환자별 vulnerability label이 RAI-readable, HLA-visible, HLA-low, CD8-excluded, barrier-high 등으로 분화
- Driver mutation은 관련되지만 label을 결정하지 않음

**Figure to use:**
- `kthyro_public_pilot/results/figures/proposal/tcga_vulnerability_axes_heatmap_dark.png`
- 또는 `kthyro_public_pilot/results/figures/proposal/figure2_tcga_therapeutic_vulnerability_landscape.png`

**Caption:**  
TCGA bulk는 환자 평균 발현 기반이므로 공간 이질성을 직접 증명하지 않는다. 다만 치료취약성 축의 patient-level separability를 지지한다.

---

### Slide 4. Patient-level subtype counts and mutation-only insufficiency

**Title:** Driver mutation만으로 치료취약성 상태를 결정할 수 없다

**One-line message:**  
BRAF group 내부에서도 7개 vulnerability label로 갈라져 mutation-only 전략의 한계를 보인다.

**Bullets:**
- BRAF group 274명: 7개 vulnerability label로 분산
- BRAF 최대 label은 Drug-delivery barrier-high 91명, 33.2%
- 즉 BRAF만 알면 66.8%는 vulnerability state가 ambiguous
- HLA-I/APM axis의 driver 설명력 eta² 0.156, immune visibility eta² 0.125
- Mutation + spatial/immune/barrier state를 함께 읽는 companion strategy 필요

**Figure to use:**
- `kthyro_public_pilot/results/public_validation_expansion/figures/tcga_mutation_not_enough_driver_label_fraction.png`
- `kthyro_public_pilot/results/public_validation_expansion/figures/tcga_mutation_not_enough_axis_variance.png`
- `kthyro_public_pilot/results/figures/proposal/tcga_label_distribution_dark.png`

**Caption:**  
Driver mutation과 vulnerability label은 유의하게 관련되지만, driver group 내부 heterogeneity가 커서 mutation-only stratification은 치료취약성 생태계를 충분히 설명하지 못한다.

---

### Slide 5. Spatial ST 16-slide niche coherence

**Title:** 치료취약성 niche는 spatial spot noise가 아니라 조직 내 구조다

**One-line message:**  
GSE250521 16개 slide 모두에서 같은 niche가 공간적으로 이웃하는 non-random coherence를 보였다.

**Bullets:**
- GSE250521: normal 4, PTC 4, locally advanced PTC 4, ATC 4
- 총 57,144 spots에 동일한 K-Thyro vulnerability score 적용
- KNN same-niche permutation: 16/16 slide z-score > 2
- Coherence z-score min/median/max: 6.09 / 22.00 / 43.16
- Spots는 biological replicate가 아니며 slide-level summary로 해석

**Figure to use:**
- `kthyro_public_pilot/results/figures/proposal/spatial_coherence_barplot_dark.png`
- `kthyro_public_pilot/results/figures/proposal/figure4_spatial_coherence_statistics.png`

**Caption:**  
Spatial coherence는 치료취약성 expression state가 조직 내에서 응집됨을 지지하지만, 임상 반응 예측이나 치료효과를 증명하지 않는다.

---

### Slide 6. Representative spatial niche maps

**Title:** RAI, HLA/APM, CAF/myeloid barrier niche가 공간적으로 분리된다

**One-line message:**  
같은 조직 안에서도 RAI-low, HLA-visible, CD8-excluded, barrier-high territory가 서로 다른 위치에 나타난다.

**Bullets:**
- Spot-level RAI differentiation, HLA/APM visibility, cytotoxicity, CAF/myeloid barrier score 계산
- Quantile threshold 기반 exploratory niche label 부여
- Representative slide에서 치료취약성 niche map 시각화
- Spatial map은 FFPE/mIHC, GeoMx/ROI, functional assay의 ROI 선정 근거

**Figure to use:**
- `kthyro_public_pilot/results/figures/proposal/spatial_representative_maps_dark.png`
- `kthyro_public_pilot/results/figures/proposal/figure3_spatial_therapeutic_niche_maps.png`
- 추가 montage: `kthyro_public_pilot/results/figures/spatial_niche_montage.png`

**Caption:**  
Spatial transcriptomics는 mRNA expression territory를 보여준다. HLA peptide presentation, iodine uptake, drug penetration은 별도 단백질/기능실험으로 검증해야 한다.

---

### Slide 7. Integrated evidence matrix

**Title:** 5개 치료취약성 축은 다층 공개데이터에서 상보적으로 지지된다

**One-line message:**  
TCGA, spatial, external bulk, GeoMx/scRNA, perturbation resource가 서로 다른 claim boundary 안에서 같은 platform hypothesis를 지지한다.

**Bullets:**
- Public validation matrix v2: 301 tests
- strong 144, moderate 30
- External GPL570 exact K-Thyro rescoring: 206 samples, 168 tests
- Axis separability: TCGA/external/spatial에서 축들이 하나의 score로 붕괴되지 않음
- Integrated evidence는 clinical deployment가 아니라 validation design의 근거

**Figure to use:**
- `kthyro_public_pilot/results/figures/proposal/integrated_public_pilot_evidence_matrix_dark.png`
- `kthyro_public_pilot/results/public_validation_expansion/figures/public_validation_expansion_v2_support_counts.png`
- `kthyro_public_pilot/results/public_validation_expansion/figures/axis_separability_correlation_heatmaps.png`

**Caption:**  
Evidence matrix는 각 데이터 layer의 관찰/inference/가설을 분리해 보여주며, 모든 축은 prospective tissue validation이 필요하다.

---

### Slide 8. Drug/perturbation class hypotheses

**Title:** Drug 후보가 아니라 perturbation class를 검증한다

**One-line message:**  
Public drug resources는 final therapy가 아니라 RAI/HLA/barrier/drug-delivery niche를 재프로그래밍할 실험 후보군을 제안한다.

**Bullets:**
- MAPK-axis modulation: RAI redifferentiation hypothesis
- Epigenetic/IFN/APM modulation: HLA/APM visibility restoration hypothesis
- Cell-cycle stress vulnerability: aggressive/proliferative niche hypothesis
- Myeloid/CAF barrier modulation: CD8 exclusion 해소 후보축
- Drug-delivery strategy: fluorescent drug/liposome/nanoparticle penetration imaging으로 검증

**Figure to use:**
- `kthyro_public_pilot/results/figures/proposal/figure6_drug_perturbation_candidate_pilot.png`
- `kthyro_public_pilot/results/figures/drug_pilot_candidate_heatmap.png`

**Caption:**  
DepMap/PRISM는 cell-line 또는 pan-cancer resource이므로 myeloid/CAF barrier와 drug-delivery failure를 직접 검증하지 못한다. Drug section은 class-level hypothesis로만 제시한다.

---

### Slide 9. Experimental validation plan

**Title:** AI-spatial inference를 FFPE/GeoMx/기능실험으로 닫는다

**One-line message:**  
Public pilot에서 나온 niche를 실제 환자 조직과 perturbation assay로 검증하는 폐루프 플랫폼을 구축한다.

**Bullets:**
- FFPE cohort: spatial marker validation and ROI selection
- mIHC/IF core panel: PAX8/TG/TPO/NIS, HLA-I/B2M/TAP1, CD8/GZMB, CD68/CD163, ACTA2/FAP/COL1A1, PD-L1
- GeoMx/ROI: HLA-low, CD8-excluded, RAI-low, CAF/myeloid-rich territory 검증
- Fresh tissue organoid/slice: iodide uptake, HLA/APM rescue, PBMC/TIL co-culture
- Fluorescent drug/liposome/nanoparticle imaging: delivery failure proxy 검증

**Figure to use:**
- `kthyro_public_pilot/results/figures/proposal/figure7_samsung_experimental_validation_design.png`

**Caption:**  
Public-data score는 치료취약성 가설과 ROI를 제안한다. 임상 적용 가능성은 단백질 검증, 기능실험, prospective cohort로 단계적으로 확인해야 한다.

---

### Slide 10. Samsung 30억/3년 execution logic

**Title:** 30억/3년: 공간 치료취약성 생태계를 실험으로 증명한다

**One-line message:**  
공개데이터 pilot은 GO를 만들었고, 본 과제는 이를 실제 수술 전후 perturbation atlas와 기능검증 플랫폼으로 확장한다.

**Bullets:**
- Year 1: FFPE/mIHC + GeoMx/ROI로 niche atlas 구축
- Year 2: fresh tissue slice/organoid perturbation으로 RAI/HLA/barrier rescue 검증
- Year 3: AI-spatial companion model + theranostic validation package 완성
- 30억/3년은 clinical deployment가 아니라 high-risk/high-impact PoC 구축 예산
- 50억/5년 확장 시 prospective cohort와 interventional validation으로 Nature급 proof 가능

**Figure to use:**
- `kthyro_public_pilot/results/figures/proposal/figure1_public_pilot_data_flow_dark.png`
- `kthyro_public_pilot/results/figures/proposal/figure7_samsung_experimental_validation_design.png`
- `kthyro_public_pilot/results/figures/proposal/figure8_pilot_conclusion_slide.png`

**Caption:**  
삼성 과제의 핵심 산출물은 drug list가 아니라, 수술 전후 perturbation과 공간오믹스·병리·기능실험을 연결하는 AI-spatial theranostic validation platform이다.

---

## 4. Allowed Claim / Forbidden Claim

### TCGA section

**Allowed claim**

- TCGA-THCA 505명에서 K-Thyro therapeutic vulnerability axes가 환자별로 분리된다.
- Driver mutation group은 vulnerability label과 관련되지만, label을 완전히 결정하지 않는다.
- BRAF group 내부에서도 7개 vulnerability label이 존재한다.
- Bulk expression 기반 결과는 patient-level hypothesis generation에 적합하다.

**Forbidden claim**

- TCGA bulk만으로 spatial heterogeneity를 증명했다.
- TCGA score가 실제 RAI 반응, immunotherapy 반응, drug-delivery 실패를 예측한다.
- BRAF/RAS mutation이 무의미하다.
- TCGA score만으로 환자 치료전략을 결정할 수 있다.

### Spatial section

**Allowed claim**

- GSE250521 16 slide에서 K-Thyro niche label이 non-random spatial coherence를 보인다.
- Spatial transcriptomics는 치료취약성 expression territory의 존재를 지지한다.
- Spatial map은 FFPE/mIHC, GeoMx/ROI, fresh tissue assay의 ROI 설계 근거가 된다.

**Forbidden claim**

- Spatial transcriptomics만으로 clinical response를 예측했다.
- HLA/APM mRNA score가 실제 peptide presentation을 증명한다.
- RAI differentiation mRNA score가 iodine uptake를 증명한다.
- Drug-delivery failure proxy가 실제 약물 침투 실패를 증명한다.
- Spots를 독립 biological replicate처럼 사용했다.

### Drug/perturbation section

**Allowed claim**

- Public drug resources nominate perturbation classes for validation.
- MAPK, epigenetic/APM, cell-cycle, myeloid/CAF, delivery-enhancing strategy는 실험 가설이다.
- DepMap/PRISM는 thyroid line sensitivity와 pan-cancer association을 탐색하는 resource다.

**Forbidden claim**

- 특정 drug가 갑상선암 환자에서 효과적이라고 주장한다.
- JAK inhibitor가 HLA/APM을 회복한다고 주장한다.
- Mycophenolic acid를 immune-restoration therapy로 주장한다.
- Cell-line drug sensitivity로 myeloid/CAF barrier 또는 drug-delivery failure를 검증했다고 말한다.

---

## 5. 가장 강한 근거 7개

1. **TCGA 505명에서 vulnerability label이 분리됨**  
   RAI-readable, HLA-visible, HLA-low, CD8-excluded, barrier-high, RAI-low 상태가 환자별로 나뉜다.

2. **Spatial 16/16 slide에서 non-random niche coherence**  
   KNN same-niche z-score가 모든 slide에서 > 2이며, median z-score가 22.00이다.

3. **Mutation-only 설명 부족**  
   BRAF 274명 내부가 7개 label로 갈라지고, BRAF 최대 label도 33.2%에 불과하다.

4. **External exact K-Thyro validation**  
   GSE29265/GSE33630/GSE53157/GSE65144 206샘플에서 168 tests, strong/moderate 77 tests.

5. **Public validation matrix v2**  
   전체 301 tests 중 strong 144, moderate 30.

6. **Axis separability**  
   TCGA/external/spatial에서 축들이 관련은 있지만 하나의 score로 붕괴되지 않음.

7. **실험가설로 바로 연결됨**  
   RAI-low/MAPK-high, HLA/APM-low, CD8-excluded, CAF/myeloid-rich, hypoxia/vascular-low territory가 각각 검증 assay로 연결된다.

---

## 6. 가장 약한 부분과 방어 문장

### Weak point 1. Direct RAI avid/refractory validation is weak

**Issue:** GSE151179 direct RAI avid/refractory bulk signal은 강한 validation으로 쓰기 어렵다.

**Defense wording:**  
RAI score는 direct clinical predictor가 아니라 differentiation/restorability hypothesis로 사용하며, iodide uptake assay와 redifferentiation perturbation으로 검증한다.

### Weak point 2. Drug candidates are not final therapies

**Issue:** DepMap/PRISM는 drug-delivery, myeloid/CAF, immune restoration을 직접 검증하지 못한다.

**Defense wording:**  
Public drug resources nominate perturbation classes, not final therapies. 후보는 organoid/slice, mIHC, HLA rescue, nanoparticle penetration assay에서 검증한다.

### Weak point 3. HLA/APM mRNA is not peptide presentation

**Issue:** mRNA HLA/APM score는 실제 antigen presentation이나 patient-specific neoantigen presentation을 증명하지 않는다.

**Defense wording:**  
HLA/APM-low territory는 HLA-I/B2M/TAP1 protein staining, GeoMx, IFN/epigenetic rescue assay로 검증한다.

### Weak point 4. Spatial spots are not independent biological samples

**Issue:** spot 수가 많아도 biological n은 slide/patient level이다.

**Defense wording:**  
공간 통계는 slide별로 계산하고, group comparison은 slide-level summary로 해석한다.

### Weak point 5. Drug-delivery failure is proxy

**Issue:** transcriptomic barrier/hypoxia/vascular score는 실제 약물 침투를 측정하지 않는다.

**Defense wording:**  
Drug-delivery failure score는 fluorescent drug/liposome/nanoparticle distribution imaging으로 검증할 ROI nomination tool이다.

---

## 7. Experimental validation panel

### Core panel

| Axis | Markers |
|---|---|
| Tumor/thyroid/RAI | PAX8, TG, TPO, NIS/SLC5A5 |
| HLA/APM | HLA-I/HLA-ABC, B2M, TAP1 |
| T cell/cytotoxic | CD3, CD8, GZMB |
| Myeloid/TAM | CD68, CD163 |
| CAF/ECM | ACTA2/alpha-SMA, FAP, COL1A1 |
| Checkpoint | PD-L1 |
| Hypoxia/delivery proxy | CA9, VEGFA |

### Extended panel

| Axis | Markers |
|---|---|
| Thyroid differentiation | TSHR, FOXE1, NKX2-1 |
| APM regulation | NLRC5, PSMB8, PSMB9 |
| Myeloid state | MRC1, SPP1 |
| CAF/ECM barrier | POSTN, FN1 |
| Hypoxia/metabolism | GLUT1/SLC2A1 |
| Proliferation | Ki-67 |

### Functional assays

- Iodide uptake assay after MAPK-axis perturbation
- HLA-I/B2M/TAP1 rescue assay after IFN/epigenetic perturbation
- PBMC/TIL co-culture if feasible
- Fresh tissue slice/organoid perturbation
- Fluorescent drug/liposome/nanoparticle penetration imaging
- GeoMx/ROI validation for HLA-low, CD8-excluded, RAI-low, CAF/myeloid-rich territories

---

## 8. Figure path quick list

### Core proposal figures

- Data flow: `kthyro_public_pilot/results/figures/proposal/figure1_public_pilot_data_flow_dark.png`
- TCGA landscape: `kthyro_public_pilot/results/figures/proposal/figure2_tcga_therapeutic_vulnerability_landscape.png`
- TCGA heatmap dark: `kthyro_public_pilot/results/figures/proposal/tcga_vulnerability_axes_heatmap_dark.png`
- TCGA label distribution dark: `kthyro_public_pilot/results/figures/proposal/tcga_label_distribution_dark.png`
- TCGA quadrant dark: `kthyro_public_pilot/results/figures/proposal/tcga_therapeutic_quadrant_dark.png`
- Spatial maps: `kthyro_public_pilot/results/figures/proposal/spatial_representative_maps_dark.png`
- Spatial coherence: `kthyro_public_pilot/results/figures/proposal/spatial_coherence_barplot_dark.png`
- Spatial niche fraction: `kthyro_public_pilot/results/figures/proposal/spatial_niche_fraction_by_condition_dark.png`
- Evidence matrix: `kthyro_public_pilot/results/figures/proposal/integrated_public_pilot_evidence_matrix_dark.png`
- Drug perturbation pilot: `kthyro_public_pilot/results/figures/proposal/figure6_drug_perturbation_candidate_pilot.png`
- Experimental design: `kthyro_public_pilot/results/figures/proposal/figure7_samsung_experimental_validation_design.png`
- Conclusion slide: `kthyro_public_pilot/results/figures/proposal/figure8_pilot_conclusion_slide.png`

### New validation figures

- Public validation v2 support: `kthyro_public_pilot/results/public_validation_expansion/figures/public_validation_expansion_v2_support_counts.png`
- Exact K-Thyro external axis heatmap: `kthyro_public_pilot/results/public_validation_expansion/figures/exact_kthyro_external_bulk_axis_heatmap.png`
- Exact K-Thyro external boxplots: `kthyro_public_pilot/results/public_validation_expansion/figures/exact_kthyro_external_bulk_boxplots.png`
- Axis separability heatmaps: `kthyro_public_pilot/results/public_validation_expansion/figures/axis_separability_correlation_heatmaps.png`
- Mutation not enough, label fraction: `kthyro_public_pilot/results/public_validation_expansion/figures/tcga_mutation_not_enough_driver_label_fraction.png`
- Mutation not enough, axis variance: `kthyro_public_pilot/results/public_validation_expansion/figures/tcga_mutation_not_enough_axis_variance.png`

---

## 9. Data/table path quick list

### Main pilot tables

- TCGA patient scores: `kthyro_public_pilot/results/tables/tcga_thca_patient_vulnerability_scores.tsv`
- TCGA label summary: `kthyro_public_pilot/results/tables/tcga_label_summary_cleaned.tsv`
- Spatial spot scores: `kthyro_public_pilot/results/tables/spatial_spot_vulnerability_scores.tsv`
- Spatial slide summary: `kthyro_public_pilot/results/tables/spatial_slide_niche_summary.tsv`
- Spatial coherence: `kthyro_public_pilot/results/tables/spatial_coherence_summary_cleaned.tsv`
- Drug cleaned candidates: `kthyro_public_pilot/results/tables/drug_pilot_candidate_rankings_cleaned.tsv`

### Additional validation tables

- Public validation v2: `kthyro_public_pilot/results/public_validation_expansion/tables/external_public_validation_signal_matrix_v2.tsv`
- Public validation v2 summary: `kthyro_public_pilot/results/public_validation_expansion/tables/external_public_validation_layer_summary_v2.tsv`
- Exact external scores: `kthyro_public_pilot/results/public_validation_expansion/tables/exact_kthyro_external_bulk_scores.tsv`
- Exact external contrasts: `kthyro_public_pilot/results/public_validation_expansion/tables/exact_kthyro_external_bulk_contrasts.tsv`
- Axis separability summary: `kthyro_public_pilot/results/public_validation_expansion/tables/axis_separability_dataset_summary.tsv`
- Axis pairwise correlations: `kthyro_public_pilot/results/public_validation_expansion/tables/axis_separability_pairwise_correlations.tsv`
- TCGA driver-label counts: `kthyro_public_pilot/results/public_validation_expansion/tables/tcga_driver_by_vulnerability_label_counts.tsv`
- TCGA driver-label fractions: `kthyro_public_pilot/results/public_validation_expansion/tables/tcga_driver_by_vulnerability_label_fractions.tsv`
- TCGA driver diversity: `kthyro_public_pilot/results/public_validation_expansion/tables/tcga_driver_within_group_vulnerability_diversity.tsv`
- TCGA driver variance explained: `kthyro_public_pilot/results/public_validation_expansion/tables/tcga_driver_axis_variance_explained.tsv`

---

## 10. 제안서 본문용 preliminary evidence paragraph

본 공개데이터 pilot은 K-Thyro Surgical-Spatial Theranostic Perturbation Atlas의 핵심 가설을 다층적으로 지지한다. TCGA-THCA 505명 bulk expression에서 RAI 분화상태, HLA/APM 면역가시성, CD8 exclusion, myeloid/CAF barrier, drug-delivery failure proxy, aggressive dedifferentiation 축을 계산한 결과, 환자군은 RAI-readable differentiated, HLA-visible inflamed, HLA-low immune-invisible, CD8-excluded myeloid/CAF-high, drug-delivery barrier-high, RAI-low dedifferentiated 등 서로 다른 치료취약성 상태로 분리되었다. 특히 BRAF mutation group 내부에서도 7개 vulnerability label이 관찰되어, driver mutation alone으로 치료취약성 상태를 충분히 설명할 수 없음을 보였다. GSE250521 spatial transcriptomics 16개 slide, 57,144개 spot에서는 치료취약성 niche label이 모든 slide에서 비무작위적 공간 응집성을 보여, bulk 평균으로는 포착되지 않는 조직 내 치료취약성 territory가 존재할 가능성을 지지했다. 추가로 외부 GPL570 thyroid cohorts 206샘플을 동일한 K-Thyro gene set으로 재스코어링한 결과 168개 contrast test 중 strong 65개, moderate 12개가 관찰되어 주요 축의 외부 재현성을 보강하였다. 이 결과는 임상 배포용 예측모델이 아니라, FFPE/mIHC, GeoMx/ROI, fresh tissue organoid/slice perturbation, iodide uptake, HLA/APM rescue, fluorescent drug/nanoparticle penetration imaging으로 이어지는 실험 검증 설계의 근거이다.

---

## 11. PPT 제작 순서 추천

1. Title/GO decision slide에는 숫자 3개만 크게: **TCGA 505**, **spatial 16 slides**, **16/16 coherent**
2. TCGA slide는 heatmap보다 먼저 label count 또는 mutation-not-enough figure를 보여도 좋음
3. Spatial slide는 coherence barplot 다음 representative maps 순서가 좋음
4. External validation은 appendix가 아니라 main 중간에 1장 넣으면 reviewer 방어력이 올라감
5. Drug slide는 8번 이후로 배치하고, drug name보다 class를 크게 표시
6. 마지막은 experimental validation loop로 닫기

---

## 12. 최종 one-line proposal message

> **본 과제는 갑상선암을 평균 유전체 질환이 아니라, 수술 perturbation과 공간 치료취약성 niche의 조합으로 치료반응이 결정되는 동적 생태계로 재정의한다.**

