# 10-Slide Korean Copy

## 1. 공개데이터 pilot 결과: K-Thyro는 GO

메시지:  
TCGA 505명과 GSE250521 spatial 16 slide에서 치료취약성 축이 분리되고, 공간 niche가 비무작위적으로 응집된다.

Bullets:

- TCGA-THCA 505명에서 7개 vulnerability label 도출
- GSE250521 16 slide, 57,144 spots에서 치료취약성 niche mapping
- 16/16 slide에서 KNN same-niche coherence z-score > 2
- Public validation matrix v2: 301 tests, strong/moderate 174 tests
- 최종 판단: GO, 단 public-data pilot은 hypothesis generation 단계

---

## 2. 공개데이터 기반 K-Thyro 분석 흐름

메시지:  
Bulk, spatial, scRNA/GeoMx, external cohort, drug resource를 하나의 치료취약성 scoring framework로 통합했다.

Bullets:

- TCGA-THCA: patient-level vulnerability scoring
- GSE250521 spatial: spot-level niche mapping and coherence test
- External GPL570 cohorts: exact K-Thyro gene set reproducibility
- scRNA/GeoMx: cell-type and ROI-level sanity check
- DepMap/PRISM: perturbation class hypothesis nomination

---

## 3. TCGA 505명에서 치료취약성 축이 환자별로 분리된다

메시지:  
평균 bulk mutation만으로는 RAI, immune visibility, barrier, proliferation 상태를 모두 설명할 수 없다.

Bullets:

- RAI differentiation, HLA/APM, CD8/cytotoxic, myeloid/CAF, drug-delivery proxy, aggressive axis 계산
- 환자별 vulnerability label이 RAI-readable, HLA-visible, HLA-low, CD8-excluded, barrier-high 등으로 분화
- Bulk 평균에서도 치료취약성 축의 separability가 관찰됨
- 공간 검증 없이는 niche biology를 주장하지 않음

---

## 4. Driver mutation만으로 치료취약성 상태를 결정할 수 없다

메시지:  
BRAF group 내부에서도 7개 vulnerability label로 갈라져 mutation-only 전략의 한계를 보인다.

Bullets:

- BRAF group 274명: 7개 vulnerability label로 분산
- BRAF 최대 label은 Drug-delivery barrier-high 91명, 33.2%
- BRAF만 알면 66.8%는 vulnerability state가 ambiguous
- HLA-I/APM axis driver 설명력 eta² 0.156, immune visibility eta² 0.125
- Mutation + spatial/immune/barrier state를 함께 읽는 companion strategy 필요

---

## 5. 치료취약성 niche는 spatial spot noise가 아니라 조직 내 구조다

메시지:  
GSE250521 16개 slide 모두에서 같은 niche가 공간적으로 이웃하는 non-random coherence를 보였다.

Bullets:

- GSE250521: normal 4, PTC 4, locally advanced PTC 4, ATC 4
- 총 57,144 spots에 동일한 K-Thyro vulnerability score 적용
- KNN same-niche permutation: 16/16 slide z-score > 2
- Coherence z-score min/median/max: 6.09 / 22.00 / 43.16
- Spots는 biological replicate가 아니며 slide-level summary로 해석

---

## 6. RAI, HLA/APM, CAF/myeloid barrier niche가 공간적으로 분리된다

메시지:  
같은 조직 안에서도 RAI-low, HLA-visible, CD8-excluded, barrier-high territory가 서로 다른 위치에 나타난다.

Bullets:

- Spot-level RAI, HLA/APM, cytotoxicity, CAF/myeloid barrier score 계산
- Quantile threshold 기반 exploratory niche label 부여
- Representative slide에서 치료취약성 niche map 시각화
- Spatial map은 FFPE/mIHC, GeoMx/ROI, functional assay의 ROI 선정 근거

---

## 7. 5개 치료취약성 축은 다층 공개데이터에서 상보적으로 지지된다

메시지:  
TCGA, spatial, external bulk, GeoMx/scRNA, perturbation resource가 서로 다른 claim boundary 안에서 같은 platform hypothesis를 지지한다.

Bullets:

- Public validation matrix v2: 301 tests
- strong 144, moderate 30
- External GPL570 exact K-Thyro rescoring: 206 samples, 168 tests
- Axis separability: 축들이 하나의 score로 붕괴되지 않음
- Integrated evidence는 clinical deployment가 아니라 validation design의 근거

---

## 8. Drug 후보가 아니라 perturbation class를 검증한다

메시지:  
Public drug resources는 final therapy가 아니라 RAI/HLA/barrier/drug-delivery niche를 재프로그래밍할 실험 후보군을 제안한다.

Bullets:

- MAPK-axis modulation: RAI redifferentiation hypothesis
- Epigenetic/IFN/APM modulation: HLA/APM visibility restoration hypothesis
- Cell-cycle stress vulnerability: aggressive/proliferative niche hypothesis
- Myeloid/CAF barrier modulation: CD8 exclusion 해소 후보축
- Drug-delivery strategy: fluorescent drug/liposome/nanoparticle penetration imaging으로 검증

---

## 9. AI-spatial inference를 FFPE/GeoMx/기능실험으로 닫는다

메시지:  
Public pilot에서 나온 niche를 실제 환자 조직과 perturbation assay로 검증하는 폐루프 플랫폼을 구축한다.

Bullets:

- FFPE cohort: spatial marker validation and ROI selection
- mIHC/IF core panel: PAX8/TG/TPO/NIS, HLA-I/B2M/TAP1, CD8/GZMB, CD68/CD163, ACTA2/FAP/COL1A1, PD-L1
- GeoMx/ROI: HLA-low, CD8-excluded, RAI-low, CAF/myeloid-rich territory 검증
- Fresh tissue organoid/slice: iodide uptake, HLA/APM rescue, PBMC/TIL co-culture
- Fluorescent drug/liposome/nanoparticle imaging: delivery failure proxy 검증

---

## 10. 30억/3년: 공간 치료취약성 생태계를 실험으로 증명한다

메시지:  
공개데이터 pilot은 GO를 만들었고, 본 과제는 이를 실제 수술 전후 perturbation atlas와 기능검증 플랫폼으로 확장한다.

Bullets:

- Year 1: FFPE/mIHC + GeoMx/ROI로 niche atlas 구축
- Year 2: fresh tissue slice/organoid perturbation으로 RAI/HLA/barrier rescue 검증
- Year 3: AI-spatial companion model + theranostic validation package 완성
- 30억/3년은 clinical deployment가 아니라 high-risk/high-impact PoC 구축 예산
- 50억/5년 확장 시 prospective cohort와 interventional validation로 확장 가능

