# Samsung Preliminary Results Section (KR)

## 1. Public-data pilot purpose

본 pilot의 목적은 제안서 작성 전, 갑상선암 치료취약성이 평균 driver mutation이나 bulk 평균 발현만으로 설명되는지, 아니면 RAI 분화상태, HLA/APM 면역가시성, CD8 exclusion, myeloid/CAF barrier, drug-delivery failure proxy가 분리된 치료축으로 나타나는지를 공개데이터만으로 검증하는 것이었다.

## 2. Data used

- TCGA-THCA primary tumor bulk RNA-seq: 505명.
- GSE250521 spatial transcriptomics: 16개 slide, 57144개 spot.
- scRNA reference: local GSE193581 h5ad 기반 cell-type module sanity check.
- DepMap/PRISM-derived local public resource: cleaned perturbation-class candidate table.

## 3. Key result 1: TCGA therapeutic vulnerability subtypes

TCGA 505명 분석에서 vulnerability label은 {'RAI-readable differentiated': 127, 'Mixed/Other': 125, 'Drug-delivery barrier-high': 105, 'HLA-visible inflamed': 44, 'CD8-excluded myeloid/CAF-high': 42, 'HLA-low immune-invisible': 32, 'RAI-low dedifferentiated': 30}로 분리되었다. 이는 갑상선암이 단일 driver mutation 질환이 아니라 RAI, immune visibility, stromal exclusion, delivery proxy, aggressive dedifferentiation이 조합된 치료취약성 상태로 나뉠 가능성을 보여준다.

## 4. Key result 2: spatial therapeutic niches are coherent

GSE250521 spatial pilot에서 16/16개 slide가 KNN same-niche coherence z-score > 2를 보였다. 이는 spot-level score가 무작위 noise가 아니라 조직 내 공간적으로 응집된 치료취약성 territory를 형성함을 지지한다. 단, spot은 독립 환자 수가 아니므로 모든 stage/condition 비교는 slide 수준 요약으로 해석해야 한다.

## 5. Key result 3: perturbation hypotheses are testable

Drug table은 약 이름 중심이 아니라 perturbation class 중심으로 정리했다. 현재 class 분포는 {'Proliferation stress vulnerability': 47, 'HLA/APM restoration / immune visibility modulation': 20, 'RAI redifferentiation / MAPK-axis modulation': 5, 'Myeloid/CAF barrier modulation': 1, 'Drug-delivery or nanoparticle distribution validation strategy': 1}이다. 제안서에서는 final drug claim이 아니라 RAI redifferentiation, HLA/APM restoration, proliferative stress, myeloid/CAF barrier modulation, delivery imaging strategy로 제시한다.

## 6. Claim boundary

Public pilot은 hypothesis generation 근거이며 clinical deployment 근거가 아니다. Spatial transcriptomics는 peptide presentation을 직접 증명하지 않고, drug-delivery failure는 RNA proxy이며, RAI-restorable score는 iodide uptake assay 전까지 가설이다.

## 7. Why this justifies 30억/3년

30억/3년 과제는 이미 공개데이터에서 관찰된 치료취약성 축과 공간 niche를 실제 수술 검체에서 검증하는 translational closure 단계이다. 비용의 핵심은 신규 발견 탐색이 아니라 FFPE/mIHC, GeoMx ROI, fresh tissue slice/organoid perturbation, iodide uptake, HLA/APM rescue, fluorescent delivery imaging을 하나의 수술 전후 theranostic platform으로 닫는 데 있다.

## 8. Why this can expand to 50억/5년

50억/5년 확장에서는 WES/RNA/HLA typing, serial liquid biopsy, pathology foundation model, multi-region spatial profiling, organoid/slice perturbation atlas를 통합해 환자별 치료취약성 niche의 동역학과 intervention response를 추적할 수 있다. Nature급 proof는 public-data observation이 아니라 paired surgical-spatial-functional perturbation atlas에서 나온다.
