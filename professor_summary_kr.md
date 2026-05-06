# 교수님 보고용 요약

## 결론

이 advance paper는 **진행 가능**합니다. 단, 주제는 예후 논문이나 치료 표적 논문이 아니라 **BRAF/RAS-negative 갑상선암의 분자 taxonomy 논문**으로 잡아야 합니다.

추천 제목:

**Molecular taxonomy of BRAF/RAS-negative thyroid cancer reveals noncanonical driver classes and radioiodine-related differentiation states**

## 이미 확보된 근거

TCGA-THCA mutation-callable cohort에서 BRAF/RAS-negative 군은 **137/482명, 28.4%**입니다. 분석 가능한 크기입니다.

이 군 안에서 DM2/FVPTC-like 쪽에 **DICER1/EIF1AX/PPM1D** 변이가 유의하게 enrichment 됩니다.

- DM1: 1/81 = 1.2%
- DM2: 6/55 = 10.9%
- Fisher p = 0.0175

TERT promoter는 별도 회수 분석에서 TCGA-THCA **36명**이 확인되었습니다. 그중 BRAF+TERT 25명, RAS+TERT 6명, triple-negative/TERT+ 5명입니다. 다만 triple-negative TERT-only는 n=5라서 예후 결론으로 쓰면 안 되고 descriptive class로만 쓰는 것이 안전합니다.

Fusion overlay는 매우 작습니다.

- DM1: ALK fusion 1, NTRK fusion 1, RET fusion 1
- DM2: fusion 0

따라서 fusion은 "overlay"로만 제시하고, fusion-driven subtype이라고 주장하면 안 됩니다.

## Differentiation / RAI gene 축

DM2는 DM1보다 thyroid differentiation 및 RAI 관련 gene score가 높습니다.

- RAI score: DM1 7.73 vs DM2 9.19, p=4.57e-16
- TDS16: DM1 6.93 vs DM2 8.11, p=3.96e-16

Single-cell에서도 기존 결과가 좋습니다.

- GSE241184: 8-gene score와 FVPTC signature 상관 r=0.905
- GSE184362: pooled r=0.889, median patient r=0.865
- multisite GSE184362: pooled r=0.914

따라서 "RAI resistance"가 아니라 **RAI-related differentiation state**라고 표현하는 것이 안전합니다.

## 반드시 피해야 할 주장

예후 claim은 현재 데이터로 불가능합니다.

TCGA BRAF/RAS-negative subset에서 event가 너무 적습니다.

- PFI: 136명 중 9 events, HR 1.33, p=0.675
- OS: 136명 중 6 events, HR 0.83, p=0.824
- DFI/DSS도 events 3개 수준

따라서 survival/recurrence는 secondary exploratory 또는 limitation으로만 넣어야 합니다.

치료 표적 discovery도 주장하면 안 됩니다. 현재 drug-response나 functional assay가 없습니다.

## 외부 validation 후보

가장 중요한 RAI validation은 **GSE151179/GSE151180**입니다.

- GSE151179: mRNA array, 52 samples
- 39 PTC 및 13 matched non-neoplastic thyroid
- RAI-avid / RAI-refractory label 있음
- common PTC mutation/fusion도 PTC-MA로 characterization 되어 있음
- GSE151180은 miRNA subseries, 47/52 samples

Single-cell validation은 다음 두 개가 좋습니다.

- GSE184362: 158,577 cells, 11 PTC patients, 23 specimens
- GSE232237: Korean scRNA-seq, normal/PTC/ATC case 포함

Korean advanced thyroid dataset은 **EGAD00001004845**가 가장 좋지만 controlled access입니다.

- DNA sequencing: 113 advanced thyroid cancers
- RNA sequencing: 25 advanced thyroid cancers
- EGA DAC 승인 전에는 사용 불가

KCG/K-BDS는 KAP210106에 thyroid cancer/thymoma RNA-seq 444 samples가 보이지만, thyroid-only 분리와 임상 endpoint는 로그인/접근 확인이 필요합니다.

CODA는 portal은 확인되었지만 thyroid-specific accession은 이번 조사에서 확인되지 않았습니다. 현재는 blocked로 두는 것이 맞습니다.

## 추천 분석 흐름

1. TCGA에서 BRAF V600E 및 RAS hotspot 음성 샘플을 정의합니다.
2. 그 안에서 DICER1/EIF1AX/PPM1D, TERT-only, fusion/alternative MAPK, true driver-negative로 hierarchy를 만듭니다.
3. SLC5A5, TPO, TSHR, TG, PAX8, NKX2-1 기반 RAI gene score와 TDS16, 기존 8-gene score를 비교합니다.
4. GSE151179에서 RAI-avid vs refractory label로 signature를 검증합니다.
5. GSE184362/GSE232237에서 thyrocyte/malignant cell gradient를 검증합니다.
6. Survival/recurrence는 exploratory limitation으로만 넣습니다.

## 최종 판단

**Go입니다.** 다만 논문 정체성은 "새 치료 표적"이나 "RAI resistance proof"가 아니라, **BRAF/RAS-negative thyroid cancer 안의 noncanonical driver taxonomy와 differentiation/RAI-gene expression state**로 고정해야 합니다.

