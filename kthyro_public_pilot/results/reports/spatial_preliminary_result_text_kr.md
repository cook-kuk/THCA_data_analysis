# Spatial Preliminary Result Text (KR)

## 계산 내용

GSE250521 spatial transcriptomics 16개 slide, 57144개 spot에서 TCGA와 동일한 치료취약성 module score와 spot-level niche label을 계산했다. Spot은 독립 생물학적 반복으로 취급하지 않고 slide 내부에 nested된 관측치로 다루었다.

## 핵심 결과

KNN same-niche permutation 결과, 16/16개 slide에서 niche coherence z-score > 2가 관찰되었다. 조건 분포는 {'normal': 4, 'PTC': 4, 'locally_advanced_PTC': 4, 'ATC': 4}이다.


상위 coherence slide:


| sample_id         | condition            |   n_spots |   same_niche_z |   same_niche_empirical_p |
|:------------------|:---------------------|----------:|---------------:|-------------------------:|
| GSM7980870_LPTC-3 | locally_advanced_PTC |      4355 |        43.1644 |               0.00199601 |
| GSM7980864_PTC-1  | PTC                  |      2243 |        38.1139 |               0.00199601 |
| GSM7980873_ATC-2  | ATC                  |      2656 |        37.2604 |               0.00199601 |
| GSM7980869_LPTC-2 | locally_advanced_PTC |      4129 |        33.2437 |               0.00199601 |
| GSM7980871_LPTC-4 | locally_advanced_PTC |      4688 |        28.8097 |               0.00199601 |
| GSM7980866_PTC-3  | PTC                  |      4475 |        25.8125 |               0.00199601 |
| GSM7980872_ATC-1  | ATC                  |      4672 |        23.9135 |               0.00199601 |
| GSM7980875_ATC-4  | ATC                  |      2566 |        22.3509 |               0.00199601 |

## 조심스러운 해석

이 결과는 치료취약성 expression state가 조직 내에서 공간적으로 응집되어 있음을 지지하지만, 임상 반응 예측이나 실제 drug delivery를 직접 증명하지 않는다.

## 삼성 제안서 문장

공간전사체 pilot에서 치료취약성 niche는 무작위 spot noise가 아니라 조직 내에서 공간적으로 응집된 구조를 보였으며, 이는 FFPE/mIHC, GeoMx ROI, fresh tissue 기능실험으로 이어지는 수술-공간 theranostic platform의 필요성을 뒷받침한다.

## Figure caption

Public GSE250521 spatial transcriptomics reveals non-random spatial organization of therapeutic vulnerability niches across PTC, locally advanced PTC, and ATC slides. Spots are visualized for spatial structure, while inference is summarized at slide level.
