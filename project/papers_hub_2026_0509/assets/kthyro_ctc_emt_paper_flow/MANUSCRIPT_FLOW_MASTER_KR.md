# CTC-EMT-NGS In Silico Paper Flow Master

## 한 줄 정리

지금 논문은 `Serial CTC-EMT transition discovery`가 아니라 **public-data 기반 PTC CTC-EMT-NGS prior map / framework paper**로 써야 한다.

## 논문 제목 후보

**Tissue-State and Genetic Prior Mapping for Post-Thyroidectomy CTC-EMT Monitoring in Papillary Thyroid Cancer**

## 핵심 흐름

1. PTC 수술 후 residual-risk biology는 직접 읽기 어렵다.
2. Yu 2024는 serial CTC-EMT 측정 가능성을 보여준 published anchor다.
3. 하지만 public dataset에는 serial CTC-EMT + matched NGS + cfDNA + tissue context + outcome이 없다.
4. 그러므로 공개 multi-omics로 tissue-state / marker / NGS prior map을 먼저 세운다.
5. 이 prior map은 Samsung/Yu prospective cohort에서 thyroidectomy perturbation으로 테스트한다.

## 주요 숫자

- TCGA BRAF-mutant subset: 274 tumors, 7 tissue-state labels, largest state 33.2%.
- Spatial GSE250521: 16 slides, 57144 spots, same-niche z>2 in 16/16 slides.
- Yu 2024 anchor: 62 prospective PTC patients, pre-op / 2 weeks / 3 months, 87% CTC detection reported.

## 논문에서 금지할 문장

- 우리 데이터가 thyroidectomy 후 CTC-EMT transition을 증명했다.
- public pilot이 CTC shedding을 증명한다.
- early PTC에서 CTC/cfDNA NGS가 universal하게 성공한다.
- 이 결과가 ready clinical diagnostic이다.
