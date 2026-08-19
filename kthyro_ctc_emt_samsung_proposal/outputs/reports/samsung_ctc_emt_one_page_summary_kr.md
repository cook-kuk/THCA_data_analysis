# 삼성 Science 과제 1페이지 요약

## 제목

**갑상선절제술 전후 유두갑상선암 CTC-EMT 상태전이와 유전지도의 연속 해독**  
**Serial Genetic Mapping of CTC-EMT State Transitions After Thyroidectomy in Papillary Thyroid Cancer**

## 핵심 질문

갑상선절제술은 유두갑상선암에서 원발 종양을 제거하는 가장 강한 생체 내 perturbation이다. 수술 후 CTC가 단순히 감소하는지를 넘어서, **상피형(E), hybrid E/M, mesenchymal(M) CTC 상태가 어떻게 전이되는지**, 그리고 지속되는 postoperative CTC/cfDNA 신호가 절제 종양의 유전 clone과 연결되는지를 묻는다.

## 중심 가설

갑상선절제술 후 비특이적 tumor shedding은 감소해야 한다. 그럼에도 지속되는 hybrid/mesenchymal CTC와 tumor-matched genetic signal은 잔존위험 생물학을 반영할 수 있다.

## 왜 Science인가

이 과제는 ICT 플랫폼이 아니라, 수술이라는 인간 perturbation을 이용해 암세포의 circulating EMT state transition을 해석하는 생명과학 과제다. AI/ICT는 clone-state modeling, longitudinal transition modeling, tissue/CTC/cfDNA 통합해석, report generation을 위한 분석 인프라로만 사용한다.

## Preliminary Evidence

Yu 2024는 62명 PTC 전향 코호트에서 수술 전, 수술 후 2주, 3개월 blood CTC 분석을 수행했고, CTC가 87%에서 검출되며 EMT/mesenchymal phenotype이 우세하고 수술 후 CTC count가 감소함을 보였다. 그러나 이 연구는 CTC NGS clone map을 만들지 않았다.

Local public pilot은 TCGA-THCA 505명에서 driver mutation만으로 tissue state가 설명되지 않음을 보였다. BRAF 274명도 7개 vulnerability label로 갈라졌고 최대 label은 33.2%뿐이었다. GSE250521 spatial 16 slides, 57,144 spots에서는 16/16 slides가 same-niche coherence z > 2를 보였다. 이 결과는 CTC 증거가 아니라 marker selection과 tissue-state logic의 보조 근거다.

## 3년 목표

1. PTC thyroidectomy 전후 CTC E/EM/M state transition을 정의한다.
2. 절제 종양 matched tumor-normal NGS로 CTC/cfDNA 신호의 유전적 anchor를 만든다.
3. persistent EM/M CTC와 genetic concordance가 LVI, LNM, ETE, BRAF/TERT/fusion, Tg/US follow-up, recurrence-risk feature와 연결되는지 검증한다.
4. 임상진단이 아니라 future clinical validation을 정당화하는 research-grade residual-risk state framework를 만든다.

## 금지할 주장

본 과제는 CTC-EMT가 이미 재발을 예측한다고 주장하지 않는다. Public pilot이 CTC biology를 증명한다고 말하지 않는다. cfDNA/CTC NGS가 모든 early PTC에서 작동한다고 말하지 않는다. Inocras/vendor 기술을 novelty로 포장하지 않는다.
