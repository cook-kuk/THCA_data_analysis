# Paper 9 교수님 업데이트: 대사 취약성 우선순위 지도

## 핵심 결론
Paper 9는 GLS 단일 타깃의 synthetic lethality 논문이 아니라, lineage-silenced thyroid cancer에서 대사 취약성 후보를 우선순위화하는 지도로 재정의하는 것이 가장 안전합니다.

## 현재 GLS 상태
- GLS verdict: KEEP_AS_CANDIDATE (MEDIUM_MINUS).
- pan-cancer dependency 방향은 가설과 맞지만 효과 크기는 작고, thyroid-only 분석은 모델 수가 부족합니다.
- GDSC BPTES drug concordance는 현재 지지적이지 않습니다.

## 확장 방향
- GLS 단독이 아니라 glutamine transport/metabolism module, MYC glutamine-addiction module, glycolysis/OXPHOS/NAD salvage comparator를 함께 평가합니다.
- 이 프레임은 negative/weak drug evidence를 숨기지 않고, wet-lab 검증 후보를 선별하는 구조입니다.

## 사용할 수 있는 표현
- lineage-silenced thyroid cancer may exhibit glutamine-axis vulnerability candidates requiring experimental validation.
- metabolic vulnerability prioritization map.

## 피해야 할 표현
- validated synthetic lethality.
- clinical treatment recommendation.
- patient selection biomarker.
