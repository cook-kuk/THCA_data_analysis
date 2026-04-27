# U2B — Han 2023 Korean parallel narrative

## 한국어 요약 (Han SC, Park YJ et al. ENM 2023, PMID 37461149)

**원논문 (Han 2023):** TCGA-PTC 503명을 분석하여 BRAF-like (BL)
및 RAS-like (RL) 분자 아형의 진행 메커니즘이 다름을 확인. BL-PTC의 공격적 케이스는
세포외기질(ECM) 유전자 상향조절과 CAF 풍부도 증가, RL-PTC의 공격적 케이스는
면역반응 유전자 하향조절을 보임.

**본 연구 (v17 DM1/DM2):**

- DM1 (n=140) vs DM2 (n=360) marker 유전자 비교.
- DM1 상향 top50 ∩ ECM 패널: **0개** (없음)
- DM2 상향 top50 ∩ ECM 패널: **1개** (FN1)
- DM1 상향 top50 ∩ 면역 패널: **0개** (없음)
- DM2 상향 top50 ∩ 면역 패널: **0개** (없음)

**v17 Hot vs Cold 면역 Cohen's d:** +1.683 (대형 효과). Han 2023 초록은 효과 크기를
보고하지 않음.

**정렬 가능성 (DM1 ↔ BL-aggressive):** low
(점수 0/3 — ECM_DM1, IMMUNE_DM1 풍부, ECM_DM1 > ECM_DM2 차이 평가).

**주의사항:**
1. Han 2023 초록은 구체적 유전자명을 나열하지 않아 표준 ECM/면역 패널로 대리 비교.
2. 두 연구 모두 TCGA를 사용 — 외부 검증 아닌 평행 코호트 분석.
3. v17 Hot/Cold 클러스터 vs Han의 BL/RL aggressive 분류는 정의가 다름. 정성적 정렬만 가능.

## English summary

Han 2023 analyzed TCGA-PTC (N=503) and reported BL-PTC aggressive cases show
ECM upregulation + CAF enrichment, while RL-PTC aggressive cases show immune
downregulation. Our DM1/DM2 marker genes were screened against canonical ECM
and immune panels:

| Panel | DM1-up overlap | DM2-up overlap |
|---|---|---|
| ECM (20 genes) | 0 | 1 |
| Immune (20 genes) | 0 | 0 |

v17 Hot/Cold immune Cohen's d = +1.683 (no comparable Han effect size reported).

Alignment likelihood DM1 ↔ BL-aggressive: **low** (score 0/3).
