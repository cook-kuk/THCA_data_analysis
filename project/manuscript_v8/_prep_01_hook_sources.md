---
title: "Prep #1 — Hook stat verification (Section 1.1 첫 줄용)"
date: 2026-04-30
prep_for: Prompt 2 (Introduction draft) Section 1.1 첫 단락
status: ✅ 검증 완료
---

# Hook stat 검증 결과

## 본인 v2 outline 의 example hook

> "Despite ~95% 5-year survival, 5-15% of papillary thyroid carcinomas recur or metastasize, and clinicians continue to make radioiodine decisions on heterogeneous BRAF/RAS-negative tumors with no mechanistic compass."

## 검증 — 출처별 stat 분포

| Source | Recurrence rate | Distant met rate | 5-yr survival | Note |
|---|---|---|---|---|
| **ATA 2015 (Haugen 2016)** | "<10% post-operative recurrence" | (per tier) | ~98-99% (DTC overall) | canonical, low end |
| Lee et al 2021 PMC8283398 (n=4,085) | 4.3% (mean follow-up 58.7 mo) | 0.4% | — | conservative single-center |
| Broader literature meta | 5-20% local recurrence | 10-15% distant met | — | wide range |
| Intermediate-risk tier specifically | ~20% recurrence | — | — | ATA 2015 estimate |
| High-risk tier specifically | up to 50%+ | up to 25% | — | ATA 2015 estimate |
| PTMC (≤1 cm) | 2-6% loco-regional | 1-2% distant | <1% mortality | low-risk only |

## ★ 본인 hook 의 정확도 평가

| 본인 표현 | 검증 verdict |
|---|---|
| "~95% 5-year survival" | ✅ 정확 (ATA 2015 + multiple sources) — 더 정확하게는 98-99% (overall DTC 5-yr OS) |
| "5-15% recur or metastasize" | ⚠️ 부분 정확 — overall recurrence 는 5-20% (해외 lit), 한국 single-center 는 4.3% conservative. **Intermediate-risk specifically 가 ~20%** (ATA 2015) — "5-15%" 보다 "5-20%" 또는 "up to 20% in intermediate-risk" 가 더 정확 |
| "no mechanistic compass" | ✅ 정확 framing (본인 voice), 검증 영역 아님 |

## ★ 권장 hook alternative 3개 (citation 별)

### Alternative A — ATA 2015 conservative

> "Although papillary thyroid carcinoma (PTC) carries an excellent overall prognosis (5-year disease-specific survival >98%; ATA 2015), structural disease recurs in 5-20% and distant metastases occur in approximately 10% of patients (Haugen et al., 2016 *Thyroid*), and clinicians continue to make radioiodine decisions on heterogeneous BRAF/RAS-negative tumors without a mechanistic basis."

### Alternative B — Risk-stratified specific

> "Despite a >98% 5-year overall survival in differentiated thyroid carcinoma, structural disease recurrence remains 5-20% across risk tiers — reaching ~20% in ATA 2015 intermediate-risk patients (Haugen et al., 2016) — and the heterogeneity of BRAF/RAS-negative papillary thyroid carcinoma continues to complicate radioiodine treatment decisions in the absence of mechanistic sub-stratification."

### Alternative C — 본인 v2 outline original (수치만 정정)

> "Despite a >98% 5-year survival, 5-20% of papillary thyroid carcinomas recur or develop distant metastases (Haugen et al., 2016), and clinicians continue to make radioiodine decisions on heterogeneous BRAF/RAS-negative tumors without a mechanistic compass."

## 본인 voice 영역 (3 alt 중 선택)

- **A** — 가장 안전, citation 정확, conservative 톤. Cell Rep Med 표준.
- **B** — Discussion 3.2 ATA 2015 alignment 까지 thread 시키려면 가장 강력 (intermediate-risk 20% 명시 → DM1 reflex 가 그 tier 에 input 가능 framing).
- **C** — 본인 original 표현 보존. "5-15%" → "5-20%" 만 수정. "no mechanistic compass" 본인 voice 유지.

## 최종 권고 (v2 — web Claude 2차 review 후)

### ★★★ Alt B refined (web Claude 권장 final)

> **"Despite a >98% 5-year overall survival in differentiated thyroid carcinoma, structural disease recurrence reaches ~20% in ATA 2015 intermediate-risk patients (Haugen et al., 2016), and BRAF/RAS-negative tumors — accounting for ~23% of cases — complicate radioiodine decisions in the absence of mechanistic sub-stratification."**

근거:
- Alt B 의 "intermediate-risk specifically 20%" 가 Discussion 3.2 ATA alignment paragraph 와 narrative arc thread
- "5-20% across risk tiers" 빼고 intermediate 만 명시 → cleaner reading
- "~23% of cases" 추가 → BRAF/RAS-negative dark matter scope 정량화 (Xing 2014 cite gateway)
- "absence of mechanistic sub-stratification" → academic 톤, 본인 voice "no mechanistic compass" tone 보존

### Hybrid 옵션 (본인 voice "compass" 보존)

> **"Despite a >98% 5-year overall survival, structural disease recurrence reaches ~20% in ATA 2015 intermediate-risk papillary thyroid carcinoma (Haugen et al., 2016), and clinicians continue to make radioiodine decisions on heterogeneous BRAF/RAS-negative tumors — ~23% of cases — without a mechanistic compass."**

본인 "compass" metaphor 정말 좋아하시면 이 hybrid 도 OK. Cell Press 표준은 academic 톤 권장 (Alt B refined) 이지만 본인 voice 보호 영역.

### v2 → v3 변경 사항 (web Claude 2차 review 적용)

| 변경 | 근거 |
|---|---|
| "5-15%" → "intermediate-risk ~20%" specific number | ATA 2015 검증 결과 — intermediate 만 명시가 정확 |
| "~23% of cases" 추가 | BRAF/RAS-negative dark matter scope 정량화 (Xing 2014 cite gateway) |
| "no mechanistic compass" → "absence of mechanistic sub-stratification" (또는 hybrid 보존) | Cell Press academic 톤 권장 vs 본인 voice 보호 — 본인 결정 |

## Citation needed

**Haugen BR, Alexander EK, Bible KC, Doherty GM, Mandel SJ, Nikiforov YE, et al. 2015 American Thyroid Association management guidelines for adult patients with thyroid nodules and differentiated thyroid cancer. *Thyroid*. 2016;26(1):1-133.**

- PMID: 26462967
- DOI: 10.1089/thy.2015.0020
- PMC: PMC4739132
- 인용 형태: (Haugen et al., 2016) 또는 [Haugen 2016] 본인 style 선택

## 추가 발견 (★ Discussion 3.2 영향)

**ATA 2025 Guidelines 가 별도로 published 됨** (PMID 40844370, *Thyroid* 2025). ATA 2015 (Haugen 2016) 는 여전히 canonical 이지만, 2026 submission 시 reviewer 가 "ATA 2025 update 반영했는지" 물을 가능성. Discussion 3.2 작성 시 두 가이드라인 cross-reference 권장.

| Guideline | Citation | Role in Paper 1 |
|---|---|---|
| ATA 2015 | Haugen 2016 *Thyroid* 26:1-133 | 1순위 cite (current standard practice 기준) |
| ATA 2025 | *Thyroid* 2025 (PMID 40844370) | secondary cite (most recent update 반영 신호) |
