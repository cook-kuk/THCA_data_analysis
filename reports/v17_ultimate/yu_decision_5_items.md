# 유 교수님 confirm 요청 — 5개 결정 항목 (v17 ULTIMATE)

**날짜**: 2026-04-27 KST
**작성**: 국승호 (Seungho Cook)
**문서 역할**: ULTIMATE sprint 결과 정리 + 유 교수님 직접 confirm 받을 5개 항목

---

## 1. Cluster orientation 정정 (v5 → v6 ULTIMATE)

**v5 preprint 표기**: "DM1 = MAPK-active dedifferentiated, DM2 = canonical thyroid-differentiation"

**ULTIMATE 분석 결과**: 4개 독립 증거 모두 v5 표기가 **반대**임을 확인:
- U2A WHO 2022 매핑: DM2 = 85.6% BRAF-like cPTC, DM1 = 77.5% RAS-like FVPTC
- U2C Pu 2021 direction concordance: 88.9% (8/9 signature genes)
- U1D methylation: 모든 aggressive histology 100%가 met_DM2
- U1C GSE213647: cPTC > ATC dm-score (DM1 = 분화된 것이 맞음)

**v6 ULTIMATE 정정 표기**:
- **DM1 = 분화 (well-differentiated, RAS-like, FVPTC/IEFVPTC, panel-high)**
- **DM2 = 탈분화 (dedifferentiated, BRAF-like, cPTC/TCV/DHGTC, panel-low)**

**유 교수님 confirm 필요**: 이 정정 표기로 manuscript v6 ULTIMATE 진행해도 되는지?

---

## 2. De-circularization 정직 보고 framing

**문제**: v5 preprint의 8-gene panel CV AUC = 0.954는 cluster 정의에 8 gene 자체가 포함되어 있어 자기상관 (circular). 진짜 prediction 아님.

**ULTIMATE 수정**: 4개 leak-free variant 재생성, 진짜 AUC 측정:
- Variant A (TIERA67-8gene, 47 genes): **AUC 0.962** [0.940-0.979]
- Variant B (MAPK+immune+EMT, 30 genes, zero overlap): **AUC 0.925** [0.895-0.952]
- Variant D (BRS71-8gene, 26 genes): **AUC 0.957** [0.936-0.976]
- Variant C (immune-only 5 genes): AUC 0.664 (붕괴 — biology specificity 입증)

**ΔAUC vs BRAF V600E baseline**: 모든 variant +0.111 ~ +0.130

**Cover letter v5 framing 제안**:
> "Following pre-submission self-audit, we identified a circular validation issue in our v5 manuscript and corrected it through four independent leak-free cluster reconstructions. The reported AUC values now reflect honest predictive performance (range 0.92–0.96 across three leak-free constructions); the circular 0.954 figure is no longer claimed."

**유 교수님 confirm 필요**: 이 honest reframe으로 진행?

---

## 3. GSE213647 (n=632) 결과 — manuscript 어떻게 표기?

**결과**: 가장 큰 외부 cohort (n=632; 370 tumour, 262 normal)에서 8-gene panel 적용:
- 8/8 gene 모두 매칭 성공
- cPTC dm_score 0.193 > PDTC 0.103 > ATC 0.080 (DM1=분화 convention과 일치)
- 단, "aggressive vs PTC" AUC = 0.672 (TCGA in-sample 0.962 대비 낮음)
- Cross-platform / batch heterogeneity가 일부 영향

**3가지 옵션**:

a) **명시적 limitation으로 보고** (추천): "GSE213647에서 효과 크기가 TCGA in-sample 대비 감소 — cross-platform generalisability flag, 미해결 future work"
b) **별도 supplementary**: GSE213647 결과만 supp figure로 빼고 본문 forest plot에서 제외
c) **본문 strong, limitation 짧게**: 결과를 본문에 포함하되 limitation 한 문장만

**유 교수님 confirm 필요**: 옵션 a / b / c 중?

---

## 4. Park YJ outreach 보내는 timing

**제안된 outreach 자료** (`outreach/email_park_yj_KR.md`):
- 박영주 SNU 교수 (Han SC ENM 2023 senior author)
- Korean parallel work positioning + 협력 framework 3개

**Timing 옵션**:

a) **즉시 (오늘-내일)**: Manuscript v6 ULTIMATE draft + dashboard 첨부, 1-2주 내 submit 직전 협력 가능성 빠르게 sound out
b) **Submit 직후**: npj submit 후 첨부 강도 높여서 응답률 높임
c) **Yu confirm 받은 후**: Manuscript final 검토 + 첨부 정합성 100% 후 발송

**추천**: a (즉시) — Park YJ 응답 24-48h 시 cover letter / Discussion에 "Korean parallel work outreach in progress"가 더 강함

**유 교수님 confirm 필요**: 발송 timing?

---

## 5. Submit venue 최종 결정

**Scenario 분류 (v6 ULTIMATE 결과 기반)**:

| 기준 | ULTIMATE 결과 | Threshold | 통과? |
|---|---|---|---|
| Best variant ARI | 0.642 | > 0.4 | ✓ |
| Honest 8-gene CV AUC | 0.962 | > 0.85 | ✓ |
| ΔAUC vs BRAF (95% CI) | +0.113 [+0.071, +0.155] | exclude 0 | ✓ |
| External GSE76039 AUC | 0.935 | > 0.80 | ✓ |
| WHO 2022 alignment OR | 20.4 (p=2.5e-33) | > 5 | ✓ |
| Methylation cross-modality | 100% (aggressive→met_DM2) | qualitative strong | ✓ |
| GSE213647 generalisability | 0.672 (감소) | preserve | ✗ |

**6/7 통과** → **Scenario A (npj Precision Oncology) 추천**

대안:
- BMC Cancer / Frontiers Oncol: 1점 fail이면 더 안전한 venue
- Bioinformatics: methodology paper로 pivot (DIAL framework)

**추천**: **npj Precision Oncology submit** + GSE213647 limitation 정직 명시

**Submit timeline**:
- 유 교수님 confirm: 1-2일
- Park YJ outreach 응답 대기: 24-48h
- 분당서울대 outreach: 24h 이내 발송
- Manuscript v6 ULTIMATE final read-through: 2-3시간
- npj submit: **1주 이내**

**유 교수님 confirm 필요**: npj submit 진행해도 되는지?

---

## 결정 요약 표

| # | 항목 | 추천 | 유 교수님 결정 |
|---|---|---|---|
| 1 | Cluster orientation 정정 | DM1=분화/DM2=탈분화 | ☐ confirm |
| 2 | Honest reframe (0.954 → 0.92~0.96 range) | yes | ☐ confirm |
| 3 | GSE213647 limitation framing | 옵션 a (본문 + limitation) | ☐ a / ☐ b / ☐ c |
| 4 | Park YJ outreach timing | 즉시 발송 | ☐ confirm / ☐ 다른 옵션 |
| 5 | Submit venue | npj Precision Oncology | ☐ confirm / ☐ 다른 venue |

---

## 본인 (국승호) 다음 액션

유 교수님 confirm 받으면:

1. **분당서울대 outreach v2 발송** (24h 이내)
2. **Park YJ outreach 발송** (timing 결정 후)
3. **Manuscript v6 ULTIMATE final polish** + figure 정합성 검토
4. **Cover letter v5 작성** (honest reframe + Korean parallel work + Yu confirm 반영)
5. **npj submit click**

---

## 첨부 / 참고

- 본 manuscript v6 ULTIMATE: `project/submission/npj/manuscript_v6_ULTIMATE.md`
- 인터랙티브 종합 리포트: `interactive_report.html` (현재 280KB, 28 figs + 17 noteboooks)
- ULTIMATE 결과 폴더: `project/results/v17_ultimate/` (16 JSON + 8 TSV + 5 PNG)
- Outreach 메일 3종: `outreach/email_park_yj_KR.md`, `email_bundang_KR_v2.md`, `email_wetlab_5labs_KR.md`
- Reviewer defense v3: `project/submission/npj/v17p35_REVIEWER_DEFENSE_v3.md` (25 attacks)
