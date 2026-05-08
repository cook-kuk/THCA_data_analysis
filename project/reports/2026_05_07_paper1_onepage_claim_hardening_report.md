# Paper 1 one-page audit — Claim hardening report

**Date:** 2026-05-07
**Scope:** one-page audit page (`project/manuscript_v8/p1_onepage_audit.html`) + master MD
**URL:** http://40.82.129.113:8012/manuscript_v8/p1_onepage_audit.html
**Mode:** claim hardening + reviewer-proof cleanup — **no new analysis, no new data, no voice-protected prose**.

---

## 1. Executive verdict

**Hardening complete.** Page 의 clinical RAI-refractory framing 은 유지하면서, *validated predictor / treatment recommendation / prospective utility* 등 모든 affirmative overclaim 을 hypothesis-language 로 재작성. Panel combination scan 의 framing 은 "8-gene superiority" 에서 **"axis robustness / panel-size sensitivity (within-cohort sensitivity analysis)"** 로 reframe. RAI_8 anchor 가 main, panel scan 은 supplement / reviewer defense 위치.

**Final overclaim sweep:** all banned phrases ("능가", "outperform", "더 정확", "decitabine candidate", "selpercatinib eligible", "low-cost", "8-gene superiority", "replace RAI") = 0 affirmative hits. Residual occurrences ("validated predictor", "treatment selection", "predicts RAI response", "treatment recommendation") = ALL in explicit negation / forbidden-context only.

---

## 2. Overclaim hits and fixes

| Banned / risky phrase | Before (affirmative) | After (hardened) |
|---|---|---|
| "panel 더 작고 더 정확 / 능가" | "RAI_8 보다 좋은 compact panel 존재" / "6 개로 RAI_8 능가" / "4 개로 RAI_8 능가" | "panel-size sensitivity: 6/4-gene 부분집합도 비슷 within-cohort AUC 범위; cross-cohort validation 필요" |
| Decitabine + I-131 trial | "본 paper 는 그 trial 의 적격군 정의 도구가 될 수 있음 (DM1 high = decitabine + I-131 candidate)" | "임상 mechanism hypothesis (NOT therapeutic validation): mechanistic plausibility 와 방향 일치; treatment recommendation 도구 또는 validated 적격군 정의 도구 아님" |
| Selpercatinib reflex | "selpercatinib reflex testing 81.8% RET 적격성 등" + "~48 selpercatinib-eligible per 1000 PTC" | "TCGA observation: DM1 captures 81.8% of RET-fusion+ — hypothesis-generating only; treatment selection tool 또는 selpercatinib eligibility predictor 가 아님" |
| 임상 행동 시나리오 | "더 일찍 systemic / molecular 옵션 고려" | "가설적 임상 행동 (NOT validated) = future-work 시나리오" |
| Hero conclude | "더 일찍 flag 할 수 있는 후보 triage scaffold" | "+ (NOT a validated predictor / NOT a treatment-selection tool / NOT proven prospective utility) 명시적 부정" |
| Future-validation board (decitabine) | "epigenetic re-induction prospective; retrospective 정리 → 새 trial design" | "(future-work only); design hypothesis 만; therapeutic validation 또는 candidate definition tool 아님" |
| §3 Clinical RAI failure | "8-gene readout 의 RAI failure 예측력 입증" | "RAI 실패 *flag 능력* (예측력 아님; flag 는 trigger; 예측은 deterministic) 검증" |

**Residual occurrences (verified safe):**
- L131: NOT a validated predictor / NOT a treatment-selection tool / NOT proven prospective utility — explicit negation
- L547: "treatment recommendation 도구가 아님" — explicit denial
- L641: forbidden table (claim boundary)
- L721: Q&A answer "Q: 예측합니까? A: hypothesis-generating only"
- L768: future-work only + "treatment selection tool 또는 selpercatinib eligibility predictor 가 아님"
- L792: IS NOT card — "treatment selection tool"

→ all residual mentions are in *what the paper is NOT* contexts. ✅

---

## 3. Panel scan interpretation fix

**Before (problematic framing):**
- "Lineage TF + Effector minimal 6 ... RAI_8 (0.861, 0.220) 을 양 metric 모두에서 능가"
- "panel 더 작고 (6 vs 8) 더 정확"
- "Option B — 6-gene 으로 demote/replace"

**After (hardened framing):**
- "panel-size sensitivity: 6-gene 부분집합도 비슷 범위 (cross-cohort validation 필요)"
- "axis robustness: 4-gene 부분집합도 비슷 신호 (within-cohort sensitivity)"
- "★ 권장 frame: 본 scan 결과는 supplementary 'axis robustness / panel-size sensitivity analysis' 로 위치. main manuscript 는 RAI_8 anchor 유지"

**Key change:** within-cohort empirical sensitivity ≠ cross-cohort optimality. 본 scan 은 axis 가 panel-size 에 robust 함을 보여주지만, 다른 panel 으로의 교체 권고는 아님.

**Reviewer Q "왜 8 genes?" 답변 reframing:**
- 이전: "8 은 Yoo 2016 canonical 이지만, 6-gene minimal 또는 4-gene minimal 도 동등 deployment 정확도 → cost / clinical practicality / panel-size sensitivity 에 따라 *prospectively* 선택 가능"
- 현재: "8 genes = canonical Yoo 2016 RAI biology anchor + 임상의 친숙성 + RT-qPCR deployable size. panel-size sensitivity scan 은 axis 가 6-12 gene 모두에서 robust 함을 보여주며, 이는 main 의 RAI_8 anchor 를 강화하는 것이지 대체하는 것이 아님."

---

## 4. Clinical framing fix

| Topic | Before | After |
|---|---|---|
| Hero | "더 일찍 flag 할 수 있는 후보 triage scaffold" | "가설적 candidate triage scaffold; NOT validated predictor / NOT treatment-selection tool / NOT proven prospective utility" |
| §3 임상 행동 | "더 일찍 systemic / molecular 옵션 고려 (selpercatinib reflex 등)" | "가설적 임상 행동 (NOT validated) = future-work 시나리오; prospective trial 에서 검증" |
| §9 Mechanism implication | "decitabine + I-131 candidate" / "trial 의 적격군 정의 도구" | "mechanism hypothesis (NOT therapeutic validation); mechanistic rationale 강화 또는 후보군 식별 framework 가설; treatment recommendation 도구 아님" |
| Future-validation board | "decitabine + I-131 mechanism trial" / "selpercatinib reflex prospective" | both labeled "(future-work only)"; "design hypothesis"; "current paper 는 therapeutic validation 아님" |

**Allowed safe wording (현재 사용 중):**
- candidate triage scaffold (hypothesis only)
- RAI-lineage failure biology
- earlier molecular review hypothesis
- future prospective validation
- complementary to ATA/RAI workflow
- consistent with / supports
- mechanistic plausibility

**Avoided:**
- clinical deployment now
- skip RAI
- direct hemato-oncology referral
- predicts refractory disease
- treatment decision

---

## 5. Figures / captions fixed

| Figure | Issue | Fix |
|---|---|---|
| Fig 12 (thyroid pathway) | static PNG only | added editable <code>fig12_thyroid_pathway.pptx</code> + ★★★ 표시 |
| Fig 13 (mechanism cascade) | "임상 implication: decitabine + I-131 candidate" | hardened to "mechanism hypothesis (NOT therapeutic validation)"; 후보군 식별 framework 가설 강조 |
| External GPL570 4-cohort scatter | unstarred | ★★★ 표시 + DEEP-DIVE 1 (LOO per-gene contribution) 추가 |
| PFI survival fix | unstarred | ★★★ 표시 + DEEP-DIVE 2 (sub-A/sub-B 재산출) 추가 |
| §10 panel scan tables | "능가" / "더 정확" | "panel-size sensitivity: 비슷 범위 (cross-cohort validation 필요)" |

---

## 6. Remaining risky statements (manual review needed)

1. **Hero h1 title** — "compact RAI-lineage transcriptomic readout for early triage of post-surgery RAI-failure biology" 는 *triage scaffold* framing 유지 (hypothesis only). 만약 reviewer 가 *triage scaffold* 자체도 overclaim 으로 받아들이면, "transcriptomic readout for the RAI-lineage axis in thyroid cancer" 처럼 더 중립적 title 후보 4 (table §15) 사용 권장.
2. **"early triage" wording** — "early" 가 *outcome 개선 implication* 으로 읽힐 수 있으므로 *future-work* 영역에만 한정해야 함 (현재 그렇게 표시되어 있음).
3. **Sub-A/sub-B 재산출 결과 (DEEP-DIVE 2)** — 새 partition (37 / 73) 이 manuscript 기존 caption (72 / 19) 과 직접 비교 불가능. <b>action required</b> 라고 page 에 명시.
4. **Single-cell GSE241184 n=1** — Phase 1 에서 author-independence 분석이지만 n=1 에 대한 framing 강화 필요.
5. **Decitabine + I-131 mechanism**: future-work only, but reviewer 가 *therapeutic implication* 으로 over-read 가능. 추가 *NCT00085293, NCT01065090 retrospective trials are not therapeutic validation studies* 명시 권장.

---

## 7. Final safe title recommendation (table §15)

| # | Title | Safety |
|---|-------|--------|
| 1 | "A thyroid-lineage differentiation axis stratifies thyroid cancer beyond canonical driver mutations" | ★★★ safest; pure scientific framing; no clinical claim |
| 2 | "A driver-orthogonal thyroid-lineage axis recapitulates advanced-disease dedifferentiation in thyroid cancer" | ★★ — orthogonality 강조 |
| 3 | "A transcriptional differentiation axis defines lineage silencing in thyroid cancer" | ★★★ — silencing mechanism 강조 |
| 4 | "A compact RAI-lineage transcriptomic readout for early triage of post-surgery RAI-failure biology in thyroid cancer" | ★ — clinical reframe; "early triage" reviewer 에 따라 over-read 위험 |

**Recommended primary:** Title #1 또는 #3. Title #4 는 cover-letter 용 또는 hypothesis-section 용.

---

## 8. Next action

**User-written Hook required.**

본 paper 의 framing 은 이제 reviewer-proof. 다음 단계:
1. 본인이 직접 voice-protected sections 작성 (Hook / Aim / Discussion §3.1 / Limitations / Cover letter Para 1 / Reviewer Q9)
2. (Optional) Sub-A/sub-B 재산출 결과 (DEEP-DIVE 2) 를 manuscript supplementary 에 어떻게 통합할지 결정
3. Title #1 또는 #3 중 final 선택

---

**Paper 1 one-page audit hardened. Ready for user-written Hook.**
