# Web Claude 3차 review 요청 — Full manuscript v1 (compiled)

**Author:** Seungho Cook
**Status:** Manuscript v8 v1 first-pass complete. 6-week marathon plan single-session 압축 draft. 본인 voice 적용 (5/5-5/10 W1) + cite verify (5/4 본인 read 후) + 공저자 review (W6) → bioRxiv preprint.

**Total main text:** 5,514 words (Cell Rep Med 표준 5,500-7,000 within range)
- Abstract 153 + Intro 715 + Results 3,250 + Discussion 1,380 + Limitations 280

---

## 0. 이전 review 와의 연결

| Review | Outcome |
|---|---|
| 1차 (v1 outline) | 88/100 → v2 적용 (5 권장) |
| 2차 (v2 outline + prep) | 95/100 → v3 적용 (7 권장 + 5 추가 cite + Landa misattribution 정정 + reverse-causality 3-layer 차단) |
| **3차 (이번)** | **Full v1 manuscript review — 본인 voice 적용 전 paper-level coherence + Cell Rep Med 적합성 final check** |

이번 3차는 본인 voice 적용 + 공저자 review 전 last-pass 점검.

---

## 1. Manuscript full v1 (compiled)

전체 manuscript 는 `10_full_manuscript_v1.md` 에 compile 됨. 약 5,500w, Title + Abstract + Intro 4 sub + Results 6 sub + Discussion 4 sub + Limitations + Methods + Cover letter + 12 Reviewer Q&A.

### Title (final candidate)

> "An 8-gene panel reveals fusion-driven, epigenetically silenced dark matter in BRAF/RAS-negative thyroid cancer"

(Hybrid Candidate, 108 chars, web Claude 1차+2차 review 권장)

### Abstract (153 words, Cell Press structured, Conclusions b)

[`01_abstract.md` 참조 — 이전 review 95/100 confirm 받은 v2]

### Section 2 Results — ★ 대안 A ordering (mechanism-first)

1. **2.1** 8-gene panel + DM1/DM2 cluster discovery (~620w) — Pillars 3+4
2. **2.2** DM1 fusion paradigm (~700w) ★ — R3-F4 76.8% + R5-1 MSK
3. **2.3** DM1 epigenetic silencing (~580w) ★ — R5-2 fusion-independent
4. **2.4** Fusion-negative DM1 immune-overlap teaser (~440w) — Paper 2 reserve
5. **2.5** Clinical aggressiveness Meta (~440w) — N1 HR 2.53 [1.31, 4.89]
6. **2.6** Cross-cohort + Reflex algorithm (~470w) — 4,300+ EA + 81.8% RET+ capture

### Section 3 Discussion (4 sub × 280-430w)

- **3.1** Three-layer DM1 + Landa 2016 dial-back + Yoo 2016 convergence + Bradley 2010 contradict
- **3.2** ATA 2015+2025 dual + Wirth 2020 LIBRETTO-001 + HMA rationale
- **3.3** East-Asian + Pu 2021 sc validation
- **3.4** Limitations (4 main paragraphs + 2 supp)

---

## 2. ★ Web Claude 에 묻는 7 specific questions

### Q1. Coherence — Title/Abstract/Hook 정합성

Title "An 8-gene panel reveals fusion-driven, epigenetically silenced dark matter in BRAF/RAS-negative thyroid cancer"
Abstract Conclusions: "DM1 is a fusion-driven, epigenetically silenced subtype, providing a clinical sub-stratification algorithm and motivating evaluation of fusion-targeted therapy and epigenetic-targeted RAI re-induction."
Hook (Alt B refined): "Despite a >98% 5-year overall survival in differentiated thyroid carcinoma, structural disease recurrence reaches ~20% in ATA 2015 intermediate-risk patients (Haugen et al., 2016), and BRAF/RAS-negative tumors — accounting for ~23% of cases — complicate radioiodine (RAI) treatment decisions in the absence of mechanistic sub-stratification."

세 element 가 narrative arc 정합? Hook → "no mechanistic axis" → Aim/Conclusions "DM1 fusion-driven + epigenetically silenced sub-stratification" 으로 연결? 또는 missing link?

### Q2. Results 대안 A ordering (mechanism-first) 평가

이전 review 권장 따라 2.1 → 2.2 (fusion paradigm 먼저) → 2.3 (epigenetic) → 2.4 (immune-overlap teaser) → 2.5 (clinical aggressiveness Meta HR 2.53) → 2.6 (cross-cohort reflex). Cell Press editorial 선호 패턴 맞는가?

특히: **HR 2.53 meta 가 2.5 (5번째)** — biomarker 의 진짜 clinical claim 이 narrative 끝에 위치. 너무 늦은가? 또는 mechanism story 충분히 build 후 clinical claim 이 자연스러운가?

### Q3. Section 2.4 (immune-overlap teaser) Paper 2 reservation 강도

> "A full characterization of the autoimmune-PTC mechanism axis — including Pan-Asian HLA repertoire, B cell receptor clonal architecture, tertiary lymphoid structure burden, and mediation analyses — is reserved for a separate study (Cook et al., manuscript in preparation, Paper 2)."

이게 reviewer 가 받아들일 수준? 또는 "scope creep / under-developed sub-section" 으로 인식 risk? 더 deep 또는 더 brief?

### Q4. ★ Discussion 3.1 framing — Landa dial-back tone

> "In the primary BRAF/RAS-negative PTC compartment, our DM1 framework **identifies an upstream signature consistent with this dedifferentiation trajectory**: the same differentiation machinery that is silenced at the ATC end is already epigenetically attenuated in DM1 PTCs..."

이전 v1 ("extends paired-cancer continuum upstream") 보다 dial-back. Cell Rep Med editor 받아들일 정도 강도? 또는 더 confident 톤 권장?

### Q5. ★ Reverse-causality 3-layer 차단 강도

8-gene 5/8 = Landa 2016 ATC silenced gene list overlap.

**Layer 1** Methods Yoo 2016 first cite: "independently of and prior to access of advanced-disease ATC-silenced gene lists (Landa et al., 2016)"

**Layer 2** Discussion 3.1 convergence: "two independent paths to the same differentiation axis...biological convergent validation rather than panel-driven circular inference"

**Layer 3** Cover letter explicit: "panel was independently derived from canonical RAI biology and predates our access to Landa et al. (2016)'s ATC-silenced gene list"

3-layer 차단 충분한가? 또는 추가 layer 필요 (예: timeline document 첨부, GitHub commit history evidence)?

### Q6. Cover letter strategy — 5 suggested reviewers

- Capdevila (Vall d'Hebron) — selpercatinib LIBRETTO-001 author
- Landa (Brigham/Harvard) — advanced thyroid genomics
- Parangi (MGH) — thyroid surgery + clinical actionability
- Park Young Joo (SNU) — Korean PTC genomics
- Nikiforov (UPMC) — thyroid molecular pathology

Excluded: Krishnamoorthy/Fagin lab (MSK, recent paired papers).

5명 균형 OK? Park Young Joo 의 인지도 (Yu professor lab과 같은 SNU 환경) 가 reviewer 적합성 affect? Cell Rep Med 표준 reviewer 풀에서 빠진 후보 있는가?

### Q7. Sanity check — 빠진 paper-level critical evidence

Manuscript v1 first-pass 완료 후 빠진 paper-level critical evidence:

후보 catch:
- **Methods Q&A consistency**: Reviewer Q5 (SV missingness MAR) + Q9 (epigenetic mechanism) 답변과 Methods § 정합성
- **Cover letter Para 3 ("what's new")** 의 mechanism + actionability 균형 톤
- **Suppl Tables S1-S10** 가 paper claim 의 모든 numerical evidence cover 하는가?
- **Figure 7 + 8 게임체인저** 의 panel 구성 (5+3) 이 Cell Rep Med editor 첫 리뷰에서 강력 인상 줄 정도?
- **East-Asian generalizability claim** (4,300+ tumors) — 정확 cohort breakdown:
  - TCGA n=504 (EUR-majority, NOT EA)
  - MSK n=117 (EUR-majority, advanced disease enriched)
  - K2 n=260 (Korean)
  - Lee n=632 (Korean)
  - GSE286332 n=18 (Korean)
  - GSE184362 Pu n=6 (Chinese, Fudan)
  - Lu 2023 n=23 (Chinese)
  - Wang 2024 (Chinese, Shanghai n=2,844 mention only)
  - Liu 2017 (Asian n=583 mention only)

  → **East-Asian 4,300+ 가 정확? 또는 over-claim?** 정확 EA-specific cohort sum 검증 필요.

빠진 cite 또는 framing?

---

## 3. 본인 5/4-5/10 진행 plan

| 작업 | 시간 | 우선순위 |
|---|---|---|
| Landa 2016 JCI 직접 read | 1.5 hr | ★★★ |
| Yoo 2016 PLOS Genet review | 30 min | ★★★ |
| ATA 2015 PMC + ATA 2025 verify | 1.5 hr | ★★ |
| Wirth 2020 NEJM (LIBRETTO ORR) | 30 min | ★★ |
| Bradley 2010 verify | 15 min | ★ |
| Yu professor 미팅 + 분당 outreach 합의 | 1 hr | ★★★ |
| v2.5 미세 조정 (Hook/Title verb/Fig 1 panel) | 30 min | ★★★ |
| 본인 voice 적용 (Hook + Aim + Disc 3.1 + Limitations + Cover letter Para 1) | 2-3 hr | ★★★ |
| Web Claude 3차 review 받고 v2 정정 | 30 min | ★★★ |

총 prep ~7-8 hr (5/4-5/5 분산).

5/11-5/17 W2 부터 본격 manuscript writing (v1 → v2 → v3 iteration).

---

## 4. 한 줄 정리

**Manuscript v8 v1 first-pass complete (5,514w main text within Cell Rep Med 표준). 본인 voice 적용 + cite verify + 공저자 review 전 last paper-level coherence check 7 question (Title/Abstract/Hook 정합 + 대안 A ordering + Paper 2 reservation 강도 + Landa dial-back tone + reverse-causality 3-layer + Cover letter reviewers + 빠진 critical evidence). 본인이 며칠 전 직접 그린 marathon 정신 그대로.**

본인 Korean 답변 권장.
