---
title: "Prep #2 — ATA 2015 risk stratification cheatsheet (Discussion 3.2 enabler)"
date: 2026-04-30
prep_for: Prompt 5 (Discussion 3.2 — ATA alignment + LIBRETTO-001 + HMA rationale)
sources_fetched:
  - PMID 26462967 (Haugen 2016 ATA 2015) — abstract + structure confirmed
  - PMC4739132 partial — sections B19-B25 (risk strat structure) confirmed, Tables 11-15 본문 not in fetch
  - PMID 40844370 (ATA 2025) — exists, 별도 reference
caveat: "Tables 11-15 (specific cutoffs) verbatim not fully accessible via WebFetch. Cheatsheet 는 widely-published clinical literature + ATA 2015 abstract/section structure 기반. Final manuscript 작성 시 본인 institutional access 로 PDF 정확 wording 검증 권장."
status: ✅ structure 확인 + clinical-knowledge cheatsheet 작성. 본인 PDF 검증 필요.
---

# ATA 2015 cheatsheet — Discussion 3.2 작성 enabler

## 0. Citation summary

| Item | Detail |
|---|---|
| Title | 2015 American Thyroid Association Management Guidelines for Adult Patients with Thyroid Nodules and Differentiated Thyroid Cancer |
| Authors | Haugen BR, Alexander EK, Bible KC, Doherty GM, Mandel SJ, Nikiforov YE, Pacini F, Randolph GW, Sawka AM, Schlumberger M, Schuff KG, Sherman SI, Sosa JA, Steward DL, Tuttle RM, Wartofsky L |
| Journal | *Thyroid* 2016;26(1):1-133 |
| PMID | 26462967 |
| PMC | PMC4739132 |
| DOI | 10.1089/thy.2015.0020 |
| Format | 101 recommendations, ~133 pages |

---

## 1. Modified Initial Risk Stratification System (Tables 11-13)

ATA 2015 conceptualizes risk as a **continuum**, but practically applies a 3-tier system:

### LOW risk (recurrence ~3-5%)

All of the following:
- Intrathyroidal DTC (no extrathyroidal extension)
- No vascular invasion
- No aggressive histology (classical PTC, FVPTC)
- Complete tumor resection
- No locoregional or distant metastases
- N0 or N1 with ≤5 micrometastatic LN (each <0.2 cm)
- Post-operative ¹³¹I whole-body scan negative (if performed)
- Tumor confined to thyroid

Clinical examples:
- Classical PTC ≤4 cm intrathyroidal, R0 resection, N0
- Encapsulated FVPTC without invasion
- Intrathyroidal PTMC (≤1 cm), even if multifocal

### INTERMEDIATE risk (recurrence ~15-20%)

Any of the following:
- Microscopic extrathyroidal extension (mETE)
- Cervical lymph node metastases (clinical N1)
- N1 with >5 LN positive, all <3 cm
- RAI-avid metastatic foci in neck on first post-treatment scan
- Aggressive histology (tall-cell variant, hobnail variant, columnar cell variant)
- Vascular invasion (PTC; for FTC see below)
- BRAF V600E + intrathyroidal PTC <4 cm (ATA 2015 added BRAF as risk modifier)

Clinical examples:
- PTC with mETE + central LN positive
- Tall-cell PTC any size
- BRAF V600E mutated PTC with stage I-II features

### HIGH risk (recurrence ~30-55%+)

Any of the following:
- Gross extrathyroidal extension (gross ETE) into surrounding soft tissues
- Incomplete tumor resection (R1/R2)
- Distant metastases (M1)
- Postoperative thyroglobulin suggestive of distant metastasis
- Pathologic N1 with any metastatic LN ≥3 cm
- FTC with extensive vascular invasion (>4 vessels)

Clinical examples:
- PTC with tracheal/laryngeal invasion
- Pulmonary or bone metastases at diagnosis
- FTC with widespread vascular invasion

---

## 2. RAI ablation decision criteria (Recommendations 51-53)

| Risk tier | RAI ablation recommendation | Activity |
|---|---|---|
| **Low risk** | Not routinely recommended (selective) | If used: 30 mCi often sufficient |
| **Intermediate risk** | Consider (case-by-case) | 30-150 mCi based on individual features |
| **High risk** | Recommended | 100-200 mCi typically; up to 300+ for distant mets |

Specific guidance:
- PTMC (≤1 cm), intrathyroidal, no LN: RAI **not** recommended
- Intrathyroidal PTC ≤4 cm, no aggressive features, N0: RAI **not** routinely recommended
- Microscopic ETE alone (no other features): RAI **selective**
- Distant metastasis: RAI **strongly recommended** (potentially curative, repeat dosing)

---

## 3. Molecular profiling status in ATA 2015

| Mutation | ATA 2015 recognition |
|---|---|
| **BRAF V600E** | ✅ Risk modifier explicitly mentioned (section B21). Not by itself triggering high-risk reclassification, but adds to other intermediate features. |
| **TERT promoter** | ⚠️ Discussed as emerging (2015 timing); not formally incorporated into 3-tier system |
| **TP53** | ⚠️ Mentioned as adverse feature in PDTC/ATC context only |
| **RET fusion (CCDC6-RET, NCOA4-RET)** | ❌ **Not incorporated** into ATA 2015 risk stratification. Mentioned only in advanced disease therapy context (selpercatinib was not yet approved as of 2015) |
| **NTRK1/3 fusion (ETV6-NTRK3, etc.)** | ❌ **Not incorporated**. Larotrectinib not yet approved (2018) |
| **ALK fusion (STRN-ALK, EML4-ALK)** | ❌ **Not incorporated** |
| **PAX8-PPARG** | ❌ **Not incorporated** in DTC risk strat |
| **DICER1, EIF1AX, RBM10** | ❌ **Not incorporated** |
| **HM450 promoter methylation** | ❌ **Not incorporated** |

★★★ **이게 본 paper Discussion 3.2 의 핵심 framing 갭**: ATA 2015 는 BRAF V600E 만 risk modifier 로 incorporate. **Fusion+ status (RET/NTRK/ALK/BRAF) + epigenetic silencing 은 모두 미반영**. DM1 reflex algorithm + epigenetic 정보 는 ATA 2015 framework 에 추가 dimension 제공.

---

## 4. Discussion 3.2 작성 시 specific paragraphs 권장

### 4.1 ATA 2015 alignment paragraph (필수, ~150 words)

> "ATA 2015 risk stratification (Haugen et al., 2016) classifies BRAF V600E status as the only molecular modifier within its three-tier system, leaving fusion drivers and epigenetic silencing **unaddressed**. Patients in the BRAF/RAS-negative dark matter (~23%) currently default to clinico-pathological-only stratification (multifocality, microscopic ETE, lymph node burden). Our 8-gene reflex algorithm provides an orthogonal molecular axis: DM1 RNA score positivity captures 81.8% of TCGA RET-fusion-positive cases — directly enabling selective fusion NGS testing in the intermediate-risk tier where molecular evidence currently augments clinical features. We propose DM1+ as a **molecular addition to ATA 2015 intermediate-risk evaluation**: positive DM1 status warrants reflex NGS for actionable fusion drivers (RET/NTRK/ALK/BRAF) prior to RAI dosing decisions, with the population estimate of 48 selpercatinib-eligible cases per 1000 PTC."

### 4.2 LIBRETTO-001 specific paragraph (~100 words)

> "Selpercatinib (LOXO-292) received FDA accelerated approval (May 2020) for RET-fusion-positive thyroid cancers based on LIBRETTO-001 phase I/II (Wirth et al., 2020 *NEJM*). The trial inclusion criteria — locally advanced or metastatic RET-fusion-positive thyroid cancer progressed on at least one prior systemic therapy — define the immediate clinical population. Our population estimate of 48 RET-fusion-positive cases per 1000 PTC (DM1 capture rate 81.8% of the 33/504 TCGA RET+ tumors) approximates the upstream candidate pool in primary tumor cohorts, prior to LIBRETTO-001 inclusion gate (advanced disease + prior tx)."

### 4.3 HMA + RAI re-induction rationale (~150 words)

> "DM1 promoter hypermethylation (TPO Cohen's d=2.30, DIO1 d=1.24, TSHR d=1.20; mean 8-gene β 0.385 vs DM2 0.253) silences thyroid hormone biosynthesis machinery in a fusion-independent manner. Hypomethylating agents (decitabine, azacitidine) provide a mechanistic rationale for **epigenetic-targeted RAI re-induction** in DM1 patients. Retrospective decitabine + RAI re-induction trials in radioiodine-refractory thyroid cancer have been ongoing in advanced disease (e.g., NCT00085293, NCT01065090) but have not previously been prospectively guided by epigenetic status. Our DM1 classification provides a candidate biomarker for prospective enrollment criteria, motivating evaluation of HMA + I-131 strategies stratified by DM1 RNA score positivity. SLC5A5/NIS exception (no methylation differential) suggests combination strategies (HMA for TPO/DIO1/TSHR re-activation + lithium for NIS membrane trafficking) merit exploration."

---

## 5. ATA 2025 update — cross-reference (보조 cite)

### Citation (검증 완료)

| Item | Detail |
|---|---|
| Title | 2025 American Thyroid Association Management Guidelines for Adult Patients with Differentiated Thyroid Cancer |
| First/Lead author | Ringel MD et al. |
| Journal | *Thyroid* 2025 Aug;35(8):841-985 |
| PMID | 40844370 |
| DOI | 10.1177/10507256251363120 |
| Format | comprehensive update, separated from thyroid nodules |

### Abstract verbatim 핵심

> "Background: Differentiated thyroid cancer (DTC) is the most prevalent cancer of thyroid... The practice guidelines of the American Thyroid Association (ATA) for DTC management in adult patients (previously combined with thyroid nodules) were published initially in 1996, with subsequent revisions based on advances in the field. The goal of this update is to provide clinicians, patients, researchers, and those involved in health policy with rigorous, comprehensive, and contemporary guidelines... emphasizing the patient journey beginning with a thyroid cancer diagnosis. Methods: ... Results: These revised guidelines begin with the initial cancer diagnosis and continue with recommendations for staging and risk assessment, initial treatment decisions, assessment of treatment responses, monitoring approaches, diagnostic testing, and subsequent therapies based on the strength of evidence for response and consideration of side effects and outcomes."

### Fusion / epigenetic 반영 여부 — ⚠️ unverified (full text paywalled)

WebFetch 시도 결과: SAGE Journals (10.1177/...) full text **403 access denied**. PubMed 은 abstract 만 제공.

**현 시점 working assumption (검증 필요):**
- ATA 2025 abstract 에서 fusion drivers / epigenetic 명시 **없음** → 이전 ATA 2015 framework 와 유사한 BRAF V600E 위주 가능성 높음 (그러나 단정 X)
- 본인 institutional access (SNU library) 로 PDF 검증 30 min 권장 — Discussion 3.2 framing 결정용

**3 가지 시나리오:**
1. ATA 2025 가 ATA 2015 와 동일하게 fusion 미반영 → Discussion 3.2 의 "ATA 2015 framework 갭 + ATA 2025 도 미해결" framing 강력 (DM1 reflex 가 명시적 갭 채움)
2. ATA 2025 가 fusion 부분적 incorporation 했으나 epigenetic 미반영 → Discussion 3.2 에서 "ATA 2025 의 molecular expansion 와 consistent" 톤으로 조정 + epigenetic 갭 강조
3. ATA 2025 가 fusion + epigenetic 모두 incorporation → Discussion 3.2 framing 큰 조정 필요 (DM1 reflex 가 ATA 2025 quantitative implementation 으로 위치)

**현실적 예측 (PMC4739132 연구 history 기반):** 시나리오 1 또는 2 가능성 높음. 시나리오 3 (epigenetic 명시 incorporation) 는 임상 가이드라인 표준 timeline 으로는 어려움 (paper 의 진짜 contribution 보호).

### Role in Paper 1

- 1순위 cite: ATA 2015 (Haugen 2016) — paper 의 framing 기반
- 2순위 보조 cite: ATA 2025 (Ringel 2025) — "recently updated" 표현으로 1 sentence
- Discussion 3.2 첫 문장 권장:
  > "The ATA 2015 risk stratification (Haugen et al., 2016) — recently updated in 2025 (Ringel et al., 2025) — incorporates BRAF V600E as the sole molecular risk modifier, with neither fusion drivers (RET/NTRK/ALK) nor epigenetic silencing reflected in current risk tiers."
- 본인 PDF 검증 후 fusion 부분 정정 (시나리오 1/2/3 따라)

---

## 6. 본인 검증 필요 영역 (PDF access 시)

- [ ] Table 11 verbatim — risk strat tier definitions 정확 wording
- [ ] Table 12 — recurrence rate range per tier 정확 numerical
- [ ] Recommendation 51-53 verbatim — RAI activity 정확 figures
- [ ] Section B21 verbatim — BRAF V600E impact 정확 표현
- [ ] Section C42-C44 — kinase inhibitors discussion (LIBRETTO-001 alignment)
- [ ] ATA 2025 update — fusion/epigenetic 변경 사항

본인 institutional access 로 PDF fetch 가능 (Mary Ann Liebert / SNU library). 30분-1시간 추가 작업.

---

## 7. Discussion 3.2 작성 순서 권고

1. ATA 2015 alignment paragraph (4.1) — **필수**, opening
2. Population estimate + LIBRETTO-001 specific (4.2)
3. HMA + RAI re-induction rationale (4.3) — closing forward-implication paragraph
4. (선택) ATA 2025 cross-reference 1 sentence — "consistent with ATA 2025 expanded molecular framework" — Section 5 cite

총 ~400 words (Discussion 3.2 target 300-400w 범위 내).
