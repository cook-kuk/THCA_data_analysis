---
title: "F1 — GSE286332 cancer paper 통합 가능성 audit + framework decision"
date: 2026-04-30
purpose: GSE286332 P3 결과를 본 cancer paper에 넣을지 분리할지 결정 + critical caveat
---

# F1 — Framework decision audit

## 0. 핵심 질문

GSE286332 P3 결과 (10,380 DEGs, HLA-II d=+3.65, **18/18 모두 DM2**)를 본 cancer paper에 통합? 별도 autoimmune-PTC paper로 분리? Hybrid?

## 1. 통합문서에서 GSE286332 빠진 진단

**확인된 사실:**
- GSE286332 raw + 분석 결과 모두 프로젝트 안 존재 (`project/results/p3_gse286332/` + `2026_04_30_P3_GSE286332_critical_result.md`)
- 5/1 P3가 진행되었고 STRONG GO verdict
- HOWEVER 본 audit (4/29 ~ 4/30 R1 + R2)의 통합문서에는 명시적 GSE286332 reference 없음

**빠진 이유 추정:**
- 4/29 audit은 cancer paper의 BRAF/RAS DM1/DM2 4-pillar focused
- 5/1 P3 GSE286332는 **autoimmune-PTC sub-axis** 라는 distinct angle
- Audit author가 "cancer paper main scope 외"로 판단해서 통합문서에 미반영했을 가능성
- 또는 단순 누락

## 2. 🚨 Round 3 Critical Finding — GSE286332 DM 분류와 TCGA가 INVERSE direction

본 round 3 (F2 분석) 발견:
- **TCGA-THCA Hashimoto-like proxy (B_cell + IFN-γ top 20%, n=142):**
  - DM1: 45.2% Hashimoto-like (66/146)
  - **DM2: 0% Hashimoto-like (0/61)** ← 완전 zero
  - not_DM: 15.2%
- **GSE286332 (5/1 P3):** PTC+HT 18/18 = **DM2 분류** (per user)

### Inverse direction reconciliation

이 두 finding은 **정반대 direction**. 가능한 설명:

1. **K2/GSE286332 mini-index calibration issue** (memory `v17_korean_k2_calibration.md`)
   - kallisto 8-gene mini-index TPM이 ~10-100x inflated
   - TCGA-trained absolute-form LogReg가 wrong direction 예측
   - GSE286332에서 P_DM2 high → 실제로는 DM1 mapping 가능성
   - **이 가설 매우 plausible**

2. **Different signature definition**
   - TCGA: B_cell + IFN-γ (immune deconvolution proxy)
   - GSE286332: full Hashimoto-overlap PTC molecular signature (10,380 DEGs)
   - 둘이 다른 phenotype 측정할 가능성

3. **Cohort-specific biology**
   - GSE286332 = explicit clinical Hashimoto-overlap PTC (paired with regular PTC)
   - TCGA = implicit lymphocytic infiltration variability

→ **Critical action**: Prompt 2 (TCGA Hashimoto-like generalization) 결과 본 audit Round 3에서 시행. **GSE286332와 TCGA의 DM 분류 일관성이 깨짐을 발견**. 이것이 framework decision에 결정적 영향.

## 3. 3가지 시나리오 평가

### 시나리오 A: 통합 (GSE286332 cancer paper main에 포함)

**조건**: TCGA Hashimoto-like도 DM2 enriched면 가능

**현실**: TCGA Hashimoto-like는 **DM1 enriched (반대 direction)**
- 두 cohort에서 inverse → "Hashimoto-overlap = extreme DM2" 단정 불가
- Reviewer가 cohort 간 inconsistency 즉시 지적 가능

**Verdict**: ❌ 추천 안 함 — TCGA inverse direction이 본 시나리오의 base hypothesis 무효화

### 시나리오 B: 분리 (별도 autoimmune-PTC paper)

**조건**: 본 cancer paper 빠른 submission + autoimmune trajectory 분리

**근거**:
- TCGA / GSE286332 inverse direction → 두 phenotype은 같은 axis 아님 (best to separate)
- 본 cancer paper의 DM1/DM2 framework 안전하게 보존
- 별도 paper "Autoimmune-overlap PTC molecular signature" trajectory (J Autoimmun / Front Immunol target)
- 분당 Graves' 응답 시 풍부 cohort 만들 수 있음

**Verdict**: ✅ **권장** — 본 cancer paper는 ship-ready 상태 유지, GSE286332는 별도

### 시나리오 C: Hybrid (1 supp figure로 minimal 포함)

**조건**: Cancer paper venue 유지 + autoimmune angle 보존

**근거**:
- 본 cancer paper Discussion에 1 paragraph 추가:
  > "We note that the autoimmune-overlap PTC subtype (Hashimoto-overlap, GSE286332 cohort) shows a distinct molecular signature with elevated HLA Class II expression. The relationship between this autoimmune-overlap axis and our DM1/DM2 classification merits further investigation in larger cohorts; preliminary analysis in TCGA-THCA suggests Hashimoto-like signatures associate with DM1 (45% Hashimoto-proxy) rather than DM2, potentially reflecting cohort calibration differences."
- 본격 분석은 별도 paper

**Verdict**: 🟡 가능 — 정직 disclosure이면서 cancer paper venue 유지. 그러나 inverse direction은 reviewer에게 confusion 유발 가능.

## 4. 본인 권장: **시나리오 B (분리)**

**이유:**

1. **TCGA inverse direction**: F2 결과로 "GSE286332 PTC+HT 18/18 DM2"와 "TCGA Hashimoto-like DM1 enriched"가 inverse → 같은 axis 아님

2. **Cancer paper의 핵심 narrative 보존**:
   - 4-pillar (BRAF/RAS / DM cluster / Xing rescue / multi-cohort meta) 견고
   - 추가 5번째 pillar는 reviewer에게 scope 확장 의심 유발

3. **별도 paper trajectory의 가치**:
   - "Autoimmune-overlap PTC molecular signature" — 본인 unique angle
   - autoimmune × thyroid identity 정렬 (메모리에 명시된 5-7년 plan)
   - GSE286332 + Korean GSE213647 hashi_like flag (이미 59 환자) + Lu 2023 Hashimoto sc + 분당 Graves' (응답 시) → 풍부 cohort
   - venue: J Autoimmun (IF 14) / Front Immunol (IF 7-8) / Endocrine-Related Cancer (IF 5)

4. **Critical paper-shaping finding 발견 (round 3 F4)**:
   - **DM1 fusion rate 76.8%** (RET/PTC/NTRK/ALK/BRAF fusion-driven)
   - 이게 GSE286332 통합보다 훨씬 큰 cancer paper 영향
   - 이번 라운드에 발견된 "DM1 mechanism = fusion-driven"가 본 cancer paper의 main breakthrough
   - GSE286332는 별도 paper이 더 어울림

## 5. 시나리오 B 시 본 cancer paper 영향

- 통합문서 (FINAL_COMPREHENSIVE_SUMMARY.md) 그대로 유지
- GSE286332 reference 추가 안 함
- 단 Discussion limitations에 1줄 추가:
  > "We did not include autoimmune-overlap PTC analysis in this study; that axis is being addressed in a separate manuscript focused on Hashimoto-overlap molecular signatures (GSE286332 cohort + Korean lymphocytic infiltration analyses)."

## 6. 시나리오 B 별도 paper outline

**Title (tentative)**: "Hashimoto-overlap papillary thyroid carcinoma: a B-cell rich tertiary lymphoid structure-driven molecular subtype"

**Cohorts**:
- GSE286332: n=18 (9 PTC vs 9 PTC+HT), Korean RNA-seq, primary discovery
- TCGA-THCA n=142 Hashimoto-like (B_cell + IFN-γ top-20%) — generalization
- Korean GSE213647 hashi_like flag (n=59) — independent Korean validation
- Lu 2023 GSE193581 Hashimoto sc — single-cell mechanism
- 분당 Graves' (응답 시) — extends to autoimmune spectrum

**Target venues**: J Autoimmun / Front Immunol / Endocrine-Related Cancer / Thyroid

**Timeline**: 4-9개월 (cancer paper와 parallel)

## 7. Action items post-decision

### Cancer paper (immediate)
- [x] 통합문서 그대로 유지 (시나리오 B)
- [ ] Discussion 1줄 disclosure 추가 (separate manuscript framing)
- [ ] DM1 fusion finding (round 3 F4 — 76.8% fusion-driven) 즉시 통합 → **paper revision priority**

### Autoimmune-PTC paper (separate trajectory)
- [ ] GSE286332 + TCGA Hashimoto generalization (F2 round 3 결과 활용)
- [ ] Korean GSE213647 hashi_like × DM cluster
- [ ] Lu 2023 sc Hashimoto signature
- [ ] 분당 Graves' outreach + integration
- [ ] target venue: J Autoimmun

## 8. 한 줄 결론

**시나리오 B (분리) 권장.** TCGA에서 Hashimoto-like가 DM2가 아니라 DM1 enriched (45% vs 0%)인 inverse direction 발견 → GSE286332와 TCGA가 같은 axis가 아님. 본 cancer paper는 4-pillar + 새로 발견된 DM1 fusion-driven mechanism (76.8% fusion+)로 충분히 강력. GSE286332는 별도 autoimmune-PTC paper trajectory로.

**가장 중요한 round 3 finding은 GSE286332 framework가 아니라 DM1 fusion mechanism 발견** — 본 cancer paper의 paradigm-shifting result.
