# Web Claude 2차 review 요청 — Prompt 2 (Introduction draft) 시작 전 Prep work 검증

**Author:** Seungho Cook
**Status:** v2 outline 95/100 confirm 받음. Prep work (Hook stat 검증 + ATA 2015 cheatsheet + Krishnamoorthy/Landa cite 검증) 완료. Prompt 2 (Introduction 600-900 w) 진입 전 본인 confirm 4개 + critical finding 1개 web Claude 판단 요청.

---

## 0. 이전 review 와의 연결

이전 review 에서 권장한 5건 모두 v2 적용 완료:
1. Hybrid Title Cand 1 (108 chars)
2. Abstract Conclusions (b)
3. Results 대안 A ordering (2.1 → 2.3 → 2.4 → 2.2 → 2.5)
4. Section 2.4 split (2.4a epigenetic + 2.4b immune-overlap teaser, Pillar 5 Paper 2 reserve)
5. Discussion 3.1 Krishnamoorthy + 3.2 ATA 2015/LIBRETTO-001 + Limitations 4+2 + venue ladder 정직화

이전 review 끝에서 권장한 prep 3건:
- Item 1: Hook stat ("5-15% recur") 검증
- Item 2: ATA 2015 (Haugen 2016) cheatsheet 정리
- Item 3: **Krishnamoorthy 2025 Nat Comm 본인 직접 read** (framing tone 결정)

이 3건 prep 결과 — 이번 review 의 대상.

---

## 1. Prep #1 — Hook stat 검증 결과

### 본인 v2 outline example hook

> "Despite ~95% 5-year survival, 5-15% of papillary thyroid carcinomas recur or metastasize, and clinicians continue to make radioiodine decisions on heterogeneous BRAF/RAS-negative tumors with no mechanistic compass."

### 검증 결과

| 본인 표현 | Verdict |
|---|---|
| "~95% 5-year survival" | ✅ 정확 (실은 98-99%, 5-yr OS DTC overall) |
| "5-15% recur or metastasize" | ⚠️ 부분 정확 — overall 5-20% (해외 lit), 한국 single-center 4.3% conservative. **Intermediate-risk specifically ~20%** (ATA 2015) |
| "no mechanistic compass" | ✅ 본인 voice — 검증 영역 아님 |

### 3 alternative hook 작성 (citation 별)

**Alt A (ATA 2015 conservative):**
> "Although papillary thyroid carcinoma (PTC) carries an excellent overall prognosis (5-year disease-specific survival >98%; ATA 2015), structural disease recurs in 5-20% and distant metastases occur in approximately 10% of patients (Haugen et al., 2016 *Thyroid*), and clinicians continue to make radioiodine decisions on heterogeneous BRAF/RAS-negative tumors without a mechanistic basis."

**Alt B (Risk-stratified specific, ATA framing thread):**
> "Despite a >98% 5-year overall survival in differentiated thyroid carcinoma, structural disease recurrence remains 5-20% across risk tiers — reaching ~20% in ATA 2015 intermediate-risk patients (Haugen et al., 2016) — and the heterogeneity of BRAF/RAS-negative papillary thyroid carcinoma continues to complicate radioiodine treatment decisions in the absence of mechanistic sub-stratification."

**Alt C (본인 v2 voice 보존, 수치만 정정):**
> "Despite a >98% 5-year survival, 5-20% of papillary thyroid carcinomas recur or develop distant metastases (Haugen et al., 2016), and clinicians continue to make radioiodine decisions on heterogeneous BRAF/RAS-negative tumors without a mechanistic compass."

---

## 2. Prep #2 — ATA 2015 cheatsheet 핵심 finding

### Haugen 2016 *Thyroid* 26(1):1-133 (PMID 26462967)

3-tier risk system (continuum 으로 visualization):
- **Low (~3-5% recur)**: intrathyroidal DTC, no ETE, no aggressive histology, N0 or N1 ≤5 micromets <0.2cm
- **Intermediate (~15-20% recur)**: microscopic ETE, N1 with >5 LN, aggressive histology (tall-cell/hobnail/columnar), vascular invasion, **BRAF V600E + intrathyroidal**
- **High (~30-55%+ recur)**: gross ETE, distant mets, R1/R2 resection, N1 LN ≥3cm, FTC extensive vascular invasion (>4 vessels)

### ★★★ Discussion 3.2 의 핵심 framing 갭 — ATA 2015 의 molecular profiling 부재

| Mutation | ATA 2015 incorporation |
|---|---|
| BRAF V600E | ✅ 유일하게 risk modifier (section B21) |
| TERT promoter | ⚠️ 2015 timing 으로 emerging only, 미공식 |
| **RET fusion** | ❌ 미반영 (selpercatinib 미승인 시점) |
| **NTRK1/3 fusion** | ❌ 미반영 (larotrectinib 2018 승인 후) |
| **ALK fusion** | ❌ 미반영 |
| **PAX8-PPARG** | ❌ 미반영 |
| DICER1, EIF1AX, RBM10 | ❌ 미반영 |
| **HM450 promoter methylation** | ❌ 미반영 |

→ **ATA 2015 framework 가 fusion+ 와 epigenetic silencing 둘 다 미반영** = 본 paper 의 정확한 unmet need 갭. DM1 reflex algorithm + epigenetic layer 가 이 갭을 직접 채움.

### 추가 발견

- **ATA 2025 Guidelines** 별도 published (PMID 40844370, *Thyroid* 2025). 본 paper Cell Rep Med submission 2026 시 reviewer 가 "ATA 2025 update 반영했는지" 물을 risk. Discussion 3.2 에서 ATA 2015 (1순위 cite) + ATA 2025 (보조 cite) cross-reference 권장.

### 본인 검증 필요 (PDF access)

Tables 11-15 verbatim wording 은 WebFetch 에서 추출 안 됨. 본인 institutional access (SNU library / Mary Ann Liebert) 로 30분-1시간 검증 권장. Cheatsheet 의 risk strat 내용은 widely-published clinical literature 기반으로 정확하지만, manuscript 의 quoted text 는 PDF verbatim 이 안전.

---

## 3. ★★★ Prep #3 — Krishnamoorthy/Landa cite — Critical misattribution 발견

### 이전 web Claude 권장

> "Discussion 3.1 + Krishnamoorthy 2025 PDTC/ATC Nat Comm cite (paired cancer continuum framing)"

### 검증 결과 — 4가지 모두 잘못

| 항목 | 이전 표기 | 실제 검증 |
|---|---|---|
| 첫 저자 | Krishnamoorthy | **Iñigo Landa** (Krishnamoorthy GP 는 9번째 공저자) |
| 연도 | 2025 | **2016** |
| 저널 | Nat Commun | **J Clin Invest (JCI)** |
| 내용 framing | "PDTC/ATC dark matter 정의" | ✅ paired-continuum 정확히 다룸 (표현은 "PDTCs and ATCs arise from well-differentiated tumors") |

### 진짜 cite — Landa 2016 JCI

> **Landa I, Ibrahimpasic T, Boucai L, Sinha R, Knauf JA, Shah RH, Dogan S, Ricarte-Filho JC, Krishnamoorthy GP, Xu B, Schultz N, Berger MF, Sander C, Taylor BS, Ghossein R, Ganly I, Fagin JA. Genomic and transcriptomic hallmarks of poorly differentiated and anaplastic thyroid cancers. *J Clin Invest.* 2016;126(3):1052-1066. doi:10.1172/JCI85271. PMID: 26878173.**

### 본 paper 와의 정합성 — ★ 원래 framing 보다 훨씬 강력

| Landa 2016 finding (verbatim 추출) | 본 paper 8-gene/DM1 finding | Discussion 3.1 framing power |
|---|---|---|
| **"ATCs have profoundly suppressed mRNA levels for TG, TSHR, TPO, PAX8, SLC26A4, DIO1, and DUOX2 genes"** | DM1 promoter hypermethylation TPO d=2.30, DIO1 d=1.24, TSHR d=1.20, PAX8 d=0.97, TG d=0.86, FOXE1 d=0.84, NKX2-1 d=0.63 (R5-2) | ★★★ **거의 1:1 gene overlap** (TG/TSHR/TPO/PAX8/DIO1 5/8 panel + SLC26A4/DUOX2 supplementary) |
| TERT promoter stepwise: 9% PTC → 40% PDTC → 73% ATC | 본 paper R3-F4 fusion finding 은 PTC primary tumor 에서 fusion-driven dark matter 확인 | Stepwise progression 신호 — DM1 = early stage in trajectory |
| "ATCs were BRAF-like irrespective of driver mutation" | DM1 = 76.8% fusion+; BRAF transcript d=−0.04 mutation-neutral | Convergence 패턴 — DM1 cluster 도 transcript-level convergence |
| **n=117 advanced thyroid (84 PDTC + 33 ATC) — MSK-IMPACT** | 본 paper Pillar 의 MSK-IMPACT n=117 (R5-1 SV cross-validation) | **같은 cohort** — 본 paper 가 이미 기반 cohort 사용 중 |

### ★ 권장 Discussion 3.1 framing (Cell Rep Med editor 강력)

> "Landa et al. (2016) characterized the genomic landscape of advanced thyroid cancer (84 PDTC + 33 ATC), establishing that PDTC and ATC arise from well-differentiated tumors through accumulated genetic abnormalities — including TERT promoter mutations that increase stepwise from 9% in PTC to 40% in PDTC and 73% in ATC, and a profound suppression of thyroid differentiation transcripts (TG, TSHR, TPO, PAX8, SLC26A4, DIO1, DUOX2) in ATC. **Our DM1 framework extends this paired-cancer continuum upstream into the primary PTC compartment**: the same differentiation machinery silenced at the ATC end of the trajectory is already epigenetically attenuated in DM1 PTCs (TPO Cohen's d=2.30, DIO1 d=1.24, TSHR d=1.20 vs DM2; mean 8-gene β 0.385 vs 0.253), even before TERT promoter mutations clonally expand. DM1 thus identifies — within the BRAF/RAS-negative dark matter — the early epigenetic prefiguration of the dedifferentiation trajectory Landa described in advanced disease."

→ Cell Rep Med editor 가 좋아할 framing: **continuum 의 upstream end 를 본 paper 가 정의 → Landa downstream 과 paired**.

### 보조 cite 정리

| Paper | Cite role | Fit |
|---|---|---|
| **★ Landa 2016 JCI 126(3):1052-1066** | Discussion 3.1 primary | ★★★ Canonical paired-continuum cite. n=117 cohort 정확 일치. 8-gene panel 의 5/8 직접 cite. |
| Nat Commun 2025 PDTC/ATC (s41467-025-58910-3) | Discussion 3.1 secondary 보조 | ★★ 다른 group (likely Wang/Fudan), 348 thyroid + Pro-I/II/III subtypes. 본인이 D3-P4 audit 에서 다룬 Wang 2024 와 같은 group 가능 — supporting cite only |
| Krishnamoorthy 2025 JEM RBM10 (PMID 39992626) | Discussion 3.2 또는 supplementary | ★ first-author 정확 (Krishnamoorthy GP MSKCC), 단 framing 은 metastasis splicing — Discussion 3.1 paired-continuum 과 부적합. Discussion 3.2 또는 limitations forward implication 에서 broad context cite |

### 본인 직접 read 권장 (이전 web Claude 권장 그대로)

| Paper | URL | 우선순위 | 시간 |
|---|---|---|---|
| **★ Landa 2016 JCI** | https://www.jci.org/articles/view/85271 | ★★★ 필수 | 1-1.5 hr |
| Nat Commun 2025 PDTC/ATC | https://www.nature.com/articles/s41467-025-58910-3 | ★★ optional | 30 min |
| Krishnamoorthy 2025 JEM | (PMID 39992626) | ★ optional | 30 min |

---

## 4. ★ Web Claude 에 묻는 specific questions

### Q1. Hook stat 3 alt 중 어느 것?

| Option | Tone | ATA 3.2 thread |
|---|---|---|
| Alt A | conservative, safest, Cell Rep Med 표준 | weak |
| **Alt B** | risk-stratified, intermediate-risk 20% explicit | **★ strong** (Discussion 3.2 ATA alignment 까지 thread) |
| Alt C | 본인 voice 보존 (수치만 정정) | medium |

본인 voice 보호 측면 vs Cell Rep Med editor 선호 — 어느 쪽?

### Q2. ★★★ Landa 2016 JCI cite 채택 confirm

이전 review 의 "Krishnamoorthy 2025 Nat Comm cite" 권장이 misattribution 으로 확인 (4가지 fact 모두 틀림). 진짜 cite 는 **Landa 2016 JCI** (Krishnamoorthy GP 9번째 공저자, MSK-IMPACT n=117 — 본 paper 가 이미 사용 중인 cohort). 다음 confirm:

- (a) Landa 2016 JCI 채택 OK?
- (b) Discussion 3.1 권장 framing ("DM1 extends Landa's paired-continuum upstream into BRAF/RAS-negative primary PTC compartment") 의 강도 OK? 또는 더 약하게/강하게?
- (c) 8-gene panel 의 5/8 (TG/TSHR/TPO/PAX8/DIO1) 가 Landa 2016 의 ATC silenced gene list (TG, TSHR, TPO, PAX8, SLC26A4, DIO1, DUOX2) 와 직접 1:1 overlap — 이거 framing power 가 reviewer 에게 강력해 보이는가? 또는 "panel selection 이 Landa 결과에 의존?" 이라는 reverse-causality Q risk 가 있는가?

### Q3. Nat Commun 2025 PDTC/ATC paper — Discussion 3.1 보조 cite 가치

본인이 이 paper 직접 read 할 가치 있는가 (30분 추가)? 아니면 generic 한 supplementary cite 로 충분한가?

특히 본인이 D3-P4 에서 audit 한 Wang/Fudan group (Wang 2024 Endocr Connect 13(11):e240301 patient-level supplementary X / FU 데이터 X / DICER1 미언급) 와 같은 group 일 가능성 — 이게 본 paper Discussion 3.1 에서 cite 시 신뢰도 issue?

### Q4. Krishnamoorthy 2025 JEM RBM10 처리

이 paper 는 Krishnamoorthy GP first-author 이지만 framing 다름 (RBM10 splicing → metastasis). 본 paper 에서:

- (a) 완전 제외
- (b) Discussion 3.2 또는 supplementary 에서 broad context only ("metastasis-specific mechanisms in advanced thyroid cancer")
- (c) Limitations forward implication ("DM1 framework prospective integration with Krishnamoorthy 2025 RBM10 axis merits future work")

어느 옵션?

### Q5. ATA 2025 Guidelines (PMID 40844370) cross-reference

본 paper Cell Rep Med submission 2026 시 ATA 2025 update 반영해야 하는가? (1) ATA 2015 1순위 + ATA 2025 보조 cite 권장 vs (2) ATA 2015 만 cite 후 reviewer Q 받으면 답변?

### Q6. v2.5 미세 조정 — Fig 1 panel count

이전 review 에서 web Claude 권장: "Fig 1 5 panels → 4 panels (Heatmap C → S1 supp 이동)". 본인 적용 vs 5 panels 유지 — 다시 한번 final 판단.

### Q7. Sanity check — 빠진 critical evidence

이번 prep 결과 이후 (Hook stat 검증 + ATA 2015 + Landa 2016 cite 확인), 본 paper Discussion 또는 Introduction 에 빠진 critical evidence/cite 가 추가로 있는가?

후보 catch 시도:
- Pozdeyev 2018 CCR (n=779 advanced DTC + ATC, 4-cluster ATC) — Landa 보다 더 광범위 cohort, Discussion 3.1 보조 cite?
- Yoo 2016 PLoS Genet (BRS 273-gene, 본인 8-gene 의 출처 paper) — Methods 외에 Introduction 1.2 에서 hard cite?
- TCGA 2014 Cell PTC paper — Introduction 1.2 의 BRAF-like/RAS-like spectrum 정의 cite?
- Pu 2021 Nat Commun (sc PTC primary external) — Methods 에서 cite, Discussion 에서도 cite?

---

## 5. 본인 직접 read 일정 — Prompt 2 진입 전

| 작업 | 시간 | 우선순위 |
|---|---|---|
| Landa 2016 JCI 직접 read + 5-question reading guide | 1-1.5 hr | ★★★ 필수 |
| ATA 2015 PDF Tables 11-15 verbatim 검증 | 30 min - 1 hr | ★★ 권장 |
| Nat Commun 2025 PDTC/ATC paper read | 30 min | ★ optional |
| Hook 3 alt 중 1개 선택 + Title verb 결정 + Fig 1 panel count 결정 | 30 min | ★★★ 필수 |
| Web Claude 답변 받고 confirm 정리 | 30 min | ★★★ 필수 |

총 prep 마무리 ~3-4 hr → Prompt 2 (Introduction 600-900 w draft) 진입.

---

## 6. 한 줄 정리

Prep #1+#2 정상. **Prep #3 에서 Landa 2016 JCI 가 진짜 cite 로 확인 — framing power 가 원래보다 훨씬 강력 (8-gene 의 5/8 직접 1:1 overlap, 같은 MSK-IMPACT n=117 cohort).** Web Claude 7개 specific question (Q1 Hook 선택 + Q2 Landa framing 강도 + Q3-Q5 보조 cite + Q6 Fig 1 + Q7 sanity check) 답변 받으면 본인 voice 적용 + Prompt 2 진입.
