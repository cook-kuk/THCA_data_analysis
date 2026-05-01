# Web Claude review request — Paper 1 manuscript v8 Prompt 1 deliverables

**Author:** Seungho Cook (1st author)
**Corresponding:** 유형원 교수
**Target venue:** Cell Reports Medicine (1순위) → Nature Communications (reach) → npj Precision Oncology (fallback)
**What I need:** 두 번째 의견. Title 3개 후보 / Abstract 147-word / Outline 에 대해 critical review + 본인 voice 적용 전 paper-level coherence sanity check.

---

## 1. Project context (web Claude 가 모르는 부분)

### 1.1 What this paper is

8-gene panel (SLC5A5, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, DIO1 — RAI 분화 transcript axis) 으로 BRAF/RAS-negative PTC dark matter 를 sub-stratify 한 paper. DM1/DM2 cluster 발견. **DM1 이 critical**:

- **L1 genetic:** DM1 76.8% tyrosine-kinase-fusion-positive (RET/NTRK/ALK/BRAF), OR 7.41 vs DM2
- **L2 mechanism heterogeneity:** DM1 sub-A (n=72, 84.7% fusion+, age 37, less advanced) vs sub-B (n=19, 57.9% fusion+, age 51, immune-hot Hashimoto-overlap)
- **L3 epigenetic:** DM1 promoter hypermethylation TPO Cohen's d=2.30, DIO1 d=1.24, TSHR d=1.20 — fusion-independent (within-DM1 fusion+ vs - methylation NS, d=−0.36)

**Combined paradigm:** DM1 = "fusion-driven + epigenetically-silenced dark matter"

Pooled OS HR (TCGA + MSK meta) = 2.53 [1.31, 4.89], I²=0%.

### 1.2 5-Pillar evidence backbone

| # | Pillar | Evidence |
|---|---|---|
| 1 | Korean Pan-Asian HLA cohort n=874 | K2+Lee+GSE286332-PTC arcasHLA, DPB1*05:01 53.2%, Chen 2018 GD forest replication |
| 2 | GSE286332 PTC vs PTC+HT molecular dissection | 10,380 DEGs, 8-gene Cohen d=−1.60, HLA-II d=+3.65, IFN-γ FDR=2e-4 |
| 3 | Driver mRNA neutrality | BRAF mRNA d=−0.04 (NS, V600E vs WT), all driver AUC <0.61, Driver_anchor cluster ARI=0 |
| 4 | Pan-genome cluster robustness | TIERA67 ARI=0.90 ≈ pan-genome top-5000 ARI=0.92, 8-gene alone ARI=0.49 |
| 5 | Autoimmune-PTC mechanism layer | TCGA Hashimoto-like 18-30% (DM2 OR up to 5×, p=6e-10), HLA-II 140% mediation, BCR clonal + TLS d=+1.96 + AICDA up, sub-B = NBNR cluster |

### 1.3 Key clinical claims

- **Reflex testing algorithm:** DM1 RNA-score positive → reflex RET fusion NGS panel. DM1 captures 81.8% of TCGA RET-fusion-positive tumors.
- **Population estimate:** 48 selpercatinib-eligible per 1000 PTC (LIBRETTO-001 alignment)
- **HMA + RAI re-induction trial rationale** (decitabine + I-131): TPO/DIO1/TSHR promoter hypermethylation → re-activation possible
- **East-Asian generalizability:** 4,300+ EA PTC tumors (Korean GSE213647 n=632 + GSE286332 + K2 + Lee + Wang Shanghai 2,844 + Liu 583 + MSK)

### 1.4 Cohort overview

- Discovery: TCGA-THCA n=504 (with OS), WGS+RNA-seq+miRNA+HM450 methylation (n=503)+SV (n=542 cBioPortal)
- Validation: K2 PRJEB11591 n=260, Lee n=632 GSE213647, MSK-IMPACT n=117 (Landa 2016 Cell, +12 SV), GSE286332 n=18 (autoimmune-PTC), 분당 SNUH outreach 단계
- Single-cell: Pu 2021 Nat Commun GSE184362 n=6 PTC (PRIMARY), Lu 2023 GSE193581 n=23 (suppl), GSE241184 (Phase 1, n=1)

### 1.5 What's already done vs not

✅ **Done (42+ analytical angles, 31 paper-shaping findings, 67 charts):**
- TCGA + MSK meta survival
- DM1 fusion paradigm (R3-F4) + robustness (R4-1)
- DM1 promoter hypermethylation (R5-2) — paradigm-shifting epigenetic layer
- 5-Pillar evidence assembled
- arcasHLA Korean Pan-Asian forest
- GSE286332 mediation analysis (HLA-II 140%)
- TCGA Hashimoto-like signature transfer (D4-P2)
- DM1 sub-B = NBNR cluster (D6-P7)
- Korean GSE213647 replication 22.8% (D8-B)

🟡 **Not yet done (limitations):**
- Bundang prospective Korean cohort (outreach 단계, 0%)
- K2 mini-index calibration FAIL (D7-P3) — alternate evidence chain holds
- Korean fusion direct calling (raw fastq + STAR-Fusion 1-2주 deferred)
- Single-probe-per-gene methylation aggregate (multi-probe deeper analysis 추후)

---

## 2. Prompt 1 deliverable A — Title candidates

### Reference (current v4 baseline, NOT recommended)

> "An 8-gene RAI-responsiveness biomarker reveals a young-onset fusion-driven epigenetically-silenced actionable subtype within BRAF/RAS-negative papillary thyroid carcinoma"

Char count: 174 (over by ~50). Four adjectives stacked = title overload.

### ★ Candidate 1 — Mechanism-forward (99 chars)

> **"Fusion-driven, epigenetically silenced dark matter in BRAF/RAS-negative papillary thyroid carcinoma"**

| Aspect | Note |
|---|---|
| **Hook** | Dual-mechanism (fusion + methylation) in one phrase; "dark matter" cites Xing 2014 NEJM |
| **Foregrounds** | R3-F4 (76.8% fusion+) + R5-2 (TPO d=2.30 hypermethylation) — Tier-1 paper-defining findings |
| **Misses** | No panel mention; no actionability signal upfront |
| **Pros** | Cell-Press style brevity; both mechanism layers; "dark matter" reviewer-known framing |
| **Cons** | "Dark matter" jargon needs 1-line abstract definition; less actionable upfront |

### Candidate 2 — Translational-forward (110 chars)

> **"An 8-gene RAI panel enables reflex fusion testing and re-induction therapy in BRAF/RAS-negative thyroid cancer"**

| Aspect | Note |
|---|---|
| **Hook** | Clinical algorithm (DM1+ → fusion NGS + HMA/RAI re-induction); 48/1000 selpercatinib-eligible implied |
| **Foregrounds** | R4-3 (DM1 captures 81.8% RET+) + R5-2 → HMA rationale |
| **Misses** | Dark-matter framing dropped; epigenetic mechanism not visible |
| **Pros** | Clinician-friendly; panel + therapy upfront; matches Cell Rep Med "actionable" emphasis |
| **Cons** | "Re-induction therapy" forward-looking (rationale, not trial); reads as biomarker paper |

### Candidate 3 — Discovery-forward (109 chars)

> **"An 8-gene panel sub-stratifies BRAF/RAS-negative thyroid cancer into fusion-driven and immune-overlap subtypes"**

| Aspect | Note |
|---|---|
| **Hook** | DM1 vs DM2 + DM1 sub-A vs sub-B heterogeneity foregrounded |
| **Foregrounds** | DM1/DM2 split + R4-2 fusion+/− heterogeneity + Pillar 5 autoimmune-PTC |
| **Misses** | Epigenetic layer not visible; HR/clinical aggressiveness not signaled |
| **Pros** | Honest about heterogeneity; "immune-overlap" foreshadows Pillar 5; flat discovery tone |
| **Cons** | Reads as biomarker; loses meta HR 2.53 signal |

### Decision matrix (Claude code 추천)

| Criterion | Cand 1 (Mech) | Cand 2 (Trans) | Cand 3 (Disc) |
|---|---|---|---|
| Cell Rep Med fit | ★★★ | ★★★ | ★★ |
| Nat Commun reach | ★★★ | ★★ | ★★ |
| npj Precision Onc fallback | ★★ | ★★★ | ★★★ |
| Clinical actionability signal | low | high | low |
| Mechanism signal | high | low | medium |

---

## 3. Prompt 1 deliverable B — Abstract (147 words, Cell Press structured)

**Background.** BRAF- and RAS-negative papillary thyroid carcinoma (PTC), the "dark matter" representing ~23% of cases, lacks mechanistic sub-stratification, hampering radioiodine (RAI) treatment decisions.

**Methods.** We applied an 8-gene RAI-responsiveness panel (SLC5A5, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, DIO1) to TCGA-THCA (n=504), MSK-IMPACT thyroid (n=117), and Korean PTC cohorts (K2/Lee/GSE286332, n=874), with single-cell validation (Pu 2021; Lu 2023). cBioPortal structural variants (n=542) and HM450 methylation (n=503) characterized cluster heterogeneity.

**Results.** The panel resolves a DM1/DM2 split. DM1 is 76.8% tyrosine-kinase-fusion-positive (RET/NTRK/ALK/BRAF; OR 7.41 vs DM2), captures 81.8% of TCGA RET-fusion tumors, and harbors fusion-independent promoter hypermethylation of differentiation genes (TPO Cohen's d=2.30). Pooled DM1 overall-survival hazard is 2.53 [1.31, 4.89] (TCGA + MSK meta). The axis is panel-free reproducible (pan-genome top-5000 ARI 0.92 vs TIERA67 0.90) and driver-mRNA-neutral (BRAF transcript d=−0.04).

**Conclusions.** DM1 is a fusion-driven, epigenetically silenced subtype, supporting a clinical 8-gene reflex algorithm for fusion-targeted therapy and a rationale for hypomethylating-agent + RAI re-induction trials.

### Word breakdown

| Section | Words |
|---|---|
| Background | 20 |
| Methods | 42 |
| Results | 60 |
| Conclusions | 25 |
| **Total** | **147** |

### Alternative Conclusions (forward-looking phrasing 후퇴 옵션)

- (a) [현재] "...supporting a clinical 8-gene reflex algorithm for fusion-targeted therapy and a rationale for hypomethylating-agent + RAI re-induction trials."
- (b) "...providing a clinical sub-stratification algorithm and motivating evaluation of fusion-targeted therapy and epigenetic-targeted RAI re-induction."
- (c) "...with immediate utility for reflex fusion testing and forward implications for hypomethylating-agent re-induction strategies."

---

## 4. Prompt 1 deliverable C — Outline (single-page)

### Section 1 · Introduction (600-900 words)

| # | Sub | Length | Content |
|---|---|---|---|
| 1.1 | Clinical context | 150-200 w | PTC global incidence ↑; RAI ablation over-treatment concern; Bethesda III/IV indeterminate FNA unmet need; BRAF V600E status alone insufficient |
| 1.2 | Existing molecular framework | 150-200 w | TCGA 2014 BRAF-like/RAS-like spectrum; Yoo 2016 BRS 273-gene; 23% BRAF-/RAS- gap |
| 1.3 | Dark matter concept | 200-250 w | Xing 2014 NEJM; BRAF-/TERT- 23-37%; East Asian 37.8% (Korean) vs 28.4% (TCGA); existing panels insufficient |
| 1.4 | Aim + brief preview | 100-150 w | 8-gene RAI panel; DM1/DM2 axis orthogonal to BRAF/RAS framework; 3-layer DM1 pathology; multi-cohort validation + clinical actionability |

### Section 2 · Results (2,500-4,000 words, 5 sub-results × 500-800 w)

| # | Sub-result | Pillars | Key claims | Figures |
|---|---|---|---|---|
| 2.1 | 8-gene panel + DM1/DM2 cluster discovery | P3 + P4 | TIERA67 7-category curation; drivers excluded by design; 8-gene = TDS_core subset; ΔAUC 8 vs 16 = 0.013 NS; pan-genome top-5000 ARI 0.92 ≈ TIERA67 0.90; 8-gene alone 0.49; Driver_anchor ARI=−0.007 | Fig 1 (5 panels) |
| 2.2 | Clinical aggressiveness + Meta | (Tier 1) | Xing 73% rescue (180 BRAF-/TERT- → 131 sub-stratified); TCGA HR 2.30 + MSK HR 2.67; Pooled HR 2.53 [1.31, 4.89] I²=0% | Figs 2, 6 |
| 2.3 | ★ DM1 fusion paradigm | (R3-F4) | DM1 76.8% fusion+: RET 33, NTRK 10, ALK 4, BRAF 5; OR 7.41; MSK SV cross-validation; missingness MAR (chi² p=0.56) | Fig 7 A-C |
| 2.4 | ★ DM1 heterogeneity + epigenetic layer | P5 | Sub-A vs sub-B (age, fusion%, immune); promoter hypermethylation TPO d=2.30; mean β 0.385 vs 0.253; methylation fusion-independent; 3-layer model | Figs 7D, 8 |
| 2.5 | Cross-cohort + Reflex algorithm | P1 + translational | sc external (Pu 2021, Lu 2023); FFPE robust; Korean n=874 + GSE213647 n=632; DM1 captures 81.8% TCGA RET+; 48 selpercatinib-eligible/1000 PTC | Figs 3, 5 |

**대안 ordering:** 2.1 → 2.3 (fusion paradigm 먼저) → 2.4 → 2.2 → 2.5 (mechanism→clinical narrative arc 강화) — 본인 review point.

### Section 3 · Discussion (1,000-1,500 words)

| # | Sub | Length | Content |
|---|---|---|---|
| 3.1 | Three-layer DM1 pathology synthesis | 300-400 w | L1 genetic + L2 heterogeneity + L3 epigenetic; compare to Krishnamoorthy 2025 PDTC/ATC Nat Comm |
| 3.2 | Clinical actionability + HMA + RAI re-induction | 300-400 w | 57% DM1 actionable; 81.8% RET+ capture; 48/1000 selpercatinib-eligible; HMA + RAI re-induction trial rationale |
| 3.3 | East-Asian generalizability | 200-300 w | 4,300+ EA tumors; Korean DM 37.8% vs TCGA 28.4%; K2 NBNR mixed phenotype |
| 3.4 | Limitations | 200-300 w | 28-point honest disclosure → 6-paragraph condensed |

### Figures (Cell Press main 8 + Suppl 6)

| Fig | Title | Highlight |
|---|---|---|
| 1 | 8-gene panel definition + DM1/DM2 discovery | Sankey + UMAP + heatmap + driver neutrality + ARI bar |
| 2 | Clinical aggressiveness within Xing rescue | Sankey 73% + KM |
| 3 | Single-cell validation | Pu 6 patients + Lu thyrocyte-intrinsic |
| 4 | Mutation × TERT × outcome | 8-cell decomposition + KM |
| 5 | Korean validation + FFPE | GSE213647 + K2 NBNR + FFPE robust |
| 6 | META forest | TCGA + MSK pooled 2.53 |
| **7 ★** | DM1 mechanism paradigm | Sub-A/B silhouette + fusion enrichment + partner stack + phenotype + reflex 81.8% |
| **8 ★** | Epigenetic layer (NEW R5-2) | Per-gene β heatmap + mean β bar + fusion-independent + HMA schematic |

---

## 5. 본인 (paper 저자) 의 voice 보호 영역

이 부분은 본인이 직접 작성:

- Title 1개 선택 (또는 본인 hybrid)
- Abstract 본인 voice 적용
- Intro 1.4 hook 첫 단락
- Discussion 3.1 첫 paragraph (paper 의 진짜 mechanism story)
- Discussion 3.4 Limitations (정직 disclosure 정신)
- Cover letter 첫 줄

Claude code 는 draft + section structure + Methods 표현 + Reference 정리 도구.

---

## 6. ★ Web Claude 에 묻는 specific questions

너무 generic 한 "review this" 보다는 specific 답을 받고 싶음:

### Q1. Title 선택
세 후보 중 **Cell Reports Medicine 1순위 + Nature Communications reach** 전략에 가장 fit 한 것?
본인이 hybrid 도 가능 — "Mechanism-forward + 8-gene panel mention" combo 같은 새 hybrid 후보를 제안할 의향?
**판단 근거 + 1-2줄 alternative phrasing 도 같이.**

### Q2. Abstract 강도
- "Re-induction trials" 가 너무 forward-looking 인가? Reviewer 가 reject 할 risk?
- Conclusions (a)/(b)/(c) 중 어느 것이 venue-appropriate?
- 5-Pillar 가 "panel-free reproducible + driver-mRNA-neutral + Korean replication" 으로 weaving 되어 있는데, **explicit 5 pillar 명시** vs **implicit weaving** — 어느 쪽이 Cell Rep Med abstract 표준?
- "dark matter" 단어 abstract 에 있어도 OK? (Background 첫 줄)

### Q3. Outline narrative arc
- Results 5 sub-section 순서:
  - 현재: 2.1 panel → 2.2 clinical agg → 2.3 fusion → 2.4 epi → 2.5 cross
  - 대안 A: 2.1 → 2.3 → 2.4 → 2.2 → 2.5 (mechanism-first → clinical-last)
  - 대안 B: 2.1 → 2.2 → 2.3 → 2.4 → 2.5 (현재, biomarker-first → mechanism → clinical)
- 어느 ordering 이 reviewer flow 좋은가?

### Q4. Pillar 5 (autoimmune-PTC) 위치
- 현재 outline 에서 Pillar 5 (HLA-II 140% mediation, BCR clonal, sub-B NBNR) 가 Results 2.4 에 통합. 단독 Section 으로 분리해서 6 sub-results (2.1-2.6) 로 가야 하는가?
- 단독 분리 시: 8-gene cancer paper 의 정체성이 흐려져서 reviewer 가 "scope creep" 으로 인식할 risk.
- 통합 유지 시: Pillar 5 autoimmune story 가 묻힐 risk.
- Web Claude 의 panel 균형 판단?

### Q5. Limitation framing
28-point honest disclosure → 6-paragraph condensed 가 Cell Rep Med 적절한가?
- 너무 많으면 paper weakness 누적
- 너무 적으면 reviewer 가 "limitations 너무 가볍다" 비판

권장 paragraph 수 + 어느 limitation 을 main 에서 빼고 supp 으로 옮길 수 있는지.

### Q6. Figure 7 + 8 의 game-changer 강도
- Fig 7 (DM1 mechanism, 5 panels): R3-F4 fusion 76.8% + R4-2 sub-A/B + R4-3 reflex 81.8%
- Fig 8 (epigenetic, 4 panels): R5-2 TPO β 2.30 + 8-gene β bar + fusion-independent + HMA schematic
- 이 두 figure 가 Cell Rep Med editor desk-reject 면제 + Nature Communications reach 가능 강도?
- 어느 panel 이 약하면 빼고 supplementary?

### Q7. Venue ladder strategic call
- Cell Rep Med 1순위 → Nat Commun reach → npj Precision Oncology fallback 의 ladder 가 이 evidence base (5-Pillar + 31 findings + 4,300+ EA + meta HR 2.53) 에서 합리적인가?
- 더 reach 가능 (Nat Med? Cancer Discov?) or 더 conservative (JCI Insight?) 권장?
- Bundang 분당 prospective 0% 가 최대 venue 결정 risk — bundang 없이 Cell Rep Med 통과 확률 객관적 판단.

### Q8. Sanity check — 빠진 내용
이 outline 에 빠진 critical evidence / framing / clinical claim 이 있는가? (예: Krishnamoorthy 2025 paired cancer paradigm, ATA 2015 guideline alignment, FDA companion diagnostic pathway 등)

---

**본인 Korean 으로 답변 받으면 좋고, English 도 OK.**

이 review 받은 후 → Prompt 2 (Introduction 600-900 words draft) 진행.
