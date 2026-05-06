---
title: "Paper 1 manuscript v8 — Abstract (Cell Press structured, ~150 words)"
date: 2026-04-30
revised: 2026-04-30 (v2 — web Claude review 통합)
author: Seungho Cook
target_venue: Cell Reports Medicine (1순위) → JCI Insight + Nat Commun dual reach → npj Precision Oncology (fallback)
word_count: 153 (target 150 ± 5)
status: v2 — Conclusions (b) 채택, Methods 40w / Results 66w 재분배, "dark matter" 따옴표 제거
---

# Abstract (Cell Press structured, 153 words)

**Background.** BRAF- and RAS-negative papillary thyroid carcinoma (PTC), the dark matter representing ~23% of cases, lacks mechanistic sub-stratification, hampering radioiodine (RAI) treatment decisions.

**Methods.** We applied an 8-gene RAI-responsiveness panel (SLC5A5, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, DIO1) to TCGA-THCA (n=504), MSK-IMPACT (n=117), and Korean cohorts (n=874). External single-cell validation (Pu 2021; Lu 2023), cBioPortal structural variants (n=542), and HM450 methylation (n=503) characterized cluster heterogeneity.

**Results.** The panel resolves a DM1/DM2 split. DM1 is 76.8% tyrosine-kinase-fusion-positive (RET/NTRK/ALK/BRAF; OR 7.41 vs DM2), captures 81.8% of TCGA RET-fusion tumors, and harbors fusion-independent promoter hypermethylation of differentiation genes (TPO Cohen's d=2.30; mean 8-gene β 0.385 vs 0.253). Pooled DM1 overall-survival hazard is 2.53 [1.31, 4.89] (TCGA + MSK meta). Candidate-pool restriction does not artificially impose this structure (pan-genome top-5000 ARI 0.92 vs TIERA67 0.90); BRAF transcript is mutation-status-neutral (d=−0.04).

**Conclusions.** DM1 is a fusion-driven, epigenetically silenced subtype, providing a clinical sub-stratification algorithm and motivating evaluation of fusion-targeted therapy and epigenetic-targeted RAI re-induction.

---

## Word-count breakdown (v2)

| Section | Words v1 | Words v2 | Cell Press 표준 |
|---|---|---|---|
| Background | 20 | 22 | 20-30 |
| Methods | 42 | 40 | 35-40 ★ web Claude 권장 |
| Results | 60 | 66 | 60-70 ★ web Claude 권장 (R5-2 강도 보강) |
| Conclusions | 25 | 25 | 20-30 |
| **Total** | **147** | **153** | **150 ± 10** |

## v1 → v2 변경 사항

| 변경 | 근거 (web Claude) |
|---|---|
| `"dark matter"` → `dark matter` (no quotes) | Academic abstract 에서 따옴표 처리 awkward. Xing 2014 NEJM cite 가 reviewer 익숙. |
| Conclusions (a) `"rationale for HMA + RAI re-induction trials"` → (b) `"motivating evaluation of fusion-targeted therapy and epigenetic-targeted RAI re-induction"` | (a) 가 너무 forward-looking — reviewer 가 "trial 결과 없는데 trial rationale 만 결론?" 의심 risk. (b) 는 정직 + venue-safe. |
| Methods 42w → 40w (Korean cohort 명시 단순화) | Cell Press Methods 표준 35-40w. K2/Lee/GSE286332 cohort 명세는 STAR Methods 에서 풀음. |
| Results 60w → 66w (R5-2 epigenetic 강도 보강 + Methods-level neutrality 명시) | "mean 8-gene β 0.385 vs 0.253" 추가 — TPO d=2.30 단독 claim 보다 정량 contrast 강함. "Candidate-pool restriction does not artificially impose this structure" 추가 — driver-mRNA-neutrality + pan-genome reproducibility 가 reviewer 의 "panel trivial?" Q 사전 차단 (web Claude Q8 catch). |

## Numerical claims — sourced

| Claim | Source |
|---|---|
| BRAF/RAS-neg PTC ~23% | TCGA 2014 Cell + Xing 2014 NEJM |
| TCGA-THCA n=504 (with OS) | Source v4 §9 Cohorts |
| MSK-IMPACT n=117 | Landa 2016 Cell |
| Korean cohorts n=874 | D4-P1: K2(235) + Lee(630) + GSE286332-PTC(9) |
| cBioPortal SV n=542 | R3-F4 (542/557 SV-tested) |
| HM450 methylation n=503 | R5-2 |
| DM1 fusion 76.8% (63/82) | R3-F4 |
| OR 7.41 vs DM2 | R3-F4 |
| DM1 captures 81.8% RET+ | R4-3 (27/33) |
| TPO Cohen's d=2.30 | R5-2 |
| Mean 8-gene β 0.385 vs 0.253 | R5-2 (DM1 vs DM2) |
| Pooled HR 2.53 [1.31, 4.89] | R2-N1 (TCGA + MSK meta) |
| Pan-genome top-5000 ARI 0.92 | P4 |
| TIERA67 ARI 0.90 | P4 |
| BRAF transcript d=−0.04 | P1 (V600E vs WT) |

## 본인 review focus (v2)

1. **"Dark matter" 단어** — 따옴표 제거 적용. 본인 호불호 한 번 더 검토.
2. **Conclusions (b)** — "motivating evaluation of fusion-targeted therapy and epigenetic-targeted RAI re-induction" — 정직 톤. 본인 voice 적용 시 "evaluation" → "investigation" 또는 "exploration" 도 고려.
3. **5-Pillar implicit weaving 유지** — explicit "5-Pillar evidence" 명시 안 함 (web Claude 권장). Pillar 1 (Korean cohort) + 3 (driver neutrality) + 4 (pan-genome) 자연스럽게 woven. **Companion immune-overlap axis 는 abstract 에서 의도적 보류 — Paper 2 backbone 으로 reserve.**
4. **R5-2 강도** — TPO d=2.30 단독에서 "mean 8-gene β 0.385 vs 0.253" 정량 contrast 추가로 강화.

## Alternative Conclusions (이전 옵션 보존)

- (a) [폐기] "...supporting a clinical 8-gene reflex algorithm for fusion-targeted therapy and a rationale for hypomethylating-agent + RAI re-induction trials." — forward-looking risk
- **(b) [채택] "...providing a clinical sub-stratification algorithm and motivating evaluation of fusion-targeted therapy and epigenetic-targeted RAI re-induction."**
- (c) [예비] "...with immediate utility for reflex fusion testing and forward implications for hypomethylating-agent re-induction strategies." — npj Prec Onco 시 톤 다운 옵션
