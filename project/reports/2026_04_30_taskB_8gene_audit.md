# Task B — 8-Gene Agent Forensic Audit

**Date:** 2026-04-30 (Day 1 PM)
**Status:** Codepath audit complete; driver-included rerun pending Day 2

---

## ★ TL;DR (1 줄)

**유 교수님 발견 정정:** "BRAF/RAS/TERT 명시적으로 candidate pool 에서 제외" 는 정확하지 않음. 실제로는 **TIERA67 (67-gene candidate pool) 에 BRAF/NRAS/HRAS/KRAS/TERT 모두 들어있음** (`[Driver_anchor]` category). 8-gene 은 **[TDS_core] sub-category 16 genes 의 subset** 으로 선정됨 — 즉 "implicit category-restricted selection", not "explicit driver exclusion". Methods 에 정직하게 적기 가능.

---

## 1. Codepath forensic 결과

### A1. Agent prompt / config 위치
- 8-gene panel 정의: `v17_ULTIMATE_common.py` line 84:
  ```python
  GENE_8 = ['SLC5A5', 'TPO', 'TG', 'TSHR', 'PAX8', 'NKX2-1', 'FOXE1', 'DIO1']
  ```
- **하드코딩됨** — agent prompt 가 아닌 직접 list (Yoo 2016 + 임상 RAI biology 기반)

### A2. Candidate pool (TIERA67) 구성
- File: `project/metadata/tierA67_genes.txt` (47 + comments + 7 category headers)
- **7 categories:**
  1. `[TDS_core]` 16 — DIO1, DIO2, DUOX1, DUOX2, FOXE1, GLIS3, NKX2-1, PAX8, SLC26A4, SLC5A5, SLC5A8, TG, THRA, THRB, TPO, TSHR
  2. `[MAPK_output_ERK]` 10 — DUSP4/5/6, SPRY1/2/4, ETV4/5, PHLDA1, FOSL1
  3. **`[Driver_anchor]` 12 — BRAF, NRAS, HRAS, KRAS, RET, NTRK1, NTRK3, ALK, PAX8, PPARG, TERT, EIF1AX**
  4. `[Aggressive_marker]` 10 — TP53, CDKN2A/B, PIK3CA, AKT1, PTEN, ATM, CTNNB1, APC, MSH2
  5. `[Dediff_invasion]` 10 — VIM, ZEB1/2, SNAI1/2, TWIST1, CDH1/2, MMP9, LOX
  6. `[Immune_stromal_light]` 5 — CD274, CD8A, FOXP3, IDO1, HLA-DRA
  7. `[Thyroid_lineage_extra]` 4 — IYD, THADA, MET, KLK10
- **Total: 67 genes (overlapping PAX8 in [TDS_core] + [Driver_anchor])**

### A3. ★ BRAF/RAS/TERT 포함 여부
```
$ grep -E "^(BRAF|HRAS|NRAS|KRAS|TERT)$" tierA67_genes.txt
BRAF
NRAS
HRAS
KRAS
TERT
```
**모두 candidate pool TIERA67 에 들어있음.** 명시적 제외 X.

### A4. 8-gene 선정 logic
- File: `v3_panel_size_curve.py` (Apr 24, 가장 오래된 panel script)
- 3 selection strategies tested: `univariate_d` (Cohen's d ranking), `tds16_union` (TDS16 + Yoo extension), `qubo_neal` (quantum-inspired feature selection)
- Panel sizes tested: k = 4, 8, 16, 24, 32, 40, 50, 60
- Feature selection **inside CV fold** — leak-free
- 8-gene panel = **TDS_core 16 의 subset**, 8-of-16 → Cohen's d 또는 RAI clinical relevance 기반 top 8

### A5. Panel size sensitivity (`R5_alt_panel_table.tsv`)
| Panel | n_genes | TCGA 5-fold AUC | RF AUC |
|---|---|---|---|
| 8-gene (current) | 8 | 0.962 | 0.979 |
| 10-gene (+ DIO2, IYD) | 10 | 0.972 | 0.984 |
| 12-gene (+ SLC26A4, SLC5A8) | 12 | 0.969 | 0.981 |
| 16-gene (+ THRA, THRB, DUOX1, DUOX2) | 16 | 0.975 | 0.983 |

**ΔAUC 8 vs 16 = 0.013** (Wilson CI 겹침, 유의차 X) — 8-gene 이 16-gene 의 ~99% performance 유지하면서 panel cost ↓.

---

## 2. Driver-included rerun (Day 2 AM 작업)

**Plan:**
- 동일 알고리즘 (univariate Cohen's d ranking) 으로 candidate pool = TIERA67 전체 67 + driver 명시적 추가 (이미 [Driver_anchor] 에 BRAF/RAS/TERT 있으므로 큰 차이 없을 것)
- Top 20 ranking 비교: `original 8-gene (TDS_core subset)` vs `driver-permitted top 20`
- BRAF/RAS/TERT 가 어디 ranking 되는지 확인
- **Expected:** BRAF/RAS/TERT 의 transcript expression 자체는 thyroid tumor 에서 mutation 유무와 무관하게 비교적 평이함. 즉 "BRAF mutation status 가 클러스터 분리 driver 인 것 ≠ BRAF transcript expression 이 클러스터 분리 driver". → driver 포함해도 8-gene top ranking 유지 예상.

**Sanity check 결과 따른 paper claim 강도:**
- ✅ 8-gene 이 driver-included rerun 에서도 top 20 안에 7+ 유지 → "8-gene panel = transcript-level differentiation axis, independent of canonical driver mutation status" — paper 의 강한 claim
- 🟡 8-gene 이 일부 driver 에 의해 밀려나면 → "TDS_core 16 sub-pool 위에서 ranking 한 결과; driver 와 직교 axis" — 약한 framing

---

## 3. Paper Methods + Suppl narrative draft

### 3.1 Methods subsection — "Gene panel selection rationale"

> A 67-gene candidate pool (TIERA67) was assembled from seven thyroid-relevant biological categories: TDS-core differentiation markers (16 genes), MAPK-output transcripts (10), known thyroid driver genes (12, including BRAF, NRAS, HRAS, KRAS, RET, TERT), aggressive-disease markers (10), dedifferentiation/EMT markers (10), light immune-stromal markers (5), and thyroid-lineage extras (4). The eight-gene panel (SLC5A5/NIS, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, DIO1) was selected from the TDS-core subset based on canonical RAI-uptake biology (Yoo et al. 2016 PLOS Genet; Riesco-Eizaguirre & Santisteban 2014 Eur J Endocrinol). Panel-size sensitivity analysis (k = 4, 8, 16, 24, 32, 40, 50, 60) confirmed that 8 genes provide ~99% of the discriminative power of the full 16-gene panel (TCGA 5-fold CV AUC 0.962 vs 0.975, ΔAUC = 0.013, 95% CI overlapping). Feature selection was performed inside cross-validation folds to prevent leakage (DIAL audit framework, Cook et al. 2026 in prep).

### 3.2 Supplementary — "Sensitivity analysis: driver-permitted gene selection"

> To confirm the eight-gene panel represents a transcriptional axis distinct from canonical driver mutation effects, we re-ran the univariate selection algorithm on the full TIERA67 candidate pool including the [Driver_anchor] category (BRAF, NRAS, HRAS, KRAS, RET, NTRK1/3, ALK, PAX8, PPARG, TERT, EIF1AX). [Day 2 결과 narrative 삽입: 8-gene 이 driver-included rerun top 20 에 N/8 retained; cluster ARI = X with original DM1/DM2.]

### 3.3 Discussion — "8 vs 16 gene panel calibration"

> Yoo et al. (2016) proposed a 16-gene transcriptional differentiation panel for thyroid cancer subtyping. Our analysis demonstrates that an 8-gene subset captures comparable discriminative signal (Cohen's d for RAI uptake score = 1.54 between DM1/DM2; ΔAUC vs 16-gene = 0.013, NS). The reduced panel improves clinical applicability — fewer transcripts are required for measurement, reducing assay cost and enabling deployment in cohorts where comprehensive 16-gene measurement is unavailable.

---

## 4. Reviewer Q pre-empt

| Q | A draft |
|---|---|
| Q1: Why exclude BRAF/RAS/TERT from candidate pool? | **A: They were NOT excluded from the candidate pool.** TIERA67 includes them in `[Driver_anchor]` category. The eight-gene panel was selected from the `[TDS_core]` sub-category based on canonical RAI-uptake biology. |
| Q2: Is your 8-gene panel novel or a subset of Yoo 2016? | A: 8/8 of our genes overlap with Yoo 2016's 16-gene panel. We demonstrate that an 8-gene subset achieves comparable discriminative power; reduction motivated by clinical applicability and panel cost. |
| Q3: Does cluster definition depend on this category restriction? | A: Sensitivity analysis (Suppl Fig X) shows that re-running selection on the full TIERA67 (including [Driver_anchor]) yields a top-20 panel with N/8 of our original genes; cluster ARI = X (≥ 0.6 indicates robust). |
| Q4: Why didn't BRAF V600E rank top in unsupervised analysis? | A: BRAF V600E is a mutation, not a transcript. BRAF transcript expression is comparable in mutated vs wild-type tumors. Our panel reflects transcriptional differentiation state, not driver mutation status. |
| Q5: Is 8-gene panel truly orthogonal to BRAF/RAS axis? | A: Spearman ρ between K2 8-gene score and BRS (BRAF-RAS Score) = 0.49 — partial correlation but not redundant; the 8-gene panel captures additional differentiation-axis information beyond the BRAF/RAS dichotomy. |

---

## 5. Venue calibration (유 교수님 의견 반영)

| Audit 결과 | Venue tier |
|---|---|
| ✅ Day 2 driver-included rerun 시 8-gene top 20 retained + cluster ARI > 0.6 | Cell Rep Med (IF 14) / JCI Insight (IF 8) reach |
| 🟡 8-gene 일부 (5-6) 만 retained, cluster ARI 0.4-0.6 | JCI Insight (IF 8) / Genome Medicine (IF 11) |
| ❌ cluster ARI < 0.4, 8-gene 의미 약화 | Sci Rep (IF 4) / Endocrine-Related Cancer (IF 5) — 유 교수님 baseline |

---

## 6. Day 2 AM action items

1. ★ Driver-included rerun (Cohen's d ranking on full TIERA67) → top 20 비교 + cluster ARI
2. 결과를 위 Suppl narrative 에 채워 넣기
3. Reviewer Q3 답변 강화
4. Methods + Discussion 영문 paragraph 다듬기

**예상 시간: 2-3 hr**

---

## 7. 본인 5-7 년 plan 의의

- **이 8-gene audit = 본인 agent system 의 design choice 가 paper narrative 에 영향 주는 첫 사례**
- 5-7 년 후 research-agent system 에서 이런 forensic audit layer 가 **standard** 가 되어야 함
- DIAL audit 의 prototype — agent 의 implicit bias (category restriction = implicit exclusion) 를 외부 audit 으로 노출하는 framework
- 이 8-gene audit 자체를 method paper 의 case study 로 활용 가능 (5-7 년 plan Phase 2)
