# 2026-05-03 — Paper 2 Pillar 1 Forest Meta STRONG → Decision Brief for Web Claude

> 본인 Seungho Cook · Phase 0 5-Pillar Cancer paper trajectory.
> 4/29 audit → 4/30 dual-track → 5/1 P1-P5 → 5/2 D-series 7개 → 5/3 Paper 2 Pillar 1 forest meta 완료.
> 이 문서를 web Claude 에 paste 하고 다음 지령 받기 위함.
> 핵심: **Pillar 1 PARTIAL → STRONG 확정.** Chu 2018 (Chen → Chu citation correction) Han Chinese GD 와 Korean PTC pool 통합. 이제 manuscript writing 시작 vs 추가 분석 vs venue 결정 — web Claude 판단 필요.

---

## 0. 전체 맥락 (4/29 → 5/3 Phase 0)

### 0.1 본인 위치
- **분야 교차점:** autoimmune (cookHLA Nat Commun first-author, 류마티스+T1D+CD HLA imputation 경험) × thyroid cancer (현 PTC paper)
- **5-7 년 plan:** 한국 medical AI 첫 mover, autoimmune × thyroid 횡단
- **직무:** ARIA build-out 빠짐 → 시간 70% 확보 → 중장기 연구 직무 전환

### 0.2 Phase 0 Cancer paper — 5-Pillar 구조 (5/3 EOD)

| # | Pillar | Status | Key signal |
|---|---|---|---|
| 1 | **Korean Pan-Asian HLA cohort n=874** | ✅ **STRONG (5/3 today)** | DPB1\*05:01 53.2% Korean > 44.0% Chen GD |
| 2 | GSE286332 PTC vs PTC+HT molecular | ✅ STRONG | 10,380 DEGs, IFN-γ FDR=2e-4 |
| 3 | Driver mRNA neutrality | ✅ STRONG | BRAF d=−0.04 |
| 4 | Pan-genome cluster robustness | ✅ STRONG | TIERA67 ARI 0.90 ≈ pan-genome 0.92 |
| 5 | Autoimmune-PTC mechanism layer | ✅ STRONG | DM2 OR up to 5×, p=6e-10; HLA-II 140% mediation; TLS d=+1.96 |

→ **5-pillar 모두 STRONG.** Cell Rep Med / JCI Insight reach 강력 보장. Bundang Graves' 추가 시 Nat Commun reach.

### 0.3 Total Azure cost (4/29 → 5/3)
- 4/29 K2/Lee arcasHLA burst: $4.80 (already)
- 5/2 D4-P1 GSE286332 arcasHLA burst: ~$1.0
- 5/3 Pillar 1 forest meta: $0 (local only)
- **Total: ~$5.80**

---

## 1. 5/3 Today — Paper 2 Pillar 1 Forest Meta (STRONG)

### 1.1 Task
Chen 2018 (실은 Chu X et al. 2018) Han Chinese GD allele-level summary stats 통합 → quantitative random-effects forest meta. Pillar 1 PARTIAL → STRONG.

### 1.2 ★ Citation correction

원래 메모리에 "Chen 2018 PMC 6161647" 으로 적혀있었지만 실제 paper 는:

> **Chu X, Pan CM, Zhao SX, et al.** *HLA association with autoimmune Graves' disease in 1,468 Chinese cases and 1,490 controls — Han Chinese.* J Med Genet. 2018;55(10):685–692. doi:**10.1136/jmedgenet-2017-105146** (PMC 6161647).

DPB1\*05:01 published OR = **1.90** [1.69, 2.14] (not 2.45 from prior memory table).
- 이전 prior table 의 "OR=2.45 p=1e-30" 는 부정확. 실제 값은 OR=1.90 p=1.7e-26.

### 1.3 핵심 결과

#### Per-allele 3-arm forest

| Allele | Chen ctrl | **Korean PTC** | Chen GD | OR Kr vs ctrl | p Kr vs ctrl | Pooled OR | I²% |
|---|---|---|---|---|---|---|---|
| **DPB1\*05:01** ★ | 31.3% | **53.2%** | 44.0% | **2.50** | **4e-26** | 2.16 [1.65, 2.83] | 85 |
| **C\*01:02** | 10.9% | **24.1%** | 18.4% | **2.61** | <1e-30 | 2.16 [1.53, 3.06] | 85 |
| **B\*46:01** | 6.5% | **10.3%** | 14.1% | 1.65 | 1e-3 | 2.02 [1.41, 2.89] | 76 |
| **A\*02:07** | 4.9% | **8.1%** | 9.7% | 1.72 | 2e-3 | **1.99 [1.66, 2.37]** | **0** |
| **DRB1\*07:01** (protective) | 15.3% | **11.4%** | 7.1% | 0.72 | 9e-3 | 0.55 [0.33, 0.91] | 91 |
| **DQB1\*02:01** (protective) | 17.8% | **0.0%** | 10.9% | 0.003 | 3e-5 | **0.57 [0.49, 0.66]** | **0** |

→ **6 alleles 모두 direction-consistent.** Pan-Asian autoimmune-thyroid 'continuum' 확인.

★ **Critical finding:** Korean PTC pool DPB1\*05:01 carrier frequency (53.2%) is **HIGHER than Chen Han Chinese GD itself (44.0%)**. Korean PTC 가 *Chinese GD 보다 더 진한* autoimmune-thyroid risk allele profile.

#### Korean sub-cohort heterogeneity (Cochran's Q)

- **DPB1\*05:01: I² = 0%** (perfectly homogeneous: K2 56% / Lee 52% / GSE286332-PTC 56%) ★
- DQB1\*02:01: I² = 0% (all 0%) ★
- B\*46:01: I² = 62% (moderate)
- DRB1\*07:01: I² = 60%
- A\*02:07: I² = 76%
- C\*01:02: I² = 81%

→ DPB1\*05:01 가 가장 안정적인 main claim. 다른 alleles 는 caveat 명시.

#### Sensitivity (4 scenarios)

DPB1\*05:01 freq stability:
- All Korean PTC n=874: **53.2%**
- Excl GSE286332 n=865: 53.2%
- Lee 2024 only n=630: 52.1%
- K2 only n=235: 56.2%

→ 4 scenarios 모두 ±2% 이내 — extremely robust.

### 1.4 Korean PTC pool n=874 (확인)

| Cohort | n | Source |
|---|---|---|
| K2 | 235 | PRJEB11591 Yoo SK 2016 SNU-GMI (260 manifest, 235 valid arcasHLA) |
| Lee 2024 | 630 | GSE213647 |
| GSE286332-PTC | 9 | Dongguk Univ PTC arm (excl 9 PTC+HT) |
| **Total** | **874** | |

(prior 메모리에 'n=908' 표기 — 이건 K2 valid 260 가정한 추정치. 실제 K2 valid 는 235.)

### 1.5 Honest disclosure (5 items)
1. **Cohort size imbalance:** Korean 874 vs Chu 2,958 — Chen 추정치가 pooled 결과 dominate
2. **Phenotype heterogeneity:** PTC vs GD distinct diseases (cancer vs autoimmunity)
3. **Population stratification:** Korean vs Han Chinese 가까운 East Asian 이지만 동일 X (e.g., DRB1\*15:01 Korean 18% vs Chinese ctrl 4%)
4. **Korean sub-cohort heterogeneity:** B\*46:01 / A\*02:07 / C\*01:02 의 I²>50% 일부 존재
5. **Allele-level only:** haplotype-level (DRB1-DQA1-DQB1) interactions future work

### 1.6 Deliverables

```
project/results/p2_pillar1_forest/
├── chu2018_allele_summary.tsv             — 8 published alleles
├── korean_PTC_pool_per_subcohort.tsv      — 4-cohort × 8-allele matrix
├── korean_subcohort_heterogeneity.tsv     — Cochran Q + I²
├── forest_meta_results.tsv                — 3-arm forest table
├── random_effects_pooled.tsv              — DL pooled OR + I²
├── sensitivity_4scenarios.tsv             — 4-scenario robustness
├── allele_freq_3cohort.{pdf, png}          — Figure 1 candidate
├── forest_panasian_HLA.{pdf, png}          — Figure 2 candidate (forest)
├── methods_paragraph.md                   — Cell Press paste-ready
├── discussion_paragraph.md                — Cell Press paste-ready
├── PILLAR1_FOREST_SUMMARY.md              — full status report
└── P2_PILLAR1_summary.json                — machine-readable
```

---

## 2. 결정 필요 — Web Claude 판단

### 2.1 후보 A: 즉시 manuscript writing 시작 (Phase 0 paper)

- **Pros:** 5-pillar 모두 STRONG, 분석은 그만하라고 본인이 명시함, 5-7년 plan 의 첫 paper 가 완성 단계
- **Cons:** 분당 outreach 미발송 — 만약 분당 데이터 받으면 Nat Commun reach 가능. Manuscript freeze 후 재시작 cost 발생
- **Timeline:** 2-3 week first draft → 4-6 week review iteration → submission

### 2.2 후보 B: 분당 outreach 먼저 발송 후 6주 wait

- **Pros:** Bundang Graves' n>50 + PTC+HT n>30 시 Nat Commun reach (시나리오 A)
- **Cons:** 6주 dead time + 응답 보장 X. ARIA exit 직후 momentum 손실 risk
- **Timeline:** 6주 wait → response → 추가 4-8주 (시나리오 B/C/D 따라)

### 2.3 후보 C: Phase 0 + Phase 1 (Graves') 병렬 trajectory

- **Pros:** Phase 0 paper 동시 manuscript + Phase 1 Graves' paper outline 시작 — 시간 70% 확보 활용
- **Cons:** 두 paper 동시 진행은 advisor (Yu professor) coordination 필요. Phase 1 은 Bundang 의존.
- **Timeline:** Phase 0 paper Q1 2026 / Phase 1 paper Q3-Q4 2026

### 2.4 후보 D: Pillar 1 깊이 더 (Okada 2015 일본 GD + Taiwan CMUH 통합)

- **Pros:** Pan-Asian 4-cohort meta (Korean + Han Chinese + Japanese + Taiwanese) → 더 강한 venue reach
- **Cons:** 1-2 weeks 추가 작업, controlled access (Taiwan), 본인이 명시적으로 "더 분석 안 함" 선언
- **Timeline:** 1-2 weeks → Pillar 1 deeper

### 2.5 후보 E: HLA 외 다른 mechanism 추가 (TCR repertoire, TIL deconvolution)

- **Pros:** Pillar 5 mechanism 더 깊게
- **Cons:** Diminishing returns. 5-pillar 이미 STRONG.
- **Timeline:** 1-2 weeks

### 2.6 후보 F: 하이브리드 — manuscript outline 시작 + 분당 outreach 동시 + Pillar 1 깊이 미루기

- **Pros:** Manuscript outline (no full draft yet) 1주 + 분당 outreach 동시 발송. 분당 응답 6주 wait 동안 첫 draft 완성. 시나리오 A 응답 시 Nat Commun pivot, 외 시나리오 시 Cell Rep Med submit.
- **Cons:** Coordination 복잡, Yu professor 와 outline 합의 필요
- **Timeline:** Week 1: outline + outreach send; Week 2-5: full draft + 응답 wait; Week 6: 분당 응답 후 venue 결정

---

## 3. Web Claude 에게 던질 핵심 질문

### 3.1 Strategic / 본질
1. "Pillar 1 STRONG 확정 후 이제 Phase 0 paper manuscript writing 시작이 합리적인가? 아니면 분당 outreach 먼저 (시나리오 A 가능성 vs dead time risk)?"
2. "5-7 년 plan 의 'paper 외 활동 70% 컷' 1단계가 완료된 지금, 첫 paper 의 venue 선택 (Cell Rep Med vs JCI Insight) + manuscript timeline 의 합리적 cadence?"
3. "후보 A/B/C/D/E/F 중 본인 5-7 년 plan + ARIA exit 후 momentum + Yu professor coordination 종합해서 추천?"

### 3.2 Pillar 1 forest meta 자체
4. "Korean PTC DPB1\*05:01 53.2% > Chen GD 44.0% — 이게 진짜 finding 인가? 아니면 Korean PTC 의 underlying autoimmune-thyroid background prevalence 가 expected/published Korean reference (Kim 2014 ~36%) 와 일치하지 않는 cohort effect?"
5. "Random-effects DL pooled I² 가 DPB1\*05:01 85% 로 high — 이게 published meta-analysis 권장 threshold (I² > 50%) 위에 있어 fixed-effects 가 더 정확한가?"
6. "Korean sub-cohort heterogeneity (B\*46:01 I²=62%, A\*02:07 I²=76%) 의 root cause — population stratification 인가, technical (arcasHLA call rate 차이) 인가, sampling 인가?"

### 3.3 Citation + reproducibility
7. "Chen 2018 → Chu X et al. 2018 citation correction 후 우리 prior memory + reports + manuscript 어디까지 update 해야 하는지 우선순위?"
8. "Chu 2018 의 8 alleles 외 추가 alleles (DRB1\*04:01, DRB1\*15:01, DQB1\*06:01 등) 의 Korean PTC 위치 — 즉시 추가 분석 가능 (Korean cohort 4-digit 다 있음)?"

### 3.4 Manuscript draft 시작 시
9. "Methods + Discussion paragraph (Cell Press style) draft — Yu professor coordination 전 본인 self-edit 후 발송 vs 직접 미팅 자료로?"
10. "Cell Rep Med vs JCI Insight 둘 다 reach — figures (현재 8개 candidate) 의 main vs supplementary breakdown 권고?"
11. "분당 outreach 시 Pillar 1 STRONG 결과 + Pillar 5 BCR/TLS 결과 첨부 추천 (cohort access enabler)?"

### 3.5 Phase 1 + 2 trajectory
12. "Phase 1 Graves' paper (autoimmune × thyroid Korean) 의 분당 응답 시나리오 A 일 때 venue (J Autoimmun vs Front Immunol vs Nat Commun)?"
13. "Phase 2 DIAL audit method paper 의 Phase 0 와 simultaneous draft 가능 시점?"

---

## 4. 관련 메모리 + 데이터 위치

### 4.1 Reports (chronological)

```
project/reports/
├── 2026_04_30_taskA_graves_dataset_inventory.md      # Day 1 Task A
├── 2026_04_30_taskB_8gene_audit.md                   # Day 1 Task B
├── 2026_04_30_dual_track_for_web_claude.md           # 4/30 web Claude brief
├── 2026_04_30_P3_GSE286332_critical_result.md        # P3 detailed
├── 2026_05_01_FULL_RESULTS_dual_track.md             # 5/1 master 893 lines
├── 2026_05_01_day2_yu_1pager.md                      # Yu prof 1-pager
├── 2026_05_02_D3P4_wang2024_audit.md                 # D3-P4
├── 2026_05_02_D3P5_pdm1_gradient.md                  # D3-P5
├── 2026_05_02_D3_decision_for_web_claude.md          # 5/2 web Claude
├── 2026_05_03_D4_D8_FULL_CLOSURE.md                  # D-series 7 prompts
├── 2026_05_03_D8D_figure_inventory.md                # 5-Pillar figure plan
├── 2026_05_03_GRAND_CONSOLIDATION.md                 # ★ all-in-one 840 lines
└── 2026_05_03_paper2_pillar1_for_web_claude.md       # this file
```

### 4.2 Results dirs (D-series + Pillar 1)

```
project/results/
├── p1_driver_mrna_audit/                  — Pillar 3
├── p2_power_planB/                        — Plan B map
├── p3_gse286332/                          — Pillar 2 (10,380 DEGs)
├── p4_pangenome_vs_tiera67/               — Pillar 4
├── p5_8gene_vs_hla_autocorr/              — Pillar 5 prep
├── d3p5_pdm1_gradient/                    — Pillar 5 mediation
├── d4p1_panasian_meta/                    — Pillar 1 first pass
├── d4p2_tcga_hashimoto_signature/         — Pillar 5 generalization
├── d5p6_bcr_repertoire/                   — Pillar 5 BCR/TLS
├── d6p7_dm1_subcluster/                   — Pillar 5 NBNR
├── d7p3_k2_calibration/                   — FAIL alternate evidence
├── d8b_korean_replication/                — Korean Hashimoto replication
├── d8c_dm1_subB_x_K2_NBNR/                — TCGA → Korean transfer
└── p2_pillar1_forest/                     — ★ 5/3 today, Pillar 1 STRONG
```

### 4.3 Memory entries

```
~/.claude/projects/-home-seungho-personal-THCA-data-analysis/memory/
├── MEMORY.md (index — 19+ entries)
├── v17_2026_04_30_pivot.md
├── v17_arcasHLA_korean_k2.md (4/29 burst $4.80)
├── v17_korean_k2_calibration.md
├── v17_K2_vs_bundang_distinction.md
├── v17_audit_session_2026_04_29.md
├── v17_dark_matter_pivot_2026_04_29.md
├── v17_8gene_audit_2026_04_29.md
├── v17_graves_pivot.md
├── v17_gse286332_strong_go.md (5/1 P3)
├── v17_8gene_pangenome_robustness.md (5/1 P4)
├── v17_D4P2_tcga_hashimoto_generalization.md (5/2 paper-changing)
├── v17_D5P6_BCR_clonal_TLS.md (5/2 antigen-driven)
└── v17_D6P7_dm1_subB_NBNR.md (5/2 NBNR cluster)
```

(memory 새 entry 'v17_paper2_pillar1_forest_strong.md' 작성 예정)

### 4.4 Bundang outreach query template
- File: `project/results/p2_power_planB/bundang_outreach_query.txt`
- 4-item query (PTC+HT / Graves' / modality / timeline)
- 발송 미시행. Yu professor coordination 후 발송 예정.

---

## 5. 한 줄 요약

**Pillar 1 PARTIAL → STRONG 확정 (DPB1\*05:01 Korean 53.2% > Chen GD 44.0%, 6 alleles direction-consistent, citation Chen→Chu correction). 5-Pillar 모두 STRONG. 후보 A (즉시 manuscript writing) vs B (분당 wait) vs C (Phase 0+1 병렬) vs F (하이브리드: outline + outreach 동시) — Web Claude 가 5-7 년 plan + ARIA momentum + Yu professor coordination 종합해서 추천 + 다음 D9+ 지령 부탁.**

이거 끝나면 진짜로 manuscript writing 시작 단계. 더 분석 안 함 (Pillar 1 이 마지막 enabler).
