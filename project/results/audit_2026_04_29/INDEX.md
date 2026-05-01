---
title: "rThyroid 마무리 분석 — 전체 결과 인덱스"
date: 2026-04-29
session: A → J 일괄 실행 (PM, Dark Matter reframe 후)
status: A/B/D/E/F/G/H/J completed; C delivered as coverage matrix + within-TCGA AUC; I plan-only
---

# 한 줄 요약

**8-gene paper는 차단되지 않고, 오히려 reframe 후 강화됨.** A audit는 의도된 design, B paradox는 group definition artifact (BRAF+TERT+ 25/36이 진짜 worst), D는 K2≠Bundang 정정, E는 MSK refractory enrichment 정량화 (chi² p=6.6e-131), F는 thyrocyte-intrinsic 확인, G는 monotonic Normal→PTC→PDTC→ATC 회복 (GSE213647), H는 FFPE 강건성 (KS p=0.44), C는 P8 가성비 (ΔAUC +0.007 vs P16, but cohort N 1518 vs 630), J는 Wang = Chen 2024 PMID 39235852.

# 결과 디렉토리

```
project/results/audit_2026_04_29/
├── INDEX.md                                    ← 이 파일
├── audit_report_8gene.md                       ← A: 포렌식 audit
├── meeting_brief_pm.md                         ← 미팅 1-pager (오전 작성)
├── 4way/
│   ├── 4way_audit.md                           ← B writeup
│   ├── 8cell_crosstab.tsv
│   ├── tert_subgroup_breakdown.tsv
│   ├── forest_HR_8cell.tsv
│   ├── summary.json
│   └── figure_4way_revalidation.pdf/png        ← B 4-panel figure
├── msk_bias/
│   ├── caveat_paragraph.md                     ← E paper-ready text
│   ├── histology_comparison.tsv
│   ├── msk_bias_summary.json
│   └── msk_bias_panel.pdf/png                  ← E 2-panel figure
├── ffpe_qc/
│   ├── ffpe_qc_summary.json                    ← H summary
│   └── figure_ffpe_qc.pdf/png                  ← H 3-panel figure
├── sc_wrapup/
│   ├── sc_wrapup_summary.json                  ← F summary
│   ├── celltype_DM_score.tsv
│   ├── histology_DM_score.tsv
│   ├── intra_tumor_heterogeneity.tsv
│   └── figure_sc_wrapup.pdf/png                ← F 4-panel figure
├── trajectory/
│   ├── trajectory_summary.tsv                  ← G summary
│   ├── trajectory_data.tsv
│   ├── trajectory_summary.json
│   └── figure_trajectory.pdf/png               ← G figure
├── robustness/
│   ├── gene_coverage_matrix.tsv                ← C coverage
│   ├── tcga_panel_auc.json
│   ├── robustness_summary.json
│   └── figure_panel_coverage.pdf/png           ← C 2-panel figure
└── nrg1_plan.md                                ← I deferred plan
```

# 결과별 핵심 (manuscript 반영용)

## A — 포렌식 audit (BLOCKER 해제)

- 8-gene은 RandomForest feature importance ranking 결과, candidate pool은 TIERA67 67-gene 큐레이션에서 `Driver_anchor` 카테고리 (BRAF, TERT, NRAS 등) 명시적 제외 (`rerun_v2.py:198`)
- Manuscript v6 title 자체가 "RAI-responsiveness biomarker" — RAI framing 이미 disclose
- Leak-free 재검증 (R1-B 30-gene MAPK+immune+EMT panel, 8-gene과 zero overlap)이 이미 존재. 8-gene은 새 cluster를 AUC 0.925로 예측 (ΔAUC +0.130 vs BRAF)
- **Action:** Methods 한 문장 reframe, "unsupervised genome-wide" → "RandomForest-ranked from curated 55-gene pool"

## B — TERT × BRAF × RAS 4-way (paradox 해소)

- 8-cell breakdown: TERT+ 36명 중 25명(69%)이 BRAF+, 4명(11%)만이 OTHER_TERT+ ("triple-neg-otherwise")
- "TERT+ triple-neg-otherwise" subgroup HR=6.90, but **CI=[0.009, 34.09]** ← 4 환자 small-N artifact
- 진짜 worst: BRAF+TERT+ (n=25, e=4, HR=3.04, p=0.04)
- **Action:** Figure 4 caption + supplementary 8-cell breakdown 표 추가

## D — K2 ≠ Bundang (cohort identity 정정)

- K2 = PRJEB11591 = Yoo 2016 SNU-GMI public (n=260)
- Bundang = outreach 단계 (2026-04-27 email v2), 데이터 없음
- **Action:** 미팅 + manuscript figure caption에서 정확히 호명

## E — MSK enrichment bias

- MSK 117명: 0% PTC, 28.2% ATC, 71.8% PDTC vs TCGA 94% PTC
- chi² p=6.6e-131, age MW p=8.5e-12 (median 61 vs 46)
- **Action:** Methods + Discussion에 caveat paragraph 추가 (file: caveat_paragraph.md)

## F — Single-cell wrap-up (thyrocyte-intrinsic 확인)

- Lu 2023 GSE193581 (n=67,678 cells, GSE184362 not located in project)
- DM_score는 Epithelial (n=706, median=1.36) + Malignant (n=14,624, median=-0.04)에서만 유의미
- T cell, B cell, Myeloid 등 immune/stromal 8-gene 발현 없음 → **8-gene은 thyrocyte-intrinsic**
- Inter-sample heterogeneity 측정 가능
- **Action:** Figure 5 panel, "thyrocyte-intrinsic differentiation signature" claim 강화

## G — PTC → PDTC → ATC trajectory

- GSE213647 (Korean Kim cohort): Normal=+0.50 → PTC=-0.51 → PDTC=-0.78 → ATC=-1.99 ← **monotonic**
- TCGA + K2: PTC only (PDTC/ATC 없음)
- GSE76039: PDTC + ATC만 (37 samples), small N
- **Action:** 4-cohort layered trajectory figure (figure_trajectory.pdf)

## H — FFPE robustness

- GSE213647 within-study FFPE (n=80, TruSeq) vs FF (n=169, same kit)
- panel_z KS p=0.44, MW p=0.75 → **차이 없음**
- 같은 kit 안에서 fixation effect는 8-gene signature를 흔들지 못함
- **Action:** Methods에서 "FFPE-compatible" claim 정량적 근거로 인용
- **주의:** "262 NGS FFPE" 중 GEO에 있는 건 80개. 나머지 182개는 unpublished. 미팅에서 분명히 짚을 것.

## C — Panel size robustness (P8 가성비)

- TCGA AUC for BRAF-like classification: P8=0.875 vs P16=0.882 (ΔAUC=+0.007)
- Cohort applicability: P8 N=1518 (3 RNA-seq cohorts), P10/12 N=630 (TCGA+MSK only, 둘다 mutation 필요)
- P10/12 추가하려면 MSK refractory population (E의 caveat 적용)
- **Action:** Discussion subsection "Why an 8-gene panel" — gene coverage figure 인용

## J — Wang 인용 (PMID 39235852)

- Chen XF / Wang YL 2024 Endocrine Connections, n=2,844 Shanghai NGS
- BRAF 71% / RAS 4% / TERT 3% — TCGA 일관
- **Action:** Discussion에 East-Asian comparator landscape로 인용

## I — NRG1 (defer)

- Korean germline 데이터 없음 → 별도 trajectory plan만 존재
- 분당/KoGES germline 도착 후 3-month side-paper로 진행
- 메인 paper에 영향 없음

# Manuscript 반영 우선순위

1. **즉시 (오늘 밤)**: A reframe 한 문장, B 8-cell figure, D K2 caption 정정
2. **revision round 1**: E caveat paragraph, F sc subsection, J Discussion citation
3. **revision round 2**: G trajectory figure (PDTC/ATC 데이터 더 모이면), C 가성비 sub
4. **별도 paper**: I NRG1 defer

# 재현

전체 분석:
```bash
cd /home/seungho/personal/THCA_data_analysis
python project/notebooks_or_scripts/v17_4way_revalidation.py    # B
python project/notebooks_or_scripts/v17_4way_figure.py          # B figure
python project/notebooks_or_scripts/v17_msk_bias_doc.py         # E
python project/notebooks_or_scripts/v17_h_ffpe_qc.py            # H
python project/notebooks_or_scripts/v17_f_sc_wrapup.py          # F
python project/notebooks_or_scripts/v17_g_trajectory.py         # G
python project/notebooks_or_scripts/v17_c_robustness.py         # C
```
