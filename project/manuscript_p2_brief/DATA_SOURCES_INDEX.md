---
title: "Paper 2 advisor brief — DATA SOURCES INDEX"
date: 2026-05-01
purpose: "p2_advisor_discussion.html 의 21 figures + 5 tables 를 underlying TSV/JSON 파일로 1:1 mapping. Advisor 가 source 검증 가능."
verified: "2026-05-01 — 54 hard-coded data points (Figures 2/2b/2c/2d) all match source TSV (Pillar I)."
---

# Paper 2 advisor brief — DATA SOURCES INDEX

`project/manuscript_p2_brief/p2_advisor_discussion.html` 의 모든 numeric claim → underlying file mapping.

## A. Figures (21개, 6 paper-defining ★★★)

### A.1 Pillar I — Korean Pan-Asian HLA forest (NEW 2026-05-03)

| Fig | Title | Source file | Verified |
|---|---|---|---|
| Fig 1 | arcasHLA top 5 carrier freq | `results/p2_pillar1_forest/forest_meta_results.tsv` | ✅ |
| **Fig 2 ★★★** | 6-allele 3-arm forest (Korean PTC vs Chu ctrl vs Chu GD) | `results/p2_pillar1_forest/forest_meta_results.tsv` + `chu2018_allele_summary.tsv` | ✅ 18/18 cells |
| Fig 2b | Random-effects pooled OR + I² | `results/p2_pillar1_forest/random_effects_pooled.tsv` | ✅ 6/6 |
| Fig 2c | Korean sub-cohort heterogeneity (K2 / Lee / GSE286332) | `results/p2_pillar1_forest/korean_subcohort_heterogeneity.tsv` | ✅ 18/18 cells |
| Fig 2d | 4-scenario sensitivity DPB1*05:01 | `results/p2_pillar1_forest/sensitivity_4scenarios.tsv` | ✅ |

### A.2 Pillar II — GSE286332 PTC vs PTC+HT in-cohort

| Fig | Title | Source file |
|---|---|---|
| Fig 3 | 8-gene RAI Cohen's d per-gene | `results/p3_gse286332/P3_summary.json` (`panel_8gene_compare`) + `8gene_per_gene_compare.tsv` |
| Fig 4 | GSEA top pathways NES | `results/p3_gse286332/gsea_MSigDB_Hallmark_2020.tsv` + `gsea_KEGG_2021_Human.tsv` + `gsea_Reactome_2022.tsv` |
| Fig 5 | Top 10 up DEGs | `results/p3_gse286332/P3_summary.json` (`deg_top10_up`) + `deg_ptcht_vs_ptc.tsv` |
| **Fig 6 ★★★** | HLA-I / HLA-II module Cohen's d (+3.65 anchor) | `results/p3_gse286332/hla_module_scores.tsv` |
| Fig 7 | 2-cohort meta forest 8-gene RAI | `results/p5_8gene_vs_hla_autocorr/meta_3cohort.tsv` |

### A.3 Pillar III — Cross-cohort generalization

| Fig | Title | Source file |
|---|---|---|
| Fig 8 | Hashimoto-like prevalence × 4 thresholds × 2 cohorts | `results/d4p2_tcga_hashimoto_signature/D4P2_summary.json` + `results/d8b_korean_replication/D8B_summary.json` |
| **Fig 9 ★★★** | TCGA DM1/DM2 × Hashimoto-like cross-tab (OR=0.20, p=6.4e−10) | `results/d4p2_tcga_hashimoto_signature/D4P2_summary.json` (`crosstabs.hashi_top30`) |
| Fig 10 | Stromal+immune residualized OR (robustness) | `results/d4p2_tcga_hashimoto_signature/D4P2_summary.json` (`crosstabs.hashi_resid_otsu`) |
| Fig 11 | HLA-II saturation analysis (DM1/DM2 × Hashimoto subset) | `results/d4p2_tcga_hashimoto_signature/D4P2_summary.json` (HLA-II residualization) |

### A.4 Pillar IV — Antigen-driven specificity

| Fig | Title | Source file |
|---|---|---|
| **Fig 12 ★★★** | TLS / IGHV / IGKV clonality + entropy | `results/d5p6_bcr_repertoire/group_comparison.tsv` |
| Fig 13 | Spearman correlation matrix (8-gene × IGHV × TLS × HLA-II) | `results/d5p6_bcr_repertoire/spearman_matrix.tsv` |

### A.5 Pillar V — DM1 sub-B = NBNR cluster

| Fig | Title | Source file |
|---|---|---|
| **Fig 14 ★★★** | DM1 sub-A vs sub-B mutation landscape | `results/d6p7_dm1_subcluster/dm1_subcluster_labels.tsv` |
| Fig 15 | Korean GSE213647 sub-B-like rate | `results/d8c_korean_subB/korean_subB_score.tsv` |

### A.6 § 3.6 Mediation

| Fig | Title | Source file |
|---|---|---|
| **Fig 16 ★★★** | Mediation % bootstrap (HLA-II 140%) | `results/d3p5_mediation/mediation_results.json` |
| Fig 17 | Within-PTC severity gradient (g8 vs P_DM1) | `results/d3p5_mediation/within_ptc_severity.json` |

### A.7 § 6 Timeline

| Fig | Title | Source |
|---|---|---|
| Fig 18 | 6-week parallel Gantt — Paper 1 writing × Paper 2 gap analysis | 본 brief § 6.0 plan (no external file) |

## B. Tables (5개, NEW 2026-05-03 격상)

| Table | Title | Source file |
|---|---|---|
| Table 1 | Mediation summary bootstrap 5000 iter | `results/d3p5_mediation/mediation_results.json` |
| Table 1.0 (NEW) | Korean sub-cohort heterogeneity | `results/p2_pillar1_forest/korean_subcohort_heterogeneity.tsv` |
| Table 1.1 (NEW) | 4-scenario sensitivity 6 alleles | `results/p2_pillar1_forest/sensitivity_4scenarios.tsv` |
| Table 1.2 (NEW) | Pillar I 5-row 한계 disclosure | (수동 정리 — pre-emptive reviewer 답변용) |
| Table 2 | Paper 2 writing schedule | 본 brief § 6.1 |
| Table 2.0 (NEW) | P5 8-gene × HLA-II residualization (TCGA + GSE286332) | `results/p5_8gene_vs_hla_autocorr/P5_summary.json` |

## C. Paste-ready paragraphs (Cell Press style)

### C.1 Pillar 1 forest meta — Methods + Discussion

| Paragraph | Source file |
|---|---|
| § 3.1.6 Methods (Pillar 1 forest meta) | `results/p2_pillar1_forest/methods_paragraph.md` |
| § 3.1.7 Discussion (Pillar 1 forest meta) | `results/p2_pillar1_forest/discussion_paragraph.md` |

## D. Cohort raw data

| Cohort | Description | n | Source |
|---|---|---|---|
| TCGA-THCA | Discovery cohort | 500 | TCGA pan-cancer (`results/v17_tcga_*/` 각종) |
| GSE286332 | Korean Dongguk Univ Lim 2025 | 18 (9 PTC + 9 PTC+HT) | `data/GSE286332/` (BioProject PRJNA1208932) |
| GSE213647 | Lee 2024 Korean replication | 632 | `data/GSE213647/` |
| K2 / PRJEB11591 | Yoo 2016 SNU-GMI | 235 typeable (260 manifest) | `data/k2_korean/` |
| Chu 2018 reference | Han Chinese GD vs ctrl | 1,468 GD / 1,490 ctrl | `results/p2_pillar1_forest/chu2018_allele_summary.tsv` (published values verified PMC 6161647) |

## E. Verification methodology

본 index 의 ✅ 표시 항목은 다음 sanity check 으로 검증:

```python
# project/manuscript_p2_brief/p2_advisor_discussion.html
# 의 Plotly 'const KOREAN_F2', 'const POOL_OR' 등 hard-coded array
# 를 source TSV 의 해당 column 과 0.5% tolerance 으로 비교

forest_meta_results.tsv     → Fig 2 KOREAN_F2 / CHU_CTRL_F2 / CHU_GD_F2  → 18/18 cells ✅
random_effects_pooled.tsv   → Fig 2b POOL_OR                              → 6/6 ✅
korean_subcohort_heterogeneity.tsv → Fig 2c SUB_K2 / SUB_LEE / SUB_286   → 18/18 cells ✅
```

검증 script: `project/manuscript_p2_brief/PROMPT_DECISION_LOG.md` § 4.3 참조.

## F. 2026-05-01 통합 변경 추적 (v4 → brief)

| 변경 | 위치 | 출처 |
|---|---|---|
| Pillar I PARTIAL → STRONG | Hero, exec summary, pillar grid, § 3.1 전체 | `results/p2_pillar1_forest/` 신규 자료 |
| Citation Chen → Chu 2018 | 모든 occurrence (Fig 2, 2b, 2c, 2d, Methods, Discussion) | Web fetch PMC 6161647 verified 2026-05-03 |
| 8-gene × HLA-II autocorr (Table 2.0) | § 3.6.1 | `results/p5_8gene_vs_hla_autocorr/P5_summary.json` |
| Q9–Q12 advisor questions | § 7 | 신규 (citation propagation, venue ladder, DQB1 typing, haplotype) |
| Gap 3 RESOLVED | § 4 | 동일 source — 작업 완료 |
| Scenario B 확률 50% → 55% | § 5 | Pillar I STRONG 격상 reflect |

## G. 다른 산출물에서 source 사용 시

| 사용 케이스 | 권장 path |
|---|---|
| Paper 2 manuscript (voice-protected, user 직접 작성) | `results/p2_pillar1_forest/methods_paragraph.md` + `discussion_paragraph.md` 통째로 paste 가능 |
| Cover letter Pillar I 한 줄 sentence | "DPB1*05:01 53.2% in n=874 Korean PTC vs 31.3% in Chu Han Chinese ctrl, OR 2.50 [2.10, 2.97], pooled OR 2.16 [1.65, 2.83]." |
| Reviewer Q1 (Pillar I quantitative claim) | § 3.1 전체 + § 3.1.8 Honest disclosure 5-row |
| Supplementary table S5 (arcasHLA forest) | `forest_meta_results.tsv` rename → `Table_S5_panasian_forest.tsv` |
| Supplementary table S6 (sub-cohort heterogeneity) | `korean_subcohort_heterogeneity.tsv` |
| Supplementary table S7 (4-scenario sensitivity) | `sensitivity_4scenarios.tsv` |
| Supplementary table S9 (per-sample 4-digit alleles) | `results/d4p1_panasian_meta/korean_PTC_pool_n908.tsv` (NOT deprecated, KEEP) |

## H. ⚠️ Deprecated / superseded 자료 사용 금지

| 사용 금지 path | 대신 사용 |
|---|---|
| `results/d4p1_panasian_meta/D4P1_summary.json` | `results/p2_pillar1_forest/P2_PILLAR1_summary.json` |
| `results/d4p1_panasian_meta/panasian_forest_meta.tsv` | `results/p2_pillar1_forest/forest_meta_results.tsv` |
| `results/d4p1_panasian_meta/chen2018_han_chinese_GD_published.tsv` | `results/p2_pillar1_forest/chu2018_allele_summary.tsv` |
| `notebooks_or_scripts/v17_D4P1_forest_meta.py` | `notebooks_or_scripts/v17_paper2_pillar1_forest.py` |

→ DEPRECATED 사유 + 자세한 분류: `results/d4p1_panasian_meta/DEPRECATED.md`

---

*Generated 2026-05-01 by audit of `project/manuscript_p2_brief/p2_advisor_discussion.html`. 54 hard-coded numeric assertions (Pillar I figures) verified against source TSV. Other 17 figures rely on source files listed above.*
