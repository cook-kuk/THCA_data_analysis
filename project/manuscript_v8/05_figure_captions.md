---
title: "Paper 1 manuscript v8 — Figure 1-8 captions + Suppl S1-S9 (draft v1)"
date: 2026-05-01
status: v1 draft — Cell Press style. 본인 panel 그림 작성 완료 후 caption 정확 표현 검증.
target_format: bold title + (A)/(B)/.. panel descriptions + statistical notes (n, test, p)
---

# Main Figures (8 main, Cell Press 4-6 panels each)

## Figure 1. The 8-gene RAI-responsiveness panel resolves a DM1/DM2 cluster within papillary thyroid carcinoma. (4 panels v3)

(A) Sankey flow diagram showing the curation path from the TIERA67 67-gene candidate pool (seven thyroid-relevant biological categories) to the final 8-gene panel (TDS-core sub-category: SLC5A5, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, DIO1). (B) UMAP embedding of TCGA-THCA primary tumors (n = 504) on the eight-dimensional 8-gene transcript space, colored by KMeans-derived DM1 (red) versus DM2 (blue) cluster assignment. (C) Driver mRNA neutrality. Box plots of BRAF transcript expression in V600E carriers (n = 273) versus wild-type tumors (n = 182), Cohen's d = −0.044, Mann-Whitney p = 0.57. Single-feature classification AUCs for DM1 versus DM2 are shown for BRAF (0.602), TERT (0.578), KRAS (0.525), NRAS (0.521), and HRAS (0.500). (D) Pan-genome cluster reproducibility. Adjusted Rand index (ARI) of unsupervised KMeans partitions versus the 8-gene-driven DM1/DM2 reference: 8-gene panel alone (0.49), TIERA67 (0.90), pan-genome top-5000 by median absolute deviation (0.92), and Driver_anchor genes alone (−0.007). Hypergeometric enrichment p of TIERA67 within pan-genome top 100 = 3 × 10⁻⁴.

Statistical tests: KMeans k = 2; ARI computed against DM1/DM2 reference; Cohen's d with pooled standard deviation; Mann-Whitney U for transcript distributions.

---

## Figure 2. Clinical aggressiveness within Xing 2014 dark matter. (3 panels)

(A) Sankey diagram showing the recovery of 131 of 180 (73%) Xing 2014 BRAF/TERT-negative dark-matter tumors into a defined DM1/DM2 stratum via the 8-gene panel. (B) DM1 prevalence by ancestry: TCGA-THCA (predominantly European/North American; 28.4% BRAF/RAS-negative dark matter; n = 504) versus Korean cohort (37.8% dark matter; n = 874 K2/Lee/GSE286332 pool). (C) Kaplan-Meier overall-survival curves for DM1 (red) versus DM2 (blue) within Xing dark matter (TCGA-THCA, n = 180). Log-rank test reported.

Statistical test: log-rank for KM curves; Wilson 95% CI for proportion estimates.

---

## Figure 3. Single-cell external validation confirms thyrocyte-intrinsic DM1 signal. (3 panels)

(A) Per-patient Spearman correlation heatmap of DM1 score in patient-matched tumor versus normal thyrocytes from GSE184362 (Pu et al., 2021; n = 6 PTC patients, Fudan University), with all 6 patients showing r between 0.798 and 0.886 (p < 10⁻¹⁰). (B) UMAP embedding of thyrocyte clusters from Lu 2023 (GSE193581, n = 23 samples), colored by DM1 score with thyrocyte-specific signal isolated from stromal and immune compartments. (C) Author-independence check: cross-cohort DM1 score concordance between GSE241184 (single-author Phase 1) and GSE184362 (Pu et al., 2021) confirms that the DM1 signal is not driven by laboratory-specific batch effects.

Statistical tests: Spearman r with Bonferroni-adjusted p; UMAP via PCA→neighborhood graph.

---

## Figure 4. Mutation, TERT promoter, and outcome stratification. (3 panels)

(A) Stacked bar plot of 8-cell decomposition (BRAF mutated × TERT mutated × DM cluster) across TCGA-THCA tumors, showing the distribution of patients into the eight defined molecular strata. (B) Kaplan-Meier overall-survival curves for BRAF V600E + TERT promoter mutated (BRAF_TERT+) tumors versus the BRAF/TERT-negative dark matter stratum stratified by DM cluster. (C) The 4-patient OTHER_TERT+ caveat box: TCGA contains only 4 BRAF-wildtype + TERT promoter-mutated cases, limiting any small-N inference about this stratum.

Statistical tests: log-rank for KM; Cox proportional hazards.

---

## Figure 5. Korean cohort validation and FFPE compatibility. (3 panels)

(A) HLA-II module Cohen's d between DM1 and DM2 in Korean GSE213647 (Lee et al., n = 632), with replication of d ≈ 0.95 — comparable to the TCGA-THCA estimate. (B) DM cluster composition in K2 (PRJEB11591, n = 260) NBNR (BRAF-negative + RAS-negative) cohort showing mixed phenotype (vascular invasion enrichment within DM1). (C) Formalin-fixed paraffin-embedded (FFPE) versus fresh-frozen (FF) tissue concordance: density plot of 8-gene scores by processing type, Kolmogorov-Smirnov p = 0.44 (no detectable distributional shift).

Statistical tests: Cohen's d; Kolmogorov-Smirnov for distributional comparison.

---

## Figure 6. Pooled meta-analysis of DM1 overall-survival hazard. (3 panels)

(A) Forest plot of Cox proportional-hazards estimates for DM1 versus DM2 overall-survival risk, by cohort: TCGA-THCA (HR 2.30, 95% CI 0.77–6.88; n = 504, 16 events), MSK-IMPACT thyroid (HR 2.67, 95% CI 1.17–6.10; n = 117, 38 events). Pooled DerSimonian-Laird random-effects estimate: HR 2.53 (95% CI 1.31–4.89), Cochran I² = 0%. (B) I² heterogeneity panel showing no detectable between-cohort heterogeneity (I² = 0%). (C) Sub-cluster Cox analysis for DM1 sub-A versus sub-B (TCGA n = 91, 4 events): DM1 sub-B HR 0.31 (95% CI 0.07–1.39, p = 0.13) — directionally consistent with the immune-overlap phenotype but underpowered as a standalone claim.

Statistical tests: Cox proportional hazards (R `survival`, Python `lifelines`); DerSimonian-Laird random-effects meta-analysis.

---

## Figure 7. DM1 is a fusion-driven dark matter subtype with mechanistic heterogeneity. (5 panels) ★

(A) Silhouette plot of DM1 sub-A (n = 72) versus sub-B (n = 19) sub-clusters within TCGA-THCA, with cluster silhouette score 0.584 supporting the two-population structure. (B) Stacked bar of fusion-positivity by DM cluster: DM1 76.8% (63/82), DM2 30.9%, Fisher OR 7.41 (95% CI 4.38–12.55, p = 1.9 × 10⁻¹³). (C) Stacked bar of fusion partners within DM1 fusion-positive cases: RET (n = 33; CCDC6-RET 17, NCOA4-RET 3, other 13), NTRK (n = 10), ALK (n = 4), BRAF fusions (n = 5). (D) DM1 sub-A versus sub-B phenotype panel: age (37.3 vs 51.3 years; Cohen's d = −0.82, MW p = 0.004), stage III/IV (15.3% vs 44.4%; OR 0.23, p = 0.020), CD8/IFN-γ/checkpoint signatures (Cohen's d = −0.5 to −0.6 vs sub-B). (E) DM1 capture rate of TCGA RET-fusion-positive tumors: 27 of 33 (81.8%), supporting a clinical reflex testing algorithm.

Statistical tests: KMeans silhouette; Fisher exact for OR; Cohen's d; Mann-Whitney U; sensitivity analyses for SV missingness (chi² p = 0.56 MAR).

---

## Figure 8. DM1 epigenetically silences thyroid differentiation machinery. (3 panels v3) ★

(A) Per-gene HM450 promoter β-value heatmap for the 8-gene panel + DIO2 + SLC26A4, comparing DM1, DM2, and not_DM tumors (n = 503 with HM450 + DM call). DM1 vs DM2 per-gene Cohen's d: TPO (2.30, p = 1.9 × 10⁻¹⁸), DIO1 (1.24, p = 6.5 × 10⁻¹¹), TSHR (1.20, p = 9.8 × 10⁻¹²), PAX8 (0.97, p = 4.5 × 10⁻⁸), TG (0.86, p = 2.2 × 10⁻⁶), FOXE1 (0.84, p = 1.0 × 10⁻⁵), NKX2-1 (0.63, p = 8.9 × 10⁻⁷), SLC5A5 (0.22, p = 0.42, NS). (B) Mean 8-gene panel β-value bar plot by DM cluster: DM1 = 0.385, DM2 = 0.253, not_DM = 0.356 — corresponding to 52% higher DM1 promoter methylation versus DM2. (C) Methylation fusion-independence: within DM1, fusion-positive (n = 63) versus fusion-negative (n = 19) tumors show equivalent mean panel β (Cohen's d = −0.36, MW p = 0.31, NS), supporting epigenetic silencing as a mechanism layer parallel to fusion drivers.

Statistical tests: Cohen's d (pooled SD); Mann-Whitney U for per-gene β; Fisher exact for fusion × DM cross-tab.

---

# Supplementary Figures (S1-S9)

**S1.** 8-gene panel heatmap sorted by P_DM1 score, comparing TCGA versus K2 cohort distributions (Fig 1 C 이동, v3).

**S2.** 8-gene panel similarity matrix from external single-cell datasets (Pu 2021, Lu 2023, GSE241184).

**S3.** Pan-genome cluster ARI ladder spanning panel sizes 8, 16, 67 (TIERA67), 200, 1000, and 5000 by median absolute deviation.

**S4.** HLA-II residualization analysis. DM1 vs DM2 Cohen's d for 8-gene panel score before residualization (raw d = 1.78), after residualization on HLA-II module (d = 1.00), after residualization on generic immune signature (d = 1.50), and after dual residualization on HLA-II + immune (d = 0.87).

**S5.** arcasHLA Korean Pan-Asian three-arm forest plot. Korean PTC pool (n = 874) vs GSE286332 PTC+HT (n = 9) vs Chu et al. 2018 Han Chinese Graves' disease cohort. (Pillar 1, primary role in Paper 2.)

**S6.** B cell receptor (BCR) clonal architecture and tertiary lymphoid structure (TLS) heatmap in GSE286332 PTC+HT versus PTC. (Paper 2 main; Paper 1 supp only.)

**S7.** DM1 sub-B × Korean K2 NBNR signature transfer in GSE213647 (Lee et al., n = 632), with sub-B-like rate 47.2% (GMM) to 52.5% (Otsu). (Paper 2 main.)

**S8.** Hypomethylating-agent + radioiodine re-induction schematic + literature meta. Decitabine + I-131 retrospective trial summary (NCT00085293, NCT01065090) and the rationale for prospective trial design stratified by DM1 status. (Fig 8 D 이동, v3.)

**S9.** K2 mini-index calibration FAIL diagnostic + alternate evidence chain. Per-gene inflation factors (4.9–12.5× across panel genes); 4-metric direction check (raw, within-sample-z, per-gene-z, per-gene-rank); alternate evidence (DM call distribution 94.6% DM2, HLA arm n = 874 Korean PTC pool).

---

# Figure 작성 우선순위 (Cell Press)

| Priority | Figure | Status | Source code |
|---|---|---|---|
| ★★★ | Fig 7 (DM1 mechanism, 5 panels) | Build needed — combines R3-F4 + R4-2 + R4-3 | `v17_audit_F2_F3_F4.py` + `v17_audit_R4_all.py` + new |
| ★★★ | Fig 8 (Epigenetic, 3 panels v3) | Build needed — R5-2 paradigm | `v17_audit_R5_all.py` + new |
| ★★★ | Fig 6 (Meta forest, 3 panels) | Build needed — N1 meta + R4-2 sub-cluster Cox | `v17_D4P1_forest_meta.py` + new |
| ★★ | Fig 1 (Panel + cluster, 4 panels v3) | Build needed | `v17_8gene_figs_v2.py` + `p4_pangenome_vs_tiera67.py` |
| ★★ | Fig 2 (Xing rescue + KM) | Build needed | `v17_4way_figure.py` + new |
| ★ | Fig 3 (sc validation) | Existing figs available | sc figs + new author-independence panel |
| ★ | Fig 4 (BRAF×TERT 8-cell) | Existing figs available | `v17_quad_tert_stack.html` + new |
| ★ | Fig 5 (Korean + FFPE) | Existing figs available | `v17_KOREAN_K2_v260_figure.py` + new |

## 본인 voice 영역

- [ ] 각 figure caption 첫 줄 (figure title) — 본인 voice 검증
- [ ] Statistical notes 표기 일관성 — Mann-Whitney U vs MW vs ranksum 결정
- [ ] Sub-figure cross-reference (Methods, Suppl Table) — 본인 검증

## 다음 step

1. Figure code build (W3-W4) — Plotly/matplotlib script per figure
2. 본인 figure preview → caption verbatim 정정
3. Cover letter Figure 7+8 game-changer claim 정합 (W5)
