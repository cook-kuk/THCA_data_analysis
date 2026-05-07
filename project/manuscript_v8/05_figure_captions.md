---
title: "Paper 1 manuscript v8 — Figure 1-8 captions + Suppl S1-S9"
date: 2026-05-04
status: clean draft
target_format: bold title + (A)/(B)/.. panel descriptions + statistical notes (n, test, p)
---

# Main Figures (8 main, Cell Press 4-6 panels each)

## Figure 1. The 8-gene RAI-responsiveness panel resolves a DM1/DM2 cluster within papillary thyroid carcinoma. (4 panels v3)

(A) Sankey flow diagram showing the curation path from the TIERA67 67-gene candidate pool (seven thyroid-relevant biological categories) to the final 8-gene panel (TDS-core sub-category: SLC5A5, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, DIO1). (B) UMAP embedding of TCGA-THCA primary tumors (n = 504) on the eight-dimensional 8-gene transcript space, colored by KMeans-derived DM1 (red) versus DM2 (blue) cluster assignment. (C) Driver mRNA neutrality. Box plots of BRAF transcript expression in V600E carriers (n = 273) versus wild-type tumors (n = 182), Cohen's d = −0.044, Mann-Whitney p = 0.57. Single-feature classification AUCs for DM1 versus DM2 are shown for BRAF (0.602), TERT (0.578), KRAS (0.525), NRAS (0.521), and HRAS (0.500). (D) Pan-genome cluster reproducibility. Adjusted Rand index (ARI) of unsupervised KMeans partitions versus the 8-gene-driven DM1/DM2 reference: 8-gene panel alone (0.49), TIERA67 (0.90), pan-genome top-5000 by median absolute deviation (0.92), and Driver_anchor genes alone (−0.007). Hypergeometric enrichment p of TIERA67 within pan-genome top 100 = 3 × 10⁻⁴.

Statistical tests: KMeans k = 2; ARI computed against DM1/DM2 reference; Cohen's d with pooled standard deviation; Mann-Whitney U for transcript distributions.

---

## Figure 2. Clinical aggressiveness within Xing 2014 dark matter. (3 panels)

(A) Sankey diagram showing the recovery of 131 of 180 (73%) Xing 2014 BRAF/TERT-negative dark-matter tumors into a defined DM1/DM2 stratum via the 8-gene panel. (B) DM1 prevalence by ancestry: TCGA-THCA (predominantly European/North American; 28.4% BRAF/RAS-negative dark matter; n = 504) versus Korean cohort (37.8% dark matter; n = 865 K2 + Lee pool — GSE286332-PTC dropped to preserve Paper 2 boundary). (C) Kaplan-Meier overall-survival curves for DM1 (red) versus DM2 (blue) within Xing dark matter (TCGA-THCA, n = 180). Log-rank test reported.

Statistical test: log-rank for KM curves; Wilson 95% CI for proportion estimates.

---

## Figure 3. Single-cell external validation confirms thyrocyte-intrinsic DM1 signal. (3 panels)

(A) Per-patient Spearman correlation heatmap of DM1 score in patient-matched tumor versus normal thyrocytes from GSE184362 (Pu et al., 2021), with per-patient r values ranging from 0.798 to 0.886 (p < 10⁻¹⁰). (B) UMAP embedding of thyrocyte clusters from Lu 2023 (GSE193581, n = 23 samples), colored by DM1 score with thyrocyte-specific signal isolated from stromal and immune compartments. (C) Author-independence check: cross-cohort DM1 score concordance between GSE241184 (single-author Phase 1) and GSE184362 (Pu et al., 2021) confirms that the DM1 signal is not driven by laboratory-specific batch effects.

Statistical tests: Spearman r with Bonferroni-adjusted p; UMAP via PCA→neighborhood graph.

---

## Figure 4. Mutation, TERT promoter, and outcome stratification. (3 panels)

(A) Stacked bar plot of 8-cell decomposition (BRAF mutated × TERT mutated × DM cluster) across TCGA-THCA tumors, showing the distribution of patients into the eight defined molecular strata. (B) Kaplan-Meier overall-survival curves for BRAF V600E + TERT promoter mutated (BRAF_TERT+) tumors versus the BRAF/TERT-negative dark matter stratum stratified by DM cluster. (C) The 4-patient OTHER_TERT+ caveat box: TCGA contains only 4 BRAF-wildtype + TERT promoter-mutated cases, limiting any small-N inference about this stratum.

Statistical tests: log-rank for KM; Cox proportional hazards.

---

## Figure 5. Korean cohort validation and FFPE compatibility. (3 panels)

(A) Independent Korean cohort validation. Distribution of 8-gene DM scores in GSE213647 (Lee et al., n = 632), showing preservation of the DM1/DM2 score boundary in an external East-Asian cohort. (B) DM cluster composition in K2 (PRJEB11591, n = 260) NBNR (BRAF-negative + RAS-negative) cohort showing mixed phenotype (vascular invasion enrichment within DM1). (C) Formalin-fixed paraffin-embedded (FFPE) versus fresh-frozen (FF) tissue concordance: density plot of 8-gene scores by processing type, Kolmogorov-Smirnov p = 0.44 (no detectable distributional shift).

Statistical tests: Cohen's d; Kolmogorov-Smirnov for distributional comparison.

---

## Figure 6. Pooled meta-analysis of DM1 overall-survival hazard. (3 panels)

(A) Forest plot of Cox proportional-hazards estimates for DM1 versus DM2 overall-survival risk, by cohort: TCGA-THCA (HR 2.30, 95% CI 0.77–6.88; n = 504, 16 events), MSK-IMPACT thyroid (HR 2.67, 95% CI 1.17–6.10; n = 117, 38 events). Pooled DerSimonian-Laird random-effects estimate: HR 2.53 (95% CI 1.31–4.89), Cochran I² = 0%. (B) I² heterogeneity panel showing no detectable between-cohort heterogeneity (I² = 0%). (C) Sub-cluster Cox analysis for DM1 sub-A versus sub-B (TCGA n = 91, 4 events): DM1 sub-B HR 0.31 (95% CI 0.07–1.39, p = 0.13) — directionally consistent with the immune-overlap phenotype but underpowered as a standalone claim.

Statistical tests: Cox proportional hazards (R `survival`, Python `lifelines`); DerSimonian-Laird random-effects meta-analysis.

---

## Figure 7. DM1 is a fusion-driven dark matter subtype. (4 panels — v3 → v4 2026-05-08, panel D demoted to S6 per audit P1-7)

(A) Silhouette plot of DM1 sub-A (n = 72) versus sub-B (n = 19) sub-clusters within TCGA-THCA, with cluster silhouette score 0.584 supporting the two-population structure. (B) Stacked bar of fusion-positivity by DM cluster: DM1 76.8% (63/82), DM2 30.9%, Fisher OR 7.41 (95% CI 4.38–12.55, p = 1.9 × 10⁻¹³). (C) Stacked bar of fusion partners within DM1 fusion-positive cases: RET (n = 33; CCDC6-RET 17, NCOA4-RET 3, other 13), NTRK (n = 10), ALK (n = 4), BRAF fusions (n = 5). (D) DM1 capture rate of TCGA RET-fusion-positive tumors: 27 of 33 (81.8%), supporting a clinical reflex testing algorithm.

<em>Panel demotion (v4, 2026-05-08): Original panel D (DM1 sub-A vs sub-B phenotype: age, stage, CD8/IFN-γ/checkpoint) moved to Supplementary Figure S6 to keep main Figure 7 mechanism-focused (4 panels) and reserve sub-B mechanistic claims for the companion Paper 2; Cell Press main-figure count thereby reduces from 8 to 7. Original panel E becomes new panel D.</em>

Statistical tests: KMeans silhouette; Fisher exact for OR; Cohen's d; Mann-Whitney U; sensitivity analyses for SV missingness (chi² p = 0.56 MAR).

---

## Figure 8. DM1 epigenetically silences thyroid differentiation machinery. (2 panels — v3 → v4 2026-05-08, panel C demoted to S5b per audit P1-5)

(A) Per-gene HM450 promoter β-value heatmap for the 8-gene panel + DIO2 + SLC26A4, comparing DM1, DM2, and not_DM tumors (n = 503 with HM450 + DM call). DM1 vs DM2 per-gene Cohen's d: TPO (2.30, p = 1.9 × 10⁻¹⁸), DIO1 (1.24, p = 6.5 × 10⁻¹¹), TSHR (1.20, p = 9.8 × 10⁻¹²), PAX8 (0.97, p = 4.5 × 10⁻⁸), TG (0.86, p = 2.2 × 10⁻⁶), FOXE1 (0.84, p = 1.0 × 10⁻⁵), NKX2-1 (0.63, p = 8.9 × 10⁻⁷), SLC5A5 (0.22, p = 0.42, NS). (B) Mean 8-gene panel β-value bar plot by DM cluster: DM1 = 0.385, DM2 = 0.253, not_DM = 0.356 — corresponding to 52% higher DM1 promoter methylation versus DM2.

<em>Panel demotion (v4, 2026-05-08): Original panel C (within-DM1 fusion+ vs fusion− methylation, n = 63 vs 19, Cohen's d = −0.36, NS p = 0.31) moved to Supplementary Figure S5b. Rationale per audit P1-5: a non-significant Cohen's d on n = 19 fusion-negative subgroup is power-limited (insufficient to robustly support a "fusion-independent" claim as a main panel) and creates a reviewer-attack surface; the supplementary placement preserves the data while clarifying the boundary of the fusion-independence interpretation, which is now framed in Limitations §3.4 as "consistent with but not formally proving fusion-independent epigenetic silencing".</em>

Statistical tests: Cohen's d (pooled SD); Mann-Whitney U for per-gene β; Fisher exact for fusion × DM cross-tab.

---

# Supplementary Figures (S1-S9)

**S1.** 8-gene panel heatmap sorted by P_DM1 score, comparing TCGA versus K2 cohort distributions (Fig 1 C 이동, v3).

**S2.** 8-gene panel similarity matrix from external single-cell datasets (Pu 2021, Lu 2023, GSE241184).

**S3.** Pan-genome cluster ARI ladder spanning panel sizes 8, 16, 67 (TIERA67), 200, 1000, and 5000 by median absolute deviation.

**S4.** Immune-residualization analysis. DM1 vs DM2 Cohen's d for the 8-gene panel score before residualization (raw d = 1.78), after residualization on an antigen-presentation module (d = 1.00), after residualization on a generic immune signature (d = 1.50), and after dual residualization on both covariates (d = 0.87).

**S5.** Structural-variant missingness sensitivity analyses. Best-case, worst-case, and case-control-matched imputations for the 15 TCGA tumors lacking SV-tested status, showing preserved DM1-versus-DM2 fusion enrichment across scenarios.

**S5b.** *(v4 2026-05-08, demoted from main Figure 8C per audit P1-5.)* Methylation fusion-independence within DM1: fusion-positive (n = 63) versus fusion-negative (n = 19) tumors show equivalent mean 8-gene panel β-value (Cohen's d = −0.36, Mann-Whitney p = 0.31, NS). Interpretation: directionally consistent with a fusion-independent epigenetic silencing layer, but the small fusion-negative subgroup limits formal proof; the result is reported as supportive of, rather than definitive evidence for, parallel mechanism. See Limitations §3.4 for power discussion.

**S6.** *(v4 2026-05-08, contains both: original S6 sub-B detail + main Figure 7D phenotype panel demoted per audit P1-7.)* DM1 sub-A versus sub-B fusion-negative phenotype detail plots. (A) Age (sub-A 37.3 vs sub-B 51.3 years; Cohen's d = −0.82, MW p = 0.004). (B) Stage III/IV (15.3% vs 44.4%; OR 0.23, p = 0.020). (C) CD8/IFN-γ/checkpoint signatures (Cohen's d = −0.5 to −0.6 vs sub-B). (D) Hashimoto-like prevalence (sub-A 3.6% vs sub-B 12.5%; trend). Provided as supporting detail; deeper mechanistic dissection of the immune-overlap sub-B phenotype is reserved for the companion paper (Paper 2).

**SX.** *(NEW 2026-05-08, multi-method bulk cell-type deconvolution.)* Multi-method bulk cell-type deconvolution and canonical RAI-score residualization in TCGA-THCA. (A) Residualization Cohen's d grid for DM1 vs DM2 across four deconvolution methods (NNLS, Ridge-NNLS [L2 α=1], LR-clip, nu-SVR [CIBERSORT-style 3-ν]) and five regression stages (raw / + stromal / + immune / + epithelial-only / + all 8 fractions). Reference: Lu 2023 (GSE193581) `author_celltype` 67,678 cells × 8 cell types in 1,898 HVG ∩ TCGA gene symbols. Bulk: TCGA-THCA log2(TPM+1) with canonical labels `dm_like` (DM1_like 403, DM2_like 110) and canonical 8-gene RAI score `rai_score_recalc`. (B) Effect-retention ratio after residualization on all 8 cell-type fractions: NNLS / Ridge-NNLS 24% (sparse-weight artifact with T cell collapse), LR-clip 70% (over-permissive clip + renorm), nu-SVR 47% (methodologically aligned with canonical S4 immune-residualization 56% retention). (C) Per-cell-type DM1 vs DM2 fraction Cohen's d across the four methods, demonstrating direction-consistent enrichment of Malignant cell (d ≈ +1.08 to +1.35) and Myeloid cell (d ≈ +0.92 to +1.14) in DM1, and depletion of Epithelial cell (d ≈ −1.12 to −2.33; Lu 2023 cluster ≈ benign/normal thyroid epithelium, indicating sample-purity confound) and Endothelial cell. T cell DM1 vs DM2 d = −0.84 is uniquely resolved by nu-SVR (the three NNLS-family methods collapse T cell to 0). (D) nu-SVR primary mean cell-type fractions across TCGA-THCA n = 572. (E) Methodology and caveats panel. Together, panels A–C demonstrate that the bulk-level 8-gene RAI signal is partially explained by sample-level cell-type composition variation, but a non-trivial cohort-level signature persists after full residualization (47% retained, nu-SVR primary). The patient-level resolution comes from the Pu 2021 single-cell analysis (Figure 3; r = 0.798–0.886 between patient-matched tumor and normal thyrocyte 8-gene scores), which controls for sample purity within paired thyrocyte clusters. Source code, per-method fraction tables, and grids: `project/results/p_deconv_2026_05_08/`.

**S7.** Cross-cohort DM score portability in Korean validation cohorts. Score-distribution overlays and threshold-portability checks across K2, Lee/GSE213647, and the small Korean reference cohort used for calibration.

**S8.** Hypomethylating-agent + radioiodine re-induction schematic + literature meta. Decitabine + I-131 retrospective trial summary (NCT00085293, NCT01065090) and the rationale for prospective trial design stratified by DM1 status. (Fig 8 D 이동, v3.)

**S9.** K2 mini-index calibration diagnostic + alternate evidence chain. Per-gene inflation factors (4.9–12.5× across panel genes); 4-metric direction check (raw, within-sample-z, per-gene-z, per-gene-rank); alternate evidence from score-distribution portability and DM-call consistency.

---

# Additional supplementary figures (v17 audit + Dark matter phase 1/2 + Landa 2016)

These figures complement the main Fig 1-8 + S1-S9 set with sanity-check (v17 audit phase, 2026-04-29), single-cell foundation + multisite validation (Dark matter phase 1/2), and Discussion §3.1 cite-save evidence (Landa 2016 GSE76039 heatmap). Full panel-by-panel descriptions with verification markers are maintained in the parallel detail file `05_supp_figure_captions_v17_dm.md`.

## Part A — v17 Audit phase sanity-check

**SA1.** 4-way revalidation of the 8-gene DM1/DM2 cluster across BRAF V600E × TERT promoter mutational strata. (A-D) Stratum-specific Cox HR + KM curves; small-N caveat for BRAF−/TERT+ (n=4).

**SA2.** FFPE versus fresh-frozen tissue compatibility QC. (A) 8-gene panel score distribution by tissue type, Kolmogorov-Smirnov p=0.44 (no detectable shift). (B-C) Per-gene paired analysis + DM-call concordance.

**SA3.** Panel-size sensitivity. (A) 5-fold cross-validated AUC across 8/10/12/16-gene variants; ΔAUC 8 vs 16 = 0.013, NS. (B-C) Cluster-call concordance matrix + ARI ladder.

**SA4.** Single-cell wrap-up — thyrocyte-intrinsic 8-gene signal across GSE184362 (Pu et al., 2021), GSE193581 (Lu 2023), GSE241184 (Phase 1). (A-B) Per-cohort summary + per-patient r 0.798–0.886, all p < 10⁻¹⁰. (C) Author-independence cross-cohort concordance.

**SA5.** Differentiation trajectory along the 8-gene panel score in TCGA-THCA. (A) Pseudotime ordering of primary tumors. (B-D) Driver mutation distribution, differentiation gene expression, and DM cluster overlay along the trajectory.

**SA6.** MSK-IMPACT bias panel — advanced-disease cohort (Landa et al., 2016) labeled explicitly. (A) PDTC n=84 + ATC n=33 composition. (B-D) Mutation-frequency, stage, and DM-prevalence comparison versus TCGA-THCA primary tumors.

## Part B — Dark matter phase 1 single-cell UMAP foundations

**SB1.** Single-cell UMAP overview of GSE184362 (Pu et al., 2021; n=6 PTC patients, Fudan University). (A) Embedding colored by patient identity. (B) Cell-type annotation (thyrocyte / immune / stromal / endothelial). (C) Tumor versus adjacent normal labeling.

**SB2.** Single-cell UMAP — molecular signatures. (A) 8-gene panel score per cell (continuous gradient). (B) HLA-II module signature (Paper 1 residualization control only; deeper HLA-II analysis = Paper 2 territory). (C) Canonical thyroid differentiation transcripts (TG, TPO, TSHR averaged). (D) Proliferation signature (MKI67, TOP2A).

**SB3.** Thyrocyte-restricted UMAP. (A) Subset to KRT8+ KRT19+ EPCAM+ cells. (B) 8-gene panel score gradient on thyrocyte UMAP. (C) Tumor versus adjacent-normal labeling. (D) Per-patient thyrocyte distribution along the 8-gene score axis.

## Part C — Dark matter phase 2 multisite + per-patient validation

**SC1.** Per-patient Spearman r forest plot — tumor versus adjacent-normal thyrocyte 8-gene scores in GSE184362 (n=6 patients; per-patient r 0.798–0.886; Bonferroni-adjusted p < 10⁻¹⁰; per-patient n_cells annotated).

**SC2.** Multisite trajectory — DM1/DM2 score reproducibility across TCGA-THCA (n=504), MSK-IMPACT advanced-disease (n=117), K2 / PRJEB11591 (n=260), and Lee / GSE213647 (n=632).

**SC3.** Multisite single-cell UMAP — joint cell embedding across GSE184362 (Pu 2021), GSE193581 (Lu 2023), and GSE241184 (Phase 1) after batch-corrected integration. Dataset-independent gradient for the 8-gene panel score.

**SC4.** GSE184362 (Pu 2021) summary table — per-patient metadata, n_cells per condition, per-patient Spearman r, DM1/DM2 prediction.

**SC5.** External-validation pooled scatter — per-patient pseudobulk 8-gene score in tumor versus adjacent-normal thyrocytes across GSE184362 + Lu 2023.

**SC6.** Pooled scatter — DM1 probability versus 8-gene panel score across external single-cell cohorts; per-cohort markers + pooled regression with 95% CI.

**SC7.** K2 (PRJEB11591) versus Yoo 2016 reference — 8-gene score concordance + K2 mini-index calibration mismatch diagnostic (R4-4 audit; per-gene inflation factors 4.9–12.5×; within-sample-centered profile restores TCGA direction).

## Part D — Landa 2016 evidence (Discussion §3.1 cite save)

**SD1.** Landa et al. (2016) GSE76039 — differentiation transcript suppression in advanced thyroid cancer. (A) Per-sample heatmap of seven canonical differentiation transcripts (TG, TSHR, TPO, PAX8, SLC26A4, DIO1, DUOX2) across PDTC + ATC samples from the GSE76039 transcriptome subset. (B) 5-of-8 overlap with this paper's panel (TG, TSHR, TPO, PAX8, DIO1). (C) Convergence framing — Yoo et al. (2016) panel origin and Landa et al. (2016) advanced-disease list reach the same differentiation core via independent paths (Discussion §3.1 reverse-causality framing).

Cite: Landa I, Ibrahimpasic T, Boucai L, Sinha R, Knauf JA, Shah RH, Dogan S, Ricarte-Filho JC, Krishnamoorthy GP, Xu B, Schultz N, Berger MF, Sander C, Taylor BS, Ghossein R, Ganly I, Fagin JA. Genomic and transcriptomic hallmarks of poorly differentiated and anaplastic thyroid cancers. *J Clin Invest.* 2016;126(3):1052–1066. doi:10.1172/JCI85271. PMID 26878173. PMC4767360.

For full panel-by-panel descriptions, exact statistical-test specifications, and the 17-item verification queue (per-figure ⚠ markers), see `05_supp_figure_captions_v17_dm.md`.

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
