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

<em>Mechanism layer cross-reference (2026-05-09).</em> Upstream-effector candidates and panel-canonicality checks are reported as supplementary support: MAPK-pathway output × HM450 mean 8-gene β methylation co-variation (BRAF V600E β = 0.37 ≈ RET fusion β = 0.39 vs RAS β = 0.27; Supplementary Figure SX panels H–J / `v5A_per_driver_class.tsv`); MAPK output × 8-gene panel z RNA Spearman ρ = −0.291 (TCGA n = 572) / −0.395 (Lee n = 632); DM1 sub-A vs sub-B two-axis convergence (sub-A MAPK-active d = +1.79, sub-B HT-active d = +0.62) both reaching panel silencing (panel d = +0.07 NS); and reproduction on the canonical Yoo 2014 TDS-16 (TDS-16 ≈ Panel-8: ΔAUC for MAPK-high classification +0.007 / +0.012 NS, sub-A vs sub-B TDS-16 d = +0.03 NS) — Supplementary Figure SX_v13 / `v13_*.tsv`. The mechanism layer is reported as supportive, not causally proven; the Q9 framing in `09_reviewer_qa.md` retains the "motivates rather than confirms" boundary.

<em>Reserve extension cross-reference (2026-05-09 v15).</em> Supplementary Figure SX_v15 adds K2 score-only overlay (246/260 K2 samples DM2-like using the TCGA-centered 8-gene classifier), a GSE250521 spatial stress test (direct per-spot MAPK × Panel association is positive, not anti-correlated; QC-partial pooled ρ = +0.059), and PRISM/DepMap drug-vulnerability support (7/11 FDR < 0.05 PRISM hits are canonical MAPK-axis inhibitors; AZD-0364 MEK d = −0.594, FDR = 5.9 × 10⁻⁷). These panels are reviewer-reserve support only and do not change the main Figure 8 panel count.

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

**SX.** *(UPDATED 2026-05-08 v5, multi-method deconvolution + compositional axes; cumulative scope n = 1,241 bulk samples = TCGA 572 + Lee/GSE213647 632 + GSE76039 37.)* Multi-method bulk cell-type deconvolution and compositional axes of the DM1↔DM2 phenotype. **(A) Multi-method residualization grid.** Cohen's d (DM1 − DM2) for the canonical 8-gene RAI score (`rai_score_recalc`) under five residualization stages (raw / + stromal [Endo + Fibro] / + immune [T + Myel + B + NK] / + Epithelial-only / + all 8 fractions) across four deconvolution methods (NNLS, Ridge-NNLS [L2 α=1], LR-clip, nu-SVR [CIBERSORT-style; 3 ν 0.25/0.5/0.75, lowest-RMSE]). Reference: Lu 2023 (GSE193581) `author_celltype` 67,678 cells × 8 cell types in 1,898 HVG ∩ TCGA gene symbols. Bulk: TCGA-THCA log2(TPM+1), n = 572, canonical labels `dm_like` (DM1_like 403 / DM2_like 110). **(B) Effect-retention ratio after full residualization:** NNLS / Ridge-NNLS 24% (sparse-weight artifact with T-cell collapse), LR-clip 70%, nu-SVR 47% (methodologically aligned with canonical S4 immune-residualization 56% retention). **(C) Per-cell-type DM1 vs DM2 fraction Cohen's d** across the four methods. Direction-consistent: DM1 enriched for Malignant (d ≈ +1.08 to +1.35) and Myeloid (d ≈ +0.92 to +1.14); depleted for Epithelial (d ≈ −1.12 to −2.33; benign/normal thyroid epithelium proxy) and Endothelial. T cell DM1 vs DM2 d = −0.84 uniquely resolved by nu-SVR. **(D) nu-SVR primary mean cell-type fractions** in TCGA-THCA (n = 572). **(E) Methodology and caveats panel.** **(F) Cross-cohort direction consistency.** TCGA nu-SVR Cohen's d (DM1 − DM2) versus Lee/GSE213647 (n = 632) Spearman ρ vs `panel_z` (sign-flipped to align with DM1 direction). Direction-consistent for 7/8 cell types; T-cell discordance reflects binary-vs-continuous score frame. **(G) Within-DM1 fusion+ vs fusion− cell-type fraction Cohen's d** (n = 74 vs 391; kinase fusion: RET/NTRK/ALK/BRAF/PAX8/PPARG aggregated from cBioPortal SV table). All |d| ≤ 0.42 — within-DM1 fusion+/− tumors share near-identical compositions, direct support for the fusion-independent epigenetic silencing claim (Figure 8 / §2.5a). **(H) Per-driver-class TCGA-THCA cell-type composition** (BRAF V600E n = 280 / RAS-mutant n = 54 / driver-negative n = 179). RAS-mutant tumors are Epithelial-cluster–enriched (d_RAS−BRAF = +1.52) and Malignant-cluster–depleted (d = −1.66) versus BRAF V600E, consistent with Landa 2016 / Paper 1 §2.1 differentiation framing. **(I) HM450 8-gene mean β × cell-type fraction Spearman ρ heatmap** (TCGA n ≈ 484 paired): mean β positively tracks Myeloid (ρ = +0.39), B cell (+0.27), Fibroblast (+0.23) and Malignant (+0.30); negatively tracks T cell (ρ = −0.42) and Epithelial (−0.31), connecting the Round-4 methylation layer (DM1 vs DM2 mean-β d = −1.75) to compartment composition. **(J) DM1 sub-A vs sub-B cell-type fraction Cohen's d** (sub-A n = 84, sub-B n = 56). Sub-A is Malignant-cell rich (d = +1.22, p = 2.5 × 10⁻⁹) and Epithelial-poor (d = −1.26, p = 5.2 × 10⁻¹⁰); immune compartment differences are not significant. The sub-A/B split is therefore a tumor-purity-vs-thyrocyte split, not immune-hot vs immune-cold — Paper-2 boundary marker (sub-B = fusion-/mutation-negative Hashimoto-overlap retains thyrocyte identity). **(K) Cell-type composition pseudotime along the canonical 8-gene score.** TCGA-THCA (n = 513, score = `rai_score_recalc`) and Lee/GSE213647 (n = 632, score = `panel_z`) samples were ranked by score, binned into 10 deciles; the mean per-decile cell-type fraction (nu-SVR against Lu 2023) defines a 10-step pseudotime trajectory. Four compartments — Malignant (↓), Epithelial (↑), Myeloid (↓), Endothelial (↑) — show monotonic decile-level Spearman ρ ≥ |0.95| in **both** cohorts, defining a reproducible compositional pseudotime independent of cohort, scoring scheme, or sample size. Source code, per-method fraction tables, panel-level TSVs, and the full v5 composite (panels A–K) at `project/results/p_deconv_2026_05_08/` (build: `plot_deconv_v5_composite.py`; per-panel: `run_deconv_v5_{abc, d_final, e_trajectory}.py`; output figure: `Fig_SX_deconvolution_v5_composite.{png, pdf}`).

**SX_v13.** *(NEW 2026-05-09; Paper 1 Fig 8 mechanism support; cumulative scope n = 1,204 = TCGA 572 + Lee/GSE213647 632.)* The MAPK→thyroid-silencing axis operates indistinguishably on the deployable 8-gene panel and the canonical Yoo 2014 TDS-16. **(A) Cross-cohort MAPK output × thyroid-score Spearman ρ.** TCGA-THCA n = 572 and Lee/GSE213647 n = 632 z-mean MAPK output (DUSP4/5/6, SPRY2/4, ETV4/5, PHLDA1, CCND1; n = 9 genes) versus 8-gene Panel (deployable: DIO1/FOXE1/NKX2-1/PAX8/SLC5A5/TG/TPO/TSHR), TDS-16 (Yoo 2014 canonical: Panel + DIO2/DUOX1/DUOX2/GLIS3/SLC26A4/SLC5A8/THRA/THRB), and TDS−panel (the 8 disjoint TDS-16 genes). Panel-8 ρ = −0.291 (p = 1.3 × 10⁻¹²) / −0.395 (p = 4.7 × 10⁻²⁵); TDS-16 ρ = −0.306 / −0.435; TDS−panel ρ = −0.310 / −0.449. **(B) Per-driver-class TCGA score means** (BRAF V600E n = 344 / RAS-mutant n = 61 / RET fusion n = 43 / NTRK fusion n = 10 / driver-negative n = 103) for MAPK output, Panel-8, TDS-16, and TDS−panel z-scores; the Panel-8 and TDS-16 traces are essentially superimposable, with BRAF V600E vs RAS-mutant Cohen's d Panel-8 = −1.615 versus TDS-16 = −1.616 (identical to three significant figures). **(C) DM1 sub-A vs sub-B Cohen's d** (sub-A n = 93, sub-B n = 62; labels from `d6p7_dm1_subcluster/dm1_subcluster_labels.tsv`) for five metrics: MAPK d = +1.79 (Mann-Whitney p = 4.1 × 10⁻¹⁷), HT signature d = −0.12 (NS), Panel-8 d = +0.07 (NS), TDS-16 d = +0.03 (NS), TDS−panel d = −0.03 (NS). The v12 two-axis convergence (MAPK-active sub-A or HT-active sub-B, both reach 8-gene silencing) reproduces on the full TDS-16. **(D) Per-gene MAPK output × thyroid-gene Spearman ρ heatmap** for the 16 TDS-16 genes in TCGA-THCA and Lee; ★ marks the 8-gene Panel members. Eight-gene members occupy median rank within the 16-gene heatmap (not selectively top-extreme: strongest TCGA negatives = SLC5A8 −0.516 [non-panel], DIO2 −0.481 [non-panel], TPO −0.465 [panel]; NKX2-1 is the consistent positive-ρ outlier in both cohorts, +0.307 / +0.425, a known panel-mean-absorbed caveat). **(E) MAPK-decile pseudotime.** Samples ranked by MAPK output and binned into 10 deciles; mean Panel-8, TDS-16, and TDS−panel score per decile in both cohorts. The two panels trace identical monotonic descent (decile-rank ρ Panel = −0.600 / −0.648, TDS-16 = −0.600 / −0.818). **(F) ROC-AUC for MAPK-high vs MAPK-low classifier** using Panel-8 versus TDS-16 versus TDS−panel score (sign-flipped, low score = high MAPK). Panel-8 AUC = 0.623 (TCGA) / 0.728 (Lee); TDS-16 AUC = 0.631 / 0.740; ΔAUC TDS-16 − Panel-8 = +0.007 (TCGA) / +0.012 (Lee), both non-significant — independently reproducing the existing manuscript p1_onepage_audit ΔAUC = 0.013 NS (DM1/DM2 task) on the orthogonal MAPK-classification task. *Methods.* z-score = within-cohort z-mean of the named gene set; gene lists 100% recovered in both cohorts (Lee via `F1_gene_recovery_mapping.tsv` Ensembl→symbol). Driver class from cBioPortal SV/MAF anchors (`audit_2026_04_30/round3/cbio_sv_thca.tsv`). DM1 sub-A/sub-B labels from `d6p7_dm1_subcluster`. Source code, panel-level TSVs, and the 6-panel composite at `project/results/p_deconv_2026_05_08/` (build: `plot_v13_tds16_panel_mapk.py`; pipeline: `run_v13_tds16_panel_mapk.py`; output figure: `Fig_SX_v13_TDS16_MAPK.{png, pdf}`).

**SX_v14.** *(NEW 2026-05-09 evening; Paper 1 Fig 8 / Q14 cross-cohort generalizability lock; cumulative scope n = 1,287 = TCGA 572 + Lee 632 + GSE126698 28 + GSE286332 18 + GSE76039 37.)* Cross-cohort forest of MAPK output × thyroid-score Spearman ρ. Three score panels (8-gene Panel deployable / TDS-16 Yoo 2014 canonical / TDS−panel 8 disjoint genes) × five cohorts. Squares = per-cohort Spearman ρ (size proportional to √n); horizontal bars = 95% CI (Fisher z-transform); diamond = pooled fixed-effect Fisher-z ρ. **Pooled MAPK × Panel-8 ρ = −0.327, 95% CI [−0.376, −0.278], n_total = 1,287.** Cochran Q = 14.77 (df = 4, p = 0.005), I² = 72.9% — heterogeneity is biologically expected and predicted by the v12 two-axis convergence model (see Q9 / Q14 in `09_reviewer_qa.md`): GSE286332 (Korean PTC vs PTC+HT, n = 18) is dominated by the HT route where v12 predicts MAPK decouples (ρ = −0.040 NS), and GSE76039 (Landa 2016 PDTC + ATC, n = 37) shows panel saturation at the dedifferentiated end of the axis (ρ = +0.125 NS). The two well-differentiated primary cohorts (TCGA + Lee, n = 1,204) carry the entire signal (ρ = −0.291 / −0.395, both p ≪ 10⁻¹²). Source code, panel TSV, and the 3-panel forest figure at `project/results/p_deconv_2026_05_08/` (build: `plot_v14_cross_cohort_forest.py`; pipeline: `run_v14_cross_cohort_forest.py`; output figure: `Fig_SX_v14_cross_cohort_forest.{png,pdf}`; per-cohort table: `v14_cross_cohort_forest.tsv`).

**SX_v15.** *(NEW 2026-05-09 evening; Paper 1 Fig 8 reviewer-reserve extensions; main Figure 8 panel count unchanged.)* K2 score-only overlay, GSE250521 spatial stress test, and PRISM/DepMap drug-vulnerability overlay. **(A) K2 mini-index overlay.** A TCGA-trained centered-profile 8-gene classifier gives TCGA 5-fold CV AUC = 0.960 and classifies 246/260 K2 PRJEB11591 samples as DM2-like (median p_DM2 = 0.978). K2 lacks MAPK-output genes, so this panel is score-distribution evidence only and is not included in the MAPK × Panel forest. **(B) Spatial Lu 2023 GSE250521 MAPK × Panel test.** Raw h5ad spots were log1p(CP10K)-normalized, within-slide z-scored, and summarized as MAPK-9 and Panel-8 scores. Across 12 tumor-region slides (PTC/LPTC/ATC; 43,180 spots), the fixed-effect per-slide MAPK × Panel association is positive rather than negative (ρ = +0.326 [0.318, 0.335]); QC-partial ρ after total-count, gene-count, and mitochondrial-fraction adjustment is much smaller (+0.059 [0.049, 0.068]). This panel is a non-confirmatory stress test, not a spatial mechanism lock. **(C) PRISM overlay.** Among 1,518 PRISM drugs, 11 are FDR < 0.05 and DM1-high selective; 7/11 are canonical MAPK-axis inhibitors (MEK/RAF/ERK). The top hit is AZD-0364 (MEK; Cohen's d = −0.594, FDR = 5.9 × 10⁻⁷). Across 690 nonmissing cell lines, DM1 score anti-correlates with canonical MAPK-inhibitor mean LFC (Spearman ρ = −0.236, p = 3.3 × 10⁻¹⁰), and high-vs-low DM1 bins differ by d = −0.617. **(D) DepMap overlay.** Top DM1-high dependencies are MYC (d = −0.499, p = 1.0 × 10⁻¹¹) and NAMPT (d = −0.445, p = 1.5 × 10⁻⁹); thyroid transcription factors are not the primary therapeutic lever. PRISM/DepMap supports MAPK-axis vulnerability but does not measure restoration of thyroid-gene expression after inhibitor treatment. Source code and outputs: `run_v15_mechanism_extensions.py`, `v15A_*`, `v15B_*`, `v15C_*`, `v15_mechanism_extensions_summary.json`, and `Fig_SX_v15_mechanism_extensions.{png,pdf}` in `project/results/p_deconv_2026_05_08/`.

**SX_v16.** *(NEW 2026-05-09 evening; reviewer-reserve diagnostic stress test.)* Spatial rescue failure and PRISM claim-boundary lock. **(A) GSE250521 MAPK × Panel adjustment grid.** Across 57,997 spatial spots and 43,180 tumor-stage spots, raw MAPK × Panel association is positive; QC and cell-state partialing attenuate but do not flip the direction. Tumor-slide pooled ρ = +0.326 [0.318, 0.335] raw, +0.058 [0.048, 0.068] after QC adjustment, and +0.034 [0.024, 0.044] after QC plus epithelial/CAF/EMT/hypoxia/proliferation adjustment. **(B) Epithelial-score quartile grid.** Stage × within-slide epithelial-score quartiles show positive or near-zero MAPK × Panel relationships; no robust anti-correlation emerges. **(C) PRISM top-k enrichment sensitivity.** The correct claim is top 7/7 canonical MAPK-axis inhibitors, 7/11 FDR < 0.05 canonical MAPK-axis inhibitors (hypergeometric p = 3.25 × 10⁻¹⁵), and 8/15 top-15 canonical MAPK-axis inhibitors; the earlier shorthand "top 15 all MAPK" is not used. **(D) PRISM lineage sensitivity.** DM1 score × canonical MAPK-inhibitor mean LFC remains negative when thyroid cell lines are excluded (ρ = −0.240, p = 2.27 × 10⁻¹⁰); thyroid-only cell-line correlation is not interpretable because the available thyroid lines share a constant DM1 score. Source code and outputs: `run_v16_spatial_prism_diagnostics.py`, `v16_spatial_*`, `v16_prism_*`, `v16_spatial_prism_diagnostics_summary.json`, and `Fig_SX_v16_spatial_prism_diagnostics.{png,pdf}` in `project/results/p_deconv_2026_05_08/`.

**SX_v17.** *(NEW 2026-05-09 evening; spatial signal decomposition.)* Why the GSE250521 spatial MAPK × Panel stress test fails as mechanism confirmation. **(A) Adjustment-family ladder.** In tumor epithelial-top50 spots, raw MAPK × Panel ρ = +0.208, but full covariate adjustment reduces the linear-residual ρ to +0.052 and full spatial+detection adjustment reduces it to +0.007. **(B) Covariate R².** Detection/depth covariates explain much of both modules (median tumor epithelial R²: MAPK = 0.809; Panel-8 = 0.683). **(C) Depth/detection structure.** MAPK-detection and Panel-detection are positively correlated (median ρ = +0.383), consistent with broad Visium detection co-localization. **(D) MAPK submodule tests.** DUSP/SPRY, ETV/PHLDA1, no-CCND1, and CCND1-only variants do not recover a strong full-adjusted anti-correlation. **(E) Gene-pair sign consistency.** The most negative pair is CCND1 × TPO (median ρ = −0.096), but gene-pair FDR is non-significant (FDR = 0.694). **(F) Random-module null.** Observed MAPK × Panel raw median correlation is near the random-module median (percentile = 0.545), and full-residual observed correlation is below the random-module median (percentile = 0.343). These diagnostics support treating GSE250521 as a spatial caveat rather than positive Fig 8 mechanism evidence. Source code and outputs: `run_v17_spatial_signal_decomposition.py`, `v17_spatial_*`, `v17_spatial_signal_decomposition_summary.json`, and `Fig_SX_v17_spatial_signal_decomposition.{png,pdf}`.

**SX_v18.** *(NEW 2026-05-09 evening; spatial lag and pocket closure.)* Final GSE250521 spatial rescue test. **(A) Same-spot and neighborhood-lag MAPK × Panel correlations** in tumor epithelial-top50 spots. Raw same-spot median ρ = +0.110; raw KNN 1-6 neighborhood-lag median ρ = +0.090; raw KNN 7-18 median ρ = +0.071. Full spatial+detection residual same-spot median ρ = −0.011; residual KNN 1-6 median ρ = +0.033; residual KNN 7-18 median ρ = +0.010. **(B) MAPK-high/Panel-low anti-pocket enrichment.** Raw anti-pocket enrichment = 0.882 (OR = 0.693), so anti-pockets are depleted rather than enriched; residual anti-pocket enrichment = 1.010 (OR = 1.018), approximately independence-level. **(C) Concordant MAPK-high/Panel-high pockets** are not strongly enriched (raw enrichment = 1.002; residual enrichment = 0.982), consistent with broad co-detection rather than focal anti-silencing domains. **(D) Anti-pocket covariates.** The strongest raw anti-pocket covariate depletion is Panel detection (Cohen's d = −1.094), indicating low panel-detection spots rather than robust biological MAPK-high/Panel-low neighborhoods. These panels close the last spatial escape hatch: GSE250521 is retained as a Visium-resolution/detection caveat, not Fig 8 mechanism support. Source code and outputs: `run_v18_spatial_lag_pockets.py`, `v18_spatial_lag_*`, `v18_spatial_pocket_*`, `v18_spatial_lag_pockets_summary.json`, and `Fig_SX_v18_spatial_lag_pockets.{png,pdf}`.

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

---

# Narrative tightening notes (2026-05-11)

The manuscript reads best when the main figures are treated as a clean claim ladder rather than a list of analyses.

## Main-text priority order

1. **Figure 1** should carry the compact panel-defining claim only.
2. **Figure 2** should handle clinical relevance only.
3. **Figure 3** should establish single-cell intrinsic validation.
4. **Figure 4** should stop at driver-context stratification and not overreach mechanistically.
5. **Figure 5** should be the portability / assay-compatibility proof.
6. **Figure 6** should remain the survival-risk summary.
7. **Figure 7** should carry the fusion-mechanism support.
8. **Figure 8** should carry the epigenetic support and explicitly stay correlative.

## What should stay out of the main figures

- calibration mismatch diagnostics
- spatial stress-test caveats
- niche reserve extensions
- method-comparison stress tests
- any counterexample that is useful but not claim-building

## What makes the paper stronger

- fewer claims per figure
- more separation between discovery and mechanism
- cleaner boundary language in captions
- fewer supplementary panels inside the main-figure storyline
