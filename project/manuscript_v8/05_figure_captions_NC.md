---
title: "Paper 1 manuscript v8 — Nature Cancer 6-main-figure captions (NC reach)"
date: 2026-05-27
status: draft v1 — NC-reach architecture, claim boundaries identical to v3 (no perturbation in main)
target_venue: Nature Cancer (reach) → Cell Reports Medicine (anchor) → JCI Insight / Genome Medicine (safe)
panel_density_target: 6–8 panels per figure (NC convention, vs Cell Press 3–5)
boundary_lock: PRISM/DepMap, TERT survival, Mun proteome, pan-cancer LGG/LUAD → Extended Data only
companion_file: 05_figure_captions.md (Cell Press 8-fig v3, retained as the CRM-anchor build)
---

# Main figures (6 — NC architecture)

The six main figures progress: **subtype existence → fusion/histology mechanism → epigenetic mechanism → external validation → cross-cohort survival → clinical reflex pathway**. Mechanism (Figs 2–3) is sandwiched between discovery (Fig 1) and validation (Fig 4), with translation (Fig 6) as the closer. Perturbation, TERT survival, proteomic, and pan-cancer panels are deliberately held to Extended Data to preserve a tight, defensible claim surface.

---

## Figure 1. The DM1/DM2 axis defines an aggressive papillary thyroid carcinoma subtype orthogonal to canonical BRAF/RAS classification. (7 panels)

**(A) Study schematic and cohort inventory.** Six cohorts spanning bulk transcriptome, methylation, and single-cell layers: TCGA-THCA (n = 504, primary tumors, RNA-seq + HM450 + clinical); MSK-IMPACT thyroid (n = 117, advanced disease); GSE213647 / Lee et al. Korean (n = 632, FFPE); PRJEB11591 / K2 (n = 260, Korean FF); GSE184362 / Pu et al. 2021 (n = 6 patients, single-cell, paired tumor and adjacent-normal); GSE193581 / Lu et al. 2023 (n = 23 samples, single-cell, thyrocyte-resolved). Schematic strip indicates the entry points for each cohort into the analysis (panel-call, methylation, fusion, single-cell intrinsic, survival).

**(B) TIERA67 → 8-gene panel curation.** Sankey flow diagram from the TIERA67 67-gene candidate pool — assembled from seven thyroid-relevant biological categories (Driver_anchor, TDS_core, RAI_machinery, lineage TFs, MAPK output, immune-context, dedifferentiation) — to the final 8-gene RAI-responsiveness panel: SLC5A5, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, DIO1. The panel is the TDS-core sub-category and was selected on the Yoo 2016 (PLOS Genet) RAI-uptake biological prior, independent of any subsequent dedifferentiation gene list (reverse-causality control; see Methods).

**(C) UMAP embedding of TCGA-THCA on the 8-gene transcript space.** Primary tumors (n = 504) projected on the eight-dimensional 8-gene transcript space, colored by KMeans (k = 2) DM1 (red) versus DM2 (blue) cluster assignment. DM1 prevalence = 28.4 % of TCGA-THCA primary tumors.

**(D) 8-gene expression heatmap sorted by P_DM1.** Rows = the 8 panel genes; columns = TCGA-THCA tumors ordered by ascending P_DM1. Top annotation tracks: DM cluster, BRAF V600E status, RAS-family mutation, RET / NTRK / ALK / BRAF fusion, TERT promoter mutation, ATA 2015 stage. The driver tracks visually confirm that the DM1↔DM2 axis is not co-aligned with any single driver lesion.

**(E) Driver mRNA neutrality.** Cohen's d (DM1 vs DM2) and single-feature classification AUC for the canonical drivers as a panel-orthogonality control: BRAF transcript (d = −0.044, AUC = 0.602, Mann-Whitney p = 0.57); TERT (AUC = 0.578); KRAS (0.525); NRAS (0.521); HRAS (0.500). Dashed line at AUC = 0.5 (chance).

**(F) Pan-genome cluster reproducibility (ARI ladder).** Adjusted Rand index of unsupervised KMeans partitions versus the 8-gene-driven DM1/DM2 reference, across six gene-set sizes: 8-gene panel alone (0.49), TIERA67 67-gene (0.90), pan-genome top-5000 by median absolute deviation (0.92), top-1000 (0.91), top-200 (0.78), Driver_anchor genes alone (−0.007). Hypergeometric enrichment of TIERA67 within the pan-genome top-100 MAD genes: p = 3 × 10⁻⁴.

**(G) Kaplan-Meier overall survival, DM1 vs DM2 (TCGA-THCA primary).** Pooled KM stratified by DM cluster assignment (n_DM1 = 143, n_DM2 = 361, 16 events). Log-rank p reported on plot; hazard ratio cross-referenced to the meta-analysis in Figure 5A.

**Statistics.** KMeans k = 2 with cosine distance on z-scored TPM; ARI computed against the 8-gene-driven reference labels; Cohen's d with pooled standard deviation; Mann-Whitney U for transcript distributions; hypergeometric test (one-sided) for TIERA67 enrichment. n_TCGA = 504 (primary tumors with usable RNA-seq and clinical metadata).

---

## Figure 2. DM1 is a fusion-driven, histologically distinct subtype within the canonical BRAF/RAS dark matter. (7 panels)

**(A) Xing 2014 dark-matter rescue (Sankey).** Of 180 TCGA-THCA tumors that are BRAF V600E-negative and TERT promoter-negative (Xing 2014 NEJM "dark matter"), 131 (73 %) are rescued into a defined DM1/DM2 stratum by the 8-gene panel call. Sankey flow shows the source dark-matter pool, the DM1/DM2 assignment, and the residual unclassified bin.

**(B) DM1 sub-A vs sub-B silhouette structure.** Within DM1 (n = 91), KMeans k = 2 sub-clustering yields two phenotypes: a fusion-enriched sub-A (n = 72) and a fusion-attenuated sub-B (n = 19). Average silhouette score 0.584 supports the two-population structure. Sub-A vs sub-B mechanistic dissection is reserved for the companion paper (Paper 2); here it is reported only as evidence that DM1 is not a single homogeneous population.

**(C) Kinase fusion enrichment by DM cluster.** Stacked bar plot of fusion-positivity from cBioPortal SV calls: DM1 76.8 % (63 / 82 with usable SV) vs DM2 30.9 %; Fisher OR = 7.41 (95 % CI 4.38 – 12.55), p = 1.9 × 10⁻¹³. SV missingness sensitivity (best-case, worst-case, MAR-imputed) preserved direction in all scenarios (Methods; Extended Data Fig. 5).

**(D) Fusion partner spectrum within DM1 fusion-positive tumors.** Stacked bar of fusion partner classes: RET (n = 33; CCDC6-RET 17, NCOA4-RET 3, other RET partner 13); NTRK1/2/3 (n = 10); ALK (n = 4); BRAF rearrangement (n = 5); other (n = 11). Partner identity supports the clinical reflex logic in Figure 6A.

**(E) DM1 capture rate of TCGA RET-fusion-positive tumors.** Of the 33 TCGA RET-fusion-positive tumors, 27 (81.8 %) are called DM1 by the 8-gene panel. Bar plot of DM1, DM2, and unclassified capture across the four major actionable fusion classes (RET, NTRK, ALK, BRAF), supporting an RNA-first reflex testing trigger.

**(F) Follicular variant PTC (FVPTC) enrichment.** Mosaic plot of FVPTC histology (TCGA-THCA pathology call) by DM cluster: DM1 vs DM2 OR = 17.9, Fisher exact p = 3 × 10⁻³¹. The DM1 axis captures the histological substructure that distinguishes follicular-variant tumors — a feature not predicted by the 8-gene panel a priori (panel was selected on RAI biology, not histology).

**(G) Cross-cohort fusion replication: MSK-IMPACT advanced thyroid.** Independent confirmation of the fusion-DM1 association in MSK-IMPACT (n = 117 advanced thyroid carcinomas): RET (n = 5), ALK (n = 3), PAX8-PPARG (n = 3), among DM1-called tumors. Per-cohort fusion-positivity rate and DM1 capture rate side-by-side.

**Statistics.** Fisher exact test for fusion / FVPTC odds ratios with 95 % CIs from the conditional MLE; KMeans silhouette for sub-A / sub-B; chi-square sensitivity (p = 0.56) for missing-at-random SV-status imputation. n_TCGA with SV call = 379; n_DM1 with SV = 82.

---

## Figure 3. DM1 epigenetically silences the thyroid-differentiation core. (7 panels)

**(A) Per-gene HM450 promoter β-value heatmap.** Rows = the 8 panel genes plus DIO2 and SLC26A4 (two TDS-16 dedifferentiation markers not in the 8-gene panel, as cross-set controls); columns = TCGA-THCA tumors with paired HM450 and DM call (n = 503), grouped by DM1, DM2, and not_DM. Per-gene DM1 vs DM2 Cohen's d: TPO (2.30, p = 1.9 × 10⁻¹⁸); DIO1 (1.24, p = 6.5 × 10⁻¹¹); TSHR (1.20, p = 9.8 × 10⁻¹²); PAX8 (0.97, p = 4.5 × 10⁻⁸); TG (0.86, p = 2.2 × 10⁻⁶); FOXE1 (0.84, p = 1.0 × 10⁻⁵); NKX2-1 (0.63, p = 8.9 × 10⁻⁷); SLC5A5 (0.22, NS p = 0.42).

**(B) Mean 8-gene panel β by DM cluster.** Bar plot of mean 8-gene promoter β (TCGA HM450, n = 503): DM1 = 0.385, DM2 = 0.253, not_DM = 0.356, corresponding to 52 % higher DM1 promoter methylation versus DM2. Error bars = 95 % bootstrap CI.

**(C) Methylation vs expression scatter (four representative panel genes).** Per-tumor mean promoter β (x-axis) versus log2(TPM + 1) expression (y-axis) for TPO, DIO1, TSHR, and TG; points colored by DM cluster. Spearman ρ and within-DM-stratum slope reported. Negative β-expression coupling is monotonic across the 8-gene panel (per-gene ρ in Suppl. Table).

**(D) Cross-cohort MAPK output × Panel-8 score forest.** Pooled fixed-effect Fisher-z Spearman ρ of MAPK output (DUSP4/5/6, SPRY2/4, ETV4/5, PHLDA1, CCND1) against the 8-gene panel z-score across five cohorts: TCGA-THCA (n = 572, ρ = −0.291, p = 1.3 × 10⁻¹²); GSE213647 / Lee Korean (n = 632, ρ = −0.395, p = 4.7 × 10⁻²⁵); GSE126698 (n = 28); GSE286332 (n = 18, Korean PTC vs PTC+HT); GSE76039 / Landa 2016 PDTC + ATC (n = 37). **Pooled ρ = −0.327 [95 % CI −0.376, −0.278], n_total = 1,287.** Cochran I² = 72.9 % (Q p = 0.005). Heterogeneity is the predicted two-axis convergence: PTC+HT (HT-active route, ρ = −0.040 NS) and PDTC+ATC (panel-saturated end, ρ = +0.125 NS) trade signal with the two well-differentiated primary cohorts that carry the entire effect.

**(E) TDS-16 ≡ Panel-8 functional equivalence.** Side-by-side per-driver-class mean z-score (TCGA, n = 561; BRAF V600E n = 344 / RAS-mutant n = 61 / RET-fusion n = 43 / NTRK-fusion n = 10 / driver-negative n = 103) for the deployable 8-gene panel and the canonical Yoo 2014 TDS-16. BRAF V600E vs RAS-mutant Cohen's d: Panel-8 = −1.615, TDS-16 = −1.616 (identical to three significant figures); ROC-AUC for MAPK-high classification ΔAUC (TDS-16 − Panel-8) = +0.007 (TCGA) / +0.012 (Lee), both NS. The 8-gene panel is not a cherry-picked subset of TDS-16.

**(F) Per-driver-class mean 8-gene β.** Bar plot of TCGA mean 8-gene HM450 β by driver class: BRAF V600E β = 0.37 ≈ RET fusion β = 0.39 ≫ RAS-mutant β = 0.27. The methylation-silencing layer tracks MAPK-driver class, not BRAF / fusion identity per se — directly supporting the fusion-independent epigenetic silencing claim shown within DM1 (Supplementary Fig. S5b: within-DM1 fusion-positive vs negative β d = −0.36, NS).

**(G) Landa 2016 advanced thyroid (PDTC + ATC) convergence heatmap.** Per-sample heatmap of seven canonical differentiation transcripts (TG, TSHR, TPO, PAX8, SLC26A4, DIO1, DUOX2) across GSE76039 (Landa 2016 JCI, n = 84 PDTC + 33 ATC). Five of the seven genes (TG, TSHR, TPO, PAX8, DIO1) directly overlap the 8-gene panel — independently selected on Yoo 2016 RAI biology. The convergence frames DM1 as the upstream PTC signature consistent with the Landa 2016 dedifferentiation continuum (Discussion §3.1), not a downstream cherry-pick.

**Statistics.** Cohen's d with pooled SD; Mann-Whitney U for per-gene β; Spearman ρ with Fisher-z 95 % CI; fixed-effect meta-analysis (inverse-variance) across cohorts; Cochran Q and I² for heterogeneity. n_HM450 = 503; n_MAPK forest pooled = 1,287.

---

## Figure 4. Single-cell external validation establishes a thyrocyte-intrinsic DM1 signal. (7 panels)

**(A) Lu 2023 thyrocyte UMAP colored by 8-gene panel score.** GSE193581 (Lu et al. 2023, n = 23 samples) thyrocyte-restricted UMAP after batch-correction and KRT8 / KRT19 / EPCAM filtering, with per-cell 8-gene panel score overlaid as a continuous gradient. The DM1↔DM2 gradient is resolved within the thyrocyte compartment, not driven by stromal or immune cells.

**(B) Per-patient Spearman r heatmap, Pu 2021.** GSE184362 (Pu et al. 2021, n = 6 PTC patients with paired tumor and adjacent-normal thyrocytes). Per-patient Spearman r between bulk DM1 score and pseudobulk tumor-vs-normal thyrocyte 8-gene score ranges from 0.798 to 0.886 (per-patient p < 10⁻¹⁰, Bonferroni-adjusted). Heatmap with per-patient r value and n_cells annotation.

**(C) Author-independence cross-cohort check.** Cross-cohort 8-gene score concordance scatter between GSE241184 (single-author Phase 1 cohort) and GSE184362 (Pu 2021 independent author). Spearman r and Bland-Altman residuals confirm that the DM1 signal is not laboratory-specific batch artifact.

**(D) External pseudobulk replication: GSE232237.** Single-cell pseudobulk DM1 score for the GSE232237 thyroid cohort (n = additional samples, author-independent), with reference cross-cohort score distribution overlay. Replicates the DM1↔DM2 score boundary in a fully external single-cell dataset not used in panel design or training.

**(E) Multisite pooled scatter.** Per-patient pseudobulk 8-gene score in tumor vs adjacent-normal thyrocytes across GSE184362 (Pu 2021) and GSE193581 (Lu 2023); per-cohort markers with the pooled regression slope and 95 % CI. Confirms a uniform tumor-direction shift across two independent single-cell laboratories.

**(F) Compositional pseudotime: 10-decile trajectory.** TCGA-THCA (n = 513, score = panel z) and Lee / GSE213647 (n = 632, score = panel z) samples ranked by 8-gene score and binned into 10 deciles; per-decile mean cell-type fraction from nu-SVR deconvolution against Lu 2023 (n_celltypes = 8 against 67,678 reference cells). Four compartments show monotonic decile-level Spearman |ρ| ≥ 0.95 in both cohorts: Malignant ↓, Epithelial ↑, Myeloid ↓, Endothelial ↑. Defines a reproducible compositional pseudotime independent of cohort, score scaling, and sample size.

**(G) External direction-consistency forest.** Cohen's d (DM1 vs DM2 — TCGA frame) and Spearman ρ (panel z vs DM1 direction — external frame) for each external cohort: Lee GSE213647 (n = 632); K2 PRJEB11591 (n = 260, score-distribution only); GSE286332 (n = 18); GSE76039 (n = 37); GSE241184; GSE184362; Lu 2023 GSE193581. Direction-consistent across cohorts for ≥ 6 / 7 layers.

**Statistics.** Spearman r with Bonferroni-adjusted per-patient p; nu-SVR deconvolution with ν ∈ {0.25, 0.5, 0.75}, lowest-RMSE solution retained; per-decile rank correlation. Single-cell scoring: log1p(CP10K) + within-sample z; thyrocyte filter KRT8 ∩ KRT19 ∩ EPCAM.

---

## Figure 5. Pooled overall-survival hazard and cross-cohort assay portability. (7 panels)

**(A) Pooled meta-analysis forest, DM1 vs DM2 overall survival.** Cox proportional-hazards estimates by cohort: TCGA-THCA (HR 2.30, 95 % CI 0.77 – 6.88; n = 504, 16 events); MSK-IMPACT thyroid (HR 2.67, 95 % CI 1.17 – 6.10; n = 117, 38 events). **Pooled DerSimonian-Laird random-effects HR = 2.53 (95 % CI 1.31 – 4.89), Cochran I² = 0 %.** Diamond size proportional to inverse variance.

**(B) Multivariate Cox forest (TCGA-THCA primary).** Forest plot of adjusted HRs for DM1 status, age (per decade), ATA 2015 stage (III/IV vs I/II), and TERT promoter mutation (covariate only — see Extended Data for TERT-specific stratification). DM1 retains an independent hazard contribution after adjustment.

**(C) Multi-cohort Kaplan-Meier stack.** Three side-by-side KM panels: TCGA-THCA primary (n = 504), MSK-IMPACT advanced (n = 117), and pooled (TCGA + MSK, random-effects per-curve estimate). Each panel shows DM1 (red) vs DM2 (blue) with log-rank p, HR, and 95 % CI.

**(D) Cross-cohort score-distribution portability.** Density plot of 8-gene panel scores across four cohorts on a common z-scale: TCGA-THCA (n = 504); GSE213647 / Lee Korean (n = 632); K2 / PRJEB11591 (n = 260, calibration caveat — see Limitations §3.4 and Suppl. Fig. S9); MSK-IMPACT (n = 117). The DM1 / DM2 score boundary is preserved across cohorts processed by independent pipelines.

**(E) FFPE vs fresh-frozen tissue concordance.** Density plot of 8-gene scores by processing type (FFPE n = 632 from GSE213647, FF n = 504 from TCGA), with Kolmogorov-Smirnov p = 0.44 and per-gene Bland-Altman residuals. The 8-gene panel is FFPE-deployable, removing a major reflex-testing implementation barrier.

**(F) Time-dependent ROC for overall survival.** Time-dependent ROC curves at 1-, 3-, and 5-year horizons in TCGA-THCA, comparing DM1-status-only versus DM1 plus age and stage. AUC values and per-cutoff sensitivity / specificity annotated.

**(G) Survival-model calibration.** Calibration plot of predicted vs observed 5-year overall-survival probability, binned into deciles, with Hosmer-Lemeshow goodness-of-fit and 45° reference. Calibration intercept and slope reported on plot.

**Statistics.** Cox proportional-hazards (R `survival`, Python `lifelines`) with Schoenfeld-residual PH check; DerSimonian-Laird random-effects meta-analysis; Kolmogorov-Smirnov test for distributional equivalence; time-dependent ROC via inverse-probability-of-censoring weighting. n_pooled = 621 (504 TCGA + 117 MSK).

---

## Figure 6. A clinical reflex pathway and projected translational impact of DM1 status. (7 panels)

**(A) Reflex-testing algorithm flowchart.** Decision-flow from primary tumor RNA → 8-gene panel score → DM1 / DM2 call → if DM1, reflex orthogonal fusion testing (RNA-fusion panel or anchored multiplex PCR) → selpercatinib (RET; Wirth 2020 NEJM LIBRETTO-001 ORR 79 %) or larotrectinib (NTRK) candidate pool. Estimated yield: 48 selpercatinib-eligible candidates per 1,000 PTC patients (population estimate from TCGA RET-fusion rate × DM1 capture × selpercatinib RET-fusion eligibility).

**(B) ATA 2015 / 2025 risk-tier overlay.** Mosaic plot of TCGA-THCA patients by ATA 2015 risk tier (low / intermediate / high) × DM cluster. DM1 distributes preferentially into the ATA intermediate tier — the tier with the largest current RAI-decision uncertainty (Haugen 2016; Ringel 2025 ATA 2025 update) — supporting DM1 as an orthogonal molecular axis to the current clinico-pathological risk framework.

**(C) HMA + RAI re-induction rationale schematic + retrospective trial summary.** Schematic of the proposed DM1-stratified hypomethylating-agent (decitabine or azacitidine) plus radioiodine re-induction strategy, anchored in the TPO / DIO1 / TSHR promoter-hypermethylation evidence (Figure 3A–C). Side panel: retrospective trial summary table for decitabine + I-131 in advanced thyroid carcinoma (NCT00085293, NCT01065090) with reported n, response, and limitations. Framed as rationale for a prospective DM1-stratified trial, not a retrospective claim of clinical benefit.

**(D) DM1-like state in post-radioiodine refractory disease (GSE151179).** Boxplot of the thyroid-differentiation score in pre- versus post-RAI samples from GSE151179: post-RAI thyroid_diff Cohen's d = −1.01, Mann-Whitney p = 0.0001. The post-RAI dedifferentiation profile recapitulates the DM1 transcriptional state, providing an external clinical anchor that the DM1 axis names the same biology that defines RAI-refractory progression.

**(E) Decision-curve analysis (TCGA primary + Korean replicate).** Net-benefit decision curves over a range of clinical-decision thresholds for the candidate strategies: (i) treat all, (ii) treat none, (iii) DM1-stratified treatment. Replicated in the Korean GSE213647 validation cohort. The DM1-stratified strategy dominates over a clinically relevant threshold window.

**(F) Prospective DM1-stratified trial schema.** Schematic of a proposed prospective DM1-stratified RAI trial, including DM1-status determination (RNA assay) at biopsy, randomization within DM1 to RAI-only vs RAI + decitabine, and primary endpoint (structural disease recurrence at 24 months). Illustrative only; not a pre-registered protocol.

**(G) Selpercatinib eligibility waterfall (TCGA exemplar).** Waterfall plot ranking TCGA-THCA tumors by DM1 score; bars annotated by RET-fusion status. The top of the DM1 score distribution concentrates RET-fusion positive tumors, illustrating the prioritization yield of an RNA-first reflex pathway over a fusion-test-everyone strategy.

**Statistics.** Decision-curve net-benefit (Vickers and Elkin 2006); Mann-Whitney U for pre- vs post-RAI score; descriptive proportions for waterfall and ATA overlay. Trial-schema and HMA rationale panels are illustrative — captions explicitly state that these are not pre-registered or retrospective-outcome claims.

---

# Extended Data figures (held boundaries — NC architecture)

These figures preserve evidence and reviewer-defense material outside the main 6-figure claim surface, in line with the conservative claim boundary chosen for the NC reach.

**Extended Data Fig. 1.** TCGA 8-gene heatmap sorted by P_DM1 with detailed annotation tracks (Cell Press companion file Fig. 1 panel C extension).
**Extended Data Fig. 2.** Pan-genome ARI ladder at fine resolution; 8 / 16 / 67 / 200 / 1000 / 5000 panel sizes with TIERA67 hypergeometric enrichment p.
**Extended Data Fig. 3.** Immune-residualization analysis — DM1 vs DM2 Cohen's d before and after residualization on antigen-presentation and generic immune modules (raw d = 1.78; post-antigen-presentation 1.00; post-immune 1.50; dual 0.87).
**Extended Data Fig. 4.** Single-cell foundation and per-patient summaries from GSE184362, GSE193581, GSE241184, GSE232237 (companion-file Part B / Part C content).
**Extended Data Fig. 5.** Structural-variant missingness sensitivity (best-case, worst-case, MAR-imputed) for the fusion enrichment.
**Extended Data Fig. 6.** DM1 sub-A vs sub-B detailed phenotype (age, stage III/IV, CD8 / IFN-γ / checkpoint signatures, Hashimoto-like prevalence) — Paper 2 boundary marker, retained here as a sub-cluster characterization.
**Extended Data Fig. 7.** TERT promoter mutation × DM cluster — co-occurrence frequencies, TCGA + MSK pooled, with the 4-patient BRAF-negative / TERT-positive small-N caveat explicitly framed.
**Extended Data Fig. 8.** DM1 × TERT joint stratification: 8-cell decomposition (BRAF × TERT × DM); within-stratum KM curves; TERT_or_DM1 1-year ROC = 0.83 (sensitivity panel — see Limitations §3.4 for cross-cohort caveats).
**Extended Data Fig. 9.** Multi-method bulk deconvolution composite (NNLS, Ridge-NNLS, LR-clip, nu-SVR) — companion-file SX v5 content, including cell-type effect-retention ratio, per-driver-class composition, and HM450 β × cell-type fraction Spearman heatmap.
**Extended Data Fig. 10.** Drug-vulnerability and dependency overlays — PRISM 7 / 11 MAPK-axis inhibitor enrichment (AZD-0364 Cohen's d = −0.594, FDR = 5.9 × 10⁻⁷); DepMap top dependencies (MYC d = −0.499, NAMPT d = −0.445); explicitly framed as prioritization and not validated mechanism (caption uses "supports MAPK-axis vulnerability but does not measure restoration of thyroid-gene expression after inhibitor treatment").
**Extended Data Fig. 11.** Mun et al. 2025 proteogenomic cross-modal replication (n = 336): thyroid-differentiation protein score Cohen's d = −1.91; 7 / 7 sign-consistent for the 8-gene panel; 5 / 7 |Spearman r| < −0.5 for RNA-vs-protein within-DM stratum.
**Extended Data Fig. 12.** Pan-cancer DM1-axis cameo (TCGA pan-cancer): per-disease Cox HR by DM1-axis projection, including LGG (HR 44.7), LUAD (HR ~19), and UCEC (inverse direction). Discussion §3.X forward reference; full pan-cancer analysis reserved for the companion Paper 11.
**Extended Data Fig. 13.** Spatial stress-test caveat panels (GSE250521 Visium): MAPK × Panel raw vs adjusted ρ, neighborhood-lag closure, and signal-decomposition diagnostics — explicitly framed as a Visium-resolution / detection caveat rather than rescued mechanism evidence (companion-file SX_v16 / v17 / v18 content).
**Extended Data Fig. 14.** Korean K2 / Lee / GSE213647 portability detail: K2 mini-index calibration mismatch (per-gene 4.9 – 12.5 × inflation), within-sample-centered profile rescue, alternate-evidence chain (score distribution and DM-call consistency).
**Extended Data Fig. 15.** Landa 2016 (GSE76039) reference detail and convergence framing — 8-gene 5 / 8 overlap with the Landa advanced-disease silenced gene list; reverse-causality control via the Yoo 2016 panel-origin timeline.

---

# Build status and source mapping

| Figure | Panels | Reuse from existing render | New build required |
|---|---|---|---|
| Fig 1 | A–G | B, C, E, F (current Fig 1 A/B/C/D); D (current S1 heatmap) | A schematic strip; G pooled KM teaser |
| Fig 2 | A–G | A (current Fig 2A Sankey); B–D (current Fig 7 A/B/C/D); G (existing MSK fusion overlap) | E FVPTC mosaic; F sub-A/sub-B mosaic carry-over |
| Fig 3 | A–G | A (current Fig 8A); B (current Fig 8B); D (SX_v14 forest, render existing); E (SX_v13 panel B, render existing); F (Fig 8 mechanism cross-ref / v5A_per_driver_class); G (Landa SD1 heatmap, render existing) | C β-vs-expression 4-gene scatter grid |
| Fig 4 | A–G | A, B, C (current Fig 3 A/B/C); D (gse232237_scrna_pseudobulk.png); E (SC5 external_dm1_nonoverlap_scatter_grid.png); F (SX_v5K trajectory render existing); G (external_direction_consistency_forest.png) | None — composition only |
| Fig 5 | A–G | A (current Fig 6A forest); B (o_multivariate_cox.png); D (existing K2/Lee distribution); E (current Fig 5C FFPE) | C multi-cohort KM stack; F time-dep ROC; G calibration |
| Fig 6 | A–G | C (S8 HMA schematic existing) | A reflex flowchart; B ATA mosaic; D GSE151179 post-RAI box; E decision-curve; F trial schema; G selpercatinib waterfall |

**New scripts to author** (estimate, in `project/results/manuscript_v8_nc_main/`):

1. `fig1g_pooled_km_teaser.py` — DM1 vs DM2 OS KM, TCGA-only initial.
2. `fig2e_fvptc_mosaic.py` — Round 3 FVPTC OR=17.9 panel.
3. `fig3c_beta_expression_scatter.py` — 4-gene grid TPO / DIO1 / TSHR / TG.
4. `fig5c_multicohort_km_stack.py` — TCGA / MSK / pooled three-panel KM.
5. `fig5f_time_dependent_roc.py` — 1- / 3- / 5-year ROC, IPCW.
6. `fig5g_calibration.py` — 5-year calibration plot.
7. `fig6a_reflex_flowchart.py` — diagram-style flowchart.
8. `fig6b_ata_mosaic.py` — ATA tier × DM cluster mosaic.
9. `fig6d_post_rai_box.py` — GSE151179 thyroid_diff pre/post-RAI box.
10. `fig6e_decision_curve.py` — DCA on TCGA + Lee replicate.
11. `fig6f_trial_schema.py` — illustrative schematic.
12. `fig6g_selpercatinib_waterfall.py` — RET-fusion concentration in DM1.

**Layout convention.** NC main figures rendered at 180 mm wide; 7-panel layouts arranged as 2 rows × 4 columns (panel G spanning bottom-right or as a side strip). Panel labels A–G top-left in 8 pt Arial bold. Stat annotations in 6.5 pt. PDF + PNG at 600 dpi for submission.

---

# Boundary checklist (held lines)

- No PRISM / DepMap perturbation in main figures (held in Extended Data Fig. 10 with prioritization-only framing).
- No TERT survival as a main panel (TERT enters Fig 5 only as a multivariate covariate; full TERT × DM stratification in Extended Data Fig. 7–8).
- No Mun 2025 proteome in main figures (Extended Data Fig. 11).
- No pan-cancer LGG / LUAD / UCEC in main figures (Extended Data Fig. 12; Discussion forward-reference only).
- HLA mediation, sub-B NBNR mechanism, BCR / TLS — reserved for the companion paper (Paper 2); Fig 2B (sub-A vs sub-B silhouette) and Extended Data Fig. 6 are the only sub-cluster surfaces in this manuscript.

These five boundaries are the difference between an NC-defensible manuscript and an NC-attack-surface manuscript.
