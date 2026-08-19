---
title: "Nature Communications-style full manuscript draft"
date: 2026-07-08
author: Codex scaffold for Seungho Cook
target_journal: Nature Communications
status: full manuscript-form scaffold; protected voice slots preserved
source_context:
  - project/manuscript_v8/NATURE_CANCER_MANUSCRIPT_ARCHITECTURE_2026_07_08.md
  - project/manuscript_v8/04_results.md
  - project/manuscript_v8/05_figure_captions_NC.md
  - project/manuscript_v8/07_star_methods.md
  - project/manuscript_v8/DM1_MASTER_BRIEF_FOR_GPT_AND_MEETING_2026_06_25.md
protected_sections:
  - Introduction opening hook
  - Introduction final aim paragraph
  - Discussion 3.1 mechanism interpretation
  - Limitations paragraph block
  - Cover letter paragraph 1
  - Reviewer Q9 prose
---

# A thyroid-lineage state marks iodine-handling failure in papillary thyroid cancer

Seungho Cook^1,2,* and Hyeong-Won Yu^1,2,*

^1 Department of Internal Medicine, Seoul National University Bundang Hospital, Seongnam, Republic of Korea  
^2 Seoul National University College of Medicine, Seoul, Republic of Korea  

*Correspondence: kukshomr@gmail.com; [VERIFY institutional corresponding email]

## Abstract

Radioiodine remains central to the management of differentiated thyroid cancer, yet current molecular frameworks do not directly measure the thyroid-lineage machinery required for iodine uptake and organification. Here we define a compact eight-gene thyroid differentiation and iodine-handling axis comprising **TG, TPO, TSHR, SLC5A5, DIO1, PAX8, NKX2-1** and **FOXE1**, and apply it across bulk transcriptomic, methylation, proteomic, single-cell and clinical thyroid cancer cohorts. In TCGA-THCA, this axis resolves iodine-handling-low (DM1) and iodine-handling-high (DM2) states that are recovered by pan-genome clustering (adjusted Rand index 0.92), not by driver-anchor genes alone (adjusted Rand index -0.007). DM1 is enriched for kinase fusions (76.8%; odds ratio 7.41 versus DM2), captures 81.8% of RET-fusion-positive tumours, and shows promoter hypermethylation of thyroid differentiation machinery (TPO Cohen's d = 2.30; mean eight-gene beta 0.385 versus 0.253). The axis is reproduced across external expression cohorts, Korean cohorts, single-cell thyrocytes and proteomic data, and is associated with overall survival in a TCGA plus MSK meta-analysis (hazard ratio 2.53, 95% confidence interval 1.31-4.89). These findings position the eight-gene axis as a compact readout of thyroid-lineage collapse and a candidate framework for prospective radioiodine harm-avoidance studies.

## Introduction

**[AUTHOR HOOK SLOT - protected voice section. Insert the opening clinical hook here. Recommended direction: start from the clinical paradox that most differentiated thyroid cancers have excellent survival, yet a subset receives repeated radioiodine despite molecularly poor iodine-handling biology.]**

Papillary thyroid carcinoma is commonly interpreted through canonical driver classes, including BRAF-like and RAS-like molecular programmes. These frameworks have clarified tumour evolution and risk stratification, but they do not directly quantify whether a tumour retains the thyroid-cell functions required for radioiodine uptake, organification and hormone synthesis. This distinction is clinically important because radioiodine response is constrained by cellular differentiation state, not by driver identity alone.

The thyroid differentiation score and related lineage markers have previously shown that advanced thyroid cancers lose expression of canonical thyroid genes. However, most existing signatures are either broad transcriptomic scores or mechanistic descriptors rather than compact, clinically translatable readouts. A practical molecular state classifier should satisfy three conditions: it should be anchored in thyroid biology, reproduce across platforms and populations, and retain a clear boundary between retrospective association and prospective clinical utility.

We therefore treated iodine handling as a circuit rather than as a single gene. The panel includes five effector genes directly involved in thyroid hormone production and iodine handling (**TG, TPO, TSHR, SLC5A5/NIS, DIO1**) and three thyroid-lineage transcription factors (**PAX8, NKX2-1/TTF-1, FOXE1/TTF-2**). The premise is that a coordinated decrease across this circuit reports a tumour state in which radioiodine uptake and processing are biologically compromised.

**[AUTHOR AIM SLOT - protected voice section. Insert the final introduction paragraph here. Factual ingredients available: TCGA-THCA discovery, 8-gene axis, DM1/DM2, fusion/MAPK biology, HM450 methylation, external validation, RAI harm avoidance.]**

## Results

### A compact thyroid-lineage axis resolves iodine-handling states in TCGA-THCA

We first assembled a thyroid-relevant candidate space designed to separate lineage-state information from canonical driver identity. The 67-gene TIERA67 candidate pool included thyroid differentiation, iodine-handling, MAPK-output, driver-anchor, aggressive-disease, dedifferentiation and immune-context categories. From this pool, we selected an eight-gene thyroid differentiation and iodine-handling panel: **SLC5A5, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1** and **DIO1** (Fig. 1a,b). The panel was chosen on thyroid biology rather than on genome-wide optimization, with the driver-anchor category retained only as a leakage-control comparator.

Unsupervised clustering of TCGA-THCA primary tumours on the eight-gene transcript space resolved two states, which we refer to operationally as DM1 and DM2 (Fig. 1c,d). DM1 represented an iodine-handling-low state, whereas DM2 retained higher thyroid differentiation and iodine-handling expression. In TCGA-THCA, DM1 accounted for 28.4% of primary tumours.

We next asked whether the eight-gene panel artificially imposed a cluster boundary or instead captured a broader transcriptomic state. Pan-genome clustering using the top 5000 genes by median absolute deviation reproduced the eight-gene partition with an adjusted Rand index (ARI) of 0.92, comparable to the full TIERA67 pool (ARI = 0.90). In contrast, driver-anchor genes alone produced no meaningful recovery of the same structure (ARI = -0.007; Fig. 1e). The axis was also not explained by driver transcript abundance: BRAF expression was indistinguishable between BRAF V600E and wild-type tumours (Cohen's d = -0.044, Mann-Whitney p = 0.57), and single-feature AUCs for BRAF, TERT, KRAS, NRAS and HRAS were near chance (Fig. 1f).

These analyses indicate that the eight-gene panel is not merely a small classifier trained to recover driver status. Rather, it is a compact readout of a broader thyroid-lineage programme that remains visible when the transcriptome is considered without panel restriction.

### The iodine-handling-low state recovers clinically relevant dark matter

Canonical driver frameworks leave a clinically important BRAF- and RAS-negative compartment incompletely stratified. Applying the eight-gene axis to this dark-matter space recovered 131 of 180 BRAF/TERT-negative TCGA tumours into a defined DM1 or DM2 state, corresponding to 73% molecular recovery of the previously unresolved compartment (Fig. 2a).

The axis also tracked population-level differences relevant to thyroid cancer epidemiology. DM1 prevalence was 28.4% in TCGA-THCA and increased to 37.8% in Korean cohorts, consistent with the higher burden of driver-negative or non-BRAF-like thyroid cancer in East Asian populations (Fig. 2b). This enrichment supports the clinical relevance of an axis that is not anchored to BRAF status alone.

We then evaluated outcome association. In TCGA-THCA, DM1 was associated with higher overall-survival hazard than DM2, although the low event rate in primary thyroid cancer limited precision. We therefore combined TCGA-THCA with the MSK-IMPACT thyroid cohort, an advanced-disease-enriched cohort with higher event density. The TCGA hazard ratio was 2.30 (95% CI 0.77-6.88), and the MSK hazard ratio was 2.67 (95% CI 1.17-6.10). Random-effects meta-analysis yielded a pooled hazard ratio of 2.53 (95% CI 1.31-4.89) with no detectable between-cohort heterogeneity (I^2 = 0%; Fig. 5a). These data support DM1 as a retrospectively aggressive state while preserving the boundary that prospective outcome and treatment-response validation remain required.

### DM1 is enriched for kinase-fusion biology

To define the genetic context of DM1, we integrated TCGA-THCA structural-variant annotations from cBioPortal. Among tumours with usable structural-variant data, DM1 was 76.8% kinase-fusion positive (63 of 82), compared with 30.9% in DM2 (Fisher odds ratio 7.41, 95% CI 4.38-12.55, p = 1.9 x 10^-13; Fig. 2c). The fusion spectrum included RET, NTRK, ALK and BRAF rearrangements, with RET fusions forming the largest class (Fig. 2d).

This enrichment was robust to missing structural-variant status. Structural-variant missingness was independent of DM status (chi-square p = 0.56), and best-case, worst-case and matched-imputation sensitivity analyses preserved the direction and magnitude of the fusion enrichment. Independent support was observed in MSK-IMPACT thyroid cancers, where RET, ALK and PAX8-PPARG fusions were concentrated among DM1-called tumours despite the more limited structural-variant ascertainment in that cohort.

The fusion enrichment has immediate translational logic. Of 33 TCGA RET-fusion-positive tumours, 27 were classified as DM1, corresponding to an 81.8% RET-fusion capture rate (Fig. 2e). Thus, an iodine-handling-low RNA state could serve as an upstream trigger for reflex fusion testing, especially in settings where universal comprehensive fusion testing is not yet routine.

### DM1 carries promoter methylation of thyroid differentiation machinery

The high fusion rate in DM1 does not by itself explain the coordinated repression of iodine-handling and thyroid-lineage genes. We therefore examined Illumina HumanMethylation450 promoter methylation in TCGA-THCA tumours with paired methylation and DM calls (n = 503). DM1 tumours showed higher mean eight-gene promoter beta values than DM2 tumours (0.385 versus 0.253), corresponding to a 52% increase in mean panel methylation (Fig. 3a,b).

Per-gene methylation differences were strongest for thyroid hormone biosynthesis and differentiation genes. TPO showed the largest effect (Cohen's d = 2.30, p = 1.9 x 10^-18), followed by DIO1 (d = 1.24), TSHR (d = 1.20), PAX8 (d = 0.97), TG (d = 0.86), FOXE1 (d = 0.84) and NKX2-1 (d = 0.63). SLC5A5/NIS was the exception (d = 0.22, p = 0.42), consistent with the known complexity of NIS regulation beyond promoter methylation alone.

We also tested whether the eight-gene axis is a cherry-picked subset of canonical thyroid differentiation biology. Across TCGA and Lee/GSE213647, MAPK-output scores were inversely associated with both the deployable eight-gene panel and the canonical Yoo TDS-16. In five-cohort meta-analysis, the pooled MAPK-output versus Panel-8 Spearman correlation was -0.327 (95% CI -0.376 to -0.278; n = 1,287). BRAF V600E versus RAS-mutant differentiation effects were essentially identical for the eight-gene panel and TDS-16 (Cohen's d = -1.615 versus -1.616), and the TDS-16 added only minimal, non-significant AUC increments for MAPK-high classification (Fig. 3d,e).

Finally, the panel converged with an independent advanced-thyroid-cancer lineage-loss programme. Five of the eight panel genes overlapped the Landa 2016 PDTC/ATC differentiation-suppressed gene set (**TG, TSHR, TPO, PAX8, DIO1**), although the panel itself was selected from thyroid iodine-handling biology rather than from advanced-disease expression data (Fig. 3g). This convergence supports a shared biological axis of thyroid differentiation collapse.

### External validation establishes a multi-platform thyroid-lineage state

We next evaluated whether the DM1/DM2 axis generalized beyond TCGA. Across 19 external validation entries spanning expression, methylation, proteomic and single-cell contexts, the master cross-cohort forest showed a mean Cohen's d of 2.81 and median Cohen's d of 2.37. A per-gene by cohort by contrast matrix showed 80 of 80 cells moving in the expected direction, arguing against dependence on any single gene, cohort or platform.

In the Lee 2024 Korean PTC cohort (GSE213647, n = 632), within-cohort DM1/DM2 separation was strong (Cohen's d = 5.93). In four GPL570 microarray cohorts, the axis reproduced in the expected direction in all cohorts, with Spearman correlations <= -0.84. At the protein level, the Mun 2025 proteogenomic thyroid cohort (n = 336) showed 7 of 7 available panel proteins direction-consistent with the RNA-defined thyroid differentiation state.

Single-cell analyses supported a thyrocyte-intrinsic interpretation. In GSE184362 (Pu et al. 2021), patient-matched tumour and adjacent-normal thyrocyte scores showed per-patient Spearman correlations of 0.798 to 0.886, each with Bonferroni-adjusted p < 10^-10. In the Lu 2023 single-cell cohort (GSE193581), the eight-gene gradient remained visible after restricting to KRT8/KRT19/EPCAM-positive thyrocyte-lineage cells, arguing against stromal or immune contamination as the primary source of the signal (Fig. 4a,b).

The axis was also compatible with clinical tissue processing. FFPE and fresh-frozen score distributions were not detectably shifted (Kolmogorov-Smirnov p = 0.44), supporting the feasibility of retrospective pathology-specimen deployment and prospective clinical assay development (Fig. 5e).

### A translational pathway for radioiodine harm avoidance

The most defensible clinical position for the eight-gene axis is not immediate treatment selection, but harm avoidance. In this framing, an iodine-handling-low tumour state identifies patients in whom repeated high-dose RAI may be biologically less plausible and in whom additional molecular work-up may be valuable.

The axis aligns with external radioiodine-refractory biology. In GSE151179, post-RAI refractory tumours shifted toward a DM1-like thyroid-differentiation state (Cohen's d approximately -1.0; Mann-Whitney p approximately 10^-4), indicating that the DM1 programme resembles the transcriptional state observed after clinical RAI failure (Fig. 6d). Because this comparison is retrospective and public-cohort based, we interpret it as a biological anchor rather than as a validated response-prediction model.

A practical clinical pathway would begin with RNA-based assessment of the eight-gene panel in the primary tumour or archival FFPE tissue. A DM1 call would trigger reflex orthogonal fusion testing for RET, NTRK, ALK and BRAF rearrangements and would identify candidates for prospective evaluation of epigenetic re-differentiation strategies. The TSO500-only subset is insufficient as a standalone approximation (3-gene AUC 0.76), whereas a compact add-on including TPO and DIO1 improves performance toward the full panel (AUC 0.940). These observations inform clinical assay design but do not replace prospective validation.

## Discussion

### 3.1 Mechanistic interpretation

**[AUTHOR DISCUSSION 3.1 SLOT - protected voice section. Do not replace with AI prose. Suggested factual ingredients: DM1 as iodine-handling-low lineage state; fusion/MAPK enrichment; HM450 promoter methylation; TDS-16 equivalence; Landa 2016 convergence; "motivates rather than confirms" causal boundary.]**

### Clinical implications

The immediate clinical value of the eight-gene axis lies in its ability to report thyroid-lineage function using a compact assay surface. Current management decisions often integrate histology, stage, ATA risk tier and selected driver mutations, but these variables do not directly measure whether iodine-handling machinery remains transcriptionally intact. The DM1/DM2 framework adds a functional layer: whether the tumour retains the lineage programme required for radioiodine uptake and organification.

This distinction matters because differentiated thyroid cancer is often long-lived. For many patients, the harm of ineffective repeated high-dose RAI is not only treatment delay but cumulative toxicity. A lineage-state assay that identifies iodine-handling-low tumours could support earlier reflex molecular testing and more deliberate discussion of RAI benefit, especially in the ATA intermediate-risk setting where treatment decisions are frequently uncertain.

The strongest near-term application is therefore a reflex-testing pathway rather than a direct treatment-selection biomarker. DM1 positivity enriches for actionable kinase fusions and captures most RET-fusion-positive tumours in TCGA-THCA. This supports the concept of RNA-state-guided orthogonal fusion testing, particularly where comprehensive fusion testing is unavailable or resource-constrained.

### Generalizability and translational readiness

The axis is supported by several forms of generalizability: pan-genome recovery within TCGA, external Korean expression cohorts, microarray cohorts, proteomic directionality, single-cell thyrocyte restriction and FFPE compatibility. This breadth is important because a compact panel is only useful if it remains stable across preprocessing, ancestry, platform and tissue-processing conditions.

Nevertheless, the final clinical step remains prospective. Public cohorts can establish state biology, cross-platform reproducibility and retrospective association, but they cannot determine whether treatment decisions based on the state improve outcomes. The next validation layer should combine institutional FFPE or RNA testing, annotated RAI dose and response, toxicity endpoints, recurrence, structural disease progression and orthogonal IHC assessment of key lineage proteins.

### Limitations

**[AUTHOR LIMITATIONS SLOT - protected voice section. Suggested factual ingredients: retrospective public cohorts; TCGA low event rate; MSK advanced-disease enrichment; methylation/MAPK correlative not causal; n=19 fusion-negative DM1 methylation comparison is underpowered; K2 calibration mismatch; spatial Visium stress test non-confirmatory; prospective SNUBH/IHC validation pending.]**

## Methods

### Cohorts and data sources

TCGA-THCA RNA-seq, clinical, mutation, structural-variant and HM450 methylation data were obtained from TCGA and cBioPortal. The MSK-IMPACT thyroid cohort was used as an advanced-disease validation cohort. Korean validation cohorts included PRJEB11591/K2 and GSE213647/Lee. Single-cell validation used GSE184362/Pu 2021 and GSE193581/Lu 2023, with additional external single-cell support from GSE241184 and GSE232237 where indicated. Cohort-level sample sizes are reported in each figure caption and Supplementary Table 1.

### Eight-gene panel definition

The eight-gene panel comprised **TG, TPO, TSHR, SLC5A5, DIO1, PAX8, NKX2-1** and **FOXE1**. These genes were selected from the thyroid differentiation and iodine-handling component of a 67-gene thyroid-relevant candidate pool. Driver-anchor genes were retained as controls but were not used to define the final panel. Analyses used official HGNC symbols; clinical aliases were reported as SLC5A5/NIS, NKX2-1/TTF-1 and FOXE1/TTF-2.

### Score construction and clustering

Expression values were log-transformed and z-standardized within cohort. KMeans clustering with k = 2 was applied to the eight-dimensional panel expression matrix in TCGA-THCA to define DM1 and DM2 states. For cross-cohort transfer, panel scores were computed as within-cohort z-mean summaries, with additional within-sample-centered profile scoring for cohorts requiring calibration adjustment. Cluster reproducibility was evaluated using adjusted Rand index against panel-defined DM1/DM2 labels.

### Structural-variant analysis

Structural variants were retrieved from cBioPortal for TCGA-THCA and MSK-IMPACT. Fusion classes included RET, NTRK, ALK, BRAF and PAX8-PPARG rearrangements. Fusion enrichment was tested using Fisher exact tests. Structural-variant missingness was assessed by chi-square testing and sensitivity analyses under best-case, worst-case and matched-imputation scenarios.

### Methylation analysis

HM450 promoter beta values were retrieved for TCGA-THCA tumours with paired DM calls. Per-gene and mean eight-gene promoter beta values were compared between DM1 and DM2 using Cohen's d and Mann-Whitney U tests. MAPK-output scores were computed using DUSP4/5/6, SPRY2/4, ETV4/5, PHLDA1 and CCND1. Cross-cohort MAPK-output versus thyroid-score correlations were pooled using Fisher-z transformation.

### Single-cell processing

Single-cell datasets were normalized using log1p(CP10K) and scored using within-sample z-standardized panel expression. Thyrocyte-lineage cells were defined using KRT8, KRT19 and EPCAM positivity when raw annotations permitted. Patient-level pseudo-bulk scores were generated by averaging panel expression within thyrocyte-lineage cells. Per-patient correlations and cross-cohort score concordance were evaluated using Spearman correlation.

### Survival analysis

Overall-survival analyses used Cox proportional-hazards models and Kaplan-Meier curves. TCGA-THCA and MSK-IMPACT hazard ratios were pooled using DerSimonian-Laird random-effects meta-analysis. Between-cohort heterogeneity was quantified using Cochran Q and I^2. Multivariable models included age, stage and selected driver covariates where sample size allowed.

### Statistical analysis

Cohen's d was computed using pooled standard deviation. Continuous variables were compared using two-sided Mann-Whitney U tests unless otherwise specified. Proportions were compared using Fisher exact tests. Multiple-testing correction used Benjamini-Hochberg false-discovery-rate control where applicable. Confidence intervals for correlations were estimated using Fisher-z transformation.

## Data availability

All source datasets are publicly available from TCGA, cBioPortal, GEO or ENA under the accessions listed in Supplementary Table 1. Processed per-sample score tables, methylation summaries, structural-variant joins and meta-analysis inputs will be deposited with the public code release at submission.

## Code availability

Analysis and figure-generation scripts are maintained in the project repository and will be archived with a DOI at submission. Current working scripts are located under `project/results/manuscript_v8_nc_main/`, `project/results/p_deconv_2026_05_08/` and related manuscript analysis directories.

## Acknowledgements

[AUTHOR SLOT - add funding, institutional support, collaborator acknowledgements and data-source acknowledgements.]

## Author contributions

S.C. conceived the analysis, curated datasets, performed computational analyses, generated figures and wrote the manuscript draft. H.-W.Y. supervised the clinical framing and interpretation. [VERIFY AND EXPAND BEFORE SUBMISSION.]

## Competing interests

The authors declare no competing interests. [VERIFY BEFORE SUBMISSION.]

## Figure legends

### Figure 1. The eight-gene thyroid-lineage axis resolves iodine-handling states.

(a) Study overview and cohort inventory. (b) Curation of the eight-gene panel from thyroid differentiation and iodine-handling biology. (c) UMAP of TCGA-THCA tumours in eight-gene transcript space, coloured by DM1 and DM2 state. (d) Heatmap of panel genes ordered by DM1 probability. (e) ARI ladder comparing panel, TIERA67, pan-genome and driver-anchor clustering. (f) Driver mRNA neutrality and single-feature classification AUCs. (g) TCGA overall-survival Kaplan-Meier curve by DM state.

### Figure 2. DM1 is enriched for kinase-fusion biology.

(a) Recovery of Xing dark-matter tumours into DM1/DM2 states. (b) DM1 prevalence across TCGA and Korean cohorts. (c) Fusion enrichment by DM state. (d) Fusion partner spectrum within DM1. (e) DM1 capture of RET-fusion-positive tumours. (f) Histological enrichment of follicular-variant PTC. (g) MSK-IMPACT fusion replication.

### Figure 3. DM1 shows epigenetic silencing of thyroid differentiation machinery.

(a) HM450 promoter beta heatmap for panel genes and TDS controls. (b) Mean eight-gene beta by DM state. (c) Methylation-expression coupling for representative thyroid genes. (d) Cross-cohort MAPK-output versus Panel-8 score forest. (e) Functional equivalence of Panel-8 and TDS-16. (f) Per-driver-class mean eight-gene beta. (g) Landa 2016 advanced-thyroid-cancer convergence heatmap.

### Figure 4. Single-cell validation establishes a thyrocyte-intrinsic state.

(a) Lu 2023 thyrocyte UMAP coloured by eight-gene score. (b) Pu 2021 per-patient tumour and adjacent-normal thyrocyte correlations. (c) Author-independence cross-cohort check. (d) External single-cell pseudo-bulk replication. (e) Multisite pooled tumour versus normal pseudo-bulk scatter. (f) Compositional pseudotime across TCGA and Lee cohorts. (g) External direction-consistency forest.

### Figure 5. DM1 is associated with survival and assay portability.

(a) TCGA plus MSK random-effects overall-survival meta-analysis. (b) Multivariable Cox model. (c) Multi-cohort Kaplan-Meier stack. (d) Cross-cohort score-distribution portability. (e) FFPE versus fresh-frozen concordance. (f) Time-dependent ROC. (g) Survival-model calibration.

### Figure 6. Translational pathway for radioiodine harm avoidance.

(a) RNA-state-guided reflex fusion-testing algorithm. (b) ATA risk-tier overlay. (c) Hypomethylating-agent plus RAI re-induction rationale. (d) DM1-like state in post-RAI refractory disease. (e) Decision-curve analysis. (f) Prospective DM1-stratified trial schema. (g) RET-fusion concentration across the DM1 score distribution.

