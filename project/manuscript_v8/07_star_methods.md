---
title: "Paper 1 manuscript v8 — STAR Methods (draft v1)"
date: 2026-05-04
author: Seungho Cook
target_format: Cell Press STAR Methods (unlimited length, reproducibility-focused)
status: clean draft
---

# STAR Methods (draft v1)

---

## Key Resources Table

(structured Cell Press format — populated below as plain table; final manuscript: per-section table)

| Reagent / Resource | Source | Identifier |
|---|---|---|
| **Biological samples — bulk** | | |
| TCGA-THCA RNA-seq + WGS + HM450 + clinical | The Cancer Genome Atlas | dbGaP phs000178; cBioPortal `thca_tcga_pub` |
| MSK-IMPACT thyroid (Landa 2016 cohort) | Landa et al., 2016 | cBioPortal `thca_mskcc_2016` |
| K2 / PRJEB11591 Korean PTC RNA-seq (Yoo 2016) | Yoo et al., 2016 | ENA PRJEB11591 |
| Lee Korean PTC cohort | Lee et al., GEO | GSE213647 |
| GSE286332 Korean PTC reference arm | Macrogen / Dongguk Univ | GEO GSE286332 |
| GSE184362 Pu 2021 single-cell PTC | Pu et al., 2021 | GEO GSE184362 |
| GSE193581 Lu 2023 single-cell PTC | Lu et al., 2023 | GEO GSE193581 |
| GSE241184 Phase 1 single-cell PTC | (Phase 1) | GEO GSE241184 |
| **Software** | | |
| pyDESeq2 (DEG analysis) | Snakemake/lab | https://github.com/owkin/PyDESeq2 |
| scikit-learn (KMeans, LogReg, AUC) | Pedregosa et al. | https://scikit-learn.org |
| lifelines (Cox PH, log-rank) | Davidson-Pilon | https://lifelines.readthedocs.io |
| STAR aligner | Dobin et al. | https://github.com/alexdobin/STAR |
| kallisto (pseudo-alignment) | Bray et al. | https://pachterlab.github.io/kallisto |
| GENCODE v44 reference | EBI | https://www.gencodegenes.org/human/release_44.html |
| **Deposited data** | | |
| cBioPortal SV (RET/NTRK/ALK/BRAF fusions) | cBioPortal API | `thca_tcga_pub` study, SV endpoint |
| HM450 promoter methylation (Illumina HumanMethylation450) | cBioPortal API | `thca_tcga` legacy study |
| **Source code (this paper)** | | |
| 8-gene panel + DM cluster pipeline | Cook et al., this paper | [TODO: insert public code repository URL on submission]; [TODO: insert Zenodo DOI on submission] |

---

## Resource availability

**Lead contact.** Further information and requests for resources should be directed to the lead contact, Seungho Cook (kukshomr@gmail.com).

**Materials availability.** This study did not generate new unique reagents. All analyses were performed on publicly accessible datasets (TCGA, GEO, ENA, cBioPortal). Korean cohort access (K2 / PRJEB11591, GSE213647, GSE286332 reference arm) is available via the listed repositories.

**Data and code availability.** Public source data are available from TCGA, GEO, ENA, and cBioPortal under the identifiers listed above. Analysis code and figure-generation scripts will be released in a public repository together with an archival DOI at submission or acceptance ([TODO: insert public code repository URL on submission]; [TODO: insert Zenodo DOI on submission]). Intermediate data tables used in the manuscript, including per-sample DM scores, fusion annotations, methylation summaries, and meta-analysis inputs, are provided through Supplementary Tables S1-S10.

---

## Experimental model and study participant details

This study uses publicly available genomic and transcriptomic data from previously published cohorts:

- **TCGA-THCA** (n = 504 with overall survival annotation; 513 primary tumors total). Discovery cohort. Publicly available via The Cancer Genome Atlas (Cancer Genome Atlas Research Network, 2014).
- **MSK-IMPACT thyroid** (n = 117; advanced disease, mostly PDTC + ATC). Validation cohort. (Landa et al., 2016).
- **K2 / PRJEB11591** (n = 235; primary Korean PTC, post-QC). Validation cohort. (Yoo et al., 2016).
- **Lee / GSE213647** (n = 630; Korean PTC, post-QC). Validation cohort.
- **GSE286332 reference arm** (n = 9 Korean PTC). Small external Korean reference set retained here for calibration and score-portability description only; <em>dropped from the Paper 1 Korean validation total per the 2026-05-07 P1-2 audit. The canonical Paper 1 Korean validation cohort = K2 (n = 235) + Lee et al. (n = 630) = n = 865. GSE286332 PTC samples are retained in Paper 2 scope (GSE286332 PTC vs PTC+HT main cohort, n = 18).</em>
- **GSE184362 Pu 2021** (n = 7 PTC patients; single-cell). External validation.
- **GSE193581 Lu 2023** (n = 23 single-cell samples). External validation.
- **GSE241184** (n = 1; Phase 1 single-cell). Internal pilot.

All studies were originally approved by the respective institutional review boards. The present analysis used only de-identified public data and is exempt from additional IRB review.

---

## Method details

### Cohort assembly and clinical metadata harmonization

Bulk RNA-seq quantifications were obtained as log2(TPM + 1) (TCGA) and log2(FPKM + 1) (Korean cohorts and the GSE286332 reference arm) and processed through unified gene-level filtering (≥10 reads in ≥30% of samples; GENCODE v44 protein-coding annotation). Clinical metadata (age, sex, stage, vital status, time-to-event) were harmonized from cBioPortal and source publications. Per-cohort missingness was tabulated (Supplementary Table S2) and addressed in sensitivity analyses.

### 8-gene panel selection

The 67-gene candidate pool (TIERA67) was assembled from seven thyroid-relevant biological categories (Supplementary Table S1): TDS-core differentiation markers (16 genes), MAPK-output transcripts (10), thyroid driver genes (12, including BRAF, NRAS, HRAS, KRAS, RET, NTRK1, NTRK3, ALK, PAX8, PPARG, TERT, EIF1AX), aggressive-disease markers (10), dedifferentiation/EMT markers (10), light immune-stromal markers (5), and thyroid-lineage extras (4). The 8-gene panel (SLC5A5/NIS, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, DIO1) was selected from the TDS-core sub-category based on canonical RAI-uptake biology (Yoo et al., 2016 *PLOS Genet*) — independently of and prior to access of the Landa et al. (2016) ATC-silenced gene list. This panel design therefore predates exposure to advanced-disease ATC transcriptomic outputs and represents an independent path to the same differentiation axis.

Panel size sensitivity was tested across n = 8, 10, 12, and 16 gene variants (Supplementary Figure S2). The 16-gene full TDS-core yielded TCGA 5-fold AUC of 0.975 versus the 8-gene panel's 0.962 (ΔAUC = 0.013, NS); the 10-gene and 12-gene intermediate panels yielded 0.972 and 0.969, respectively. The 8-gene panel was retained for clinical interpretability.

### DM1/DM2 cluster definition

Bulk RNA expression of the 8-gene panel was z-standardized within cohort, and KMeans (k = 2, scikit-learn default initialization) was applied to assign DM1 versus DM2 cluster identity. The cluster boundary was confirmed by silhouette analysis (mean silhouette score > 0.4 within both clusters). For per-sample classification probability (P(DM1)), a logistic regression classifier was trained on the TCGA DM1/DM2 partition using the within-sample-centered profile of the 8-gene panel (mean panel score subtracted), yielding cohort-portable class probabilities (cross-validated AUC = 0.968 in TCGA; Methods, Supplementary Figure S3).

### Single-cell external validation

For GSE184362 (Pu et al., 2021), per-patient pseudo-bulk DM1 scores were computed by averaging thyrocyte-marker-positive (KRT8, KRT19, EPCAM) cell profiles per patient and per condition (tumor versus adjacent normal). Per-patient Spearman correlation between tumor and normal scores was computed. Lu 2023 (GSE193581) was processed with the same DM1-scoring framework together with stromal and immune contamination control.

### Survival analysis and meta-analysis

Cox proportional hazards regression was performed using `lifelines` (Davidson-Pilon, Python). Pooled meta-analysis was performed using random-effects DerSimonian-Laird estimation; between-cohort heterogeneity was assessed via Cochran I². Sub-cluster Cox analyses (DM1 sub-A vs sub-B) used patient-level event data within TCGA-THCA primary tumors (Supplementary Table S8).

### cBioPortal API access (SV + methylation)

Structural variants (RET, NTRK1/3, ALK, BRAF, PAX8-PPARG, others) were retrieved from cBioPortal `thca_tcga_pub` study via REST API (endpoint `/structural-variant/fetch`); SV-tested status was captured per tumor (n = 542 / 557 = 97.3%). HM450 promoter β-values for the 8-gene panel + DIO2 + SLC26A4 were retrieved from cBioPortal `thca_tcga` legacy study (n = 503 with HM450 + DM call); per-gene aggregate β values were computed from cBioPortal merged probe-per-gene format.

### Korean cohort processing

Yoo 2016 K2 (PRJEB11591): kallisto pseudo-alignment to GENCODE v44 transcriptome, transcript-to-gene aggregation, log2(TPM + 1) normalization. The 8-gene mini-index calibration mismatch (R4-4) was identified and addressed via within-sample-centered profile classification (Memory cross-ref: `v17_korean_k2_calibration.md`). Lee (GSE213647) and the GSE286332 reference arm used pre-computed FPKM tables for score-portability checks.

### Software and statistical environment

All analyses were performed in Python 3.11 with scikit-learn 1.5, pandas 2.2, numpy 1.26, scipy 1.13, lifelines 0.27, and pyDESeq2 0.4. Exact environment specifications and figure scripts will accompany the public repository release.

---

## Quantification and statistical analysis

### Cohen's d (pooled SD)

For all between-cluster comparisons, Cohen's d was computed using pooled standard deviation:
d = (μ₁ − μ₂) / σ_pooled, where σ_pooled = √[((n₁ − 1)·σ₁² + (n₂ − 1)·σ₂²) / (n₁ + n₂ − 2)].

### Mann-Whitney U and Fisher exact

Continuous distributions were compared via two-sided Mann-Whitney U (`scipy.stats.mannwhitneyu`). Discrete proportions were compared via two-sided Fisher exact test (`scipy.stats.fisher_exact`).

### Multiple testing correction

For per-gene tests across the 8-gene panel and Supplementary Tables, Benjamini-Hochberg false-discovery-rate (FDR) correction was applied at q = 0.05.

### Bootstrap confidence intervals

For mediation analysis (R²-based) and hazard ratio meta-analysis, 5,000 bootstrap iterations with patient-level resampling were used to estimate 95% confidence intervals.

### Immune-residualization analysis

To test whether DM1 represented a generic immune-infiltration artifact, the 8-gene panel score was residualized against predefined stromal and immune covariates, and residualized effect sizes were compared with the raw DM1-versus-DM2 contrast. This analysis was used only as a specificity check and not as a primary discovery endpoint.

### Bulk cell-type deconvolution (multi-method)

For per-cell-type granularity beyond the immune-residualization analysis, bulk RNA-seq cell-type fractions for TCGA-THCA primary tumors (n = 572) were estimated against pseudobulk profiles built from the Lu 2023 single-cell reference (GSE193581, `author_celltype` annotations, 67,678 cells, eight cell types: B / Endothelial / Epithelial / Fibroblast / Malignant / Myeloid / NK / T; Lu et al., 2023; processed `.h5ad`). The reference matrix was the per-cell-type mean of expression values restricted to the 1,898 highly variable genes that overlap between Lu 2023 and TCGA gene symbols (95% retention of HVG). Four deconvolution methods were compared: (i) non-negative least squares (NNLS, `scipy.optimize.nnls`); (ii) ridge-regularized NNLS (L2 penalty α=1, augmented design matrix); (iii) ordinary least squares with non-negative clipping and renormalization; and (iv) nu-SVR (CIBERSORT-style; `sklearn.svm.NuSVR(kernel='linear')` with three ν values 0.25/0.5/0.75, lowest-RMSE selection, weights clipped to non-negative and renormalized to sum-to-one). Per-method per-sample weights were normalized to sum-to-one cell-type fractions. The canonical 8-gene RAI score (`rai_score_recalc` from the v17p3 `A2_dm_score_full_cohort.tsv` table; n = 513 with `dm_like` ∈ {DM1_like, DM2_like}, 403 DM1 / 110 DM2) was regressed (OLS) on cell-type fraction subsets — stromal (Endothelial + Fibroblast); immune (T + Myeloid + B + NK); Epithelial-only (purity-confound proxy); all eight fractions — and DM1 vs DM2 Cohen's d was recomputed on residuals. Effect retention after full residualization on all eight cell-type fractions was 24% (NNLS), 24% (Ridge-NNLS), 70% (LR-clip), and 47% (nu-SVR); the nu-SVR retention is methodologically aligned with the canonical immune-residualization analysis (56% retention, see above). T-cell fraction collapse to zero in NNLS-family methods is a known sparse-weight artifact under correlated immune signatures and was resolved by nu-SVR (T cell mean fraction 11.3%, DM1 vs DM2 d = −0.84). Direction-consistent across all four methods: DM1 enriched for Malignant cell (d ≈ +1.1 to +1.4) and Myeloid cell (d ≈ +0.9 to +1.1) fractions, depleted for Epithelial cell (d ≈ −1.1 to −2.3) and Endothelial cell fractions (Supplementary Figure SX). Pseudobulk source code, per-method fraction tables, and the residualization grid are provided at `project/results/p_deconv_2026_05_08/`.

**Per-driver-class composition.** TCGA nu-SVR cell-type fractions were merged with the cBioPortal `v3_anchor_6class` driver call (BRAF V600E n = 280, RAS-mutant n = 54, driver-negative n = 179; `fusion_calls_per_sample.tsv`) and Cohen's d was computed per cell type for each pairwise contrast. RAS-mutant tumors retain Epithelial-cluster identity (d_RAS−BRAF = +1.52) while BRAF V600E tumors are Malignant-cell rich (d = −1.66) (Supplementary Figure SX panel H).

**Methylation × composition.** Per-sample HM450 mean 8-gene β values (`r5_2_sample_methylation_8gene.tsv`, TCGA-THCA n ≈ 484) were correlated (Spearman) with per-sample cell-type fractions for each cell type. The per-gene 8-gene heatmap reports gene × cell-type Spearman ρ; positive tracking with Myeloid (ρ = +0.39), Malignant (+0.30), B (+0.27), and Fibroblast (+0.23); negative tracking with T cell (ρ = −0.42) and Epithelial (−0.31) (Supplementary Figure SX panel I).

**DM1 sub-A vs sub-B teaser.** Sub-cluster labels from `d6p7_dm1_subcluster/dm1_subcluster_labels.tsv` (sub-A n = 84, sub-B n = 56) were intersected with cell-type fractions; Cohen's d and Mann-Whitney U two-sided p were reported per cell type. Sub-A is Malignant-cell rich (d = +1.22, p = 2.5 × 10⁻⁹) and Epithelial-poor (d = −1.26, p = 5.2 × 10⁻¹⁰); immune compartment differences are not significant — defining the sub-A/B split as a tumor-purity-vs-thyrocyte split rather than immune-hot vs immune-cold (Paper-2 boundary marker; Supplementary Figure SX panel J).

**Pseudotime trajectory.** Samples were ranked by canonical 8-gene score (TCGA `rai_score_recalc` n = 513; Lee/GSE213647 `panel_z` n = 630), binned into 10 deciles, and mean per-decile cell-type fraction was computed. Decile-level Spearman ρ between mean score and mean fraction quantifies monotonic trajectory; four compartments — Malignant (↓), Epithelial (↑), Myeloid (↓), Endothelial (↑) — show |ρ| ≥ 0.95 in both cohorts (Supplementary Figure SX panel K).

**Full-transcriptome reference robustness.** Pu 2021 raw counts (33,694 genes × 66,015 cells) were subsampled to 5,000 cells balanced across 7 patients (seed = 42); each cell was assigned a Lu 2023 `author_celltype` label by maximum cosine similarity over the Lu HVG ∩ Pu intersection (1,898 genes after Ensembl → symbol conversion). A full-transcriptome pseudobulk per cell type was constructed (per-cell-type log-normalized mean, library-size-corrected, n_genes = 33,694; `Pu_pseudobulk_full`). TCGA bulk was re-deconvolved by NNLS over the 21,369-gene Pu × TCGA intersection (`scipy.optimize.nnls`, sum-to-one normalization). Per-cell-type Cohen's d (DM1 − DM2) was compared to the v2 Lu HVG nu-SVR primary result; all four informative axes (Malignant, Epithelial, Myeloid, Endothelial) sign-match between the two references with Pu full NNLS magnitudes ≥ Lu HVG. NNLS sparse-collapse zeros the B / Fibroblast / NK / T compartments in the full-transcriptome regime; these are interrogated by the v2 nu-SVR Lu HVG primary instead.

### Random-effects meta-analysis

Cox-derived log-hazard-ratios and standard errors from TCGA-THCA and MSK-IMPACT were combined using random-effects DerSimonian-Laird estimation (`statsmodels.stats.meta_analysis`). Cochran I² was reported for between-cohort heterogeneity.

---

## Additional resources

- **GENCODE v44 reference annotation**: https://www.gencodegenes.org/human/release_44.html
- **TIERA67 gene definition**: Supplementary Table S1 (also available in `metadata/tierA67_genes.txt` in the source code repository)
- **Statistical analysis notebook**: All quantification code is reproducible from the source code repository, with per-figure script paths documented in the README.

---
