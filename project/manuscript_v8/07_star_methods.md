---
title: "Paper 1 manuscript v8 — STAR Methods (draft v1)"
date: 2026-05-01
author: Seungho Cook
target_format: Cell Press STAR Methods (unlimited length, reproducibility-focused)
status: v1 draft — 본인 voice 적용 영역 명시. Code release timeline + 정확 software versions verify.
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
| GSE286332 Korean PTC vs PTC+HT | Macrogen / Dongguk Univ | GEO GSE286332 |
| GSE184362 Pu 2021 single-cell PTC | Pu et al., 2021 | GEO GSE184362 |
| GSE193581 Lu 2023 single-cell PTC | Lu et al., 2023 | GEO GSE193581 |
| GSE241184 Phase 1 single-cell PTC | (Phase 1) | GEO GSE241184 |
| **Software** | | |
| pyDESeq2 (DEG analysis) | Snakemake/lab | https://github.com/owkin/PyDESeq2 |
| scikit-learn (KMeans, LogReg, AUC) | Pedregosa et al. | https://scikit-learn.org |
| lifelines (Cox PH, log-rank) | Davidson-Pilon | https://lifelines.readthedocs.io |
| arcasHLA (HLA imputation from RNA-seq) | Orenbuch et al. | https://github.com/RabadanLab/arcasHLA |
| STAR aligner | Dobin et al. | https://github.com/alexdobin/STAR |
| kallisto (pseudo-alignment) | Bray et al. | https://pachterlab.github.io/kallisto |
| GENCODE v44 reference | EBI | https://www.gencodegenes.org/human/release_44.html |
| IMGT/HLA v3.55 | EBI | https://www.ebi.ac.uk/ipd/imgt/hla |
| **Deposited data** | | |
| cBioPortal SV (RET/NTRK/ALK/BRAF fusions) | cBioPortal API | `thca_tcga_pub` study, SV endpoint |
| HM450 promoter methylation (Illumina HumanMethylation450) | cBioPortal API | `thca_tcga` legacy study |
| **Source code (this paper)** | | |
| 8-gene panel + DM cluster pipeline | Cook et al., this paper | `github.com/seungho-cook/v17_THCA_paper1` (release v1.0 + Zenodo DOI on acceptance) |
| arcasHLA Korean Pan-Asian forest | Cook et al., this paper | (same repo, branch `arcasHLA`) |

---

## Resource availability

**Lead contact.** Further information and requests for resources should be directed to the lead contact, Seungho Cook (kukshomr@gmail.com).

**Materials availability.** This study did not generate new unique reagents. All analyses were performed on publicly accessible datasets (TCGA, GEO, ENA, cBioPortal). Korean cohort access (K2 / PRJEB11591, GSE213647, GSE286332) is available via the listed repositories.

**Data and code availability.** All source code is available at https://github.com/seungho-cook/v17_THCA_paper1 (release v1.0 will be tagged on manuscript acceptance) and archived on Zenodo (DOI to be assigned). Intermediate data tables (TCGA per-sample DM scores, fusion partner annotations, methylation β-values per panel gene, Cox meta-analysis inputs) are deposited as Supplementary Tables S1–S10.

---

## Experimental model and study participant details

This study uses publicly available genomic and transcriptomic data from previously published cohorts:

- **TCGA-THCA** (n = 504 with overall survival annotation; 513 primary tumors total). Discovery cohort. Publicly available via The Cancer Genome Atlas (Cancer Genome Atlas Research Network, 2014).
- **MSK-IMPACT thyroid** (n = 117; advanced disease, mostly PDTC + ATC). Validation cohort. (Landa et al., 2016).
- **K2 / PRJEB11591** (n = 260; primary Korean PTC). Validation cohort. (Yoo et al., 2016).
- **Lee / GSE213647** (n = 632; Korean PTC). Validation cohort.
- **GSE286332** (n = 18; 9 Korean PTC + 9 Korean PTC+HT). Validation cohort.
- **GSE184362 Pu 2021** (n = 6 Fudan University PTC patients; single-cell). External validation.
- **GSE193581 Lu 2023** (n = 23 single-cell samples). External validation.
- **GSE241184** (n = 1; Phase 1 single-cell). Internal pilot.

All studies were originally approved by the respective institutional review boards. The present analysis used only de-identified public data and is exempt from additional IRB review.

---

## Method details

### Cohort assembly and clinical metadata harmonization

Bulk RNA-seq quantifications were obtained as log2(TPM + 1) (TCGA) and log2(FPKM + 1) (Korean cohorts and GSE286332) and processed through unified gene-level filtering (≥10 reads in ≥30% of samples; GENCODE v44 protein-coding annotation). Clinical metadata (age, sex, stage, vital status, time-to-event) were harmonized from cBioPortal and source publications. Per-cohort missingness was tabulated (Supplementary Table S2) and addressed in sensitivity analyses.

### 8-gene panel selection

The 67-gene candidate pool (TIERA67) was assembled from seven thyroid-relevant biological categories (Supplementary Table S1): TDS-core differentiation markers (16 genes), MAPK-output transcripts (10), thyroid driver genes (12, including BRAF, NRAS, HRAS, KRAS, RET, NTRK1, NTRK3, ALK, PAX8, PPARG, TERT, EIF1AX), aggressive-disease markers (10), dedifferentiation/EMT markers (10), light immune-stromal markers (5), and thyroid-lineage extras (4). The 8-gene panel (SLC5A5/NIS, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, DIO1) was selected from the TDS-core sub-category based on canonical RAI-uptake biology (Yoo et al., 2016 *PLOS Genet*) — independently of and prior to access of the Landa et al. (2016) ATC-silenced gene list. This panel design therefore predates exposure to advanced-disease ATC transcriptomic outputs and represents an independent path to the same differentiation axis.

Panel size sensitivity was tested across n = 8, 10, 12, and 16 gene variants (Supplementary Figure S2). The 16-gene full TDS-core yielded TCGA 5-fold AUC of 0.975 versus the 8-gene panel's 0.962 (ΔAUC = 0.013, NS); the 10-gene and 12-gene intermediate panels yielded 0.972 and 0.969, respectively. The 8-gene panel was retained for clinical interpretability.

### DM1/DM2 cluster definition

Bulk RNA expression of the 8-gene panel was z-standardized within cohort, and KMeans (k = 2, scikit-learn default initialization) was applied to assign DM1 versus DM2 cluster identity. The cluster boundary was confirmed by silhouette analysis (mean silhouette score > 0.4 within both clusters). For per-sample classification probability (P(DM1)), a logistic regression classifier was trained on the TCGA DM1/DM2 partition using the within-sample-centered profile of the 8-gene panel (mean panel score subtracted), yielding cohort-portable class probabilities (cross-validated AUC = 0.968 in TCGA; Methods, Supplementary Figure S3).

### Single-cell external validation

For GSE184362 (Pu et al., 2021), per-patient pseudo-bulk DM1 scores were computed by averaging thyrocyte-marker-positive (KRT8, KRT19, EPCAM) cell profiles per patient and per condition (tumor vs adjacent normal). Per-patient Spearman correlation between tumor and normal scores was computed. Lu 2023 (GSE193581) was similarly processed with stromal/immune contamination control.

### Survival analysis and meta-analysis

Cox proportional hazards regression was performed using `lifelines` (Davidson-Pilon, Python). Pooled meta-analysis was performed using random-effects DerSimonian-Laird estimation; between-cohort heterogeneity was assessed via Cochran I². Sub-cluster Cox analyses (DM1 sub-A vs sub-B) used patient-level event data within TCGA-THCA primary tumors (Supplementary Table S8).

### cBioPortal API access (SV + methylation)

Structural variants (RET, NTRK1/3, ALK, BRAF, PAX8-PPARG, others) were retrieved from cBioPortal `thca_tcga_pub` study via REST API (endpoint `/structural-variant/fetch`); SV-tested status was captured per tumor (n = 542 / 557 = 97.3%). HM450 promoter β-values for the 8-gene panel + DIO2 + SLC26A4 were retrieved from cBioPortal `thca_tcga` legacy study (n = 503 with HM450 + DM call); per-gene aggregate β values were computed from cBioPortal merged probe-per-gene format.

### Korean cohort processing

Yoo 2016 K2 (PRJEB11591): kallisto pseudo-alignment to GENCODE v44 transcriptome, transcript-to-gene aggregation, log2(TPM + 1) normalization. The 8-gene mini-index calibration mismatch (R4-4) was identified and addressed via within-sample-centered profile classification (Memory cross-ref: `v17_korean_k2_calibration.md`). Lee (GSE213647) and GSE286332 used pre-computed FPKM tables. arcasHLA HLA-allele imputation was performed on PRJEB11591 (n = 260) and GSE286332 (n = 18); 4-digit allele frequencies were aggregated for Korean PTC pool (n = 874).

### Software and statistical environment

All analyses were performed in Python 3.11 with scikit-learn 1.5, pandas 2.2, numpy 1.26, scipy 1.13, lifelines 0.27, and pyDESeq2 0.4. Code is available at https://github.com/seungho-cook/v17_THCA_paper1.

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

### Mediation analysis (Baron-Kenny)

For autoimmune-PTC mediation (P3 / D3-P5 framework), Baron-Kenny mediation was performed with HLA-II module as candidate mediator. Indirect effect (a × b) and percent-mediated were computed; percentile bootstrap 95% CI was used. (Pillar 5, Paper 2 reserve.)

### Random-effects meta-analysis

Cox-derived log-hazard-ratios and standard errors from TCGA-THCA and MSK-IMPACT were combined using random-effects DerSimonian-Laird estimation (`statsmodels.stats.meta_analysis`). Cochran I² was reported for between-cohort heterogeneity.

---

## Additional resources

- **GENCODE v44 reference annotation**: https://www.gencodegenes.org/human/release_44.html
- **IMGT/HLA v3.55 reference alleles**: https://www.ebi.ac.uk/ipd/imgt/hla
- **TIERA67 gene definition**: Supplementary Table S1 (also available in `metadata/tierA67_genes.txt` in the source code repository)
- **Statistical analysis notebook**: All quantification code is reproducible from the source code repository, with per-figure script paths documented in the README.

---

# 본인 voice 적용 영역

- [ ] **8-gene panel selection** paragraph — independent design + Yoo 2016 first cite (reverse-causality 차단 Layer 1)
- [ ] **Code release timeline** — paper publish 시 release vs preprint 시 release 결정
- [ ] **Korean cohort processing** — K2 mini-index calibration mismatch acknowledgment (R4-4)
- [ ] **GitHub repo URL** — `seungho-cook/v17_THCA_paper1` 본인 username 검증 + 실제 repo 생성 timeline

# 다음 step

1. Software version 정확 verify (5/4)
2. Code release timeline 본인 결정
3. Korean cohort processing methodology 본인 PDF read (Yoo 2016) 후 정확 표현 검증
4. Prompt 7 (Cover letter + Reviewer Q&A) 진행
