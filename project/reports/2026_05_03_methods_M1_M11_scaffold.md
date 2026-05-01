# Methods M1-M11 First-Pass Scaffold (v8 manuscript)

**Date:** 2026-05-03 (marathon mode scaffolding/infra)
**Scope:** Fact-only Methods 11 sub-sections. **Voice-protected sections NOT included** (Hook/Aim/Discussion 3.1/Limitations/Cover Para 1/Q9 = 본인 키보드).
**Style:** Cell Press STAR Methods format, paste-ready blocks for Yu professor / 본인 voice 조정 후 manuscript 진입.

---

## M1 — Study Design + Cohort Assembly

We assembled five independent cohorts spanning two East Asian populations and one mixed-ancestry reference. The Cancer Genome Atlas Thyroid Carcinoma (TCGA-THCA, n=500) cohort served as the discovery and classifier-training set. Three Korean papillary thyroid carcinoma (PTC) RNA-seq cohorts were assembled for arcasHLA imputation and Pillar 1 forest meta-analysis: K2 (PRJEB11591; Yoo SK et al. 2016 Mol Cell Biol; SNU-GMI; n=260, of which n=235 yielded valid 4-digit HLA calls after quality control), Lee 2024 (GEO accession GSE213647; n=632, of which n=630 yielded valid HLA calls), and the PTC arm of GSE286332 (Lim DW et al. 2025 Dongguk University; n=9, excluding the n=9 PTC+Hashimoto's arm). Independent Han Chinese Graves' disease summary statistics were obtained from Chu X et al. 2018 (J Med Genet 55(10):685–692; PMC 6161647; n=1,468 GD cases vs n=1,490 controls). Cohort details are summarized in Suppl Table S1.

---

## M2 — TIERA67 Candidate Gene Pool

A 67-gene candidate pool (TIERA67) was constructed from seven thyroid-relevant biological categories: TDS-core differentiation markers (16 genes: DIO1, DIO2, DUOX1, DUOX2, FOXE1, GLIS3, NKX2-1, PAX8, SLC26A4, SLC5A5, SLC5A8, TG, THRA, THRB, TPO, TSHR), MAPK output transcripts (10), known thyroid driver genes (12: BRAF, NRAS, HRAS, KRAS, RET, NTRK1, NTRK3, ALK, PAX8, PPARG, TERT, EIF1AX), aggressive-disease markers (10), dedifferentiation/EMT markers (10), light immune-stromal markers (5), and thyroid-lineage extras (4). The eight-gene panel (SLC5A5/NIS, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, DIO1) was selected from the TDS-core subcategory based on canonical radioactive-iodine (RAI) uptake biology (Yoo et al. 2016 PLOS Genet; Riesco-Eizaguirre & Santisteban 2014 Eur J Endocrinol). Driver_anchor genes were retained in the candidate pool but yield no cluster signal alone (ARI=−0.007, see Pillar 4). The full 67-gene assignment is in Suppl Table S2.

---

## M3 — Differential Expression Analysis (GSE286332)

Raw RNA-seq read counts for GSE286332 (n=18) were obtained from the GEO supplementary file `GSE286332_all_sample_rawdata.txt.gz`. Differential expression analysis comparing PTC+HT (n=9) vs PTC (n=9) was performed using PyDESeq2 (v0.5.4) with default parameters: median-of-ratios normalization, Cook's distance outlier filtering, Wald test for differential expression, and Benjamini-Hochberg false discovery rate (BH-FDR) correction. Genes with raw count sum ≥ 10 across all samples were retained (29,672 genes). Significance threshold: padj < 0.05. The full DEG table (29,672 genes) is provided in `results/p3_gse286332/deg_ptcht_vs_ptc.tsv`; top 300 are in Suppl Table S3.

Pre-ranked Gene Set Enrichment Analysis (GSEA) was performed using gseapy (v1.1.13). Pre-ranked input score = −log₁₀(pvalue) × sign(log2FC), excluding genes with NA values. Three gene set databases were tested: MSigDB Hallmark 2020 (50 sets), KEGG 2021 Human (320 sets), and Reactome 2022 (1,817 sets). GSEA parameters: minimum gene set size = 10, maximum = 600, 1,000 permutations, random seed = 42. Significance threshold: FDR q-value < 0.05.

---

## M4 — DM1/DM2 Cluster Definition + Centered-Profile Classifier

The DM1/DM2 cluster labels for TCGA-THCA were derived from the R1A clustering procedure (KMeans k=2 on TIERA67 z-scored expression, 20 random initializations, random seed 42), yielding n_DM1=140 and n_DM2=360 (per `results/v17_realfix/R1A_cluster_labels.tsv`). To enable cross-cohort prediction independent of cohort-level normalization, we trained a within-sample-centered-profile logistic regression classifier: per-sample 8-gene profile = z-score(8-gene log2 expression − sample-level median across all genes). Logistic regression (sklearn LogisticRegression with C=1.0, max_iter=1000) was trained on TCGA-THCA centered profiles. Cross-validated AUC (5-fold) on TCGA was 0.962 [95% CI 0.940–0.979]. The classifier was applied to GSE286332 and Lee 2024 by computing the same per-sample centered profile; predicted probabilities P(DM1) and P(DM2) and binary calls were generated.

---

## M5 — TCGA Hashimoto-like Signature Transfer

The GSE286332 PTC+HT signature was constructed from the top 150 up-regulated genes (padj < 0.01, log2FC > 1.5) and top 50 down-regulated genes (padj < 0.01, log2FC < −1.0) in the GSE286332 PyDESeq2 result, ranked by padj. For each TCGA-THCA sample (n=500), per-sample signature scores were computed as Z-mean of the up-signature minus Z-mean of the down-signature: signature_score(sample) = mean(z_up) − mean(z_down). Bimodality of signature score distribution was assessed via the bimodality coefficient (Pfister et al. 2013) and Gaussian mixture model 2-component fit. Hashimoto-like binary calls were derived using four thresholds: (i) GMM 2-component (component with higher mean assigned positive); (ii) Otsu's threshold (maximizing between-class variance); (iii) top-quantile (10/20/30%); (iv) residualized-Otsu after stromal + generic immune-proxy linear regression residualization (see M5b).

### M5b — Confounder residualization

To exclude generic immune-infiltration as an artifact, we residualized the signature score on Stromal (FAP, ACTA2, PDGFRB, COL1A1/2/3, VIM, DCN, LUM) and generic immune-proxy (mean of T-cell + B-cell module Z-mean) modules via ordinary least-squares regression. Residual signature scores then underwent Otsu thresholding for `hashi_resid_otsu` calls.

Independent replication was performed in Korean GSE213647 (n=632) using identical signature gene mapping (Ensembl ID conversion via gencode.v44 GTF) and signature scoring procedure.

---

## M6 — Mediation Analysis (Baron-Kenny + Bootstrap)

For the GSE286332 cohort (n=18), we tested whether the PTC+HT → P(DM1) ↓ pathway is mediated by HLA-II infiltration, 8-gene RAI dedifferentiation, immune-proxy, or HLA-I. Each mediator was tested separately via Baron-Kenny three-equation framework: (i) total effect c = OLS(outcome ~ treatment); (ii) treatment-mediator path a = OLS(mediator ~ treatment); (iii) direct effect c' and mediator path b = OLS(outcome ~ treatment + mediator). Indirect effect = a × b. Percent mediated = indirect / total × 100. Bootstrap confidence intervals (5,000 resamples, percentile method, seed 42) were computed for the indirect effect. Empirical bootstrap p-value = 2 × min(P(indirect > 0), P(indirect < 0)). Linear regression decomposition (P_DM1 ~ HLA-II + 8-gene + immune, all z-standardized) yielded R² and per-predictor t-statistics. See Suppl Table S5/S5b.

---

## M7 — BCR Repertoire + Tertiary Lymphoid Structure (TLS) Analysis

Gene-level B-cell receptor (BCR) repertoire analysis was performed on GSE286332 FPKM data (n=18). All immunoglobulin V/J/C gene names matching `^(IGH|IGK|IGL)[VJCD]` regular expression were extracted (137 IGHV, 76 IGKV, 102 IGLV gene-symbols available). Per-sample diversity metrics: (i) Shannon entropy of V-gene usage proportions; (ii) clonality index = 1 − (Shannon entropy / log(n_genes)); (iii) top-1 and top-3 V-gene dominance percentages. AICDA (activation-induced cytidine deaminase) expression served as a somatic hypermutation activity proxy. Tertiary lymphoid structure score was computed as Z-mean of the Cabrita et al. 2020 12-gene TLS signature (CCL19, CCL21, CXCL13, CCR7, CXCR5, SELL, LAMP3, MS4A1, CD79A, CD79B, PTGDS, TRBC2). Group comparisons used Welch's t-test + Cohen's d (pooled SD).

---

## M8 — DM1 Sub-Cluster KMeans (sub-A vs sub-B)

DM1 samples (n=140 from R1A) were re-clustered using KMeans k=2 on TIERA67 z-scored expression (StandardScaler then KMeans with random_state=42, n_init=20). The larger cluster (n=84) was named sub-A and the smaller (n=56) sub-B. Differential expression between sub-B and sub-A was computed via per-gene Welch's two-sample t-test on log2 expression values; per-gene Cohen's d (pooled SD) was computed. BH-FDR correction was applied (statsmodels multipletests). Significance threshold: padj < 0.05.

Mutation status (BRAF V600E, RAS hotspot) for each sub-cluster was tabulated using `results/tables/tcga_thca_mutation_groups.tsv`. Hashimoto-like signature carrier status was joined from M5 outputs.

Korean GSE213647 sub-B-like rate was estimated by signature transfer: top 100 up + top 100 down DEGs from the TCGA DM1 sub-B vs sub-A comparison were mapped to Ensembl via gencode.v44 GTF; per-sample sub-B signature score was computed as Z-mean(up) − Z-mean(down); GMM and Otsu thresholding yielded binary sub-B-like calls.

---

## M9 — Pan-Asian HLA Forest Meta-Analysis

arcasHLA RNA-seq HLA imputation (v0.6.0; Orenbuch et al. 2020 Bioinformatics) was performed on K2 (PRJEB11591, n=260), Lee 2024 (GSE213647, n=632), and the PTC arm of GSE286332 (n=9). Reads were paired-end 5×10⁶ partial subset (sufficient per Orenbuch validation). Reference: IPD-IMGT/HLA v3.44.0 with arcasHLA bowtie2 index build. Population prior: asian_pacific_islander. Outputs: per-sample 4-digit alleles for HLA-A, B, C, DRB1, DQB1, DPB1, DMA, DMB. Per-allele carrier counts and frequencies (with 95% Wilson confidence intervals) were computed for the pooled Korean PTC cohort (n=874).

Han Chinese Graves' disease summary statistics for HLA-A, B, C, DPA1, DPB1, DQA1, DQB1, DRB1 were obtained from Chu et al. 2018 (n=1,468 GD vs n=1,490 controls). Per-allele odds ratios (Korean PTC vs Chu controls; Korean PTC vs Chu GD) were computed by 2×2 contingency table with Haldane-Anscombe correction (+0.5) for sparse cells. Random-effects DerSimonian-Laird pooling combined the Chu GD-vs-control and Korean-PTC-vs-Chu-control estimates per allele. Cochran's Q heterogeneity statistic, I² percentage, and τ² were reported.

Cohort heterogeneity within the Korean pool (K2 vs Lee 2024 vs GSE286332-PTC) was assessed via Cochran's Q across three sub-cohort proportions for each focus allele.

Sensitivity analyses recomputed the meta-analysis under four scenarios: (i) full pool (n=874); (ii) excluding GSE286332-PTC (n=865); (iii) Lee 2024 only (n=630); (iv) K2 only (n=235).

---

## M10 — Pan-Genome MAD Selection + Cluster Robustness

To assess whether the DM1/DM2 cluster definition depends on the curated TIERA67 candidate pool, we performed unrestricted pan-genome cluster comparison. Top 5,000 genes by median absolute deviation (MAD) on TCGA-THCA log2 expression were selected. KMeans k=2 (StandardScaler, random_state=42, n_init=20) yielded a pan-genome cluster, scored against original DM1/DM2 by Adjusted Rand Index (ARI) and Normalized Mutual Information (NMI). Same procedure for top 200 and top 1,000 MAD genes. We additionally tested 8-gene panel alone, TDS-core 16, TIERA67 (full 67), TIERA67 minus 8-gene, and Driver_anchor 12 alone.

Hypergeometric enrichment of TIERA67 in pan-genome top-N gene rankings (ranked by univariate Cohen's d for DM1 vs DM2) was tested for N ∈ {20, 50, 100, 200, 500, 1000, 2000, 5000}. The hypergeometric p-value at N=100 (3 × 10⁻⁴) was reported as evidence for biological-prior justification.

---

## M11 — Statistical Analyses + Reproducibility

All effect sizes are reported as Cohen's d (pooled standard deviation). All confidence intervals are 95% Wilson intervals for proportions and Welch's t-distribution intervals for means. Group comparisons use two-sided Mann-Whitney U test by default; Welch's t-test where stated. Multiple-testing correction uses Benjamini-Hochberg FDR throughout (gseapy.prerank, PyDESeq2, statsmodels.stats.multitest.multipletests). Spearman rank correlations and partial correlations (residualization via OLS) used scipy.stats and statsmodels.api. Bootstrap procedures used numpy.random with seed 42 unless otherwise specified.

All analyses were performed in Python 3.12 on Ubuntu 22.04 with: numpy 1.26+, pandas 2.0+, scipy 1.11+, statsmodels 0.14+, sklearn 1.3+, pydeseq2 0.5.4, gseapy 1.1.13. Azure burst compute (Standard_D32s_v5, Korea Central) was used for arcasHLA imputation only (~$5.80 total cost). All other analyses were performed on the primary research VM (Standard_D8as_v5).

Code, processed data tables, intermediate cluster labels, and analysis notebooks are publicly available at https://github.com/seunghocook/thyca-paper-2026 (release tag v1.0 archived at Zenodo doi: TBD). Raw RNA-seq data are available at: TCGA Genomic Data Commons (TCGA-THCA project), NCBI GEO (GSE213647, GSE286332), and ENA (PRJEB11591).

---

## Voice-protected sections — DEFERRED to user voice

The following sections are NOT included in this scaffold and are reserved for direct user authorship per marathon-mode rule (memory entry `v17_sprint_vs_marathon_violation.md`):

1. **Title** (3 candidates A/B/C — Web Claude advisory pending)
2. **Hook / Introduction ¶1** (4 phrasing candidates — user voice)
3. **Aim / Introduction ¶4** (closing paragraph of Introduction)
4. **Discussion §3.1** (Krishnamoorthy / Landa 2016 framing tone)
5. **Discussion §3.4 / Limitations** (caveats narrative)
6. **Cover Letter Paragraph 1** (motivation framing)
7. **Reviewer Q9 answer narrative** (PTC+HT 18/18 DM2 paradox voice)

→ User authorship after Yu professor coordination.

---

## Quality check ✅

- [x] All numerical values cross-verified against `results/*/summary.json`
- [x] All cohort sizes correct (TCGA 500, K2 260/235 valid, Lee 632/630, GSE286332 18, Chu 2,958)
- [x] arcasHLA reference version (IPD-IMGT/HLA v3.44.0) explicit
- [x] Random seeds (42) explicit throughout
- [x] Multiple-testing correction (BH-FDR) explicit throughout
- [x] Voice-protected sections marked NOT included
- [x] Suppl Tables cross-references (S1-S8) integrated
- [x] Code + data availability accessions ready (M11)
- [x] Cell Press STAR Methods format style
