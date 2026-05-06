# External expression validation — CANONICAL SINGLE-MD REPORT

**Date:** 2026-05-04
**Paper scope:** Paper 1 (DM1 molecular dark matter)
**Question:** Does the driver-orthogonal RAI-lineage / DM1-differentiation axis reproduce in independent thyroid expression cohorts?

**Canonical decision:** use this file as the single report for downstream review/commit. Treat the other `2026_05_04_external_expression_*.md` files as intermediate scratch outputs unless a specific audit trail is needed.

This file consolidates everything from the sweep:
- Access check (Step 1)
- Sweep configuration & forbidden-action audit
- Used-cohort sample summary
- Gene-panel coverage
- Per-dataset score tests (Mann–Whitney + Cohen's d) — full numbers
- Direction-consistency forest table — full numbers
- Spearman correlations — full numbers
- Figures
- Commit proposal (NOT executed)

Intermediate source reports (not needed if committing only one md):
- `project/reports/2026_05_04_external_expression_access_check.md`
- `project/reports/2026_05_04_external_expression_validation_report.md`
- `project/reports/2026_05_04_external_expression_commit_proposal.md`
- `project/results/p_external_expression_validation/` (data + figures + scripts)

---

## 0. TL;DR

**REPLICATED — direction-consistent lineage silencing in 4 independent GPL570 cohorts (n=205 used).**

| Cohort | n | Headline (advanced vs PTC/normal) | DM1_like ↔ NONOVERLAP Spearman ρ |
|--------|---|----------------------------------|-----------------------------------|
| GSE33630 | 105 (11 ATC / 49 PTC / 45 normal) | RAI_8 ATC vs normal **d=−5.48** (p=3.5e-7) | **−0.93** (p=9e-46) |
| GSE65144 | 25 (12 ATC / 13 normal) | RAI_8 ATC vs normal **d=−2.79** (p=4.0e-5) | **−0.94** (p=7e-12) |
| GSE29265 | 49 (9 ATC / 20 PTC / 20 paired-normal) | RAI_8 ATC vs normal d=−2.29 (p=7.5e-4) | −0.85 (p=8e-15) |
| GSE53157 | 26 after pool drop (5 PDTC / 7 PTC / 8 FVPTC / 4 FTC / 2 normal) | PDTC vs DTC d=−1.05 (p=0.09) — **direction-consistent, underpowered** | −0.84 (p=9e-8) |

**Dropped:** GSE126698 (Series Matrix metadata-only; raw=SRA-only; alignment forbidden).
**Held (per user spec, not downloaded):** GSE53072, GSE60542, GSE120177.

Forbidden actions audit (all ✅): no extra GEO search, no FASTQ alignment, no CEL bulk, no GPU/RunPod, no TCGA WSI/methylation, no DepMap/CCLE/PRISM, no Paper 3/4 touch, no manuscript prose, no voice-protected section drafted, no commit performed.

---

## 1. Access check (Step 1)

GEO landing pages were fetched for the user-supplied 8-dataset list only. **No additional GEO search/download performed.**

| # | Tier | GEO | Title (short) | Platform | Modality | n | Histology groups (n) | Series Matrix | Per-sample expr | Raw | Verdict |
|---|------|-----|---------------|----------|----------|---|----------------------|---------------|-----------------|-----|---------|
| 1 | T1 | GSE126698 | IGF2BP1 ATC marker | GPL15456 (Illumina HiScanSQ) | RNA-seq | 28 | ATC 10 / PTC 6 / FTC 6 / Normal 6 | Yes (87 lines) | **No — DE summary only** | SRA (PRJNA523137 / SRP186236) | **DROP at access** — alignment forbidden by spec |
| 2 | T2 | GSE33630 | Normal vs PTC vs ATC | GPL570 | Affy U133+2 | 105 | ATC 11 / PTC 49 / Normal 45 | Yes | Yes | CEL TAR (849 MB) | **GO** |
| 3 | T2 | GSE29265 | Sporadic vs Chernobyl PTC + ATC | GPL570 | Affy U133+2 | 49 | ATC 9 / PTC 20 / Paired-normal 20 | Yes | Yes | CEL TAR (213 MB) | **GO** |
| 4 | T2 | GSE65144 | ATC vs matched/unmatched normal | GPL570 | Affy U133+2 | 25 | ATC 12 / Normal 13 | Yes | Yes | CEL TAR (107 MB) | **GO** |
| 5 | T2 | GSE53157 | PDTC progression series | GPL570 | Affy U133+2 | 27 | PDTC 5 / cPTC 7 / fvPTC 8 / FTC 4 / Normal 2 / pool 1 | Yes | Yes | CEL TAR (126 MB) | **GO** (drop pool) |
| 6 | T3 hold | GSE53072 | Tiny ATC vs normal | GPL6244 | Array | 9 | — | — | — | — | **HOLD** (tiny + non-GPL570) |
| 7 | T3 hold | GSE60542 | PTC primary vs nodal mets | GPL570 | Affy U133+2 | 92 | — | — | — | — | **HOLD** (lymphoid confounding) |
| 8 | T3 hold | GSE120177 | CDK7/THZ1 in ATC cell lines | GPL23227 (BGISEQ-500) | RNA-seq + ChIP | 12 | — | — | — | — | **HOLD** (Paper 9 scope; perturbation not histology) |

GSE126698 verification: Series Matrix is 87 lines, no `series_matrix_table_begin/end`-flanked expression block; supplementary files are `GSE126698_DE_Thyroid_totalRNA_*.csv.gz` (DE summary tables, not per-sample counts). Raw data = SRA only. Per spec, raw FASTQ alignment is forbidden → marked FAIL/DROP.

---

## 2. Sweep configuration

### 2.1 Gene panels

| Panel | Genes |
|-------|-------|
| RAI_8 | SLC5A5, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, DIO1 |
| THYROID_NONOVERLAP | SLC26A4, IYD, DUOX1, DUOX2, TFF3, HHEX, GLIS3, DIO2 |
| TDS_TF | FOXE1, NKX2-1, PAX8, HHEX |
| MECHANISM (STAT3_AP1_DNMT) | STAT3, FOSL1, JUNB, DNMT1, DNMT3B (TACSTD2 reported separately) |
| OPTIONAL_TARGETS | NAMPT, KCNN4, LYN, OSMR, IL6R, SRC, FYN, ATR, CHEK2, MYC, GLS, LDHA |

Aliases applied: `NKX2-1` ⇄ `TITF1`, `FOXE1` ⇄ `TTF2/FKHL15`, `SLC5A5` ⇄ `NIS` (none required for these 4 series — canonical symbols all present).

### 2.2 Pipeline (per-dataset, no pooling)

1. Parse Series Matrix → metadata + probe × sample matrix.
2. GPL570 SOFT annot.gz (`Aug 09 2016`, 45,118 probes mapped) → probe → symbol(s) (split on `///`).
3. Probe → gene collapse: **mean** across probes mapping to a symbol → 22,836 genes per dataset.
4. Auto-scale check (max < 50) → no re-log2 applied (GPL570 series matrices are RMA log2).
5. Within-dataset z-score per gene.
6. Panel scores = mean of available z's (RAI_8, NONOVERLAP, TDS_TF, STAT3_AP1_DNMT). DM1_like = −RAI_8.
7. Mann–Whitney U + Cohen's d per contrast; Spearman correlations across DM1_like ↔ {NONOVERLAP, TF_collapse, STAT3_AP1_DNMT, TACSTD2}.
8. Direction-consistency forest: advanced (ATC ∪ PDTC) vs DTC (PTC ∪ FVPTC ∪ FTC) and vs normal.

### 2.3 Forbidden-action audit

| Constraint | Status |
|------------|--------|
| No GEO accessions outside the 8 supplied | ✅ |
| No raw FASTQ alignment | ✅ |
| No raw CEL bulk processing | ✅ |
| No RunPod / GPU | ✅ |
| No TCGA WSI touch | ✅ |
| No H&E-DM1 retry | ✅ |
| No TCGA methylation touch | ✅ |
| No DepMap / CCLE / PRISM touch | ✅ |
| No Paper 3 / Paper 4 work | ✅ |
| No manuscript prose written | ✅ |
| No voice-protected section drafted (Hook / Aim / Discussion / Limitations / Cover / Q9) | ✅ |
| No commit performed | ✅ |

---

## 3. Sample summary (used cohorts, n=205)

| Dataset | n_total | ATC | PDTC | PTC | FVPTC | FTC | normal | dropped | Histology source |
|---------|---------|-----|------|-----|-------|-----|--------|---------|------------------|
| GSE33630 | 105 | 11 | — | 49 | — | — | 45 | 0 | `Sample_source_name_ch1` ("non-tumor") + `Sample_characteristics_ch1` ("anaplastic"/"papillary") |
| GSE29265 | 49 | 9 | — | 20 | — | — | 20 | 0 | `Sample_title` ("Anaplastic …" / "Papillary …" / "Patient-matched non-tumor …") |
| GSE65144 | 25 | 12 | — | — | — | — | 13 | 0 | `Sample_source_name_ch1` ("Anaplastic Thyroid Carcinoma" / "Normal …") |
| GSE53157 | 27 | — | 5 | 7 | 8 | 4 | 2 | 1 commercial pool | `Sample_source_name_ch1` (`fresh-frozen <hist>`) |
| **Total used** | **205** | 32 | 5 | 76 | 8 | 4 | 80 | 1 pool | — |

Counts match GEO landing-page declarations exactly. Full per-sample table: `project/results/p_external_expression_validation/sample_metadata.tsv` (206 rows + header).

---

## 4. Gene-panel coverage

Every panel gene was found in every used dataset after probe → gene collapse.

| Panel | GSE33630 | GSE29265 | GSE65144 | GSE53157 |
|-------|----------|----------|----------|----------|
| RAI_8 (8) | 8/8 | 8/8 | 8/8 | 8/8 |
| THYROID_NONOVERLAP (8) | 8/8 | 8/8 | 8/8 | 8/8 |
| TDS_TF (4) | 4/4 | 4/4 | 4/4 | 4/4 |
| MECHANISM (6) | 6/6 | 6/6 | 6/6 | 6/6 |
| OPTIONAL_TARGETS (12) | 12/12 | 12/12 | 12/12 | 12/12 |

Per-gene matched-symbol audit → `external_gene_coverage.tsv` (152 rows). QC heatmap → `external_dataset_qc_heatmap.png`.

---

## 5. Per-dataset score tests — FULL numbers

`external_score_tests.tsv` (Mann–Whitney p, Cohen's d). Cohen's d sign convention: positive = higher in group1.

| dataset | contrast | n1 | n2 | RAI_8 d | RAI_8 p | DM1_like d | DM1_like p | NONOVERLAP d | NONOVERLAP p | TDS_like d | TDS_like p | TF_collapse d | TF_collapse p | STAT3_AP1_DNMT d | STAT3_AP1_DNMT p | TACSTD2 d | TACSTD2 p |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GSE33630 | ATC_vs_PTC | 11 | 49 | **−3.906** | 6.80e-07 | +3.906 | 6.80e-07 | **−4.216** | 4.57e-07 | −4.202 | 3.38e-07 | −4.202 | 3.38e-07 | +1.187 | 5.60e-03 | −1.584 | 1.02e-03 |
| GSE33630 | ATC_vs_normal | 11 | 45 | **−5.481** | 3.51e-07 | +5.481 | 3.51e-07 | **−6.723** | 3.51e-07 | −4.808 | 3.91e-07 | −4.808 | 3.91e-07 | +1.655 | 2.42e-04 | +0.044 | 7.57e-01 |
| GSE33630 | PTC_vs_normal | 49 | 45 | −2.253 | 3.50e-13 | +2.253 | 3.50e-13 | −1.978 | 8.07e-13 | −1.400 | 7.36e-09 | −1.400 | 7.36e-09 | +0.563 | 1.08e-03 | +2.139 | 8.22e-12 |
| GSE29265 | ATC_vs_PTC | 9 | 20 | **−1.565** | 1.96e-02 | +1.565 | 1.96e-02 | −1.515 | 1.17e-02 | −2.209 | 1.46e-03 | −2.209 | 1.46e-03 | +1.195 | 1.33e-02 | −2.062 | 3.11e-04 |
| GSE29265 | ATC_vs_normal | 9 | 20 | −2.294 | 7.50e-04 | +2.294 | 7.50e-04 | −2.663 | 4.55e-05 | −2.630 | 1.48e-04 | −2.630 | 1.48e-04 | +1.141 | 1.52e-02 | −0.861 | 8.89e-03 |
| GSE29265 | PTC_vs_normal | 20 | 20 | −1.161 | 7.58e-04 | +1.161 | 7.58e-04 | −1.271 | 4.16e-04 | −0.664 | 2.56e-02 | −0.664 | 2.56e-02 | −0.031 | 6.75e-01 | +1.230 | 8.36e-04 |
| GSE65144 | ATC_vs_normal | 12 | 13 | **−2.792** | 4.01e-05 | +2.792 | 4.01e-05 | **−3.531** | 3.17e-05 | −2.934 | 3.17e-05 | −2.934 | 3.17e-05 | +1.200 | 7.09e-03 | −0.044 | 9.78e-01 |
| GSE53157 | PDTC_vs_PTC | 5 | 7 | −0.279 | 1.00 | +0.279 | 1.00 | +0.155 | 0.876 | −0.199 | 1.00 | −0.199 | 1.00 | −1.080 | 0.149 | −2.693 | 5.05e-03 |
| GSE53157 | PTC_vs_normal | 7 | 2 | −1.912 | 0.222 | +1.912 | 0.222 | −3.287 | 5.56e-02 | −0.672 | 1.00 | −0.672 | 1.00 | +4.204 | 5.56e-02 | +3.334 | 5.56e-02 |
| GSE53157 | PDTC_vs_normal | 5 | 2 | −1.425 | 0.190 | +1.425 | 0.190 | −1.696 | 9.52e-02 | −0.519 | 0.857 | −0.519 | 0.857 | +3.317 | 9.52e-02 | −0.421 | 0.381 |
| GSE53157 | FVPTC_vs_PTC | 8 | 7 | +1.504 | 2.89e-02 | −1.504 | 2.89e-02 | +2.533 | 1.24e-03 | +1.594 | 9.32e-03 | +1.594 | 9.32e-03 | −1.496 | 4.01e-02 | −2.363 | 2.18e-03 |
| GSE53157 | FTC_vs_normal | 4 | 2 | +0.743 | 0.533 | −0.743 | 0.533 | −0.415 | 1.00 | +1.162 | 0.267 | +1.162 | 0.267 | +1.995 | 0.133 | −2.234 | 0.133 |

Bold = headline replication contrast.

---

## 6. Direction-consistency forest — FULL table

`external_direction_consistency.tsv`. Sign convention: positive d ⇒ higher in advanced group (ATC ∪ PDTC).

| dataset | contrast | score | cohen_d | p | n1 | n2 |
|---|---|---|---|---|---|---|
| GSE33630 | advanced_vs_DTC | RAI_8_score | −3.906 | 6.80e-07 | 11 | 49 |
| GSE33630 | advanced_vs_DTC | DM1_like_score | +3.906 | 6.80e-07 | 11 | 49 |
| GSE33630 | advanced_vs_DTC | THYROID_NONOVERLAP_score | −4.216 | 4.57e-07 | 11 | 49 |
| GSE33630 | advanced_vs_DTC | TDS_like_score | −4.202 | 3.38e-07 | 11 | 49 |
| GSE33630 | advanced_vs_DTC | TF_collapse_score | −4.202 | 3.38e-07 | 11 | 49 |
| GSE33630 | advanced_vs_DTC | STAT3_AP1_DNMT_score | +1.187 | 5.60e-03 | 11 | 49 |
| GSE33630 | advanced_vs_DTC | TACSTD2_z | −1.584 | 1.02e-03 | 11 | 49 |
| GSE33630 | advanced_vs_normal | RAI_8_score | −5.481 | 3.51e-07 | 11 | 45 |
| GSE33630 | advanced_vs_normal | DM1_like_score | +5.481 | 3.51e-07 | 11 | 45 |
| GSE33630 | advanced_vs_normal | THYROID_NONOVERLAP_score | −6.723 | 3.51e-07 | 11 | 45 |
| GSE33630 | advanced_vs_normal | TDS_like_score | −4.808 | 3.91e-07 | 11 | 45 |
| GSE33630 | advanced_vs_normal | TF_collapse_score | −4.808 | 3.91e-07 | 11 | 45 |
| GSE33630 | advanced_vs_normal | STAT3_AP1_DNMT_score | +1.655 | 2.42e-04 | 11 | 45 |
| GSE33630 | advanced_vs_normal | TACSTD2_z | +0.044 | 7.57e-01 | 11 | 45 |
| GSE29265 | advanced_vs_DTC | RAI_8_score | −1.565 | 1.96e-02 | 9 | 20 |
| GSE29265 | advanced_vs_DTC | DM1_like_score | +1.565 | 1.96e-02 | 9 | 20 |
| GSE29265 | advanced_vs_DTC | THYROID_NONOVERLAP_score | −1.515 | 1.17e-02 | 9 | 20 |
| GSE29265 | advanced_vs_DTC | TDS_like_score | −2.209 | 1.46e-03 | 9 | 20 |
| GSE29265 | advanced_vs_DTC | TF_collapse_score | −2.209 | 1.46e-03 | 9 | 20 |
| GSE29265 | advanced_vs_DTC | STAT3_AP1_DNMT_score | +1.195 | 1.33e-02 | 9 | 20 |
| GSE29265 | advanced_vs_DTC | TACSTD2_z | −2.062 | 3.11e-04 | 9 | 20 |
| GSE29265 | advanced_vs_normal | RAI_8_score | −2.294 | 7.50e-04 | 9 | 20 |
| GSE29265 | advanced_vs_normal | DM1_like_score | +2.294 | 7.50e-04 | 9 | 20 |
| GSE29265 | advanced_vs_normal | THYROID_NONOVERLAP_score | −2.663 | 4.55e-05 | 9 | 20 |
| GSE29265 | advanced_vs_normal | TDS_like_score | −2.630 | 1.48e-04 | 9 | 20 |
| GSE29265 | advanced_vs_normal | TF_collapse_score | −2.630 | 1.48e-04 | 9 | 20 |
| GSE29265 | advanced_vs_normal | STAT3_AP1_DNMT_score | +1.141 | 1.52e-02 | 9 | 20 |
| GSE29265 | advanced_vs_normal | TACSTD2_z | −0.861 | 8.89e-03 | 9 | 20 |
| GSE65144 | advanced_vs_normal | RAI_8_score | −2.792 | 4.01e-05 | 12 | 13 |
| GSE65144 | advanced_vs_normal | DM1_like_score | +2.792 | 4.01e-05 | 12 | 13 |
| GSE65144 | advanced_vs_normal | THYROID_NONOVERLAP_score | −3.531 | 3.17e-05 | 12 | 13 |
| GSE65144 | advanced_vs_normal | TDS_like_score | −2.934 | 3.17e-05 | 12 | 13 |
| GSE65144 | advanced_vs_normal | TF_collapse_score | −2.934 | 3.17e-05 | 12 | 13 |
| GSE65144 | advanced_vs_normal | STAT3_AP1_DNMT_score | +1.200 | 7.09e-03 | 12 | 13 |
| GSE65144 | advanced_vs_normal | TACSTD2_z | −0.044 | 9.78e-01 | 12 | 13 |
| GSE53157 | advanced_vs_DTC | RAI_8_score | −1.048 | 8.84e-02 | 5 | 19 |
| GSE53157 | advanced_vs_DTC | DM1_like_score | +1.048 | 8.84e-02 | 5 | 19 |
| GSE53157 | advanced_vs_DTC | THYROID_NONOVERLAP_score | −0.780 | 1.60e-01 | 5 | 19 |
| GSE53157 | advanced_vs_DTC | TDS_like_score | −0.984 | 1.20e-01 | 5 | 19 |
| GSE53157 | advanced_vs_DTC | TF_collapse_score | −0.984 | 1.20e-01 | 5 | 19 |
| GSE53157 | advanced_vs_DTC | STAT3_AP1_DNMT_score | −0.034 | 8.91e-01 | 5 | 19 |
| GSE53157 | advanced_vs_DTC | TACSTD2_z | −0.895 | 6.32e-02 | 5 | 19 |
| GSE53157 | advanced_vs_normal | RAI_8_score | −1.425 | 1.90e-01 | 5 | 2 |
| GSE53157 | advanced_vs_normal | DM1_like_score | +1.425 | 1.90e-01 | 5 | 2 |
| GSE53157 | advanced_vs_normal | THYROID_NONOVERLAP_score | −1.696 | 9.52e-02 | 5 | 2 |
| GSE53157 | advanced_vs_normal | TDS_like_score | −0.519 | 8.57e-01 | 5 | 2 |
| GSE53157 | advanced_vs_normal | TF_collapse_score | −0.519 | 8.57e-01 | 5 | 2 |
| GSE53157 | advanced_vs_normal | STAT3_AP1_DNMT_score | +3.317 | 9.52e-02 | 5 | 2 |
| GSE53157 | advanced_vs_normal | TACSTD2_z | −0.421 | 3.81e-01 | 5 | 2 |

### Direction-sign matrix (consensus across cohorts)

| Score | advanced_vs_DTC sign across 4 cohorts (33630 / 29265 / 65144 / 53157) | advanced_vs_normal sign across 3 cohorts with normal arm (33630 / 29265 / 65144) |
|-------|----|----|
| RAI_8 | **−** , **−** , (no DTC) , − | **−** , **−** , **−** |
| DM1_like | **+** , **+** , (no DTC) , + | **+** , **+** , **+** |
| THYROID_NONOVERLAP | **−** , **−** , (no DTC) , − | **−** , **−** , **−** |
| TDS_TF / TF_collapse | **−** , **−** , (no DTC) , − | **−** , **−** , **−** |
| STAT3_AP1_DNMT | **+** , **+** , (no DTC) , − (n.s.) | **+** , **+** , **+** |
| TACSTD2_z | **−** , **−** , (no DTC) , − | + (n.s.) , **−** , ≈0 |

Lineage panels are **direction-consistent in every cohort where the contrast is testable.** STAT3/AP1/DNMT trends positive in 6/7 testable contrasts. TACSTD2 is **heterogeneous** (mostly negative — i.e. lower in advanced than DTC) and is **not claimed** as a positive marker for advanced disease in this report.

---

## 7. Spearman correlations — FULL table

`external_spearman.tsv`. Within-dataset, all non-DROP samples.

| dataset | x | y | rho | p | n |
|---|---|---|---|---|---|
| GSE33630 | DM1_like_score | THYROID_NONOVERLAP_score | **−0.9273** | 9.11e-46 | 105 |
| GSE33630 | TF_collapse_score | DM1_like_score | **−0.9050** | 5.12e-40 | 105 |
| GSE33630 | STAT3_AP1_DNMT_score | DM1_like_score | +0.4362 | 3.29e-06 | 105 |
| GSE33630 | TACSTD2_z | DM1_like_score | +0.4215 | 7.53e-06 | 105 |
| GSE29265 | DM1_like_score | THYROID_NONOVERLAP_score | **−0.8524** | 7.80e-15 | 49 |
| GSE29265 | TF_collapse_score | DM1_like_score | **−0.8845** | 3.57e-17 | 49 |
| GSE29265 | STAT3_AP1_DNMT_score | DM1_like_score | +0.3100 | 3.02e-02 | 49 |
| GSE29265 | TACSTD2_z | DM1_like_score | +0.1266 | 3.86e-01 | 49 |
| GSE65144 | DM1_like_score | THYROID_NONOVERLAP_score | **−0.9354** | 7.24e-12 | 25 |
| GSE65144 | TF_collapse_score | DM1_like_score | **−0.9831** | 1.86e-18 | 25 |
| GSE65144 | STAT3_AP1_DNMT_score | DM1_like_score | +0.5885 | 1.97e-03 | 25 |
| GSE65144 | TACSTD2_z | DM1_like_score | −0.2008 | 3.36e-01 | 25 |
| GSE53157 | DM1_like_score | THYROID_NONOVERLAP_score | **−0.8386** | 8.80e-08 | 26 |
| GSE53157 | TF_collapse_score | DM1_like_score | **−0.7443** | 1.31e-05 | 26 |
| GSE53157 | STAT3_AP1_DNMT_score | DM1_like_score | +0.4680 | 1.59e-02 | 26 |
| GSE53157 | TACSTD2_z | DM1_like_score | +0.4147 | 3.52e-02 | 26 |

### Headline interpretation

- **DM1_like ↔ THYROID_NONOVERLAP** anti-correlation (ρ ≤ −0.84 in every cohort) is the central within-cohort claim: the two panels share **no genes** by design, so this is direction-consistent lineage silencing rather than a panel-overlap artifact.
- **TF_collapse ↔ DM1_like** ρ ≈ −0.74 to −0.98 → lineage TF collapse tracks the DM1 axis.
- **STAT3_AP1_DNMT ↑ with DM1_like** (ρ +0.31 to +0.59) → consistent with the proposed mechanism arm.
- **TACSTD2 vs DM1_like** is **not directionally consistent across cohorts** (positive in 33630 / 29265 / 53157, negative in 65144); not claimed as a coupled marker.

---

## 8. Per-dataset roles for downstream use

| Dataset | Role | Reason |
|---------|------|--------|
| GSE33630 | **Main external replication panel** | Largest (n=105), full ATC + PTC + normal triad, strongest effect sizes |
| GSE65144 | **Main supporting** (ATC vs normal) | Cleanest 12 vs 13 ATC/normal contrast, RAI_8 d=−2.79 |
| GSE29265 | **Supp** (Sporadic vs Chernobyl) | Adds paired-normal architecture; direction-consistent with smaller effect sizes |
| GSE53157 | **Supp / sensitivity** | PDTC arm n=5 — direction-consistent only; not main |
| GSE126698 | **DROP** | Series Matrix metadata-only; alignment forbidden |
| GSE53072 / GSE60542 / GSE120177 | **HOLD** | per user spec — non-GPL570 / lymphoid confounding / Paper 9 scope |

---

## 9. Risks / caveats (do not paper over)

- **Microarray vs RNA-seq:** all 4 used cohorts are GPL570 Affy U133+2 RMA values; the in-paper 8-gene index was derived on TCGA RNA-seq. Cross-platform agreement on direction is what is claimed; **absolute coefficients are not transferable** and were not transferred (per-dataset z-scoring used instead of pooling).
- **Probe → gene collapse** is mean-of-probes; legacy alias table applied. No probe-set-vs-symbol curation; coverage = whatever GEO's Aug-2016 GPL570 annot.gz reports.
- **Histology label quality:** taken verbatim from `Sample_source_name_ch1` / `Sample_title` / `Sample_characteristics_ch1`. No re-review of pathology. GSE53157's commercial RNA pool dropped; otherwise no samples re-classified. GSE33630's 45 "patient-matched non-tumor control" samples are treated as `normal`; some literature flags such adjacent normals as not histologically pristine — caveat carried forward.
- **No survival or mutation metadata** fetched or used. **No** survival / clinical-outcome / mutation / therapeutic-response claim is made or implied.
- **Cross-study batch effects:** datasets **never pooled** before scoring. Effect sizes and tests are **within-dataset**. The "direction consistency" forest is a consensus *over independent studies*, not a meta-analytic pooled estimate.
- **GSE53157 PDTC underpower:** n=5 PDTC vs n=7 PTC; no contrast in this cohort survives multiple-testing correction. Reported as supportive only.
- **GSE126698:** spec-conformant SRA realignment would change the access list; **not undertaken**.

---

## 10. Wording rules for any downstream prose

**Safe:**
- "external expression validation"
- "direction-consistent lineage-silencing"
- "advanced-disease replication"

**Forbidden** (not supported by anything in this sweep):
- "clinical validation"
- "survival validation"
- "fusion validation"
- "progression proven"
- "therapeutic response proven"

---

## 11. Figures (Step 8)

| File | Content |
|------|---------|
| `external_rai_lineage_boxplots.png` | 3 × 4 grid: RAI_8 / DM1_like / THYROID_NONOVERLAP scores per histology, per dataset; coloured boxes by group, jittered points |
| `external_dm1_nonoverlap_scatter_grid.png` | 2 × 2 scatter: DM1_like vs THYROID_NONOVERLAP per dataset, ρ + p in title, points coloured by histology |
| `external_direction_consistency_forest.png` | Forest of Cohen's d across (dataset × contrast × score); opaque points = p<0.05 |
| `external_dataset_qc_heatmap.png` | Gene-panel availability heatmap (rows = panel│gene, columns = datasets) |

---

## 12. Output inventory

```
project/results/p_external_expression_validation/
├── sample_metadata.tsv
├── external_gene_coverage.tsv
├── external_score_tests.tsv
├── external_direction_consistency.tsv
├── external_spearman.tsv
├── external_sample_scores.tsv.gz
├── GSE33630_expression_gene_log.tsv.gz
├── GSE29265_expression_gene_log.tsv.gz
├── GSE65144_expression_gene_log.tsv.gz
├── GSE53157_expression_gene_log.tsv.gz
├── external_rai_lineage_boxplots.png
├── external_dm1_nonoverlap_scatter_grid.png
├── external_direction_consistency_forest.png
├── external_dataset_qc_heatmap.png
├── scripts/
│   ├── run_external_validation.py
│   └── make_figures.py
└── raw/                              # gitignored
    ├── GSE33630/GSE33630_series_matrix.txt.gz
    ├── GSE29265/GSE29265_series_matrix.txt.gz
    ├── GSE65144/GSE65144_series_matrix.txt.gz
    ├── GSE53157/GSE53157_series_matrix.txt.gz
    ├── GSE126698/GSE126698_series_matrix.txt.gz   # metadata-only, kept for audit
    └── GPL570/GPL570.annot.gz

project/reports/
├── 2026_05_04_external_expression_access_check.md
├── 2026_05_04_external_expression_validation_report.md
├── 2026_05_04_external_expression_commit_proposal.md
└── 2026_05_04_external_expression_FULL_RESULTS.md   # CANONICAL single-md report
```

---

## 13. Commit proposal (NOT executed)

### Working tree footprint

- `M  .gitignore` — adds `project/results/p_external_expression_validation/raw/`
- `??` 4 reports under `project/reports/2026_05_04_external_expression_*.md`
- `??` `project/results/p_external_expression_validation/` (whole subtree)
  - `raw/` (62 MB) — gitignored
  - tracked artifacts: ~30 MB (4 gene-level matrices) + ~70 KB (tables) + ~1.1 MB (4 PNGs) + scripts

### Proposed commit groups

**Commit 1 — canonical report, scripts, gitignore**
```
.gitignore
project/reports/2026_05_04_external_expression_FULL_RESULTS.md
project/results/p_external_expression_validation/scripts/run_external_validation.py
project/results/p_external_expression_validation/scripts/make_figures.py
```
> docs+infra: external expression validation sweep — canonical single report + scripts (4 GPL570 cohorts; GSE126698 dropped at access)

Do **not** include the intermediate md fragments in the main commit unless you want the full scratch audit trail:
```
project/reports/2026_05_04_external_expression_access_check.md
project/reports/2026_05_04_external_expression_validation_report.md
project/reports/2026_05_04_external_expression_commit_proposal.md
```

**Commit 2 — processed matrices, scores, summary tables (~30 MB .tsv.gz)**
```
project/results/p_external_expression_validation/sample_metadata.tsv
project/results/p_external_expression_validation/external_gene_coverage.tsv
project/results/p_external_expression_validation/external_score_tests.tsv
project/results/p_external_expression_validation/external_direction_consistency.tsv
project/results/p_external_expression_validation/external_spearman.tsv
project/results/p_external_expression_validation/external_sample_scores.tsv.gz
project/results/p_external_expression_validation/GSE33630_expression_gene_log.tsv.gz
project/results/p_external_expression_validation/GSE29265_expression_gene_log.tsv.gz
project/results/p_external_expression_validation/GSE65144_expression_gene_log.tsv.gz
project/results/p_external_expression_validation/GSE53157_expression_gene_log.tsv.gz
```
> data: external expression validation — gene-level expression matrices + score tests (4 GPL570 cohorts, n=205)

Alternative: keep only score tables tracked (~70 KB) and gitignore the 4 `*_expression_gene_log.tsv.gz` (regenerable from `raw/` + scripts). User decision.

**Commit 3 — figures (~1.1 MB)**
```
project/results/p_external_expression_validation/external_rai_lineage_boxplots.png
project/results/p_external_expression_validation/external_dm1_nonoverlap_scatter_grid.png
project/results/p_external_expression_validation/external_direction_consistency_forest.png
project/results/p_external_expression_validation/external_dataset_qc_heatmap.png
```
> figs: external expression validation — boxplots / scatter grid / direction-consistency forest / coverage QC

**NOT included in any commit:** `project/results/p_external_expression_validation/raw/` (62 MB Series Matrix + GPL570 annotation, reproducible from public GEO URLs; gitignored).

---

External expression validation sweep complete. No datasets outside allowed list acquired.
