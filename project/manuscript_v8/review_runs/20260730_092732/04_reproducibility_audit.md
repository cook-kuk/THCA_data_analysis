---
title: Reproducibility audit
audit_date: 2026-07-30
auditor: manuscript-audit-agent (reproducibility auditor mode)
---

# 04 — Reproducibility audit

## 1. Cohort accession numbers

### Present and complete

| Cohort | Accession stated | Location |
|---|---|---|
| TCGA-THCA RNA-seq, WGS, HM450, clinical | dbGaP phs000178; cBioPortal `thca_tcga_pub` | STAR Methods Key Resources Table |
| MSK-IMPACT thyroid (Landa 2016) | cBioPortal `thca_mskcc_2016` | STAR Methods |
| K2 / PRJEB11591 (Yoo 2016) | ENA PRJEB11591 | STAR Methods |
| Lee Korean PTC | GEO GSE213647 | STAR Methods and v2 Methods |
| Pu 2021 single-cell | GEO GSE184362 | STAR Methods and v2 Methods |
| Lu 2023 single-cell | GEO GSE193581 | STAR Methods and v2 Methods |
| GSE241184 (Phase 1 single-cell) | GEO GSE241184 | STAR Methods |
| GSE151179 (post-RAI refractory) | listed in v2 text | v2 Methods line 160 |
| GSE286332 (Korean PTC reference) | GEO GSE286332 | STAR Methods |

### Present in Results text but missing from v2 Methods accession listing

| Cohort | Accession | Issue |
|---|---|---|
| GPL570 microarray cohorts | GSE29265, GSE33630, GSE65144 listed in v2 Methods (line 160) | Only 3 accessions listed but text says "four GPL570 microarray cohorts" (v2 line 88 and line 148). A fourth GPL570 cohort accession is missing. |
| GSE126698 | Listed in 04_results.md and figure captions (09_reviewer_qa.md Q14) | Not listed in v2 Methods cohort section; missing from cohort inventory |
| GSE232237 | Listed in v2 Methods (line 160) "additional single-cell support" | Accession present but no per-cohort n or details |
| Mun 2025 proteogenomic cohort (n=336) | Not stated in v2 Methods; only cited as "Mun 2025" | No public accession number given |
| Landa 2016 dataset (GSE76039) | Referenced in STAR Methods Key Resources and figure captions | GSE76039 accession implied by "Landa 2016" but not explicitly stated in v2 Methods |

---

## 2. Clustering seed, k, and preprocessing

### Specified

- k=2 stated explicitly: "KMeans clustering with k = 2 was applied to the eight-dimensional panel expression matrix" (v2 Methods line 168)
- Within-cohort z-standardization stated: "log-transformed and z-standardized within cohort" (v2 Methods line 168)
- Python scikit-learn KMeans used: stated in STAR Methods (line 32)
- scikit-learn "default initialization" stated in STAR Methods (line 84) — this is k-means++ by default

### Not specified

- **Random seed:** The clustering random seed is not stated. scikit-learn KMeans default uses `random_state=None`, which means results are not reproducible across runs. The manuscript must specify `random_state=42` or equivalent and document it.
- **Number of initializations (n_init):** Not stated. scikit-learn default is n_init=10; this should be documented.
- **Silhouette threshold:** STAR Methods (line 84) says "confirmed by silhouette analysis (mean silhouette score > 0.4 within both clusters)" — this is a post-hoc confirmation threshold, not a pre-specified criterion. The observed silhouette scores per cluster are not reported.
- **Logistic regression for P(DM1):** STAR Methods (line 84) states cross-validated AUC=0.968, but the logistic regression hyperparameters (C, solver, penalty) are not specified.

---

## 3. Meta-analysis method specification

**Specified:** DerSimonian-Laird random-effects (STAR Methods line 147; v2 Methods line 204; `statsmodels.stats.meta_analysis` stated in STAR Methods).

**Not specified:**
- Software version for statsmodels (only scikit-learn 1.5 and lifelines 0.27 are explicitly version-listed; statsmodels version missing from STAR Methods environment table)
- Exact function call for DerSimonian-Laird in statsmodels (the function name varies across versions)
- Whether log-HR and SE were derived from Cox model output or a separate function

---

## 4. ATA risk-tier mapping

The v2 draft and figure captions reference ATA 2015/2025 risk-tier overlay (Fig. 6B in NC caption build, Fig. 5D in reviewer-first plan). The DM1_PROJECT_STATE.md (line 75) lists "ATA intermediate-risk mosaic or replacement panel, only if the exact risk definitions and denominators are recoverable" as a remaining computational task.

**Issue:** The ATA 2015 risk-tier mapping rules (low/intermediate/high) are not defined in the Methods. The specific clinical variables used (size, lymph node status, extrathyroidal extension, distant metastasis, histology) and the TCGA field names used to operationalize these rules are not stated. This is required for reproducibility.

---

## 5. Figure source data files

The NC checklist requires source data for each main figure. The manuscript text (v2 Methods, Data Availability section, line 216) states: "Processed per-sample score tables, methylation summaries, structural-variant joins, interaction Cox model outputs and Monte Carlo simulation seeds will be deposited with the public code release at submission."

**Assessment:** Source data deposition is committed but not yet completed. No source-data manifest mapping each figure panel to a specific file currently exists in the manuscript. A figure-to-file mapping is required per NC policy.

**Files available in `project/results/manuscript_v8_nc_main/`:**
- `master_tcga.tsv` — appears to be master TCGA table
- `cohort_scorecard.tsv` — likely score summary
- Individual figure scripts exist (fig1_discovery.py through fig6_reflex.py)
- Fig 7 and Fig 8 scripts: `deep_analysis_5x_parallel.py`, `ihc_3plex_killer_figure.py`, `gene_reduction_v2.py`, `a2_replication_lee_k2.py`

**Missing:** A manifest document (Supplementary Table format or README) that maps each manuscript figure panel (e.g., Fig. 1c = UMAP from master_tcga.tsv, Fig. 3b = methylation from r5_2_sample_methylation_8gene.tsv) to its source file and generating script.

---

## 6. Code repository and DOI

**Stated intention:** Code will be archived with Zenodo DOI at submission.

**Current status:**
- STAR Methods Key Resources Table (line 39): "[TODO: insert public code repository URL on submission]; [TODO: insert Zenodo DOI on submission]"
- STAR Methods Data and Code Availability (line 49): same two TODO placeholders
- No public repository URL currently exists in the manuscript

**Blockers:** Two TODO placeholders remain in STAR Methods. These must be filled before submission. The NC checklist explicitly prohibits "TBD," "XXX," "<DOI>," or private local paths in the final package.

---

## 7. All TODO / TBD / PLACEHOLDER found in manuscript files

### In `NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md`

| Line | Text |
|---|---|
| 36 | `[VERIFY institutional corresponding email]` |
| 44 | `[AUTHOR HOOK SLOT — protected voice section...]` |
| 52 | `[AUTHOR AIM SLOT — protected voice section...]` |
| 134 | `[AUTHOR DISCUSSION 3.1 SLOT — protected voice section...]` |
| 154 | `[AUTHOR LIMITATIONS SLOT — protected voice section...]` |
| 224 | `[AUTHOR SLOT — add funding, institutional support...]` |

### In `07_star_methods.md`

| Line | Text |
|---|---|
| 39 | `[TODO: insert public code repository URL on submission]` |
| 39 | `[TODO: insert Zenodo DOI on submission]` |
| 49 | `[TODO: insert public code repository URL on submission]` (repeated) |
| 49 | `[TODO: insert Zenodo DOI on submission]` (repeated) |

### In `08_cover_letter.md` (not audited in full, but known from VOICE_PROTECTED_SLOTS.md)

- Cover letter paragraph 1 voice slot unfilled

### In `09_reviewer_qa.md`

- Q9 prose is voice-protected and not finalized

---

## 8. Environment specification completeness

STAR Methods (line 104) states: "All analyses were performed in Python 3.11 with scikit-learn 1.5, pandas 2.2, numpy 1.26, scipy 1.13, lifelines 0.27, and pyDESeq2 0.4."

**Missing versions:**
- `statsmodels` (used for DerSimonian-Laird meta-analysis)
- `matplotlib` and `seaborn` (used for figures)
- STAR aligner version (listed in Key Resources but version not specified)
- kallisto version (listed in Key Resources but version not specified)

**Missing specification:**
- GENCODE v44 release date / chromosome build (GRCh38 assumed but not stated)
- HM450 probe-to-gene aggregation method (max promoter beta? mean over TSS-200 and 1stExon probes? Not specified)
- cBioPortal SV endpoint URL and access date

---

## 9. Critical reproducibility gaps summary

| Gap | Severity |
|---|---|
| Random seed for KMeans clustering not specified | MAJOR |
| Fourth GPL570 cohort accession missing | MAJOR |
| GSE126698 not listed in v2 Methods cohort section | MAJOR |
| Mun 2025 proteogenomic accession not stated | MAJOR |
| 2× TODO placeholders for Zenodo DOI and repo URL | MAJOR |
| ATA risk-tier operationalization rules absent | MAJOR |
| HM450 probe aggregation method not specified | MAJOR |
| No figure-to-source-data manifest | MAJOR |
| statsmodels version missing | MINOR |
| Fig 3C not yet built (missing figure) | MAJOR |
| Fig 3D noted as needing rebuild | MAJOR |
