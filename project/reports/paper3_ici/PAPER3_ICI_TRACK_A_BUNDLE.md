<!-- FROZEN 2026-05-04 — design bundle locked. Track B blocked until Paper 1 bioRxiv + Paper 2 Task A/B/C closure + explicit user command. Do not edit. -->


# Paper 3 ICI — Track A Bundle

**Working title:** HLA loss, neoantigen architecture, and immune ecotypes define ICI vulnerability in molecularly dark thyroid cancer
**Author:** Seungho Cook
**Track:** A — design only. No data download, no analysis execution.
**Bundle date:** 2026-05-04

---

## Table of contents

1. [Dataset Registry](#1-dataset-registry)
2. [Signature Registry](#2-signature-registry)
3. [scRNA Reference Atlas Plan](#3-scrna-reference-atlas-plan)
4. [HLA & Neoantigen Feasibility](#4-hla--neoantigen-feasibility)
5. [DIAL Audit Plan](#5-dial-audit-plan)
6. [Figure Plan](#6-figure-plan)
7. [Go / No-Go Verdict](#7-go--no-go-verdict)
8. [12-Week Execution Plan](#8-12-week-execution-plan)

---

## 1. Dataset Registry

**Working title:** HLA loss, neoantigen architecture, and immune ecotypes define ICI vulnerability in molecularly dark thyroid cancer
**Track:** A — design only. No download, no analysis. All accessions marked `to_verify` until Track B kickoff.
**Author:** Seungho Cook
**Status date:** 2026-05-04

**Claim guard:** Datasets are catalogued for ICI vulnerability / immunogenomic prioritization, never for "ICI response prediction" in thyroid unless thyroid ICI-treated raw RNA-seq is acquired.

---

## 0. Verification policy

For each accession:
- **`to_verify`** — accession string is from prior literature / user input; existence, version, file types, and access conditions must be checked at Track B kickoff via GEO/SRA/EGA/dbGaP web UI (not by Claude during Track A).
- **`internal`** — already on local mirror or processed for Paper 1/2; reuse read-only with paper-boundary discipline.
- **`controlled`** — known to be controlled-access (dbGaP / EGA / institutional). Application timing must be added to the 12-week plan.
- **`tier_x`** — see priority tiers below.

Priority tiers:
- **Tier 1** — required for headline figure. Track B blocks if unavailable.
- **Tier 2** — replication / cross-platform. Track B proceeds if at least 50% of Tier 2 are usable.
- **Tier 3** — supplementary / sensitivity. Loss tolerable.

---

## 1. Bulk RNA-seq ± WES — discovery & dedifferentiation cohorts

| Accession | Modality | n (claimed) | Subtype | Role | Tier | Status | Notes |
|---|---|---|---|---|---|---|---|
| TCGA-THCA | RNA-seq + WES + clinical | ~500 | PTC, FVPTC, TCV; ATC limited | Primary discovery: ecotype NMF, HLA LOH, neoantigen, dark-matter stratification | 1 | internal | Driver/fusion calls already curated in Paper 1 ETL. WES BAM available via GDC; LOHHLA needs raw BAM not just MAF. |
| GSE76039 | RNA-seq | ~37 (PDTC + ATC; Landa 2016) | PDTC, ATC | Dedifferentiation axis anchor | 1 | to_verify | Landa 2016 JCI (cite already saved in `v17_landa2016_cite_save`). Small n — pair with TCGA aggressive subset for power. |
| GSE33630 | microarray (Affy U133 Plus 2) | ~105 | PTC, ATC, normal | Cross-platform replication | 2 | to_verify | Microarray, so signature scoring only — no neoantigen, no HLA LOH. |
| GSE65144 | microarray | ~13 | ATC + normal | ATC enrichment | 2 | to_verify | Tiny n; suitable as ATC-vs-normal contrast in figure inset. |
| GSE54958 | microarray | small | PTC stages | Stage-axis sensitivity | 3 | to_verify | |
| GSE53157 | microarray | ~27 | PTC | Replication | 3 | to_verify | |
| GSE29265 | microarray | ~49 | PTC, ATC, normal | Cross-platform | 2 | to_verify | |
| Yoo SK et al. 2019 *Nat Commun* | RNA-seq | Korean ATC + advanced DTC | ATC, advanced DTC | Asian dedifferentiation generalization | 1 | to_verify | Accession not yet pinned; check supp; may be EGA-controlled. If controlled, add EGA application to Wk 1. |
| PRJKA210106 | RNA-seq | n=282 Korean PTC | PTC | Asian baseline | 2 | to_verify | Cross-paper care — already touched by Paper 1/2 lineage. Use for HLA-typing reference, not for ecotype discovery, to keep Paper 3 distinct. |
| PRJEB11591 (Yoo 2016 SNU-GMI) | RNA-seq | n=260 Korean PTC | PTC | HLA-typing reference (arcasHLA already done) | 3 | internal (read-only) | Per `v17_arcasHLA_korean_k2`. Do not re-run; do not duplicate Paper 1/2 analyses. Use only as Korean HLA-frequency anchor. |

---

## 2. Single-cell RNA-seq — atlas construction

| Accession | Platform | n samples (claimed) | Subtype focus | Role | Tier | Status | Notes |
|---|---|---|---|---|---|---|---|
| GSE184362 | 10x 3' | thyroid tumor + normal | PTC + ATC | Primary atlas — epithelial + immune + stromal | 1 | to_verify | Per literature, ATC/PDTC content. Verify n PDTC/ATC vs PTC. |
| GSE193581 | 10x | PTC | PTC | T cell / myeloid sub-states | 1 | to_verify | |
| GSE232237 | 10x | thyroid | mixed | Atlas integration | 2 | to_verify | |
| GSE191288 | 10x | thyroid | PTC | Replication | 2 | to_verify | |
| GSE148673 | 10x | thyroid | mixed | Replication | 3 | to_verify | |
| Han et al. 2024 *JCI Insight* | 10x | thyroid | dedifferentiation | Atlas | 2 | to_verify | Accession to be pulled from supp at Track B kickoff. |

**Atlas-side risks:** ATC/PDTC scRNA samples are scarce. Combined ATC + PDTC scRNA cells across all six accessions are likely <5K — viable for state annotation by reference projection but not for de novo trajectory inference. The atlas plan in `paper3_ici_scRNA_reference_atlas_summary.md` accounts for this.

---

## 3. Spatial transcriptomics — TLS / niche localization

| Accession | Platform | n (claimed) | Subtype | Role | Tier | Status | Notes |
|---|---|---|---|---|---|---|---|
| GSE250521 | Visium | DM1-spatial validation set | PTC + DM1 | TLS / B-cell niche localization (already used in Paper 1/2 supplement) | 2 | internal (read-only) | Paper-boundary care. Paper 3 may reference niche maps but must NOT duplicate Paper 1's DM1 spatial analysis. |
| GSE230424 | Visium | thyroid | mixed | Spatial replication | 3 | to_verify | |
| GSE248205 | Visium | thyroid | mixed | Spatial replication | 3 | to_verify | |
| Ning / Liao / Zheng spatial | Visium / VisiumHD | thyroid | mixed | Asian spatial | 3 | to_verify | Accessions to be located via PubMed → supp at Track B kickoff. |

---

## 4. Pan-cancer ICI reference — DIAL audit + transfer

ICI-treated raw RNA-seq cohorts used for direction-invariance testing of the signature suite. Thyroid ICI-treated raw RNA-seq is *not* expected to be available, hence the transfer-only framing.

| Accession | Disease | n (claimed) | ICI agent | Role | Tier | Status | Notes |
|---|---|---|---|---|---|---|---|
| IMvigor210 (EGAS00001002556 / SRP172670) | UC | ~298 | atezolizumab | DIAL anchor 1 | 1 | controlled / partial-public | Pre-treatment RNA-seq + RECIST. Public package available via `IMvigor210CoreBiologies` R data drop. |
| Hugo GSE78220 | melanoma | ~28 | anti-PD-1 | DIAL anchor 2 | 1 | to_verify | Small but heavily reused in pan-cancer ICI papers; pre-treatment. |
| Riaz GSE91061 | melanoma | ~109 | nivolumab | DIAL anchor 3 | 1 | to_verify | Pre + on-treatment; use pre only for DIAL. |
| Gide melanoma (PRJEB23709) | melanoma | ~73 | anti-PD-1 ± anti-CTLA-4 | DIAL anchor 4 | 2 | to_verify | |
| Liu melanoma (dbGaP phs000452) | melanoma | ~121 | anti-PD-1 | DIAL anchor 5 | 2 | controlled | dbGaP application required. If timeline blocks, drop to Tier 3. |
| Kim gastric (PRJEB25780) | GC | ~45 | pembrolizumab | DIAL extension | 3 | to_verify | Pan-cancer breadth. |
| Cho NSCLC | NSCLC | varies | anti-PD-1 / anti-PD-L1 | DIAL extension | 3 | to_verify | Specific accession to be pinned at Track B kickoff. |

**Thyroid ICI gap (key risk).** No ICI-treated thyroid raw RNA-seq cohort is reliably public as of design date. Implications:
- Headline claim must be "ICI vulnerability / readiness", not "response predictor".
- Validation strategy = pan-cancer transfer (DIAL audit identifies which signatures retain direction in non-thyroid → apply only those to thyroid).
- A single-arm thyroid ICI report (small, e.g. lenvatinib + pembrolizumab in ATC) may exist as supplementary tables only; treat as case-level anecdote, not validation.

---

## 5. Internal / Paper-1 / Paper-2 reuse — boundary discipline

Reusable read-only artifacts (no re-analysis allowed during Track B unless explicitly opened):

| Path | Origin | Allowed reuse |
|---|---|---|
| `project/results/00_qc/` | Paper 1 ETL | Sample-level QC tables — read-only reference. |
| `project/results/01_spatial_score/` | Paper 1/2 spatial pipeline | Spatial niche maps for TLS visualization only. Do not recompute. |
| `project/results/03_pathology_poc/` | Paper 2 H&E projection | Pathology features as comparator — *not* as Paper 3 input. |
| PRJEB11591 arcasHLA results | per `v17_arcasHLA_korean_k2` | Korean HLA frequency anchor. Do not re-run. |
| TCGA-THCA driver/fusion calls | Paper 1 ETL | Dark-matter (BRAF/RAS-negative, fusion-negative) stratification. |

**Boundary rules (binding for Paper 3 prose):**
- Paper 1 (DM1 molecular) — Paper 3 may *cite* the dark-matter definition. May not duplicate DM1 cluster analyses.
- Paper 2 (H&E-DM1 / HT-overlap PTC) — Paper 3 may *cite* the immune-rich PTC+HT phenotype as one comparator subgroup. May not duplicate TLS/AICDA/BCR analyses.
- Paper 4 (GD HLA backlog) — Paper 3 is cancer-side; reference at most at the population HLA-frequency level.

---

## 6. Verification checklist (to execute at Track B kickoff, NOT now)

1. For each `to_verify` accession: open GEO/SRA/EGA/ENA UI; confirm (a) records exist, (b) raw FASTQ or processed counts available, (c) clinical/sample metadata accessible, (d) license/access mode.
2. For each `controlled` accession: identify required application (dbGaP, EGA, institutional), determine 4–8 week lag, place in 12-week plan accordingly.
3. Resolve Yoo SK 2019 Korean ATC accession; if EGA-controlled, decide whether to apply or drop to Tier 2 with TCGA aggressive subset substitution.
4. Resolve Han 2024 JCI Insight accession from paper supp.
5. Resolve Ning / Liao / Zheng spatial accessions.
6. Resolve Cho NSCLC ICI accession.

Verification must be logged into `paper3_ici_dataset_registry_VERIFIED.tsv` (TSV form) at kickoff.

---

## 7. Out-of-scope datasets (explicit exclusions)

- ICI-treated cohorts in non-relevant diseases (e.g., RCC, CRC) unless DIAL audit specifically benefits — added only via amendment.
- Non-thyroid scRNA atlases used as reference (e.g., pan-cancer T-cell atlases like Zheng 2021, Andreatta 2022) — these are *reference annotations*, not data inputs, and will be cited as method dependencies in `paper3_ici_scRNA_reference_atlas_summary.md`.
- Animal models — out of scope for this paper.
- Patient-derived organoids without paired tumor RNA-seq — out of scope.

---

Track A completed. No marathon violation.

---

## 2. Signature Registry

**Working title:** HLA loss, neoantigen architecture, and immune ecotypes define ICI vulnerability in molecularly dark thyroid cancer
**Track:** A — design only. No scoring runs. Gene members are listed for reference; final operational gene lists will be locked at Track B kickoff after intersection with each cohort's expression matrix.
**Author:** Seungho Cook
**Status date:** 2026-05-04

**Claim guard:** Signatures are computed for ecotype discovery and DIAL audit, not as standalone "ICI response classifiers". Interpretation in thyroid is conditional on DIAL audit results (`paper3_ici_DIAL_audit_plan.md`).

---

## 0. Scoring conventions (locked at design)

- **Bulk scoring method:** ssGSEA (GSVA package) as default. Singscore as sensitivity (rank-based, sample-independent). Both reported in Supp.
- **Single-cell scoring method:** UCell (rank-based, robust to dropout) as default. AUCell as sensitivity.
- **Normalization input:** log2(TPM+1) for bulk RNA-seq; log-normalized counts (Seurat / scanpy default) for scRNA. Microarray cohorts use RMA-normalized log-intensity.
- **Score sign convention:** higher = more of the named feature. Direction-flips revealed by DIAL audit are documented per signature, never silently corrected.
- **Cohort-level z-score:** all signatures z-scored within cohort before cross-cohort comparison, to avoid platform-mean confounding. Pan-cohort scaling only applied after batch evaluation.
- **Validation across platforms:** any signature dropped if <50% gene members measurable on a microarray cohort; flagged in `paper3_ici_signature_registry_VERIFIED.tsv` at kickoff.

---

## 1. T-cell inflammation / IFNγ axis

| Signature ID | Source | Members (representative) | Role | Notes |
|---|---|---|---|---|
| `IFNG_AYERS6` | Ayers et al. *JCI* 2017 | IFNG, CXCL9, CXCL10, IDO1, HLA-DRA, STAT1 | Core IFNγ axis | Used in pembrolizumab IFNγ-signature studies. |
| `IFNG_AYERS18` (TIS, Tumor Inflammation Signature) | Ayers 2017 | 18 genes incl. CD8A, GZMK, CD27, LAG3, TIGIT, CXCL9, etc. | T-cell-inflamed phenotype | Dako commercial as nCounter; we re-implement as ssGSEA. |
| `CYTOLYTIC_ROONEY` | Rooney et al. *Cell* 2015 | GZMA, PRF1 (geometric mean) | Cytolytic activity | 2-gene; very platform-stable. |
| `EFFECTOR_T` | Tirosh / Sade-Feldman composite | CD8A, CD8B, GZMB, GZMK, IFNG, NKG7, PRF1 | Effector CD8 | Bulk + scRNA. |
| `EXHAUSTED_T` | Tirosh / Wherry composite | PDCD1, CTLA4, LAG3, TIGIT, HAVCR2, TOX, TOX2, EOMES | Exhausted CD8 | Co-occurs with IFNG axis; DIAL-flag candidate. |

---

## 2. Antigen presentation machinery

| Signature ID | Source | Members | Role | Notes |
|---|---|---|---|---|
| `MHC1_CORE` | composite | HLA-A, HLA-B, HLA-C, B2M, TAP1, TAP2, TAPBP, NLRC5, ERAP1 | Class I machinery | Loss = immune escape. |
| `MHC2_CORE` | composite | HLA-DRA, HLA-DRB1, HLA-DPA1, HLA-DPB1, HLA-DQA1, HLA-DQB1, CIITA | Class II machinery | Tumor-intrinsic Class II expression in PTC is documented (HT-overlap signal in Paper 2). Important DIAL candidate. |
| `IFNG_RESPONSE_HALLMARK` | MSigDB Hallmark | hallmark IFNG_RESPONSE | Upstream context | Used as IFNγ-response sanity check. |

---

## 3. TLS / B-cell / plasma niche

| Signature ID | Source | Members | Role | Notes |
|---|---|---|---|---|
| `TLS_CABRITA9` | Cabrita et al. *Nature* 2020 | CCL19, CCL21, CXCL13, CCR7, CD79B, MS4A1, BCL6, LAMP3, SELL | TLS hallmark | Core TLS signature. |
| `TLS_MEYLAN` | Meylan et al. *Immunity* 2022 | extended TLS module incl. germinal-center markers | TLS extended | Cross-validate Cabrita. |
| `B_CELL_CORE` | composite | CD19, MS4A1, CD79A, CD79B, BANK1 | B-cell mass | Distinguish TLS-driving from infiltrating B-cells. |
| `PLASMA_CELL` | composite | MZB1, JCHAIN, XBP1, IGHA1, IGHG1 | Plasma niche | Antibody-class info available only when IGH constant chains expressed. |
| `CXCL13_AXIS` | composite | CXCL13, CXCR5 | CXCL13 axis | Single-gene CXCL13 also reported separately. |

**Paper 2 boundary note.** AICDA / IGHV clonality / detailed BCR architecture are *Paper 2* analyses (HT-overlap PTC). Paper 3 limits itself to TLS-presence scoring and CXCL13 axis; deep BCR repertoire is *not* re-analyzed here.

---

## 4. Myeloid suppression / TAM

| Signature ID | Source | Members | Role | Notes |
|---|---|---|---|---|
| `MYELOID_CORE` | composite | CD68, CD163, CD14, CSF1R, MRC1 | Total myeloid mass | |
| `M2_TAM` | composite | CD163, MRC1, MARCO, MSR1, IL10, TGFB1 | M2-like TAM | DIAL-flag candidate (high in inflamed thyroid via thyroid-specific macrophage biology). |
| `M1_TAM` | composite | NOS2, IL1B, CXCL10, IDO1 (overlaps with IFNγ axis) | M1-like TAM | Overlap risk with IFNγ axis — must be deconvolved. |
| `MDSC_LIKE` | composite | S100A8, S100A9, ARG1, CD33, IDO1, NOS2 | Myeloid-derived suppressor-like | Bulk-level only; scRNA distinguishes better. |
| `NEUTROPHIL` | composite | FCGR3B, CSF3R, FPR1, S100A8, S100A9 | Neutrophil infiltrate | Optional. |

---

## 5. Treg / suppressive T

| Signature ID | Source | Members | Role | Notes |
|---|---|---|---|---|
| `TREG_CORE` | composite | FOXP3, IL2RA, IKZF2, CTLA4, TNFRSF18 | Tregs | Confounded with activation markers in inflamed tissue. |
| `TREG_RATIO_TO_CD8` | derived | (TREG / EFFECTOR_T) | Treg dominance | Ratio metric used in IMPRES-style audit. |

---

## 6. Pan-cancer composite ICI scores (DIAL-audit targets)

These are the signatures most prone to direction-flip in non-trained tissues (here: thyroid). DIAL audit (`paper3_ici_DIAL_audit_plan.md`) specifically tests each of these.

| Signature ID | Source | Method | DIAL risk | Notes |
|---|---|---|---|---|
| `TIDE_LIKE` | Jiang et al. *Nat Med* 2018 | T-cell dysfunction × exclusion composite | High | TIDE proper requires their server; we implement a "TIDE-like" composite using cytotoxic T × M2 + Treg + MDSC inhibition modules. Frame as TIDE-like, not TIDE proper. |
| `IMPRES` | Auslander et al. *Nat Med* 2018 | 15 pairwise gene-pair rules | High | Implemented exactly per published rules; scored as 0–15. |
| `TIS_18` | Ayers 2017 | 18-gene ssGSEA | Medium | T-cell-inflamed phenotype. |
| `IFNG6` | Ayers 2017 | 6-gene ssGSEA | Medium | |
| `CYTOLYTIC` | Rooney 2015 | geometric mean GZMA × PRF1 | Low | Stable across tissues empirically. |
| `EXHAUSTION_INDEX` | composite | EXHAUSTED_T − EFFECTOR_T | High | DIAL audit flag — exhaustion can appear high in tumors that are *non-inflamed*, breaking the assumption. |

---

## 7. Tumor-intrinsic axes

| Signature ID | Source | Members | Role | Notes |
|---|---|---|---|---|
| `BRS` (BRAF-RAS Score) | TCGA-THCA Cell 2014 | 71 genes; pos = BRAF-like, neg = RAS-like | Driver-axis | Required for dark-matter framing; already in Paper 1 ETL. |
| `ERK_SCORE` | Pratilas / Landa | MAPK output genes | Pathway activation | |
| `TDS` (Thyroid Differentiation Score) | TCGA-THCA Cell 2014 | 16 thyroid-differentiation genes | Dedifferentiation axis | Drops sharply in PDTC/ATC. Headline axis for dedifferentiation. |
| `THYROID_LINEAGE` | composite | TG, TPO, TSHR, NIS (SLC5A5), DUOX1, DUOX2, FOXE1, NKX2-1, PAX8 | Tumor-cell lineage identity | Use TPO IHC equivalent for pathology cross-mapping. |

**TDS is the dedifferentiation anchor.** All ecotype results are stratified by TDS tertiles in headline figures.

---

## 8. Single-cell-derived signatures (defined in atlas, projected to bulk)

These are *defined* during scRNA atlas (`paper3_ici_scRNA_reference_atlas_summary.md`) and *registered* here for bulk projection. Members will be populated at Track B Wk 5–6, not now.

| Signature ID | Source state | Method | Status |
|---|---|---|---|
| `scTLS_BCELL_GC` | Germinal-center B-cell state from scRNA atlas | top 50 marker genes | placeholder — populate Wk 5–6 |
| `scCXCL13_TFH` | CXCL13+ Tfh-like CD4 state | top 50 markers | placeholder |
| `scCXCL13_CD8` | CXCL13+ exhausted CD8 (tumor-reactive proxy) | top 50 markers | placeholder |
| `scTAM_INFLAM` | Inflammatory CXCL9/10-high TAM | top 50 markers | placeholder |
| `scTAM_M2_LIKE` | SPP1+/MRC1+ TAM | top 50 markers | placeholder |
| `scCAF_INFLAM` | Inflammatory CAF (CXCL12/IL6) | top 50 markers | placeholder |
| `scDEDIFF_EPI` | Dedifferentiated epithelial state | top 50 markers | placeholder |

---

## 9. Final integrated ICI vulnerability score (composition only — not yet weighted)

Defined here as a registry; weight calibration occurs in Module E during Track B.

```
ICI_VULN_RAW = w1·IFNG18 + w2·MHC1_CORE + w3·MHC2_CORE
             + w4·TLS_CABRITA9 + w5·CXCL13_AXIS
             + w6·NEOANTIGEN_LOAD (Module C) − w7·HLA_LOH (Module C)
             − w8·M2_TAM − w9·MDSC_LIKE − w10·TREG_CORE
             + w11·(1 − TDS)   [dedifferentiation lifts vulnerability]
```

Weights (w1..w11) are not chosen during Track A. Calibration options for Track B:
- (a) equal weighting (sanity baseline)
- (b) DIAL-pruned (drop signatures that fail audit, equal weight on survivors)
- (c) supervised against pan-cancer ICI response in transferred space (only on signatures DIAL passes)
- (d) Cox-supervised against thyroid aggressive-disease outcome (PFS/OS) as a *correlate*, NOT an ICI-response endpoint

Report all four; declare (b) as primary for headline, (c)/(d) as sensitivity. The framing in prose remains "vulnerability / readiness", not "response".

---

## 10. Out-of-scope signatures (explicit exclusions)

- TMB as a single-feature predictor (TMB is a *covariate*, not a signature; documented separately).
- Tumor purity, stromal score (ESTIMATE) — used as covariates, not features.
- HLA-allele specific metrics beyond LOH and supertype — defer to Paper 4 (GD HLA) territory.
- Hashimoto-overlap immune profile — that is Paper 2's analytic claim. Paper 3 may *use* it as a comparator group label only.

---

Track A completed. No marathon violation.

---

## 3. scRNA Reference Atlas Plan

**Working title:** HLA loss, neoantigen architecture, and immune ecotypes define ICI vulnerability in molecularly dark thyroid cancer
**Track:** A — design only. No integration runs. No Seurat / scanpy / scVI / scANVI / Harmony invocation.
**Author:** Seungho Cook
**Status date:** 2026-05-04

**Claim guard:** Atlas defines cell states for *signature derivation* and bulk projection. State frequencies in patients are descriptive, not predictive of ICI response in thyroid until DIAL audit + (where possible) thyroid ICI raw RNA-seq supports such a claim.

---

## 1. Cohorts & expected scale

Per dataset registry (`paper3_ici_dataset_registry.md` §2):

| Accession | Subtype skew | Expected cells (rough order) | Atlas role |
|---|---|---|---|
| GSE184362 | PTC + ATC | 30–80K | Anchor — broadest cell-type coverage |
| GSE193581 | PTC | 20–60K | Immune sub-states |
| GSE232237 | mixed | 20–50K | Integration validation |
| GSE191288 | PTC | 10–40K | Replication |
| GSE148673 | mixed | 10–30K | Replication |
| Han 2024 JCI Insight | dedifferentiation | 10–40K | Dedifferentiation anchor |

**Total estimate:** 100–300K cells after QC. Of these, PDTC + ATC fraction is expected to be small (<5K) — a critical caveat.

**Atlas-side risk register:**
- **R1.** PDTC/ATC scRNA cells <5K total → de novo dedifferentiation trajectory inference is underpowered. Mitigation: project pseudobulk-derived dedifferentiation signature from `paper3_ici_signature_registry` Module 7 onto epithelial cells; do not run Slingshot/Monocle3 on ATC alone.
- **R2.** Sample-level batch effects dominate disease-level effects in thyroid scRNA. Mitigation: scVI with `batch_key=sample_id`, treat patient-level as donor-level; never use disease as integration covariate.
- **R3.** Thyrocyte epithelial cells often contaminated by adjacent normal in surgical specimens. Mitigation: copy-number-based malignant-cell calling (inferCNV / CopyKAT) before downstream epithelial-state analysis.

---

## 2. Reference atlas dependencies (citation only — no runs)

Cell-type annotation will lean on published reference atlases via **transfer learning, not de novo discovery**. Each is referenced read-only — no re-derivation.

| Reference | Use | Tool |
|---|---|---|
| Zheng et al. 2021 (pan-cancer T-cell atlas) | T-cell sub-state labels | scANVI/celltypist |
| Andreatta et al. 2022 ProjecTILs | T-cell sub-state QC | ProjecTILs |
| Mulder et al. 2021 (myeloid pan-tissue) | Myeloid sub-state labels | celltypist |
| Lambrechts / Salcher pan-cancer NSCLC | Stromal CAF sub-states | celltypist |
| HCA thyroid (if available) | Thyrocyte normal reference | scANVI |

If HCA thyroid is unavailable, normal-thyroid cells from each cohort serve as the within-study normal reference.

---

## 3. Integration strategy

**Default:** scVI → scANVI (label-aware fine-tune for cell types) → Harmony as cross-check.

Sequence:
1. Per-cohort QC (scanpy default + `scrublet` doublet removal).
2. Concatenate raw counts; subset to top-3000 highly variable genes per cohort, take union (~5–8K HVGs).
3. scVI training — `n_layers=2`, `n_latent=30`, `batch_key=sample_id`. Single seed first; replicate with 3 seeds for stability.
4. scANVI fine-tune with reference cell-type labels from celltypist seeded annotation.
5. Harmony on the same HVG matrix as alternative integration; UMAP-level concordance reported.
6. Cluster (Leiden, multiple resolutions 0.3 / 0.6 / 1.0); pick resolution by silhouette × biological coherence.

**Acceptance bars:**
- kBET batch metric — improvement vs raw PCA must exceed 0.3.
- Marker-gene coherence — top 20 cluster markers must contain ≥3 known canonical markers for the assigned coarse cell type.
- LISI / iLISI — sample-mixing within cell-type bounds.

---

## 4. Cell-type taxonomy (locked at design)

**Coarse (Level 1) — cluster-then-annotate:**
- Epithelial (thyrocyte normal / malignant)
- T cell (CD4 / CD8 / Treg)
- B cell + plasma
- Myeloid (mono / macro / DC)
- Mast cell
- Endothelial
- Fibroblast / CAF
- Pericyte / smooth muscle

**Fine (Level 2) — within each Level 1, define operational sub-states:**

T-cell:
- Naive / central memory CD8 (TCF7+ LEF1+)
- Effector memory CD8 (GZMK+)
- Effector / tissue-resident CD8 (GZMB+ PRF1+ ITGAE+)
- Exhausted CD8 (PDCD1+ HAVCR2+ TOX+ — possibly CXCL13+ subset = `scCXCL13_CD8` — tumor-reactive proxy)
- CD4 Th1 (CXCR3+ TBX21+)
- CD4 Tfh-like (CXCL13+ ICOS+ BCL6+ — `scCXCL13_TFH`)
- Treg (FOXP3+)
- Proliferating T (MKI67+)

B-cell + plasma:
- Naive B (IGHD+)
- Memory B (CD27+)
- Germinal-center B (BCL6+ AICDA+ — `scTLS_BCELL_GC`) — *referenced for atlas; deeper BCR analysis is Paper 2 territory*
- Plasma cell (MZB1+ JCHAIN+)

Myeloid:
- Classical mono (S100A8+ S100A9+)
- Non-classical mono (FCGR3A+)
- Inflammatory TAM (CXCL9+ CXCL10+ — `scTAM_INFLAM`)
- M2-like / TREM2-high TAM (SPP1+ MRC1+ TREM2+ — `scTAM_M2_LIKE`)
- cDC1 (CLEC9A+ XCR1+)
- cDC2 (CD1C+)
- pDC (LILRA4+)
- LAMP3+ mature DC (mregDC)
- Mast (TPSAB1+)

CAF:
- Inflammatory CAF (CXCL12+ IL6+ — `scCAF_INFLAM`)
- Myo-CAF (ACTA2+)
- Antigen-presenting CAF (HLA-DR+) — flag as candidate Class II source distinct from myeloid

Epithelial (post-malignancy call):
- Normal thyrocyte
- Malignant — well-differentiated (TG+ TPO+)
- Malignant — dedifferentiated (TG-low, TDS-low — `scDEDIFF_EPI`)
- Malignant — IFN-responsive (HLA-DR-high, ISG-high)

---

## 5. Derived single-cell signatures → bulk projection

For each Level-2 sub-state, derive a signature by:
1. Wilcoxon DE of that sub-state vs all other cells of the same Level-1 type.
2. Take top 50 markers with adj. p < 0.01 and logFC > 0.5.
3. Filter against thyroid normal contamination markers (TG, TPO, TSHR top genes excluded).
4. Filter against Y-chromosome / mitochondrial / ribosomal genes.
5. Lock the gene list into `paper3_ici_signature_registry.md` §8 (currently placeholders).

Bulk projection method:
- ssGSEA on bulk cohorts using locked gene lists.
- Validation — predicted state fraction (CIBERSORTx with sc-derived signature matrix) correlated with ssGSEA score across samples; r > 0.5 required.

---

## 6. Headline scRNA outputs

- **Atlas UMAP** — Level 1 + Level 2 cluster overlay, dataset-source overlay, malignancy overlay.
- **State-frequency table** — per patient, fraction of each Level 2 state, separated by PTC vs PDTC vs ATC and by BRAF/RAS/fusion/dark-matter group.
- **Tumor-reactive CD8 (CXCL13+ exhausted) frequency** — per patient, by dedifferentiation tertile (TDS).
- **TLS-like B/Tfh co-occurrence** — per patient, evidence for organized B/Tfh niches; cross-validated against spatial cohorts.
- **TAM polarization landscape** — `scTAM_INFLAM` vs `scTAM_M2_LIKE` per patient.
- **Locked sc-derived signatures** populated into bulk registry.

---

## 7. What this atlas does *not* do

- **Does not** infer ICI response trajectories. State frequencies are descriptive.
- **Does not** re-analyze Paper 2's GSE286332 PTC+HT BCR / TLS / AICDA story. That dataset may appear as one cohort in the atlas, but state-level B/plasma analysis defers to Paper 2.
- **Does not** define DM1 / DM2 molecular clusters — that is Paper 1.
- **Does not** establish neoantigen-T-cell mapping at single-cell level (would require TCR-seq paired with WES; not consistently available across cohorts).

---

## 8. Outputs (filenames to be produced in Track B)

- `project/results/paper3_ici/scrna_atlas/atlas_obs.parquet` — cell-level metadata.
- `project/results/paper3_ici/scrna_atlas/atlas_markers.tsv` — Level 2 markers.
- `project/results/paper3_ici/scrna_atlas/sc_signatures_locked.json` — gene lists for §5.
- `project/results/paper3_ici/scrna_atlas/state_frequency_per_patient.tsv` — for ecotype linkage.
- `project/results/paper3_ici/scrna_atlas/figures/` — UMAPs, heatmaps.

---

Track A completed. No marathon violation.

---

## 4. HLA & Neoantigen Feasibility

**Working title:** HLA loss, neoantigen architecture, and immune ecotypes define ICI vulnerability in molecularly dark thyroid cancer
**Track:** A — design only. No HLA typing runs, no LOHHLA execution, no NetMHCpan invocation. arcasHLA results from `v17_arcasHLA_korean_k2` are referenced read-only and **not** re-run.
**Author:** Seungho Cook
**Status date:** 2026-05-04

**Claim guard:** HLA loss + neoantigen load are reported as *immunogenomic features*, not as standalone "ICI response signals". Combined into the integrated readiness score (`paper3_ici_signature_registry.md` §9) only after DIAL audit.

---

## 1. Per-cohort feasibility matrix

| Cohort | RNA-seq | WES / WGS | HLA Class I | HLA Class II | LOHHLA-feasible | Neoantigen-feasible | Notes |
|---|---|---|---|---|---|---|---|
| TCGA-THCA | yes | yes (BAM via GDC) | OptiType / Polysolver from BAM | arcasHLA from RNA BAM | yes (paired tumor-normal BAM) | yes (MAF + RNA expression filter) | Primary neoantigen + HLA LOH cohort. |
| GSE76039 PDTC/ATC (Landa 2016) | yes | yes (Landa 2016 published WES) | OptiType from RNA BAM (no germline normal) | arcasHLA | partial — paired normal availability per-sample needed | yes if MAF + paired normal | Verify normal availability; without it, somatic calls less reliable. |
| GSE33630 / GSE65144 / GSE54958 / GSE53157 / GSE29265 | microarray | no | no | no | no | no | Microarray — no HLA, no neoantigen. Excluded from Module C. |
| Yoo SK 2019 Korean ATC | yes (RNA-seq) | likely yes (matched WES) | arcasHLA | arcasHLA | conditional on access mode | yes | EGA-controlled candidate; verify at Track B kickoff. |
| PRJKA210106 Korean PTC n=282 | yes (RNA-seq) | unclear | arcasHLA | arcasHLA | likely no (no paired normal) | no (no MAF) | HLA frequency anchor only; cross-paper boundary care. |
| PRJEB11591 Korean PTC n=260 | yes | no | arcasHLA already done | arcasHLA already done | no | no | Read-only reference per `v17_arcasHLA_korean_k2`. **Do not re-run.** |
| Pan-cancer ICI cohorts (DIAL audit set) | yes | partial | per-cohort published HLA | per-cohort | per-cohort | per-cohort | Used as DIAL anchors, not as thyroid HLA evidence. |

**Effective Module C cohort = TCGA-THCA + GSE76039 + (Yoo 2019 if accessible).**

---

## 2. Tool selection (locked at design)

| Step | Primary | Sensitivity / fallback | Why |
|---|---|---|---|
| HLA Class I from WES/WGS BAM | OptiType (paired-end exome) | Polysolver | OptiType is high-accuracy on Class I; Polysolver published for cancer cohorts. |
| HLA Class II from WES | xHLA / HLA*LA | — | WES Class II coverage is uneven; require ≥30× over HLA-DRB1, HLA-DPB1, HLA-DQB1 exons. |
| HLA from RNA-seq | arcasHLA | seq2HLA | arcasHLA already validated by us on PRJEB11591 (`v17_arcasHLA_korean_k2`). Default for cohorts without WES. |
| HLA SNP-imputed 4-digit | cookHLA (Cook et al. *Nat Commun*) | — | Reserved for cohorts with SNP arrays. **Cookhla scope is Paper 4 (GD HLA backlog).** Do not introduce here unless absolutely necessary. |
| HLA LOH | LOHHLA | HLAthena (rare) | Requires paired tumor-normal BAM, allele-specific. |
| Somatic mutation calling | Mutect2 (single-sample mode if no normal); paired Mutect2 if normal exists | Strelka2 | Single-sample mode flagged; PoN required. |
| Variant filtering | dbSNP common-variant filter + gnomAD AF<1e-4 | — | |
| Neoantigen prediction (MHC-I) | NetMHCpan 4.1 | MHCflurry 2.0 | Standard. |
| Neoantigen prediction (MHC-II) | NetMHCIIpan 4.0 | — | Class II is noisier; report binders ≤500 nM with annotation that values are less reliable. |
| Neoantigen RNA-expression filter | TPM ≥1 in tumor RNA-seq for the source gene | — | Required to avoid silent-allele false positives. |
| Self-similarity filter | NetMHCpan against UniProt human reference; remove peptides identical to self | — | |
| Pipeline orchestration | nf-core/epitopeprediction or pVACseq | — | Pick at Wk 1 of Track B based on cohort BAM availability. |

---

## 3. Computational scope (estimate, NOT executed)

Per cohort, single-machine cost estimate (Azure NC-series or local HPC):

| Module | TCGA-THCA (~500) | GSE76039 (~37) | Yoo 2019 (~? aggressive) |
|---|---|---|---|
| HLA typing (RNA-seq + WES) | ~30 CPU-hr | ~3 CPU-hr | ~5 CPU-hr |
| LOHHLA (paired tumor-normal) | ~50 CPU-hr | ~5 CPU-hr | ~5 CPU-hr |
| Mutect2 (paired) | ~200 CPU-hr if BAMs already aligned | ~20 CPU-hr | ~30 CPU-hr |
| NetMHCpan + RNA-expression filter | ~10 CPU-hr | ~1 CPU-hr | ~2 CPU-hr |

**Total Module C estimate:** ~300–400 CPU-hr on TCGA-THCA dominated by Mutect2 if re-running from BAM. If re-using GDC MC3 MAF, Mutect2 cost drops to near-zero and Module C becomes ~50 CPU-hr.

**Decision for Track B:** Reuse GDC MC3 MAF for TCGA somatic calls. Run LOHHLA fresh because MC3 does not provide HLA LOH calls. Run OptiType/Polysolver fresh because per-sample HLA calls may not be packaged.

---

## 4. Headline Module C analyses (defined now, executed later)

1. **HLA typing & frequency** — per cohort, allele frequencies at 4-digit Class I/II. Compare to Korean baseline (PRJEB11591) and TCGA-THCA pan-Asian/Caucasian breakdown.
2. **HLA LOH frequency** — TCGA-THCA paired BAM via LOHHLA. Stratify by BRAF / RAS / fusion / dark-matter / TDS tertile.
3. **Neoantigen load** — per-sample SNV-derived neoantigen count (Class I, Class II), expression-filtered, self-filtered.
4. **Clonal vs subclonal neoantigens** — using PyClone-VI or sciClone CCF estimates if available.
5. **Driver-derived neoantigens** — restrict to BRAF V600E, RAS hotspots, TERT promoter mutations (non-coding so neoantigen unlikely; reported separately), TP53 hotspots — per known immunogenicity priors.
6. **HLA-intact + neoantigen-positive subgroup** — defined as no HLA LOH × neoantigen load ≥ cohort median. Cross-tabulate with bulk immune ecotype (Module A) and dedifferentiation tertile.
7. **PTC vs PDTC vs ATC neoantigen burden** — expected: ATC > PDTC > PTC in TMB; test whether HLA presentation (LOH frequency) co-evolves.

---

## 5. Risks & mitigations

- **R1: TCGA-THCA paired-normal BAMs are slow to retrieve.** Mitigation — request GDC dbGaP-controlled access at Track B Wk 1; in the interim, design and dry-run pipeline against open-access TCGA tutorial BAM.
- **R2: GSE76039 lacks paired normals for some samples.** Mitigation — drop those from LOHHLA + neoantigen; report HLA typing + somatic-MAF-based neoantigen only with caveat.
- **R3: Yoo 2019 may be EGA-controlled.** Mitigation — apply at Wk 1; if not granted by Wk 8, frame as planned validation, drop from primary Module C.
- **R4: Class II prediction is unreliable.** Mitigation — report Class II as descriptive, highlight Class I in headline; Class II analyses framed as exploratory.
- **R5: TMB in PTC is low.** Many PTC samples will have <5 SNV-derived neoantigens. Statistical handling — model as Poisson; report rate, not naive count comparisons.
- **R6: HLA LOH calls noisy in low-purity samples.** Mitigation — purity ≥ 0.4 inclusion gate; sensitivity at 0.3 reported.
- **R7: cookHLA is the user's first-author tool but scope is Paper 4.** Mitigation — do not introduce cookHLA here. Class I from OptiType/Polysolver/arcasHLA only.

---

## 6. Cross-paper boundary (binding)

- **Paper 1 (DM1 molecular)** — Paper 1 does not include HLA typing or neoantigen prediction. Module C is fully owned by Paper 3.
- **Paper 2 (H&E-DM1 / HT-overlap PTC)** — Paper 2's Pillar 1 forest is HLA *susceptibility* allele frequency, Korean PTC vs Korean baseline. Different question (susceptibility ≠ tumor HLA loss). Paper 3 Module C analyzes *somatic* HLA loss in tumors and *neoantigen* presentation; explicitly distinct from Paper 2 Pillar 1.
- **Paper 4 (GD HLA backlog)** — cookHLA SNP imputation belongs to Paper 4. Paper 3 does not run cookHLA.

---

## 7. Outputs (filenames to be produced in Track B Module C)

- `project/results/paper3_ici/hla/hla_class1_calls.tsv`
- `project/results/paper3_ici/hla/hla_class2_calls.tsv`
- `project/results/paper3_ici/hla/lohhla_calls.tsv`
- `project/results/paper3_ici/neo/neoantigens_class1.tsv`
- `project/results/paper3_ici/neo/neoantigens_class2.tsv`
- `project/results/paper3_ici/integrated/hla_intact_neoag_pos_subgroup.tsv`

---

## 8. Module C → ICI vulnerability score linkage

- `HLA_INTACTNESS = 1 − HLA_LOH_indicator × allele_lost_fraction`
- `NEO_LOAD_LOG = log10(neoantigen_count_class1 + 1)` (expression-filtered)
- `NEO_PRESENTABLE = NEO_LOAD_LOG × HLA_INTACTNESS`

These three feed Module E's integrated score with weight calibration deferred to Track B.

---

Track A completed. No marathon violation.

---

## 5. DIAL Audit Plan

**Working title:** HLA loss, neoantigen architecture, and immune ecotypes define ICI vulnerability in molecularly dark thyroid cancer
**Track:** A — design only. No DIAL execution. No bootstrap, no permutation, no Cox fitting.
**Author:** Seungho Cook
**Status date:** 2026-05-04

**Claim guard:** DIAL audit is the *gating step* that decides which pan-cancer ICI signatures are allowed to enter Paper 3's integrated readiness score in thyroid. Without DIAL, the integrated score has no validity argument.

---

## 1. Context — DIAL framework lineage

- DIAL = Direction-Invariance ALgorithm. Originated in `v52_lodo_finding`: under proper LODO ComBat, all 5 THCA classifiers returned DIAL=0.000 / true biology, falsifying the v5.1 DIAL=0.494 flip artifact.
- Codified into `v18_agentic_research framework` (5 patterns + 6 components).
- Applied in Paper 3 to **pan-cancer ICI signatures** rather than to THCA classifiers — same algebra, different inputs.

**Definition for Paper 3 use.** A signature S is *direction-invariant* across context set C if its sign of association with the chosen endpoint (ICI response or T-cell-inflamed phenotype) is preserved across all leave-one-out subsets of C, with effect size > minimum. Failure modes:
- **Flip** — sign reverses in ≥1 leave-one-out fold.
- **Collapse** — effect size shrinks to noise (|t| < 1) in ≥1 fold.
- **Tissue-flip** — sign holds within trained tissue (e.g., melanoma) but reverses in held-out tissue (here: thyroid).

Only signatures that pass DIAL on the pan-cancer reference set are eligible to enter Paper 3's integrated score in thyroid.

---

## 2. Reference set & endpoints

**Pan-cancer ICI reference (per `paper3_ici_dataset_registry.md` §4):**
- IMvigor210 (UC, atezolizumab) — Tier 1
- Hugo GSE78220 (melanoma, anti-PD-1) — Tier 1
- Riaz GSE91061 (melanoma, nivolumab; pre-treatment only) — Tier 1
- Gide melanoma (PRJEB23709) — Tier 2
- Liu melanoma (dbGaP) — Tier 2 if access granted
- Kim gastric (PRJEB25780) — Tier 3
- Cho NSCLC — Tier 3

**Endpoints:**
- Primary — RECIST objective response (CR/PR vs SD/PD).
- Secondary — PFS (Cox HR).
- Tertiary — T-cell-inflamed phenotype (TIS-positive vs negative) as a *biological* anchor independent of clinical response.

**Held-out tissue = thyroid:**
- TCGA-THCA aggressive subset (PDTC/ATC + Stage IV PTC + dedifferentiated) as proxy for "ICI-eligible thyroid".
- GSE76039 (PDTC + ATC) as second held-out cohort.
- No clinical ICI response endpoint in thyroid — endpoint in held-out is *T-cell-inflamed phenotype* (TIS) as proxy. Direction-invariance from response → TIS is checked first within reference cohorts; only signatures invariant under both are claimed.

---

## 3. Audit pipeline (defined; NOT executed)

```
For each signature S in pan-cancer registry:
    For each cohort c in reference set:
        Compute S on c (ssGSEA / IMPRES rules / Cytolytic GM)
        Fit logistic regression: response ~ S + clinical_covariates
        Record beta_S, se, p

    Direction-invariance check 1 (within-reference, response endpoint):
        sign_consistent = all(sign(beta_S, c) == majority_sign across c)
        effect_size_min  = min(|beta_S, c| > threshold)
        DIAL_response    = sign_consistent AND effect_size_min

    Direction-invariance check 2 (within-reference, TIS endpoint):
        Repeat with TIS_18 score as endpoint
        DIAL_TIS = sign_consistent AND effect_size_min

    Tissue-transfer check (held-out thyroid, TIS endpoint):
        Compute S on TCGA-THCA aggressive + GSE76039
        Fit linear regression: TIS ~ S
        DIAL_thyroid_TIS = (sign(beta_S, thyroid) == majority sign in reference)

    DIAL verdict:
        PASS if DIAL_response AND DIAL_TIS AND DIAL_thyroid_TIS
        FLIP if any check has reversed sign
        COLLAPSE if effect size below threshold
        AMBIGUOUS otherwise

    Record verdict in paper3_ici_DIAL_results.tsv
```

Bootstrap — 1000 resamples per cohort to attach CI to beta_S. Permutation — 1000 label-shuffled fits to attach empirical p.

---

## 4. Signatures audited

From `paper3_ici_signature_registry.md` (binding list):

| Signature ID | DIAL prior risk | Why this risk |
|---|---|---|
| `IFNG_AYERS6` | Low | Robust across tissues in published literature. |
| `TIS_18` | Low | T-cell-inflamed phenotype well-established. |
| `CYTOLYTIC` | Low | 2-gene; very stable. |
| `MHC1_CORE` | Low | Mechanistic — direction predictable. |
| `MHC2_CORE` | **Medium** | Tumor-intrinsic Class II expression in thyroid (HT-overlap signal in Paper 2) may flip thyroid direction. |
| `TLS_CABRITA9` | Low | Established in melanoma; expected to hold in thyroid. |
| `TIDE_LIKE` | **High** | Composite includes M2 + Treg + MDSC; some components flip in thyroid. |
| `IMPRES` | **High** | 15-pair rule trained against melanoma; pairwise rules may not transfer. |
| `EXHAUSTION_INDEX` | **High** | Exhaustion can co-occur with non-inflamed in thyroid; semantic flip likely. |
| `M2_TAM` | **Medium** | Thyroid-specific macrophage biology may invert "M2 = bad" rule. |
| `MDSC_LIKE` | **Medium** | Confounded with neutrophil signal in thyroid. |
| `TREG_CORE` | Medium | Treg in thyroid can mark inflamed-but-controlled state. |
| `EFFECTOR_T` | Low | Stable. |
| `CXCL13_AXIS` | Low | Tumor-reactive proxy stable. |

---

## 5. Statistical thresholds

- Sign-consistency required across ≥6 of 7 reference cohorts (allow 1 outlier).
- Effect-size threshold |beta_S| > 0.2 (logistic, standardized predictor).
- Empirical permutation p < 0.05.
- Tissue-transfer check |beta_S, thyroid TIS regression| ≥ 0.1 with same sign as reference majority.
- Multiple-testing correction — Benjamini-Hochberg across the signature set; report q-values.

---

## 6. Pre-specified failure handling

- **All composite ICI scores fail (TIDE_LIKE, IMPRES, EXHAUSTION_INDEX) →** report this as the headline DIAL finding ("pan-cancer composites do not transfer to thyroid"), use only the surviving signatures (likely IFNG_AYERS6, TIS_18, CYTOLYTIC, MHC1_CORE, TLS_CABRITA9, EFFECTOR_T, CXCL13_AXIS) in Module E. This is a publication-strengthening result, not a setback.
- **MHC2_CORE flips in thyroid →** investigate whether thyrocyte-intrinsic Class II expression (Hashimoto-like) is the driver, citing Paper 2. Report as biology-driven flip with mechanistic interpretation.
- **All signatures pass with same direction →** weakest publication outcome (no thyroid-specific story). Pivot prose toward "thyroid behaves canonically" framing; lean harder on Module C HLA LOH + neoantigen architecture for the headline.

---

## 7. Thyroid-adjusted readiness score

Once DIAL audit closes, define:

```
ICI_VULN_DIAL = sum over signatures S where DIAL(S) == PASS of
                (sign_majority_S * w_S * z_S)
              − sum over Module C terms (HLA loss, neoantigen presented)
```

Weights for surviving signatures = equal (default) and DIAL-effect-size-weighted (sensitivity).

This score is reported at PTC vs PDTC vs ATC, BRAF / RAS / fusion / dark-matter, and TDS-tertile strata.

---

## 8. What this audit does *not* do

- **Does not** establish that thyroid patients with high `ICI_VULN_DIAL` will respond to ICI. That requires thyroid ICI-treated raw RNA-seq (which we do not have). Frame as readiness *prioritization*, not response prediction.
- **Does not** retrofit thresholds after seeing thyroid results — DIAL verdict locks before thyroid scoring.
- **Does not** replace TIDE proper — we audit TIDE-like composite and report TIDE proper (via TIDE server) as supplementary at Track B.

---

## 9. Outputs (filenames to be produced in Track B Module D)

- `project/results/paper3_ici/dial/per_cohort_betas.tsv`
- `project/results/paper3_ici/dial/dial_verdict.tsv` (per-signature PASS/FLIP/COLLAPSE/AMBIGUOUS)
- `project/results/paper3_ici/dial/sign_consistency_heatmap.pdf`
- `project/results/paper3_ici/dial/thyroid_transfer_betas.tsv`
- `project/results/paper3_ici/dial/audit_report.md`

---

## 10. Sequence within 12-week plan

DIAL audit runs Wk 9, after bulk ecotype (Wk 3–4), scRNA atlas (Wk 5–6), and HLA/neoantigen (Wk 7–8). Pan-cancer ICI cohorts must be acquired by Wk 7. If a cohort lags, downgrade DIAL to remaining cohorts and document.

---

Track A completed. No marathon violation.

---

## 6. Figure Plan

**Working title:** HLA loss, neoantigen architecture, and immune ecotypes define ICI vulnerability in molecularly dark thyroid cancer
**Track:** A — design only. All figures are layout schematics. No real-data rendering. Captions are placeholder; final captions written during Wk 11.
**Author:** Seungho Cook
**Status date:** 2026-05-04

**Claim guard:** Each figure caption ends with a thyroid-ICI-vulnerability framing, never an "ICI response prediction" claim. Discrepancy between visualization and claim is a stop-the-press reviewer trigger; designer enforces this at draft.

---

## 1. Main figures (target = 6, ≤7)

### Figure 1 — Cohort landscape & dark-matter framing

Panels:
- **1A** Schematic — PTC → PDTC → ATC dedifferentiation continuum; BRAF / RAS / fusion / dark-matter color key; immune-rich (HT-overlap, dark matter immune) overlay.
- **1B** Sample count matrix across cohorts × subtype (TCGA-THCA, GSE76039, others) with modality availability (RNA / WES / scRNA / spatial).
- **1C** TDS distribution by subtype + driver.
- **1D** Dark-matter (BRAF/RAS/fusion-negative) prevalence by subtype.
- **1E** Cross-paper diagram showing Paper 1 / Paper 2 / Paper 3 / Paper 4 scope boundaries — defensive against reviewer scope-confusion.

Headline: dark matter exists across the differentiation spectrum and is the analytical target.

### Figure 2 — Bulk immune ecotype discovery (Module A)

- **2A** Heatmap of standardized signature scores across pooled bulk samples (rows = signatures from registry §1–§5, cols = samples), NMF-derived ecotype labels.
- **2B** NMF rank selection (cophenetic correlation, dispersion).
- **2C** Ecotype × subtype (PTC / PDTC / ATC) Sankey or stacked bar.
- **2D** Ecotype × driver group (BRAF / RAS / fusion / dark matter) stacked bar with chi-square.
- **2E** Ecotype × TDS tertile.
- **2F** Survival KM (PFS or DSS) by ecotype within aggressive subset (TCGA-THCA stage III/IV + PDTC + ATC).

Headline: thyroid immune ecotypes are non-uniform across dedifferentiation and dark-matter axes.

### Figure 3 — scRNA atlas & sub-state map (Module B)

- **3A** Integrated UMAP, Level 1 cell types.
- **3B** Level 2 sub-state UMAP for T cell / myeloid / B cell.
- **3C** Marker dot-plot for Level 2 states (top 5 markers per state).
- **3D** State frequency per patient × bulk-ecotype (link to Fig 2) — uses pseudobulk → ecotype mapping.
- **3E** CXCL13+ exhausted CD8 (`scCXCL13_CD8`) frequency vs TDS tertile.
- **3F** TLS-like B/Tfh co-occurrence per patient.

Headline: dedifferentiated and dark-matter subgroups carry distinct sub-state architectures, with CXCL13+ exhausted CD8 enriched at the dedifferentiated end.

### Figure 4 — HLA loss & neoantigen architecture (Module C)

- **4A** HLA Class I / II calls — allele-frequency comparison (TCGA-THCA Caucasian vs Asian; vs PRJEB11591 Korean baseline).
- **4B** HLA LOH frequency — overall, then stratified by BRAF / RAS / fusion / dark matter and PTC vs PDTC vs ATC.
- **4C** Neoantigen burden (Class I) — Poisson-modeled rate by subtype × driver.
- **4D** Quadrant plot — HLA-intactness (x) vs neoantigen presented (y), colored by ecotype (Fig 2).
- **4E** Cohort labeled "HLA-intact + neo-positive + dedifferentiated" — sample list / count by cohort.

Headline: an HLA-intact + neoantigen-positive + dedifferentiated subgroup exists and overlaps the inflamed ecotype.

### Figure 5 — DIAL audit (Module D)

- **5A** Pan-cancer reference matrix — signatures (rows) × cohorts (cols) signed beta heatmap.
- **5B** Sign-consistency bar — fraction of cohorts agreeing with majority sign per signature.
- **5C** Effect-size (|beta|) per signature per cohort with bootstrap CI.
- **5D** Tissue-transfer plot — reference majority beta vs thyroid TIS-regression beta; quadrant coloring (PASS / FLIP / COLLAPSE / AMBIGUOUS).
- **5E** DIAL verdict table — per-signature PASS/FLIP/COLLAPSE/AMBIGUOUS.

Headline: pan-cancer composites partially fail to transfer to thyroid; thyroid-adjusted readiness score must restrict to DIAL-passing signatures.

### Figure 6 — Integrated ICI readiness score & subgroup map (Module E)

- **6A** Score distribution across PTC / PDTC / ATC and dark-matter.
- **6B** Score components — stacked contributions per sample.
- **6C** Top decile readiness samples — characterize: HLA-intact + neo+ + CXCL13/TLS+ + dedifferentiated overlap.
- **6D** Validation panel — pan-cancer transfer of the thyroid-adjusted score (predict response in IMvigor210 / Hugo / Riaz, sanity check).
- **6E** Indirect spatial validation — TLS niche + CXCL13 in spatial cohort overlay (read-only from Paper 1/2 spatial pipeline).
- **6F** Schematic — proposed ICI-readiness candidate group definition. Caption explicitly states *prioritization hypothesis*, not *response predictor*.

Headline: a thyroid-adjusted readiness score identifies a candidate ICI-vulnerability subgroup at the intersection of dedifferentiation, HLA-intactness, neoantigen presentation, and TLS/CXCL13.

---

## 2. Supplementary figures

| ID | Content |
|---|---|
| S1 | Per-cohort QC — read counts, mapping rate, gene detection. |
| S2 | scRNA per-cohort UMAP before integration; integration metrics (kBET, LISI). |
| S3 | NMF rank stability; alternative consensus clustering. |
| S4 | Sensitivity to ssGSEA vs Singscore; sensitivity to UCell vs AUCell. |
| S5 | Microarray-cohort signature replication (subset where genes available). |
| S6 | LOHHLA QC — purity, BAF, allele-fraction plots per representative sample. |
| S7 | Neoantigen pipeline waterfall — per-step sample loss; Class I vs Class II yield. |
| S8 | DIAL bootstrap distributions per signature. |
| S9 | Cox sensitivity — score vs PFS/DSS within aggressive subset, full multivariate. |
| S10 | TIDE-server output (proper TIDE) vs our TIDE-like composite — concordance. |
| S11 | Spatial overlay (read-only from Paper 1/2) showing TLS / CXCL13 localization. |
| S12 | Alternative score weightings (a/b/c/d in registry §9) — concordance heatmap. |
| S13 | Pan-Asian HLA frequency anchor (PRJEB11591) — read-only — context only. |
| S14 | Limitations table — claims, evidence level, generalizability scope. |

Target: ≤14 supp figures. Trim at Wk 11.

---

## 3. Tables

| ID | Content |
|---|---|
| T1 | Cohort summary (subset of dataset registry). |
| T2 | Signature registry (subset of `paper3_ici_signature_registry.md`). |
| T3 | DIAL verdict per signature. |
| T4 | HLA-intact + neoantigen-positive + dedifferentiated subgroup roster. |
| T5 | Integrated readiness score per sample (anonymized cohort IDs). |
| ST1 | Full per-cohort betas, CI, p (DIAL). |
| ST2 | Full HLA Class I/II calls. |
| ST3 | Full neoantigen list (filtered). |
| ST4 | Pan-cancer transfer validation per dataset. |

---

## 4. Visual conventions (locked at design)

- Color palette — colorblind-safe (Okabe-Ito 8-color extended). Subtype: PTC=blue, PDTC=orange, ATC=red. Driver: BRAF=teal, RAS=purple, fusion=green, dark-matter=grey.
- Ecotype labels follow E1–E5 numbering (NMF rank-dependent at Track B).
- Headline figures rendered at 300 dpi, vector .pdf; supplements may be raster .png.
- Score visualizations use viridis / cividis only; never rainbow.
- Statistical annotations — exact p-values to 2 sig figs; q-values for multiple-tested panels; effect sizes (Cohen's d / log2 FC / beta) reported alongside p.
- Sample sizes annotated in every panel caption.

---

## 5. Paper-boundary defensive framing in figure captions

Each caption that references a comparator group includes a one-line scope guard:
- "The Hashimoto-overlap PTC phenotype shown for reference is analyzed in detail in Cook et al. (Paper 2)."
- "The DM1 molecular subtype labeling is from Cook et al. (Paper 1)."
- "Korean GD HLA susceptibility (Cook et al., Paper 4) shares background frequencies but is mechanism-distinct."

This prevents reviewer scope-confusion and aligns with cross-paper boundary discipline in `v19_paper3_ici_track_a.md` and `v18_paper2_HT_isolated.md`.

---

## 6. What is *not* in the figure plan

- No survival curve for ICI response in thyroid — no thyroid ICI-treated cohort.
- No nomogram of "ICI response probability" — would over-claim.
- No clinical decision tree or per-patient recommendation — premature.
- No DM1/DM2 cluster figures — Paper 1 territory.
- No PTC+HT BCR/TLS detailed figure — Paper 2 territory.
- No GD HLA forest figure — Paper 4 territory.

---

Track A completed. No marathon violation.

---

## 7. Go / No-Go Verdict

**Working title:** HLA loss, neoantigen architecture, and immune ecotypes define ICI vulnerability in molecularly dark thyroid cancer
**Track:** A — design only. Verdict below is **conditional GO** for Track B execution at the post-marathon entry point.
**Author:** Seungho Cook
**Status date:** 2026-05-04

**Claim guard:** Below verdict assumes Paper 3 is framed as ICI vulnerability / readiness / immunogenomic prioritization. If the user later wants to claim "ICI response prediction" without thyroid ICI-treated raw RNA-seq, verdict drops to NO-GO and must be re-decided.

---

## 1. Module-level verdicts

| Module | Verdict | Rationale | Critical dependency |
|---|---|---|---|
| **A — Bulk immune ecotype** | **GO** | TCGA-THCA + GSE76039 + microarray cohorts give ≥600 samples after pooling; ssGSEA + NMF mature; thyroid ecotypes literature-supported. | None blocking. |
| **B — scRNA atlas** | **CONDITIONAL GO** | 6 candidate cohorts, total ~100–300K cells. PDTC/ATC scRNA <5K — atlas viable for state annotation, not for de novo dedifferentiation trajectory. Mitigation framed in atlas plan §1. | At least 3 of 6 cohorts must be accessible at Wk 5. |
| **C — HLA & neoantigen** | **GO (degraded scope)** | TCGA-THCA paired BAM + GDC MC3 MAF available; LOHHLA + NetMHCpan well-established; GSE76039 paired-normal access uncertain per sample; Yoo 2019 may be EGA-controlled. | TCGA-THCA dbGaP controlled-access for paired BAM by Wk 1 of Track B. |
| **D — DIAL audit** | **GO** | Pan-cancer ICI cohorts mostly accessible (IMvigor210 R-package; Hugo / Riaz / Gide GEO; Liu / Cho controlled). Statistical machinery already in `v18_agentic_research framework`. | At least 3 of 5 Tier 1/2 melanoma + UC cohorts accessible by Wk 7. |
| **E — Integrated readiness score** | **GO conditional on D** | Score is well-defined; weighting calibration deferred. | DIAL audit must produce ≥4 PASS signatures for the score to be meaningful. |

**Overall verdict — CONDITIONAL GO.** Paper 3 is feasible. The conditions are dataset-access events that map to the 12-week plan's Wk 1 / Wk 5 / Wk 7 gates, not analytical uncertainties.

---

## 2. Hard gates (must close before any Track B compute spend)

| Gate | Condition | Action if not met |
|---|---|---|
| **G1** | Paper 1 bioRxiv submitted | Track B blocked. Marathon mode preserved per `v17_marathon_mode_post_pillar1`. |
| **G2** | Paper 2 Task A/B/C closure | Track B blocked. |
| **G3** | TCGA dbGaP access (or confirmation that GDC open access is sufficient for our scope) | Module C drops to "RNA-seq HLA only via arcasHLA + MAF-based neoantigen using MC3 SNV calls without LOHHLA". Headline figure 4B (HLA LOH) downgraded to "future-work" if G3 fails. |
| **G4** | At least 3 of 6 scRNA cohorts accessible | Module B drops to "single-cohort sub-state map" with reduced replication; figure 3 panels narrow. |
| **G5** | At least 3 of 5 Tier 1/2 pan-cancer ICI cohorts accessible | Module D narrows; if <3, DIAL is reported as plan-only and Module E uses unaudited signatures with explicit caveat. |
| **G6** | User explicit "Track B 시작" — not "고고" / "다 해줘" | Per `v17_sprint_vs_marathon_violation.md`. |

---

## 3. Risks that DO NOT trigger NO-GO

- Yoo 2019 inaccessible — replace with TCGA aggressive subset; degrade Asian-validation claim, do not block.
- Cho NSCLC / Kim gastric inaccessible — drop, retain Tier 1 melanoma + UC anchors, do not block.
- HCA thyroid reference unavailable — use within-study normal cells, do not block.
- Class II HLA calls noisy — frame as exploratory, do not block.

---

## 4. Risks that WOULD trigger NO-GO (kill switches)

| Trigger | Reason |
|---|---|
| K1 — User attempts to claim "ICI response predictor in thyroid" without thyroid ICI-treated raw RNA-seq | Mis-claim invalidates paper. Hard stop. |
| K2 — Paper 1 ship slips past late August 2026 such that marathon mode is still active | Track B must remain blocked; no execution until ship. |
| K3 — DIAL audit shows zero PASS signatures | Module E has no inputs. Pivot Paper 3 to a feasibility-only / negative-result paper, but this is a different paper. |
| K4 — Paper-boundary contamination — Paper 3 prose pulls Paper 2's PTC+HT TLS / AICDA / BCR claims as Paper 3's own findings | Reviewer scope-rejection certain. Hard stop until rewritten. |
| K5 — TCGA-THCA driver-call disagreement vs Paper 1 ETL produces inconsistent dark-matter rosters | Recompute against Paper 1 source-of-truth; if cannot reconcile, halt Module C/E until resolved. |
| K6 — User pivots priority again (e.g., back to GD or to a new Paper 5) | Re-decide stack; do not silently continue. |

---

## 5. Decision dependencies on Paper 1 / Paper 2

- Paper 1 ship (target 2026-06-13) frees marathon constraint per `v17_marathon_mode_post_pillar1`.
- Paper 2 Task A/B/C closure clarifies which PTC+HT analyses are Paper 2 and which are off-limits to Paper 3.
- Yu professor agreement (per cross-paper boundary discipline in `v19_paper3_ici_track_a.md`) before Paper 3 enters compute.

---

## 6. Headline summary

Paper 3 (ICI vulnerability in molecularly dark thyroid cancer) is **feasible in design** with **CONDITIONAL GO** for Track B. The five modules are mature; the dependencies are dataset-access events and the Paper 1/2 gating events, not analytical uncertainties.

The only framing in which the verdict reverses to NO-GO is over-claiming "ICI response predictor in thyroid" without thyroid ICI-treated raw RNA-seq. Track A's claim guard makes this explicit and binds the Track B prose.

**Recommendation:** lock Track A artifacts as the paper's design freeze; do not begin Track B compute until G1–G6 close.

---

Track A completed. No marathon violation.

---

## 8. 12-Week Execution Plan

**Working title:** HLA loss, neoantigen architecture, and immune ecotypes define ICI vulnerability in molecularly dark thyroid cancer
**Track:** A (this document) — design only. Plan below describes Track B execution; Track B starts only after entry gates G1–G6 close.
**Author:** Seungho Cook
**Status date:** 2026-05-04

**Claim guard:** All weeks below assume "ICI vulnerability / readiness / immunogenomic prioritization" framing. If thyroid ICI-treated raw RNA-seq is acquired during the plan, framing may upgrade — only after re-decision with user, never silently.

---

## 0. Entry preconditions (G1–G6)

Plan does not start until:
- G1 Paper 1 bioRxiv submitted (target 2026-06-13)
- G2 Paper 2 Task A/B/C closed
- G3 TCGA dbGaP access decided
- G4 ≥3 of 6 scRNA cohorts accessible
- G5 ≥3 of 5 Tier 1/2 pan-cancer ICI cohorts accessible
- G6 User explicit "Track B 시작"

If G1/G2 not closed, the plan re-baselines its Wk 1 day to the post-ship date.

---

## 1. Week-by-week schedule

### Wk 1 — Dataset acquisition & ETL scaffolding

- Verify all `to_verify` accessions in `paper3_ici_dataset_registry.md` §1–§4. Produce `paper3_ici_dataset_registry_VERIFIED.tsv`.
- Submit dbGaP / EGA applications in parallel (Liu melanoma, Yoo 2019 if EGA, TCGA paired-BAM controlled if needed).
- Stand up `project/results/paper3_ici/` directory tree per registry / atlas / hla / neo / dial / integrated outputs.
- Reuse Paper 1 ETL outputs as read-only inputs (driver / fusion calls, BRS, TDS).
- Lock signature gene lists to actually-measurable subset per cohort. Produce `paper3_ici_signature_registry_VERIFIED.tsv`.

Deliverables: verified registries, directory scaffolding, dbGaP/EGA applications submitted.

### Wk 2 — Bulk preprocessing & cohort harmonization

- Re-pull / re-normalize bulk cohorts to a common log2(TPM+1) baseline (RNA-seq) and RMA log-intensity (microarray). Reuse Paper 1 pipelines where applicable; do not re-run if no scope difference.
- ComBat / ComBat-seq batch correction per cohort family with subtype as protected covariate; LODO sensitivity per `v52_lodo_finding`.
- QC pass: read count, mapping rate, dropout, gene detection.
- Compute all bulk signatures from `paper3_ici_signature_registry.md` §1–§7 on harmonized matrix; lock `paper3_ici_bulk_signature_scores.parquet`.

Deliverables: harmonized bulk matrix, QC tables, all-signature score matrix.

### Wk 3 — Bulk immune ecotype discovery (Module A part 1)

- NMF rank selection (K = 3..8) on signature-z-score matrix; cophenetic correlation, dispersion. Pick K.
- Consensus clustering as alternative; compare ARI with NMF.
- Stability: 100 bootstrap subsamples; per-sample assignment confidence.
- Headline figure 2A/B drafts.

Deliverables: ecotype labels, stability metrics, fig 2A/B drafts.

### Wk 4 — Bulk ecotype × clinical/molecular axes (Module A part 2)

- Cross-tabulate ecotype × subtype (PTC/PDTC/ATC), driver (BRAF/RAS/fusion/dark-matter), TDS tertile.
- Survival in aggressive subset (TCGA stage III/IV + GSE76039 PDTC/ATC) with Cox.
- Sensitivity to scoring method (ssGSEA vs Singscore) and to NMF seed.
- Lock `paper3_ici_bulk_immune_ecotype_summary.tsv`.
- Fig 2C/D/E/F drafts.

Deliverables: ecotype summary table, fig 2 first complete draft.

### Wk 5 — scRNA atlas integration (Module B part 1)

- Per-cohort QC, doublet removal, normalization.
- HVG selection, scVI training, scANVI fine-tune, Harmony cross-check.
- kBET / LISI / iLISI metrics; integration acceptance per `paper3_ici_scRNA_reference_atlas_summary.md` §3.
- Coarse Level 1 cluster + annotation (celltypist seeded).

Deliverables: integrated AnnData, integration QC, Level 1 annotations.

### Wk 6 — scRNA sub-state mapping & signature lock (Module B part 2)

- Level 2 sub-state annotation per cell type.
- DE per sub-state; lock `sc_signatures_locked.json` per registry §8.
- inferCNV / CopyKAT for malignant epithelial calling.
- Per-patient state frequency table (`state_frequency_per_patient.tsv`).
- Bulk projection of sc-derived signatures via ssGSEA + CIBERSORTx cross-check.
- Fig 3 drafts.

Deliverables: locked sc signatures, per-patient state frequencies, fig 3 first complete draft.

### Wk 7 — HLA typing (Module C part 1)

- TCGA-THCA WES BAM via GDC controlled access — OptiType + Polysolver Class I; xHLA / HLA*LA Class II.
- arcasHLA on RNA-seq for cohorts without WES.
- Lock `hla_class1_calls.tsv`, `hla_class2_calls.tsv`.
- Allele frequency comparison vs Korean baseline (PRJEB11591) and TCGA pan-Asian/Caucasian breakdown.

Deliverables: locked HLA call tables, allele-frequency tables.

### Wk 8 — HLA LOH + Neoantigen prediction (Module C part 2)

- LOHHLA on TCGA-THCA paired tumor-normal BAM with purity ≥0.4 inclusion.
- Mutect2 paired re-run only if MC3 MAF inadequate; otherwise reuse MC3.
- NetMHCpan 4.1 Class I + NetMHCIIpan 4.0 Class II prediction with RNA-expression filter (TPM ≥1 in tumor).
- Self-similarity filter against UniProt human reference.
- Lock `lohhla_calls.tsv`, `neoantigens_class1.tsv`, `neoantigens_class2.tsv`.
- HLA-intactness and neoantigen-presented composite per `paper3_ici_HLA_neoantigen_feasibility.md` §8.
- Fig 4 drafts.

Deliverables: HLA LOH + neoantigen tables, HLA-intact + neo+ subgroup roster, fig 4 first complete draft.

### Wk 9 — DIAL audit (Module D)

- Pan-cancer ICI cohort signature scoring on the same registry.
- Logistic regression (response endpoint) and TIS-regression (T-cell-inflamed endpoint) per cohort.
- Bootstrap × 1000, permutation × 1000 for sign + effect size CI.
- Tissue-transfer check on TCGA aggressive subset + GSE76039 (TIS endpoint).
- Lock `dial_verdict.tsv`.
- Fig 5 drafts.

Deliverables: DIAL verdict per signature, fig 5 first complete draft.

### Wk 10 — Integrated readiness score (Module E)

- Define `ICI_VULN_DIAL` using DIAL-passing signatures only (default = equal weight on survivors).
- Sensitivity weightings (a/b/c/d) per signature registry §9.
- Stratify score by subtype, driver, TDS, ecotype, HLA-intact + neo+ subgroup.
- Pan-cancer transfer validation — predict response in IMvigor210 / Hugo / Riaz with the *thyroid-adjusted* score (sanity).
- Spatial overlay (read-only from Paper 1/2 spatial pipeline).
- Lock `integrated_readiness_score.tsv`.
- Fig 6 drafts.

Deliverables: readiness score per sample, fig 6 first complete draft.

### Wk 11 — Manuscript drafting (figures locked)

- Figure captions finalized; supplementary figures S1–S14 drafted.
- Methods written from M-tier scaffolding (mirroring Paper 1's `methods_M1_M11_scaffold.md` style).
- Results sections following figure order: F1 cohort → F2 ecotype → F3 sc atlas → F4 HLA/neo → F5 DIAL → F6 integrated score.
- Discussion drafted by user (voice-protected; per `v17_sprint_vs_marathon_violation.md` no Claude generation in voice-protected sections; Claude provides scaffolding).
- Limitations section explicit on thyroid ICI raw RNA-seq gap.
- References built on shared `2026_05_03_references.bib` extended with ICI literature.

Deliverables: manuscript v1 with locked figures + methods + results; user-drafted discussion outline.

### Wk 12 — Internal review & revision

- User-led discussion completion.
- Internal QA: claim guard sweep — flag any "response predictor" language; rewrite to "vulnerability / readiness".
- Cross-paper boundary sweep — flag any prose that overlaps Paper 1 (DM1 cluster) or Paper 2 (PTC+HT BCR/TLS) or Paper 4 (GD HLA forest).
- Cohort and accession verification refresh — every accession in dataset registry confirmed final.
- Reviewer Q anticipation list (mirroring Paper 1's `2026_05_03_reviewer_QA_consolidated.md`).
- Cover letter draft (mirroring Paper 1's `2026_05_03_cover_letter_assembly.md`).

Deliverables: manuscript v2 ship-ready, reviewer Q list, cover letter draft.

---

## 2. Critical-path & parallelism

Critical path: Wk 1 (data) → Wk 2 (preprocess) → Wk 3–4 (Module A) → Wk 7–8 (Module C) → Wk 9 (Module D) → Wk 10 (Module E) → Wk 11 (manuscript).

Parallel branches:
- Wk 5–6 (Module B scRNA) runs alongside Wk 3–4 / Wk 7–8 because it does not block Module A or C.
- DIAL pan-cancer ICI ETL can begin Wk 5 in parallel with Module B.
- Manuscript scaffolding (Methods / supp tables) can begin Wk 8 in parallel with Module D.

---

## 3. Compute budget (rough estimate)

| Module | Compute | Notes |
|---|---|---|
| Bulk preprocess + signatures | low (~50 CPU-hr) | mostly normalization. |
| NMF + bootstrap | low | |
| scRNA atlas integration | medium (1–2 GPU-day on A100; or CPU-only ~2–3 days) | scVI on 100–300K cells. |
| LOHHLA + Mutect2 | medium-high (~200–400 CPU-hr) | biggest variable. |
| NetMHCpan / NetMHCIIpan | low-medium | per-allele × per-peptide. |
| DIAL audit | low | bootstrap dominates. |

Total estimate: ~1500–2500 CPU-hr + ~2 GPU-days. Single Azure NC-series burst over 1–2 weeks suffices; the `v17_arcasHLA_korean_k2` $4.80 burst pattern generalizes.

---

## 4. Risk & contingency triggers

| Trigger | Contingency |
|---|---|
| dbGaP TCGA paired BAM not granted by Wk 4 | Drop LOHHLA from Module C; report HLA call frequencies + MAF-based neoantigen only; mark fig 4B as future-work. |
| <3 scRNA cohorts accessible | Run Module B on best single cohort; fig 3 narrows; sc signature lock proceeds with caveat. |
| <3 ICI reference cohorts accessible | DIAL audit narrows; fig 5 reduced; if <2, defer DIAL to follow-up paper, drop fig 5, retain Module E with unaudited signatures + explicit caveat. |
| MC3 MAF disagreement with Paper 1 driver calls | Halt Module C until reconciled with Paper 1 ETL. |
| Paper 1 ship slips past late August 2026 | Push Track B start; do NOT begin compute before Paper 1 ship. |
| User wants to claim ICI response prediction in thyroid mid-plan | Hard stop; re-decide framing; do not silently proceed. |

---

## 5. Manuscript venue ladder (target ladder, decided at Wk 10)

- **Reach** — *Nature Cancer* / *Cancer Cell* — only if DIAL produces a thyroid-specific failure-mode story AND HLA-intact + neo+ + dedifferentiated subgroup is well-defined AND spatial cohort validates TLS/CXCL13 niche.
- **Default base** — *Nature Communications* / *Cell Reports Medicine* / *JCI Insight* — the registered design supports this tier even with degraded scope.
- **Safe** — *npj Precision Oncology* / *Genome Medicine* — fallback if multiple G-gates fail.

Avoid Frontiers/MDPI per `v17_graves_pivot`. ASCO / SITC abstract spin-off optional after submission.

---

## 6. Deliverables checklist (Track B end of Wk 12)

- [ ] `paper3_ici_dataset_registry_VERIFIED.tsv`
- [ ] `paper3_ici_signature_registry_VERIFIED.tsv`
- [ ] `paper3_ici_bulk_immune_ecotype_summary.tsv`
- [ ] `scrna_atlas/atlas_obs.parquet`, `atlas_markers.tsv`, `sc_signatures_locked.json`, `state_frequency_per_patient.tsv`
- [ ] `hla/hla_class1_calls.tsv`, `hla_class2_calls.tsv`, `lohhla_calls.tsv`
- [ ] `neo/neoantigens_class1.tsv`, `neoantigens_class2.tsv`
- [ ] `dial/dial_verdict.tsv`, `audit_report.md`
- [ ] `integrated/integrated_readiness_score.tsv`, `hla_intact_neoag_pos_subgroup.tsv`
- [ ] Manuscript v2 + supp + cover letter draft

---

## 7. Track A → Track B handoff

This 12-week plan, the registries, the atlas plan, the HLA/neoantigen feasibility, the DIAL plan, the figure plan, and the go/no-go verdict are the complete design freeze. Track B begins by re-reading these documents, executing G1–G6 closure, and starting Wk 1.

---

Track A completed. No marathon violation.
