# Paper 3 ICI — Signature Registry

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
