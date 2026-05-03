# Paper 3 ICI — Figure Plan

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
