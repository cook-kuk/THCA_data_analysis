---
title: DM1 canonical project state
date: 2026-07-30
status: wet-lab stopped; in-silico completion and Nature Communications submission target
branch: paper9-perturbation-extension-20260506
---

# Canonical project state

## Study identity

**Working title**  
*A thyroid-lineage state predicts radioiodine-refractoriness in BRAF V600E-mutant papillary thyroid cancer*

**Primary target**  
Nature Communications.

**Latest known draft**  
`NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md`

**Figure architecture note**  
`REVIEWER_FIRST_NC_REBUILD_2026_07_08.md`

## Clinical motivation

The intended clinical motivation is to reduce futile, repeated high-dose radioiodine exposure in patients unlikely to benefit. This manuscript does **not** establish a treatment rule and does **not** justify withholding RAI. The manuscript must frame this as a candidate biological state requiring pretreatment and prospective validation.

## Eight-gene panel

| Gene | Clinical alias | Functional category |
|---|---|---|
| TG | Thyroglobulin | hormone precursor and storage |
| TPO | Thyroperoxidase | iodine oxidation and organification |
| TSHR | TSH receptor | TSH signalling |
| SLC5A5 | NIS | membrane iodide transport |
| DIO1 | Type 1 deiodinase | thyroid hormone metabolism |
| PAX8 | — | thyroid-lineage transcription factor |
| NKX2-1 | TTF-1 | thyroid-lineage transcription factor |
| FOXE1 | TTF-2 | thyroid-lineage transcription factor |

Current labels:

- DM1: iodine-handling-low / lineage-low / RAI-refractory-like state.
- DM2: iodine-handling-high / lineage-high / RAI-responsive-like state.

Use “RAI-refractory-like” only as a biological alignment, not a measured clinical phenotype in discovery data.

## Wet-lab status

As of 2026-07-30, internal validation is not proceeding.

- TSO500 is DNA-oriented and did not provide usable RNA expression validation.
- Routine RNA-seq validation was not feasible.
- Research IHC antibodies raised reliability concerns.
- TG, PAX8, and NKX2-1 are clinically used antibodies, but no internal retrospective IHC validation was completed.

Therefore, phrases such as “SNUBH validation pending,” “validation underway,” or “will be added during review” are inaccurate unless the author explicitly reopens the experiment.

## Completed analysis domains

1. TCGA-THCA discovery and driver-orthogonal state analysis.
2. Driver-negative “dark matter” rescue.
3. Overall survival and progression-free interval analyses.
4. Kinase-fusion enrichment.
5. HM450 methylation analysis.
6. Multi-cohort external transcriptomic, proteomic, and single-cell concordance.
7. Post-RAI-refractory alignment using GSE151179.
8. In-silico three-marker IHC proxy and power simulation.

## Remaining computational work

- Fig. 3C: paired HM450 beta versus RNA expression scatter for TPO, DIO1, TSHR, and TG.
- Fig. 3D: driver-class methylation visualization.
- ATA intermediate-risk mosaic or replacement panel, only if the exact risk definitions and denominators are recoverable.
- Resolve the complete main-figure architecture and numbering.
- Assemble final figures, captions, and source-data manifest.

## Known narrative conflicts that must be resolved

1. Five-main-figure reviewer-first plan versus v2 multi-figure structure.
2. References to Fig. 7 and Fig. 8 without a clearly documented complete sequence.
3. “Predictive biomarker” wording versus a genotype-stratified PFI interaction without a treatment interaction.
4. “Internal validation pending” versus wet-lab termination.
5. Direct RAI prediction language versus absence of pretreatment uptake or response measurement.
6. “Nineteen cohorts,” “fourteen forest entries,” and other validation counts must be defined consistently.
7. The direction of HR=0.66 must be verified from model coding before saying DM1 or DM2 is protective.
8. “Extended Data” terminology should be checked against current Nature Communications requirements and replaced with Supplementary Figures if needed.

## Manuscript files known from project notes

```text
manuscript_v8/
├── NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md
├── NATURE_COMMUNICATIONS_FULL_DRAFT_2026_07_08.md
├── REVIEWER_FIRST_NC_REBUILD_2026_07_08.md
├── NATURE_CANCER_MANUSCRIPT_ARCHITECTURE_2026_07_08.md
├── 01_abstract.md
├── 04_results.md
├── 05_figure_captions_NC.md
├── 06_discussion.md
├── 07_star_methods.md
├── 08_cover_letter.md
├── 09_reviewer_qa.md
├── 13_supplementary_tables.md
├── figures/
└── dm1_story_web/public/figures/
```

Claude must inventory the real repository rather than assuming every path still exists.
