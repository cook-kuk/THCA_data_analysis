---
title: "Reviewer-first Nature Communications rebuild — Paper 1 DM1"
date: 2026-07-08
author: Codex reviewer-mode scaffold
status: implementation blueprint; no voice-protected prose generated
source_context:
  - project/manuscript_v8/EDITORIAL_REVIEW_RUTHLESS_CUT_2026_06_05.md
  - project/manuscript_v8/MANUSCRIPT_REORG_6FIG_PLAN_2026_06_05.md
  - project/manuscript_v8/05_figure_captions_NC.md
  - project/papers_hub_2026_0509/paper1_nature_cancer_board.html
---

# Reviewer-First Nature Communications Rebuild

## Executive Verdict

The paper should not be presented as:

> "Here is an 8-gene panel that predicts RAI failure."

That version gets attacked immediately: panel selection, retrospective endpoints, no wet-lab causality, too many cohorts, over-translation.

The paper should be presented as:

> "A hidden thyroid-lineage state exists in papillary thyroid cancer, is not explained by canonical drivers, converges with known dedifferentiation biology, has an epigenetic correlate, reproduces across independent cohorts, and is positioned for prospective RAI harm-avoidance validation."

This is the version a Nature Communications reviewer can follow.

## Reviewer Attack Map

| Reviewer attack | Main figure that must answer it | Required answer |
|---|---|---|
| "Your 8 genes created the cluster." | Figure 1 | Pan-genome recovery ARI = 0.92 and driver-only ARI = -0.007. |
| "This is just BRAF/RAS/fusion biology." | Figure 1 + Figure 3 | Driver transcript AUC near chance; methylation and DM state cross driver classes. |
| "What does DM1 biologically mean?" | Figure 2 | It is thyroid-lineage / iodine-handling loss, converging with TDS-16 and Landa PDTC/ATC silencing. |
| "Methylation is correlative." | Figure 3 + limitations | State "epigenetic correlate"; do not claim causality. |
| "TCGA-only cherry-pick." | Figure 4 | Independent Korean PTC, advanced thyroid, proteome, and compressed forest. |
| "Clinical claim is too strong." | Figure 5 | Retrospective survival, post-RAI state, ATA uncertainty, FFPE deployability; prospective validation boundary. |
| "Too many side analyses." | Extended Data only | Spatial, PRISM, pan-cancer, sub-A/B, single-cell expansion, decision curves stay out of main. |

## Main Claim Ladder

1. **Existence:** DM1/DM2 is a real transcriptomic state, not an artifact of eight genes.
2. **Biology:** DM1 means loss of thyroid lineage and iodine-handling machinery.
3. **Mechanism correlate:** DM1 shows promoter hypermethylation of lineage genes and is not reducible to one driver.
4. **Validation:** DM1 recurs across independent cohorts and modalities.
5. **Clinical relevance:** DM1 aligns with reduced differentiation, RAI-failure biology and adverse outcome, while remaining a prospective-validation candidate.

Every main figure must map to exactly one rung.

## New Five-Figure Architecture

### Figure 1 — State Discovery

**Title:** A driver-orthogonal thyroid-lineage axis separates papillary thyroid cancer into two reproducible states

**Reviewer job:** Defeat "panel artifact" and "driver proxy" before the reader reaches biology.

| Panel | Content | Current asset / rebuild instruction |
|---|---|---|
| 1A | Cohort and design strip: TCGA discovery, external validation, methylation, clinical layers | Rebuild simple schematic. Do not list 19 cohorts. |
| 1B | Ordered eight-gene heatmap sorted by P_DM1 with driver/stage tracks | Use heatmap from `Fig1_discovery_axis` or rebuild from master table. |
| 1C | Driver-neutrality bar: BRAF/TERT/KRAS/NRAS/HRAS AUC near chance | Use `F3_driver_neutrality.png`. |
| 1D | Pan-genome ARI ladder: eight-gene, TIERA67, top-5000, driver-only | Use `F4_pangenome_robustness.png`; simplify to four bars. |
| 1E | DM1 prevalence by driver class: BRAF V600E 0.7%, BRAF/RAS-negative 49%, RAS+ 96% | Rebuild stacked bar. |

**Delete from main:** UMAP, PCA, sub-A/B silhouette, 256-subset grid, quantum kernel.

**Figure legend claim:** The eight-gene axis is a compact readout of a broader driver-orthogonal state.

### Figure 2 — Biological Meaning

**Title:** DM1 represents coordinated loss of thyroid lineage and radioiodine machinery

**Reviewer job:** Explain what the state means biologically.

| Panel | Content | Current asset / rebuild instruction |
|---|---|---|
| 2A | Panel gene expression DM1 vs DM2: TF3 + effector5 | Rebuild clean box/violin plot. |
| 2B | Thyroid TF / RAI pathway enrichment | Use or rebuild from `q1_gsea_hallmark.png` and `q5_thyroid_tf_network.png`; simplify. |
| 2C | TDS-16 convergence / Panel-8 is the compact form of canonical differentiation loss | Use v13 outputs; rebuild as one convergence panel. |
| 2D | Landa 2016 PDTC/ATC overlap: 5/8 genes | Use `gse76039_mechanism_heatmap.png`; add overlap annotation. |
| 2E | Proteomic confirmation: Mun 2025 direction consistency | Use ED11/Mun panel as main if clean; otherwise rebuild. |

**Delete from main:** detailed single-cell, deconvolution, pseudotime, all pathway walls.

**Figure legend claim:** DM1 is not just a cluster; it is a thyroid differentiation-loss state.

### Figure 3 — Mechanism Correlate

**Title:** DM1 is associated with promoter hypermethylation of thyroid-lineage genes across driver contexts

**Reviewer job:** Give mechanistic depth without overclaiming causality.

| Panel | Content | Current asset / rebuild instruction |
|---|---|---|
| 3A | HM450 promoter beta heatmap, 8 genes + controls | Use `Fig3_epigenetic` / `Fig8_epigenetic`; simplify labels. |
| 3B | Mean eight-gene beta: DM1 0.385 vs DM2 0.253 | Use existing Fig3/Fig8 panel. |
| 3C | Beta-expression scatter for TPO/DIO1/TSHR/TG | Rebuild; this is currently the missing "mechanism bridge" panel. |
| 3D | Per-driver-class beta: BRAF V600E approx RET fusion much greater than RAS | Rebuild from v5A/per-driver class table. |
| 3E | Caption-only boundary box: "epigenetic correlate, not causal proof" | Add small inset, not a schematic therapy claim. |

**Delete from main:** PRISM, DepMap, HMA trial schematic, spatial Visium, pan-cancer.

**Figure legend claim:** Methylation supports a silencing model but does not prove causality.

### Figure 4 — External Validation

**Title:** The DM1 axis recurs across independent PTC, advanced-thyroid and proteomic cohorts

**Reviewer job:** Defeat TCGA cherry-pick with a small, interpretable validation panel.

| Panel | Content | Current asset / rebuild instruction |
|---|---|---|
| 4A | Lee 2024 Korean PTC, n=632, d=5.93 | Use external validation atlas or rebuild single cohort panel. |
| 4B | Landa 2016 PDTC/ATC aggressive-end alignment | Use GSE76039 lineage plot/heatmap. |
| 4C | Mun 2025 proteomics, n=336, 7/7 direction | Use ED11/Mun panel. |
| 4D | Compressed four-entry forest: TCGA, Lee, Landa, Mun | Rebuild from master forest, not the 14-entry version. |
| 4E | Full 14-cohort validation statement as text inset | "Full forest and 80/80 matrix in Supplement." |

**Delete from main:** K2 calibration-heavy panel, GPL570 four-cohort details, K-Thyro, 80-cell matrix, full external atlas.

**Figure legend claim:** Independent cohorts reproduce the state; the full validation battery is supplementary.

### Figure 5 — Clinical Relevance and Deployability

**Title:** DM1 is associated with adverse outcome, RAI-failure biology and FFPE deployability

**Reviewer job:** Show why anyone should care clinically while avoiding treatment-selection overclaim.

| Panel | Content | Current asset / rebuild instruction |
|---|---|---|
| 5A | Pooled OS forest: TCGA + MSK HR 2.53 [1.31, 4.89] | Use survival forest. |
| 5B | Multivariate Cox: DM1 adjusted with age/stage/TERT where possible | Use `o_multivariate_cox.png`. |
| 5C | Post-RAI refractory alignment: GSE151179 d approx -1.0 | Use/rebuild boxplot; include n and caveat. |
| 5D | ATA intermediate-risk enrichment / decision uncertainty | Rebuild mosaic. |
| 5E | FFPE vs fresh frozen concordance, KS p=0.44 | Use Fig5 portability panel. |

**Delete from main:** selpercatinib waterfall, prospective trial schema, decision curve, time-dependent ROC, calibration, HMA schematic.

**Figure legend claim:** Clinical relevance is retrospective and deployment-ready enough for prospective validation, not for immediate treatment selection.

## Extended Data Set

Keep only six Extended Data figures:

| ED | Content | Reason |
|---|---|---|
| ED1 | Full ordered heatmap + all annotation tracks | Detail for Fig. 1. |
| ED2 | Panel-size / leave-one-out / ARI sensitivity | Defense against "why eight genes." |
| ED3 | Single-cell thyrocyte-intrinsic validation | Defense against stromal/immune confound. |
| ED4 | Within-DM1 fusion-positive vs fusion-negative methylation | Defense for fusion-independent methylation, but small n. |
| ED5 | Full cross-cohort validation atlas + 80/80 matrix | Defense against cherry-pick. |
| ED6 | Survival diagnostics: KM, ROC, calibration, proportional-hazards checks | Clinical statistics support. |

Everything else becomes supplementary table/text or is removed from the submission package.

## Main Text Results Order

1. **A compact thyroid-lineage axis identifies a driver-orthogonal state in TCGA-THCA**  
   Figure 1

2. **DM1 represents coordinated loss of thyroid lineage and iodine-handling machinery**  
   Figure 2

3. **DM1 is associated with promoter hypermethylation of lineage genes across driver contexts**  
   Figure 3

4. **The DM1 axis is reproduced across independent PTC, advanced-thyroid and proteomic cohorts**  
   Figure 4

5. **DM1 is linked to adverse outcome, RAI-failure biology and FFPE deployability**  
   Figure 5

## Final Reviewer-Grade Central Claim

> We identify a driver-orthogonal thyroid-lineage state in papillary thyroid cancer that is recovered by pan-genome expression structure, corresponds to coordinated iodine-handling loss, has an epigenetic silencing correlate, recurs across independent cohorts and modalities, and is associated with retrospective RAI-failure biology and adverse outcome.

This is strong.

The following claim is not allowed:

> DM1 predicts RAI response and should guide treatment.

That is prospective-validation territory.

## Immediate Implementation Checklist

| Priority | Work item | Output |
|---|---|---|
| P0 | Stop using current 6-figure article as final structure | Replace web/PDF with reviewer-first five-figure version. |
| P0 | Rebuild Figure 1 from heatmap + driver-neutrality + ARI + driver-prevalence | One clean discovery figure. |
| P0 | Rebuild Figure 2 around biology, not validation | Gene expression, pathway, TDS, Landa, proteome. |
| P0 | Rebuild Figure 3 around methylation only | HM450, mean beta, beta-expression, per-driver beta. |
| P0 | Rebuild Figure 4 as compressed validation | Lee, Landa, Mun, four-entry forest. |
| P0 | Rebuild Figure 5 as clinical/deployability | OS, Cox, post-RAI, ATA, FFPE. |
| P1 | Convert all removed analyses into ED/Supp index | Reviewer defense without main-text clutter. |
| P1 | Regenerate Nature Communications article HTML/PDF after figures are rebuilt | Web deliverable. |

