# Manuscript strategy — RAI Response Genomics Atlas

## Working titles (rank in increasing risk)

1. **"A Thyroid Differentiation-Silencing Axis Identifies Radioiodine-Refractory Thyroid Cancer Across Independent Cohorts"** — safest. Claim is *identify*, not *predict*. Plural cohorts emphasized.

2. **"An Eight-Gene Iodide-Handling Panel Reveals a Molecular Gray Zone of Radioiodine Failure in Thyroid Cancer"** — sharper framing on the gray zone (Mu 2024 4-class data). Slightly more committed to the panel as the primary instrument.

3. **"Integrated Transcriptomic Evidence for Differentiation-Linked Radioiodine Failure in Thyroid Cancer"** — most conservative; deemphasizes the panel and emphasizes the integration.

Working choice for first draft: **title 1**, with title 2 saved as Cover Letter framing.

## Core claims (ranked safest → boldest)

1. The 8-gene panel is *associated with* RAI avidity in multiple Tier-2 cohorts (GSE151179 + GSE299988).
2. Driver-stratified analysis shows the panel adds information beyond BRAF/RAS/TERT status.
3. The 4-class uptake patterns of Mu 2024 are explained by a gradient of panel silencing — a **molecular gray zone**, not a binary.
4. Redifferentiation-treatment cohorts (Tier 5) show panel score rises with restored uptake, supporting mechanistic plausibility.
5. *Risk-stratification readout* framing, not clinical biomarker claim.

Bold (reserve): a single integrated panel score that **summarizes** TDS-16 with parsimony advantage. Hold this for review response, not headline.

## Figure plan (manuscript layout)

### Figure 1 — Biological model

Schematic of RAI as a multi-gate program: blood → NIS uptake → lineage TF maintenance → organification (TPO/DUOX) → storage (TG) → trafficking/colloid → radiation killing. Each gate annotated with the 8 panel genes that read it out. Cited references: TCGA-THCA Cell 2014; Landa 2016 JCI; Ho 2013 NEJM (selumetinib redifferentiation precedent).

### Figure 2 — Discovery cohort: TCGA-THCA + integrated bulk

Panels (a–d):

- (a) panel score distribution across drivers (BRAF / RAS / driver-neg) in TCGA-THCA (R17 L1).
- (b) cross-cohort forest of panel effect size (TCGA + Lee 2024 + Mun 2025 protein) — re-use Paper 1's 14-entry master forest.
- (c) HM450 β heatmap by panel zone (R17 L8).
- (d) per-zone 8-gene profile heatmap (R17 L7).

### Figure 3 — Label-anchored validation (NEW IN THIS WORKSPACE)

Panels (a–d):

- (a) GSE151179: panel score boxplot RAI-avid vs RAI-refractory; Wilcoxon p; per-sample-type breakdown.
- (b) GSE151179: ROC curve; AUC with 95% CI.
- (c) GSE299988: replication boxplot + ROC (caveat: small n).
- (d) Mu 2024 JCEM: mutation-frequency stacked bar by 4 uptake classes; arrow showing predicted panel direction.

### Figure 4 — Genomic context and molecular gray zone

Panels (a–d):

- (a) Driver × panel-zone 2D heatmap (BRAF / RAS / driver-neg × silenced / preserved). R17 L1 directly.
- (b) TERT+ enrichment per panel zone (R17 L5; OR=2.34, p=0.016).
- (c) PFI Cox forest per panel zone (R17 L6).
- (d) "Molecular gray zone" diagram — Mu 2024 4-class mapped to panel score continuum.

### Figure 5 — Redifferentiation and therapeutic translation

Panels (a–c):

- (a) Selumetinib (Ho 2013) — panel-gene RNA increase from pre to post (if data accessible).
- (b) MERAIODE / dabrafenib+trametinib post-treatment panel score.
- (c) Conceptual decision diagram: high-silencing patient → consider MAPKi redifferentiation before RAI; low-silencing patient → proceed to RAI directly.

### Supplement

- scRNA cellular substrate (Lu 2023, R17 L2) — 38.3% ATC dark-matter cells.
- Korean PTC+HT cross-ethnic replication (GSE286332, R17 L3).
- Proteomic corroboration (Mun 2025, R17 L4) — Fisher OR=8.54, p=2.7e-15.
- Methylation deep-dive (R17 L8).
- Leave-one-gene-out + random-panel permutation.

## Target journals (manuscript-strategy)

- **Tier A (NCx-floor reach)**: *Nature Communications* — if Tier-1 Boucai + Mu access succeed.
- **Tier B (NCx floor, more achievable)**: *Cell Reports Medicine* — accepts integrative discovery papers with substantial mechanistic depth.
- **Tier C (fallback)**: *JCI Insight*, *Thyroid* (the field's home journal), *npj Precision Oncology*.

## Submission gating

Before submission:

- Tier 1 (Boucai) or Tier 2 (Mu 4-class) external label-anchored validation locked.
- Leave-one-gene-out + random-panel null analyses run and reported.
- Driver-stratified analysis (BRAF / RAS / driver-neg) computed in every dataset.
- Voice-protected sections (Discussion opening + Limitations + Q9 + Cover ¶1) written by author keyboard.
