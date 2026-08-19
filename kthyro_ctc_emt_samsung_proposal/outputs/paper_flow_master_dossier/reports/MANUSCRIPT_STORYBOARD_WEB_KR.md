# Manuscript Storyboard Web

Purpose: 논문 본문 흐름에 맞춰 figure, data used, results scaffold, claim boundary를 한 줄로 고정한다.

Core boundary: 현재 public-data output은 CTC transition discovery가 아니라 tissue-state/marker/genetic-anchor prior framework다.

## Main Fig. 1. Study Logic And Data Provenance
- Results heading: A staged public-data prior map for a future serial CTC-EMT perturbation study
- Data used: TCGA-THCA public pilot; GSE250521 spatial summary; scRNA/proteomic marker summaries; Yu 2024 published CTC feasibility; local proposal reports.
- Write-now paragraph: We first separated the evidence layers that can be analyzed now from the prospective serial blood experiment that remains to be performed. Public TCGA, spatial, single-cell, and proteomic resources support tissue-state and marker-prior logic; Yu 2024 supports feasibility of serial CTC-EMT measurement; hospital serial CTC/cfDNA/tissue NGS data are required for the actual transition test.
- Claim boundary: This is a framework/provenance figure. It does not claim that public data measured postoperative CTC transition.
- Reviewer answer: The study is honest because it explicitly labels which evidence is public prior, which is published feasibility, and which is future prospective validation.

## Main Fig. 2. Driver Mutation Is Not Sufficient To Define State
- Results heading: BRAF-mutant PTC remains state-heterogeneous in public tumor data
- Data used: TCGA-THCA bulk tumor RNA/mutation/clinical outputs; BRAF-only vulnerability label distribution; driver/state variance models.
- Write-now paragraph: In the TCGA-THCA public pilot context, BRAF-mutant tumors distributed across multiple state/vulnerability labels, with the largest BRAF label accounting for only 33.2% of the BRAF-mutant subset. This supports the need to interpret postoperative blood signals with phenotype-state information rather than driver genotype alone.
- Claim boundary: This justifies genotype-plus-state interpretation. It is not evidence that CTC EMT states predict recurrence.
- Reviewer answer: The result is narrow: mutation alone is insufficient as a state descriptor in tumor tissue; no blood biology is inferred directly.

## Main Fig. 3. Spatial Tissue-State Organization
- Results heading: Public spatial data support non-random tissue-state organization
- Data used: GSE250521 spatial transcriptomics summary; 16 slides; 57,144 spots; spatial coherence, adjacency, hotspot, and interface outputs.
- Write-now paragraph: Spatial thyroid tissue data showed organized tissue-state structure rather than random mixing. In the current package, the spatial layer is used to support ROI and marker-selection logic for future CTC-EMT studies, not to prove the origin of circulating cells.
- Claim boundary: Spatial organization is a tissue-context prior only. It does not identify the shedding source of postoperative CTCs.
- Reviewer answer: The spatial analysis is deliberately positioned as marker/ROI rationale, not as direct CTC evidence.

## Main Fig. 4. Marker Module Support Across Data Layers
- Results heading: Candidate CTC-EMT marker modules have interpretable tissue and molecular context
- Data used: PTC scRNA marker-context summary; bulk proteomics proxy. The supplementary planning matrix (PF10) is shown separately and is not part of this figure.
- Write-now paragraph: Single-cell and proteomic summaries provide supporting context for epithelial, mesenchymal, thyroid-lineage, stress/survival, and immune-evasion marker modules. These data support rational panel construction for future CTC-EMT phenotyping and NGS anchoring. A separate supplementary planning matrix (PF10) records our hand-coded judgement of how useful each public/published evidence layer is for prior construction; that matrix is a planning score, not a data analysis result.
- Claim boundary: This is marker prioritization, not a validated clinical CTC assay and not a completed blood-based discovery. PF10 (supplementary) is a manuscript-planning judgement, not a cross-layer measurement.
- Reviewer answer: The marker panel is hypothesis-generating but biologically organized across independent public layers; the supplementary PF10 planning matrix records reviewer-visible judgement and does not introduce new data.

## Supplementary Planning Matrix. PF10 - Marker Module Planning Support Scores (not data-derived)
- Results heading: Supplementary, reviewer-visible record of manuscript-planning judgement on marker-module priors
- Data used: Hand-coded planning matrix marker_module_cross_layer_support_matrix.tsv; reference layers cited as context only (TCGA-THCA tissue [no CTC], GSE250521 spatial, scRNA module summary, bulk proteomics proxy, Yu 2024 published anchor summary).
- Write-now paragraph: As a supplementary planning matrix, we hand-coded (0/1/2) how useful each public/published evidence layer was for constructing each marker-module prior. These integers are manuscript-planning judgements, not measurements or validated statistics. The TCGA column records judgement about TCGA-THCA tumor tissue context only (TCGA contains no CTC, blood, or liquid biopsy data). The Yu 2024 column scores the published Yu 2024 paper's reported feasibility summary (62 PTC, 87% CTC detection); no unpublished or raw Yu CTC data was used.
- Claim boundary: These are planning scores, not validation metrics, not assay-performance values, and not a cross-layer measurement. TCGA does not contain CTC data and was never analyzed as such here.
- Reviewer answer: We disclose this matrix as a supplementary planning artifact so reviewers can see the prior-construction judgement explicitly, and so it cannot be misread as a data-derived cross-layer analysis.

## Main Fig. 5. Phenotype-To-Genotype CTC-EMT-NGS Framework
- Results heading: A matched genetic anchor is required to interpret postoperative blood signals
- Data used: Marker-to-NGS map; CTC-EMT-NGS prior panel; claim-boundary and Samsung proposal logic.
- Write-now paragraph: We define a research framework in which epithelial, hybrid E/M, and mesenchymal CTC phenotypes are interpreted together with matched tumor-normal NGS and serial cfDNA. The genetic anchor is essential because CTC phenotype alone may be nonspecific in early PTC.
- Claim boundary: The framework does not assume cfDNA or CTC-enriched NGS will work in all early PTC patients.
- Reviewer answer: This figure lowers overclaim risk by making genotype anchoring a requirement rather than an afterthought.

## Main Fig. 6. Missing Dataset And Prospective Unlock
- Results heading: No public dataset currently tests the serial postoperative CTC-EMT genetic-map question
- Data used: Dataset-gap bridge; future data unlock map; claim boundary matrix; hospital CRF table.
- Write-now paragraph: We found no public dataset that jointly provides serial pre/post-thyroidectomy blood, CTC EMT phenotype transition, matched tissue-normal NGS, cfDNA or CTC-enriched NGS, and postoperative follow-up. This absence defines the scientific need for the prospective Samsung/Yu cohort.
- Claim boundary: The missing-dataset result justifies the prospective study; it is not a substitute for the prospective study.
- Reviewer answer: The proposal is needed because the exact joint dataset required for the central biological question does not exist publicly.

## Supplement/Defense. Reviewer Defense And Execution Control
- Results heading: Defense materials define the safe manuscript lane
- Data used: Reviewer risk heatmap; publishability decision table; execution board; analysis control package.
- Write-now paragraph: The supplementary defense layer maps expected reviewer attacks to figures, sets publishability boundaries, and preserves execution control through key-number, statistical-plan, and kill-criteria tables.
- Claim boundary: These are control/defense artifacts, not additional biological discoveries.
- Reviewer answer: The defense package is included to prevent claim creep and to make the paper auditable.
