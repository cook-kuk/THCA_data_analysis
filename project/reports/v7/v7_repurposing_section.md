# v7. Computational drug repurposing for novel THCA targets

## Methods

We assembled an end-to-end repurposing workflow around the eight top-novelty druggable targets identified in our THCA differential-expression and druggability screen (TACSTD2, TMPRSS4, PLEKHA6, CYP1B1, LDLR, GABRB2, B3GNT3, PTPRE). Drug-target evidence was triangulated from three public resources: (i) Open Targets Platform `knownDrugs` (GraphQL v4) after resolving each gene symbol to its Ensembl stable ID via the `search` query; (ii) DrugCentral `drug.target.interaction.tsv.gz` (cached locally, filtered to HGNC matches of the eight symbols); and (iii) ChEMBL `activity.json` with a `pchembl_value >= 6` threshold followed by a molecule-level fetch restricted to `max_phase >= 1` (i.e., at least phase-1 clinical development). Six commercially available thyroid-cancer tyrosine kinase inhibitors (lenvatinib, sorafenib, cabozantinib, vandetanib, selpercatinib, pralsetinib) were excluded so the pipeline surfaces genuinely repurposable assets.

Each (target, drug) record was scored with a weighted rubric: `score = 0.30 * (max_phase / 4) + 0.20 * potency + 0.15 * thyroid_safety + 0.15 * oral_bioavailability + 0.10 * off_patent + 0.10 * mechanism_fit`, with a 10x downweight for withdrawn drugs. Potency was mapped from pChEMBL bins (>=8 = 1.0, >=7 = 0.7, >=6 = 0.4, else 0.2; missing = 0.3). Thyroid safety was a text-derived heuristic flag on QT or hepatotoxicity keywords. Off-patent status used first-approval year cutoffs (<=2005, 2006-2015, >2015). Mechanism fit favored inhibitors/antagonists/blockers (all eight targets are overexpressed in THCA) and penalized agonists/activators. Where multiple sources reported the same drug-target pair we merged phase tags (maximum), pChEMBL values (maximum), and mechanism strings (union).

For the ten highest-scoring unique drugs we layered (a) ClinicalTrials.gov v2 study counts with a thyroid-indication flag, (b) ChEMBL-derived pharmacokinetic tags (oral/parenteral, first-approval year, withdrawal), and (c) a PubMed E-utilities count of papers since 2020 matching `<drug> AND (thyroid OR cancer)`. Combination hypotheses were enumerated across all C(8,2)=28 target pairs where both members had at least one candidate with score >=0.4; pathway overlap was drawn live from Reactome's UniProt-to-pathway mapping with a curated fallback for offline runs. A synergy score combined pathway overlap, FDA-approval status of both partners, and mechanism-class diversity.

## Top-10 repurposing candidates

- **cannabidiol** (target GABRB2, phase 4.0, score 0.720) — unspecified; ClinicalTrials.gov studies=50 (thyroid=0); PubMed 2020+ hits=1098
- **quercetin** (target CYP1B1, phase 3.0, score 0.640) — unspecified; ClinicalTrials.gov studies=50 (thyroid=0); PubMed 2020+ hits=5896
- **cannabinol** (target CYP1B1, phase 3.0, score 0.580) — unspecified; ClinicalTrials.gov studies=17 (thyroid=0); PubMed 2020+ hits=72
- **resveratrol** (target CYP1B1, phase 3.0, score 0.580) — unspecified; ClinicalTrials.gov studies=50 (thyroid=0); PubMed 2020+ hits=5550
- **luteolin** (target CYP1B1, phase 2.0, score 0.565) — unspecified; ClinicalTrials.gov studies=24 (thyroid=0); PubMed 2020+ hits=1589
- **etizolam** (target GABRB2, phase 1.0, score 0.550) — unspecified; ClinicalTrials.gov studies=1 (thyroid=0); PubMed 2020+ hits=5
- **flunitrazepam** (target GABRB2, phase 1.0, score 0.550) — unspecified; ClinicalTrials.gov studies=2 (thyroid=0); PubMed 2020+ hits=67
- **nitrazepam** (target GABRB2, phase 1.0, score 0.550) — unspecified; ClinicalTrials.gov studies=17 (thyroid=0); PubMed 2020+ hits=12
- **midazolam** (target GABRB2, phase 1.0, score 0.550) — unspecified; ClinicalTrials.gov studies=50 (thyroid=0); PubMed 2020+ hits=1088
- **triazolam** (target GABRB2, phase 1.0, score 0.550) — unspecified; ClinicalTrials.gov studies=18 (thyroid=0); PubMed 2020+ hits=21

## Combination hypotheses

- **TACSTD2:sacituzumab govitecan** + **CYP1B1:quercetin** (synergy=0.50, pathway overlap=0)
- **TACSTD2:sacituzumab govitecan** + **GABRB2:cannabidiol** (synergy=0.50, pathway overlap=0)
- **CYP1B1:quercetin** + **GABRB2:cannabidiol** (synergy=0.50, pathway overlap=0)

A complete 8x8 visualization of synergy scores (upper triangle) and Reactome pathway overlap counts (lower triangle) is provided as an interactive heatmap (`combination_heatmap.html`).

## Summary statistics

- Total drug-target interaction rows after deduplication: **83**
- FDA-approved (phase 4) hits: **2**
- Clinical phase 2-3 hits: **7**
- Rows with repurposing score > 0.5: **20**

## Figures and tables referenced

- Figure v7-1: interactive repurposing table (`v7_repurposing_ranked.tsv`).
- Figure v7-2: 8x8 combination heatmap (`combination_heatmap.html`).
- Table v7-A: `v7_combination_hypotheses.tsv`.

## Limitations

The analysis is entirely computational and public-data-driven: (1) bioactivity measurements in ChEMBL span heterogeneous assays, so the potency term is a coarse binning rather than calibrated Ki/IC50; (2) our thyroid-safety term is a keyword heuristic and cannot replace FDA label review, DDI profiling, or in vivo QT/hepatotoxicity work; (3) Open Targets phase tags occasionally lag ClinicalTrials.gov, which is why STEP 5 re-queries CT.gov; (4) Reactome overlap is a necessary but not sufficient proxy for pharmacodynamic synergy, and all combination hypotheses require wet-lab validation before translational claims; (5) the eight targets are transcriptionally overexpressed in THCA but proteomic/surface-exposure confirmation is a separate workstream. Future work should layer PDX efficacy signals, patient-derived-organoid sensitivity data, and structure-based DTI refinement for the novel targets without co-crystal structures.
