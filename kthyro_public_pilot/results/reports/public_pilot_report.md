# K-Thyro Public-Data Pilot Report

## Key Findings

- TCGA-THCA bulk expression was scored at patient level for RAI differentiation, HLA/APM immune visibility, cytotoxicity, myeloid/CAF barrier, hypoxia/drug-delivery proxy, proliferation, and derived vulnerability axes.
- Patient-level vulnerability labels were generated for 505 TCGA primary tumor samples. Label counts: {'RAI-readable differentiated': 127, 'Mixed/Other': 125, 'Drug-delivery barrier-high': 105, 'HLA-visible inflamed': 44, 'CD8-excluded myeloid/CAF-high': 42, 'HLA-low immune-invisible': 32}.
- Public GSE250521 spatial transcriptomics was scored at spot level for 16 slides, producing spatial niche maps and slide-level summaries.
- KNN same-niche permutation analysis found 16 slides with niche coherence z-score > 2. Spots are treated as nested within slides; no fake large-N biological claims are made.
- Drug/perturbation candidates were ranked from local DepMap/PRISM-derived public resources. These nominate validation hypotheses only.

## GO/NO-GO

# GO/NO-GO Decision

Decision: **GO**

Reason: 5/5 main axes have TCGA support and spatial same-niche coherence is non-random in cancer slides.

Testable hypotheses:
- MEK/BRAF/RET/NTRK-axis inhibition can restore iodide uptake in RAI-low/MAPK-high tumor niches.
- HLA/APM-low PAX8+ territories can be rescued by IFN/epigenetic perturbation and validated by HLA-I/B2M/TAP1 protein.
- CAF/myeloid-rich regions spatially exclude CD8/GZMB and can be nominated for barrier-modulating combinations.
- ECM/hypoxia-high vascular-low territories predict poor fluorescent drug/nanoparticle penetration ex vivo.

Claim boundary: this is hypothesis-generating public-data evidence, not clinical deployment evidence.


## Claim Boundaries

- Public pilot supports hypothesis generation, not clinical deployment.
- TCGA bulk cannot prove spatial heterogeneity.
- Spatial transcriptomics supports spatial expression states, not direct peptide presentation or patient-specific HLA presentation.
- Drug-delivery failure score is a proxy unless validated by fluorescent drug/nanoparticle imaging.
- RAI-restorable score is a hypothesis unless iodide uptake assay validates it.
- HLA/APM-low state requires protein validation by mIHC/IF or GeoMx.

## Proposed Experimental Validation Panel

- Tumor/thyroid: PAX8, TG, TPO, NIS/SLC5A5.
- HLA/APM: HLA-I, B2M, TAP1.
- T cells: CD8, GZMB.
- Myeloid: CD68, CD163.
- CAF/ECM: ACTA2, FAP, COL1A1.
- Checkpoint: PD-L1.
- Hypoxia/drug delivery: CA9, VEGFA plus fluorescent drug/liposome/nanoparticle distribution imaging.

## What To Ask Hospital/Pathology Collaborators

- Can FFPE blocks be selected across PTC, locally advanced PTC, RAI-refractory disease, PDTC/ATC, and matched normal/adjacent thyroid?
- Are pre/post-operative blood, RAI treatment history, recurrence, and response annotations available?
- Can fresh tissue be routed to slice/organoid perturbation and fluorescent delivery imaging within hours of surgery?
- Can pathology annotate tumor, stromal, lymphoid, necrotic/hypoxic, and invasive-front ROIs for GeoMx/mIHC?
- Are paired WES/RNA-seq/HLA typing feasible for a subset to avoid overinterpreting RNA-only immune visibility?

## Candidate Drug/Perturbation Hypotheses

- AZD7762 DPC-000855: Proliferation stress (A: thyroid line sensitivity + mechanistic relevance).
- JAK3_7406 DPC-007447: HLA/APM restoration / IFN-epigenetic (A: thyroid line sensitivity + mechanistic relevance).
- JAK3_7406 DPC-007447: HLA/APM restoration / IFN-epigenetic (A: thyroid line sensitivity + mechanistic relevance).
- JAK3_7406 DPC-007447: HLA/APM restoration / IFN-epigenetic (A: thyroid line sensitivity + mechanistic relevance).
- ALISERTIB DPC-000428: Proliferation stress (A: thyroid line sensitivity + mechanistic relevance).
- BAY-HDAC11_4 DPC-007395: HLA/APM restoration / IFN-epigenetic (A: thyroid line sensitivity + mechanistic relevance).
- BAY-HDAC11_4 DPC-007395: HLA/APM restoration / IFN-epigenetic (A: thyroid line sensitivity + mechanistic relevance).
- MYCOPHENOLIC ACID DPC-007456: Proliferation stress (A: thyroid line sensitivity + mechanistic relevance).
- SELUMETINIB:NAVITOCLAX (8:1 MOL/MOL) DPC-005903: RAI redifferentiation / MAPK (A: thyroid line sensitivity + mechanistic relevance).
- MOMELOTINIB DPC-004334: HLA/APM restoration / IFN-epigenetic (A: thyroid line sensitivity + mechanistic relevance).

## Main Output Files

- `results/tables/tcga_thca_patient_vulnerability_scores.tsv`
- `results/tables/spatial_spot_vulnerability_scores.tsv`
- `results/tables/spatial_slide_niche_summary.tsv`
- `results/tables/drug_pilot_candidate_rankings.tsv`
- `results/tables/integrated_public_pilot_evidence_matrix.tsv`
- `results/figures/proposal`
