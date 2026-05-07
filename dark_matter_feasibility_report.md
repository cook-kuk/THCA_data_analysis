# Feasibility report: BRAF/RAS-negative thyroid cancer molecular taxonomy

Working title: **Molecular taxonomy of BRAF/RAS-negative thyroid cancer reveals noncanonical driver classes and radioiodine-related differentiation states**

Date: 2026-05-06

Scope note: This report uses existing local outputs first. It does not edit Paper 1. It treats survival, recurrence, RAI resistance, and therapeutic target discovery as limited or proposed unless directly supported by local outputs or public dataset metadata.

## Executive verdict

**Go, with a restrained molecular-taxonomy frame.** Existing outputs support a feasible paper centered on BRAF/RAS-negative thyroid cancer prevalence, noncanonical driver enrichment, TERT-only/co-mutation subclasses, sparse fusion overlay, and a thyroid differentiation/RAI-gene expression axis. The evidence does **not** support a primary survival paper, a proven RAI-resistance paper, or a therapeutic target discovery paper.

Best defensible claim:

> In TCGA-THCA, BRAF/RAS-negative thyroid cancers are common enough for analysis and contain separable molecular classes, including DICER1/EIF1AX/PPM1D-enriched FVPTC-like tumors, TERT-only high-risk candidates, fusion/alternative MAPK candidates, and a large true driver-negative residual class. Their expression state tracks thyroid differentiation and RAI-handling genes, but clinical RAI-response validation must be performed externally.

## What is already proven locally

### 1. TCGA BRAF/RAS-negative count

Primary clean denominator:

- `project/results/dark_matter_phase1/step1_to_6_summary.json`
- TCGA-THCA mutation-callable cohort: **482 tumors**
- BRAF V600E: **291**
- RAS hotspot: **55**
- BRAF/RAS-negative "Dark Matter": **137/482 = 28.42%**
- Clustered BRAF/RAS-negative cases: **136**, with DM1 **81** and DM2 **55**

Broader v17/TERT-integrated denominator:

- `project/results/v17/tables/quad_group_summary.tsv`
- v17 sample master after broader integration: **513 rows**
- Triple-negative group: **178**
- This is useful for extended bookkeeping, but the report should use **137/482** as the discovery denominator because it comes from the explicitly mutation-callable TCGA cohort.

### 2. DICER1 / EIF1AX / PPM1D enrichment

Local evidence:

- `project/results/dark_matter_phase1/step1_to_6_summary.json`
- `project/results/dark_matter_phase1/step3_dm_cluster_alt_driver_enrichment.tsv`

Within BRAF/RAS-negative TCGA:

| Driver group | DM1 | DM2 | Frequency contrast | Fisher p |
|---|---:|---:|---:|---:|
| DICER1 / EIF1AX / PPM1D | 1/81 | 6/55 | 1.2% vs 10.9% | 0.0175 |

Interpretation: DM2 is enriched for DICER1/EIF1AX/PPM1D. This supports a **noncanonical driver class** claim. It does not prove functional dependency or drug sensitivity.

Korean bulk validation already in Phase 2:

- `project/results/dark_matter_phase2/phase2_D1_verdict_AND_paper_outline.md`
- Yoo 2016/K2 parsed cohort: BRAF/RAS-negative **68/180 = 37.8%**
- DICER1 + EIF1AX in K2 BRAF/RAS-negative cases: **7/68 = 10.3%**, close to TCGA DM2 rate **10.9%**

### 3. TERT-only and co-mutation classes

Local evidence:

- `project/results/v17_tert_recovery/v2/FINAL_extended_audit.md`
- `project/results/v17_tert_recovery/v3/v3_summary.json`
- `project/results/dark_matter_phase2/p2d_driver_class_phenotype.tsv`

TERT recovery:

- Recovered TCGA-THCA TERT promoter mutated patients: **36**
- v3 quint groups:
  - BRAF V600E / TERT-: **255**
  - RAS / TERT-: **48**
  - TripleNeg / TERT-: **173**
  - BRAF V600E / TERT+: **25**
  - RAS / TERT+: **6**
  - TripleNeg / TERT+: **5**

Within the 482-patient Phase 2 driver-class map:

| Class | n | DM1 | DM2 | PFI events / n |
|---|---:|---:|---:|---:|
| Class4_DICER1_EIF1AX | 7 | 1 | 6 | 1/7 |
| Class5_TERT_only | 5 | 4 | 1 | 2/5 |
| Class6_True_driver_neg | 125 | 76 | 48 | 6/125 |

Interpretation: TERT-only cases are a small but important subclass. Because n=5, they should be shown descriptively or as a hypothesis-generating stratum, not as a definitive prognostic class.

### 4. Fusion overlay counts

Local evidence:

- `project/results/v17/tables/dark_matter_fusion_overlay.tsv`
- `project/results/v17/tables/driver_landscape_v17_summary.tsv`

In v17 dark-matter cluster overlay:

| Cluster | ALK_fusion | NTRK_fusion | RET_fusion | DICER1/EIF1AX/PPM1D | TP53 | unknown |
|---|---:|---:|---:|---:|---:|---:|
| DM1 | 1 | 1 | 1 | 1 | 0 | 103 |
| DM2 | 0 | 0 | 0 | 6 | 1 | 61 |

Interpretation: the overlay finds rare kinase-fusion candidates, mostly in DM1, but counts are too small and fusion ascertainment is incomplete. This supports a **fusion-overlay table**, not a fusion-driven subtype claim.

### 5. Thyroid differentiation score and 8-gene/RAI behavior

Local bulk evidence:

- `project/results/dark_matter_phase1/step2_rai_pivot.json`
- `project/results/v17/tables/dark_matter_cluster_clinical.tsv`

Within BRAF/RAS-negative TCGA:

| Metric | DM1 mean | DM2 mean | p value | Effect |
|---|---:|---:|---:|---:|
| RAI uptake gene score | 7.73 | 9.19 | 4.57e-16 | Cohen d magnitude 1.54 |
| TDS16 | 6.93 | 8.11 | 3.96e-16 | strong separation |

v17 broader dark-matter cohort:

| Cluster | n | TDS16 mean | RAI score mean | BRAF-like | RAS-like | unknown |
|---|---:|---:|---:|---:|---:|---:|
| DM1 | 109 | 6.81 | 7.58 | 86 | 18 | 5 |
| DM2 | 69 | 8.13 | 9.20 | 25 | 39 | 5 |

Single-cell behavior already shown:

- `project/results/dark_matter_phase1/phase1_FINAL_verdict.md`
- GSE241184: 8-gene score vs FVPTC signature in tumor thyrocytes **r = 0.905**
- `project/results/dark_matter_phase2/p2a2_summary.json`
- GSE184362: pooled 8-gene vs FVPTC signature **r = 0.889**, median patient r **0.865**
- `project/results/dark_matter_phase2/p2a2_multisite_summary.json`
- GSE184362 multisite: pooled r **0.914**
- `project/results/audit_2026_04_29/sc_wrapup/sc_wrapup_summary.json`
- GSE193581: signature is mainly thyrocyte-like/malignant-cell intrinsic, not immune/stromal.

Interpretation: the expression-state claim is strong: DM2 is more thyroid-differentiated and RAI-gene-high. The clinical RAI-response claim remains proposed until tested against RAI-avid/refractory labels.

### 6. Survival and recurrence limitations

Local evidence:

- `project/results/dark_matter_phase1/step2_pfi_results.json`
- `project/results/dark_matter_phase1/phase1_FINAL_verdict.md`
- `project/results/dark_matter_phase2/web/data/p10_enrich_summary.json`

Within BRAF/RAS-negative TCGA:

| Endpoint | n | Events | HR DM2 vs DM1 | 95% CI | p |
|---|---:|---:|---:|---:|---:|
| PFI | 136 | 9 | 1.33 | 0.35-4.95 | 0.675 |
| DFI | 103 | 3 | 3.16 | 0.29-34.98 | 0.348 |
| DSS | 134 | 3 | 0.73 | 0.07-8.06 | 0.797 |
| OS | 136 | 6 | 0.83 | 0.15-4.51 | 0.824 |

Power output:

- For HR 2.0, about **68 events** would be needed for 80% power.
- The local BRAF/RAS-negative PFI analysis has **9 events**.

Interpretation: survival/recurrence can be shown only as secondary exploratory context or as a reason the paper is not a prognostic paper.

## Public dataset registry summary

Full registry: `thyroid_public_dataset_registry.tsv`

High-value datasets for this paper:

1. **TCGA-THCA / GDC**
   - Current GDC project metadata: 507 cases with transcriptome, methylation, copy-number, SNV, SV categories.
   - Supports discovery and multi-omics annotation, but not RAI response.
   - Source: https://portal.gdc.cancer.gov/projects/TCGA-THCA and https://gdc.cancer.gov/about-data/publications/thca_2014

2. **cBioPortal TCGA THCA Firehose and PanCancer Atlas**
   - Firehose: all 516, sequenced 405, RNA-seq V2 509, methylation 511.
   - PanCancer Atlas: all 500, sequenced 490, RNA-seq 498, complete 480, structural variants 135.
   - Useful for cross-checking mutation, CNA, expression, methylation, and clinical fields.
   - Source: https://www.cbioportal.org/

3. **GSE151179 / GSE151180**
   - Best direct public RAI-response validation series.
   - GSE151179: mRNA array, 52 samples from 39 PTCs plus 13 matched non-neoplastic thyroids.
   - GSE151180: miRNA array subseries, 47/52 samples via SuperSeries GSE151181.
   - Includes RAI uptake/response stratification and PTC driver mutation/fusion characterization by PTC-MA.
   - Source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151179 and https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151181

4. **GSE184362**
   - scRNA-seq, 158,577 cells from 11 PTC patients and 23 specimens, including primary tumors, paratumors, lymph nodes, and RAI-refractory distant metastases.
   - Best public cellular-gradient validation.
   - Source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184362

5. **GSE232237**
   - Korean scRNA-seq thyroid cancer dataset; GEO summary lists 3 normal, 7 PTC, and 5 ATC cases, with 12 GEO samples.
   - Useful for Korean cellular-gradient/dedifferentiation validation.
   - Source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE232237

Blocked or limited:

- **KCG/K-BDS**: KCG has a thyroid/thymoma RNA-seq project (KAP210106; 444 samples, 2,392 data items), but thyroid-only extraction and clinical endpoints need login/access confirmation.
- **EGAD00001004845**: highly relevant Korean advanced thyroid dataset, but controlled access.
- **CODA**: portal exists for clinical/omics sharing, but no thyroid-specific CODA accession was verified in public search during this pass.
- **NCC/KCCR/NEST registry**: excellent epidemiology/treatment metadata, no molecular data and no driver/RAI-response endpoint suitable for this paper's molecular claims.

## Proposed analysis plan

### A. Discovery cohort definition

Use TCGA-THCA as the discovery cohort with the clean mutation-callable denominator.

1. Start from TCGA-THCA tumor samples with matched expression and mutation calls.
2. Exclude samples with BRAF V600E or canonical RAS hotspot mutations.
3. Keep non-V600E BRAF alterations separate unless explicitly classified as BRAF-like or RAS-like by prior curated rules.
4. Record TERT promoter status from the recovered TERT integration table.
5. Record alternative MAPK/fusion calls only where supported by curated local calls or public structural-variant evidence.

Primary count to report: **137/482 BRAF/RAS-negative TCGA tumors**.

### B. Subclassify BRAF/RAS-negative cases

Use a hierarchical class system:

1. **DICER1/EIF1AX/PPM1D class**
   - Classify if any DICER1, EIF1AX, or PPM1D mutation is present.
   - Report enrichment in DM2/FVPTC-like state.

2. **TERT class**
   - Separate TERT-only from BRAF+TERT and RAS+TERT in broader maps.
   - In the strict BRAF/RAS-negative cohort, treat TripleNeg/TERT+ as a small descriptive class.

3. **Fusion/alternative MAPK class**
   - RET, NTRK, ALK fusions and other curated kinase alterations.
   - Report as sparse overlay, not as a proven subtype.

4. **True driver-negative residual**
   - No BRAF V600E, no RAS hotspot, no DICER1/EIF1AX/PPM1D, no TERT promoter, no supported fusion/alternative MAPK call.
   - This is the unresolved discovery population.

### C. Compare differentiation and RAI-gene expression

Primary expression endpoints:

- TDS16
- Existing 8-gene score
- RAI-gene score using at minimum:
  - SLC5A5
  - TPO
  - TSHR
  - TG
  - PAX8
  - NKX2-1

Analyses:

1. Compare scores across BRAF/RAS-negative subclasses.
2. Plot class-wise distributions with effect sizes and confidence intervals.
3. Report individual gene behavior for the six RAI/thyroid-lineage genes.
4. Use "RAI-related differentiation state" language, not "RAI resistance" unless external response labels support it.

### D. Validate RAI-related signature in GSE151179/GSE151180

Use GSE151179 as the direct public RAI clinical validation dataset.

1. Recompute the six-gene RAI score and 8-gene score in mRNA array data.
2. Compare RAI-avid versus RAI-refractory labels from sample metadata.
3. If labels are lesion-level or timepoint-specific, preserve that structure and do not collapse without a rule.
4. Use GSE151180 miRNA as optional mechanism/context, not as a required driver of the main claim.
5. Report AUC/effect size as validation of a signature, not proof that TCGA subclasses are clinically RAI resistant.

### E. Validate cellular gradient in GSE184362 and GSE232237

GSE184362:

1. Use thyrocyte/malignant epithelial cells only.
2. Score 8-gene, FVPTC-like, cPTC-like, and six-gene RAI modules.
3. Test per-patient and pooled correlations.
4. Use specimen labels to compare primary tumor, paratumor, LN, and RAI-refractory distant metastasis descriptively.

GSE232237:

1. Use Korean scRNA-seq as an independent dedifferentiation/cell-state validation.
2. Compare normal/PTC/ATC cell-state gradients when metadata supports these labels.
3. Do not infer RAI response because GEO metadata does not provide that endpoint.

### F. Survival/recurrence as exploratory only

1. Keep TCGA PFI/OS/DFI/DSS as a limitations figure or supplementary table.
2. State that event scarcity prevents reliable survival inference.
3. Do not use survival as a go/no-go gate for the taxonomy paper.
4. Any recurrence analysis in external public microarray datasets should be treated as hypothesis-generating unless recurrence definitions and follow-up are explicit.

## Blocked items

| Item | Status | Reason |
|---|---|---|
| Clinical RAI resistance in TCGA subclasses | Blocked | TCGA does not provide direct RAI avid/refractory response labels suitable for this claim. |
| Survival/prognosis headline | Blocked | BRAF/RAS-negative TCGA has only 9 PFI and 6 OS events in local analysis. |
| Therapeutic target discovery | Blocked | No drug-response or functional dependency data in the requested evidence set. |
| Full fusion-driven taxonomy | Partially blocked | Local fusion overlay is sparse and incomplete; external fusion-rich cohort needed. |
| Korean advanced cohort reanalysis | Blocked pending access | EGAD00001004845 requires EGA DAC approval. |
| CODA thyroid-specific use | Blocked pending accession | Public search verified CODA as a portal but not a thyroid-specific CODA accession. |
| NCC/KCCR individual-level use | Blocked or metadata-only | Registry supports epidemiology/treatment context, not molecular subclass validation. |

## Figure and table plan

1. Cohort definition: TCGA BRAF/RAS-negative count and subclass flowchart.
2. Driver subclass map: DICER1/EIF1AX/PPM1D, TERT, fusion/alternative MAPK, true driver-negative.
3. Expression state: TDS16 and six RAI genes across subclasses.
4. GSE151179 validation: RAI-avid versus RAI-refractory signature behavior.
5. GSE184362/GSE232237 validation: single-cell gradient and per-patient correlations.
6. Limitations table: survival/recurrence event scarcity and unavailable endpoints.

## Claim boundary

Allowed:

- BRAF/RAS-negative thyroid cancers are a substantial TCGA subgroup.
- DICER1/EIF1AX/PPM1D alterations are enriched in the differentiated/FVPTC-like subgroup.
- TERT-only and fusion/alternative MAPK cases can be represented as small descriptive classes.
- The expression axis tracks thyroid differentiation and RAI-handling genes.
- Public datasets exist to validate RAI-related signatures and cellular gradients.

Not allowed:

- "This subclass predicts survival."
- "This subclass is RAI resistant" without GSE151179/151180 validation.
- "We discovered therapeutic targets" without drug-response or functional data.
- "Fusion-driven dark matter subtype" from the current sparse overlay.

