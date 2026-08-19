# 00 — Executive readiness, four new public-data candidates

**Date** 2026-08-06 · total downloaded this round **7.1 MB** (budget was 5 GB soft / 10 GB hard)
· **no DECISION raw FASTQ was downloaded** · no manuscript file was modified · no commit made.

## One-sentence conclusion

The DECISION trial RNA-seq is the largest advanced radioiodine-refractory cohort in existence
and its per-sample treatment arm is fully recoverable, but **no processed expression matrix
exists anywhere in public**, so the panel cannot be scored without processing 0.98 TB of raw
reads — and even then the per-patient survival endpoints are not machine-readable.

## Readiness grades

| Dataset | Grade | One-line reason |
|---|---|---|
| **DECISION PRJNA563018** (n = 125) | **B · METADATA_READY_RAW_TOO_LARGE** | treatment arm, histology, age, sex, QC all recovered per sample; no processed matrix anywhere; per-patient PFS/RECIST not machine-readable |
| **Zheng 2024 metabolomics** (20 / 14) | **D · MANUAL_DOWNLOAD_REQUIRED** | supplement obtained but contains methods and a chemistry table only; the sample-level matrix sits in a Baidu repository requiring a browser |
| **ERRITI redifferentiation** (n = 20) | **C · ENDPOINT_REFERENCE_ONLY** | full supplement obtained; it is a gene-panel list, eligibility criteria and figure legends — no per-patient table, no molecular matrix |
| **Wang acetoacetate metabolomics** (24 / 18) | **D · EXCLUDE** | no per-patient data published; metabolomics cannot test the panel; **non-independent from Zheng — same laboratory and corresponding author** |
| **Redifferentiation trials pooled** (9 arms, 105 pts) | **A · ANALYSED** | the one positive of the round; Tier-2 uptake endpoint extracted from the ITOG 2025 statement and pooled |

## What was actually downloaded

| File | Size | Content |
|---|---|---|
| `decision_prjna563018/metadata/decision_runinfo.csv` | 62 KB | SRA RunInfo, 125 runs |
| `decision_prjna563018/metadata/decision_ena.tsv` | 14 KB | ENA report with read counts and FASTQ byte sizes |
| `decision_prjna563018/metadata/decision_biosample_metadata.tsv` | — | **125 BioSamples with per-sample attributes, merged** |
| `decision_prjna563018/source/decision_mct2020_supplementary.pdf` | 5.1 MB | AACR supplement, 8 pages, via figshare article 22506703 |
| `zheng_rair_metabolomics_2024/source/pmc_supp.zip` | 1.3 MB | Europe PMC supplementary bundle |
| `zheng_rair_metabolomics_2024/extracted/41598_2024_61067_MOESM1_ESM.docx` | 699 KB | supplementary methods + one chemistry table |
| `erriti_redifferentiation/source/erriti_39940925.docx` | 35 KB | ERRITI Supplemental Table S1 + eligibility criteria |
| `erriti_redifferentiation/source/erriti_39940928.docx` | 13 KB | ERRITI figure legends |

SHA-256 manifests written for all three new datasets.

---

## DECISION PRJNA563018 — what is and is not recoverable

**Submitter:** Vall d'Hebron Institute of Oncology, registered 2019-08-30.
**Paper:** Mol Cancer Ther 2020;19(1):312, PMID 31540966 — **not open access, no PMCID**.
Supplement obtained through the AACR figshare mirror (article 22506703), the same route that
worked for Boucai 2023 and ERRITI.

### Recovered per sample, from BioSample attributes (all 125)

| Field | Distribution |
|---|---|
| **Treatment arm** | **sorafenib (BAY 43-9006 400 MG BID) 68 · placebo 57** |
| Histology | papillary 73 · follicular 33 · **poorly differentiated 18** · non-diagnostic 1 |
| Sex | male 65 · female 60 |
| Age | median 63 (32–87) |
| Library barcodes | `isolate` = AB####, `GENOMIC.LAB.SAMPLE.ID` = B15/###, library = L15/#### |
| Provider | Jaume Capdevila, Vall d'Hebron, Barcelona |

The arm split reproduces the published 68/57 exactly. **This is the crucial recovery: a
57-patient untreated comparator arm is precisely what TCGA lacked** (2 events among
non-radioiodine-treated patients made the treatment-by-score interaction unestimable there).

### Recovered from the supplement, aggregate only

Table 2 suppl compares the RNA-seq subset (n = 125) with the full DECISION trial (n = 417):
ECOG 0 in 64%, BRAF mutated 33 (26.4%), RAS mutated 20 (16%), and expression profiles
**BRAF-like 56 · RAS-like 28 · NoBRaL 41**. Figure 2 suppl gives subset survival:
sorafenib median PFS 10.8 months versus placebo 3.7 months, HR 0.42 (95% CI 0.26–0.67),
P < 0.001.

### Two hard blockers

**1. Supplementary Table 1 is five pages of raster images at 87–93 ppi**, and its columns are
QC (GOOD/BAD), CNAG library barcode, genomic lab sample ID, library ID, clinical trial number,
subject number and treatment group. Every one of those except QC is already in BioSample.
Worse, **the Subject Number column has been destroyed by spreadsheet formatting** — the
nine-digit trial subject IDs render as `1E+08`, `1,2E+08`, `1,4E+08`. The key that would link
these samples to trial clinical data no longer exists in the published file.

**Consequence: per-patient PFS time, PFS event, RECIST response, and expression-cluster
assignment are not machine-readable and cannot be reconstructed from public files.**

**2. No processed expression matrix exists.** Searched and confirmed absent: GEO (zero hits
for the BioProject and for the study terms), ArrayExpress/BioStudies, Zenodo, Figshare beyond
the single supplement PDF, GitHub (zero repositories), and the BioProject's own linked
resources. The only deposited data is raw reads.

**Raw volume: 0.98 TB across 125 runs** (median ~66 M reads, ~10 Gbase per sample).

### Pilot cost estimate, not executed

| Quantity | Estimate |
|---|---|
| Per-sample FASTQ | ~7.8 GB (paired, two files) |
| 3-sample pilot download | ~23 GB |
| 5-sample pilot download | ~39 GB |
| Full 125-sample download | **~980 GB** |
| Quantification approach | pseudoalignment (kallisto/salmon) against gencode — the eight panel genes only need transcript-level abundance, not alignment |
| Compute, per sample | roughly 10–20 min on 4 threads for pseudoalignment |
| Full-cohort compute | roughly 25–40 CPU-hours |
| Storage after quantification | trivial (abundance tables only); raw can be deleted per sample as processed |

A streaming design — download one sample, quantify, delete, move on — would need only ~8 GB
of transient disk rather than 980 GB. **This is technically feasible and is the honest next
option, but it should not start until the clinical-endpoint problem is solved**, because
without PFS the analysis reduces to comparing score distributions across arms, which answers
nothing.

### What DECISION could and could not answer

**Could** (if expression and endpoints were both available): the distribution of the
differentiation state in established advanced refractory disease; a pure prognostic test in
the placebo arm; a treatment-by-score interaction for sorafenib benefit; expression state
versus DNA driver class; behaviour by histology including the 18 poorly differentiated cases.

**Could not**, under any circumstance: who becomes refractory in the first place, whether the
panel separates refractory from radioiodine-avid disease, or anything about radioiodine
benefit — **every patient in this trial was already refractory at entry.**

### Minimum author request

Written to `rai-response-genomics-atlas/docs/decision_minimum_data_request.md`. The essential
items are the normalised transcript matrix, the sample-ID-to-Supplementary-Table-1 mapping,
and per-patient PFS time and event. Without the third, the matrix alone is of limited use.

---

## Zheng 2024 metabolomics — manual step required

*Scientific Reports* 2024, PMC11079026. RAIR 20 versus non-RAIR 14, serum LC-MS, refractory
status assigned by the four ATA 2015 criteria.

The Europe PMC supplementary bundle was retrieved in full. Its single supplementary document
(699 KB) contains sample-preparation methods, QC-stability description, and one chemistry
table on tyrosine/MIT/DIT iodination — **no sample-level metabolite matrix and no per-patient
clinical table.**

One design feature is worth recording: **serum was collected before the first surgery**, so
refractory status was assigned afterwards. That makes it a genuinely prospective-in-design
discovery cohort, unlike most cross-sectional comparisons — but it also means the sampling
timepoint precedes the exposure by a long and variable interval, which must be audited before
any effect is quoted.

The raw matrix is stated to be in a Baidu repository with a password given in the paper.
That requires a browser and is logged as a manual action, not attempted automatically.

**Even if obtained, this cannot test the eight-gene panel** — it is serum metabolites, not
tumour transcripts. Its role is a metabolomic layer for the atlas and possible pathway
convergence with the Liu proteomics.

---

## ERRITI — endpoint reference, nothing more

*Clin Cancer Res* 2022, ccr-22-0437. Both supplementary documents were downloaded and read.
They contain the solid-tumour sequencing panel gene list, the full eligibility criteria, and
figure legends. **There is no per-patient data table**; the individual uptake, dosimetry and
thyroglobulin values live inside main-text Figure 2 as graphics.

Two details are worth keeping. The trial's mutation panel **includes SLC5A5 (exons 2, 3, 6, 11)
and all of TSHR** — one of the very few clinical panels that sequences the iodide machinery.
And the trial measured the full endpoint ladder that our project has been arguing for:
¹²³I scintigraphy uptake, ¹²⁴I PET dosimetry with absorbed dose, ¹³¹I therapy, thyroglobulin
trajectory, RECIST, and FDG-PET.

**Role: the reference definition of a radioiodine endpoint hierarchy**, against which the
vaguer endpoints in TCGA (`treatment_best_response`), GSE151179 (uptake at the metastatic
site) and GSE138042 (refractory versus sensitive) can be positioned. That comparison table is
the deliverable, not an analysis.

---

## Correction to an external lead

`10.1038/s41467-025-61788-w` was suggested as a DECISION re-analysis with a Source Data
workbook. It is not. It is a non-small-cell lung cancer study on prolonging response to EGFR
inhibition. **It must not be cited as a DECISION source, a thyroid clinical-label source, or
a processed-matrix source.** Logged in `audit/rai_integration_20260806/17_EXTERNAL_SOURCE_CORRECTIONS.tsv`.

---

## Redifferentiation trials — the one positive result of the round

Not on the original candidate list. It came from inverting the question the project had been
asking all day: instead of *who becomes refractory*, **who can be brought back**.

The ITOG 2025 consensus statement (*Lancet Diabetes Endocrinol*, PMC13011886) Table 3 reports
restored radioiodine uptake by driver genotype across published redifferentiation trials.
Extracted from the Europe PMC full-text XML and pooled on the logit scale
(DerSimonian–Laird, Haldane–Anscombe correction):

| Genotype | Trials | Restored / evaluable | Crude | Pooled random-effects | I² |
|---|---|---|---|---|---|
| BRAF V600E | 6 | 43 / 65 | 66.2% | **62.1% [45.9, 76.0]** | 36% |
| RAS | 3 | 33 / 40 | 82.5% | **78.7% [52.7, 92.5]** | 47% |

**RAS versus BRAF V600E: OR = 2.41, Fisher exact P = 0.077** on crude counts — numerically
higher in RAS, **not statistically significant**.

**Why it matters.** It bounds the driver-axis null rather than contradicting it. The pooled
driver-versus-refractoriness meta-analysis (n = 294, OR = 0.75, P = 0.53, I² = 63%) asked
whether drivers discriminate *who becomes* refractory; they do not. This asks whether drivers
relate to *whether refractoriness can be reversed*; they appear to. Different question,
different estimand, no contradiction — and it sharpens where the eight-gene expression axis
should be aimed: at the question drivers fail to answer, not the one they partly do.

**Limits, all of which must travel with the number.** Arm-level aggregate, not patient-level —
no adjustment, no bootstrap, no confounding assessment. Drugs are heterogeneous within
genotype (selumetinib, dabrafenib, trametinib, combinations, vemurafenib + anti-ErbB3). The
definition of "restored uptake" varies across trials — I-124 PET dosimetry, I-123 scan, or a
qualitative read — which is precisely the reason the endpoint hierarchy in
`03_ERRITI_ENDPOINT_AUDIT.md` was built. And publication bias plausibly favours successful
redifferentiation trials: NCT00085293 (decitabine, **0/12 restored**) has no publication we
could find.

Ledger rows REDIFF-01 through REDIFF-03 in `audit/rai_integration_20260806/02_EVIDENCE_LEDGER.tsv`.

## The single recommended next action

**Ask the DECISION authors for the normalised transcript matrix and per-patient PFS.**
Not the raw reads — those we can process ourselves if we must. The irreplaceable item is the
subject-level endpoint mapping, which was destroyed by spreadsheet formatting in the published
supplement and exists nowhere else. Jaume Capdevila (Vall d'Hebron) is named as the
biomaterial provider on all 125 BioSamples and is the natural contact.

Everything else this round is either endpoint reference material or requires a browser step.
