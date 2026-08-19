# 01 — DECISION trial (PRJNA563018) readiness assessment

**Grade B · METADATA_READY_RAW_TOO_LARGE.**
Half of what we needed was recovered. The other half was destroyed in publication.

## What DECISION is

The phase 3 randomised trial of sorafenib versus placebo in progressive, locally advanced or
metastatic **radioiodine-refractory** differentiated thyroid carcinoma. The RNA-seq substudy
of 125 patients was deposited by Vall d'Hebron Institute of Oncology on 2019-08-30.

Paper: Mol Cancer Ther 2020;19(1):312, **PMID 31540966 — paywalled, no PMCID**. The
supplement was obtained through the AACR figshare mirror (article 22506703), the same route
that worked for Boucai 2023 and ERRITI.

## Why it was worth chasing

**TCGA gave us no untreated comparator.** Among TCGA-THCA patients who did not receive
radioiodine there were **2 events**, which makes the treatment-by-score interaction —
the only thing that could turn a prognostic association into a predictive one —
formally unestimable. That is the single structural weakness in our TCGA result
(OR = 0.221, P = 0.019 for new tumour event).

DECISION has **57 placebo-arm patients**. That is exactly the missing comparator, and it is
why this dataset was graded and pursued ahead of everything else this round.

## Recovered, per sample, all 125

From BioSample attributes:

| Field | Distribution |
|---|---|
| **Treatment arm** | **sorafenib (BAY 43-9006 400 MG BID) 68 · placebo 57** — reproduces the published split exactly |
| Histology | papillary 73 · follicular 33 · **poorly differentiated 18** · non-diagnostic 1 |
| Sex | male 65 · female 60 |
| Age | median 63 (range 32–87) |
| Library identifiers | `isolate` = AB####, `GENOMIC.LAB.SAMPLE.ID` = B15/###, library = L15/#### |
| Biomaterial provider | Jaume Capdevila, Vall d'Hebron, Barcelona |

From the supplement, **aggregate only**: ECOG 0 in 64%; BRAF mutated 33 (26.4%);
RAS mutated 20 (16%); expression profiles **BRAF-like 56 · RAS-like 28 · NoBRaL 41**;
substudy survival sorafenib median PFS 10.8 versus placebo 3.7 months,
HR 0.42 (95% CI 0.26–0.67), P < 0.001.

## Blocker 1 — the patient key was destroyed by spreadsheet formatting

Supplementary Table 1 is **five pages of raster images at 87–93 ppi**. Its columns are QC
(GOOD/BAD), CNAG library barcode, genomic lab sample ID, library ID, clinical trial number,
subject number and treatment group. Every column except QC is already in BioSample, so the
table adds almost nothing — **except the subject number, which is the key to the trial's
clinical database.**

That column renders as `1E+08`, `1,2E+08`, `1,4E+08`. The nine-digit trial subject IDs were
converted to scientific notation by a spreadsheet before the figure was rendered, and the
information is gone from the published file.

**Consequence: per-patient PFS time, PFS event, RECIST response and expression-cluster
assignment cannot be reconstructed from public files.** This is not a paywall problem or a
format problem. The data no longer exists in the published artefact.

## Blocker 2 — no processed expression matrix exists anywhere

Searched and confirmed absent:

| Repository | Result |
|---|---|
| GEO | zero hits for the BioProject and for the study terms |
| ArrayExpress / BioStudies | absent |
| Zenodo | absent |
| Figshare | only the single supplement PDF |
| GitHub | zero repositories |
| BioProject linked resources | raw reads only |

**The only deposited data is raw reads: 0.98 TB across 125 runs** (median ~66 M reads,
~10 Gbase per sample).

## Cost of processing it ourselves — estimated, not executed

| Quantity | Estimate |
|---|---|
| Per-sample FASTQ | ~7.8 GB (paired) |
| 3-sample pilot | ~23 GB |
| Full 125-sample download | **~980 GB** |
| Quantification | pseudoalignment (kallisto/salmon) against gencode — the eight panel genes need only transcript-level abundance |
| Compute per sample | ~10–20 min on 4 threads |
| Full cohort | ~25–40 CPU-hours |
| Transient disk with streaming (download → quantify → delete) | **~8 GB, not 980 GB** |

This is technically feasible on our infrastructure. It is nonetheless **not started**, and
should not be, until blocker 1 is resolved.

## What DECISION could and could not answer

**Could**, if both matrix and endpoints were available: the distribution of the
differentiation state in established advanced refractory disease; a pure prognostic test
within the placebo arm; a **treatment-by-score interaction for sorafenib benefit**;
expression state versus DNA driver class; behaviour in the 18 poorly differentiated cases.

**Could not**, under any circumstance: who becomes refractory in the first place; whether
the panel separates refractory from radioiodine-avid disease; anything about radioiodine
benefit. **Every patient in this trial was already refractory at entry.** Note the estimand
carefully — the untreated comparator here is untreated *with sorafenib*, not untreated with
radioiodine. It repairs the sorafenib interaction, not the radioiodine one.

## Minimum author request

Written to `rai-response-genomics-atlas/docs/decision_minimum_data_request.md`, not sent.
Three items, in order of irreplaceability:

1. **Per-patient PFS time and event** — exists nowhere else, cannot be reconstructed.
2. **The sample-ID-to-Supplementary-Table-1 mapping** — the key destroyed by formatting.
3. The normalised transcript matrix — useful, but we could generate our own from raw reads.

Without item 1, item 3 alone reduces the analysis to comparing score distributions between
arms, which answers nothing.
