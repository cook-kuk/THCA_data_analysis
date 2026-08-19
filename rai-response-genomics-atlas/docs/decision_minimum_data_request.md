# Minimum data request — DECISION trial RNA-seq (PRJNA563018)

**Drafted 2026-08-06. NOT SENT.**

Paper: Identification of Expression Profiles Defining Distinct Prognostic Subsets of
Radioactive-Iodine Refractory Differentiated Thyroid Cancer from the DECISION Trial.
*Mol Cancer Ther* 2020;19(1):312, PMID 31540966. Submitter of PRJNA563018: Vall d'Hebron
Institute of Oncology. Biomaterial provider named on all 125 BioSamples:
**Dr Jaume Capdevila**, Medical Oncology, Vall d'Hebron University Hospital, Barcelona.

## Why a request is needed at all

We already recovered from public sources, per sample: treatment arm (sorafenib 68 / placebo
57), histology (papillary 73, follicular 33, poorly differentiated 18), age, sex, library
barcodes, and QC status. Raw reads are public (0.98 TB) and we can quantify them ourselves.

The blocker is different. **Supplementary Table 1 as published renders the trial Subject
Number in scientific notation** — every value appears as `1E+08`, `1,2E+08` or `1,4E+08`,
so the nine-digit subject identifiers were lost to spreadsheet formatting before publication.
That column is the only bridge between the sequenced samples and the trial's clinical
endpoints, and it exists nowhere else in public.

Without it, even after processing the full 0.98 TB we could compare score distributions
between arms but could not run a single survival model.

## What we would ask for, in priority order

| # | Item | Why | Substitutable? |
|---|---|---|---|
| 1 | **Per-patient PFS time and event** keyed to the sample identifiers already in BioSample (`isolate` AB####, or `GENOMIC.LAB.SAMPLE.ID` B15/###) | The only irreplaceable item. Enables the placebo-arm prognostic test and the treatment-by-score interaction. | No |
| 2 | Sample-ID to Supplementary-Table-1 mapping with an intact subject number | Alternative route to item 1 if the trial database can be joined | No |
| 3 | Normalised transcript-level expression matrix (the ~11,106-transcript analysis matrix) | Saves ~980 GB of download and ~30 CPU-hours; we can generate an equivalent ourselves if refused | Yes, by processing raw reads |
| 4 | Per-patient BRAF and RAS status | Published only in aggregate (33 and 20) | Partly, from raw reads |
| 5 | Per-patient BRAF-like / RAS-like / NoBRaL assignment | Lets us compare our axis against theirs in the same patients | Partly, recomputable |
| 6 | RECIST best response | Secondary endpoint | No |

Items 3–5 are recoverable by our own processing. **Items 1, 2 and 6 are not.**

## The question we would state

Whether a thyroid-lineage differentiation score measured in archival tumour is prognostic in
established radioiodine-refractory disease, tested first in the **placebo arm alone** so that
the estimate is uncontaminated by treatment effect, and second as a treatment-by-score
interaction for sorafenib benefit.

We would state plainly what this cannot address: every DECISION patient was refractory at
entry, so nothing in this cohort speaks to who becomes refractory, to separating refractory
from avid disease, or to radioiodine benefit.

## What we would offer

Co-authorship on any resulting work; our analysis code and the TCGA and Zhang results in
advance so the investigators can judge whether the question is worth their time; and analysis
restricted to the stated hypothesis with no attempt at re-identification.

## Before sending

Substitute the audited effect sizes — the structural-outcome odds ratios are Firth-penalised
(0.221 and 0.224), and the initial-response figure is the first-course rule (d = 0.00,
95% CI −0.45 to +0.46, P = 0.91). Do not send a version quoting the earlier maximum-likelihood
values.
