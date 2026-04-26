# Task C — 4-cohort THCA LODO extension: label-extraction failure log

_v8.1 rigor upgrade. Attempt to extend THCA LODO from 2 cohorts
(TCGA-THCA + GSE27155) to 4 cohorts by incorporating GSE33630 and
GSE29265. Generated 2026-04-24._

## Bottom line

**Cannot extend.** Neither GSE33630 nor GSE29265 provides the joint
BRAF-and-RAS mutation annotation required by the v5.1 binary
BRAF-vs-RAS classification task. 2-fold LODO remains the correct
upper bound for THCA under honest labels; a 4-fold LODO would require
either (i) manual curation from the source publications or (ii)
expression-signature-based label inference (circular for DIAL).

## What's in the two candidate cohorts

Parsed from `GSE29265_metadata.tsv` and `GSE33630_metadata.tsv`
(characteristics_ch1 column):

| Cohort   | Total n | PTC  | ATC | Normal | BRAF annotation | RAS annotation |
|----------|--------:|-----:|----:|-------:|-----------------|----------------|
| GSE29265 |      49 |   20 |   9 |     20 | partial (PTC only, BRAF+/−/NA) | **none** |
| GSE33630 |     105 |   49 |  11 |     45 | **none**        | **none**       |

GSE29265 PTC label distribution: 8 BRAF+, 10 BRAF−, 2 NA. Even the
BRAF side has ≤ 20 samples; the RAS side is structurally absent from
the GEO phenotype record.

`characteristics_ch1` full-text grep (case-insensitive) for
`ras|kras|hras|nras` returned **0 matches in either cohort**.

## Why this is not a pipeline bug

v5.1's DIAL task is a supervised BRAF-vs-RAS binary classification.
LODO requires that every test fold carries both BRAF and RAS
positive samples. Adding GSE29265 (BRAF+/BRAF− but no RAS) or
GSE33630 (no mutation labels at all) to the training pool would
either:

1. Force invalid folds where the test cohort has only one class
   (GSE29265 as a fold would hold out 18 PTCs of which 0 are RAS+),
   or
2. Require inferring RAS status from expression signatures (for
   example, the TCGA-2014 BRS BRAF-RAS Score, or the Landa 2016
   "RAS-like" signature), which are themselves expression-based
   classifiers and would bias DIAL by introducing label-feature
   circularity.

Neither option yields a clean LODO extension.

## What *could* rescue a 4-fold LODO (documented for v8.2, not run)

1. **Manual curation from source publications.**
   - GSE33630: Tomás *et al.* (2012) *Endocrine-Related Cancer* —
     original paper may report per-sample BRAF status in supplementary.
   - GSE29265: Giordano *et al.* (2009/2011) — supplementary table
     may carry the RAS genotype.
   Estimated effort: 1–2 days of manual spreadsheet work per cohort.

2. **BRS-based expression inference** (circular for DIAL — use only
   for *validation* not for DIAL computation).
   - Define BRAF-like / RAS-like classes from pre-specified BRS genes
     from Landa *et al.* 2016 (non-overlapping with v5.1 biomarkers).
   - Validate on TCGA-THCA where both labels are known, then propagate
     to GSE33630/GSE29265.
   - Use the propagated labels for *downstream biomarker checks* only;
     do **not** feed into DIAL computation for those cohorts (it
     would circularly make DIAL = 0 for BRS genes).

3. **Orthogonal mutation datasets.**
   - dbGaP / cBioPortal may have the original Giordano cohort's
     mutation calls that can be cross-referenced to the GSE
     sample IDs via the patient identifier.

## What this means for the v5.1 paper

- LODO fold count stays at **2** (TCGA-THCA vs GSE27155).
- The specificity claim for THCA is therefore bounded by the fact
  that both cohorts share a papillary-dominant histology and that
  GSE27155's 13/41 RAS samples are the *only* non-TCGA RAS pool
  available in public BRAF-vs-RAS-annotated thyroid data.
- A meaningful 3-fold or 4-fold LODO refresh is a v8.2 task blocked
  on manual curation.

## Artefact

- `cohort_label_attempt.tsv` — per-sample extracted fields
  (tissue, braf, braf_delph, ras_status_available) across both cohorts
- `label_extraction_failure_log.md` (this file)

## References

- Tomás G *et al.* A general method to derive robust
  organ-specific gene expression-based differentiation indices:
  application to thyroid cancer diagnostic.
  *Endocrine-Related Cancer* 2012.
- Giordano TJ *et al.* Thyroid gene expression profiling.
  *Endocrine Pathology* 2009/2011.
- Landa I *et al.* Genomic and transcriptomic hallmarks of
  poorly-differentiated and anaplastic thyroid cancer.
  *Cell* 169:803 (2016). (BRS)
- Cancer Genome Atlas Research Network. Integrated genomic
  characterization of papillary thyroid carcinoma.
  *Cell* 159:676 (2014). (TCGA-THCA BRS score)
