# 02 — Zheng 2024 serum metabolomics readiness assessment

**Grade D · MANUAL_DOWNLOAD_REQUIRED**, and low priority even if unblocked.

## Identification

Zheng W, Tang X, Dong J, Feng J, Chen M, Zhu X.
*Metabolomic screening of radioiodine refractory thyroid cancer patients and the underlying
chemical mechanism of iodine resistance.*
**Sci Rep** 2024, doi 10.1038/s41598-024-61067-6 · PMID 38719979 · PMC11079026.
Zhejiang Cancer Hospital, Hangzhou. Ethics IRB-2020-438 Ke.

Design: serum LC-MS metabolomics, **20 RAIR versus 14 non-RAIR**, refractory status assigned
by the four ATA 2015 criteria.

## What was obtained and what it contains

The Europe PMC supplementary bundle was retrieved in full (1.3 MB zip). Its single
supplementary document — `41598_2024_61067_MOESM1_ESM.docx`, 699 KB — contains
**sample-preparation methods, a QC-stability description, and one chemistry table** on
tyrosine/MIT/DIT iodination.

**It contains no sample-level metabolite matrix and no per-patient clinical table.**

The raw matrix is stated to sit in a **Baidu repository with a password given in the paper
text**. That requires a browser session; no automated bypass was attempted, and it is logged
as manual action 3 in `06_MANUAL_ACTIONS.md`.

## One design feature worth recording

**Serum was collected before the first surgery**, so refractory status was assigned
afterwards. That makes this genuinely prospective in design — a discovery cohort rather than
the cross-sectional refractory-versus-avid comparison that most of our sources are.

But the same feature is also its main threat: the sampling timepoint precedes the exposure by
a long and variable interval. Any effect quoted from this cohort would need that interval
audited, and neither the paper nor the supplement reports it.

## The ceiling on what it could ever do for us

**It cannot test the eight-gene panel.** These are serum metabolites, not tumour transcripts.
No amount of access changes that.

Its only possible roles are:

1. a **metabolic layer** for the atlas, alongside the Liu proteomics;
2. a check for **pathway convergence** — whether the metabolic signature points at the same
   MAPK-differentiation axis our transcriptional work implicates.

Both are supporting material, not evidence for a claim.

## Non-independence warning

Zheng 2024 and Wang 2024 (*The Oncologist*, PMID 38760956, tissue metabolomics, 24 RAIR / 18
non-RAIR) are from the **same institution, the same laboratory, and share the corresponding
author Zhu Xin and the author Tang Xi**. Both use the 2015 ATA criteria over the same period.
Patient overlap is not disclosed in either paper and should be assumed until an author says
otherwise.

**They must never be counted as two independent cohorts.** Full comparison in
`04_WANG_ACAC_SCREENING.md`.

If the Baidu matrix is pursued, ask Zhu Xin for the Wang tissue matrix in the same message —
it is one request, not two.

## Verdict

`ZHENG_2024 · HOLD · browser-gated; wrong modality to test the panel; non-independent from
WANG_2024. Lowest priority of the three manual actions.`
