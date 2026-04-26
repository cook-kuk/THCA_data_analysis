# Cover Letter — Bioinformatics (OUP)

2026-04-25

Dear Editor,

We are pleased to submit our manuscript, *"DIAL: a direction-invariant audit metric reveals thyroid-specific BRAF/RAS classifier failure under linear batch correction"*, for consideration as an **Original Paper** in *Bioinformatics*.

## Motivation

Public cancer expression cohorts (TCGA combined with GEO) are routinely harmonised via empirical-Bayes batch correction (ComBat) before training molecular subtype classifiers. When the label and cohort directions share a linear subspace, however, covariate-preserving ComBat can silently invert a classifier's decision axis — a failure mode that is invisible to standard cross-validated AUC. The resulting classifier looks perfect in-cohort and generalises backwards under held-out-cohort evaluation. Given the prevalence of multi-cohort reanalysis in cancer genomics, and the growing use of public-data-trained classifiers in translational settings, we believe this failure mode deserves a named diagnostic and a body of specificity evidence. This is the subject of our paper.

## Key contributions

1. **A new metric, DIAL** (Direction-Invariant AUC Leakage), that detects label-axis inversion under linear batch correction, together with a formal label-flip lemma characterising the conditions under which it fires.
2. **Specificity evidence across five cancer types**: DIAL flags thyroid BRAF/RAS under LODO in 4 of 5 classifier families while returning zero across 20 of 20 classifier–cancer combinations in SKCM, LGG, LUAD, and COAD. Random-effects meta-analysis gives *m*(THCA) = 1.00 with *I²* > 97% under four of five classifiers.
3. **A mechanistic explanation**: the flip is driven by cohort–subtype non-orthogonality in the harmonised THCA matrix (within-subtype PC1 KS statistic = 1.00), and it is eliminated by FastRNA-style per-cohort centering — a constructive mitigation rather than a generic alarm.
4. **Cross-paradigm robustness**: the THCA flip reproduces under a variational quantum classifier and an RBF-kernel SVC, ruling out a linear-model-specific artefact. Gene-level DIAL fires; pathway-level (MSigDB Hallmark) DIAL is zero, identifying DIAL as a gene-level early-warning that complements pathway-level analyses.

## Significance for the Bioinformatics readership

The paper delivers a practical, testable diagnostic for a failure mode that arises whenever a multi-cohort expression study combines two platforms with unequal class priors. The fix — a BUHMBOX-style concept check and FastRNA centering — is cheap, reproducible, and directly implementable by any group currently running ComBat on a TCGA + GEO pairing. We expect the metric and the specificity framework to be useful outside thyroid carcinoma to any subtype classification task where label and cohort are near-collinear.

## Recommended reviewers

- **Jeffrey T. Leek** (Johns Hopkins; Fred Hutchinson Cancer Center) — SVA/ComBat lineage, batch-effects methodology.
- **W. Evan Johnson** (Rutgers) — author of the ComBat and ComBat-seq estimators.
- **Buhm Han** (Seoul National University) — METASOFT, BUHMBOX, FastRNA; statistical genetics.
- **Christina Kendziorski** (University of Wisconsin–Madison) — single-cell and bulk batch correction.
- **Wei Sun** (University of North Carolina / Fred Hutchinson) — expression batch effect statistics.

## Administrative

- **Competing interests**: none declared.
- **Funding**: this work received no external funding.
- **Data availability**: all inputs are public (TCGA MC3 via GDC and GEO accessions GSE27155, GSE22153, GSE16011, GSE31210, GSE39582); per-cohort harmonised matrices, DIAL tables, figures, and a Dockerfile reproducing the pipeline end-to-end are available at <http://40.82.129.113:8012/>. The interactive dashboard hosts every figure at gene and pathway level.
- **Ethics**: no primary patient data were collected; all cohorts are de-identified public repositories.
- **ORCID IDs**: to be finalised at submission.
- **Corresponding author**: to be confirmed (user email: `kukshomr@gmail.com`).

Thank you for your consideration. We look forward to a constructive review.

Sincerely,

Seungho Kuk (corresponding author)
Independent research, Seoul, Republic of Korea
<kukshomr@gmail.com>
