# HLA ↔ Cancer Separation Rules

**Status:** Active boundary contract for the THCA manuscript program (Papers 1, 2, 3, 4).
**Owner:** Seungho Cook
**Last updated:** 2026-05-08
**Authority:** This document supersedes ad-hoc decisions in conversation. Violations must be repaired in code/manuscript or explicitly waived in writing here.

---

## 0. Scientific rationale

HLA allele frequencies and cancer outcomes are independently meaningful axes that **must not be casually joined**. Doing so produces three classes of error:

1. **Confounded ancestry signal.** HLA frequencies vary sharply by ancestry; cancer cohorts are ancestry-imbalanced. A naive HLA-vs-outcome test in a mixed cohort recovers ancestry, not biology.
2. **Tumor/germline conflation.** Tumor RNA-seq and tumor BAMs are not reliable substrates for germline HLA typing (LOH, expression bias, allele dropout). Reporting a "patient HLA allele" from tumor data is a methods error, even if mechanically possible.
3. **Causal-direction errors.** HLA susceptibility loci for autoimmunity (HT, GD, AITD) are not the same axis as HLA-mediated tumor immune escape. The genetic models, the cell types, and the inference targets are different. Mixing them inside one paper is a category error that reviewers and ethics committees will (correctly) reject.

The separation below is therefore **not a stylistic preference** — it is a methodological boundary. Audit scripts in `scripts/boundary_audit/` enforce it mechanically.

---

## 1. Paper 1 — Cancer / DM1 / DM2 / RAI lineage

### 1.1 Allowed in Paper 1

- Thyroid cancer **DM1 / DM2 expression** signatures (from TCGA-THCA, MSK Landa 2016, Korean K2 + Lee external validation cohort, GSE-series cancer datasets).
- **RAI-lineage gene** silencing/recovery analyses (TG, TPO, TSHR, NIS/SLC5A5, DUOX1/2, IYD, etc.).
- **Fusion enrichment** in DM1 / BRAF-RAS-negative compartment (RET, NTRK1/3, ALK, others).
- **BRAF / RAS-negative compartment** characterization.
- **Thyroid differentiation gene** expression (e.g., PAX8 module).
- **Promoter methylation** subtype analyses on tumor-only data.
- **External cancer validation** (Korean external n=865, MSK PDTC/ATC, GSE-series tumor cohorts).
- **Single-cell / spatial cancer context** for tumor cells, stromal cells, and TME (cell-type frequency, expression modules, niches).
- **Clinical framework for reflex fusion testing** and prospective epigenetic / RAI re-induction evaluation.
- **HLA-I / HLA-II as transcriptomic gene-expression module signatures** *only* when used as a residualization or contextual control with an explicit caption that says "Paper 1 residualization control only; deeper HLA analysis = Paper 2 territory." This is allowed because gene-expression module ≠ allele genotype. The audit allowlist permits this exact framing.

### 1.2 Forbidden in Paper 1

- **HLA allele association** (any DRB1/DQB1/DPB1/DQA1/DPA1 allele frequency claim).
- **MHC allele association.**
- **HLA imputation from cancer data** (TCGA-THCA, MSK, Korean K2/Lee, any cancer scRNA-seq, any cancer spatial dataset, any cancer WES/WGS BAMs).
- **HLA typing from cancer cohorts** even if labeled "exploratory."
- **Hashimoto / Graves / AITD HLA claims** of any kind.
- **Autoimmune susceptibility claim** anywhere in the Paper 1 storyline.
- **Carrier frequency or allele frequency comparisons** between cancer and population baselines.
- **HLA-based cancer risk prediction.**
- **HLA-based cancer prognosis** (HLA × OS / DSS / PFI / recurrence / RAI response, etc.).
- **HLA-based patient selection** (e.g., "HLA-X+ patients should receive Y").
- Importing Paper 2 / Paper 4 HLA results into Paper 1 narrative.

---

## 2. Paper 2 — Hashimoto / HT-overlap PTC HLA exploratory

(Per `v18_paper2_HT_isolated`: Paper 2 is HT-overlap PTC ONLY. Pillar I v2 = Korean PTC vs Korean baseline (AFND South Korea pool). GD/Chu 2018 forest is Paper 4 reserve.)

### 2.1 Allowed in Paper 2

- **Exploratory prioritization** of candidate HLA alleles from Korean AITD literature.
- **Hypothesis generation** for future validation.
- **Professor / collaborator discussion material.**
- **Validation-roadmap design** (sample-size calc, sequencing platform choice, consent template).
- **Korean AITD HLA literature curation** (Shin 2019, Cho 2011, Park 2005, Jang 2011, Baek 2021, Kwak 2014/KoGES context).
- **External HLA population baseline comparison** (AFND South Korea pool, AFND East Asian, AFND global, Korean 5,802 reference, IPD-IMGT/HLA nomenclature).
- **Allele frequency baseline comparison with caveats** (must distinguish allele frequency vs carrier frequency in every table).
- **Candidate HLA allele ranking for future validation.**

### 2.2 Forbidden in Paper 2

- **Validated cancer susceptibility claim** ("HLA-X causes thyroid cancer").
- **Causal thyroid cancer risk claim.**
- **Clinical risk prediction.**
- **Patient selection** (any HLA-based "give X to patient with allele Y").
- **Cancer survival association** (OS / DSS / PFI / recurrence / DSF).
- **RAI response association.**
- **DM1 / DM2 association.**
- **Fusion / BRAF / RAS / TERT association.**
- **Stage / lymph-node / distant-metastasis association.**
- **Pan-cancer outcome data to validate HLA-autoimmune biology.**
- Importing Paper 1 DM1 / RAI / survival claims into HLA interpretation.
- Treating HT/PTC bridge datasets as if they validate HLA-cancer association.

---

## 3. Bridge quarantine — HT + PTC expression / single-cell / spatial

### 3.1 Allowed in bridge

- **HT vs PTC vs HT+PTC expression comparison.**
- **Immune-context description** (B / T / plasma cell infiltration, TLS context, AICDA, IGHV clonality).
- **HLA-I / HLA-II *gene-expression* module signatures** as a transcriptomic immune-context proxy. Must be labeled "transcriptomic context only — not allele genotype" wherever shown.
- **Hypothesis generation only.** Bridge results frame open questions; they do not close them.

### 3.2 Forbidden in bridge

- **HLA allele association** (no allele-level genotype claim from bridge data, ever).
- **Cancer susceptibility claim** built on bridge alone.
- **Causality** language from bridge.
- **Clinical risk** language from bridge.
- **Prognosis / survival** language from bridge.
- **Patient selection** language from bridge.

---

## 4. Restricted / rejected

- Datasets behind **dbGaP / EGA controlled access** are catalogued but never downloaded without credentials.
- Datasets with **unclear phenotype labels** are flagged `NEEDS_MANUAL_REVIEW` and excluded from manifests until resolved.
- Datasets that **mix tumor and germline HLA without explicit consent / access** are rejected for current manuscript use.

---

## 5. Audit enforcement

Three audit scripts in `scripts/boundary_audit/` mechanically enforce these rules. They read `config/boundary_allowlist.yml` for permitted exceptions (boundary-warning sentences, separation-rule docs, audit reports).

- `audit_paper1_no_hla_contamination.py` — exits non-zero if HLA / Hashimoto / Graves / AITD terms appear in Paper 1 outside allowlisted contexts.
- `audit_paper2_no_cancer_claim_contamination.py` — flags cancer-outcome terms in Paper 2 for manual review.
- `audit_bridge_quarantine.py` — enforces "no allele-level HLA association tables, no cancer survival × HLA tables; immune-context summaries only."

Any modification to allowlists requires a written justification line in this file.

---

## 6. Operational reminders

1. **Never** use TCGA-THCA / MSK / Korean K2 BAMs for HLA imputation in service of any Paper 1 claim.
2. **Never** join HLA allele tables to cancer outcomes without an explicit Paper 4 sign-off.
3. **Always** distinguish allele frequency from carrier frequency in Paper 2 / reference manifests.
4. **Keep Korean external validation = n=865** unless Seungho explicitly reverts.
5. Paper 1 bioRxiv submission must complete before any Paper 2 / Paper 4 HLA-cancer crossover analysis is even discussed in writing.

---

## 7. Change log

- 2026-05-08 — initial scaffold written by Claude under Seungho's data-acquisition + boundary-hardening directive.
