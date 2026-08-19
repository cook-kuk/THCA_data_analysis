---
title: DM1 audit-locked quantitative results
date: 2026-07-30
rule: values are immutable until verified against primary outputs; interpretation and reference direction may remain unresolved
---

# Audit-locked results

These values came from the project status supplied by the author. They are a **starting ledger**, not a substitute for primary output verification.

## Discovery — TCGA-THCA

| Claim ID | Quantity | Locked value | Required verification before final prose |
|---|---|---:|---|
| D-01 | Discovery cohort n | 504 | exclusion flow and unique-patient count |
| D-02 | DM1 proportion | 28.4% | denominator and clustering assignment |
| D-03 | Pan-genome ARI | 0.92 | compared partitions and feature universe |
| D-04 | Driver-only ARI | -0.007 | adjusted Rand implementation and seed stability |
| D-05 | BRAF mRNA effect size | d=-0.044 | group coding and standardized-effect convention |
| D-06 | BRAF mRNA p-value | 0.57 | test type and two-sidedness |
| D-07 | Driver single-feature AUC | approximately 0.5 | individual values and uncertainty |

## Driver and dark-matter analyses

| Claim ID | Quantity | Locked value | Required verification |
|---|---|---:|---|
| DM-01 | BRAF/RAS-negative cases rescued | 131/180 (73%) | definition of “rescue” and denominator |
| DM-02 | BRAF-positive DM1 proportion | 0.7% | exact BRAF definition and cluster coding |
| DM-03 | RAS-mutant DM1 proportion | 96.4% | mutation list and cluster coding |
| DM-04 | BRAF-/RAS- DM1 proportion | 49.1% | missing-driver handling |
| DM-05 | BRAF-/RAS- DM2 proportion | 50.9% | missing-driver handling |

## Survival

| Claim ID | Analysis | Locked value | Critical semantic checks |
|---|---|---:|---|
| S-01 | TCGA OS | HR 2.30; 95% CI 0.77–6.88 | event coding, reference group, covariates, time origin |
| S-02 | MSK OS | HR 2.67; 95% CI 1.17–6.10 | advanced-disease enrichment, reference group, covariates |
| S-03 | Pooled OS | HR 2.53; 95% CI 1.31–4.89; I²=0% | fixed/random method, log-HR SE, heterogeneity calculation |
| S-04 | BRAF-positive PFI | n=287; HR 0.66; p=0.013 | **must verify which state is numerator/reference** |
| S-05 | DM1 × BRAF interaction | interaction HR 0.47; p=0.022 | formula, coding, covariates, scale, multiplicity |
| S-06 | BRAF-negative PFI | HR 1.19; p=0.52 | reference group and power |
| S-07 | RAS-positive PFI | HR 1.04; p=0.92 | reference group and power |

Do not describe S-04 as “DM1 is protective” or “DM1 is adverse” until source code and model summary establish direction. Do not call S-05 a treatment-predictive interaction.

## Structural variants and kinase fusions

| Claim ID | Quantity | Locked value | Required verification |
|---|---|---:|---|
| FUS-01 | DM1 fusion positivity | 63/82 (76.8%) | available-SV denominator and fusion inclusion list |
| FUS-02 | DM1 versus DM2 fusion OR | 7.41; 95% CI 4.38–12.55; p=1.9e-13 | contingency table and test |
| FUS-03 | RET-fusion cases assigned DM1 | 27/33 (81.8%) | RET fusion definition |
| FUS-04 | SV missingness test | chi-square p=0.56 | tested variables and MAR wording |

A non-significant missingness association does not prove MAR. Prefer “we found no evidence that availability differed by the tested group variable” rather than “MAR confirmed.”

## Methylation

| Claim ID | Quantity | Locked value | Required verification |
|---|---|---:|---|
| M-01 | Paired HM450 cohort n | 503 | patient matching and probe QC |
| M-02 | TPO Cohen's d | 2.30 | direction and beta aggregation |
| M-03 | TPO p-value | 1.9e-18 | test and multiplicity |
| M-04 | DIO1 Cohen's d | 1.24 | direction |
| M-05 | TSHR Cohen's d | 1.20 | direction |
| M-06 | PAX8 Cohen's d | 0.97 | direction |
| M-07 | TG Cohen's d | 0.86 | direction |
| M-08 | FOXE1 Cohen's d | 0.84 | direction |
| M-09 | NKX2-1 Cohen's d | 0.63 | direction |
| M-10 | SLC5A5 Cohen's d | 0.22; p=0.42 | promoter definition and probe coverage |
| M-11 | Mean eight-gene beta | DM1 0.385 vs DM2 0.253; +52% | formula for percentage change and paired/unpaired test |
| M-12 | Driver-class beta | BRAF 0.37; RET 0.39; RAS 0.27 | exact n and uncertainty by class |

Use “hypermethylation correlated with lineage loss” or “epigenetic correlate.” Do not write “methylation silences the state” unless longitudinal or perturbational evidence exists.

## External concordance and validation

| Claim ID | Quantity | Locked value | Required verification |
|---|---|---:|---|
| V-01 | Master forest entries | 14 | list each dataset and independent unit |
| V-02 | Mean Cohen's d | 2.81 | weighting convention |
| V-03 | Median Cohen's d | 2.37 | entry list |
| V-04 | Minimum d | all at least 1.56 | sign convention |
| V-05 | Direction consistency | 80/80 cells | define cell: gene × cohort or other unit |
| V-06 | Lee 2024 Korean PTC | n=632; d=5.93 | platform, FFPE status, leakage risk |
| V-07 | Lee replication A2-R | n=370; Mann–Whitney p=2.7e-7 | independence from discovery and group definitions |
| V-08 | GPL570 cohorts | four cohorts, total n=205; all rho ≤ -0.84 | variables correlated and sign convention |
| V-09 | Mun 2025 proteomics | n=336; 7/7 proteins sign-consistent | missing eighth protein and statistical support |
| V-10 | Pu 2021 single-cell | per-patient r=0.798–0.886; Bonferroni p<1e-10 | patient n, pseudobulk definition, non-independence |
| V-11 | Lu 2023 single-cell | thyrocyte-gated gradient retained | quantitative effect and patient-level analysis |
| V-12 | Landa 2016 overlap | 5/8 genes: TG, TSHR, TPO, PAX8, DIO1 | exact comparison and direction |
| V-13 | FFPE vs fresh-frozen | KS p=0.44 | compared distributions and interpretation |

Avoid treating “no significant difference” as equivalence unless an equivalence test was performed.

## RAI-refractory alignment

| Claim ID | Analysis | Locked value | Required verification |
|---|---|---:|---|
| RAI-01 | GSE151179 alignment | Cohen's d approximately -1.0 | exact sign and state-score definition |
| RAI-02 | Group difference | Mann–Whitney p approximately 1e-4 | exact p-value, groups, n |

This supports transcriptional concordance with post-RAI-refractory tumors. It does not establish pretreatment prediction of uptake or response.

## Three-marker in-silico proxy

| Claim ID | Analysis | Locked value | Required verification |
|---|---|---:|---|
| IHC-01 | TG + PAX8 + NKX2-1 in BRAF-positive PFI | log-rank p=1.7e-4 | thresholding and leakage |
| IHC-02 | TG + PAX8 + NKX2-1 in BRAF-positive PFI | Cox HR 0.65; 95% CI 0.47–0.90 | reference direction and covariates |
| IHC-03 | TSO500 three-gene proxy | p=0.30 | why included and platform limitation |
| IHC-04 | Power simulation | n=200; >90% power | effect-size source, prevalence, event rate, alpha, simulation code |

Call this an “in-silico three-marker expression proxy motivated by clinically available antibodies.” Do not call it a validated IHC panel.
