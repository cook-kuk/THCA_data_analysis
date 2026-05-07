# Paper 2 HLA — claim boundary v2 (per-allele honesty audit)
**Date:** 2026-05-06 (v2)
**Status:** exploratory candidate-prioritization map; *not* a case-control association.

This document defines, per allele, exactly what we are allowed to say after the strengthened metric-sensitivity audit, and exactly what is forbidden.

---

## v2 core statement (read first)

> **The S0 carrier-vs-allele forest initially showed 6/6 significant signals, but strict S3 allele-vs-allele harmonization retains only two depletion candidates: C*01:02 and DQB1*02:01. No enrichment candidate survives strict metric harmonization.**

DQB1*02:01 is further qualified by the 2026-05-06 zero-cell QC audit (`2026_05_06_paper2_dqb1_zero_cell_qc_report.md`): the depletion is most likely a typing-resolution artifact, and the claim is suspended pending orthogonal NGS-typer validation.

## v2 claim grade table (per-allele)

| Allele | Claim grade (v2) | S0 primary q | S3 direct q | S0 direction | S3 direction | Next validation priority |
|---|---|---|---|---|---|---|
| A*02:07 | **S0_only_exploratory** | 1.1×10⁻⁵ | 0.43 (n.s.) | enrichment | enrichment | medium |
| B*46:01 | **S0_only_exploratory** | 1.1×10⁻⁵ | 0.80 (n.s.) | enrichment | enrichment | medium |
| C*01:02 | **survives_strict_metric_but_direction_flips** | 4.3×10⁻⁴ | **7.4×10⁻⁴ ✓** | enrichment | **depletion** | **high** |
| DPB1*05:01 | **source_sensitive_artifact_likely** | 9.9×10⁻¹⁴ | 0.80 (n.s.) | enrichment | depletion | medium |
| DQB1*02:01 | **survives_strict_metric_but_zero_cell_caution** | 2.5×10⁻⁶ | **5.3×10⁻⁸ ✓** | depletion | depletion (zero-cell, suspended by QC audit) | **high** |
| DRB1*07:01 | **S0_only_exploratory** | 4.3×10⁻⁴ | 0.43 (n.s.) | enrichment | depletion | medium |

Reference TSV: `project/results/p2_pillar1_forest_v2/paper2_hla_claim_grade_v2.tsv`.

---

## 0. Universal boundary (all 6 alleles)

| Allowed | Forbidden |
|---|---|
| "Exploratory carrier-vs-allele forest with metric mismatch noted." | "Korean PTC HLA association." |
| "Direction-stability across S0/S1/S2/S3 metric scenarios reported." | "Causal HLA susceptibility." |
| "BH-q within six candidates reported per scenario." | "Genome-wide multiple-test-corrected association." |
| "Subcohort heterogeneity flagged where χ² p < 0.05." | "Population-level allele-frequency association." |
| "Matched-control NGS validation needed before any association claim." | "Clinical risk prediction / patient selection." |
| "Direct allele-vs-allele Fisher (S3) shows ____ direction with q = ____." | "Carrier frequency and allele frequency are interchangeable." |

---

## 1. Per-allele claim boundary

Format per row: **observation** → *allowed claim* → forbidden claim.

### 1.1 A*02:07
- **Observation:** S0 OR = 2.43 (q = 1.1e-5); S3 direct OR = 1.26 (q = 0.43, n.s.). PTC AF = 4.40%, baseline AF = 3.51%. Subcohort χ² p = 0.016 (Lee2024 dominates 66% of carriers).
- **Allowed:** "A*02:07 carrier frequency in Korean PTC is modestly higher than baseline allele frequency, but the apparent enrichment does not survive strict allele-vs-allele harmonization (S3 q = 0.43). Subcohort heterogeneity is present."
- **Forbidden:** "A*02:07 is enriched in Korean PTC." / "A*02:07 is a Korean PTC risk allele." / Any clinical or causal language.

### 1.2 B*46:01
- **Observation:** S0 OR = 2.16 (q = 1.1e-5); S3 direct OR = 1.05 (q = 0.80, n.s.). PTC AF = 5.32%, baseline AF = 5.06%.
- **Allowed:** "B*46:01 carrier-vs-allele forest signal does not survive strict allele-vs-allele harmonization (S3 q = 0.80, OR ≈ 1.05)."
- **Forbidden:** Any enrichment / association claim.

### 1.3 C*01:02
- **Observation:** S0 OR = 1.47 (q = 4.3e-4); S1 OR = 0.66; S2 OR = 0.69; **S3 direct OR = 0.68 (q = 7.4e-4, depletion).** PTC AF = 12.89%, baseline AF = 17.81%. Subcohort χ² p = 0.006 (heterogeneity flagged).
- **Allowed:** "C*01:02 shows direction-flip from S0 enrichment to S3 depletion; the depletion direction survives BH-q < 0.05 within six candidates. This is the strongest *exploratory* candidate to validate, and the validation hypothesis should be framed as **C*01:02 PTC depletion**, not enrichment. Subcohort heterogeneity (χ² p = 0.006) requires attention."
- **Forbidden:** "C*01:02 enrichment in Korean PTC" *(this was the metric-mismatch artifact)* / "C*01:02 is a Korean PTC risk allele" / Any causal claim. The depletion direction is also exploratory until matched-control NGS validation.

### 1.4 DPB1*05:01
- **Observation:** S0 OR = 1.96 (q = 9.9e-14); S3 direct OR = 0.98 (q = 0.80, n.s.). PTC AF = 36.20%, baseline AF = 36.67% (essentially equal). Baseline source = AFND South Korea pool, *not* IN2015. Jung 2023 reports DPB1*05:01:01 = 35.1% (consistent with baseline).
- **Allowed:** "DPB1*05:01 1차 forest signal of OR ≈ 2 reduces to OR ≈ 1.0 under strict allele-vs-allele comparison. The original signal is consistent with a metric-mismatch artifact compounded by source heterogeneity in the AFND pool baseline."
- **Forbidden:** "DPB1*05:01 enrichment in Korean PTC" / "DPB1*05:01 is associated with Korean PTC" / Any HLA-DP causal claim.

### 1.5 DQB1*02:01
- **Observation:** PTC carriers = 0 / 874 (callable 631). Baseline AF = 2.12%. S0 Haldane OR = 0.026 (q = 2.5e-6); **S3 direct OR = 0.018 (q = 5.3e-8, depletion).** Direction stable across all 4 scenarios (zero observation persists).
- **Allowed:** "DQB1*02:01 shows zero PTC carriers across the n = 874 pool callable for DQB1 (n = 631). The depletion signal is direction-stable across four metric scenarios and survives BH-q < 0.05. **However**, zero-observation may reflect (a) true biological depletion, (b) typing pipeline failure to call this specific allele, or (c) resolution-collapse error. Single-cohort zero observation cannot stand-alone."
- **Forbidden:** "DQB1*02:01 protective in Korean PTC" / "DQB1*02:01 absence is a Korean PTC risk factor" / Any causal interpretation. Independent matched-control NGS validation that reproduces the zero-or-near-zero rate is the minimum gate.

### 1.6 DRB1*07:01
- **Observation:** S0 OR = 1.73 (q = 4.3e-4); S3 direct OR = 0.84 (q = 0.43, n.s.). PTC AF = 5.92%, baseline AF = 6.93%.
- **Allowed:** "DRB1*07:01 carrier-vs-allele forest signal does not survive strict allele-vs-allele harmonization (S3 OR = 0.84, q = 0.43)."
- **Forbidden:** "DRB1*07:01 enrichment in Korean PTC" / Any association/causal claim.

---

## 2. Aggregate claim grade — six-allele forest

| Strict-metric outcome | Count |
|---|---|
| Survives BH-q < 0.05 in S3 (depletion direction) | **2** (C*01:02, DQB1*02:01 with zero-cell caveat) |
| Survives BH-q < 0.05 in S3 (enrichment direction) | **0** |
| S0-only signal (does not survive S3) | 4 (A*02:07, B*46:01, DPB1*05:01, DRB1*07:01) |

**Allowed aggregate claim:** "After strict metric harmonization, the six-allele HLA forest yields an exploratory candidate-prioritization map: zero alleles support an enrichment claim, two alleles show depletion direction surviving within-six BH-q < 0.05, and matched-control NGS validation is required for any association statement."

**Forbidden aggregate claim:** "Six HLA alleles are differentially distributed in Korean PTC." / "We identify Korean PTC HLA risk alleles." / "HLA susceptibility profile of Korean PTC."

---

## 3. Single sentence at maximum permitted strength

> "The six-allele HLA forest provides an exploratory candidate-prioritization map suggesting Korean PTC may have immune-genetic background signals requiring matched-control validation."

This is the only association-adjacent sentence permitted at the current evidence grade. Any sentence stronger than this is over-claim.

---

## 4. What it would take to graduate from this boundary

See `paper2_hla_next_validation_plan.md` for the full validation roadmap. Minimum gate to graduate: matched-control NGS HLA typing (independent Korean healthy population, same pipeline as the PTC pool), allele-level dosage harmonization on both sides, ≥ 2 independent PTC cohorts, pre-registered hypothesis and analysis plan.

---

**End of claim-boundary document. The boundary is binding for all Paper 2 HLA-related communication (manuscript drafts, slides, web pages, advisor briefings) until the validation gate is passed.**
