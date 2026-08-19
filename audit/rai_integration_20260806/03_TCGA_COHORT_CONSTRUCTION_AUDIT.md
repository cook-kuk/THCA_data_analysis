# AUDIT 03 — TCGA radioiodine cohort construction

Generated 2026-08-06 by `rai-response-genomics-atlas/scripts/audit/01_tcga_cohort_and_index_course.py`

## Provenance of `treatment_best_response`

The field sits in the radiation table keyed by `bcr_radiation_barcode`, its human label is `measure_of_response` (CDE 2857291), and every mCi course carries a start day. Decisively, the value **varies within a patient across courses** and tracks course timing and treatment site:

| patient | course | dose (mCi) | site | response |
|---|---|---|---|---|
| TCGA-FK-A3S3 | day 162 | 106.4 | primary | Stable Disease |
| TCGA-FK-A3S3 | day 581 | 197.8 | primary | Complete Response |
| TCGA-FE-A3PA | day 91 | 217 | primary | Partial Response |
| TCGA-FE-A3PA | day 378 | 217 | distant | Complete Response |

A patient-level best-response field would be constant across a patient's rows. It is not. The field is therefore course-specific.

**Consequence for the first-pass analysis.** That analysis aggregated to the *best* response across courses, which discards the earlier, worse response (Stable Disease for TCGA-FK-A3S3) and biases the cohort toward complete response. For the question 'does the score relate to the initial response to radioiodine', the **first course** is the correct index.

## Cohort ledger

| quantity | value |
|---|---|
| radiation course rows (all units) | 293 |
| unique patients with any radiation record | 274 |
| mCi courses (= radioiodine) | 249 |
| Gy / cGy courses (= external beam) | 17 |
| courses with unit not recorded | 27 |
| unique patients with >=1 mCi course | 237 |
| patients with >1 mCi course | 12 |
| max mCi courses per patient | 2 |
| mCi patients who also received Gy/cGy | 6 |
| mCi courses with an evaluable response | 176 |
| mCi courses with start day recorded | 249 |
| patients whose courses DISAGREE on response | 2 |
| patients in R17 panel table | 500 |
| mCi patients matched to a panel score | 235 |

## Index-course sensitivity

| index_rule                        |   n_total |   n_nonCR |   n_CR |   cohens_d |   ci_low |   ci_high |    p_mwu |   spearman_rho_ordinal |   p_ordinal |
|:----------------------------------|----------:|----------:|-------:|-----------:|---------:|----------:|---------:|-----------------------:|------------:|
| first course (earliest start day) |       167 |        26 |    141 |      -0    |   -0.445 |     0.457 | 0.91038  |                -0.0193 |    0.804426 |
| last course                       |       167 |        24 |    143 |      -0.03 |   -0.505 |     0.459 | 0.730518 |                -0.0358 |    0.645756 |
| highest-dose course               |       167 |        24 |    143 |      -0.03 |   -0.505 |     0.459 | 0.730518 |                -0.0358 |    0.645756 |
| best response across courses      |       167 |        24 |    143 |      -0.03 |   -0.505 |     0.459 | 0.730518 |                -0.0358 |    0.645756 |
| worst response across courses     |       167 |        26 |    141 |      -0    |   -0.445 |     0.457 | 0.91038  |                -0.0193 |    0.804426 |

## Statistical framing correction

The first pass wrote that the cohort 'detects d >= 0.63 at 80% power, so clinically meaningful effects are excluded'. That conflates a minimum-detectable-effect calculation with an equivalence test. No equivalence margin was pre-specified. The defensible statement is that the confidence interval excludes effects larger than its own bounds, and that no association was detected — not that no effect exists.
