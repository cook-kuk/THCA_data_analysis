# Siraj 2022 (KFSHRC, n = 158 PTC) — driver class vs radioiodine refractoriness

**Date** 2026-08-06 · **Script** `scripts/siraj2022_driver_vs_rai_refractory_2026_08_06.py`
**Source** Siraj AK et al., *Cancers* 2022;14(6):1584, doi 10.3390/cancers14061584.
Supplementary Tables S1 and S3 are open (CC BY) and carry per-patient rows. Local copy at
`/data/rai_atlas/external/siraj2022/SupplementaryTables.xlsx`.

This is the largest cohort we can use **today** without any data-access application:
158 papillary thyroid carcinomas with adjudicated radioiodine status
(**Refractory 66 / Avid 92**), cumulative radioiodine activity, thyroglobulin before and
6 months after treatment, progression-free survival, and 4,788 per-patient somatic
mutations. The matched whole-exome data (EGAS00001001788) still requires an ICGC DACO
application; the supplement does not.

## Cohort behaves exactly as it should — internal validity checks pass

| Check | Refractory | Avid | P |
|---|---|---|---|
| Cumulative radioiodine activity (mCi), median | **149.5** | 100.0 | 4.7 × 10⁻¹⁰ |
| Thyroglobulin 6 months after radioiodine (µg/L), median | **3.5** (n = 42) | 0.5 (n = 32) | 6.1 × 10⁻⁵ |
| Progression events | 53 / 66 | 1 / 92 | — |

Refractory patients received half again as much radioiodine and retained sevenfold higher
thyroglobulin. The labels are real and behave as a nuclear-medicine physician would expect.

## Primary result: driver class does not predict refractoriness

| Driver class | n | Refractory | % | OR vs rest | P |
|---|---|---|---|---|---|
| BRAF V600E | 81 | 32 | 40% | 0.83 | 0.63 |
| **BRAF/RAS-negative** | 56 | 26 | **46%** | 1.34 | 0.40 |
| RAS hotspot | 21 | 8 | 38% | 0.84 | 0.81 |

Overall χ² across the three classes: P = 0.67. The BRAF/RAS-negative compartment carries a
numerically higher refractory fraction, in the hypothesised direction, but the difference is
well inside noise at this sample size.

Progression-free survival by driver class: BRAF/RAS-negative HR = 1.67 (95% CI 0.97–2.86),
P = 0.063 — a borderline trend toward worse progression-free survival, not a result.

## What this means for the manuscript

It does **not** refute the framework, and it should not be reported as if it did. Our claim
concerns a transcriptional differentiation state *within* the BRAF/RAS-negative compartment,
not the compartment itself; this cohort has no transcriptome and therefore cannot score the
panel. What it does establish, in an independent advanced-disease cohort with genuine
refractoriness labels, is the negative premise the paper rests on:

> **Driver mutation class alone does not identify radioiodine-refractory disease
> (P = 0.67, n = 158).** That is precisely why a differentiation-axis readout is needed
> rather than a driver-axis one.

Used that way this is a supporting citation with our own re-analysis behind it, which is
stronger than citing the paper's conclusions.

## Do not use the progression-free survival contrast

Cox regression gives refractoriness HR = 116 (95% CI 16–844), P = 2.5 × 10⁻⁶. This number is
meaningless as an outcome finding: the cohort's refractoriness definition **includes
structural progression** among its seven criteria, so the exposure and the outcome share a
component. Reporting it would be circular, and a reviewer familiar with the Siraj criteria
would catch it immediately. It is retained in the output table only with an explicit
circularity flag.

## Files

- `results/tables/siraj2022_driver_vs_rai_2026_08_06.tsv`
- `results/tables/siraj2022_patient_level_2026_08_06.tsv` (158 rows, reusable)
- `results/figures/figure_siraj2022_driver_rai_2026_08_06.{png,pdf}`

Driver calls were made from Table S3 by requiring a non-synonymous coding change: BRAF
V600E; RAS = NRAS/HRAS/KRAS at codons 12, 13 or 61; everything else classed BRAF/RAS-negative.
TERT promoter status is not recoverable from the exonic annotation in Table S3 and was not used.
