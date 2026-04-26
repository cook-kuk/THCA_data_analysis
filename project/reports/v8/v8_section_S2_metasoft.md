> ⚠️ **v5.2 SUPERSEDED NOTICE (2026-04-25).** This robustness section was designed to stress-test the v5.1 THCA DIAL flip (DIAL=0.494, batch_entangled in 4/5). Under proper LODO ComBat the underlying flip disappears (v5.2: DIAL=0.000 / true_biology in 5/5; see reports/v5p2/v5p2_critical_assessment.md). Conclusions below are conditional on the v5.1 result they reference.

# S2. Random-Effects Meta-Analysis of Per-Cancer DIAL (v8 Task 2)

## Background

v5.1 measured the Delta-AUC Information Leak (DIAL) for five classifiers
on five TCGA-derived cohorts (THCA, SKCM, LGG, LUAD, COAD). THCA was the
only cohort where linear/tree classifiers flipped labels under
adversarial batch balancing (4/5 classifiers, DIAL 0.323 - 0.494). To
formalize that this is "THCA-specific" rather than pan-cancer
heterogeneity with a heavy tail, we applied random-effects (RE)
meta-analysis in the style of METASOFT (Han & Eskin, 2011, AJHG; Han,
2012, AJHG; Lee et al., 2017, Bioinformatics).

## Methods

Per-classifier, per-cancer DIAL is the effect size theta_i (i=1..5).
SE(DIAL) is approximated by Hanley-McNeil SE(AUC_post) using class
counts from v5p1_harmonization.tsv; DIAL is a shifted/flipped function
of AUC_post so this is valid to first order. SE is floored at 1e-4 to
avoid divide-by-zero under perfect separation. Between-study tau^2 uses
DerSimonian-Laird; Cochran's Q is tested against chi^2 with k-1=4 df;
I^2 = max(0, (Q-df)/Q).

Han m-values use a Gaussian two-component mixture with flat 0.5 prior:
P(data | effect) at the RE pooled mean with variance SE_i^2 + tau^2 vs
P(data | no effect) at zero with variance SE_i^2. This is the analytic
analog of the Han 2012 MSP posterior and matches what METASOFT reports
per input cohort.

## Results

| classifier        | mean  | tau^2  | Q       | p(Q)   | I^2   | m_THCA | m_SKCM | m_LGG | m_LUAD | m_COAD |
|-------------------|------:|-------:|--------:|-------:|------:|-------:|-------:|------:|-------:|-------:|
| LogReg_l2         | 0.099 | 0.1141 | 4538.84 | <1e-16 | 99.9% | 1.00   | 0.10   | 0.02  | 0.04   | 0.03   |
| RandomForest      | 0.065 | 0.0166 |  184.10 | <1e-16 | 97.8% | 1.00   | 0.22   | 0.06  | 0.16   | 0.10   |
| XGBoost           | 0.001 | 0.0000 |    0.12 | 0.998  |  0.0% | 0.50   | 0.50   | 0.50  | 0.50   | 0.50   |
| LogReg_elasticnet | 0.099 | 0.1104 | 4471.88 | <1e-16 | 99.9% | 1.00   | 0.10   | 0.02  | 0.04   | 0.04   |
| GradientBoosting  | 0.076 | 0.0262 |  186.57 | <1e-16 | 97.9% | 1.00   | 0.18   | 0.07  | 0.28   | 0.12   |

Four of five classifiers reject homogeneity overwhelmingly (I^2 > 97%,
p(Q) < 1e-16). In every rejected case THCA has m >= 0.9999, while the
other cohorts fall below or near the conventional Han m <= 0.1 "no
effect" band (Han & Eskin, 2011) - exactly the pattern METASOFT was
designed to flag: one outlier study driving the pooled effect despite
large tau^2.

XGBoost is the lone exception. Q is non-significant (p = 0.998), I^2=0,
and all five m-values collapse to the 0.5 prior because per-cohort
thetas are consistent with zero. This is a negative control: the
meta-analysis refuses to localize an effect when there is none,
consistent with XGBoost never flipping in v5.1.

## Interpretation

RE meta-analysis converts "THCA looks different" into a quantitative
claim: under 4/5 classifiers, posterior P(DIAL effect | data) is
essentially one for THCA and below 0.3 (usually below 0.1) for every
other cohort. With near-ceiling I^2 and the vanishing XGBoost control,
the DIAL signal is THCA-specific and classifier-dependent, not a
pan-cancer phenomenon with a long tail. This motivates S3
(cross-platform replication) and S4 (mechanistic causal analysis)
focused on THCA.

## References

- Han B, Eskin E. Random-effects model for meta-analysis of GWAS. Am J
  Hum Genet 2011;88:586-598.
- Han B, Eskin E. Interpreting meta-analyses of GWAS (m-values / MSP).
  Am J Hum Genet 2012;90:49-60.
- Lee CH, Cook S, Lee JS, Han B. METASOFT comparison of meta-analysis
  estimators. Bioinformatics 2017.
