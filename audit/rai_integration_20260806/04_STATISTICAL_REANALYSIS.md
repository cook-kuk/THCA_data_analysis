# AUDIT 04 — rare-event reanalysis of the structural outcomes

Generated 2026-08-06 by `scripts/audit/03_rare_event_structural_outcomes.py`

## Why the first-pass odds ratios cannot stand as reported

The first pass fitted ordinary maximum-likelihood logistic regression with up to six covariates to outcomes with 9 and 13 events — about 1.5 to 2 events per parameter. Small-sample ML logistic coefficients are biased away from the null, so the reported odds ratios near 0.19 are expected to be too extreme.

Every model below is Firth-penalised, which removes that bias and remains defined under separation. Covariates are added one at a time; the three-covariate model is marked exploratory rather than presented as primary.

| endpoint                                      | model                                         |   n |   events |    epv |       or_ |      ci_low |    ci_high |         p | method        |
|:----------------------------------------------|:----------------------------------------------|----:|---------:|-------:|----------:|------------:|-----------:|----------:|:--------------|
| new tumour event after initial treatment      | Tier 1 · Cohen's d (patient bootstrap)        | 145 |        9 | nan    | -0.837553 |  -1.3747    |  -0.247422 | 0.0048805 | nonparametric |
| new tumour event after initial treatment      | Tier 1 · rank-biserial correlation            | 145 |        9 | nan    | -0.562092 | nan         | nan        | 0.0048805 | nonparametric |
| new tumour event after initial treatment      | Tier 2 · Firth, score only                    | 145 |        9 |   9    |  0.221161 |   0.0624042 |   0.783795 | 0.0194199 | Firth         |
| new tumour event after initial treatment      | Tier 3 · Firth + stage                        | 145 |        9 |   4.5  |  0.239175 |   0.0689686 |   0.829433 | 0.0241486 | Firth         |
| new tumour event after initial treatment      | Tier 3 · Firth + purity                       | 145 |        9 |   4.5  |  0.20199  |   0.0544222 |   0.749692 | 0.0168225 | Firth         |
| new tumour event after initial treatment      | Tier 3 · Firth + stage + purity (exploratory) | 145 |        9 |   3    |  0.219439 |   0.0603556 |   0.797832 | 0.0212819 | Firth         |
| persistent disease within 3 months of surgery | Tier 1 · Cohen's d (patient bootstrap)        |  72 |       13 | nan    | -0.781851 |  -1.29429   |  -0.311077 | 0.0118015 | nonparametric |
| persistent disease within 3 months of surgery | Tier 1 · rank-biserial correlation            |  72 |       13 | nan    | -0.449804 | nan         | nan        | 0.0118015 | nonparametric |
| persistent disease within 3 months of surgery | Tier 2 · Firth, score only                    |  72 |       13 |  13    |  0.224491 |   0.0610427 |   0.825593 | 0.0245475 | Firth         |
| persistent disease within 3 months of surgery | Tier 3 · Firth + stage                        |  72 |       13 |   6.5  |  0.238041 |   0.0650529 |   0.87104  | 0.0301126 | Firth         |
| persistent disease within 3 months of surgery | Tier 3 · Firth + purity                       |  72 |       13 |   6.5  |  0.223132 |   0.0600414 |   0.829229 | 0.0251175 | Firth         |
| persistent disease within 3 months of surgery | Tier 3 · Firth + stage + purity (exploratory) |  72 |       13 |   4.33 |  0.236591 |   0.0640735 |   0.873611 | 0.0305614 | Firth         |

## Stability

| endpoint                                      |   n_events |   loo_or_min |   loo_or_max | loo_all_below_1   |   boot_median_or |   boot_ci_low |   boot_ci_high |   boot_frac_or_below_1 |   n_boot_used |
|:----------------------------------------------|-----------:|-------------:|-------------:|:------------------|-----------------:|--------------:|---------------:|-----------------------:|--------------:|
| new tumour event after initial treatment      |          9 |        0.121 |        0.282 | True              |            0.208 |         0.033 |          0.774 |                  0.988 |          1997 |
| persistent disease within 3 months of surgery |         13 |        0.156 |        0.274 | True              |            0.218 |         0.032 |          0.631 |                  0.999 |          2000 |

`loo_all_below_1` asks whether the direction survives dropping any single event; `boot_frac_or_below_1` is the fraction of 2,000 patient-level bootstrap refits whose odds ratio stays below one. Neither is a p-value — they describe how much the estimate depends on individual patients.
