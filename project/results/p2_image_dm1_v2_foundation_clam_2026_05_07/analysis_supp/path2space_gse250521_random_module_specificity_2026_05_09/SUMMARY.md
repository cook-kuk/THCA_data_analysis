# GSE250521 Random-Module Specificity

## Verdict

This is the GSE250521 specificity counterpart to the GSE230424 caveat control.

## Headline Results

| Mode | Observed UNI rho | Random-module percentile | Empirical p | Smoothness-expected rho | Smoothness-adjusted residual | Smoothness residual percentile |
|---|---:|---:|---:|---:|---:|---:|
| Raw smoothed target | 0.418 | 89.6% | 0.1076 | 0.491 | -0.073 | 9.2% |
| Coord+QC residual target | 0.213 | 100.0% | 0.0040 | 0.168 | 0.045 | 90.4% |

## Interpretation

Raw DM1/RAI is image-predictable, but it is not formally exceptional after random-module or smoothness adjustment. The stricter coord+QC residual target is the useful layer: it remains above all 250 matched random modules before smoothness adjustment and stays directionally high after smoothness adjustment, although the adjusted upper-tail p is 0.0996 rather than conventionally significant.

Safe wording: GSE250521 supports a residual image-aligned DM1/RAI component beyond coordinate/QC structure, but raw panel predictability is partly a spatial-smoothness/tissue-state property.

## Files

- `GSE250521_RANDOM_MODULE_SPECIFICITY_REPORT.md`
- `GSE250521_RANDOM_MODULE_SPECIFICITY_SUMMARY.json`
- `gse250521_random_module_specificity_results.tsv`
- `gse250521_spatial_autocorr_adjusted_specificity.tsv`
- `fig_gse250521_random_module_specificity.png`
- `fig_gse250521_random_module_specificity.pdf`
