# GSE230424 random-module specificity controls

Date: 2026-05-09

Purpose: test whether the GSE230424 H&E-to-DM1/low-RAI signal is simply a generic spatially smooth gene-module effect, and whether the 8-gene panel signal collapses when any one gene is removed.

## Key results

- Spots / samples: 15,489 spots across 4 Visium slides.
- Matched null: 250 expression/detection-matched random 8-gene modules.
- Actual DM1/low-RAI raw H&E sample-centered rho: 0.644.
- Raw random-module percentile: 99.2%; empirical upper-tail p = 0.0120.
- Actual DM1/low-RAI coord+QC residual-target H&E rho: 0.232.
- Residual-target random-module percentile: 100.0%; empirical upper-tail p = 0.0040.
- Panel leave-one-gene-out raw rho range: 0.628 to 0.654.
- Panel leave-one-gene-out coord+QC residual-target rho range: 0.220 to 0.238.

## Interpretation

This strengthens the Paper 2 claim. The raw H&E signal is partly a broad tissue-state phenomenon because matched random modules are also image-predictable, but the actual DM1/RAI axis sits in the extreme upper tail. After coordinate+QC residualization, DM1/RAI remains above all 250 matched random modules. Leave-one-out stability argues against a single-gene panel artifact.

## Files

- `GSE230424_RANDOM_MODULE_SPECIFICITY_REPORT.md`
- `GSE230424_RANDOM_MODULE_SPECIFICITY_SUMMARY.json`
- `gse230424_random_module_specificity_results.tsv`
- `gse230424_random_matched_modules.tsv`
- `gse230424_gene_expression_metrics.tsv`
- `gse230424_random_module_spot_scores.tsv.gz`
- `gse230424_random_module_actual_predictions.tsv.gz`
- `fig_gse230424_random_module_specificity.png`
- `fig_gse230424_random_module_specificity.pdf`
