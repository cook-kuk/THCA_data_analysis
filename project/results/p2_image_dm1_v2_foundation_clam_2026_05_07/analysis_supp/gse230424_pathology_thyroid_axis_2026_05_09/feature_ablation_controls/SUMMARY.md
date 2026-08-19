# GSE230424 feature-family ablation controls

Date: 2026-05-09

Purpose: test whether the GSE230424 H&E-to-DM1/low-RAI signal is driven by one simple feature family, especially RGB/stain summaries, or whether larger tissue morphology retains signal after coordinate and ST-QC adjustment.

## Key results

- Spots / samples: 15,489 spots across 4 Visium slides.
- Raw DM1/low-RAI all-H&E sample-centered rho: 0.644.
- Best raw feature family: radius96_only, rho 0.645.
- Coord+QC residual-target all-H&E sample-centered rho: 0.232.
- Best residual feature family: radius96_only, rho 0.246.
- Residual RGB-distribution-only rho: 0.143.
- Residual density/texture-only rho: 0.174.

## Interpretation

The GSE230424 DM1/RAI signal is not explained only by a single RGB or stain-summary family. Large-patch morphology carries the strongest residual signal. The evidence strengthens Paper 2 as an image-to-spatial-RNA reinforcement analysis, while preserving the explicit QC/stain caveat because raw QC-only models remain strong in the parent analysis.

## Files

- `GSE230424_FEATURE_ABLATION_REPORT.md`
- `GSE230424_FEATURE_ABLATION_SUMMARY.json`
- `gse230424_feature_ablation_model_comparison.tsv`
- `gse230424_feature_ablation_predictions.tsv.gz`
- `fig_gse230424_feature_ablation_controls.png`
- `fig_gse230424_feature_ablation_controls.pdf`
