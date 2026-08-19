# Path2Space-inspired bootstrap uncertainty

Date: 2026-05-09

Purpose: add descriptive uncertainty intervals for the main Paper 2 image-to-spatial-RNA evidence layers using 2,000 stratified bootstrap replicates.

## Key results

- GSE250521 UNI DM1/RAI smoothed spot-centered rho: 0.440; 95% CI 0.408 to 0.475.
- GSE250521 UNI DM1/RAI coordinate-domain centered rho: 0.581; 95% CI 0.381 to 0.729.
- GSE230424 H&E DM1/low-RAI smoothed spot-centered rho: 0.644; 95% CI 0.634 to 0.655.
- GSE230424 H&E DM1/low-RAI coordinate-domain centered rho: 0.910; 95% CI 0.788 to 0.956.
- GSE230424 H&E DM1/low-RAI coord+QC residual-target rho: 0.232; 95% CI 0.217 to 0.248.

## Interpretation

The main image-to-spatial-RNA evidence layers remain above zero under stratified bootstrap uncertainty. These intervals support reviewer-facing stability claims, while preserving the boundary that bootstrap uncertainty is descriptive and does not replace independent cohort validation.

## Files

- `PATH2SPACE_BOOTSTRAP_UNCERTAINTY_REPORT.md`
- `PATH2SPACE_BOOTSTRAP_UNCERTAINTY_SUMMARY.json`
- `path2space_bootstrap_uncertainty_summary.tsv`
- `path2space_bootstrap_uncertainty_replicates.tsv.gz`
- `fig_path2space_bootstrap_uncertainty.png`
- `fig_path2space_bootstrap_uncertainty.pdf`
