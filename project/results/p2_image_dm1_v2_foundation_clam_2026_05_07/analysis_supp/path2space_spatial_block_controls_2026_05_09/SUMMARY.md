# Path2Space-inspired spatial-block controls

Date: 2026-05-09

Purpose: test whether image-to-DM1/RAI correlation and hotspot localization remain after removing broad coordinate-domain structure within each slide/sample.

## Key results

- GSE250521 UNI DM1/RAI spatial-block-centered rho: 0.380.
- GSE250521 top-decile hotspot enrichment under block-stratified permutation: p = 0.001.
- GSE230424 H&E DM1/low-RAI spatial-block-centered rho: 0.497.
- GSE230424 top-decile hotspot enrichment under block-stratified permutation: p = 0.001.
- GSE230424 H&E DM1/low-RAI coord+QC residual spatial-block-centered rho: 0.198.
- GSE230424 residual-target top-decile hotspot enrichment under block-stratified permutation: p = 0.001.

## Interpretation

The signal is not explained only by broad coordinate-domain structure. Block-centering reduces effect sizes, as expected, but the primary and strict residual-target layers remain positive. The block-stratified hotspot null preserves broad predicted hotspot density within spatial blocks, yet observed hotspot overlap remains enriched.

## Files

- `PATH2SPACE_SPATIAL_BLOCK_CONTROLS_REPORT.md`
- `PATH2SPACE_SPATIAL_BLOCK_CONTROLS_SUMMARY.json`
- `path2space_spatial_block_controls_summary.tsv`
- `path2space_spatial_block_controls_by_group.tsv`
- `fig_path2space_spatial_block_controls.png`
- `fig_path2space_spatial_block_controls.pdf`
