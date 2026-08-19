# Path2Space-inspired nearest-hotspot distance controls

Date: 2026-05-09

Purpose: test whether predicted top-decile high-DM1/RAI spots are spatially proximal to observed high-DM1/RAI hotspots, allowing small localization error, under sample- and spatial-block-preserving nulls.

## Key results

- GSE250521 UNI DM1/RAI median nearest observed-hotspot distance: 1.803 grid spacings.
- GSE250521 block-null mean median distance: 2.566 grid spacings; block-stratified permutation p = 0.001.
- GSE230424 H&E DM1/low-RAI median nearest observed-hotspot distance: 1.000 grid spacing.
- GSE230424 block-null mean median distance: 1.001 grid spacings; block-stratified permutation p = 0.998.
- GSE230424 H&E DM1/low-RAI coord+QC residual median nearest observed-hotspot distance: 2.236 grid spacings.
- GSE230424 residual block-null mean median distance: 1.421 grid spacings; block-stratified permutation p = 1.000.

## Interpretation

This is a mixed, caveat-generating control. GSE250521 passes the spatial-proximity test: predicted hotspots sit closer to observed hotspots than expected under block-preserving nulls. GSE230424 remains strong by exact hotspot overlap and block-centered correlation, but nearest-distance is not supportive after preserving block density, likely because dense spatial blocks make random observed hotspots close to predicted hotspots. Use this as an honest caveat, not as a headline support layer.

## Files

- `PATH2SPACE_HOTSPOT_DISTANCE_CONTROLS_REPORT.md`
- `PATH2SPACE_HOTSPOT_DISTANCE_CONTROLS_SUMMARY.json`
- `path2space_hotspot_distance_summary.tsv`
- `path2space_hotspot_distance_by_group.tsv`
- `path2space_hotspot_distance_null.tsv.gz`
- `fig_path2space_hotspot_distance_controls.png`
- `fig_path2space_hotspot_distance_controls.pdf`
