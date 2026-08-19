# Path2Space-inspired hotspot concordance

Date: 2026-05-09

Purpose: test whether predicted high-DM1/RAI regions recover observed high-DM1/RAI spatial hotspots within each slide/sample, beyond continuous Spearman correlation.

## Key results

- GSE250521 UNI DM1/RAI exact top-decile precision: 0.350.
- GSE250521 top-decile precision lift over expected: 3.50x; Fisher OR 6.92; permutation p = 0.001.
- GSE230424 H&E DM1/low-RAI exact top-decile precision: 0.333.
- GSE230424 top-decile precision lift over expected: 3.32x; Fisher OR 6.22; permutation p = 0.001.
- GSE230424 H&E DM1/low-RAI coord+QC residual exact top-decile precision: 0.163.
- GSE230424 residual-target top-decile precision lift over expected: 1.63x; Fisher OR 1.90; permutation p = 0.001.

## Interpretation

The hotspot concordance layer supports spatial localization, not only rank correlation. Predicted top-decile high-DM1/RAI regions are enriched for observed high-DM1/RAI regions in both GSE250521 and GSE230424. The strict coordinate+QC residual target is weaker, as expected, but still above independence.

## Files

- `PATH2SPACE_HOTSPOT_CONCORDANCE_REPORT.md`
- `PATH2SPACE_HOTSPOT_CONCORDANCE_SUMMARY.json`
- `path2space_hotspot_concordance_summary.tsv`
- `path2space_hotspot_concordance_by_group.tsv`
- `fig_path2space_hotspot_concordance.png`
- `fig_path2space_hotspot_concordance.pdf`
