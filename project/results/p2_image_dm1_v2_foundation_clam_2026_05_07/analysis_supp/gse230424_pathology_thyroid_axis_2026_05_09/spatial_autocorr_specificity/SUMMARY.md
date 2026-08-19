# GSE230424 Spatial-Autocorrelation Specificity Caveat

## Verdict

This is a caveat layer, not a support layer.

The previous matched random-module test showed that DM1/RAI was unusually H&E-predictable among 250 expression/detection-matched random 8-gene modules. This follow-up asks whether that advantage remains after accounting for target spatial smoothness.

It does not.

## Headline Results

| Mode | Observed H&E rho | Expected rho from smoothness | Smoothness-adjusted residual | Percentile vs random modules | Empirical upper-tail p |
|---|---:|---:|---:|---:|---:|
| Raw smoothed target | 0.644 | 0.765 | -0.121 | 1.6% | 0.984 |
| Coord+QC residual target | 0.232 | 0.422 | -0.190 | 0.0% | 1.000 |

## Interpretation

The control regresses matched random-module H&E predictability on neighbor autocorrelation, coordinate-block variance fraction, and target standard deviation. Under that model, DM1/RAI is not above the random-module residual distribution. Its raw and strict residual-target H&E effects remain positive, but the earlier random-module specificity advantage is largely explained by spatial smoothness/block structure.

Safe wording:

- Use GSE230424 as external image-to-spatial-RNA support with QC and spatial-smoothness caveats.
- Do not claim that GSE230424 proves molecular specificity independent of target smoothness.
- Keep the GSE250521 + GSE230424 combined story as Paper 2 support, not Paper 1 MAPK mechanism rescue.

## Files

- `GSE230424_SPATIAL_AUTOCORR_SPECIFICITY_REPORT.md`
- `GSE230424_SPATIAL_AUTOCORR_SPECIFICITY_SUMMARY.json`
- `gse230424_spatial_autocorr_adjusted_specificity.tsv`
- `gse230424_spatial_autocorr_target_metrics.tsv`
- `fig_gse230424_spatial_autocorr_specificity.png`
- `fig_gse230424_spatial_autocorr_specificity.pdf`
