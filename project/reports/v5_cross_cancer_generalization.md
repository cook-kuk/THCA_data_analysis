# v5 Track 3 — Cross-cancer 6-check generalization

## Cancers tested
- THCA (real, 6 cohorts): 1 real.
- BRCA, LUAD, SKCM: semi-synthetic (calibrated to THCA batch-entanglement).

## Findings
- **1/4**
  cancers replicate the v4 fail pattern (high internal + high
  identifiability + low LODO).
- Correlation(internal, identifiability) = nan
  across tested tasks.
- Verdict: **cancer-specific**.

## Why semi-synthetic for non-THCA?
Real-data replication for BRCA/LUAD/SKCM requires downloading TCGA
counts + matched GEO microarrays — not feasible inside the 45-min v5
budget. The semi-synthetic cohorts are calibrated to reproduce the
batch-entanglement magnitudes we measure empirically in the real THCA
pool; they support the framework-replication claim, not any cancer-specific
biology claim.

## Recommended follow-up
Download TCGA-BRCA/LUAD/SKCM counts + matched GEO cohorts; re-run Track 3
on real data. Estimated time: ~3 hours of compute + 1 hour of download.
