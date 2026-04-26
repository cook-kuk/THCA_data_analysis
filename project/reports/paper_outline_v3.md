# THCA v3 — paper outline

## Working title

*Leakage-aware compact gene panels for BRAF-like vs RAS-like thyroid cancer: an honest retrospective benchmark across TCGA-THCA + 4 GEO cohorts.*

## Key narrative

1. TierA67-style panels reach AUC ~1.0 on TCGA — but removing MAPK-output genes that are essentially re-labelings of BRAF-like status drops it to **0.924**.
2. On external cohorts (GSE27155, GSE126698, GSE213647), raw transfer is substantially worse than internal CV; isotonic recalibration partially closes the gap.
3. Panel-size curve plateaus around k≈**16** genes for the QUBO-neal strategy; beyond that, more features don't help.
4. In a synthetic Bethesda cohort (prev=0.10–0.30), the best model at Se≥0.95 offers **meaningful (but prevalence-sensitive)** reduction in surgeries vs treat-all — but this is strictly a decision-support prototype.

## Suggested figures

- Fig 1: Leakage curve + MAPK-ablation bar + permutation null histogram.
- Fig 2: Dataset-identifiability OvR AUC + LODO AUC grouped bars.
- Fig 3: Panel-k curve (3 strategies × 3 cohorts) + NPV heatmap.
- Fig 4: Three-class ROC + confusion + SHAP top-15 per class.
- Fig 5: KM curves for molecular subtype + TDS tertile + Cox forest.
- Fig 6: Bethesda operating curves + decision curve + surgery-reduction bar.

## Methods to emphasise

- Feature selection done INSIDE each CV fold (no leakage).
- Isotonic recalibration on external 5-fold splits.
- Permutation null (1000 shuffles) for empirical p-values.
- Bootstrap 95% CIs (1000 iterations).
- Explicit honesty audit page reporting the uncomfortable numbers.

## What we are NOT claiming

- No diagnostic test, no clinical validity.
- Not a prospective study.
- Not a Bethesda-III/IV triage replacement.
