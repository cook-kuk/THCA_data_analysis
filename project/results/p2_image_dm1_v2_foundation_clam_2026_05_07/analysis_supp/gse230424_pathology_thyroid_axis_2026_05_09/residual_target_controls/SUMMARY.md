# GSE230424 residual-target controls

Strict follow-up to the GSE230424 H&E-to-DM1/RAI analysis.

Question: after removing coordinate and ST-QC structure from the observed spatial RNA target within each slide, can H&E tile features directly predict the remaining residual target?

Answer: yes for DM1/low-RAI, but at a smaller and appropriately conservative effect size. The coordinate+QC-residualized `DM1_low_RAI_score` target is predicted from H&E with sample-centered LOSO Spearman rho 0.232, median sample rho 0.195, and within-sample permutation p = 0.001. Coordinate-only, QC-only, and coord+QC baselines are near null on this residual target.

Interpretation: the raw GSE230424 DM1/RAI H&E signal is strongly QC/tissue-density aligned, but a smaller H&E-aligned component remains after coordinate+QC removal. Use as Paper 2 image-to-spatial-RNA support with explicit QC caveat; do not use as Paper 1 mechanism evidence.

Primary artifacts:
- `GSE230424_RESIDUAL_TARGET_REPORT.md`
- `GSE230424_RESIDUAL_TARGET_SUMMARY.json`
- `gse230424_residual_target_model_comparison.tsv`
- `gse230424_residual_target_permutation.tsv`
- `fig_gse230424_residual_target_controls.png`
