# GSE230424 pathology-to-thyroid-axis screen

Path2Space-style, label-free external thyroid Visium control using GSE230424 raw H&E, matrix, feature, barcode, and tissue-position files.

Key result: 15,489 spots across 4 slides. Supplementary Table S5 maps P1/P2 to PTC+HT and P3/P4 to HT, although disease-group size is only n=2 per group. The top axis is `DM1_low_RAI_score`: H&E tile features reach sample-centered LOSO Spearman rho 0.644 and domain-centered rho 0.910. Within-sample permutation p = 0.001.

Critical caveat: QC-only is stronger for the same axis (sample-centered rho 0.732; coord+QC 0.723). After coordinate+QC residualization, H&E alignment remains but is smaller (rho 0.283). A stricter residual-target control gives rho 0.232. Disease labels are now recovered from Supplementary Table S5, but n=2 versus n=2 is too small for robust disease-group validation. This remains Paper 2 support, not Paper 1 mechanism evidence.

Primary artifacts:
- `GSE230424_PATHOLOGY_THYROID_REPORT.md`
- `GSE230424_PATHOLOGY_THYROID_SUMMARY.json`
- `gse230424_model_comparison.tsv`
- `gse230424_residual_alignment.tsv`
- `gse230424_domain_summary.tsv`
- `fig_gse230424_pathology_thyroid_axis.png`
- `fig_gse230424_spatial_maps.png`
