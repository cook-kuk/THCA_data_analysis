# Pathology scout 2 - UNI image score as thyroid dedifferentiation morphology axis

## Verdict

- Verdict: **RAW_ONLY**
- Raw Spearman with RNA dedifferentiation proxy: 0.527 (p=4.19e-05)
- After DM-label residualization: 0.184 (p=0.182)
- After full residualization: 0.107 (p=0.44)
- Top-vs-bottom dediff tertile AUC: 0.836

## Interpretation

The raw image-score association with dedifferentiation is real, but it does not survive conservative residualization. Treat it as explanation for the image-DM1 model, not as a new independent finding.

## Output files

- `pathology_tds_axis_merged.tsv`
- `pathology_tds_axis_spearman.tsv`
- `pathology_tds_axis_auc.tsv`
- `pathology_tds_axis_covariate_r2.tsv`
- `fig_pathology_tds_axis.png`
