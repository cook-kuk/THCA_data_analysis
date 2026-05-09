# Pathology finding scout - spatial-DM1 morphology transfer

## Verdict

- Verdict: **SPATIAL_ONLY**
- Spatial leave-slide-out Spearman rho: 0.301
- Spatial median per-slide Spearman rho: 0.317
- TCGA DM1 vs DM2 AUC from spatial-trained morphology score: 0.550
- TSS-residualized TCGA DM1 vs DM2 AUC: 0.503
- Spearman versus UNI LOTO DM1 probability: 0.175

## Candidate finding

H&E tiles carry within-spatial-cohort DM1 transcriptomic information, but the transfer to TCGA DM1 labels is too weak or confounded for a new finding.

## Methods

- Fit a PCA(32)+ridge model from GSE250521 UNI tile embeddings to spot-level DM1_like_score, evaluated by leave-one-slide-out prediction.
- Refit the same model on all spatial spots, projected it to TCGA UNI tile features, and summarized each WSI by mean, q90, and high-tile fractions.
- Tested TCGA associations with DM group, molecular subtype, histology, TDS, and conservative residualization against TSS/sex/histology/molecular subtype.

## Key files

- `PATHOLOGY_SCOUT_SUMMARY.json`
- `spatial_spot_loso_predictions.tsv`
- `spatial_slide_summary.tsv`
- `tcga_slide_pathology_scores.tsv`
- `tcga_association_tests.tsv`
- `tcga_confound_residual_tests.tsv`
- `fig_pathology_scout_spatial_loso.png`
- `fig_pathology_scout_tcga_projection.png`
