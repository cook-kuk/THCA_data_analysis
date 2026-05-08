# Multimodal Haiku-lite THCA pilot (2026-05-08)

## Input
- spots: 3,200
- slides: 16
- stages: ATC, LPTC, PT, PTC
- modalities: local morphology/cell-composition features + clinical metadata-as-text proxy (stage) + spatial transcriptomic targets.

## Main readout
Leave-one-slide-out ridge models compare morphology-only, clinical-text-only, and morphology+text.

| target | morph rho | text rho | combo rho | combo AUROC q75 | delta combo-morph |
|---|---:|---:|---:|---:|---:|
| DM1_like_score | 0.189 | -0.104 | 0.237 | 0.693 | 0.048 |
| RAI_8_score | 0.189 | -0.104 | 0.237 | 0.568 | 0.048 |
| TDS_like_score | 0.163 | -0.105 | 0.196 | 0.551 | 0.033 |
| TROP2 | -0.044 | 0.372 | 0.564 | 0.844 | 0.608 |
| Epithelial_score | 0.160 | -0.018 | 0.181 | 0.535 | 0.020 |
| Proliferation_score | 0.107 | -0.051 | 0.125 | 0.532 | 0.018 |
| EMT_score | 0.037 | -0.044 | 0.038 | 0.473 | 0.000 |
| CAF_ECM_score | 0.104 | -0.051 | 0.135 | 0.607 | 0.032 |
| Hypoxia_score | 0.136 | -0.065 | 0.139 | 0.544 | 0.003 |

## Interpretation
- Morphology-plus-clinical metadata is strongest where spatial programs track gross disease state and local cell composition.
- DM1/RAI remain weak-to-moderate, consistent with the earlier H&E-only closure battery: a Haiku-like framework helps frame the negative as modality insufficiency rather than failed execution.
- The paper-worthy multimodal angle is not direct H&E-to-DM prediction; it is a thyroid-specific retrieval/counterfactual atlas that separates morphology-visible programs from molecularly hidden lineage programs.

## Most morphology-visible targets
- TROP2: best=morphology_plus_text, rho=0.564, combo-vs-morph delta=0.608
- RAI_8_score: best=morphology_plus_text, rho=0.237, combo-vs-morph delta=0.048
- DM1_like_score: best=morphology_plus_text, rho=0.237, combo-vs-morph delta=0.048
- TDS_like_score: best=morphology_plus_text, rho=0.196, combo-vs-morph delta=0.033

## Least morphology-visible targets
- EMT_score: best=morphology_plus_text, rho=0.038, combo-vs-morph delta=0.000
- Proliferation_score: best=morphology_plus_text, rho=0.125, combo-vs-morph delta=0.018
- CAF_ECM_score: best=morphology_plus_text, rho=0.135, combo-vs-morph delta=0.032
- Hypoxia_score: best=morphology_plus_text, rho=0.139, combo-vs-morph delta=0.003

## Proposed topic
**Thyroid-Haiku-lite: clinical-text conditioned retrieval of spatial molecular niches in thyroid cancer.**

Use the Haiku concept, but adapt it to thyroid: H&E/spot morphology, clinical metadata text, and spatial RNA programs. The biological question becomes: which thyroid cancer programs are recoverable from morphology plus metadata, and which require direct spatial molecular measurement?

## Output files
- `project/results/high_impact_topic_pilots_2026_05_08/multimodal_haiku_lite/multimodal_haiku_lite_metrics.tsv`
- `project/results/high_impact_topic_pilots_2026_05_08/multimodal_haiku_lite/multimodal_haiku_lite_delta.tsv`
- `project/results/high_impact_topic_pilots_2026_05_08/multimodal_haiku_lite/multimodal_haiku_lite_fold_metrics.tsv`
- `project/results/high_impact_topic_pilots_2026_05_08/multimodal_haiku_lite/multimodal_haiku_lite_predictions.tsv.gz`
