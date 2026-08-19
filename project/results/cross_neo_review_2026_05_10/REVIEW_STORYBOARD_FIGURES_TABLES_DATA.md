# CROSS-Neo Review Storyboard

## One-Line Story

Personalized neoantigen vaccines are clinically credible again, but the field needs vaccine-ready prioritization: contamination-controlled, source-aware, top-k ranking of the few candidates worth manufacturing.

## Whole Manuscript Flow

1. **Clinical urgency**: recent vaccine studies make candidate selection practically important.
2. **Biological layering**: binding/presentation is necessary but not sufficient for immunogenicity.
3. **Benchmark problem**: overlap, HLA shortcut, source shift and ambiguous negatives can inflate rankings.
4. **Metric contract**: AUPRC, top-k precision, enrichment, calibration and abstention match the vaccine decision better than AUROC alone.
5. **Two-paper bridge**: the review defines the standard; CROSS-Neo is the implementation paper.

## Representative Figures

### Figure 1

File: `/home/seungho/personal/THCA_data_analysis/project/results/cross_neo_review_2026_05_10/figures/fig1_vaccine_ready_prioritization_funnel.png`

Message: candidate neoantigens pass through a funnel, but the final bottleneck is only 5-20 practical vaccine/assay slots.

Where used: opening section and graphical abstract.

### Figure 2

File: `/home/seungho/personal/THCA_data_analysis/project/results/cross_neo_review_2026_05_10/figures/fig2_predictor_taxonomy_by_biological_layer.png`

Message: NetMHCpan/MHCflurry/MixMHCpred/BigMHC-EL model presentation; PRIME/DeepImmuno/BigMHC-IM/TransPHLA/T-SCAPE/IMPROVE move toward immunogenicity or patient context; vaccine triage needs a separate evaluation layer.

Where used: predictor taxonomy section.

### Figure 3

File: `/home/seungho/personal/THCA_data_analysis/project/results/cross_neo_review_2026_05_10/figures/fig3_benchmark_failure_modes_and_controls.png`

Message: exact overlap, near overlap, HLA shortcuts, source shift, ambiguous negatives and AUROC-only reporting each need a matching control.

Where used: benchmark failure-mode section.

### Figure 4

File: `/home/seungho/personal/THCA_data_analysis/project/results/cross_neo_review_2026_05_10/figures/fig4_metric_contract_and_competitor_gauntlet.png`

Message: internal strict demo gauntlet shows how top-k/AUPRC reporting changes practical interpretation compared with public predictor-only baselines.

Where used: metrics section and optional supplementary CROSS-Neo bridge. Claim boundary: internal retrospective example, not external validation.

### Figure 5

File: `/home/seungho/personal/THCA_data_analysis/project/results/cross_neo_review_2026_05_10/figures/fig5_review_to_original_launch_flow.png`

Message: review first defines the field standard; original CROSS-Neo paper implements it.

Where used: author strategy, cover pitch, final review outlook.

## Representative Tables

### Table 1. Clinical vaccine signals

Purpose: show why the topic is timely. Use melanoma, pancreatic cancer and renal cell carcinoma neoantigen vaccine studies.

### Table 2. Predictor families

Purpose: separate binding/presentation predictors from immunogenicity, patient-context and triage tools.

### Table 3. Benchmark failure modes

Purpose: make contamination and source shift a field-wide methodological issue.

### Table 4. Vaccine-ready metrics

Purpose: explain why AUPRC, top-k precision, enrichment and abstention should be primary.

### Supplementary Table 1. Data-use map

Purpose: show exactly what data is used where and what claim each data source can support.

## Data Provenance

The full data-use map is in `review_data_use_map.tsv`. Critical boundaries:

- Public predictor scores are used as literature context or benchmark comparators, not clean CROSS-Neo training features.
- CROSS-Neo product/demo metrics are internal retrospective examples only.
- Clinical vaccine papers motivate the problem; they do not validate CROSS-Neo.
- Public overlap/source-shift audits are used to motivate benchmark standards, not to attack public predictors.

## Editor Pitch Version

This is not another algorithm review. It is a translational evaluation standard for vaccine-ready neoantigen prioritization.
