# CROSS-Neo-TCR Quick Subset Comparison

## Claim Boundary

This is a fast diagnostic pilot. It is not a final TCR-aware SOTA benchmark and does not transfer public TCR labels into the main CROSS-Neo task.

## Pilot A: Main CROSS-Neo Labels With TCR Evidence Features

Question: on existing CROSS-Neo rows, does adding TCR registry evidence improve ranking compared with pMHC-only features?

| split_name                       | model_name                                    |   n_mean |   n_pos_mean |   prevalence_mean |   auprc_mean |   auroc_mean |   precision_at_20_mean |   enrichment_at_20_mean |
|:---------------------------------|:----------------------------------------------|---------:|-------------:|------------------:|-------------:|-------------:|-----------------------:|------------------------:|
| exact_peptide_hla_holdout        | tcr_evidence_exsource_only                    |     17.8 |          4.2 |          0.235948 |     0.790139 |     0.915302 |               0.235948 |                1        |
| exact_peptide_hla_holdout        | tcr_evidence_only                             |     17.8 |          4.2 |          0.235948 |     0.790139 |     0.915302 |               0.235948 |                1        |
| exact_peptide_hla_holdout        | pmhc_counterfactual_plus_tcr                  |     17.8 |          4.2 |          0.235948 |     0.756226 |     0.883242 |               0.235948 |                1        |
| exact_peptide_hla_holdout        | pmhc_counterfactual_plus_tcr_exsource         |     17.8 |          4.2 |          0.235948 |     0.756226 |     0.883242 |               0.235948 |                1        |
| exact_peptide_hla_holdout        | pmhc_multimodal_plus_tcr                      |     17.8 |          4.2 |          0.235948 |     0.737984 |     0.858132 |               0.235948 |                1        |
| exact_peptide_hla_holdout        | pmhc_multimodal_plus_tcr_exsource             |     17.8 |          4.2 |          0.235948 |     0.737984 |     0.858132 |               0.235948 |                1        |
| exact_peptide_hla_holdout        | tcr_evidence_exsource_nolabel_only            |     17.8 |          4.2 |          0.235948 |     0.706209 |     0.874918 |               0.235948 |                1        |
| exact_peptide_hla_holdout        | tcr_evidence_nolabel_only                     |     17.8 |          4.2 |          0.235948 |     0.706209 |     0.874918 |               0.235948 |                1        |
| exact_peptide_hla_holdout        | pmhc_multimodal_plus_tcr_exsource_nolabel     |     17.8 |          4.2 |          0.235948 |     0.688262 |     0.823571 |               0.235948 |                1        |
| exact_peptide_hla_holdout        | pmhc_multimodal_plus_tcr_nolabel              |     17.8 |          4.2 |          0.235948 |     0.688262 |     0.823571 |               0.235948 |                1        |
| exact_peptide_hla_holdout        | pmhc_counterfactual_plus_tcr_exsource_nolabel |     17.8 |          4.2 |          0.235948 |     0.66593  |     0.838736 |               0.235948 |                1        |
| exact_peptide_hla_holdout        | pmhc_counterfactual_plus_tcr_nolabel          |     17.8 |          4.2 |          0.235948 |     0.66593  |     0.838736 |               0.235948 |                1        |
| exact_peptide_hla_holdout        | pmhc_multimodal                               |     17.8 |          4.2 |          0.235948 |     0.622409 |     0.754451 |               0.235948 |                1        |
| exact_peptide_hla_holdout        | pmhc_counterfactual                           |     17.8 |          4.2 |          0.235948 |     0.588304 |     0.74456  |               0.235948 |                1        |
| near_peptide_cluster_holdout     | tcr_evidence_exsource_only                    |     17.8 |          4.2 |          0.235638 |     0.811667 |     0.92794  |               0.235638 |                1        |
| near_peptide_cluster_holdout     | tcr_evidence_only                             |     17.8 |          4.2 |          0.235638 |     0.811667 |     0.92794  |               0.235638 |                1        |
| near_peptide_cluster_holdout     | tcr_evidence_exsource_nolabel_only            |     17.8 |          4.2 |          0.235638 |     0.786167 |     0.91272  |               0.235638 |                1        |
| near_peptide_cluster_holdout     | tcr_evidence_nolabel_only                     |     17.8 |          4.2 |          0.235638 |     0.786167 |     0.91272  |               0.235638 |                1        |
| near_peptide_cluster_holdout     | pmhc_counterfactual_plus_tcr                  |     17.8 |          4.2 |          0.235638 |     0.721158 |     0.837912 |               0.235638 |                1        |
| near_peptide_cluster_holdout     | pmhc_counterfactual_plus_tcr_exsource         |     17.8 |          4.2 |          0.235638 |     0.721158 |     0.837912 |               0.235638 |                1        |
| near_peptide_cluster_holdout     | pmhc_multimodal_plus_tcr                      |     17.8 |          4.2 |          0.235638 |     0.647098 |     0.793242 |               0.235638 |                1        |
| near_peptide_cluster_holdout     | pmhc_multimodal_plus_tcr_exsource             |     17.8 |          4.2 |          0.235638 |     0.647098 |     0.793242 |               0.235638 |                1        |
| near_peptide_cluster_holdout     | pmhc_counterfactual_plus_tcr_exsource_nolabel |     17.8 |          4.2 |          0.235638 |     0.557737 |     0.736154 |               0.235638 |                1        |
| near_peptide_cluster_holdout     | pmhc_counterfactual_plus_tcr_nolabel          |     17.8 |          4.2 |          0.235638 |     0.557737 |     0.736154 |               0.235638 |                1        |
| near_peptide_cluster_holdout     | pmhc_multimodal_plus_tcr_exsource_nolabel     |     17.8 |          4.2 |          0.235638 |     0.507514 |     0.724615 |               0.235638 |                1        |
| near_peptide_cluster_holdout     | pmhc_multimodal_plus_tcr_nolabel              |     17.8 |          4.2 |          0.235638 |     0.507514 |     0.724615 |               0.235638 |                1        |
| near_peptide_cluster_holdout     | pmhc_multimodal                               |     17.8 |          4.2 |          0.235638 |     0.455481 |     0.662527 |               0.235638 |                1        |
| near_peptide_cluster_holdout     | pmhc_counterfactual                           |     17.8 |          4.2 |          0.235638 |     0.393079 |     0.602747 |               0.235638 |                1        |
| repeated_stratified_5x5_internal | tcr_evidence_exsource_nolabel_only            |     17.8 |          4.2 |          0.235948 |     0.773694 |     0.86878  |               0.235948 |                1        |
| repeated_stratified_5x5_internal | tcr_evidence_nolabel_only                     |     17.8 |          4.2 |          0.235948 |     0.773694 |     0.86878  |               0.235948 |                1        |
| repeated_stratified_5x5_internal | tcr_evidence_exsource_only                    |     17.8 |          4.2 |          0.235948 |     0.771783 |     0.887302 |               0.235948 |                1        |
| repeated_stratified_5x5_internal | tcr_evidence_only                             |     17.8 |          4.2 |          0.235948 |     0.771783 |     0.887302 |               0.235948 |                1        |
| repeated_stratified_5x5_internal | pmhc_counterfactual_plus_tcr                  |     17.8 |          4.2 |          0.235948 |     0.673669 |     0.842934 |               0.235948 |                1        |
| repeated_stratified_5x5_internal | pmhc_counterfactual_plus_tcr_exsource         |     17.8 |          4.2 |          0.235948 |     0.673669 |     0.842934 |               0.235948 |                1        |
| repeated_stratified_5x5_internal | pmhc_multimodal_plus_tcr                      |     17.8 |          4.2 |          0.235948 |     0.650458 |     0.830681 |               0.235948 |                1        |
| repeated_stratified_5x5_internal | pmhc_multimodal_plus_tcr_exsource             |     17.8 |          4.2 |          0.235948 |     0.650458 |     0.830681 |               0.235948 |                1        |
| repeated_stratified_5x5_internal | pmhc_multimodal_plus_tcr_exsource_nolabel     |     17.8 |          4.2 |          0.235948 |     0.598843 |     0.793099 |               0.235948 |                1        |
| repeated_stratified_5x5_internal | pmhc_multimodal_plus_tcr_nolabel              |     17.8 |          4.2 |          0.235948 |     0.598843 |     0.793099 |               0.235948 |                1        |
| repeated_stratified_5x5_internal | pmhc_counterfactual_plus_tcr_exsource_nolabel |     17.8 |          4.2 |          0.235948 |     0.56935  |     0.769451 |               0.235948 |                1        |
| repeated_stratified_5x5_internal | pmhc_counterfactual_plus_tcr_nolabel          |     17.8 |          4.2 |          0.235948 |     0.56935  |     0.769451 |               0.235948 |                1        |
| repeated_stratified_5x5_internal | pmhc_multimodal                               |     17.8 |          4.2 |          0.235948 |     0.5533   |     0.736956 |               0.235948 |                1        |
| repeated_stratified_5x5_internal | pmhc_counterfactual                           |     17.8 |          4.2 |          0.235948 |     0.45317  |     0.645516 |               0.235948 |                1        |
| source_heldout_NEPdb             | tcr_evidence_only                             |    572   |        151   |          0.263986 |     0.964049 |     0.962908 |               1        |                3.78808  |
| source_heldout_NEPdb             | pmhc_multimodal_plus_tcr                      |    572   |        151   |          0.263986 |     0.876265 |     0.919979 |               1        |                3.78808  |
| source_heldout_NEPdb             | pmhc_counterfactual_plus_tcr                  |    572   |        151   |          0.263986 |     0.847281 |     0.892876 |               1        |                3.78808  |
| source_heldout_NEPdb             | tcr_evidence_exsource_only                    |    572   |        151   |          0.263986 |     0.523278 |     0.808411 |               0.65     |                2.46225  |
| source_heldout_NEPdb             | pmhc_multimodal_plus_tcr_exsource             |    572   |        151   |          0.263986 |     0.481586 |     0.796181 |               0.4      |                1.51523  |
| source_heldout_NEPdb             | pmhc_multimodal_plus_tcr_nolabel              |    572   |        151   |          0.263986 |     0.47127  |     0.599062 |               0.95     |                3.59868  |
| source_heldout_NEPdb             | pmhc_counterfactual_plus_tcr_exsource         |    572   |        151   |          0.263986 |     0.430696 |     0.775118 |               0.2      |                0.757616 |
| source_heldout_NEPdb             | tcr_evidence_nolabel_only                     |    572   |        151   |          0.263986 |     0.411305 |     0.446996 |               1        |                3.78808  |
| source_heldout_NEPdb             | pmhc_counterfactual_plus_tcr_nolabel          |    572   |        151   |          0.263986 |     0.406595 |     0.489311 |               0.95     |                3.59868  |
| source_heldout_NEPdb             | pmhc_multimodal                               |    572   |        151   |          0.263986 |     0.30528  |     0.533907 |               0.4      |                1.51523  |
| source_heldout_NEPdb             | tcr_evidence_exsource_nolabel_only            |    572   |        151   |          0.263986 |     0.304573 |     0.475437 |               0.6      |                2.27285  |
| source_heldout_NEPdb             | pmhc_multimodal_plus_tcr_exsource_nolabel     |    572   |        151   |          0.263986 |     0.280307 |     0.538091 |               0.15     |                0.568212 |
| source_heldout_NEPdb             | pmhc_counterfactual_plus_tcr_exsource_nolabel |    572   |        151   |          0.263986 |     0.259301 |     0.508848 |               0        |                0        |
| source_heldout_NEPdb             | pmhc_counterfactual                           |    572   |        151   |          0.263986 |     0.233027 |     0.458881 |               0.2      |                0.757616 |

### Rows With Any TCR Evidence

| split_name           | model_name                                    |   n_mean |   n_pos_mean |   prevalence_mean |   auprc_mean |   auroc_mean |   precision_at_20_mean |   enrichment_at_20_mean |
|:---------------------|:----------------------------------------------|---------:|-------------:|------------------:|-------------:|-------------:|-----------------------:|------------------------:|
| source_heldout_NEPdb | tcr_evidence_only                             |      567 |          146 |          0.257496 |     0.962315 |     0.961637 |                   1    |                3.88356  |
| source_heldout_NEPdb | pmhc_multimodal_plus_tcr                      |      567 |          146 |          0.257496 |     0.86893  |     0.917239 |                   1    |                3.88356  |
| source_heldout_NEPdb | pmhc_counterfactual_plus_tcr                  |      567 |          146 |          0.257496 |     0.838578 |     0.889207 |                   1    |                3.88356  |
| source_heldout_NEPdb | tcr_evidence_exsource_only                    |      567 |          146 |          0.257496 |     0.518778 |     0.80486  |                   0.65 |                2.52432  |
| source_heldout_NEPdb | pmhc_multimodal_plus_tcr_exsource             |      567 |          146 |          0.257496 |     0.467561 |     0.792585 |                   0.4  |                1.55342  |
| source_heldout_NEPdb | pmhc_multimodal_plus_tcr_nolabel              |      567 |          146 |          0.257496 |     0.438467 |     0.585364 |                   0.8  |                3.10685  |
| source_heldout_NEPdb | pmhc_counterfactual_plus_tcr_exsource         |      567 |          146 |          0.257496 |     0.417978 |     0.770735 |                   0.2  |                0.776712 |
| source_heldout_NEPdb | tcr_evidence_nolabel_only                     |      567 |          146 |          0.257496 |     0.382054 |     0.428058 |                   1    |                3.88356  |
| source_heldout_NEPdb | pmhc_counterfactual_plus_tcr_nolabel          |      567 |          146 |          0.257496 |     0.372579 |     0.471822 |                   0.9  |                3.49521  |
| source_heldout_NEPdb | tcr_evidence_exsource_nolabel_only            |      567 |          146 |          0.257496 |     0.295325 |     0.460482 |                   0.6  |                2.33014  |
| source_heldout_NEPdb | pmhc_multimodal                               |      567 |          146 |          0.257496 |     0.293162 |     0.528796 |                   0.4  |                1.55342  |
| source_heldout_NEPdb | pmhc_multimodal_plus_tcr_exsource_nolabel     |      567 |          146 |          0.257496 |     0.253815 |     0.523623 |                   0.05 |                0.194178 |
| source_heldout_NEPdb | pmhc_counterfactual_plus_tcr_exsource_nolabel |      567 |          146 |          0.257496 |     0.244954 |     0.49481  |                   0    |                0        |
| source_heldout_NEPdb | pmhc_counterfactual                           |      567 |          146 |          0.257496 |     0.225034 |     0.457025 |                   0.15 |                0.582534 |

## Pilot B: Paired TCR-pMHC Positives Versus Shuffled TCR Decoys

Question: on paired class-I public positives, can simple sequence and peptide-TCR cross features distinguish observed TCR-pMHC pairs from same-pMHC shuffled-TCR decoys?

| split_name             | model_name     |   n_mean |   n_pos_mean |   prevalence_mean |   auprc_mean |   auroc_mean |   precision_at_20_mean |   enrichment_at_20_mean |
|:-----------------------|:---------------|---------:|-------------:|------------------:|-------------:|-------------:|-----------------------:|------------------------:|
| grouped_by_peptide_hla | pmhc_only      |     8000 |         4000 |               0.5 |     0.5      |      0.5     |                   0.58 |                    1.16 |
| grouped_by_peptide_hla | pmhc_tcr       |     8000 |         4000 |               0.5 |     0.5      |      0.5     |                   0.57 |                    1.14 |
| grouped_by_peptide_hla | tcr_only       |     8000 |         4000 |               0.5 |     0.5      |      0.5     |                   0.58 |                    1.16 |
| grouped_by_peptide_hla | pmhc_tcr_cross |     8000 |         4000 |               0.5 |     0.496714 |      0.49336 |                   0.55 |                    1.1  |

## Interpretation

- A TCR gain in Pilot A is diagnostic evidence that local TCR evidence overlaps with useful ranking signal, not proof of universal TCR recognition prediction.
- Pilot B is a decoy stress test. It checks whether the current simple TCR sequence features are non-random before running expensive structure jobs.
- Strong claims still require paired TCR alpha/beta data, strict TCR/peptide-HLA/source-heldout splits, and external validation.

## Outputs

- `cross_neo_tcr_evidence_pilot_metrics.tsv`
- `cross_neo_tcr_evidence_pilot_predictions.parquet`
- `paired_tcr_decoy_pilot_metrics.tsv`
- `paired_tcr_decoy_pilot_predictions.parquet`
- `quick_tcr_subset_compare_aggregate.tsv`
- `quick_tcr_subset_compare_key_deltas.tsv`
