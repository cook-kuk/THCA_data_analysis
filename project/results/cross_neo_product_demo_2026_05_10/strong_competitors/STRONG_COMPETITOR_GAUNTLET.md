# Strong Competitor Gauntlet

This is a business/product benchmark artifact. Public predictor scores are used only as competitors or predictor-assisted context, not as clean manuscript training features.

## Strict Product Demo Set
| method                            | method_family                 |   n_scored |    AUPRC |    AUROC |   top5_precision |   top10_precision |   top20_precision |   coverage |
|:----------------------------------|:------------------------------|-----------:|---------:|---------:|-----------------:|------------------:|------------------:|-----------:|
| Product_testset_aware_upper_bound | CROSS_Neo_product_upper_bound |         89 | 0.63373  | 0.777311 |              0.8 |               0.9 |              0.5  |          1 |
| v2_selective_exact                | CROSS_Neo_v2_internal         |         89 | 0.636885 | 0.759454 |              1   |               0.7 |              0.65 |          1 |
| v2_selective_hla                  | CROSS_Neo_v2_internal         |         89 | 0.593558 | 0.77381  |              0.8 |               0.7 |              0.55 |          1 |
| Product_fixed_pan_allele_score    | CROSS_Neo_product             |         89 | 0.5822   | 0.77451  |              0.8 |               0.7 |              0.5  |          1 |
| v2_selective_near                 | CROSS_Neo_v2_internal         |         89 | 0.518669 | 0.722689 |              0.8 |               0.7 |              0.4  |          1 |
| hard_decoy_near_meta              | CROSS_Neo_v1_internal         |         89 | 0.538574 | 0.709384 |              1   |               0.6 |              0.4  |          1 |
| hard_decoy_rule_hla               | CROSS_Neo_v1_internal         |         89 | 0.520642 | 0.723389 |              0.6 |               0.6 |              0.5  |          1 |
| MHCflurry_presentation            | public_predictor              |         89 | 0.488891 | 0.658263 |              0.6 |               0.6 |              0.45 |          1 |
| MHCflurry_wave11                  | public_predictor              |         89 | 0.488891 | 0.658263 |              0.6 |               0.6 |              0.45 |          1 |
| Public_binding_mean               | public_predictor_panel        |         89 | 0.441142 | 0.644958 |              0.8 |               0.6 |              0.4  |          1 |
| TransPHLA                         | public_predictor              |         89 | 0.270541 | 0.543417 |              0.6 |               0.6 |              0.3  |          1 |
| TransPHLA_wave3                   | public_predictor              |         89 | 0.270541 | 0.543417 |              0.6 |               0.6 |              0.3  |          1 |
| v1_gated_cqk_hla                  | CROSS_Neo_v1_internal         |         89 | 0.533398 | 0.752101 |              0.8 |               0.5 |              0.45 |          1 |
| MHCflurry_affinity_inverse        | public_predictor              |         89 | 0.43597  | 0.584034 |              0.6 |               0.5 |              0.35 |          1 |
| BigMHC_EL                         | public_predictor              |         89 | 0.367573 | 0.617297 |              0.4 |               0.5 |              0.45 |          1 |

## Bottom Line
- Fixed pan-allele product score: AUPRC 0.582, top10 precision 0.700.
- Best public-only competitor in this gauntlet: MHCflurry_presentation, AUPRC 0.489, top10 precision 0.600.
- Product_testset_aware_upper_bound is internal rehearsal only and should not be used as a generalization claim.
- The right customer-facing claim is stronger practical prioritization on an internal retrospective demo set, not external validation.

## Strong Competitors Included
- BigMHC, DeepImmuno, PRIME, TransPHLA, MHCflurry, NetMHCpan, NetMHCstabpan, MHCnuggets, T-SCAPE.
- CROSS-Neo v1/v2 internal challengers, including hard-decoy, gated C+QK, selective ensemble, and ESM2/QK gate summaries.

## Strong Competitors Still Missing Locally
- IMPROVE: strong broad-scale neoepitope immunogenicity model, but no local score file was found.
- MixMHCpred standalone: indirectly represented through PRIME-style scoring, but no standalone local output was found for the strict set.
- DeepHLApan / CIImm / GraphMHC: useful literature comparators, but no local runnable score files were found.

## Files
- Score matrix: `/home/seungho/personal/THCA_data_analysis/project/results/cross_neo_product_demo_2026_05_10/strong_competitors/strong_competitor_score_matrix.tsv`
- Gauntlet metrics: `/home/seungho/personal/THCA_data_analysis/project/results/cross_neo_product_demo_2026_05_10/strong_competitors/strong_competitor_gauntlet.tsv`
- Context public benchmark: `/home/seungho/personal/THCA_data_analysis/project/results/cross_neo_product_demo_2026_05_10/strong_competitors/strong_competitor_context_public_benchmarks.tsv`
- Missing competitor manifest: `/home/seungho/personal/THCA_data_analysis/project/results/cross_neo_product_demo_2026_05_10/strong_competitors/missing_strong_competitor_manifest.tsv`
- Figure: `/home/seungho/personal/THCA_data_analysis/project/results/cross_neo_product_demo_2026_05_10/strong_competitors/figure_strong_competitor_gauntlet.png`
