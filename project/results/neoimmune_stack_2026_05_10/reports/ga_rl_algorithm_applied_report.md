# GA/RL Algorithm Applied to NeoImmune-Stack

## What was used
- GA score file: `/home/seungho/personal/THCA_data_analysis/project/results/cross_neo_kg_ga_evolutionary_ensemble_2026_05_10/kg_ga_evolved_candidate_scores.tsv`
- Algorithm: `KG_GA_evolved_controller` from knowledge-graph genetic/evolutionary architecture search.
- Related RL artifact: Darwin/RL blueprint remains a method-search blueprint, not a trained scoring model.

## Boundary
- This branch is **production/experiment-priority only**.
- It is not allowed in the clean science track because the evolved controller includes BigMHC/public predictor and impact/product components.
- Claim-safe wording: GA/RL found a retrospective candidate-prioritization controller ready for prospective assay validation.

## Result
- Integrated GA/RL best row: KG_GA_evolved_controller AUPRC=0.877, AUROC=0.910, patient_hit_rate@34=0.915.

## Feature weights
| feature          | source_column                       | concept                     | role                        | layer   |      weight | transform   |
|:-----------------|:------------------------------------|:----------------------------|:----------------------------|:--------|------------:|:------------|
| MD_control       | md_control_score_norm               | structure_control_readiness | claim_safety                | atomic  | 0.266764    | square      |
| Impact_portfolio | impact_portfolio_score              | translation_impact          | product_value               | atomic  | 0.22        | sigmoid     |
| TCR_expert       | tcr_recognition_score_norm          | tcr_recognition             | pmhc_tcr_likelihood         | atomic  | 0.179852    | sqrt        |
| Foreignness      | fitness_foreignness_proxy           | mutant_foreignness          | self_nonself_gap            | atomic  | 0.177546    | sigmoid     |
| BigMHC_IM        | bigmhc_im_score                     | public_immunogenicity_prior | presentation_immunogenicity | atomic  | 0.0906103   | identity    |
| CROSS_BMA        | bma_v2_discovery_score              | bayesian_model_average      | ensemble_stability          | atomic  | 0.0316841   | square      |
| BigMHC_EL        | bigmhc_el_score                     | public_presentation_prior   | hla_presentation            | atomic  | 0.0187562   | square      |
| CROSS_finetuned  | finetuned_experiment_priority_score | experiment_priority         | wetlab_priority             | atomic  | 0.00604986  | square      |
| CROSS_core       | crossneo_core_score_norm            | internal_ranker             | candidate_ranking           | atomic  | 0.00414097  | square      |
| CROSS_stress     | stress_guarded_discovery_score      | stress_guarded_ranker       | robust_candidate_ranking    | atomic  | 0.00286742  | square      |
| Fixed_integrated | immunogenicity_discovery_score      | fixed_integrated_score      | legacy_controller           | meta    | 0.000864304 | sqrt        |
| Fixed_claimsafe  | immunogenicity_claim_safe_score     | fixed_claim_safe_score      | legacy_claim_guard          | meta    | 0.000864304 | sigmoid     |

## Gates and synergies
| gate                                       |        value | synergy                      | left            | right       |      weight |
|:-------------------------------------------|-------------:|:-----------------------------|:----------------|:------------|------------:|
| high_leakage_penalty                       |   0.00315274 | nan                          | nan             | nan         | nan         |
| medium_leakage_penalty                     |   0.0244854  | nan                          | nan             | nan         | nan         |
| low_presentation_threshold                 |   0.247009   | nan                          | nan             | nan         | nan         |
| low_presentation_penalty                   |   0.0958028  | nan                          | nan             | nan         | nan         |
| claim_blocker_penalty                      |   0.240232   | nan                          | nan             | nan         | nan         |
| single_concept_default_weight_cap          |   0.3        | nan                          | nan             | nan         | nan         |
| single_concept_weight_cap:TCR_expert       |   0.18       | nan                          | nan             | nan         | nan         |
| single_concept_weight_cap:Impact_portfolio |   0.22       | nan                          | nan             | nan         | nan         |
| single_concept_weight_cap:Fixed_integrated |   0.2        | nan                          | nan             | nan         | nan         |
| single_concept_weight_cap:Fixed_claimsafe  |   0.2        | nan                          | nan             | nan         | nan         |
| nan                                        | nan          | presentation_x_tcr           | BigMHC_EL       | TCR_expert  |   0.232469  |
| nan                                        | nan          | immunogenicity_x_foreignness | BigMHC_IM       | Foreignness |   0         |
| nan                                        | nan          | stress_x_presentation        | CROSS_stress    | BigMHC_EL   |   0.0151821 |
| nan                                        | nan          | ranker_x_tcr                 | CROSS_core      | TCR_expert  |   0.196157  |
| nan                                        | nan          | claim_safety_x_foreignness   | MD_control      | Foreignness |   0.0328476 |
| nan                                        | nan          | claimsafe_x_public_im        | Fixed_claimsafe | BigMHC_IM   |   0         |

## Outputs
- `metrics/ga_rl_algorithm_metrics.tsv`
- `predictions/kg_ga_evolved_ranked_candidates.tsv`
- `predictions/patient_top20_candidates_kg_ga.tsv`
- `predictions/patient_top34_candidates_kg_ga.tsv`
