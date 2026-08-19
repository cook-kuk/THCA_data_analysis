# Paperability-aware Darwin-RL

Generated: 2026-05-10T23:16:54

## Bottom line

- This layer makes paper-readiness an explicit optimization target, not an afterthought.
- GA/RL should select algorithms that are accurate, validatable, explainable, ablatable and easy to defend.
- Current top candidate by paperability reward: `KG_GA_evolved` (0.922).

## Paperability reward

| reward_component              |   weight | purpose                                                              |
|:------------------------------|---------:|:---------------------------------------------------------------------|
| benchmark_lift                |     0.16 | clear quantitative improvement over public and internal comparators  |
| external_generalization       |     0.15 | held-out source/cohort/HLA/sequence-cluster performance              |
| prospective_lock_readiness    |     0.14 | score and thresholds can be frozen before wetlab or clinical readout |
| novelty_of_method             |     0.13 | not just another ensemble; contains a distinct method contribution   |
| ablation_completeness         |     0.11 | each major primitive has an on/off or replacement control            |
| failure_mode_interpretability |     0.1  | algorithm explains why hits and fails occur                          |
| comparator_coverage           |     0.08 | includes SOTA/public baselines and internal baselines                |
| figure_table_readiness        |     0.07 | naturally yields publishable figures, tables and decision matrices   |
| reviewer_risk_mitigation      |     0.04 | co-located caveats, leakage audits and claim boundaries              |
| deployment_simplicity         |     0.02 | can be run and reproduced without excessive infrastructure           |

## Candidate algorithm ranking

| algorithm        |   all_AUPRC |   all_AUROC |   all_top96_precision |   validation_like_AUPRC |   low_leakage_AUPRC |   v6_smoke_AUPRC |   benchmark_lift |   external_generalization |   prospective_lock_readiness |   novelty_of_method |   ablation_completeness |   failure_mode_interpretability |   comparator_coverage |   figure_table_readiness |   reviewer_risk_mitigation |   deployment_simplicity |   paperability_reward |
|:-----------------|------------:|------------:|----------------------:|------------------------:|--------------------:|-----------------:|-----------------:|--------------------------:|-----------------------------:|--------------------:|------------------------:|--------------------------------:|----------------------:|-------------------------:|---------------------------:|------------------------:|----------------------:|
| KG_GA_evolved    |    0.933411 |    0.943428 |              1        |                0.97775  |            0.757078 |         0.811595 |        1         |                 1         |                         0.95 |                0.95 |                    0.72 |                            0.88 |                  0.95 |                     0.95 |                       0.85 |                    0.6  |              0.9222   |
| CROSS_claimsafe  |    0.860895 |    0.863588 |              1        |                0.84648  |            0.616845 |         0.788084 |        0.746166  |                 0.804122  |                         0.85 |                0.6  |                    0.58 |                            0.82 |                  0.75 |                     0.82 |                       0.88 |                    0.72 |              0.749805 |
| CROSS_integrated |    0.878422 |    0.895868 |              1        |                0.713841 |            0.57277  |         0.790492 |        0.807518  |                 0.606201  |                         0.8  |                0.55 |                    0.55 |                            0.78 |                  0.75 |                     0.82 |                       0.72 |                    0.72 |              0.702733 |
| CROSS_stress     |    0.871905 |    0.873389 |              0.947917 |                0.650789 |            0.600589 |         0.776958 |        0.784705  |                 0.512115  |                         0.7  |                0.5  |                    0.45 |                            0.7  |                  0.6  |                     0.65 |                       0.55 |                    0.8  |              0.61637  |
| CROSS_BMA        |    0.860525 |    0.859654 |              0.927083 |                0.659761 |            0.635656 |         0.766804 |        0.744873  |                 0.525503  |                         0.7  |                0.55 |                    0.45 |                            0.65 |                  0.6  |                     0.65 |                       0.55 |                    0.75 |              0.612505 |
| CROSS_finetuned  |    0.841371 |    0.862878 |              0.875    |                0.605625 |            0.668925 |         0.681142 |        0.677825  |                 0.444723  |                         0.65 |                0.5  |                    0.42 |                            0.62 |                  0.6  |                     0.65 |                       0.55 |                    0.74 |              0.56966  |
| BigMHC_IM        |    0.671609 |    0.722291 |              0.875    |                0.333705 |            0.373655 |         0.672325 |        0.0835959 |                 0.0389681 |                         0.55 |                0.25 |                    0.35 |                            0.42 |                  0.6  |                     0.65 |                       0.55 |                    0.78 |              0.340321 |
| BigMHC_EL        |    0.647727 |    0.71921  |              0.729167 |                0.30759  |            0.306236 |         0.865599 |        0         |                 0         |                         0.55 |                0.25 |                    0.35 |                            0.42 |                  0.6  |                     0.65 |                       0.55 |                    0.78 |              0.3211   |

## Paperability primitives

| primitive                    | category               | gene                             | rl_action                       | why_it_matters                                                           | required_evidence                                                             |
|:-----------------------------|:-----------------------|:---------------------------------|:--------------------------------|:-------------------------------------------------------------------------|:------------------------------------------------------------------------------|
| novel_method_thesis          | paper_claim            | algorithm_contribution_statement | increase_method_distinctiveness | Turns a performance stack into a methods-paper claim.                    | method graph, genome definition, baseline contrast                            |
| sota_comparator_panel        | benchmark              | external_comparator_set          | add_or_update_public_baseline   | Prevents the result from looking like an internal leaderboard.           | BigMHC, PRIME, NetMHCpan, DeepImmuno, ImmunoStruct when available             |
| leakage_resistant_validation | reviewer_defense       | split_policy                     | switch_to_stricter_holdout      | Blocks the most common biomedical AI rejection risk.                     | source, patient, HLA, sequence-cluster and low-leakage splits                 |
| ablation_matrix              | method_evidence        | primitive_on_off_controls        | schedule_ablation               | Shows which part of the algorithm actually contributes.                  | KG-only, GA-only, RL-policy, guardrail, modality, fusion ablations            |
| prospective_freeze_protocol  | validation             | locked_threshold_and_score       | freeze_candidate_for_wetlab     | Converts retrospective SOTA signal into defensible prospective evidence. | hash/versioned score table, preregistered threshold, no post-hoc tuning       |
| hit_fail_interpreter         | translational_evidence | failure_mode_taxonomy            | allocate_failure_mode_assay     | Makes the wetlab result useful even when a candidate fails.              | candidate call table, endpoint unlock table, decomposition by assay arm       |
| figure_native_outputs        | paper_asset            | auto_figure_bundle               | prefer_plot_ready_artifact      | Forces the algorithm to emit artifacts that reviewers can inspect.       | benchmark plot, KG plot, reward plot, top-k plot, failure-mode plot           |
| claim_safety_layer           | reviewer_defense       | claim_boundary_cap               | add_claim_guardrail             | Lets strong results be shown without overselling clinical proof.         | claim tiers, blockers, caveats, prospective-validation status                 |
| negative_control_lane        | validation             | negative_control_design          | add_specificity_control         | Distinguishes real biological signal from dataset artifacts.             | decoy peptides, scrambled controls, unrelated target controls or null modules |
| decision_matrix              | paper_asset            | editor_reviewer_decision_table   | summarize_claim_unlock          | Compresses complex algorithm output into publishable decision logic.     | algorithm, evidence, caveat, disposition, next validation                     |

## Manuscript unlock map

| asset      | name                       | purpose                                                       | source_artifact                        |
|:-----------|:---------------------------|:--------------------------------------------------------------|:---------------------------------------|
| Figure 1   | Concept and method graph   | Show literature-to-primitive-to-GA/RL architecture            | method_kg_nodes/edges + loop schematic |
| Figure 2   | Benchmark lift             | Show KG-GA vs BigMHC/CROSS-Neo comparators                    | kg_ga_benchmark.tsv                    |
| Figure 3   | Validation-like robustness | Show source/low-leakage/generalization checks                 | kg_ga_source_benchmark.tsv             |
| Figure 4   | Guardrail ablation         | Show proxy/dominance/leakage guard effect                     | planned ablation matrix                |
| Figure 5   | 96-well interpreter        | Show actual hit/fail decomposition after wetlab               | v6 candidate_calls + endpoint_calls    |
| Table 1    | Method primitive library   | List public methods and reusable primitives                   | method_primitives.tsv                  |
| Table 2    | GA/RL genome and actions   | Define algorithm discovery search space                       | rl_ga_action_space.tsv                 |
| Table 3    | Reviewer risk register     | Leakage/proxy/overfit/claim-safety defenses                   | guardrails.tsv + paperability guards   |
| Supplement | Ablation dossier           | KG-only, GA-only, RL-policy, modality and guardrail ablations | future locked run                      |

## Claim boundary

This is a manuscript-readiness optimization scaffold. It should guide algorithm search and evidence packaging, but it must not generate protected manuscript prose or claim prospective superiority before locked wetlab/external validation.
