# GA/RL Rescue Quadrant KR

## 한 줄 결론

이 표는 GA/RL 결과를 세 갈래로 분리한다: `T1_survivor`, `recoverable_watchlist`, `hard_block_overlap`.
현재 회복 후보는 `recoverable_watchlist` 3개뿐이고, 나머지 blocked는 clean claim으로는 막혀 있다.

## Claim boundary

Allowed:

- leakage-aware rescue queue
- manual-review prioritization
- benchmark-adaptive reliability audit

Forbidden:

- clinical vaccine selection
- new SOTA predictor
- public clean comparator claim without audit

## Summary

| rescue_category       |   n_candidates |   n_pos |   mean_priority_score |   mean_rescue_priority_score |   mean_top_nn_similarity |
|:----------------------|---------------:|--------:|----------------------:|-----------------------------:|-------------------------:|
| T1_survivor           |              3 |       3 |              0.723071 |                     0.321024 |                 0.555556 |
| hard_block_overlap    |             69 |      69 |              0.25     |                     0.134647 |                 0.461413 |
| recoverable_watchlist |              3 |       3 |              0.563031 |                     0.274252 |                 0.503704 |

## Reason summary

| rescue_category       | rescue_reason_primary                            |   n_candidates |   mean_priority_score |   mean_top_nn_similarity |
|:----------------------|:-------------------------------------------------|---------------:|----------------------:|-------------------------:|
| T1_survivor           | passes BAR-Neo-NG T1 gate                        |              3 |              0.723071 |                 0.555556 |
| hard_block_overlap    | exact/near/public overlap or non-watchlist block |             69 |              0.25     |                 0.461413 |
| recoverable_watchlist | high GA score without exact/near/public overlap  |              3 |              0.563031 |                 0.503704 |

## Rescue queue

|   rescue_order | rescue_category       | candidate_id   | peptide         | hla_allele_4digit   | source_name   |   label |   ga_rl_score |   barneo_ng_score |   ga_rl_barneo_priority_score | ga_rl_barneo_decision                | exact_peptide_hla_train_overlap   | near_peptide_train_overlap   | public_tool_training_overlap_any   | split_source_heldout   | split_hla_heldout   | rescue_action                       | rescue_reason_primary                            | top_nn_candidate_id   |   top_nn_similarity |
|---------------:|:----------------------|:---------------|:----------------|:--------------------|:--------------|--------:|--------------:|------------------:|------------------------------:|:-------------------------------------|:----------------------------------|:-----------------------------|:-----------------------------------|:-----------------------|:--------------------|:------------------------------------|:-------------------------------------------------|:----------------------|--------------------:|
|              1 | T1_survivor           | CNV0_02407     | KLMNIQQKL       | HLA-A*02:01         | ITSNdb_main   |       1 |      0.740475 |         0.72      |                      0.729214 | GA_RL_hit_confirmed_T1_by_BAR_Neo_NG | False                             | False                        | False                              | ITSNdb_main            | HLA-A*02:01         | assay_design_handoff                | passes BAR-Neo-NG T1 gate                        | CNV0_01694            |            0.666667 |
|              2 | T1_survivor           | CNV0_02504     | LLVDLAEEL       | HLA-A*02:01         | ITSNdb_main   |       1 |      0.685061 |         0.72      |                      0.72     | GA_RL_hit_confirmed_T1_by_BAR_Neo_NG | False                             | False                        | False                              | ITSNdb_main            | HLA-A*02:01         | assay_design_handoff                | passes BAR-Neo-NG T1 gate                        | CNV0_00394            |            0.444444 |
|              3 | T1_survivor           | CNV0_02410     | MLGEQLFPL       | HLA-A*02:01         | ITSNdb_main   |       1 |      0.638571 |         0.72      |                      0.72     | GA_RL_hit_confirmed_T1_by_BAR_Neo_NG | False                             | False                        | False                              | ITSNdb_main            | HLA-A*02:01         | assay_design_handoff                | passes BAR-Neo-NG T1 gate                        | CNV0_00197            |            0.555556 |
|              4 | recoverable_watchlist | CNV0_02448     | ILDKVLVHL       | HLA-A*02:01         | ITSNdb_main   |       1 |      0.996333 |         0.499141  |                      0.62     | GA_RL_hit_watchlist_manual_review    | False                             | False                        | False                              | ITSNdb_main            | HLA-A*02:01         | manual_review_with_metadata_intake  | high GA score without exact/near/public overlap  | CNV0_01760            |            0.555556 |
|              5 | recoverable_watchlist | CNV0_02708     | ALDPLLLRI       | HLA-A*02:01         | ITSNdb_Val    |       1 |      0.709598 |         0.521087  |                      0.605917 | GA_RL_hit_watchlist_manual_review    | False                             | False                        | False                              | ITSNdb_Val             | HLA-A*02:01         | manual_review_with_metadata_intake  | high GA score without exact/near/public overlap  | CNV0_00209            |            0.555556 |
|              6 | recoverable_watchlist | CNV0_02405     | KMIGNHLWV       | HLA-A*02:01         | ITSNdb_main   |       1 |      0.460144 |         0.465655  |                      0.463175 | GA_RL_hit_watchlist_manual_review    | False                             | False                        | False                              | ITSNdb_main            | HLA-A*02:01         | manual_review_with_metadata_intake  | high GA score without exact/near/public overlap  | CNV0_00741            |            0.4      |
|              7 | hard_block_overlap    | CNV0_00688     | YYYGIKDLATVFF   | HLA-A*24:02         | CEDAR         |       1 |      1        |         0.122851  |                      0.25     | GA_RL_hit_blocked_from_clean_claim   | True                              | True                         | False                              | CEDAR                  | HLA-A*24:02         | do_not_promote_without_new_evidence | exact/near/public overlap or non-watchlist block | CNV0_00323            |            0.307692 |
|              8 | hard_block_overlap    | CNV0_01132     | GADGVGKSAL      | HLA-C*08:02         | NEPdb         |       1 |      1        |         0.0997747 |                      0.25     | GA_RL_hit_blocked_from_clean_claim   | True                              | True                         | False                              | NEPdb                  | HLA-C*08:02         | do_not_promote_without_new_evidence | exact/near/public overlap or non-watchlist block | CNV0_01092            |            1        |
|              9 | hard_block_overlap    | CNV0_00784     | ETVNALISDQKL    | HLA-A*68:02         | CEDAR         |       1 |      0.999559 |         0.111528  |                      0.25     | GA_RL_hit_blocked_from_clean_claim   | True                              | True                         | False                              | CEDAR                  | HLA-A*68:02         | do_not_promote_without_new_evidence | exact/near/public overlap or non-watchlist block | CNV0_00245            |            0.333333 |
|             10 | hard_block_overlap    | CNV0_00169     | AALEDTLAETEAR   | HLA-B*37:01         | CEDAR         |       1 |      0.9513   |         0.112326  |                      0.25     | GA_RL_hit_blocked_from_clean_claim   | True                              | True                         | False                              | CEDAR                  | HLA-B*37:01         | do_not_promote_without_new_evidence | exact/near/public overlap or non-watchlist block | CNV0_00120            |            0.307692 |
|             11 | hard_block_overlap    | CNV0_00804     | SRSYTSGPGSRISSS | HLA-B*27:09         | CEDAR         |       1 |      0.946788 |         0.11132   |                      0.25     | GA_RL_hit_blocked_from_clean_claim   | True                              | True                         | False                              | CEDAR                  | HLA-B*27:09         | do_not_promote_without_new_evidence | exact/near/public overlap or non-watchlist block | CNV0_00902            |            0.5625   |
|             12 | hard_block_overlap    | CNV0_00874     | FPSEYLSSHLEA    | HLA-B*54:01         | CEDAR         |       1 |      0.939424 |         0.111181  |                      0.25     | GA_RL_hit_blocked_from_clean_claim   | True                              | True                         | False                              | CEDAR                  | HLA-B*54:01         | do_not_promote_without_new_evidence | exact/near/public overlap or non-watchlist block | CNV0_00039            |            0.416667 |
|             13 | hard_block_overlap    | CNV0_00791     | EVIDKNSGGWWYV   | HLA-A*68:02         | CEDAR         |       1 |      0.91656  |         0.106931  |                      0.25     | GA_RL_hit_blocked_from_clean_claim   | True                              | True                         | False                              | CEDAR                  | HLA-A*68:02         | do_not_promote_without_new_evidence | exact/near/public overlap or non-watchlist block | CNV0_00651            |            0.384615 |
|             14 | hard_block_overlap    | CNV0_00821     | NADPASHEIW      | HLA-B*53:01         | CEDAR         |       1 |      0.911398 |         0.119311  |                      0.25     | GA_RL_hit_blocked_from_clean_claim   | True                              | True                         | False                              | CEDAR                  | HLA-B*53:01         | do_not_promote_without_new_evidence | exact/near/public overlap or non-watchlist block | CNV0_00379            |            0.4      |
|             15 | hard_block_overlap    | CNV0_00771     | YYGTGETFLYTF    | HLA-A*24:02         | CEDAR         |       1 |      0.910634 |         0.12171   |                      0.25     | GA_RL_hit_blocked_from_clean_claim   | True                              | True                         | False                              | CEDAR                  | HLA-A*24:02         | do_not_promote_without_new_evidence | exact/near/public overlap or non-watchlist block | CNV0_00133            |            0.333333 |
|             16 | hard_block_overlap    | CNV0_00245     | ATSPHLESLLK     | HLA-A*11:01         | CEDAR         |       1 |      0.908619 |         0.109555  |                      0.25     | GA_RL_hit_blocked_from_clean_claim   | True                              | True                         | False                              | CEDAR                  | HLA-A*11:01         | do_not_promote_without_new_evidence | exact/near/public overlap or non-watchlist block | CNV0_00841            |            0.454545 |
|             17 | hard_block_overlap    | CNV0_00120     | AEEEASAVSTAA    | HLA-B*45:01         | CEDAR         |       1 |      0.906665 |         0.111753  |                      0.25     | GA_RL_hit_blocked_from_clean_claim   | True                              | True                         | False                              | CEDAR                  | HLA-B*45:01         | do_not_promote_without_new_evidence | exact/near/public overlap or non-watchlist block | CNV0_00194            |            0.416667 |
|             18 | hard_block_overlap    | CNV0_00704     | MTEVISSLENANY   | HLA-A*01:01         | CEDAR         |       1 |      0.905715 |         0.120643  |                      0.25     | GA_RL_hit_blocked_from_clean_claim   | True                              | True                         | False                              | CEDAR                  | HLA-A*01:01         | do_not_promote_without_new_evidence | exact/near/public overlap or non-watchlist block | CNV0_00475            |            0.307692 |
|             19 | hard_block_overlap    | CNV0_00779     | RTLMENQHW       | HLA-B*58:01         | CEDAR         |       1 |      0.90527  |         0.118176  |                      0.25     | GA_RL_hit_blocked_from_clean_claim   | True                              | True                         | False                              | CEDAR                  | HLA-B*58:01         | do_not_promote_without_new_evidence | exact/near/public overlap or non-watchlist block | CNV0_02037            |            0.4      |
|             20 | hard_block_overlap    | CNV0_00812     | AEIDVIFKDFVNKY  | HLA-B*44:03         | CEDAR         |       1 |      0.904865 |         0.108313  |                      0.25     | GA_RL_hit_blocked_from_clean_claim   | True                              | True                         | False                              | CEDAR                  | HLA-B*44:03         | do_not_promote_without_new_evidence | exact/near/public overlap or non-watchlist block | CNV0_00268            |            0.285714 |

## Interpretation

- `T1_survivor`: already assay-design ready.
- `recoverable_watchlist`: high-GA rows without exact/near/public overlap; these are the only realistic rescue candidates.
- `hard_block_overlap`: blocked by overlap or non-watchlist behavior; do not promote without new evidence.

## Output files

- `/home/seungho/personal/THCA_data_analysis/project/results/clean_neobench_barneo_2026_05_09/ga_rl_rescue_quadrant.tsv`
- `/home/seungho/personal/THCA_data_analysis/project/results/clean_neobench_barneo_2026_05_09/ga_rl_rescue_quadrant_summary.tsv`
- `/home/seungho/personal/THCA_data_analysis/project/results/clean_neobench_barneo_2026_05_09/ga_rl_rescue_quadrant_reason_summary.tsv`
- `/home/seungho/personal/THCA_data_analysis/project/results/clean_neobench_barneo_2026_05_09/figures/fig_ga_rl_rescue_quadrant.png`
