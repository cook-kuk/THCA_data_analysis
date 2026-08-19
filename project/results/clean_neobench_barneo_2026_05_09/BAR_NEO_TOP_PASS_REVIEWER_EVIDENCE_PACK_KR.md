# BAR-Neo Top-Pass Reviewer Evidence Pack KR

## 한 줄 결론

킬 오딧을 통과한 high-impact 후보만 따로 분리했다. 현재 자동 기준으로 바로 밀 수 있는 후보는 3개이며, 모두 ITSNdb_main / HLA-A*02:01 쪽에서 low leakage, no exact/near overlap, no same-source negative pressure, low public/QK/fallback dependence를 만족한다.

## Top-pass 후보 카드

| candidate_id   |   stress_guarded_rank_global |   stress_guarded_final_review_score | source_name   | hla_allele_4digit   | peptide   | exact_overlap_clean   | near_neighbor_clean   | decoy_pressure_clean   | public_dependency_clean   | fallback_dependency_clean   | top_method            | reviewer_evidence_summary                                                                           |
|:---------------|-----------------------------:|------------------------------------:|:--------------|:--------------------|:----------|:----------------------|:----------------------|:-----------------------|:--------------------------|:----------------------------|:----------------------|:----------------------------------------------------------------------------------------------------|
| CNV0_02407     |                            1 |                            0.551702 | ITSNdb_main   | HLA-A*02:01         | KLMNIQQKL | True                  | True                  | True                   | True                      | True                        | Stack_mean_E1         | low leakage, no exact/near overlap, no same-source decoy pressure, clean internal support dominates |
| CNV0_02504     |                            2 |                            0.541321 | ITSNdb_main   | HLA-A*02:01         | LLVDLAEEL | True                  | True                  | True                   | True                      | True                        | Stack_LR_inmaster_E3a | low leakage, no exact/near overlap, no same-source decoy pressure, clean internal support dominates |
| CNV0_02410     |                            3 |                            0.511978 | ITSNdb_main   | HLA-A*02:01         | MLGEQLFPL | True                  | True                  | True                   | True                      | True                        | Stack_LR_inmaster_E3a | low leakage, no exact/near overlap, no same-source decoy pressure, clean internal support dominates |

## 후보별 method support

| candidate_id   |   support_rank_within_candidate | method_name                | method_role        |   score_for_fusion |   candidate_method_weight |   candidate_method_contribution |
|:---------------|--------------------------------:|:---------------------------|:-------------------|-------------------:|--------------------------:|--------------------------------:|
| CNV0_02407     |                               1 | Stack_mean_E1              | internal_candidate |           0.692308 |                  0.610431 |                        0.422606 |
| CNV0_02407     |                               2 | W7A_full                   | internal_candidate |           0.615385 |                  0.649666 |                        0.399794 |
| CNV0_02407     |                               3 | Stack_LR_inmaster_E3a      | internal_candidate |           0.560976 |                  0.619192 |                        0.347352 |
| CNV0_02407     |                               4 | Wave8_TCR_SelfSim_no_exact | internal_candidate |           0.888889 |                  0.366647 |                        0.325908 |
| CNV0_02407     |                               5 | W7B_stacked                | internal_candidate |           0.459016 |                  0.649494 |                        0.298128 |
| CNV0_02407     |                               6 | Wave8_TCR_motif_only       | internal_candidate |           0.818182 |                  0.351448 |                        0.287549 |
| CNV0_02407     |                               7 | kNN                        | internal_candidate |           0.666667 |                  0.41262  |                        0.27508  |
| CNV0_02410     |                               1 | Stack_LR_inmaster_E3a      | internal_candidate |           0.636364 |                  0.619192 |                        0.394032 |
| CNV0_02410     |                               2 | W7A_full                   | internal_candidate |           0.545455 |                  0.649666 |                        0.354363 |
| CNV0_02410     |                               3 | Stack_mean_E1              | internal_candidate |           0.516667 |                  0.610431 |                        0.315389 |
| CNV0_02410     |                               4 | W7B_stacked                | internal_candidate |           0.45     |                  0.649494 |                        0.292272 |
| CNV0_02410     |                               5 | Wave8_TCR_only             | internal_candidate |           0.866667 |                  0.326857 |                        0.283276 |
| CNV0_02410     |                               6 | Wave8_TCR_motif_only       | internal_candidate |           0.805556 |                  0.351448 |                        0.283111 |
| CNV0_02410     |                               7 | MultiTask_A_only           | internal_candidate |           0.625    |                  0.393212 |                        0.245758 |
| CNV0_02504     |                               1 | Stack_LR_inmaster_E3a      | internal_candidate |           0.613636 |                  0.619192 |                        0.379959 |
| CNV0_02504     |                               2 | W7A_full                   | internal_candidate |           0.540423 |                  0.649666 |                        0.351095 |
| CNV0_02504     |                               3 | Wave8_TCR_SelfSim_no_exact | internal_candidate |           0.9      |                  0.366647 |                        0.329982 |
| CNV0_02504     |                               4 | Stack_mean_E1              | internal_candidate |           0.513499 |                  0.610431 |                        0.313456 |
| CNV0_02504     |                               5 | Wave8_TCR_motif_only       | internal_candidate |           0.833333 |                  0.351448 |                        0.292874 |
| CNV0_02504     |                               6 | MultiTask_A_only           | internal_candidate |           0.702703 |                  0.393212 |                        0.276311 |
| CNV0_02504     |                               7 | MultiTask_AD               | internal_candidate |           0.686275 |                  0.381783 |                        0.262008 |

## Lab / manual-review handoff queue

| candidate_id   | handoff_tier                  |   stress_guarded_rank_global | source_name   | hla_allele_4digit   | peptide     | reviewer_kill_disposition          | handoff_action                                                            | must_verify_before_claim                                                                        |
|:---------------|:------------------------------|-----------------------------:|:--------------|:--------------------|:------------|:-----------------------------------|:--------------------------------------------------------------------------|:------------------------------------------------------------------------------------------------|
| CNV0_02407     | T1_clean_manual_review_lead   |                            1 | ITSNdb_main   | HLA-A*02:01         | KLMNIQQKL   | passes_current_reviewer_kill_audit | advance to manual biological plausibility and assay-design review         | source provenance; WT/gene/mutation; expression/clonality; patient gate; safety context         |
| CNV0_02504     | T1_clean_manual_review_lead   |                            2 | ITSNdb_main   | HLA-A*02:01         | LLVDLAEEL   | passes_current_reviewer_kill_audit | advance to manual biological plausibility and assay-design review         | source provenance; WT/gene/mutation; expression/clonality; patient gate; safety context         |
| CNV0_02410     | T1_clean_manual_review_lead   |                            3 | ITSNdb_main   | HLA-A*02:01         | MLGEQLFPL   | passes_current_reviewer_kill_audit | advance to manual biological plausibility and assay-design review         | source provenance; WT/gene/mutation; expression/clonality; patient gate; safety context         |
| CNV0_02480     | T2_resolve_before_claim       |                            5 | ITSNdb_main   | HLA-A*01:01         | ILDTAGKEEY  | manual_audit_required              | resolve near-peptide/source/protein-window caveat before any lead claim   | near-neighbor lineage; source/protein-window overlap; then standard provenance and patient gate |
| CNV0_00487     | T3_do_not_use_for_clean_claim |                           17 | CEDAR         | HLA-B*40:01         | QELNELSAISL | blocked_from_clean_claim           | exclude from clean lead list; keep only as leakage/stress-control example | exact overlap and high-leakage flags; do not advance as clean lead                              |
| CNV0_00245     | T3_do_not_use_for_clean_claim |                           18 | CEDAR         | HLA-A*11:01         | ATSPHLESLLK | blocked_from_clean_claim           | exclude from clean lead list; keep only as leakage/stress-control example | exact overlap and high-leakage flags; do not advance as clean lead                              |
| CNV0_00231     | T3_do_not_use_for_clean_claim |                           19 | CEDAR         | HLA-A*02:01         | AMFGKLMTI   | blocked_from_clean_claim           | exclude from clean lead list; keep only as leakage/stress-control example | exact overlap and high-leakage flags; do not advance as clean lead                              |
| CNV0_00156     | T3_do_not_use_for_clean_claim |                           20 | CEDAR         | HLA-A*02:01         | ALPEVLAVIQV | blocked_from_clean_claim           | exclude from clean lead list; keep only as leakage/stress-control example | exact overlap and high-leakage flags; do not advance as clean lead                              |
| CNV0_00028     | T3_do_not_use_for_clean_claim |                           21 | CEDAR         | HLA-A*24:02         | LFMNVQFLF   | blocked_from_clean_claim           | exclude from clean lead list; keep only as leakage/stress-control example | exact overlap and high-leakage flags; do not advance as clean lead                              |
| CNV0_00295     | T3_do_not_use_for_clean_claim |                           22 | CEDAR         | HLA-A*11:01         | KSFKLSGFSFK | blocked_from_clean_claim           | exclude from clean lead list; keep only as leakage/stress-control example | exact overlap and high-leakage flags; do not advance as clean lead                              |
| CNV0_00751     | T3_do_not_use_for_clean_claim |                           23 | CEDAR         | HLA-A*11:01         | GSFPENLRHLK | blocked_from_clean_claim           | exclude from clean lead list; keep only as leakage/stress-control example | exact overlap and high-leakage flags; do not advance as clean lead                              |
| CNV0_00061     | T3_do_not_use_for_clean_claim |                           24 | CEDAR         | HLA-A*02:01         | FLALIICNA   | blocked_from_clean_claim           | exclude from clean lead list; keep only as leakage/stress-control example | exact overlap and high-leakage flags; do not advance as clean lead                              |

## 다음 액션

1. peptide/HLA/source provenance 수동 확인.
2. wild-type peptide, gene/mutation, expression/clonality metadata 연결.
3. PAAD/THCA patient-gate demo에 실제 disease timing과 immune/safety context를 붙여 triage view 생성.
4. public training-corpus row-level overlap audit 전까지 public tool support는 caveated로 유지.

## Claim boundary

manual-review research lead만 허용한다. clinical vaccine selection, new SOTA predictor, external validation proven, quantum advantage claim은 금지한다.
