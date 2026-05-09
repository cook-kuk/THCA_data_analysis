# CLEAN-NeoBench Challenge Pack KR

## 한 줄 결론

이제 다음 할 일이 명확하다. aggregate leaderboard가 아니라 **실패하는 분포를 고정 challenge set으로 만들고**, 그 축에서 다시 검증해야 한다.

## Challenge axis 요약

| challenge_axis                       |   n_unique_candidates |   positive_prevalence |   mean_contextual_bma_score |   mean_confidence | recommended_split_contract                          |
|:-------------------------------------|----------------------:|----------------------:|----------------------------:|------------------:|:----------------------------------------------------|
| missed_positive_rescue_watchlist     |                   120 |              1        |                   0.0787338 |          0.149917 | HLA_heldout + source_heldout                        |
| korean_hla_focus_stress              |                   120 |              0.95     |                   0.774668  |          0.387933 | korean_hla_focus                                    |
| external_holdout_fragility           |                   120 |              0.791667 |                   0.165663  |          0.135854 | source_heldout_external                             |
| patient_gate_metadata_blocker        |                   120 |              0.966667 |                   0.781021  |          0.408078 | patient_gated_PAAD_THCA_demo                        |
| rare_hla_support_gap                 |                   120 |              0.991667 |                   0.703143  |          0.380167 | HLA_heldout + Korean_HLA_focus_if_applicable        |
| high_score_claim_blocked             |                   120 |              0.858333 |                   0.799003  |          0.38     | exact_phla_holdout + near_peptide_cluster_holdout   |
| low_prevalence_false_positive_stress |                   120 |              0        |                   0.308646  |          0.425333 | low_prevalence_heldout                              |
| public_internal_disagreement         |                    43 |              0.767442 |                   0.537812  |          0.447029 | public_overlap_audit + clean_internal_only_ablation |
| claim_safe_priority_review           |                     2 |              1        |                   0.669062  |          0.752358 | exact_phla + source_heldout + HLA_heldout           |

## 다음 실험

| priority   | experiment                               | why                                                                                  | success_metric                                                              |
|:-----------|:-----------------------------------------|:-------------------------------------------------------------------------------------|:----------------------------------------------------------------------------|
| P0         | public_training_corpus_row_overlap_audit | public pretrained tools remain caveated until row-level overlap is known             | zero unresolved public-training overlap for clean comparator status         |
| P0         | exact_near_overlap_lockdown              | high scores are currently claim-blocked by exact/near overlap risk                   | top-k performance after exact and near peptide-HLA removal                  |
| P1         | low_prevalence_topk_stress               | TESLA-like settings punish false positives more than aggregate AUPRC shows           | top10/top20 precision and false-positive pressure in low-prevalence sources |
| P1         | rare_hla_and_korean_hla_calibration      | underrepresented alleles drive both missed positives and unstable confidence         | allele-specific calibration ECE, AUPRC, positive rank                       |
| P1         | contextual_bma_clean_internal_ablation   | public tool support cannot be a clean feature until overlap audit is complete        | delta between all-expert contextual BMA and clean-contextual BMA            |
| P2         | PAAD_THCA_patient_gate_live_demo         | patient-gated score cannot be interpreted without disease/presentation/safety fields | metadata completion and gate-driven rank stability                          |
| P2         | MHC_II_separate_benchmark                | Class I and Class II must not be pooled as a single predictor                        | separate Class-II contracts and leaderboard                                 |

## 실용 판단

- high score라도 leakage high면 claim-blocked.
- low-prevalence source에서는 top-k false-positive stress가 우선.
- rare/Korean HLA는 allele-specific calibration이 필요.
- public pretrained agreement는 audit 전까지 caveated support다.
- PAAD/THCA patient gate는 실제 patient metadata 전까지 demo-only다.
