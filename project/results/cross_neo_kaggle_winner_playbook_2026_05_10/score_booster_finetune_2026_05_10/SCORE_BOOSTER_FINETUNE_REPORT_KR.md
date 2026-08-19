# CROSS-Neo score booster fine-tune report

Generated: 2026-05-10T12:12:55

## Result

Method-matrix stacker가 source-heldout에서 낮게 나와서, 기존 BAR-Neo stress score 계열을 clean source/HLA OOF로 calibration/booster fine-tune했다. Best config는 `stress_only` + `logistic_C1.0`.

## Best clean-CV score

| feature_set   | model_name    | protocol           |   n_rows |   n_pos |   auprc |   auroc |   brier |   ece |   top10_precision |   top20_precision |
|:--------------|:--------------|:-------------------|---------:|--------:|--------:|--------:|--------:|------:|------------------:|------------------:|
| stress_only   | logistic_C1.0 | source_heldout_oof |      410 |      35 |   0.603 |   0.920 |   0.159 | 0.204 |             0.900 |             0.650 |
| stress_only   | logistic_C1.0 | hla_heldout_oof    |      410 |      35 |   0.643 |   0.898 |   0.115 | 0.185 |             0.800 |             0.700 |

## Reference scores

| feature_set                                 |   auprc |   auroc |   brier |   ece |   top10_precision |   top20_precision |
|:--------------------------------------------|--------:|--------:|--------:|------:|------------------:|------------------:|
| reference_stress_guarded_final_review_score |   0.658 |   0.937 |   0.078 | 0.169 |             0.900 |             0.700 |
| reference_bma_v2_discovery_score            |   0.620 |   0.892 |   0.158 | 0.328 |             0.800 |             0.650 |
| reference_stress_guarded_discovery_score    |   0.591 |   0.914 |   0.110 | 0.234 |             0.700 |             0.800 |
| reference_stress_guarded_clean_score        |   0.547 |   0.899 |   0.109 | 0.235 |             0.700 |             0.800 |
| reference_bma_v2_claim_safe_score           |   0.532 |   0.898 |   0.121 | 0.259 |             0.700 |             0.650 |

## Top experiment-priority scores

| candidate_id   | peptide     | hla_allele_4digit   |   label | source_name   | leakage_risk_level   |   finetuned_score_booster_prob |   finetuned_experiment_priority_score |   finetuned_claim_capped_score | finetuned_score_action                  | primary_claim_blocker   |
|:---------------|:------------|:--------------------|--------:|:--------------|:---------------------|-------------------------------:|--------------------------------------:|-------------------------------:|:----------------------------------------|:------------------------|
| CNV0_00245     | ATSPHLESLLK | HLA-A*11:01         |       1 | CEDAR         | high                 |                          0.974 |                                 0.884 |                          0.250 | score_high_but_overlap_blocked          | patient context gate    |
| CNV0_00487     | QELNELSAISL | HLA-B*40:01         |       1 | CEDAR         | high                 |                          0.974 |                                 0.883 |                          0.250 | score_high_but_overlap_blocked          | patient context gate    |
| CNV0_02452     | LLDGFLATV   | HLA-A*02:01         |       1 | ITSNdb_main   | low                  |                          0.990 |                                 0.882 |                          0.420 | experiment_priority_with_claim_boundary | patient context gate    |
| CNV0_00231     | AMFGKLMTI   | HLA-A*02:01         |       1 | CEDAR         | high                 |                          0.973 |                                 0.882 |                          0.420 | score_high_but_overlap_blocked          | patient context gate    |
| CNV0_00156     | ALPEVLAVIQV | HLA-A*02:01         |       1 | CEDAR         | high                 |                          0.973 |                                 0.881 |                          0.420 | score_high_but_overlap_blocked          | patient context gate    |
| CNV0_00028     | LFMNVQFLF   | HLA-A*24:02         |       1 | CEDAR         | high                 |                          0.972 |                                 0.880 |                          0.420 | score_high_but_overlap_blocked          | patient context gate    |
| CNV0_00295     | KSFKLSGFSFK | HLA-A*11:01         |       1 | CEDAR         | high                 |                          0.970 |                                 0.879 |                          0.250 | score_high_but_overlap_blocked          | patient context gate    |
| CNV0_00751     | GSFPENLRHLK | HLA-A*11:01         |       1 | CEDAR         | high                 |                          0.970 |                                 0.878 |                          0.250 | score_high_but_overlap_blocked          | patient context gate    |
| CNV0_00061     | FLALIICNA   | HLA-A*02:01         |       1 | CEDAR         | high                 |                          0.969 |                                 0.876 |                          0.420 | score_high_but_overlap_blocked          | patient context gate    |
| CNV0_00742     | YWNEYGGGLLW | HLA-A*24:02         |       1 | CEDAR         | high                 |                          0.968 |                                 0.876 |                          0.420 | score_high_but_overlap_blocked          | patient context gate    |
| CNV0_02448     | ILDKVLVHL   | HLA-A*02:01         |       1 | ITSNdb_main   | low                  |                          0.988 |                                 0.876 |                          0.420 | experiment_priority_with_claim_boundary | patient context gate    |
| CNV0_00816     | FLSEVWNTHTL | HLA-A*02:01         |       1 | CEDAR         | high                 |                          0.967 |                                 0.874 |                          0.420 | score_high_but_overlap_blocked          | patient context gate    |

## Top claim-capped scores

| candidate_id   | peptide     | hla_allele_4digit   |   label | source_name   | leakage_risk_level   |   finetuned_score_booster_prob |   finetuned_claim_capped_score | finetuned_score_action                  | primary_claim_blocker   |
|:---------------|:------------|:--------------------|--------:|:--------------|:---------------------|-------------------------------:|-------------------------------:|:----------------------------------------|:------------------------|
| CNV0_01063     | RFLEYLPLRF  | HLA-A*24:02         |       1 | NEPdb         | high                 |                          0.420 |                          0.420 | score_support_only                      | patient context gate    |
| CNV0_01014     | SSYTGFANK   | HLA-A*11:01         |       1 | NEPdb         | high                 |                          0.439 |                          0.420 | score_support_only                      | patient context gate    |
| CNV0_00009     | SLSPALPGA   | HLA-A*02:01         |       1 | CEDAR         | high                 |                          0.960 |                          0.420 | score_high_but_overlap_blocked          | patient context gate    |
| CNV0_00050     | AMIPKDWPL   | HLA-A*02:01         |       1 | CEDAR         | high                 |                          0.962 |                          0.420 | score_high_but_overlap_blocked          | patient context gate    |
| CNV0_00541     | SLQAIQQLV   | HLA-A*02:01         |       1 | CEDAR         | high                 |                          0.962 |                          0.420 | score_high_but_overlap_blocked          | patient context gate    |
| CNV0_02450     | KELEGILLL   | HLA-B*44:03         |       1 | ITSNdb_main   | low                  |                          0.976 |                          0.420 | experiment_priority_with_claim_boundary | patient context gate    |
| CNV0_00064     | LLAGLVSLL   | HLA-A*02:01         |       0 | CEDAR         | high                 |                          0.963 |                          0.420 | score_high_but_overlap_blocked          | patient context gate    |
| CNV0_02397     | LPIQYEPVL   | HLA-B*35:03         |       0 | ITSNdb_main   | low                  |                          0.994 |                          0.420 | experiment_priority_with_claim_boundary | patient context gate    |
| CNV0_00540     | EVLETRVMER  | HLA-A*24:02         |       1 | CEDAR         | high                 |                          0.965 |                          0.420 | score_high_but_overlap_blocked          | patient context gate    |
| CNV0_00037     | MLAAPISGL   | HLA-A*02:01         |       1 | CEDAR         | high                 |                          0.963 |                          0.420 | score_high_but_overlap_blocked          | patient context gate    |
| CNV0_00088     | YAYDNFGVLGL | HLA-C*03:03         |       1 | CEDAR         | high                 |                          0.949 |                          0.420 | score_high_but_overlap_blocked          | patient context gate    |
| CNV0_02407     | KLMNIQQKL   | HLA-A*02:01         |       1 | ITSNdb_main   | low                  |                          0.999 |                          0.420 | experiment_priority_with_claim_boundary | patient context gate    |

## Feature importance

| feature                           |   importance |   abs_importance |
|:----------------------------------|-------------:|-----------------:|
| stress_guarded_final_review_score |        1.883 |            1.883 |

## Decision

1. Final score column: `finetuned_score_booster_prob`.
2. Experiment ranking column: `finetuned_experiment_priority_score`.
3. Claim-safe ranking column: `finetuned_claim_capped_score`, which keeps overlap/patient/source gates active.
4. Low-performing method-matrix fine-tune is not used as the main score; it remains a diagnostic table for cap/drop/recalibration.
