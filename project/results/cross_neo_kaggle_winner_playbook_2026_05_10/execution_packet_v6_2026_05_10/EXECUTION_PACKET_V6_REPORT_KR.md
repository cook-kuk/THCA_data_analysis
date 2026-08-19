# CROSS-Neo execution packet v6

Generated: 2026-05-10T13:51:53

## What changed

v6 turns the v5 preregistration package into a concrete execution handoff: order manifest, 96-well execution map, blank candidate/well result sheets, frozen decision worksheet, and an automatic endpoint interpreter.

## Headline metrics

- Candidates: 24
- 96-well entries: 96
- Reagent/audit order lines: 114
- Endpoint arms: 6
- Confirmatory family null-risk sum: 0.149
- Dry-run pending endpoints: 6

## Top execution priorities

| plate_v4_arm       | candidate_id   | peptide     | hla_allele_4digit   | order_batch             |   order_priority_score | claim_layer_after_unlock                      |
|:-------------------|:---------------|:------------|:--------------------|:------------------------|-----------------------:|:----------------------------------------------|
| B_mechanism_TCR_MD | CNV0_01132     | GADGVGKSAL  | HLA-C*08:02         | batch_1_claim_unlock    |                  0.737 | orthogonal TCR/MD mechanistic support         |
| B_mechanism_TCR_MD | CNV0_01151     | HMTEVVRHC   | HLA-A*02:01         | batch_1_claim_unlock    |                  0.718 | orthogonal TCR/MD mechanistic support         |
| A_clean_discovery  | CNV0_02708     | ALDPLLLRI   | HLA-A*02:01         | batch_1_claim_unlock    |                  0.612 | wetlab-supported clean antigen discovery lane |
| A_clean_discovery  | CNV0_02709     | ALPVALPSL   | HLA-A*02:01         | batch_1_claim_unlock    |                  0.602 | wetlab-supported clean antigen discovery lane |
| A_clean_discovery  | CNV0_02460     | GIVEGLITTV  | HLA-A*02:01         | batch_1_claim_unlock    |                  0.593 | wetlab-supported clean antigen discovery lane |
| A_clean_discovery  | CNV0_02448     | ILDKVLVHL   | HLA-A*02:01         | batch_1_claim_unlock    |                  0.592 | wetlab-supported clean antigen discovery lane |
| A_clean_discovery  | CNV0_02450     | KELEGILLL   | HLA-B*44:03         | batch_1_claim_unlock    |                  0.579 | wetlab-supported clean antigen discovery lane |
| A_clean_discovery  | CNV0_02452     | LLDGFLATV   | HLA-A*02:01         | batch_1_claim_unlock    |                  0.578 | wetlab-supported clean antigen discovery lane |
| C_label_rescue     | CNV0_00682     | NYRREPPHL   | HLA-C*07:02         | batch_2_risk_resolution |                  0.575 | label-noise or assay-mismatch evidence        |
| C_label_rescue     | CNV0_00711     | RREPPHLARNF | HLA-C*07:02         | batch_2_risk_resolution |                  0.574 | label-noise or assay-mismatch evidence        |

## Frozen endpoints

| plate_v4_arm          |   n_candidates |   threshold_successes |   confirmatory_threshold_successes | main_text_claim_tier     |   confirmatory_false_unlock_risk_under_null |
|:----------------------|---------------:|----------------------:|-----------------------------------:|:-------------------------|--------------------------------------------:|
| B_mechanism_TCR_MD    |              2 |                     1 |                                  2 | confirmatory_ready       |                                       0.003 |
| A_clean_discovery     |              6 |                     2 |                                  3 | confirmatory_ready       |                                       0.016 |
| C_label_rescue        |              4 |                     1 |                                  2 | confirmatory_ready       |                                       0.014 |
| E_positive_QC_control |              4 |                     3 |                                  3 | confirmatory_ready       |                                       0.027 |
| D_specificity_moat    |              4 |                     3 |                                  4 | confirmatory_with_caveat |                                       0.062 |
| F_model_boundary      |              4 |                     2 |                                  3 | confirmatory_ready       |                                       0.027 |

## Interpreter command

```bash
python scripts/interpret_cross_neo_assay_results_v6.py \
  --results /home/seungho/personal/THCA_data_analysis/project/results/cross_neo_kaggle_winner_playbook_2026_05_10/execution_packet_v6_2026_05_10/candidate_result_entry_v6.tsv \
  --endpoint-plan /home/seungho/personal/THCA_data_analysis/project/results/cross_neo_kaggle_winner_playbook_2026_05_10/preregistered_impact_v5_2026_05_10/preregistered_endpoint_plan_v5.tsv \
  --out-dir /home/seungho/personal/THCA_data_analysis/project/results/cross_neo_kaggle_winner_playbook_2026_05_10/execution_packet_v6_2026_05_10/interpreted_results_v6
```

## Claim boundary

The packet is an assay execution and interpretation scaffold. Positive controls remain QC only, source-overlap candidates remain mechanistic/support only, and patient-specific clinical selection remains locked.
