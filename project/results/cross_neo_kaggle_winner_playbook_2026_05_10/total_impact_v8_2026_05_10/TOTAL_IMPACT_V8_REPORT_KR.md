# CROSS-Neo total impact v8

Generated: 2026-05-10T21:25:17

## Headline

v8 combines the computational recovery, known-answer 96-well response board, preregistered assay endpoints, and execution packet into a single reviewer-safe impact dossier.

## Main numbers

- Known-answer top96: 91/96 (94.8%)
- Known-answer top10: 10/10
- Score recovery: 1.93x vs failed method-matrix stacker
- Clean-CV mean AUPRC/AUROC: 0.623 / 0.909
- Confirmatory null-risk sum: 0.149
- Execution packet: 24 candidates, 96 wells, 114 order/audit lines
- False-positive wells isolated: 5

## Impact stack

| impact_layer             | headline                                          | primary_metric             |   value | display       | secondary_metric                               | claim_boundary                                                     |
|:-------------------------|:--------------------------------------------------|:---------------------------|--------:|:--------------|:-----------------------------------------------|:-------------------------------------------------------------------|
| model_recovery           | score stack recovered clean-CV performance        | mean_clean_oof_auprc       |   0.623 | 0.623         | AUROC 0.909; recovery 1.93x                    | computational prioritization, not immunogenicity proof             |
| known_answer_response    | top-ranked 96-well board is visually high-yield   | top96_known_positive_hits  |  91.000 | 91/96         | precision 94.8%; top10 10/10                   | retrospective known labels, not prospective wetlab                 |
| positive_control_gallery | positive responder plate can be filled completely | known_positive_gallery     |  96.000 | 96/96         | positive-control visualization by construction | demo and control-board design only                                 |
| wetlab_translation       | candidate plan is assay-mapped                    | plate_candidates           |  24.000 | 24 candidates | 96 wells; mean unlock score 0.548              | execution plan only until assay data exists                        |
| preregistration          | endpoint rules are frozen before assay readout    | confirmatory_null_risk_sum |   0.149 | 0.149         | 6 endpoints; confirmatory mean power 0.591     | confirmatory-ready assay statistics, not pivotal clinical evidence |
| execution_packet         | handoff can be run without redesign               | order_audit_lines          | 114.000 | 114 lines     | 24 candidates; 96 wells; smoke 6/6             | lab-facing packet and interpreter only                             |
| failure_audit            | visible failures are already isolated             | false_positive_wells       |   5.000 | 5 wells       | all failure wells exported for action/audit    | failure analysis strengthens transparency                          |

## Evidence matrix

| claim_or_asset                       |   computational_score |   known_answer_visual |   preregistered_endpoint |   execution_ready |   prospective_wetlab | current_disposition                | safe_sentence                                                                       |
|:-------------------------------------|----------------------:|----------------------:|-------------------------:|------------------:|---------------------:|:-----------------------------------|:------------------------------------------------------------------------------------|
| general computational prioritization |                     1 |                     1 |                        1 |                 1 |                    0 | strong demo / prioritization claim | CROSS-Neo is a high-yield prioritization and assay-triage engine.                   |
| known-responder retrieval            |                     1 |                     1 |                        0 |                 0 |                    0 | visual demo only                   | Known-answer top-96 retrieval is 91/96 and should be labeled retrospective.         |
| clean antigen discovery lane         |                     1 |                     1 |                        1 |                 1 |                    0 | plate-ready, claim locked          | Clean discovery is ready for mutant pMHC + WT/decoy assay, but remains unvalidated. |
| TCR/MD mechanism support             |                     1 |                     1 |                        1 |                 1 |                    0 | orthogonal support track           | TCR/MD can support mechanism only after the predeclared readout passes.             |
| clinical vaccine selection           |                     0 |                     0 |                        0 |                 0 |                    0 | locked / not claimed               | No current result supports clinical vaccine selection.                              |

## Failure-action audit

| well   |   plate_rank | candidate_id   | peptide     | hla_allele_4digit   | source_name   | leakage_risk_level   | impact_lane        |   stress_guarded_discovery_score | impact_use                            | recommended_v8_action                                                         |
|:-------|-------------:|:---------------|:------------|:--------------------|:--------------|:---------------------|:-------------------|---------------------------------:|:--------------------------------------|:------------------------------------------------------------------------------|
| H2     |           16 | CNV0_00027     | GNNDVKEDP   | HLA-C*07:02         | CEDAR         | high                 | LABEL_NOISE_RESCUE |                            0.721 | label-noise/provenance audit showcase | keep as red well in visual; route to provenance audit before any rescue claim |
| E3     |           21 | CNV0_00006     | VKEDPKWEF   | HLA-C*07:02         | CEDAR         | high                 | LABEL_NOISE_RESCUE |                            0.716 | label-noise/provenance audit showcase | keep as red well in visual; route to provenance audit before any rescue claim |
| E7     |           53 | CNV0_00003     | LYGLLLEML   | HLA-A*02:01         | CEDAR         | high                 | LABEL_NOISE_RESCUE |                            0.694 | label-noise/provenance audit showcase | keep as red well in visual; route to provenance audit before any rescue claim |
| H7     |           56 | CNV0_00711     | RREPPHLARNF | HLA-C*07:02         | CEDAR         | high                 | LABEL_NOISE_RESCUE |                            0.693 | label-noise/provenance audit showcase | keep as red well in visual; route to provenance audit before any rescue claim |
| E9     |           69 | CNV0_00064     | LLAGLVSLL   | HLA-A*02:01         | CEDAR         | high                 | LABEL_NOISE_RESCUE |                            0.686 | label-noise/provenance audit showcase | keep as red well in visual; route to provenance audit before any rescue claim |

## Boundary

This dossier is intentionally split into retrospective known-answer evidence, preregistered assay design, and execution readiness. It should not be used as prospective wetlab validation or clinical vaccine-selection evidence.
