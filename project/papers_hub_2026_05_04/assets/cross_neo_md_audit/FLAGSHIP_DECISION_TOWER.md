# CROSS-Neo flagship decision tower

## Executive verdict

This is the top-most decision layer. It turns the two flagship leads into an explicit go/hold/escalate structure with branch conditions.

## Key numbers

- candidates: 13
- paired TCR evidence: 6
- physics flagships: 2
- preclinical-ready: 2
- real assay summaries: 0

## Decision tower

| lead                 | peptide    | hla_4digit   |   top_md_score |   preclinical_readiness |   physics_priority_score | preclinical_next_action                                               | physics_tier        | launch_sequence                                                                                          | endpoint_step                                                        | command_priority                                                                    | next_72h_goal                                                 | claim_boundary                                  | command_branch                                                   | current_action                                              | abort_condition                                                  | next_wetlab_step                                        | command_state   | go_no_go_state   | escalation_branch                                                           | 24h_action                                               | 72h_action                                                       | decision_owner   |
|:---------------------|:-----------|:-------------|---------------:|------------------------:|-------------------------:|:----------------------------------------------------------------------|:--------------------|:---------------------------------------------------------------------------------------------------------|:---------------------------------------------------------------------|:------------------------------------------------------------------------------------|:--------------------------------------------------------------|:------------------------------------------------|:-----------------------------------------------------------------|:------------------------------------------------------------|:-----------------------------------------------------------------|:--------------------------------------------------------|:----------------|:-----------------|:----------------------------------------------------------------------------|:---------------------------------------------------------|:-----------------------------------------------------------------|:-----------------|
| TP53_R175H_HLA_A0201 | HMTEVVRHC  | HLA-A*02:01  |       0.821789 |                0.935447 |                 0.839191 | launch WT/decoy-controlled HLA stability + multimer/activation assays | P0_FLAGSHIP_PHYSICS | endpoint energy now -> PMF if endpoint supports mutant specificity -> FEP only after WT mapping is clean | compare mutant vs WT vs decoy endpoint energies over ensemble frames | run endpoint energy now; keep PMF/FEP held unless the endpoint gate stays favorable | WT/decoy-controlled HLA stability + multimer/activation setup | prioritized execution, not immunogenicity proof | endpoint energy first; PMF/FEP held unless endpoint is favorable | prepare WT/decoy-controlled HLA stability and binding setup | WT or decoy does not stay separated from mutant on endpoint gate | multimer/activation only after presentation gate passes | RUN_ENDPOINT    | GO_ENDPOINT      | PMF/FEP hold until endpoint energy and WT/decoy separation remain favorable | prepare presentation gate controls and QC if GO_ENDPOINT | run endpoint energy; only then unlock multimer/activation branch | execution lead   |
| KRAS_G12D_HLA_C0802  | GADGVGKSAL | HLA-C*08:02  |       0.577335 |                0.796652 |                 0.656129 | launch WT/decoy-controlled HLA stability + multimer/activation assays | P0_FLAGSHIP_PHYSICS | endpoint energy now -> PMF only if TCR-pMHC model stays stable                                           | compare mutant vs WT vs decoy endpoint energies over ensemble frames | run endpoint energy now; keep PMF/FEP held unless the endpoint gate stays favorable | WT/decoy-controlled HLA stability + multimer/activation setup | prioritized execution, not immunogenicity proof | endpoint energy first; PMF only if TCR-pMHC stays stable         | prepare WT/decoy-controlled HLA stability and binding setup | pairing or TCR model becomes unstable                            | multimer/activation only after presentation gate passes | RUN_ENDPOINT    | HOLD_ENDPOINT    | do not escalate beyond endpoint physics                                     | recheck model/assay compatibility and hold physics spend | hold until a cleaner assay or structure signal appears           | execution lead   |

## Required refresh chain

```bash
python project/scripts/cross_neo_md/44_build_executive_impact_console.py
python project/scripts/cross_neo_md/43_build_decision_atlas.py
python project/scripts/cross_neo_md/42_build_integrated_decision_matrix.py
python project/scripts/cross_neo_md/41_build_physics_case_study_board.py
python project/scripts/cross_neo_md/38_build_physics_gate_package.py
python project/scripts/cross_neo_md/39_build_physics_launch_sheet.py
python project/scripts/cross_neo_md/40_build_endpoint_energy_gate.py
python project/scripts/cross_neo_md/45_build_flagship_execution_packet.py
python project/scripts/cross_neo_md/46_build_flagship_command_center.py
python project/scripts/cross_neo_md/29_build_high_impact_decision_web.py
```

## Claim boundary

This tower coordinates escalation and stop rules. It does not prove immunogenicity, activation, killing, or clinical utility.

