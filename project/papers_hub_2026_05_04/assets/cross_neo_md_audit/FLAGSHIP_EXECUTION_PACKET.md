# CROSS-Neo flagship execution packet

## Executive verdict

This packet converts the two flagship candidates into a concrete 72-hour action board. It is designed for coordination, not for claim inflation.

## Anchor numbers

- candidates: 13
- physics flagships: 2
- preclinical-ready: 2
- real assay summaries: 0

## 72h execution board

| lead                 | peptide    | hla_4digit   |   top_md_score |   preclinical_readiness |   physics_priority_score | preclinical_next_action                                               | physics_tier        | launch_sequence                                                                                          | endpoint_step                                                        | command_priority                                                                    | next_72h_goal                                                 | claim_boundary                                  |
|:---------------------|:-----------|:-------------|---------------:|------------------------:|-------------------------:|:----------------------------------------------------------------------|:--------------------|:---------------------------------------------------------------------------------------------------------|:---------------------------------------------------------------------|:------------------------------------------------------------------------------------|:--------------------------------------------------------------|:------------------------------------------------|
| TP53_R175H_HLA_A0201 | HMTEVVRHC  | HLA-A*02:01  |       0.821789 |                0.935447 |                 0.839191 | launch WT/decoy-controlled HLA stability + multimer/activation assays | P0_FLAGSHIP_PHYSICS | endpoint energy now -> PMF if endpoint supports mutant specificity -> FEP only after WT mapping is clean | compare mutant vs WT vs decoy endpoint energies over ensemble frames | run endpoint energy now; keep PMF/FEP held unless the endpoint gate stays favorable | WT/decoy-controlled HLA stability + multimer/activation setup | prioritized execution, not immunogenicity proof |
| KRAS_G12D_HLA_C0802  | GADGVGKSAL | HLA-C*08:02  |       0.577335 |                0.796652 |                 0.656129 | launch WT/decoy-controlled HLA stability + multimer/activation assays | P0_FLAGSHIP_PHYSICS | endpoint energy now -> PMF only if TCR-pMHC model stays stable                                           | compare mutant vs WT vs decoy endpoint energies over ensemble frames | run endpoint energy now; keep PMF/FEP held unless the endpoint gate stays favorable | WT/decoy-controlled HLA stability + multimer/activation setup | prioritized execution, not immunogenicity proof |

## Required commands

```bash
python project/scripts/cross_neo_md/44_build_executive_impact_console.py
python project/scripts/cross_neo_md/43_build_decision_atlas.py
python project/scripts/cross_neo_md/42_build_integrated_decision_matrix.py
python project/scripts/cross_neo_md/41_build_physics_case_study_board.py
python project/scripts/cross_neo_md/38_build_physics_gate_package.py
python project/scripts/cross_neo_md/39_build_physics_launch_sheet.py
python project/scripts/cross_neo_md/40_build_endpoint_energy_gate.py
python project/scripts/cross_neo_md/29_build_high_impact_decision_web.py
```

## Claim boundary

The packet supports execution and prioritization. It does not prove immunogenicity, activation, killing, or clinical utility.

