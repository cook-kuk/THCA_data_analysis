# CROSS-Neo executive impact console

## Verdict

The current strongest story is a two-flagship, physics-anchored, WT/decoy-controlled validation plan. It is still a prioritization and escalation system, not immunogenicity proof.

## Key numbers

- candidates: 13
- paired TCR evidence rows: 6
- physics flagships: 2
- preclinical-ready: 2
- real assay summaries: 0

## Flagship verdict

- top MD score: 0.822
- top preclinical readiness: 0.935
- top physics priority: 0.839

## Immediate actions

1. Run WT/decoy-controlled HLA stability and pHLA binding for both flagship candidates.
2. Advance multimer/activation only if the presentation gate stays favorable.
3. Use endpoint-energy comparisons before spending PMF/FEP budget.
4. Treat any real assay row as posterior update input; do not use simulated rows as labels.

## Refresh commands

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

## Forbidden claims

- immunogenicity proof
- clinical utility
- universal TCR-aware prediction
- MD as final validation

| row_id     | peptide    | hla_4digit   | source_dataset   | prediction_outcome   |   main_dl_score |   tcr_augmented_score_mean |   paired_tcr_evidence_count |   bayes_mean | md_label       |   md_score |   preclinical_readiness |   cross_neo_i_readiness | claim_state               | next_experiment                                   | physics_tier        |   physics_priority_score | launch_sequence                                                                                          | endpoint_step                                                        | case_boundary                                                    |   scenario_pressure | integrated_next_action                                                | integrated_stage           |
|:-----------|:-----------|:-------------|:-----------------|:---------------------|----------------:|---------------------------:|----------------------------:|-------------:|:---------------|-----------:|------------------------:|------------------------:|:--------------------------|:--------------------------------------------------|:--------------------|-------------------------:|:---------------------------------------------------------------------------------------------------------|:---------------------------------------------------------------------|:-----------------------------------------------------------------|--------------------:|:----------------------------------------------------------------------|:---------------------------|
| CNV0_01132 | GADGVGKSAL | HLA-C*08:02  | NEPdb            | TP                   |        0.770461 |                   0.959699 |                          13 |     0.558867 | MD_MODERATE    |   0.616608 |                0.796652 |                0.61065  | PRIOR_ONLY_NO_ASSAY_LABEL | WT/decoy-controlled HLA stability or pHLA binding | P0_FLAGSHIP_PHYSICS |                 0.656129 | endpoint energy now -> PMF only if TCR-pMHC model stays stable                                           | compare mutant vs WT vs decoy endpoint energies over ensemble frames | physics and assay support ranking only, not immunogenicity proof |            0.100662 | launch WT/decoy-controlled HLA stability + multimer/activation assays | ASSAY_AND_PHYSICS_PRIORITY |
| CNV0_01151 | HMTEVVRHC  | HLA-A*02:01  | NEPdb            | TP                   |        0.467097 |                   0.983185 |                          18 |     0.409026 | MD_VERY_STRONG |   0.821789 |                0.935447 |                0.605838 | PRIOR_ONLY_NO_ASSAY_LABEL | WT/decoy-controlled HLA stability or pHLA binding | P0_FLAGSHIP_PHYSICS |                 0.839191 | endpoint energy now -> PMF if endpoint supports mutant specificity -> FEP only after WT mapping is clean | compare mutant vs WT vs decoy endpoint energies over ensemble frames | physics and assay support ranking only, not immunogenicity proof |            0.100584 | launch WT/decoy-controlled HLA stability + multimer/activation assays | ASSAY_AND_PHYSICS_PRIORITY |
| CNV0_02398 | KLILWRGLK  | HLA-A*03:01  | ITSNdb_main      | TP                   |        0.803188 |                   0.985393 |                           0 |     0.648496 | nan            | nan        |              nan        |                0.269031 | PRIOR_ONLY_NO_ASSAY_LABEL | WT/decoy-controlled HLA stability or pHLA binding | P3_NO_PHYSICS_NOW   |                 0.385522 | cheap geometry/endpoint check only                                                                       | geometry sanity only                                                 | nan                                                              |            0.231201 | DL_PRIMARY_SURVIVOR_TCR_CONTEXT_ONLY                                  | ASSAY_AND_PHYSICS_PRIORITY |
| CNV0_02448 | ILDKVLVHL  | HLA-A*02:01  | ITSNdb_main      | TP                   |        0.773359 |                   0.96225  |                           0 |     0.63716  | nan            | nan        |              nan        |                0.260519 | PRIOR_ONLY_NO_ASSAY_LABEL | WT/decoy-controlled HLA stability or pHLA binding | P3_NO_PHYSICS_NOW   |                 0.390074 | cheap geometry/endpoint check only                                                                       | geometry sanity only                                                 | nan                                                              |            0.236014 | DL_PRIMARY_SURVIVOR_TCR_CONTEXT_ONLY                                  | ASSAY_AND_PHYSICS_PRIORITY |
| CNV0_02478 | YVDFREYEYY | HLA-A*01:01  | ITSNdb_main      | TP                   |        0.661379 |                   0.839473 |                           0 |     0.551744 | nan            | nan        |              nan        |                0.224618 | PRIOR_ONLY_NO_ASSAY_LABEL | WT/decoy-controlled HLA stability or pHLA binding | P3_NO_PHYSICS_NOW   |                 0.370838 | cheap geometry/endpoint check only                                                                       | geometry sanity only                                                 | nan                                                              |            0.256314 | DL_PRIMARY_SURVIVOR_TCR_CONTEXT_ONLY                                  | ASSAY_AND_PHYSICS_PRIORITY |
| CNV0_01118 | GADGVGKSA  | HLA-C*08:02  | NEPdb            | TP                   |        0.780925 |                   0.566884 |                           0 |     0.560605 | nan            | nan        |              nan        |                0.218542 | PRIOR_ONLY_NO_ASSAY_LABEL | WT/decoy-controlled HLA stability or pHLA binding | P3_NO_PHYSICS_NOW   |                 0.318724 | nan                                                                                                      | nan                                                                  | nan                                                              |            0.259749 | DL_TCR_DISCORDANT_AUDIT_NOT_POSITIVE                                  | ASSAY_AND_PHYSICS_PRIORITY |

