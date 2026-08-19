# CROSS-Neo Physics Gate Package

## Executive verdict

This is the last expensive computational gate before wetlab specificity testing. Use cheap interface energies first, PMF/umbrella only when TCR-plausibility matters, and FEP only for the two flagship cases after WT mapping is clean.

## Physics modules

| module                                       | physics_question                                                                                    | when_to_use                                                                                   | output                                                     | claim_boundary                                                   | cost_class   |
|:---------------------------------------------|:----------------------------------------------------------------------------------------------------|:----------------------------------------------------------------------------------------------|:-----------------------------------------------------------|:-----------------------------------------------------------------|:-------------|
| endpoint_mmgbsa_or_openmm_interaction_energy | Does the mutant preserve a more favorable interface than WT/decoy over an ensemble?                 | After short MD exists; cheap ranking before PMF/FEP.                                          | trajectory endpoint interaction energies with bootstrap CI | ranking diagnostic only; no activation or cytokine claim         | low          |
| umbrella_or_steered_md_tcr_unbinding_pmf     | Does the mutant/TCR interface remain harder to unbind than WT/decoy?                                | Only if TCR-pMHC interface is already stable and multimer/TCR-binding decision depends on it. | PMF / off-rate proxy / rupture path risk                   | mechanistic proxy only; not cytokine or killing proof            | high         |
| alchemical_fep_mutant_to_wt_delta_delta_g    | Is mutant more favorable than WT in a formally estimated free-energy sense?                         | Only for top flagship cases after WT mapping and stable 10ns controls.                        | DeltaDeltaG with uncertainty                               | does not prove immunogenicity; only supports specificity ranking | very_high    |
| structure_confidence_and_clash_crosscheck    | Is the starting pHLA/TCR-pMHC geometry physically plausible enough to justify expensive simulation? | Before any expensive physics module is launched.                                              | clash / anchor / geometry sanity score                     | presentation plausibility only                                   | very_low     |

## Candidate plan

| peptide    | hla_4digit   |   physics_priority_score | physics_tier        | recommended_physics_module                                                                                                          | recommended_next_step                                                              | cost_class        |
|:-----------|:-------------|-------------------------:|:--------------------|:------------------------------------------------------------------------------------------------------------------------------------|:-----------------------------------------------------------------------------------|:------------------|
| HMTEVVRHC  | HLA-A*02:01  |                 0.839191 | P0_FLAGSHIP_PHYSICS | endpoint_mmgbsa_or_openmm_interaction_energy + umbrella_or_steered_md_tcr_unbinding_pmf + alchemical_fep_mutant_to_wt_delta_delta_g | run endpoint energy now; hold PMF/FEP until WT/decoy control geometry is confirmed | high_to_very_high |
| GADGVGKSAL | HLA-C*08:02  |                 0.656129 | P0_FLAGSHIP_PHYSICS | endpoint_mmgbsa_or_openmm_interaction_energy + umbrella_or_steered_md_tcr_unbinding_pmf                                             | run endpoint energy now; PMF only if paired TCR model remains stable               | high              |
| ILDKVLVHL  | HLA-A*02:01  |                 0.390074 | P3_NO_PHYSICS_NOW   | structure_confidence_and_clash_crosscheck                                                                                           | hold until assay or MD uncertainty justifies more physics budget                   | very_low          |
| SYLDSGIHF  | HLA-A*24:02  |                 0.386389 | P3_NO_PHYSICS_NOW   | structure_confidence_and_clash_crosscheck                                                                                           | hold until assay or MD uncertainty justifies more physics budget                   | very_low          |
| KLILWRGLK  | HLA-A*03:01  |                 0.385522 | P3_NO_PHYSICS_NOW   | structure_confidence_and_clash_crosscheck                                                                                           | hold until assay or MD uncertainty justifies more physics budget                   | very_low          |
| YVDFREYEYY | HLA-A*01:01  |                 0.370838 | P3_NO_PHYSICS_NOW   | structure_confidence_and_clash_crosscheck                                                                                           | hold until assay or MD uncertainty justifies more physics budget                   | very_low          |
| VVGAVGVGK  | HLA-A*11:01  |                 0.364997 | P3_NO_PHYSICS_NOW   | structure_confidence_and_clash_crosscheck                                                                                           | hold until assay or MD uncertainty justifies more physics budget                   | very_low          |
| KLYGLDWAEL | HLA-A*02:01  |                 0.354926 | P3_NO_PHYSICS_NOW   | structure_confidence_and_clash_crosscheck                                                                                           | hold until assay or MD uncertainty justifies more physics budget                   | very_low          |
| ILDTAGREEY | HLA-A*01:01  |                 0.3548   | P3_NO_PHYSICS_NOW   | structure_confidence_and_clash_crosscheck                                                                                           | hold until assay or MD uncertainty justifies more physics budget                   | very_low          |
| ALDPHSGHFV | HLA-A*02:01  |                 0.343416 | P3_NO_PHYSICS_NOW   | structure_confidence_and_clash_crosscheck                                                                                           | hold until assay or MD uncertainty justifies more physics budget                   | very_low          |
| ALDPLLLRI  | HLA-A*02:01  |                 0.332248 | P3_NO_PHYSICS_NOW   | structure_confidence_and_clash_crosscheck                                                                                           | hold until assay or MD uncertainty justifies more physics budget                   | very_low          |
| KMIGNHLWV  | HLA-A*02:01  |                 0.33034  | P3_NO_PHYSICS_NOW   | structure_confidence_and_clash_crosscheck                                                                                           | hold until assay or MD uncertainty justifies more physics budget                   | very_low          |
| GADGVGKSA  | HLA-C*08:02  |                 0.318724 | P3_NO_PHYSICS_NOW   | structure_confidence_and_clash_crosscheck                                                                                           | hold until assay or MD uncertainty justifies more physics budget                   | very_low          |

## Flagship physics candidates

| peptide    | hla_4digit   |   physics_priority_score | physics_tier        | recommended_physics_module                                                                                                          | recommended_next_step                                                              |
|:-----------|:-------------|-------------------------:|:--------------------|:------------------------------------------------------------------------------------------------------------------------------------|:-----------------------------------------------------------------------------------|
| HMTEVVRHC  | HLA-A*02:01  |                 0.839191 | P0_FLAGSHIP_PHYSICS | endpoint_mmgbsa_or_openmm_interaction_energy + umbrella_or_steered_md_tcr_unbinding_pmf + alchemical_fep_mutant_to_wt_delta_delta_g | run endpoint energy now; hold PMF/FEP until WT/decoy control geometry is confirmed |
| GADGVGKSAL | HLA-C*08:02  |                 0.656129 | P0_FLAGSHIP_PHYSICS | endpoint_mmgbsa_or_openmm_interaction_energy + umbrella_or_steered_md_tcr_unbinding_pmf                                             | run endpoint energy now; PMF only if paired TCR model remains stable               |

## Budget

| peptide    | hla_4digit   | module                                       |   estimated_gpu_hours_low |   estimated_gpu_hours_high | decision_use                                        |
|:-----------|:-------------|:---------------------------------------------|--------------------------:|---------------------------:|:----------------------------------------------------|
| HMTEVVRHC  | HLA-A*02:01  | endpoint_mmgbsa_or_openmm_interaction_energy |                         6 |                         20 | cheap ranking and consistency check                 |
| HMTEVVRHC  | HLA-A*02:01  | umbrella_or_steered_md_tcr_unbinding_pmf     |                       180 |                        650 | only if multimer/TCR-binding decision depends on it |
| HMTEVVRHC  | HLA-A*02:01  | alchemical_fep_mutant_to_wt_delta_delta_g    |                       250 |                        900 | WT specificity only, after the cheap gates pass     |
| GADGVGKSAL | HLA-C*08:02  | endpoint_mmgbsa_or_openmm_interaction_energy |                         6 |                         20 | cheap ranking and consistency check                 |
| GADGVGKSAL | HLA-C*08:02  | umbrella_or_steered_md_tcr_unbinding_pmf     |                       180 |                        650 | only if multimer/TCR-binding decision depends on it |
| GADGVGKSAL | HLA-C*08:02  | alchemical_fep_mutant_to_wt_delta_delta_g    |                       250 |                        900 | WT specificity only, after the cheap gates pass     |

## Rules

| observed_pattern                               | physics_gate_response                                        |
|:-----------------------------------------------|:-------------------------------------------------------------|
| endpoint energy supports mutant > WT/decoy     | promote to PMF or FEP shortlist                              |
| PMF supports stronger mutant unbinding barrier | keep candidate alive for wetlab specificity testing          |
| FEP supports mutant > WT specificity           | use as highest-cost physics evidence before assay escalation |
| WT or decoy looks as good as mutant            | hold mutant-specific claim and recheck specificity           |
| physics disagrees with MD and assay readiness  | treat as suspicion flag, not automatic rejection             |
| no paired TCR                                  | do not spend PMF/FEP on recognition claims                   |

## Claim boundary

Physics supports rejection, ranking, and specificity sanity checks. It does not prove activation, killing, clinical utility, or universal immunogenicity.
