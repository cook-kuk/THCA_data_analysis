# Immunogenicity Simulation Escalation Report

## Executive verdict

Add simulations as a staged pre-wetlab culling layer, not as immunogenicity proof. The best immediate value is WT/decoy-aware pMHC and TCR-pMHC stability for the strict Top-13, with expensive FEP/PMF held for the top two only.

## Summary

- planned simulation/analysis modules: 81
- modules with explicit MD length: 38
- requested short-MD trajectory length: 510.0 ns
- ready/manifested/supported modules: 9
- blocked or hold modules: 21

## Highest-impact immediate actions

1. Keep the completed/synced HMTEVVRHC 10 ns evidence locked in the reports and use it as the MD-strong case study.
2. Run/monitor the already prepared diverse pack for GADGVGKSAL mutant TCR-pMHC replicates and same-HLA positive controls.
3. Curate WT peptides and anchor-preserved decoys for Top-4 candidates before spending more GPU on specificity claims.
4. Add endpoint interaction-energy/MMGBSA only after replicate trajectories exist.
5. Reserve FEP or umbrella/steered MD for GADGVGKSAL and HMTEVVRHC only if WT/decoy 10 ns controls are clean.

## Assay bridge

| wetlab_assay                                       | simulation_support                                                                                | decision_use                                                                | not_supported_claim                              |
|:---------------------------------------------------|:--------------------------------------------------------------------------------------------------|:----------------------------------------------------------------------------|:-------------------------------------------------|
| HLA binding / stability assay                      | pMHC structure crosscheck; mutant pMHC 3x10ns; WT/decoy pMHC controls                             | Cull peptides that do not remain in the groove before T-cell assays.        | Does not prove TCR recognition.                  |
| pMHC multimer / tetramer with known TCR or T cells | TCR-pMHC 3x10ns; CDR3-peptide contacts; TCR unbinding PMF only for top cases                      | Prioritize candidates whose TCR interface is persistent and peptide-facing. | Does not prove cytokine secretion or killing.    |
| IFN-gamma ELISpot / ICS                            | TCR contact persistence, mutation-site contacts, WT/decoy specificity, kinetic proofreading proxy | Interpret why a candidate might pass/fail activation assays.                | Cannot simulate complete immune-cell activation. |
| Killing assay                                      | Only upstream presentation/recognition plausibility plus WT/decoy risk                            | Run only after binding and activation screen are positive.                  | No simulation here proves tumor killing.         |

## Go/no-go rules

| stage                | go_rule                                                                                                         | no_go_or_hold                                                                                              |
|:---------------------|:----------------------------------------------------------------------------------------------------------------|:-----------------------------------------------------------------------------------------------------------|
| T0 structure         | Peptide modeled in groove; anchor contacts plausible; no severe clashes; confidence acceptable.                 | Unsupported HLA/length, failed groove placement, heavy clash, or missing WT for specificity-critical case. |
| T1 10ns pMHC MD      | Peptide RMSD/drift plateaus; native pMHC contacts and anchors persist in most replicates.                       | Peptide exits groove, anchor contacts collapse, or replicate behavior is inconsistent.                     |
| T1/T2 TCR-pMHC MD    | Persistent TCR-peptide/CDR3 contacts; mutation-site contact if biologically relevant; TCR not only MHC-focused. | Interface falls apart, contacts are MHC-only, or mutation residue is buried away from TCR in all models.   |
| WT/decoy specificity | Mutant shows stronger pMHC/TCR interface evidence than WT and anchor-preserved decoy.                           | WT/decoy looks equally good or better; self-reactivity risk requires wetlab caution.                       |
| T3 expensive FEP/PMF | Only launch if 10ns mutant/control results are stable and wetlab decision depends on specificity ranking.       | Do not spend FEP/PMF budget on unstable, uncurated, or non-TCR-available cases.                            |

## Candidate budget

| peptide    | hla_4digit   |   planned_jobs |   requested_ns |   gpu_hours_low |   gpu_hours_high |   blocked_jobs |
|:-----------|:-------------|---------------:|---------------:|----------------:|-----------------:|---------------:|
| GADGVGKSAL | HLA-C*08:02  |              4 |            120 |              72 |              144 |              0 |
| HMTEVVRHC  | HLA-A*02:01  |              4 |            120 |              72 |              144 |              2 |
| ILDKVLVHL  | HLA-A*02:01  |              2 |             60 |              36 |               72 |              1 |
| KLILWRGLK  | HLA-A*03:01  |              2 |             60 |              36 |               72 |              1 |
| ALDPHSGHFV | HLA-A*02:01  |              4 |             25 |              15 |               30 |              2 |
| ALDPLLLRI  | HLA-A*02:01  |              4 |             25 |              15 |               30 |              2 |
| KLYGLDWAEL | HLA-A*02:01  |              4 |             25 |              15 |               30 |              2 |
| VVGAVGVGK  | HLA-A*11:01  |              4 |             25 |              15 |               30 |              2 |
| GADGVGKSA  | HLA-C*08:02  |              2 |             10 |               6 |               12 |              1 |
| ILDTAGREEY | HLA-A*01:01  |              2 |             10 |               6 |               12 |              1 |
| KMIGNHLWV  | HLA-A*02:01  |              2 |             10 |               6 |               12 |              1 |
| SYLDSGIHF  | HLA-A*24:02  |              2 |             10 |               6 |               12 |              1 |
| YVDFREYEYY | HLA-A*01:01  |              2 |             10 |               6 |               12 |              1 |

## Claim boundary

Simulation can support pMHC stability, TCR-interface plausibility, and mutant-over-WT/decoy specificity hypotheses. It cannot prove IFN-gamma release, T-cell activation, killing, clinical response, or universal immunogenicity.

