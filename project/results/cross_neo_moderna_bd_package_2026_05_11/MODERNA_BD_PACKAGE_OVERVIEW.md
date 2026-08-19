# CROSS-Neo Moderna-Ready BD Package

## Executive Decision

The package is now framed as a selection-layer collaboration, not a claim that CROSS-Neo replaces Moderna's mRNA platform.

## Core Evidence

| algorithm        | algorithm_family          |   all_AUPRC |   all_AUROC |   frozen_validation_AUPRC |   low_medium_leakage_AUPRC |   v6_smoke_AUPRC | primary_disposition                                                                                    |
|:-----------------|:--------------------------|------------:|------------:|--------------------------:|---------------------------:|-----------------:|:-------------------------------------------------------------------------------------------------------|
| KG_GA_evolved    | new_ga_controller         |    0.933411 |    0.943428 |                  0.97775  |                   0.757078 |         0.811595 | CHAMPION for retrospective + validation-like + low-leakage benchmark; freeze before prospective wetlab |
| CROSS_integrated | fixed_integrated_crossneo |    0.878422 |    0.895868 |                  0.713841 |                   0.57277  |         0.790492 | Fixed integrated comparator; strong but less novel than KG-GA                                          |
| CROSS_stress     | stress_guarded_internal   |    0.871905 |    0.873389 |                  0.650789 |                   0.600589 |         0.776958 | Comparator or component                                                                                |
| CROSS_claimsafe  | claim_safe_crossneo       |    0.860895 |    0.863588 |                  0.84648  |                   0.616845 |         0.788084 | Claim-safe fallback when GA/RL proxy/product features are too risky                                    |
| CROSS_finetuned  | finetuned_internal        |    0.841371 |    0.862878 |                  0.605625 |                   0.668925 |         0.681142 | Comparator or component                                                                                |
| BigMHC_IM        | public_immunogenicity     |    0.671609 |    0.722291 |                  0.333705 |                   0.373655 |         0.672325 | External public comparator/component, not CROSS-Neo claim engine                                       |
| BigMHC_EL        | public_presentation       |    0.647727 |    0.71921  |                  0.30759  |                   0.306236 |         0.865599 | External public comparator/component, not CROSS-Neo claim engine                                       |

## Strict Reliability Context

| algorithm                  | algorithm_family              | claim_track                                                   |   neo_strict_AUPRC |   neo_strict_AUROC |   patients_evaluated | primary_disposition                                                                            |
|:---------------------------|:------------------------------|:--------------------------------------------------------------|-------------------:|-------------------:|---------------------:|:-----------------------------------------------------------------------------------------------|
| BAR-Neo_confidence         | production_reliability_ranker | strict no-overlap production                                  |           0.856903 |           0.97994  |                    1 | Best strict no-overlap production reliability baseline; compare separately from 2,715-board GA |
| BAR-Neo                    | production_reliability_ranker | strict no-overlap production                                  |           0.825318 |           0.978856 |                    1 | Best strict no-overlap production reliability baseline; compare separately from 2,715-board GA |
| RF_biophys                 | neoimmune_stack_baseline      | strict no-overlap production                                  |           0.751907 |           0.685223 |                    1 | Comparator or component                                                                        |
| MHCflurry_2.0_presentation | neoimmune_stack_baseline      | strict no-overlap production                                  |           0.725841 |           0.620698 |                    1 | Comparator or component                                                                        |
| KG_GA_evolved_controller   | new_ga_controller             | production/experiment-priority + strict no-overlap production |           0.624074 |           0.897086 |                   82 | Production GA/RL branch in NeoImmune stack; strong but not clean-track allowed                 |

## Broader Public Competitor Context

| method        |   n_context_rows |   total_scored |   datasets_with_auroc |   mean_context_AUROC |   best_context_AUROC |   min_context_AUROC | claim_use                                                                                    |
|:--------------|-----------------:|---------------:|----------------------:|---------------------:|---------------------:|--------------------:|:---------------------------------------------------------------------------------------------|
| RF_biophys    |                8 |           3816 |                     4 |             0.8914   |               0.9731 |              0.7531 | context-only, AUROC-only in this repo; do not use as locked same-board AUPRC unless rescored |
| PRIME         |                8 |           3839 |                     4 |             0.6691   |               0.8054 |              0.5673 | context-only, AUROC-only in this repo; do not use as locked same-board AUPRC unless rescored |
| ESM2_Bayesian |                1 |            311 |                     1 |             0.7841   |               0.7841 |              0.7841 | context-only, AUROC-only in this repo; do not use as locked same-board AUPRC unless rescored |
| MHCflurry     |                8 |           3578 |                     4 |             0.696225 |               0.7564 |              0.6367 | context-only, AUROC-only in this repo; do not use as locked same-board AUPRC unless rescored |
| BigMHC_IM     |                5 |           2318 |                     2 |             0.645    |               0.7477 |              0.5423 | context-only, AUROC-only in this repo; do not use as locked same-board AUPRC unless rescored |
| TransPHLA     |                8 |           3587 |                     4 |             0.681    |               0.7281 |              0.5797 | context-only, AUROC-only in this repo; do not use as locked same-board AUPRC unless rescored |
| DeepImmuno    |                1 |            319 |                     1 |             0.7019   |               0.7019 |              0.7019 | context-only, AUROC-only in this repo; do not use as locked same-board AUPRC unless rescored |
| Structure_LR  |                1 |            311 |                     1 |             0.6894   |               0.6894 |              0.6894 | context-only, AUROC-only in this repo; do not use as locked same-board AUPRC unless rescored |
| VQC           |                1 |            319 |                     1 |             0.5985   |               0.5985 |              0.5985 | context-only, AUROC-only in this repo; do not use as locked same-board AUPRC unless rescored |
| GP_quantum    |                1 |            319 |                     1 |             0.587    |               0.587  |              0.587  | context-only, AUROC-only in this repo; do not use as locked same-board AUPRC unless rescored |
| NetMHCpan_4.1 |                1 |            319 |                     1 |             0.5669   |               0.5669 |              0.5669 | context-only, AUROC-only in this repo; do not use as locked same-board AUPRC unless rescored |

## Claim Boundary

| claim                       | wording                                                                                                                    | why                                                               | do_not_say                                                 |
|:----------------------------|:---------------------------------------------------------------------------------------------------------------------------|:------------------------------------------------------------------|:-----------------------------------------------------------|
| Allowed pre-NDA             | Retrospective CROSS-Neo prioritization engine with AUPRC 0.933 on a 2,715-candidate local board                            | Supported by same-board metrics                                   | Clinically validated or Moderna-ready efficacy proof       |
| Allowed under NDA           | Frozen algorithm can be evaluated on blinded sponsor candidate pools for top-34 enrichment                                 | Validation protocol is pre-specified                              | Will outperform sponsor stack before seeing blinded result |
| Comparator boundary         | Same-board comparison beats BigMHC_IM/EL and fixed CROSS-Neo; public context includes RF_biophys/MHCflurry/PRIME/TransPHLA | Different benchmark tables have different leakage/overlap caveats | Universal SOTA over all public algorithms                  |
| Strict reliability boundary | BAR-Neo_confidence is still the strict no-overlap reliability winner in current local stack                                | Strict AUPRC 0.857 vs KG_GA_controller 0.624                      | KG-GA wins every split                                     |
| Wetlab boundary             | Assay-ready prioritization and control design                                                                              | No prospective mutant-vs-WT/decoy readout yet                     | Immunogenicity proof                                       |

## Data Room Index

| folder                 | item                                          | share_timing                             | contains                                                                   |
|:-----------------------|:----------------------------------------------|:-----------------------------------------|:---------------------------------------------------------------------------|
| 00_non_confidential    | MODERNA_NONCONFIDENTIAL_TEASER.md             | pre-NDA                                  | Value proposition, masked metrics, validation ask, no algorithm internals  |
| 01_validation_protocol | MODERNA_BLINDED_VALIDATION_PROTOCOL.md        | NDA                                      | Input schema, frozen scoring, endpoints, success criteria                  |
| 02_metrics             | moderna_same_board_comparator_metrics.tsv     | NDA                                      | KG-GA vs CROSS/BigMHC same-board metrics                                   |
| 02_metrics             | moderna_strict_reliability_metrics.tsv        | NDA                                      | BAR-Neo, RF_biophys, MHCflurry, KG-GA-controller strict no-overlap context |
| 03_claim_boundary      | moderna_claim_boundary.tsv                    | pre-NDA summary / NDA detail             | Allowed claims and blocked claims                                          |
| 04_deal                | MODERNA_OPTION_LICENSE_TERM_SHEET_SCAFFOLD.md | after technical interest                 | Option/license terms scaffold; not legal advice                            |
| 05_ip                  | INVENTION_DISCLOSURE_SCAFFOLD.md              | attorney / tech transfer before outreach | Provisional-patent disclosure skeleton and claim areas                     |
| 06_outreach            | MODERNA_OUTREACH_EMAIL_AFTER_IP.md            | after IP hygiene and contact selection   | Short pre-NDA email draft that does not disclose internals                 |

## Contact Routing

| company       | route                                     | url                                                                | why                                                                                                              | pre_nda_action                                                                                       | confidentiality_note                                                                                        |
|:--------------|:------------------------------------------|:-------------------------------------------------------------------|:-----------------------------------------------------------------------------------------------------------------|:-----------------------------------------------------------------------------------------------------|:------------------------------------------------------------------------------------------------------------|
| Merck         | Business development and licensing (BD&L) | https://www.merck.com/research/business-development-and-licensing/ | Official Merck collaboration route; includes oncology, vaccines, and enabling/platform technologies.             | Use non-confidential teaser only; request CDA/NDA and blinded validation discussion.                 | Do not send source code, feature list, candidate tables, or scoring weights before CDA/NDA.                 |
| Merck         | Corporate contact fallback                | https://www.merck.com/contact-us/                                  | Fallback for routing to BD&L or oncology external innovation if web BD intake is not available.                  | Ask for the correct BD&L contact for individualized neoantigen vaccine selection technology.         | General contact route; keep all technical details non-confidential.                                         |
| Moderna       | Official contact page                     | https://www.modernatx.com/en-US/contact-moderna                    | Official Moderna page lists corporate contact routes, HQ, media, investor, and general contact channels.         | Use only to request routing or a warm BD contact; do not disclose internals through general inboxes. | No public BD-specific email was verified here; treat as routing only, not a confidential submission portal. |
| Moderna/Merck | Warm intro / conference BD route          | N/A                                                                | Best practical route for a sensitive platform/selection-layer pitch when no public Moderna BD inbox is verified. | Send one-paragraph teaser and ask for CDA/NDA-enabled technical diligence.                           | Do not attach NDA-stage validation protocol unless CDA/NDA is in place.                                     |

