# DL-First Candidate Funnel

## Intent

This funnel runs cheap CROSS-Neo model evidence before expensive structure or MD. MD is a late audit/escalation layer, not an early rescue signal.

## Stage Counts

| stage                              |   n |
|:-----------------------------------|----:|
| 00_all_candidates                  | 649 |
| 01_supported_peptide_length        | 491 |
| 02_main_dl_broad_pass              | 280 |
| 03_bayesian_dropout_pass           |  97 |
| 04_robust_main_dl_pass             |  97 |
| 05_source_compatible_tcr_after_dl  |  90 |
| 06_paired_tcr_after_dl             |   9 |
| 07_tcr_branch_supported_after_dl   |  18 |
| 08_tcr_rescue_reviews              |   1 |
| 09_md_escalation_after_cheap_gates |   2 |
| 10_wetlab_shortlist_after_dl       |   1 |

## Decision Counts

| decision                             |   n |
|:-------------------------------------|----:|
| DL_CULL_LOW_MAIN_MODEL_SCORE         | 211 |
| DL_CULL_UNCERTAIN_OR_LOW_POSTERIOR   | 182 |
| DL_CULL_UNSUPPORTED_PEPTIDE_LENGTH   | 158 |
| DL_TCR_DISCORDANT_AUDIT_NOT_POSITIVE |  72 |
| DL_PRIMARY_SURVIVOR_TCR_CONTEXT_ONLY |  17 |
| DL_PRIMARY_SURVIVOR_MODEL_ONLY       |   7 |
| DL_PRIMARY_SURVIVOR_WITH_PAIRED_TCR  |   1 |
| TCR_RESCUE_REVIEW_AFTER_WEAK_MAIN_DL |   1 |

## Highest Priority After Cheap Gates

| row_id     | peptide    | hla_4digit   | dl_first_decision                    | next_action                                          |   dl_first_priority_score |   main_dl_score |   bayes_mean |   bayes_q05 |   bayes_q95 |   perturb_prob_gt_050 |   tcr_augmented_score_mean |   paired_tcr_evidence_count | md_label       |
|:-----------|:-----------|:-------------|:-------------------------------------|:-----------------------------------------------------|--------------------------:|----------------:|-------------:|------------:|------------:|----------------------:|---------------------------:|----------------------------:|:---------------|
| CNV0_01132 | GADGVGKSAL | HLA-C*08:02  | DL_PRIMARY_SURVIVOR_WITH_PAIRED_TCR  | prepare_wt_decoy_structure_then_md_if_controls_ready |                  0.863726 |        0.770461 |     0.558867 |    0.368107 |    0.742052 |                0.9875 |                  0.959699  |                          13 | MD_MODERATE    |
| CNV0_01151 | HMTEVVRHC  | HLA-A*02:01  | TCR_RESCUE_REVIEW_AFTER_WEAK_MAIN_DL | curate_tcr_wt_then_structure_before_md               |                  0.564203 |        0.467097 |     0.409026 |    0.241737 |    0.586356 |                0.042  |                  0.983185  |                          18 | MD_VERY_STRONG |
| CNV0_02452 | LLDGFLATV  | HLA-A*02:01  | DL_PRIMARY_SURVIVOR_TCR_CONTEXT_ONLY | curate_paired_tcr_or_template_before_md              |                  0.790093 |        0.882663 |     0.703096 |    0.590366 |    0.806007 |                1      |                  0.891245  |                           0 | nan            |
| CNV0_02450 | KELEGILLL  | HLA-B*44:03  | DL_PRIMARY_SURVIVOR_TCR_CONTEXT_ONLY | curate_paired_tcr_or_template_before_md              |                  0.789985 |        0.874081 |     0.717429 |    0.603114 |    0.820701 |                1      |                  0.792053  |                           0 | nan            |
| CNV0_02556 | SLFTLHVLGL | HLA-A*02:01  | DL_PRIMARY_SURVIVOR_TCR_CONTEXT_ONLY | curate_paired_tcr_or_template_before_md              |                  0.769857 |        0.840888 |     0.68327  |    0.569012 |    0.788667 |                1      |                  0.998847  |                           0 | nan            |
| CNV0_02398 | KLILWRGLK  | HLA-A*03:01  | DL_PRIMARY_SURVIVOR_TCR_CONTEXT_ONLY | curate_paired_tcr_or_template_before_md              |                  0.747881 |        0.803188 |     0.648496 |    0.532123 |    0.757687 |                1      |                  0.985393  |                           0 | nan            |
| CNV0_02416 | VLLRALPVL  | HLA-A*02:01  | DL_PRIMARY_SURVIVOR_TCR_CONTEXT_ONLY | curate_paired_tcr_or_template_before_md              |                  0.739079 |        0.789494 |     0.632139 |    0.515007 |    0.74288  |                1      |                  0.890772  |                           0 | nan            |
| CNV0_02448 | ILDKVLVHL  | HLA-A*02:01  | DL_PRIMARY_SURVIVOR_TCR_CONTEXT_ONLY | curate_paired_tcr_or_template_before_md              |                  0.734052 |        0.773359 |     0.63716  |    0.520245 |    0.747441 |                1      |                  0.96225   |                           0 | nan            |
| CNV0_02442 | ETVSEQSNV  | HLA-A*68:01  | DL_PRIMARY_SURVIVOR_TCR_CONTEXT_ONLY | curate_paired_tcr_or_template_before_md              |                  0.719002 |        0.742064 |     0.622809 |    0.502366 |    0.737008 |                1      |                  0.889752  |                           0 | nan            |
| CNV0_01284 | ITDVDELPI  | HLA-A*01:01  | DL_PRIMARY_SURVIVOR_TCR_CONTEXT_ONLY | curate_paired_tcr_or_template_before_md              |                  0.716742 |        0.804363 |     0.535108 |    0.343285 |    0.722326 |                0.9585 |                  0.998756  |                           0 | nan            |
| CNV0_02712 | DKESEEEVS  | HLA-C*12:03  | DL_PRIMARY_SURVIVOR_TCR_CONTEXT_ONLY | curate_paired_tcr_or_template_before_md              |                  0.714422 |        0.729675 |     0.623387 |    0.503705 |    0.736873 |                1      |                  0.845793  |                           0 | nan            |
| CNV0_02559 | FLDPDIGGL  | HLA-A*02:01  | DL_PRIMARY_SURVIVOR_TCR_CONTEXT_ONLY | curate_paired_tcr_or_template_before_md              |                  0.714156 |        0.736568 |     0.610272 |    0.49235  |    0.722861 |                1      |                  0.93625   |                           0 | nan            |
| CNV0_02424 | LADEAEVYL  | HLA-A*02:01  | DL_PRIMARY_SURVIVOR_TCR_CONTEXT_ONLY | curate_paired_tcr_or_template_before_md              |                  0.713597 |        0.737968 |     0.605316 |    0.487248 |    0.718288 |                1      |                  0.737138  |                           0 | nan            |
| CNV0_02572 | TTYSPIGEK  | HLA-A*03:01  | DL_PRIMARY_SURVIVOR_TCR_CONTEXT_ONLY | curate_paired_tcr_or_template_before_md              |                  0.694463 |        0.699815 |     0.584244 |    0.465705 |    0.698707 |                1      |                  0.595154  |                           0 | nan            |
| CNV0_01218 | ETWRETGIF  | HLA-A*26:01  | DL_PRIMARY_SURVIVOR_TCR_CONTEXT_ONLY | curate_paired_tcr_or_template_before_md              |                  0.680275 |        0.718898 |     0.517335 |    0.326305 |    0.706091 |                0.958  |                  0.955009  |                           0 | nan            |
| CNV0_02478 | YVDFREYEYY | HLA-A*01:01  | DL_PRIMARY_SURVIVOR_TCR_CONTEXT_ONLY | curate_paired_tcr_or_template_before_md              |                  0.672548 |        0.661379 |     0.551744 |    0.432914 |    0.66807  |                0.999  |                  0.839473  |                           0 | nan            |
| CNV0_01202 | GVYPMPGTQK | HLA-A*03:01  | DL_PRIMARY_SURVIVOR_TCR_CONTEXT_ONLY | curate_paired_tcr_or_template_before_md              |                  0.661081 |        0.651443 |     0.524058 |    0.334389 |    0.71063  |                0.989  |                  0.963049  |                           0 | nan            |
| CNV0_02592 | KRTNVGILK  | HLA-B*27:05  | DL_PRIMARY_SURVIVOR_TCR_CONTEXT_ONLY | curate_paired_tcr_or_template_before_md              |                  0.648717 |        0.62081  |     0.522587 |    0.403938 |    0.640143 |                0.9865 |                  0.516557  |                           0 | nan            |
| CNV0_02479 | ILDTAGREEY | HLA-A*01:01  | DL_PRIMARY_SURVIVOR_TCR_CONTEXT_ONLY | curate_paired_tcr_or_template_before_md              |                  0.642711 |        0.603585 |     0.523586 |    0.404924 |    0.641107 |                0.9885 |                  0.718606  |                           0 | nan            |
| CNV0_01079 | AQQITKTEV  | HLA-B*52:01  | DL_PRIMARY_SURVIVOR_MODEL_ONLY       | curate_tcr_wt_or_use_as_pmhc_only_candidate          |                  0.72239  |        0.90263  |     0.635048 |    0.443624 |    0.808777 |                0.998  |                  0.989682  |                           0 | nan            |
| CNV0_01011 | SEHEGSGPEL | HLA-B*38:01  | DL_PRIMARY_SURVIVOR_MODEL_ONLY       | curate_tcr_wt_or_use_as_pmhc_only_candidate          |                  0.673142 |        0.798815 |     0.589783 |    0.397134 |    0.77066  |                0.999  |                  0.967027  |                           0 | nan            |
| CNV0_02460 | GIVEGLITTV | HLA-A*02:01  | DL_PRIMARY_SURVIVOR_MODEL_ONLY       | curate_tcr_wt_or_use_as_pmhc_only_candidate          |                  0.66387  |        0.75281  |     0.626372 |    0.509007 |    0.737625 |                1      |                  0.72034   |                           0 | nan            |
| CNV0_01137 | EHEGSGPEL  | HLA-B*38:01  | DL_PRIMARY_SURVIVOR_MODEL_ONLY       | curate_tcr_wt_or_use_as_pmhc_only_candidate          |                  0.646857 |        0.909922 |     0.641665 |    0.450572 |    0.814199 |                0.9995 |                  0.0144646 |                           0 | nan            |
| CNV0_02425 | GIVEGLITT  | HLA-A*02:01  | DL_PRIMARY_SURVIVOR_MODEL_ONLY       | curate_tcr_wt_or_use_as_pmhc_only_candidate          |                  0.645259 |        0.71965  |     0.599056 |    0.480825 |    0.712495 |                1      |                  0.762416  |                           0 | nan            |
| CNV0_02565 | LLSTAEAAV  | HLA-A*02:01  | DL_PRIMARY_SURVIVOR_MODEL_ONLY       | curate_tcr_wt_or_use_as_pmhc_only_candidate          |                  0.579234 |        0.74584  |     0.617342 |    0.499647 |    0.729361 |                1      |                  1         |                         532 | nan            |

## Operating Rule

1. Use the main pMHC/CROSS-Neo ensemble to remove low-scoring or unsupported-length peptides.
2. Use Bayesian shrinkage and dropout-like perturbation to remove candidates whose probability changes too much.
3. Use TCR-aware scores only as an optional rescue/review branch for TCR-available rows.
4. Send only DL survivors or TCR-rescue reviews to WT/decoy structure preparation.
5. Run explicit-solvent MD only on the tiny subset with model support, paired/source-compatible TCR evidence, and controls.

## Claim Boundary

These are prioritization scores, not calibrated immunogenicity probabilities. A candidate that survives this funnel is a better wetlab/MD candidate, not a proven positive.
