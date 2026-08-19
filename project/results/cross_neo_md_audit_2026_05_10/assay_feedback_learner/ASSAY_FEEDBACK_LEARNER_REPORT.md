# CROSS-Neo Assay Feedback Learner

## Decision

This package makes the neoantigen workflow data-driven. CROSS-Neo produces priors; WT/decoy-controlled assay rows update posterior readiness, claim state, and the next experiment queue.

## Current Assay Data Status

- Assay feedback file: `/home/seungho/personal/THCA_data_analysis/project/data/cross_neo_assay_feedback/assay_feedback.tsv`
- Real assay summaries detected: 0
- Mode: prior-only template and active-learning plan

## Candidate Posterior Updates

| row_id     | peptide    | hla_4digit   |   prior_cross_neo_i_readiness |   posterior_immunogenicity_readiness |   posterior_uncertainty |   assay_rows_observed | claim_state               | next_experiment                                   |   active_learning_utility |
|:-----------|:-----------|:-------------|------------------------------:|-------------------------------------:|------------------------:|----------------------:|:--------------------------|:--------------------------------------------------|--------------------------:|
| CNV0_01132 | GADGVGKSAL | HLA-C*08:02  |                      0.61065  |                             0.582988 |                0.243113 |                     0 | PRIOR_ONLY_NO_ASSAY_LABEL | WT/decoy-controlled HLA stability or pHLA binding |                  1.10376  |
| CNV0_01151 | HMTEVVRHC  | HLA-A*02:01  |                      0.605838 |                             0.579379 |                0.243699 |                     0 | PRIOR_ONLY_NO_ASSAY_LABEL | WT/decoy-controlled HLA stability or pHLA binding |                  1.09954  |
| CNV0_02398 | KLILWRGLK  | HLA-A*03:01  |                      0.269031 |                             0.326773 |                0.219992 |                     0 | PRIOR_ONLY_NO_ASSAY_LABEL | WT/decoy-controlled HLA stability or pHLA binding |                  0.739023 |
| CNV0_02448 | ILDKVLVHL  | HLA-A*02:01  |                      0.260519 |                             0.32039  |                0.21774  |                     0 | PRIOR_ONLY_NO_ASSAY_LABEL | WT/decoy-controlled HLA stability or pHLA binding |                  0.728259 |
| CNV0_02478 | YVDFREYEYY | HLA-A*01:01  |                      0.224618 |                             0.293463 |                0.207343 |                     0 | PRIOR_ONLY_NO_ASSAY_LABEL | WT/decoy-controlled HLA stability or pHLA binding |                  0.68196  |
| CNV0_01118 | GADGVGKSA  | HLA-C*08:02  |                      0.218542 |                             0.288907 |                0.20544  |                     0 | PRIOR_ONLY_NO_ASSAY_LABEL | WT/decoy-controlled HLA stability or pHLA binding |                  0.673982 |
| CNV0_01439 | SYLDSGIHF  | HLA-A*24:02  |                      0.202107 |                             0.27658  |                0.200084 |                     0 | PRIOR_ONLY_NO_ASSAY_LABEL | WT/decoy-controlled HLA stability or pHLA binding |                  0.65219  |
| CNV0_02479 | ILDTAGREEY | HLA-A*01:01  |                      0.199764 |                             0.274823 |                0.199295 |                     0 | PRIOR_ONLY_NO_ASSAY_LABEL | WT/decoy-controlled HLA stability or pHLA binding |                  0.649059 |
| CNV0_01044 | ALDPHSGHFV | HLA-A*02:01  |                      0.177326 |                             0.257995 |                0.191433 |                     0 | PRIOR_ONLY_NO_ASSAY_LABEL | WT/decoy-controlled HLA stability or pHLA binding |                  0.618759 |
| CNV0_01129 | KLYGLDWAEL | HLA-A*02:01  |                      0.160567 |                             0.245425 |                0.185192 |                     0 | PRIOR_ONLY_NO_ASSAY_LABEL | WT/decoy-controlled HLA stability or pHLA binding |                  0.595759 |
| CNV0_02405 | KMIGNHLWV  | HLA-A*02:01  |                      0.154978 |                             0.241233 |                0.18304  |                     0 | PRIOR_ONLY_NO_ASSAY_LABEL | WT/decoy-controlled HLA stability or pHLA binding |                  0.588018 |
| CNV0_01193 | VVGAVGVGK  | HLA-A*11:01  |                      0.148319 |                             0.236239 |                0.18043  |                     0 | PRIOR_ONLY_NO_ASSAY_LABEL | WT/decoy-controlled HLA stability or pHLA binding |                  0.578749 |
| CNV0_02708 | ALDPLLLRI  | HLA-A*02:01  |                      0.140354 |                             0.230265 |                0.177243 |                     0 | PRIOR_ONLY_NO_ASSAY_LABEL | WT/decoy-controlled HLA stability or pHLA binding |                  0.567597 |

## Claim State Updates

| row_id     | peptide    | hla_4digit   | claim_state               | allowed_claim_now                 | required_next_gate                                |
|:-----------|:-----------|:-------------|:--------------------------|:----------------------------------|:--------------------------------------------------|
| CNV0_01132 | GADGVGKSAL | HLA-C*08:02  | PRIOR_ONLY_NO_ASSAY_LABEL | computational prioritization only | WT/decoy-controlled HLA stability or pHLA binding |
| CNV0_01151 | HMTEVVRHC  | HLA-A*02:01  | PRIOR_ONLY_NO_ASSAY_LABEL | computational prioritization only | WT/decoy-controlled HLA stability or pHLA binding |
| CNV0_02398 | KLILWRGLK  | HLA-A*03:01  | PRIOR_ONLY_NO_ASSAY_LABEL | computational prioritization only | WT/decoy-controlled HLA stability or pHLA binding |
| CNV0_02448 | ILDKVLVHL  | HLA-A*02:01  | PRIOR_ONLY_NO_ASSAY_LABEL | computational prioritization only | WT/decoy-controlled HLA stability or pHLA binding |
| CNV0_02478 | YVDFREYEYY | HLA-A*01:01  | PRIOR_ONLY_NO_ASSAY_LABEL | computational prioritization only | WT/decoy-controlled HLA stability or pHLA binding |
| CNV0_01118 | GADGVGKSA  | HLA-C*08:02  | PRIOR_ONLY_NO_ASSAY_LABEL | computational prioritization only | WT/decoy-controlled HLA stability or pHLA binding |
| CNV0_01439 | SYLDSGIHF  | HLA-A*24:02  | PRIOR_ONLY_NO_ASSAY_LABEL | computational prioritization only | WT/decoy-controlled HLA stability or pHLA binding |
| CNV0_02479 | ILDTAGREEY | HLA-A*01:01  | PRIOR_ONLY_NO_ASSAY_LABEL | computational prioritization only | WT/decoy-controlled HLA stability or pHLA binding |
| CNV0_02405 | KMIGNHLWV  | HLA-A*02:01  | PRIOR_ONLY_NO_ASSAY_LABEL | computational prioritization only | WT/decoy-controlled HLA stability or pHLA binding |
| CNV0_01044 | ALDPHSGHFV | HLA-A*02:01  | PRIOR_ONLY_NO_ASSAY_LABEL | computational prioritization only | WT/decoy-controlled HLA stability or pHLA binding |
| CNV0_01193 | VVGAVGVGK  | HLA-A*11:01  | PRIOR_ONLY_NO_ASSAY_LABEL | computational prioritization only | WT/decoy-controlled HLA stability or pHLA binding |
| CNV0_02708 | ALDPLLLRI  | HLA-A*02:01  | PRIOR_ONLY_NO_ASSAY_LABEL | computational prioritization only | WT/decoy-controlled HLA stability or pHLA binding |
| CNV0_01129 | KLYGLDWAEL | HLA-A*02:01  | PRIOR_ONLY_NO_ASSAY_LABEL | computational prioritization only | WT/decoy-controlled HLA stability or pHLA binding |

## Next Experiment Queue

| peptide    | hla_4digit   | next_experiment                                   |   active_learning_utility | utility_reason                     |
|:-----------|:-------------|:--------------------------------------------------|--------------------------:|:-----------------------------------|
| GADGVGKSAL | HLA-C*08:02  | WT/decoy-controlled HLA stability or pHLA binding |                  1.10376  | cheapest gate before T-cell assays |
| HMTEVVRHC  | HLA-A*02:01  | WT/decoy-controlled HLA stability or pHLA binding |                  1.09954  | cheapest gate before T-cell assays |
| KLILWRGLK  | HLA-A*03:01  | WT/decoy-controlled HLA stability or pHLA binding |                  0.739023 | cheapest gate before T-cell assays |
| ILDKVLVHL  | HLA-A*02:01  | WT/decoy-controlled HLA stability or pHLA binding |                  0.728259 | cheapest gate before T-cell assays |
| YVDFREYEYY | HLA-A*01:01  | WT/decoy-controlled HLA stability or pHLA binding |                  0.68196  | cheapest gate before T-cell assays |
| GADGVGKSA  | HLA-C*08:02  | WT/decoy-controlled HLA stability or pHLA binding |                  0.673982 | cheapest gate before T-cell assays |
| SYLDSGIHF  | HLA-A*24:02  | WT/decoy-controlled HLA stability or pHLA binding |                  0.65219  | cheapest gate before T-cell assays |
| ILDTAGREEY | HLA-A*01:01  | WT/decoy-controlled HLA stability or pHLA binding |                  0.649059 | cheapest gate before T-cell assays |
| ALDPHSGHFV | HLA-A*02:01  | WT/decoy-controlled HLA stability or pHLA binding |                  0.618759 | cheapest gate before T-cell assays |
| KLYGLDWAEL | HLA-A*02:01  | WT/decoy-controlled HLA stability or pHLA binding |                  0.595759 | cheapest gate before T-cell assays |
| KMIGNHLWV  | HLA-A*02:01  | WT/decoy-controlled HLA stability or pHLA binding |                  0.588018 | cheapest gate before T-cell assays |
| VVGAVGVGK  | HLA-A*11:01  | WT/decoy-controlled HLA stability or pHLA binding |                  0.578749 | cheapest gate before T-cell assays |
| ALDPLLLRI  | HLA-A*02:01  | WT/decoy-controlled HLA stability or pHLA binding |                  0.567597 | cheapest gate before T-cell assays |

## Feedback Rules

| observed_pattern                             | model_update                                              | interpretation                             |
|:---------------------------------------------|:----------------------------------------------------------|:-------------------------------------------|
| mutant positive, WT negative, decoy negative | increase posterior; allow next gate                       | mutant-specific signal                     |
| mutant positive, WT positive                 | penalize posterior; hold specificity claim                | WT cross-reactivity risk                   |
| mutant positive, decoy positive              | penalize posterior; repeat with orthogonal decoy          | nonspecific or motif artifact              |
| presentation positive, activation negative   | downgrade to presentation-only                            | presentation is not recognition/activation |
| recognition positive, activation negative    | mark as binding without productive activation             | TCR binding may be nonproductive           |
| activation positive, killing missing         | promote to killing assay                                  | activation is not functional killing       |
| killing positive with clean controls         | promote to functional specificity claim in tested context | strongest preclinical claim                |

## Required Input Schema

| column                  | type    | description                                                                   |
|:------------------------|:--------|:------------------------------------------------------------------------------|
| row_id                  | string  | CROSS-Neo candidate row_id. Required when available.                          |
| peptide                 | string  | Mutant peptide sequence.                                                      |
| hla_4digit              | string  | HLA allele, e.g. HLA-A*02:01.                                                 |
| wildtype_peptide        | string  | Matched WT peptide. Required for specificity claim.                           |
| decoy_peptide           | string  | Anchor-preserved or scrambled decoy peptide.                                  |
| assay_date              | string  | YYYY-MM-DD if available.                                                      |
| batch_id                | string  | Assay batch or plate identifier.                                              |
| assay_type              | enum    | hla_stability|multimer|tcr_reporter|elispot|ics|activation_marker|killing.    |
| antigen_type            | enum    | mutant|wildtype|decoy|irrelevant|positive_control|background.                 |
| mutant_signal           | float   | Mutant condition signal in assay units.                                       |
| wt_signal               | float   | Matched WT signal in same units.                                              |
| decoy_signal            | float   | Decoy signal in same units.                                                   |
| irrelevant_signal       | float   | Irrelevant pHLA/peptide control signal.                                       |
| background_signal       | float   | No peptide, vehicle, unstimulated, or assay background signal.                |
| positive_control_signal | float   | Positive control signal.                                                      |
| replicate_n             | integer | Number of biological or technical replicates summarized.                      |
| signal_units            | string  | MFI, spots/1e5 cells, %positive, %lysis, etc.                                 |
| pass_threshold          | float   | Pre-specified pass threshold in signal units or fold-change.                  |
| mutant_positive         | bool    | Whether mutant condition passes assay-specific threshold.                     |
| wt_positive             | bool    | Whether WT condition is positive. Positive WT blocks mutant-specific claim.   |
| decoy_positive          | bool    | Whether decoy condition is positive. Positive decoy blocks specificity claim. |
| specificity_pass        | bool    | Mutant exceeds WT/decoy/irrelevant under matched conditions.                  |
| presentation_pass       | bool    | HLA binding/stability/MS presentation pass.                                   |
| recognition_pass        | bool    | TCR/multimer/reporter recognition pass.                                       |
| activation_pass         | bool    | ELISpot/ICS/CD137/cytokine activation pass.                                   |
| killing_pass            | bool    | Target-cell killing pass.                                                     |
| notes                   | string  | Free-text caveat. Do not use this as a model label.                           |

## Boundary

The posterior is a decision-support readiness score for choosing the next experiment. It is not a clinical probability and not universal immunogenicity proof.
