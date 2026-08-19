# CROSS-Neo Assay Feedback Scenario Simulator

## Purpose

This package shows how the closed-loop learner will behave once real WT/decoy-controlled assay data arrive. These rows are explicitly simulated and must not be treated as observed labels.

## Lead Scenario Outcomes

| peptide    | hla_4digit   | scenario                                 |   prior_cross_neo_i_readiness |   scenario_posterior |   posterior_delta | claim_state                             | allowed_claim_if_real                                              | next_step_if_real                                  |
|:-----------|:-------------|:-----------------------------------------|------------------------------:|---------------------:|------------------:|:----------------------------------------|:-------------------------------------------------------------------|:---------------------------------------------------|
| GADGVGKSAL | HLA-C*08:02  | presentation_fail                        |                      0.61065  |             0.518211 |       -0.092439   | PRESENTATION_FAILED                     | no biological claim; hold or repeat presentation                   | reagent QC or deprioritize                         |
| GADGVGKSAL | HLA-C*08:02  | presentation_only                        |                      0.61065  |             0.629322 |        0.0186721  | PRESENTATION_SUPPORTED                  | presentation plausibility only                                     | mutant/WT/decoy multimer or TCR reporter           |
| GADGVGKSAL | HLA-C*08:02  | recognition_positive_activation_missing  |                      0.61065  |             0.743377 |        0.132727   | RECOGNITION_PLAUSIBLE                   | mutant-specific recognition plausibility in tested context         | ELISpot/ICS/CD137 activation assay                 |
| GADGVGKSAL | HLA-C*08:02  | activation_specific_positive             |                      0.61065  |             0.791494 |        0.180844   | ASSAY_SPECIFIC_IMMUNOGENICITY_SUPPORTED | assay-specific immunogenicity in tested context                    | matched HLA/mutation killing assay                 |
| GADGVGKSAL | HLA-C*08:02  | functional_best_case                     |                      0.61065  |             0.833195 |        0.222545   | FUNCTIONAL_SPECIFICITY_SUPPORTED        | functional specificity in tested HLA/mutation system               | independent replicate or external validation       |
| GADGVGKSAL | HLA-C*08:02  | recognition_positive_activation_negative |                      0.61065  |             0.603994 |       -0.00665641 | RECOGNITION_WITHOUT_ACTIVATION          | recognition without productive activation; no immunogenicity claim | dose response, alternate TCR context, or downgrade |
| GADGVGKSAL | HLA-C*08:02  | wt_cross_reactive                        |                      0.61065  |             0.534942 |       -0.0757081  | SPECIFICITY_RISK_HOLD_WT_POSITIVE       | no mutant-specific claim; WT cross-reactivity risk                 | hold candidate; investigate WT cross-reactivity    |
| GADGVGKSAL | HLA-C*08:02  | decoy_positive                           |                      0.61065  |             0.534942 |       -0.0757081  | SPECIFICITY_RISK_HOLD_DECOY_POSITIVE    | no specificity claim; decoy/nonspecific signal risk                | repeat with orthogonal decoy and artifact controls |
| HMTEVVRHC  | HLA-A*02:01  | presentation_fail                        |                      0.605838 |             0.515003 |       -0.090835   | PRESENTATION_FAILED                     | no biological claim; hold or repeat presentation                   | reagent QC or deprioritize                         |
| HMTEVVRHC  | HLA-A*02:01  | presentation_only                        |                      0.605838 |             0.626114 |        0.0202761  | PRESENTATION_SUPPORTED                  | presentation plausibility only                                     | mutant/WT/decoy multimer or TCR reporter           |
| HMTEVVRHC  | HLA-A*02:01  | recognition_positive_activation_missing  |                      0.605838 |             0.741156 |        0.135318   | RECOGNITION_PLAUSIBLE                   | mutant-specific recognition plausibility in tested context         | ELISpot/ICS/CD137 activation assay                 |
| HMTEVVRHC  | HLA-A*02:01  | activation_specific_positive             |                      0.605838 |             0.789689 |        0.183851   | ASSAY_SPECIFIC_IMMUNOGENICITY_SUPPORTED | assay-specific immunogenicity in tested context                    | matched HLA/mutation killing assay                 |
| HMTEVVRHC  | HLA-A*02:01  | functional_best_case                     |                      0.605838 |             0.831751 |        0.225913   | FUNCTIONAL_SPECIFICITY_SUPPORTED        | functional specificity in tested HLA/mutation system               | independent replicate or external validation       |
| HMTEVVRHC  | HLA-A*02:01  | recognition_positive_activation_negative |                      0.605838 |             0.602189 |       -0.00364889 | RECOGNITION_WITHOUT_ACTIVATION          | recognition without productive activation; no immunogenicity claim | dose response, alternate TCR context, or downgrade |
| HMTEVVRHC  | HLA-A*02:01  | wt_cross_reactive                        |                      0.605838 |             0.533423 |       -0.0724156  | SPECIFICITY_RISK_HOLD_WT_POSITIVE       | no mutant-specific claim; WT cross-reactivity risk                 | hold candidate; investigate WT cross-reactivity    |
| HMTEVVRHC  | HLA-A*02:01  | decoy_positive                           |                      0.605838 |             0.533423 |       -0.0724156  | SPECIFICITY_RISK_HOLD_DECOY_POSITIVE    | no specificity claim; decoy/nonspecific signal risk                | repeat with orthogonal decoy and artifact controls |

## Claim Transitions

| scenario                                 | claim_state                             | allowed_claim_if_real                                              |   n_candidates |   mean_posterior |   mean_delta |
|:-----------------------------------------|:----------------------------------------|:-------------------------------------------------------------------|---------------:|-----------------:|-------------:|
| activation_specific_positive             | ASSAY_SPECIFIC_IMMUNOGENICITY_SUPPORTED | assay-specific immunogenicity in tested context                    |             13 |         0.659787 |    0.400355  |
| decoy_positive                           | SPECIFICITY_RISK_HOLD_DECOY_POSITIVE    | no specificity claim; decoy/nonspecific signal risk                |             13 |         0.424031 |    0.164599  |
| functional_best_case                     | FUNCTIONAL_SPECIFICITY_SUPPORTED        | functional specificity in tested HLA/mutation system               |             13 |         0.72783  |    0.468398  |
| presentation_fail                        | PRESENTATION_FAILED                     | no biological claim; hold or repeat presentation                   |             13 |         0.284066 |    0.0246339 |
| presentation_only                        | PRESENTATION_SUPPORTED                  | presentation plausibility only                                     |             13 |         0.395177 |    0.135745  |
| recognition_positive_activation_missing  | RECOGNITION_PLAUSIBLE                   | mutant-specific recognition plausibility in tested context         |             13 |         0.581276 |    0.321844  |
| recognition_positive_activation_negative | RECOGNITION_WITHOUT_ACTIVATION          | recognition without productive activation; no immunogenicity claim |             13 |         0.472287 |    0.212855  |
| wt_cross_reactive                        | SPECIFICITY_RISK_HOLD_WT_POSITIVE       | no mutant-specific claim; WT cross-reactivity risk                 |             13 |         0.424031 |    0.164599  |

## Boundary

Use this for assay planning and reviewer explanation only. It is not a result and not a substitute for real activation/killing assays.
