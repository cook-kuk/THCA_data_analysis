# CROSS-Neo-I Immunogenicity Prediction Upgrade

## Executive Verdict

The next impact jump is to make CROSS-Neo-I: a calibrated immunogenicity-prediction layer whose endpoint is mutant-specific T-cell activation above WT/decoy controls. The current package should be framed as model-readiness and experimental design, not as a final immunogenicity probability.

## Model Layers

| layer                      | model_or_feature_family                                                            | needed_inputs                                                                          | target_signal                                                                      | claim_boundary                                                     |
|:---------------------------|:-----------------------------------------------------------------------------------|:---------------------------------------------------------------------------------------|:-----------------------------------------------------------------------------------|:-------------------------------------------------------------------|
| presentation               | BigMHC-like MHC-I presentation/immunogenicity transfer                             | mutant peptide, HLA allele, optional expression/processing features                    | peptide-HLA presentation and presentation-conditioned immunogenicity               | presentation is necessary but not sufficient for T-cell activation |
| foreignness_and_self_delta | PRIME/IMPROVE-like immunogenicity and TCR-recognition determinants                 | mutant peptide, WT peptide, anchor positions, hydrophobic/aromatic core features       | mutant-over-self recognition potential                                             | foreignness score does not prove a productive T-cell response      |
| paired_tcr_recognition     | PISTE/EPACT/TCR-HLA-antigen binding expert                                         | paired alpha/beta TCR where available, peptide, HLA                                    | TCR-pMHC recognition plausibility                                                  | TCR expert should abstain when paired TCR evidence is absent       |
| structure_md_specificity   | OpenMM pMHC/TCR-pMHC stability and mutant-vs-WT/decoy delta                        | mutant, WT, decoy pMHC/TCR-pMHC structures and trajectories                            | stable pMHC and TCR-facing interface under controls                                | MD supports structural plausibility, not cytokine release          |
| tumor_context              | expression, clonality, HLA/B2M/TAP, immune-context features                        | RNA expression, mutation clonality/VAF, HLA expression/LOH, antigen-processing context | whether the pHLA target exists in the tumor and can be seen by T cells             | tumor context is required before vaccine-response claims           |
| calibration_and_abstention | stacked calibrated ensemble with source-heldout and TCR-available abstention gates | all layer outputs plus assay labels                                                    | well-calibrated probability of mutant-specific T-cell activation in tested context | calibration must be measured on heldout/external labels            |

## Label Contract

| label_name                  | positive_definition                                                                                       | negative_definition                                                         | not_allowed_as_positive                                                 | used_for                                  |
|:----------------------------|:----------------------------------------------------------------------------------------------------------|:----------------------------------------------------------------------------|:------------------------------------------------------------------------|:------------------------------------------|
| presentation_positive       | mutant peptide is detected by MS elution or passes validated HLA binding/stability assay                  | matched assay shows no presentation/binding above background                | model-predicted binding only                                            | presentation layer                        |
| recognition_positive        | mutant pHLA binds paired TCR or T cells above WT, decoy, and irrelevant controls                          | no mutant-over-control binding under matched conditions                     | peptide-HLA stability without TCR/cell binding                          | TCR recognition layer                     |
| activation_positive         | ELISpot/ICS/cytokine response to mutant exceeds WT, decoy, and irrelevant controls with replicate support | mutant response absent, WT-like, decoy-like, or only at nonphysiologic dose | multimer binding alone                                                  | primary immunogenicity model endpoint     |
| functional_killing_positive | matched HLA/mutation target cells are killed while WT/decoy/mismatched controls are negative              | no specific killing after presentation and activation gates pass            | activation without target-cell assay                                    | functional validation endpoint            |
| unknown_or_weak             | not used as a positive class                                                                              | not used as a negative class unless a matched assay proves failure          | literature mention, model score, public TCR match without assay context | semi-supervised or abstention bucket only |

## Top Candidate Readiness

| peptide    | hla_4digit   |   presentation_proxy |   paired_tcr_proxy |   md_proxy |   wt_decoy_specificity_ready |   cross_neo_i_readiness_score | cross_neo_i_tier       |
|:-----------|:-------------|---------------------:|-------------------:|-----------:|-----------------------------:|------------------------------:|:-----------------------|
| GADGVGKSAL | HLA-C*08:02  |             0.770461 |           0.820335 |   0.616608 |                            1 |                      0.61065  | VALIDATION_P1_READY    |
| HMTEVVRHC  | HLA-A*02:01  |             0.467097 |           0.945752 |   0.821789 |                            1 |                      0.605838 | VALIDATION_P0_FLAGSHIP |
| KLILWRGLK  | HLA-A*03:01  |             0.803188 |           0.541966 |   0        |                            0 |                      0.269031 | CNI_CURATION_FIRST     |
| ILDKVLVHL  | HLA-A*02:01  |             0.773359 |           0.529238 |   0        |                            0 |                      0.260519 | CNI_CURATION_FIRST     |

## Feature Contract

| model_layer                | feature_block                                                                      | must_have_before_training                                                              | primary_label         | allowed_missingness_strategy                                 | leakage_risk                                                  |
|:---------------------------|:-----------------------------------------------------------------------------------|:---------------------------------------------------------------------------------------|:----------------------|:-------------------------------------------------------------|:--------------------------------------------------------------|
| presentation               | BigMHC-like MHC-I presentation/immunogenicity transfer                             | mutant peptide, HLA allele, optional expression/processing features                    | presentation_positive | abstain or separate expert branch; do not impute as negative | public epitope/TCR overlap and source leakage must be audited |
| foreignness_and_self_delta | PRIME/IMPROVE-like immunogenicity and TCR-recognition determinants                 | mutant peptide, WT peptide, anchor positions, hydrophobic/aromatic core features       | presentation_positive | abstain or separate expert branch; do not impute as negative | public epitope/TCR overlap and source leakage must be audited |
| paired_tcr_recognition     | PISTE/EPACT/TCR-HLA-antigen binding expert                                         | paired alpha/beta TCR where available, peptide, HLA                                    | activation_positive   | abstain or separate expert branch; do not impute as negative | public epitope/TCR overlap and source leakage must be audited |
| structure_md_specificity   | OpenMM pMHC/TCR-pMHC stability and mutant-vs-WT/decoy delta                        | mutant, WT, decoy pMHC/TCR-pMHC structures and trajectories                            | activation_positive   | abstain or separate expert branch; do not impute as negative | public epitope/TCR overlap and source leakage must be audited |
| tumor_context              | expression, clonality, HLA/B2M/TAP, immune-context features                        | RNA expression, mutation clonality/VAF, HLA expression/LOH, antigen-processing context | presentation_positive | abstain or separate expert branch; do not impute as negative | public epitope/TCR overlap and source leakage must be audited |
| calibration_and_abstention | stacked calibrated ensemble with source-heldout and TCR-available abstention gates | all layer outputs plus assay labels                                                    | activation_positive   | abstain or separate expert branch; do not impute as negative | public epitope/TCR overlap and source leakage must be audited |
| final_endpoint             | mutant-over-WT/decoy T-cell activation                                             | matched assay labels with negative controls                                            | activation_positive   | unknown bucket; never train unknown as negative              | assay batch/source split must be held out                     |

## Benchmark Contract

| benchmark                  | split                                                           | metric                                                      | success_criterion                                                         | claim_if_passed                                     |
|:---------------------------|:----------------------------------------------------------------|:------------------------------------------------------------|:--------------------------------------------------------------------------|:----------------------------------------------------|
| presentation_heldout       | source-heldout eluted ligand / binding assay                    | AUPRC, precision@k, calibration                             | improves over presentation-only baseline without public-overlap leakage   | presentation prioritization generalizes             |
| activation_heldout         | patient/source-heldout T-cell activation labels                 | precision@k, recall at fixed FP, enrichment over prevalence | strict top-k remains enriched after WT/decoy-aware filtering              | immunogenicity prioritization, not clinical utility |
| paired_tcr_expert_ablation | paired-TCR-available rows only; epitope-heldout and TCR-heldout | delta AUPRC and abstention safety                           | TCR expert helps when paired TCR exists and abstains otherwise            | TCR-aware expert is useful in supported contexts    |
| structure_md_value_add     | matched mutant/WT/decoy structural audit rows                   | false-positive reduction and explanation rate               | MD/contact features reduce WT/decoy or unstable-interface false positives | structure layer improves triage/explanation         |
| prospective_locked_panel   | pre-registered candidate panel before assay results             | hit rate, FP burden, calibration drift                      | pre-registered top candidates exceed baseline hit rate under assay gates  | prospective decision-support evidence               |

## Claim Boundary

- Allowed now: CROSS-Neo-I model contract, feature/label schema, benchmark plan, and readiness ranking.
- Not allowed now: final calibrated immunogenicity probability, vaccine efficacy, or clinical response prediction.
- Required next: WT/decoy-controlled activation labels and source-heldout benchmarks.
