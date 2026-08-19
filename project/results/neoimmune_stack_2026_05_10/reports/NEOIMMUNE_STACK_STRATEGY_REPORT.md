# NEOIMMUNE-STACK STRATEGY REPORT

**System names:** NeoImmune-Stack for business framing; CLEAN-Neo++ for paper framing.

**Core decision:** do not build a new foundation model yet. Integrate local algorithms, public frozen predictors, leakage-aware benchmarking, and patient-level top-N candidate ranking.

## Brutal executive decision
- **대박 가능성:** yes, as an integration/product strategy layer. It already unifies local algorithms, public comparators, leakage audit, and candidate handoff.
- **논문 claim:** upgraded. The integrated public tables now contain 86 real non-placeholder patient IDs (16,298 rows), so a retrospective public patient-level ranking scaffold is claimable.
- **Still not claimable:** prospective hospital utility, clinical vaccine efficacy, presentation without MS, or immunogenicity without T-cell assay.
- **Most defensible paper angle now:** leakage-aware patient-level immunogenicity ranking using public cohorts, plus a clear bridge to hospital-grade validation.
- **Must-have next data:** stronger per-patient expression/VAF/clonality, HLA typing QC, APM/HLA-LOH, and orthogonal presentation/immunogenicity assay labels.

## Scientific framing
`CLEAN-Neo++: leakage-aware integration of peptide-HLA presentation, TCR-visible immunogenicity, structure proxies, and quantum-kernel features for patient-level cancer vaccine candidate ranking.`

## Business framing
`NeoImmune-Stack: an AI operating layer for patient-specific cancer vaccine candidate prioritization.`

## Latest LLM strategy
- Use the latest available language model only as an **analysis infrastructure layer**: evidence summarization, candidate rationale drafting, contradiction detection, wet-lab question generation, and reviewer-response audit.
- Do not use LLM text output as a clean-science training feature.
- Do not let the LLM override leakage flags, labels, or wet-lab validation requirements.
- Add model name/config at runtime through `NEOIMMUNE_LLM_MODEL`; keep disabled by default until API/provider policy is selected.

## Latest GA/RL algorithm integration
The just-discovered `KG_GA_evolved_controller` has been integrated as a **production/experiment-priority branch**. It is deliberately not admitted to the clean science track because the controller includes public-predictor and product-value components.

| algorithm                |    n |   positives |    AUPRC |    AUROC |   Precision@34 |   Recall@34 |   patient_hit_rate@34 |   patient_recall@34 |   patients_evaluated | clean_track_allowed   |
|:-------------------------|-----:|------------:|---------:|---------:|---------------:|------------:|----------------------:|--------------------:|---------------------:|:----------------------|
| KG_GA_evolved_controller | 4430 |        1595 | 0.877279 | 0.9103   |       1        |   0.0213166 |              0.914634 |            0.98332  |                   82 | False                 |
| KG_GA_fixed_integrated   | 4430 |        1595 | 0.810954 | 0.866104 |       1        |   0.0213166 |              0.914634 |            0.978146 |                   82 | False                 |
| KG_GA_fixed_claimsafe    | 4430 |        1595 | 0.799992 | 0.839498 |       1        |   0.0213166 |              0.914634 |            0.975069 |                   82 | False                 |
| BigMHC_IM_inside_GA      | 4430 |        1595 | 0.644918 | 0.762679 |       0.794118 |   0.0169279 |              0.914634 |            0.981767 |                   82 | False                 |
| BigMHC_EL_inside_GA      | 4430 |        1595 | 0.528551 | 0.656181 |       0.647059 |   0.0137931 |              0.914634 |            0.96246  |                   82 | False                 |
| GA_TCR_expert_component  | 4430 |        1595 | 0.376094 | 0.512539 |       1        |   0.0213166 |              0.902439 |            0.93478  |                   82 | False                 |
| GA_MD_control_component  | 4430 |        1595 | 0.369675 | 0.507524 |       0.764706 |   0.0163009 |              0.902439 |            0.934658 |                   82 | False                 |
| GA_foreignness_component | 4430 |        1595 | 0.360045 | 0.5      |       0.470588 |   0.0100313 |              0.902439 |            0.934698 |                   82 | False                 |

- GA/RL patient top-34 rows: 879
- GA/RL patient coverage: 82 patients
- Mapping modes: `{'peptide_hla_source_patient_mapped': 856, 'direct_candidate_id': 23}`
- Claim-safe wording: GA/RL found a retrospective candidate-prioritization controller ready for prospective assay validation, not a clinically validated vaccine selection model.

## Competitor-facing imNEO public-source reconstruction
This is the narrow comparison requested against the public-source space visible from CG Invites/imNEO patent/presentation context. It is not the proprietary imNEO peptide list and not a reverse-engineered company model.

| model                      | track                    | clean_track_allowed   |   n |   positives |   real_patients | patient_level_claimable   |    AUPRC |    AUROC |   Precision@34 |   Recall@34 |
|:---------------------------|:-------------------------|:----------------------|----:|------------:|----------------:|:--------------------------|---------:|---------:|---------------:|------------:|
| ESM2 Bayesian              | local_clean_ablation     | True                  | 180 |          71 |               0 | False                     | 0.946193 | 0.962657 |       1        |    0.478873 |
| NMI W7B-ESM2-QK            | confident_clean_NMI      | True                  | 180 |          71 |               0 | False                     | 0.943834 | 0.964983 |       0.970588 |    0.464789 |
| NMI equal top5             | confident_clean_NMI      | True                  | 183 |          74 |               0 | False                     | 0.90417  | 0.938259 |       0.911765 |    0.418919 |
| W7B stacked                | local_clean_ablation     | True                  | 183 |          74 |               0 | False                     | 0.891753 | 0.927845 |       0.970588 |    0.445946 |
| NMI TCR-QK-Structure       | confident_clean_NMI      | True                  | 180 |          71 |               0 | False                     | 0.857525 | 0.931645 |       0.823529 |    0.394366 |
| NMI all-branch mean        | confident_clean_NMI      | True                  | 183 |          74 |               0 | False                     | 0.853879 | 0.932184 |       0.882353 |    0.405405 |
| BigMHC-IM                  | frozen_public_comparator | False                 | 183 |          74 |               0 | False                     | 0.825844 | 0.871188 |       0.852941 |    0.391892 |
| W7A QK only                | local_clean_ablation     | True                  | 183 |          74 |               0 | False                     | 0.787618 | 0.843913 |       0.911765 |    0.418919 |
| Structure LR               | local_clean_ablation     | True                  | 180 |          71 |               0 | False                     | 0.632964 | 0.71986  |       0.617647 |    0.295775 |
| MHCflurry presentation     | frozen_public_comparator | False                 | 183 |          74 |               0 | False                     | 0.56056  | 0.643566 |       0.617647 |    0.283784 |
| GP quantum                 | local_clean_ablation     | True                  | 183 |          74 |               0 | False                     | 0.544386 | 0.605381 |       0.558824 |    0.256757 |
| Wave8 TCR/SelfSim          | local_clean_ablation     | True                  | 183 |          74 |               0 | False                     | 0.543171 | 0.592363 |       0.676471 |    0.310811 |
| Wave8 TCR/SelfSim no-exact | local_clean_ablation     | True                  | 183 |          74 |               0 | False                     | 0.508372 | 0.579593 |       0.676471 |    0.310811 |
| PRIME                      | frozen_public_comparator | False                 | 183 |          74 |               0 | False                     | 0.502332 | 0.576246 |       0.470588 |    0.216216 |

Coverage of disclosed public sources in the local canonical table:
| source_name    | status                                       |   present_rows |   positive_labels |   nmi_locked_rows |   public_comparator_rows | claim_boundary                                                                 |
|:---------------|:---------------------------------------------|---------------:|------------------:|------------------:|-------------------------:|:-------------------------------------------------------------------------------|
| NEPdb          | present                                      |          16073 |               315 |               183 |                      486 | public-source reconstruction only; not CG Invites proprietary peptide list     |
| dbPepNeo       | present as dbPepNeo2                         |            523 |               523 |                 0 |                      467 | public-source reconstruction only; counts need not match patent training split |
| PRIME          | predictor score present, source table absent |              0 |                 0 |                 0 |                        0 | can compare PRIME score, cannot reconstruct PRIME training rows                |
| INeo-Epp       | absent                                       |              0 |                 0 |                 0 |                        0 | missing from current local canonical table                                     |
| TESLA          | present                                      |           1830 |                82 |                 0 |                      915 | patient-matched public benchmark; not the proprietary imNEO selected set       |
| IEDB t-cell DB | absent as explicit source                    |              0 |                 0 |                 0 |                        0 | do not silently relabel CEDAR as IEDB                                          |
| McPAS-TCR      | absent                                       |              0 |                 0 |                 0 |                        0 | requires additional ingestion                                                  |
| VDJdb          | absent                                       |              0 |                 0 |                 0 |                        0 | requires additional ingestion                                                  |

Brutal boundary: NMI beats frozen public comparators on the current NEPdb common head-to-head subset, but patient-level top-N is not claimable there because real patient IDs are absent/collapsed. TESLA and dbPepNeo need NMI locked-branch coverage before broad competitor claims.

## Source transfer gap
| source_family   |   rows |   positives | nmi_best_model      |   nmi_best_auprc | public_best_model          |   public_best_auprc |   delta_nmi_minus_public | patient_claimable   |
|:----------------|-------:|------------:|:--------------------|-----------------:|:---------------------------|--------------------:|-------------------------:|:--------------------|
| NEPdb           |  16073 |         315 | NMI_w7b_esm2_qk     |         0.943834 | BigMHC_IM                  |            0.762643 |                0.181191  | True                |
| CEDAR/IMPROVE   |   3345 |        2932 | nan                 |       nan        | MHCflurry_2.0_presentation |            0.972641 |              nan         | False               |
| TESLA           |   1830 |          82 | nan                 |       nan        | BigMHC_IM                  |            0.29631  |              nan         | True                |
| dbPepNeo        |    523 |         523 | nan                 |       nan        | nan                        |          nan        |              nan         | False               |
| ITSNdb          |    319 |         136 | NMI_all_branch_mean |         0.726581 | MHCflurry_2.0_presentation |            0.657881 |                0.0687006 | False               |

Readout: NEPdb is already positive for NMI vs public comparators; TESLA/dbPepNeo remain coverage-limited and must not be oversold.

## Source rescue queue
| dataset_source                           |   rows |   positives |   top_priority |   median_priority |   mean_clean |   mean_public_best |   mean_delta |
|:-----------------------------------------|-------:|------------:|---------------:|------------------:|-------------:|-------------------:|-------------:|
| CEDAR                                    |    863 |         816 |       0.8539   |          0.660238 |     0.874479 |           0.837567 |    0.0434579 |
| dbPepNeo2_MHCI                           |    344 |         344 |       0.83892  |          0.656049 |     0.892926 |           0.745615 |    0.147311  |
| TESLA_mmc4                               |    273 |          63 |       0.752232 |          0.6237   |     0.938864 |           0.923305 |    0.0156216 |
| TESLA_mmc7_validation                    |     41 |           8 |       0.685411 |          0.507685 |     0.888607 |           0.90787  |   -0.0189757 |
| dbPepNeo2_MHCII                          |      9 |           9 |       0.609821 |          0.525121 |     0.823037 |         nan        |  nan         |
| IMPROVE_Neoepitopes_CEDAR_benchmark_data |     40 |          40 |       0.487411 |          0.445179 |     0.736836 |         nan        |  nan         |

Interpretation: use the rescue queue for wet-lab discussion and failure-case triage. CEDAR and dbPepNeo rows dominate the current priority list because they have the strongest clean score and disagreement signal, while TESLA_mmc7_validation remains weaker and should stay in the control bucket unless new evidence appears.

## Collaborator handoff
Positive rows:
| candidate_id                                  | dataset_source                           | peptide_mut     | hla_allele          |   label_immunogenicity |   public_best |   clean_science_score |   production_stack_score |   local_rescue_delta_vs_public |   model_disagreement_score |   rescue_priority |
|:----------------------------------------------|:-----------------------------------------|:----------------|:--------------------|-----------------------:|--------------:|----------------------:|-------------------------:|-------------------------------:|---------------------------:|------------------:|
| CNV0_00446                                    | CEDAR                                    | QALILKIA        | HLA-B*15:01         |                      1 |    0.184108   |              0.994486 |                 1        |                      0.810378  |                   0.295371 |          0.8539   |
| dbPepNeo2_MHC-I_HC_neoantigens_412            | dbPepNeo2_MHCI                           | CMGGMNLR        | HLA-A*31:01         |                      1 |    0.293273   |              0.971436 |                 0.999077 |                      0.678162  |                   0.442174 |          0.83892  |
| CNV0_00216                                    | CEDAR                                    | AGRKLALKTIDW    | HLA-A*33:03         |                      1 |    0.00502039 |              0.904681 |                 1        |                      0.899661  |                   0.342682 |          0.838441 |
| CNV0_00813                                    | CEDAR                                    | FHVEEEGK        | HLA-A*02:01         |                      1 |    0.0755507  |              0.935973 |                 1        |                      0.860423  |                   0.299589 |          0.838211 |
| dbPepNeo2_MHC-I_Fusion_neoantigens_2          | dbPepNeo2_MHCI                           | DKESEEEVS       | HLA-C*04:01         |                      1 |    0.12777    |              0.942842 |                 0.986348 |                      0.815072  |                   0.353584 |          0.837601 |
| dbPepNeo2_MHC-I_HC_neoantigens_621            | dbPepNeo2_MHCI                           | FPKKIQMLA       | HLA-B*58:01         |                      1 |    0.319245   |              0.968642 |                 0.996544 |                      0.649397  |                   0.346807 |          0.817098 |
| CNV0_01495                                    | TESLA_mmc4                               | VRINTARPV       | HLA-C*06:02         |                      1 |    0.729424   |              0.982384 |                 0.999863 |                      0.25296   |                   0.257634 |          0.731283 |
| CNV0_01560                                    | TESLA_mmc4                               | KLRDEISLAK      | HLA-A*03:01         |                      1 |    0.912019   |              0.983159 |                 0.999919 |                      0.0711401 |                   0.259594 |          0.695573 |
| CNV0_01787                                    | TESLA_mmc4                               | FLGSLLILV       | HLA-A*02:01         |                      1 |    0.966542   |              0.987543 |                 0.999654 |                      0.0210018 |                   0.262912 |          0.687962 |
| CNV0_02098                                    | TESLA_mmc7_validation                    | YLDGKVVDY       | HLA-A*01:01         |                      1 |    0.988214   |              0.936935 |                 0.99999  |                     -0.0512795 |                   0.274638 |          0.662814 |
| 301f5e64ca51db57                              | TESLA_mmc7_validation                    | YLDGKVVDY       | HLA-A*01:01         |                      1 |  nan          |              0.936935 |                 0.99965  |                    nan         |                   0.271037 |          0.662206 |
| CNV0_02146                                    | TESLA_mmc7_validation                    | ALDPHSGHFV      | HLA-A*02:01         |                      1 |    0.966209   |              0.928341 |                 0.998964 |                     -0.0378679 |                   0.275569 |          0.658882 |
| dbPepNeo2_MHC-II_HC_neoantigens_29            | dbPepNeo2_MHCII                          | FGLLGNILLVI     | DPA1*0202/DPB1*0301 |                      1 |  nan          |              0.943054 |                 0.927233 |                    nan         |                 nan        |          0.609821 |
| dbPepNeo2_MHC-II_HC_neoantigens_35            | dbPepNeo2_MHCII                          | MVMGVLKQAFDVLVL | HLA-DRB1*07:01      |                      1 |  nan          |              0.899896 |                 0.919731 |                    nan         |                 nan        |          0.5889   |
| dbPepNeo2_MHC-II_HC_neoantigens_27            | dbPepNeo2_MHCII                          | IEILRNLYHEIIV   | HLA-DRB4*01:01      |                      1 |  nan          |              0.877656 |                 0.811589 |                    nan         |                 nan        |          0.557263 |
| IMPROVE_Neoepitopes_CEDAR_benchmark_data_1581 | IMPROVE_Neoepitopes_CEDAR_benchmark_data | HLABDEERIPVL    | B18:01	DEERIPVL	0                     |                      1 |  nan          |              0.800482 |                 0.63597  |                    nan         |                 nan        |          0.487411 |
| IMPROVE_Neoepitopes_CEDAR_benchmark_data_2435 | IMPROVE_Neoepitopes_CEDAR_benchmark_data | HLABEPGDGGIY    | B35:01	EPGDGGIY	0                     |                      1 |  nan          |              0.759266 |                 0.666696 |                    nan         |                 nan        |          0.475009 |
| IMPROVE_Neoepitopes_CEDAR_benchmark_data_1618 | IMPROVE_Neoepitopes_CEDAR_benchmark_data | HLABVENQKHSL    | B08:01	VENQKHSL	0                     |                      1 |  nan          |              0.767702 |                 0.630931 |                    nan         |                 nan        |          0.471652 |

Negative controls:
| candidate_id   | dataset_source        | peptide_mut   | hla_allele   |   label_immunogenicity |   public_best |   clean_science_score |   production_stack_score |   local_rescue_delta_vs_public |   model_disagreement_score |   rescue_priority |
|:---------------|:----------------------|:--------------|:-------------|-----------------------:|--------------:|----------------------:|-------------------------:|-------------------------------:|---------------------------:|------------------:|
| CNV0_00496     | CEDAR                 | AIKTSPKA      | HLA-A*02:01  |                      0 |      0.184014 |              0.963627 |                 1        |                      0.779614  |                   0.278336 |          0.831305 |
| CNV0_00011     | CEDAR                 | VPLCDLLLEM    | HLA-A*02:01  |                      0 |      0.465682 |              0.979087 |                 1        |                      0.513405  |                   0.267322 |          0.783369 |
| CNV0_01799     | TESLA_mmc4            | FMENIMEM      | HLA-A*02:01  |                      0 |      0.648125 |              0.989387 |                 0.999986 |                      0.341262  |                   0.258388 |          0.752232 |
| CNV0_01645     | TESLA_mmc4            | RRHIEIRDK     | HLA-B*27:05  |                      0 |      0.778379 |              0.988231 |                 0.999718 |                      0.209852  |                   0.293993 |          0.730717 |
| CNV0_02086     | TESLA_mmc7_validation | ILTDIEEKV     | HLA-A*02:01  |                      0 |      0.993748 |              0.982701 |                 0.999995 |                     -0.0110463 |                   0.287975 |          0.685411 |
| CNV0_02189     | TESLA_mmc7_validation | ALSPVIPHI     | HLA-A*02:01  |                      0 |      1        |              0.982648 |                 0.991554 |                     -0.0173525 |                   0.292727 |          0.684411 |

## 1. What existing public models are usable now?
| model_name                 | type                                  | runnable   | license_note                    | role                                   |
|:---------------------------|:--------------------------------------|:-----------|:--------------------------------|:---------------------------------------|
| NetMHCpan_4.1_EL           | HLA ligand elution predictor          | False      | license/manual install          | presentation comparator/frozen feature |
| NetMHCpan_4.1_BA           | binding affinity predictor            | False      | license/manual install          | binding comparator/frozen feature      |
| MHCflurry_2.0_presentation | presentation predictor                | True       | public package/artifact         | presentation comparator/frozen feature |
| MHCflurry_2.0_affinity     | binding predictor                     | True       | public package/artifact         | binding comparator/frozen feature      |
| BigMHC_EL                  | EL predictor                          | False      | external model/artifact pending | presentation comparator/frozen feature |
| BigMHC_IM                  | immunogenicity predictor              | True       | public model/artifact           | public comparator/frozen feature       |
| PRIME                      | immunogenicity/presentation predictor | True       | license sensitive/artifact      | public comparator/frozen feature       |
| PRIME2.1                   | updated PRIME family                  | False      | license/install check required  | pending comparator                     |

Existing artifacts were found for several public comparators, especially MHCflurry, BigMHC-IM, PRIME, and some NetMHCpan-formatted outputs. NetMHCpan/PRIME family tools remain license-sensitive and must be treated as frozen comparators/features, not clean training features.

## 2. Which local algorithms already exist in this repo?
| model_name                        | type                                | role                                | runnable   | clean_track_allowed   | production_track_allowed   |
|:----------------------------------|:------------------------------------|:------------------------------------|:-----------|:----------------------|:---------------------------|
| Structure_LR                      | structure_proxy/logistic            | clean baseline                      | True       | True                  | True                       |
| Wave8_TCR_SelfSim_full            | TCR/self-similarity                 | immunogenicity branch               | True       | True                  | True                       |
| ESM2_Bayesian                     | protein language embedding Bayesian | sequence representation branch      | True       | True                  | True                       |
| GP_quantum                        | Gaussian process/quantum features   | experimental branch                 | True       | True                  | True                       |
| VQC                               | variational quantum classifier      | experimental branch                 | False      | True                  | True                       |
| W7A_QK_only                       | quantum-kernel-only                 | clean algorithm branch              | True       | True                  | True                       |
| W7A_full                          | multi-feature clean model           | clean algorithm branch              | True       | True                  | True                       |
| W7B_stacked                       | stacked local ensemble              | clean stack if trained leakage-safe | True       | True                  | True                       |
| quantum_kernel_no_anchor_gamma1.0 | quantum kernel                      | must survive source-heldout         | True       | True                  | True                       |
| quantum_kernel_compact_gamma0.5   | quantum kernel                      | compact comparator                  | True       | True                  | True                       |
| Stack_mean                        | score stack                         | simple clean ensemble               | True       | True                  | True                       |
| Stack_median                      | score stack                         | robust clean ensemble               | True       | True                  | True                       |
| Stack_LR                          | logistic stack                      | clean stack if leakage-safe         | True       | True                  | True                       |
| BAR-Neo                           | patient-gated score                 | business/practical comparator       | True       | False                 | True                       |
| BAR-Neo-X                         | reviewer-aware ensemble             | production support                  | True       | False                 | True                       |

## 3. Which local algorithm contributes most under strict no-leakage evaluation?
| model                      | split                   |   n |   positives |    AUPRC |    AUROC |   patient_hit_rate@20 |   patient_recall@20 |
|:---------------------------|:------------------------|----:|------------:|---------:|---------:|----------------------:|--------------------:|
| W7B_stacked                | existing_local_artifact | 319 |         136 | 0.718157 | 0.782847 |              0.747126 |            0.744385 |
| ESM2_Bayesian              | existing_local_artifact | 311 |         130 | 0.704536 | 0.783638 |              0.747126 |            0.744379 |
| W7A_QK_only                | existing_local_artifact | 319 |         136 | 0.676748 | 0.721593 |              0.747126 |            0.744389 |
| W7A_full                   | existing_local_artifact | 319 |         136 | 0.66154  | 0.744556 |              0.747126 |            0.744385 |
| Structure_LR               | existing_local_artifact | 311 |         130 | 0.612358 | 0.682363 |              0.747126 |            0.744385 |
| Stack_LR_inmaster_E3a      | existing_local_artifact | 319 |         136 | 0.611342 | 0.749297 |              0.747126 |            0.744362 |
| Stack_mean_E1              | existing_local_artifact | 319 |         136 | 0.599588 | 0.679866 |              0.747126 |            0.744379 |
| Wave8_TCR_SelfSim_no_exact | existing_local_artifact | 319 |         136 | 0.575483 | 0.602961 |              0.747126 |            0.744385 |
| Wave8_TCR_motif_only       | existing_local_artifact | 319 |         136 | 0.564921 | 0.594403 |              0.747126 |            0.744379 |
| Wave8_TCR_SelfSim_full     | existing_local_artifact | 319 |         136 | 0.563012 | 0.604649 |              0.747126 |            0.744376 |

Brutal boundary: the integrated run can rank available local artifacts, but a publishable clean-algorithm claim requires the row to be explicitly strict/no-overlap/source-heldout and free of public predictor scores.

## 4. Does quantum_kernel_no_anchor_gamma1.0 survive source-heldout validation?
Best integrated quantum-labeled row is `W7A_QK_only` with AUPRC=0.6767484970911705; this is still artifact-level evidence unless the row's split is strict/source-heldout.

Historical JSON summaries show strong random/CV signal can degrade under leave-source-out conditions; therefore source-heldout survival is a must-pass gate, not an assumption.

## 5. Does Wave8_TCR_SelfSim_full improve immunogenicity ranking beyond HLA presentation?
Wave8/TCR-self-similarity evidence exists as `Wave8_TCR_SelfSim_no_exact` with AUPRC=0.5754831854430665; interpret only within its stated split.

This is the most important scientific branch because vaccine ranking bottleneck is immunogenicity, not binding alone. But it must be shown against presentation-only comparators under patient/source leakage controls.

## 6. Does production stacking beat BigMHC-IM / PRIME / MHCflurry / NetMHCpan?
| model                                             | split             |    n |   positives |    AUPRC |    AUROC |   Precision@20 |   Recall@20 |   patient_hit_rate@20 |
|:--------------------------------------------------|:------------------|-----:|------------:|---------:|---------:|---------------:|------------:|----------------------:|
| BAR-Neo                                           | existing_artifact | 2715 |        1179 | 0.977163 | 0.979698 |           1    |   0.0169635 |                     1 |
| BAR-Neo_patient_gated                             | existing_artifact | 2715 |        1179 | 0.97663  | 0.97946  |           1    |   0.0169635 |                     1 |
| BAR-Neo_confidence                                | existing_artifact | 2715 |        1179 | 0.961523 | 0.967173 |           1    |   0.0169635 |                     1 |
| BAR-Neo-X                                         | existing_artifact | 2715 |        1179 | 0.940181 | 0.94522  |           0.95 |   0.0161154 |                     1 |
| stress_guarded_discovery                          | existing_artifact | 2532 |        1105 | 0.901213 | 0.924117 |           0.95 |   0.0171946 |                     1 |
| source_balanced_rf_train_prior_calibrated         | existing_artifact | 2396 |        1043 | 0.894791 | 0.913431 |           1    |   0.0191755 |                     1 |
| source_balanced_plus_pu_rf_train_prior_calibrated | existing_artifact | 2396 |        1043 | 0.887599 | 0.911438 |           0.95 |   0.0182167 |                     1 |
| pu_weighted_rf_train_prior_calibrated             | existing_artifact | 2396 |        1043 | 0.885673 | 0.905549 |           0.95 |   0.0182167 |                     1 |
| clean_contextual_bma                              | existing_artifact | 2532 |        1105 | 0.884398 | 0.928542 |           1    |   0.0180995 |                     1 |
| RF_biophys                                        | existing_artifact | 2987 |        1561 | 0.853225 | 0.824752 |           1    |   0.0128123 |                     1 |
| RF_biophys                                        | existing_artifact | 2987 |        1561 | 0.853225 | 0.824752 |           1    |   0.0128123 |                     1 |
| stress_guarded_claim_safe                         | existing_artifact | 2532 |        1105 | 0.825612 | 0.889418 |           0.8  |   0.0144796 |                     1 |

Boundary: if production stack wins, the claim is practical integration, not a new clean predictor. Public predictor scores are frozen features and may carry training-corpus contamination risk.

## 7. Which patients have top-20 rescued candidates?
|   patient_id | candidate_id     | dataset_source   |   gene | peptide_mut   | hla_allele   |   production_stack_score |   clean_science_score |   model_disagreement_score |
|-------------:|:-----------------|:-----------------|-------:|:--------------|:-------------|-------------------------:|----------------------:|---------------------------:|
|            1 | 17328e005d288a20 | TESLA_mmc4       |    nan | LLYDHPLKV     | HLA-A*02:01  |                 1        |              0.983066 |                   0.257792 |
|            1 | 4b73849004bae360 | TESLA_mmc4       |    nan | LLAPLIATL     | HLA-A*02:01  |                 0.999999 |              0.994297 |                   0.247534 |
|            1 | 97be0faff546a239 | TESLA_mmc4       |    nan | GMLTNIWTV     | HLA-A*02:01  |                 0.999999 |              0.958213 |                   0.269693 |
|            1 | 5353c42c81b9040e | TESLA_mmc4       |    nan | FLDPDLTNI     | HLA-A*02:01  |                 0.999999 |              0.962062 |                   0.230146 |
|            1 | 0319abbf1a657b93 | TESLA_mmc4       |    nan | DTIDVSKLNR    | HLA-A*68:01  |                 0.999999 |              0.96024  |                   0.263162 |
|            1 | e6f8e72178c66507 | TESLA_mmc4       |    nan | FLGSLLILV     | HLA-A*02:01  |                 0.999994 |              0.987543 |                   0.253662 |
|            1 | 2fade8382f40ff2a | TESLA_mmc4       |    nan | KLINSQINL     | HLA-A*02:01  |                 0.999993 |              0.970922 |                   0.228654 |
|            1 | edc2ce155c59f22e | TESLA_mmc4       |    nan | EVCAPVIKR     | HLA-A*68:01  |                 0.999981 |              0.987649 |                   0.226161 |
|            1 | 949afa345f7a981b | TESLA_mmc4       |    nan | EAFAASAVGSR   | HLA-A*68:01  |                 0.999948 |              0.94644  |                   0.234253 |
|            1 | c69263bdb3df475a | TESLA_mmc4       |    nan | YLYHRVDVI     | HLA-A*02:01  |                 0.999872 |              0.973108 |                   0.210394 |
|            1 | a7c0ccc5322680dd | TESLA_mmc4       |    nan | EVVLFHPNYNINR | HLA-A*68:01  |                 0.999871 |              0.746795 |                   0.248777 |
|            1 | 8f72bfee1aa42358 | TESLA_mmc4       |    nan | KAWENFPNV     | HLA-A*02:01  |                 0.99986  |              0.920202 |                   0.213983 |
|            1 | 44bbf3a2a2d4edf2 | TESLA_mmc4       |    nan | FMENIMEM      | HLA-A*02:01  |                 0.999826 |              0.989387 |                   0.192987 |
|            1 | e6386c2e81c132aa | TESLA_mmc4       |    nan | EIIPQCIAR     | HLA-A*68:01  |                 0.999798 |              0.982384 |                   0.209872 |
|            1 | 5fa19f7eade9208c | TESLA_mmc4       |    nan | RVYDALNLL     | HLA-A*02:01  |                 0.999495 |              0.984764 |                   0.184425 |
|            1 | 2fb2f814a086c496 | TESLA_mmc4       |    nan | SLHDLTDGV     | HLA-A*02:01  |                 0.999368 |              0.977277 |                   0.195231 |
|            1 | 2ae5b59de49726f0 | TESLA_mmc4       |    nan | ALSDGPIDAV    | HLA-A*02:01  |                 0.999223 |              0.980666 |                   0.198775 |
|            1 | 20395f83a0eab2af | TESLA_mmc4       |    nan | FLGSLLILVV    | HLA-A*02:01  |                 0.999177 |              0.98379  |                   0.213646 |
|            1 | dc9110087349078b | TESLA_mmc4       |    nan | FVLASLCLYV    | HLA-A*02:01  |                 0.999112 |              0.957797 |                   0.226167 |
|            1 | 5dda089ac583b061 | TESLA_mmc4       |    nan | FLNCDIMLGV    | HLA-A*02:01  |                 0.999063 |              0.964898 |                   0.202158 |

## 8. Which candidates should be discussed with wet-lab collaborators?
Prioritize top-20 candidates with high production score, adequate expression/VAF, low leakage risk, model agreement, and a positive local/TCR/structure branch. Exclude or down-rank HLA-LOH/APM-loss and high-leakage candidates unless the goal is methods debugging.

## 9. What is the next minimal experiment to validate?
- Assemble a small patient-level candidate set with matched tumor expression, VAF/clonality, HLA typing, and APM/HLA-LOH metadata.
- Run frozen public predictors and local CLEAN-Neo++ branches.
- Select top 20 plus disagreement controls per patient.
- Validate presentation with MS where feasible and immunogenicity with T-cell assay where labels are claimed.

## 10. What should be the paper/사업 framing?
- Paper: leakage-aware patient-level immunogenicity ranking; endpoint is patient-level Recall@20/hit rate, secondary AUPRC under strict no-overlap/source-heldout splits.
- Business: operating layer that standardizes public tools, local algorithms, leakage audit, candidate explanation, and wet-lab handoff.

## Leakage summary
| leakage_control                 |   n_flagged |   denominator |   fraction |
|:--------------------------------|------------:|--------------:|-----------:|
| existing_overlap_flags          |        2327 |         22090 |   0.105342 |
| source_protein_window_available |           0 |         22090 |   0        |

## Data scale
- Canonical candidate rows: 22,090
- Sources: 9
- Patients: 87
- Real non-placeholder patient rows: 16,298
- Real non-placeholder unique patients: 86

## Existing public leakage context from hub
- Hub total rows: 119278
- Hub leakage-free rows: 19817
- Hub test-set safety counts: `{'TRAINING_OVERLAP': 62376, 'PREDICTED_ONLY': 30983, 'HELD_OUT': 16440, 'PARTIAL_OVERLAP': 4937, 'EXTERNAL_TEST': 3377, 'UNKNOWN': 1156, 'DEMO_ONLY': 9}`

## Existing representation context
- QX 5-fold CV summary: `{'mean_AUROC': 0.8540767891531009, 'mean_AUPRC': 0.898906717776924, 'std_AUROC': 0.007738139324958291, 'folds': [{'fold': 1, 'AUROC': 0.8618767183286631, 'AUPRC': 0.9012989943111807, 'brier': 0.15121999413786275, 'n_test': 898}, {'fold': 2, 'AUROC': 0.8449506862760618, 'AUPRC': 0.9055882585654622, 'brier': 0.15926806518027978, 'n_test': 898}, {'fold': 3, 'AUROC': 0.855828772850928, 'AUPRC': 0.8956529611255388, 'brier': 0.15415601642630722, 'n_test': 897}, {'fold': 4, 'AUROC': 0.8451429539223677, 'AUPRC': 0.8848737846569634, 'brier': 0.15748067658115947, 'n_test': 897}, {'fold': 5, 'AUROC': 0.8625848143874841, 'AUPRC': 0.9071195902254744, 'brier': 0.14792490565160335, 'n_test': 897}]}`
- ESM2 full summary: `{'n': 4487, 'positive_rate': 0.6229106307109428, 'mean_AUROC': 0.8317490963264499, 'std_AUROC': 0.019295472050992638, 'folds': [{'fold': 1, 'AUROC': 0.8511379887177376, 'AUPRC': 0.9016335760565364}, {'fold': 2, 'AUROC': 0.8113149798681801, 'AUPRC': 0.8703042492770271}, {'fold': 3, 'AUROC': 0.8421367403753532, 'AUPRC': 0.8815890641310129}, {'fold': 4, 'AUROC': 0.8056705232293508, 'AUPRC': 0.850223101807252}, {'fold': 5, 'AUROC': 0.8484852494416275, 'AUPRC': 0.8941027793281718}]}`

## Command examples
```bash
OUT=project/results/neoimmune_stack_2026_05_10
source /home/seungho/personal/THCA_data_analysis/.venv/bin/activate
python project/scripts/neoimmune_stack/00_repo_audit.py --outdir $OUT
python project/scripts/neoimmune_stack/01_build_model_registry.py --outdir $OUT
python project/scripts/neoimmune_stack/02_build_canonical_candidate_table.py --outdir $OUT
python project/scripts/neoimmune_stack/03_run_external_adapters.py --outdir $OUT
python project/scripts/neoimmune_stack/04_collect_local_model_outputs.py --outdir $OUT
python project/scripts/neoimmune_stack/05_leakage_audit.py --outdir $OUT
python project/scripts/neoimmune_stack/06_train_clean_ranker.py --outdir $OUT
python project/scripts/neoimmune_stack/07_train_production_stack.py --outdir $OUT
python project/scripts/neoimmune_stack/08_evaluate_patient_topn.py --outdir $OUT
python project/scripts/neoimmune_stack/09_failure_case_audit.py --outdir $OUT
python project/scripts/neoimmune_stack/17_apply_ga_rl_algorithm.py --outdir $OUT
python project/scripts/neoimmune_stack/18_train_nmi_clean_method.py --outdir $OUT
python project/scripts/neoimmune_stack/19_build_imneo_public_reconstruction.py --outdir $OUT
python project/scripts/neoimmune_stack/20_source_transfer_gap_audit.py --outdir $OUT
python project/scripts/neoimmune_stack/21_build_source_rescue_queue.py --outdir $OUT
python project/scripts/neoimmune_stack/22_build_collaborator_handoff.py --outdir $OUT
python project/scripts/neoimmune_stack/10_generate_strategy_report.py --outdir $OUT
```

## Claim boundaries
- No clinical vaccine efficacy claim.
- No antigen-presentation claim without MS immunopeptidomics or equivalent evidence.
- No immunogenicity claim without T-cell assay labels.
- HLA binding is not immunogenicity.
- Public predictor features are excluded from clean-science training.
- Leakage problems are surfaced, not hidden.
