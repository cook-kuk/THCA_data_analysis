# CROSS-Neo impact portfolio v3

Generated: 2026-05-10T13:15:42

## Impact move

Single-score story에서 빠져나와 **two-lane portfolio**로 올렸다: clean BAR-Neo lead는 generalizable antigen discovery 축, TCR/MD lead는 mechanistic wetlab 축이다. 두 축은 같은 claim으로 섞지 않고 plate v3에서 같이 검증한다.

## Score recovery

- method-matrix stacker best OOF AUPRC: 0.323
- score booster source OOF AUPRC: 0.603
- score booster HLA OOF AUPRC: 0.643
- score booster mean AUROC: 0.909

## Lane counts

```json
{
  "SUPPORTING_RESERVE": 1957,
  "ASSAY_POSITIVE_CONTROL_OVERLAP_BLOCKED": 627,
  "TCR_STRUCTURE_DISAGREEMENT_AUDIT": 73,
  "LABEL_NOISE_RESCUE": 32,
  "CLEAN_BARNEO_LEAD": 17,
  "HARD_NEGATIVE_OR_FALSE_NEGATIVE_AUDIT": 7,
  "TCR_MD_STRUCTURAL_LEAD": 2
}
```

## Top clean BAR-Neo leads

| candidate_id   | peptide    | hla_allele_4digit   |   label | source_name   | leakage_risk_level   |   finetuned_score_booster_prob |   finetuned_experiment_priority_score |   impact_portfolio_score | claim_boundary_v3                                                               |
|:---------------|:-----------|:--------------------|--------:|:--------------|:---------------------|-------------------------------:|--------------------------------------:|-------------------------:|:--------------------------------------------------------------------------------|
| CNV0_02452     | LLDGFLATV  | HLA-A*02:01         |       1 | ITSNdb_main   | low                  |                          0.990 |                                 0.882 |                    0.468 | clean BAR-Neo experiment priority; patient metadata still blocks clinical claim |
| CNV0_02448     | ILDKVLVHL  | HLA-A*02:01         |       1 | ITSNdb_main   | low                  |                          0.988 |                                 0.876 |                    0.561 | clean BAR-Neo experiment priority; patient metadata still blocks clinical claim |
| CNV0_02504     | LLVDLAEEL  | HLA-A*02:01         |       1 | ITSNdb_main   | low                  |                          0.999 |                                 0.872 |                    0.452 | clean BAR-Neo experiment priority; patient metadata still blocks clinical claim |
| CNV0_02407     | KLMNIQQKL  | HLA-A*02:01         |       1 | ITSNdb_main   | low                  |                          0.999 |                                 0.871 |                    0.450 | clean BAR-Neo experiment priority; patient metadata still blocks clinical claim |
| CNV0_02450     | KELEGILLL  | HLA-B*44:03         |       1 | ITSNdb_main   | low                  |                          0.976 |                                 0.868 |                    0.471 | clean BAR-Neo experiment priority; patient metadata still blocks clinical claim |
| CNV0_02480     | ILDTAGKEEY | HLA-A*01:01         |       1 | ITSNdb_main   | medium               |                          0.995 |                                 0.864 |                    0.430 | clean BAR-Neo experiment priority; patient metadata still blocks clinical claim |
| CNV0_02410     | MLGEQLFPL  | HLA-A*02:01         |       1 | ITSNdb_main   | low                  |                          0.998 |                                 0.861 |                    0.441 | clean BAR-Neo experiment priority; patient metadata still blocks clinical claim |
| CNV0_02503     | RLSDFSEQL  | HLA-A*02:01         |       1 | ITSNdb_main   | low                  |                          0.997 |                                 0.859 |                    0.439 | clean BAR-Neo experiment priority; patient metadata still blocks clinical claim |
| CNV0_02424     | LADEAEVYL  | HLA-A*02:01         |       1 | ITSNdb_main   | low                  |                          0.982 |                                 0.855 |                    0.443 | clean BAR-Neo experiment priority; patient metadata still blocks clinical claim |
| CNV0_02409     | IILVAVPHV  | HLA-A*02:01         |       1 | ITSNdb_main   | low                  |                          0.995 |                                 0.852 |                    0.433 | clean BAR-Neo experiment priority; patient metadata still blocks clinical claim |
| CNV0_02408     | FLYNLLTRV  | HLA-A*02:01         |       1 | ITSNdb_main   | low                  |                          0.991 |                                 0.850 |                    0.436 | clean BAR-Neo experiment priority; patient metadata still blocks clinical claim |
| CNV0_02509     | SLLRSLENV  | HLA-A*02:01         |       1 | ITSNdb_main   | low                  |                          0.991 |                                 0.846 |                    0.430 | clean BAR-Neo experiment priority; patient metadata still blocks clinical claim |

## TCR/MD structural leads

| candidate_id   | peptide    | hla_allele_4digit   |   label | source_name   |   tcr_md_integrated_score |   tcr_augmented_score_mean | md_label                |   md_score |   control_readiness_score | claim_boundary_v3                                                        |
|:---------------|:-----------|:--------------------|--------:|:--------------|--------------------------:|---------------------------:|:------------------------|-----------:|--------------------------:|:-------------------------------------------------------------------------|
| CNV0_01132     | GADGVGKSAL | HLA-C*08:02         |       1 | NEPdb         |                     0.783 |                      0.960 | MD_MODERATE             |      0.557 |                     0.533 | TCR/MD diagnostic wetlab priority; not clean benchmark or clinical claim |
| CNV0_01151     | HMTEVVRHC  | HLA-A*02:01         |       1 | NEPdb         |                     0.777 |                      0.983 | MD_INSUFFICIENT_RUNTIME |      0.666 |                     0.478 | TCR/MD diagnostic wetlab priority; not clean benchmark or clinical claim |

## Label-noise rescue leads

| candidate_id   | peptide     | hla_allele_4digit   |   label | source_name   |   label_noise_priority_score |   finetuned_score_booster_prob | label_noise_bucket                    | recommended_next_step_v3                               |
|:---------------|:------------|:--------------------|--------:|:--------------|-----------------------------:|-------------------------------:|:--------------------------------------|:-------------------------------------------------------|
| CNV0_00711     | RREPPHLARNF | HLA-C*07:02         |       0 | CEDAR         |                        0.739 |                          0.922 | possible_false_negative_or_assay_miss | manual provenance audit + low-cost pMHC binding rescue |
| CNV0_00682     | NYRREPPHL   | HLA-C*07:02         |       0 | CEDAR         |                        0.736 |                          0.908 | possible_false_negative_or_assay_miss | manual provenance audit + low-cost pMHC binding rescue |
| CNV0_00674     | FKHSLSVSL   | HLA-C*07:02         |       0 | CEDAR         |                        0.734 |                          0.910 | possible_false_negative_or_assay_miss | manual provenance audit + low-cost pMHC binding rescue |
| CNV0_00006     | VKEDPKWEF   | HLA-C*07:02         |       0 | CEDAR         |                        0.732 |                          0.935 | possible_false_negative_or_assay_miss | manual provenance audit + low-cost pMHC binding rescue |
| CNV0_00027     | GNNDVKEDP   | HLA-C*07:02         |       0 | CEDAR         |                        0.731 |                          0.932 | possible_false_negative_or_assay_miss | manual provenance audit + low-cost pMHC binding rescue |
| CNV0_00058     | LLTQDLVQEK  | HLA-A*11:01         |       0 | CEDAR         |                        0.678 |                          0.916 | possible_false_negative_or_assay_miss | manual provenance audit + low-cost pMHC binding rescue |
| CNV0_00047     | NQFNLHELK   | HLA-A*11:01         |       0 | CEDAR         |                        0.678 |                          0.946 | possible_false_negative_or_assay_miss | manual provenance audit + low-cost pMHC binding rescue |
| CNV0_00067     | NYPSLTPQAF  | HLA-A*24:02         |       0 | CEDAR         |                        0.678 |                          0.941 | possible_false_negative_or_assay_miss | manual provenance audit + low-cost pMHC binding rescue |

## Plate v3

|   plate_v3_slot | plate_v3_block                | candidate_id   | peptide        | hla_allele_4digit   |   label | impact_lane                            |   impact_portfolio_score | plate_v3_assay_bundle                                       | claim_boundary_v3                                                               |
|----------------:|:------------------------------|:---------------|:---------------|:--------------------|--------:|:---------------------------------------|-------------------------:|:------------------------------------------------------------|:--------------------------------------------------------------------------------|
|               1 | A_clean_generalizable_antigen | CNV0_02452     | LLDGFLATV      | HLA-A*02:01         |       1 | CLEAN_BARNEO_LEAD                      |                    0.468 | mutant_vs_WT_decoy_pMHC_binding; TCR expansion if binder    | clean BAR-Neo experiment priority; patient metadata still blocks clinical claim |
|               2 | A_clean_generalizable_antigen | CNV0_02448     | ILDKVLVHL      | HLA-A*02:01         |       1 | CLEAN_BARNEO_LEAD                      |                    0.561 | mutant_vs_WT_decoy_pMHC_binding; TCR expansion if binder    | clean BAR-Neo experiment priority; patient metadata still blocks clinical claim |
|               3 | A_clean_generalizable_antigen | CNV0_02504     | LLVDLAEEL      | HLA-A*02:01         |       1 | CLEAN_BARNEO_LEAD                      |                    0.452 | mutant_vs_WT_decoy_pMHC_binding; TCR expansion if binder    | clean BAR-Neo experiment priority; patient metadata still blocks clinical claim |
|               4 | A_clean_generalizable_antigen | CNV0_02407     | KLMNIQQKL      | HLA-A*02:01         |       1 | CLEAN_BARNEO_LEAD                      |                    0.450 | mutant_vs_WT_decoy_pMHC_binding; TCR expansion if binder    | clean BAR-Neo experiment priority; patient metadata still blocks clinical claim |
|               5 | A_clean_generalizable_antigen | CNV0_02450     | KELEGILLL      | HLA-B*44:03         |       1 | CLEAN_BARNEO_LEAD                      |                    0.471 | mutant_vs_WT_decoy_pMHC_binding; TCR expansion if binder    | clean BAR-Neo experiment priority; patient metadata still blocks clinical claim |
|               6 | A_clean_generalizable_antigen | CNV0_02480     | ILDTAGKEEY     | HLA-A*01:01         |       1 | CLEAN_BARNEO_LEAD                      |                    0.430 | mutant_vs_WT_decoy_pMHC_binding; TCR expansion if binder    | clean BAR-Neo experiment priority; patient metadata still blocks clinical claim |
|               7 | B_TCR_MD_mechanistic          | CNV0_01132     | GADGVGKSAL     | HLA-C*08:02         |       1 | TCR_MD_STRUCTURAL_LEAD                 |                    0.579 | mutant_vs_WT_decoy_pMHC_binding + TCR_tetramer + activation | TCR/MD diagnostic wetlab priority; not clean benchmark or clinical claim        |
|               8 | B_TCR_MD_mechanistic          | CNV0_01151     | HMTEVVRHC      | HLA-A*02:01         |       1 | TCR_MD_STRUCTURAL_LEAD                 |                    0.451 | mutant_vs_WT_decoy_pMHC_binding + TCR_tetramer + activation | TCR/MD diagnostic wetlab priority; not clean benchmark or clinical claim        |
|               9 | C_false_negative_rescue       | CNV0_00711     | RREPPHLARNF    | HLA-C*07:02         |       0 | LABEL_NOISE_RESCUE                     |                    0.540 | provenance_audit + low_cost_pMHC_binding_rescue             | rescue/audit case; could expose false-negative label or assay mismatch          |
|              10 | C_false_negative_rescue       | CNV0_00682     | NYRREPPHL      | HLA-C*07:02         |       0 | LABEL_NOISE_RESCUE                     |                    0.534 | provenance_audit + low_cost_pMHC_binding_rescue             | rescue/audit case; could expose false-negative label or assay mismatch          |
|              11 | C_false_negative_rescue       | CNV0_00674     | FKHSLSVSL      | HLA-C*07:02         |       0 | LABEL_NOISE_RESCUE                     |                    0.534 | provenance_audit + low_cost_pMHC_binding_rescue             | rescue/audit case; could expose false-negative label or assay mismatch          |
|              12 | C_false_negative_rescue       | CNV0_00006     | VKEDPKWEF      | HLA-C*07:02         |       0 | LABEL_NOISE_RESCUE                     |                    0.545 | provenance_audit + low_cost_pMHC_binding_rescue             | rescue/audit case; could expose false-negative label or assay mismatch          |
|              13 | D_hard_negative_specificity   | CNV0_02397     | LPIQYEPVL      | HLA-B*35:03         |       0 | HARD_NEGATIVE_OR_FALSE_NEGATIVE_AUDIT  |                    0.462 | matched_negative_specificity_control                        | hard-negative or false-negative audit; use to sharpen specificity               |
|              14 | D_hard_negative_specificity   | CNV0_02400     | ETSKQVTRW      | HLA-A*25:01         |       0 | HARD_NEGATIVE_OR_FALSE_NEGATIVE_AUDIT  |                    0.455 | matched_negative_specificity_control                        | hard-negative or false-negative audit; use to sharpen specificity               |
|              15 | D_hard_negative_specificity   | CNV0_02420     | KLANPLPYT      | HLA-A*02:01         |       0 | HARD_NEGATIVE_OR_FALSE_NEGATIVE_AUDIT  |                    0.428 | matched_negative_specificity_control                        | hard-negative or false-negative audit; use to sharpen specificity               |
|              16 | D_hard_negative_specificity   | CNV0_02566     | KIFNFYPRK      | HLA-A*03:01         |       0 | HARD_NEGATIVE_OR_FALSE_NEGATIVE_AUDIT  |                    0.416 | matched_negative_specificity_control                        | hard-negative or false-negative audit; use to sharpen specificity               |
|              17 | E_assay_positive_control      | CNV0_00245     | ATSPHLESLLK    | HLA-A*11:01         |       1 | ASSAY_POSITIVE_CONTROL_OVERLAP_BLOCKED |                    0.512 | assay_QC_positive_control_only                              | positive assay control only; overlap/leakage blocks benchmark claim             |
|              18 | E_assay_positive_control      | CNV0_00487     | QELNELSAISL    | HLA-B*40:01         |       1 | ASSAY_POSITIVE_CONTROL_OVERLAP_BLOCKED |                    0.512 | assay_QC_positive_control_only                              | positive assay control only; overlap/leakage blocks benchmark claim             |
|              19 | E_assay_positive_control      | CNV0_00231     | AMFGKLMTI      | HLA-A*02:01         |       1 | ASSAY_POSITIVE_CONTROL_OVERLAP_BLOCKED |                    0.519 | assay_QC_positive_control_only                              | positive assay control only; overlap/leakage blocks benchmark claim             |
|              20 | E_assay_positive_control      | CNV0_00156     | ALPEVLAVIQV    | HLA-A*02:01         |       1 | ASSAY_POSITIVE_CONTROL_OVERLAP_BLOCKED |                    0.518 | assay_QC_positive_control_only                              | positive assay control only; overlap/leakage blocks benchmark claim             |
|              21 | F_TCR_structure_audit         | CNV0_01049     | LSSVSFFLV      | HLA-A*11:01         |       1 | TCR_STRUCTURE_DISAGREEMENT_AUDIT       |                    0.264 | TCR_specificity_review + repeat_MD_or_pMHC_binding_triage   | TCR/structure disagreement audit; no positive claim until resolved              |
|              22 | F_TCR_structure_audit         | CNV0_01383     | PVLSGCSPRCLLQG | HLA-B*15:01         |       0 | TCR_STRUCTURE_DISAGREEMENT_AUDIT       |                    0.177 | TCR_specificity_review + repeat_MD_or_pMHC_binding_triage   | TCR/structure disagreement audit; no positive claim until resolved              |
|              23 | F_TCR_structure_audit         | CNV0_01153     | DAMDAKKRRQQNV  | HLA-C*07:01         |       0 | TCR_STRUCTURE_DISAGREEMENT_AUDIT       |                    0.162 | TCR_specificity_review + repeat_MD_or_pMHC_binding_triage   | TCR/structure disagreement audit; no positive claim until resolved              |
|              24 | F_TCR_structure_audit         | CNV0_01154     | DAMDAKKRRQQNV  | HLA-C*01:02         |       0 | TCR_STRUCTURE_DISAGREEMENT_AUDIT       |                    0.167 | TCR_specificity_review + repeat_MD_or_pMHC_binding_triage   | TCR/structure disagreement audit; no positive claim until resolved              |

## Decision

1. Impact framing: not a single immunogenicity predictor; this is a controlled vaccine-candidate portfolio with orthogonal evidence lanes.
2. Wetlab plate v3 now has discovery leads, structural/TCR leads, rescue cases, hard negatives, and positive assay controls.
3. Claim boundary remains intact: no clinical vaccine-selection claim without patient metadata and external/wetlab validation.

## Files

- `impact_candidate_portfolio_v3.tsv`
- `wetlab_plate_v3_impact_design.tsv`
- `impact_portfolio_v3_summary.json`
- `figures/fig1_clean_cv_score_recovery.png`
- `figures/fig2_two_axis_impact_portfolio.png`
- `figures/fig3_plate_v3_composition.png`
- `figures/fig4_claim_boundary_lane_heatmap.png`
