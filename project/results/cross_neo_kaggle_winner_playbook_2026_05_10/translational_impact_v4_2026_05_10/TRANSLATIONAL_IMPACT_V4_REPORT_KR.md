# CROSS-Neo translational impact v4

Generated: 2026-05-10T13:23:08

## Impact move

v4 converts the portfolio into a **value-of-information wetlab plan**. The key question is no longer just which peptide scores highest; it is which assay result unlocks the strongest claim while preserving overlap, patient-context, and clinical-use boundaries.

## Headline metrics

- Plate v4 candidates: 24
- 96-well assay map rows: 96
- Mean plate claim-unlock score: 0.548
- Mean expected information gain: 0.489
- Score recovery vs method-matrix stacker: 1.93x

## Plate v4

|   plate_v4_slot | plate_v4_arm          | candidate_id   | peptide     | hla_allele_4digit   |   label | impact_lane                            |   v4_claim_unlock_score |   v4_expected_information_gain | v4_claim_unlock                        | plate_v4_primary_readout                   | claim_boundary_v3                                                               |
|----------------:|:----------------------|:---------------|:------------|:--------------------|--------:|:---------------------------------------|------------------------:|-------------------------------:|:---------------------------------------|:-------------------------------------------|:--------------------------------------------------------------------------------|
|               1 | A_clean_discovery     | CNV0_02448     | ILDKVLVHL   | HLA-A*02:01         |       1 | CLEAN_BARNEO_LEAD                      |                   0.603 |                          0.505 | generalizable antigen discovery        | mutant-vs-WT/decoy pMHC binding            | clean BAR-Neo experiment priority; patient metadata still blocks clinical claim |
|               2 | A_clean_discovery     | CNV0_02708     | ALDPLLLRI   | HLA-A*02:01         |       1 | CLEAN_BARNEO_LEAD                      |                   0.602 |                          0.579 | generalizable antigen discovery        | mutant-vs-WT/decoy pMHC binding            | clean BAR-Neo experiment priority; patient metadata still blocks clinical claim |
|               3 | A_clean_discovery     | CNV0_02709     | ALPVALPSL   | HLA-A*02:01         |       1 | CLEAN_BARNEO_LEAD                      |                   0.579 |                          0.574 | generalizable antigen discovery        | mutant-vs-WT/decoy pMHC binding            | clean BAR-Neo experiment priority; patient metadata still blocks clinical claim |
|               4 | A_clean_discovery     | CNV0_02450     | KELEGILLL   | HLA-B*44:03         |       1 | CLEAN_BARNEO_LEAD                      |                   0.568 |                          0.512 | generalizable antigen discovery        | mutant-vs-WT/decoy pMHC binding            | clean BAR-Neo experiment priority; patient metadata still blocks clinical claim |
|               5 | A_clean_discovery     | CNV0_02460     | GIVEGLITTV  | HLA-A*02:01         |       1 | CLEAN_BARNEO_LEAD                      |                   0.568 |                          0.556 | generalizable antigen discovery        | mutant-vs-WT/decoy pMHC binding            | clean BAR-Neo experiment priority; patient metadata still blocks clinical claim |
|               6 | A_clean_discovery     | CNV0_02452     | LLDGFLATV   | HLA-A*02:01         |       1 | CLEAN_BARNEO_LEAD                      |                   0.565 |                          0.499 | generalizable antigen discovery        | mutant-vs-WT/decoy pMHC binding            | clean BAR-Neo experiment priority; patient metadata still blocks clinical claim |
|               7 | B_mechanism_TCR_MD    | CNV0_01132     | GADGVGKSAL  | HLA-C*08:02         |       1 | TCR_MD_STRUCTURAL_LEAD                 |                   0.746 |                          0.781 | mechanistic TCR/MD support             | pMHC binding + TCR tetramer/activation     | TCR/MD diagnostic wetlab priority; not clean benchmark or clinical claim        |
|               8 | B_mechanism_TCR_MD    | CNV0_01151     | HMTEVVRHC   | HLA-A*02:01         |       1 | TCR_MD_STRUCTURAL_LEAD                 |                   0.695 |                          0.773 | mechanistic TCR/MD support             | pMHC binding + TCR tetramer/activation     | TCR/MD diagnostic wetlab priority; not clean benchmark or clinical claim        |
|               9 | C_label_rescue        | CNV0_00027     | GNNDVKEDP   | HLA-C*07:02         |       0 | LABEL_NOISE_RESCUE                     |                   0.580 |                          0.509 | false-negative rescue / assay mismatch | pMHC binding rescue after provenance audit | rescue/audit case; could expose false-negative label or assay mismatch          |
|              10 | C_label_rescue        | CNV0_00006     | VKEDPKWEF   | HLA-C*07:02         |       0 | LABEL_NOISE_RESCUE                     |                   0.580 |                          0.508 | false-negative rescue / assay mismatch | pMHC binding rescue after provenance audit | rescue/audit case; could expose false-negative label or assay mismatch          |
|              11 | C_label_rescue        | CNV0_00711     | RREPPHLARNF | HLA-C*07:02         |       0 | LABEL_NOISE_RESCUE                     |                   0.579 |                          0.512 | false-negative rescue / assay mismatch | pMHC binding rescue after provenance audit | rescue/audit case; could expose false-negative label or assay mismatch          |
|              12 | C_label_rescue        | CNV0_00682     | NYRREPPHL   | HLA-C*07:02         |       0 | LABEL_NOISE_RESCUE                     |                   0.579 |                          0.518 | false-negative rescue / assay mismatch | pMHC binding rescue after provenance audit | rescue/audit case; could expose false-negative label or assay mismatch          |
|              13 | D_specificity_moat    | CNV0_02397     | LPIQYEPVL   | HLA-B*35:03         |       0 | HARD_NEGATIVE_OR_FALSE_NEGATIVE_AUDIT  |                   0.544 |                          0.506 | specificity moat                       | matched negative binding specificity       | hard-negative or false-negative audit; use to sharpen specificity               |
|              14 | D_specificity_moat    | CNV0_02400     | ETSKQVTRW   | HLA-A*25:01         |       0 | HARD_NEGATIVE_OR_FALSE_NEGATIVE_AUDIT  |                   0.542 |                          0.507 | specificity moat                       | matched negative binding specificity       | hard-negative or false-negative audit; use to sharpen specificity               |
|              15 | D_specificity_moat    | CNV0_02420     | KLANPLPYT   | HLA-A*02:01         |       0 | HARD_NEGATIVE_OR_FALSE_NEGATIVE_AUDIT  |                   0.537 |                          0.516 | specificity moat                       | matched negative binding specificity       | hard-negative or false-negative audit; use to sharpen specificity               |
|              16 | D_specificity_moat    | CNV0_02458     | IIGAGPAEV   | HLA-A*02:01         |       0 | HARD_NEGATIVE_OR_FALSE_NEGATIVE_AUDIT  |                   0.532 |                          0.522 | specificity moat                       | matched negative binding specificity       | hard-negative or false-negative audit; use to sharpen specificity               |
|              17 | E_positive_QC_control | CNV0_02473     | KEFEDDIINW  | HLA-B*44:03         |       1 | ASSAY_POSITIVE_CONTROL_OVERLAP_BLOCKED |                   0.522 |                          0.396 | assay QC positive control              | assay dynamic range/QC                     | positive assay control only; overlap/leakage blocks benchmark claim             |
|              18 | E_positive_QC_control | CNV0_00088     | YAYDNFGVLGL | HLA-C*03:03         |       1 | ASSAY_POSITIVE_CONTROL_OVERLAP_BLOCKED |                   0.521 |                          0.373 | assay QC positive control              | assay dynamic range/QC                     | positive assay control only; overlap/leakage blocks benchmark claim             |
|              19 | E_positive_QC_control | CNV0_00366     | QSPASLSSL   | HLA-C*01:02         |       1 | ASSAY_POSITIVE_CONTROL_OVERLAP_BLOCKED |                   0.521 |                          0.393 | assay QC positive control              | assay dynamic range/QC                     | positive assay control only; overlap/leakage blocks benchmark claim             |
|              20 | E_positive_QC_control | CNV0_00374     | NYPGFLTIL   | HLA-C*07:02         |       1 | ASSAY_POSITIVE_CONTROL_OVERLAP_BLOCKED |                   0.521 |                          0.393 | assay QC positive control              | assay dynamic range/QC                     | positive assay control only; overlap/leakage blocks benchmark claim             |
|              21 | F_model_boundary      | CNV0_01023     | FLDEFMEGV   | HLA-A*02:01         |       1 | TCR_STRUCTURE_DISAGREEMENT_AUDIT       |                   0.459 |                          0.363 | model-boundary clarification           | repeat MD/TCR review plus pMHC triage      | TCR/structure disagreement audit; no positive claim until resolved              |
|              22 | F_model_boundary      | CNV0_01251     | YLDELIRNT   | HLA-A*02:01         |       0 | TCR_STRUCTURE_DISAGREEMENT_AUDIT       |                   0.412 |                          0.334 | model-boundary clarification           | repeat MD/TCR review plus pMHC triage      | TCR/structure disagreement audit; no positive claim until resolved              |
|              23 | F_model_boundary      | CNV0_01099     | FRSGLDSYV   | HLA-C*06:02         |       1 | TCR_STRUCTURE_DISAGREEMENT_AUDIT       |                   0.410 |                          0.313 | model-boundary clarification           | repeat MD/TCR review plus pMHC triage      | TCR/structure disagreement audit; no positive claim until resolved              |
|              24 | F_model_boundary      | CNV0_01049     | LSSVSFFLV   | HLA-A*11:01         |       1 | TCR_STRUCTURE_DISAGREEMENT_AUDIT       |                   0.382 |                          0.292 | model-boundary clarification           | repeat MD/TCR review plus pMHC triage      | TCR/structure disagreement audit; no positive claim until resolved              |

## Claim unlock ladder

|   ladder_step | claim_layer                  | current_status   | evidence_now                                                     | unlock_rule                                                          | claim_after_unlock                                               | claim_boundary                                             |
|--------------:|:-----------------------------|:-----------------|:-----------------------------------------------------------------|:---------------------------------------------------------------------|:-----------------------------------------------------------------|:-----------------------------------------------------------|
|             1 | computational score recovery | achieved         | score booster source/HLA clean OOF mean AUPRC 0.623, AUROC 0.909 | already generated and verified                                       | validated computational prioritization layer                     | benchmark/prioritization only                              |
|             2 | clean antigen discovery      | plate-ready      | 6 clean BAR-Neo leads on plate v4                                | >=2 clean leads show mutant pMHC binding with WT/decoy specificity   | wetlab-supported generalizable antigen discovery lane            | no clinical vaccine-selection claim                        |
|             3 | TCR/MD mechanism             | plate-ready      | 2 TCR/MD structural leads on plate v4                            | >=1 structural lead passes pMHC binding plus TCR tetramer/activation | orthogonal mechanistic support for candidate prioritization      | not clean benchmark novelty because source overlap remains |
|             4 | label-noise rescue           | audit-ready      | 4 high-score negative-label rescue cases                         | >=1 rescue case passes provenance audit and pMHC binding             | label-noise/assay-mismatch evidence explaining benchmark ceiling | do not relabel without provenance review                   |
|             5 | specificity moat             | control-ready    | 4 hard negative/specificity cases plus positive controls         | hard negatives stay weak while positive controls pass                | specificity-aware assay and model-boundary defense               | assay control evidence, not patient outcome evidence       |
|             6 | clinical/patient relevance   | locked           | patient metadata gate remains limiting                           | complete patient context plus external/wetlab validation             | patient-level prioritization hypothesis                          | clinical vaccine-selection claim remains unavailable now   |

## Assay decision rules

| decision_node      | go_rule                                                             | no_go_rule                                                | action_if_go                                             | action_if_no_go                                        |
|:-------------------|:--------------------------------------------------------------------|:----------------------------------------------------------|:---------------------------------------------------------|:-------------------------------------------------------|
| assay_qc           | >=3/4 positive-control wells positive and blank/decoy controls weak | positive controls fail or widespread nonspecific binding  | interpret candidate arms                                 | repeat assay before candidate interpretation           |
| clean_BARNeo_leads | >=2/6 clean leads pass mutant binding and WT/decoy specificity      | 0/6 pass or all fail specificity                          | promote clean antigen discovery lane to wetlab-supported | downgrade clean booster score to computational-only    |
| TCR_MD_mechanism   | >=1/2 structural leads passes binding plus TCR readout              | binding or TCR readout fails in both structural leads     | promote orthogonal mechanism support                     | keep TCR/MD lane diagnostic-only                       |
| label_noise_rescue | >=1/4 rescue cases passes provenance audit and binding              | no rescue case passes after provenance audit              | use as label-noise explanation for benchmark ceiling     | remove rescue claim and keep as error-analysis reserve |
| specificity_moat   | hard negatives remain weak while positives pass                     | hard negatives bind strongly at similar rate as positives | strengthen specificity and control narrative             | flag model/assay specificity failure                   |

## Decision

1. The next impact layer is a 96-well experiment design, not another unrestricted model sprint.
2. The manuscript-safe frame is a controlled translational prioritization portfolio.
3. Clinical vaccine-selection remains locked until patient metadata and external/wetlab evidence are complete.

## Files

- `translational_value_portfolio_v4.tsv`
- `wetlab_plate_v4_value_of_information.tsv`
- `assay_96well_map_v4.tsv`
- `claim_unlock_ladder_v4.tsv`
- `assay_decision_rules_v4.tsv`
- `TRANSLATIONAL_IMPACT_V4_REPORT_KR.md`
