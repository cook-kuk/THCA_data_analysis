# HLA validation and screening meta update

## Executive verdict
- Paper 2 matched-control validation is the correct next step, but a final population case-control claim still requires independent Korean healthy/control carrier-frequency or individual-level HLA genotype data.
- Local GSE213647 adjacent-normal arcasHLA calls can support a technical sanity check only; they are cancer-patient normal tissues, not independent healthy controls.
- Paper 4 screening Pan-Asian GD meta-analysis was executed only for extractable Chu/Shin/Chen source rows. It is not the final full-source meta-analysis.

## Files created
- `project/results/p2_pillar1_forest_v2/paper2_matched_control_validation_gate.tsv`
- `project/papers_hub_2026_05_04/assets/paper2_hla/F16_paper2_matched_control_validation_gate.png`
- `project/results/paper4_gd_hla/paper4_screening_panasian_meta_source_rows.tsv`
- `project/results/paper4_gd_hla/paper4_screening_panasian_meta.tsv`
- `project/results/paper4_gd_hla/paper4_meta_readiness_by_allele.tsv`
- `project/papers_hub_2026_05_04/assets/paper4_hla/P4_F15_screening_panasian_meta_forest.png`
- `project/papers_hub_2026_05_04/assets/paper4_hla/P4_F16_meta_readiness_by_allele.png`

## Paper 4 screening meta results

| allele | k | pooled_or | ci_low | ci_high | p_random | tau2 | q | i2_pct | sources | meta_status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A*02:07 | 1 | 2.1 | 1.701 | 2.592 | 2.07e-12 | nan | nan | nan | Chu 2018 | single_source_anchor_only |
| B*46:01 | 3 | 2.038 | 1.208 | 3.437 | 0.007623 | 0.1599 | 17.86 | 88.8 | Chu 2018; Shin 2019; Chen 2011 | screening_random_effects |
| C*01:02 | 2 | 1.846 | 1.592 | 2.141 | 4.639e-16 | 0 | 0.4831 | 0 | Chu 2018; Shin 2019 | screening_random_effects |
| DPB1*05:01 | 3 | 2.168 | 1.675 | 2.806 | 4.144e-09 | 0.02736 | 4.696 | 57.41 | Chu 2018; Shin 2019; Chen 2011 | screening_random_effects |
| DQB1*02:01 | 1 | 0.57 | 0.4911 | 0.6615 | 2.31e-13 | nan | nan | nan | Chu 2018 | single_source_anchor_only |
| DRB1*07:01 | 1 | 0.43 | 0.3613 | 0.5118 | 2.49e-21 | nan | nan | nan | Chu 2018 | single_source_anchor_only |

## Caveats
- Shin 2019 and Chen 2011 rows use corrected p-values to approximate SE where exact count/CI rows are not locally extracted.
- A*02:07, DQB1*02:01, and DRB1*07:01 remain Chu-only anchors in this update.
- Liao 2022, Park 2005, Ueda 2014, and older Japanese/Hong Kong/Taiwan/Korean studies still need source table extraction before final Paper 4 meta-analysis.
