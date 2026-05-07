# Paper 2 Kim 2014 Korean HLA reference-panel validation

## Executive verdict
- Kim et al. Korean HLA Reference Panel v1.0 is usable as a Korean control-like carrier-frequency baseline candidate for Paper 2.
- The panel contains 413 unrelated Korean subjects and covers all six target loci: HLA-A, B, C, DRB1, DPB1, DQB1.
- Unlike the previous In/AFND baseline comparison, this validation is metric-matched: PTC carrier frequency vs Kim reference-panel carrier frequency.
- Caveat: the panel is an imputation reference panel, not a phenotype-matched healthy-control cohort. Use as validation/sensitivity unless Yu approves primary use.

## Result table

| allele_4digit | ptc_carriers | ptc_n_individuals | kim_carriers | kim_n_individuals | kim_carrier_frequency | kim_allele_frequency_from_bgl | or_haldane_for_ci | ci95_low | ci95_high | fisher_exact_p | interpretation_status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A*02:07 | 71 | 874 | 26 | 413 | 0.06295 | 0.03269 | 1.316 | 0.8264 | 2.096 | 0.2603 | metric_matched_reference_panel_validation_candidate |
| B*46:01 | 90 | 874 | 45 | 413 | 0.109 | 0.05811 | 0.9388 | 0.6428 | 1.371 | 0.7702 | metric_matched_reference_panel_validation_candidate |
| C*01:02 | 211 | 874 | 127 | 413 | 0.3075 | 0.1646 | 0.7167 | 0.5525 | 0.9296 | 0.0145 | metric_matched_reference_panel_validation_candidate |
| DPB1*05:01 | 465 | 874 | 251 | 413 | 0.6077 | 0.3644 | 0.7338 | 0.5783 | 0.931 | 0.01159 | metric_matched_reference_panel_validation_candidate |
| DQB1*02:01 | 0 | 874 | 20 | 413 | 0.04843 | 0.02421 | 0.01097 | 0.0006621 | 0.1819 | 9.747e-11 | metric_matched_reference_panel_validation_candidate |
| DRB1*07:01 | 100 | 874 | 58 | 413 | 0.1404 | 0.07143 | 0.7908 | 0.5589 | 1.119 | 0.2027 | metric_matched_reference_panel_validation_candidate |

## Files
- `project/results/p2_pillar1_forest_v2/paper2_kim2014_reference_panel_carrier_counts.tsv`
- `project/results/p2_pillar1_forest_v2/paper2_kim2014_reference_panel_validation.tsv`
- `project/papers_hub_2026_05_04/assets/paper2_hla/F17_kim2014_reference_panel_carrier_validation_forest.png`
- `project/papers_hub_2026_05_04/assets/paper2_hla/F18_kim2014_metric_matched_carrier_inputs.png`
