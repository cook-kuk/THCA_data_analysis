---
title: "Paper 2 HLA 6-allele exploratory statistics"
date: 2026-05-06
status: "COMPLETE - Paper 2 exploratory OR/Fisher/forest only"
---

# Executive verdict

The 6-allele Paper 2 HLA exploratory OR/Fisher/forest has been executed after approval.
The analysis keeps the metric mismatch explicit: Korean PTC uses carrier counts over individuals, while Korean baseline uses published allele frequencies reconstructed as allele counts over 2n.
No carrier-to-allele or allele-to-carrier conversion was performed.
No Paper 4 GD meta-analysis was executed.

# Primary exploratory table

| Allele | PTC carriers/n | Baseline alleles/2n | OR, Haldane CI | Fisher exact p | Baseline source |
|---|---:|---:|---:|---:|---|
| A*02:07 | 71/874 | 43/1226 | 2.43 (1.65-3.59) | 5.58e-06 | IN2015 |
| B*46:01 | 90/874 | 62/1226 | 2.16 (1.54-3.02) | 7.38e-06 | IN2015 |
| C*01:02 | 211/874 | 218/1226 | 1.47 (1.19-1.82) | 0.000435 | IN2015 |
| DPB1*05:01 | 465/874 | 499/1360 | 1.96 (1.65-2.33) | 1.65e-14 | AFND_SOUTH_KOREA_EXISTING_V2 |
| DQB1*02:01 | 0/874 | 26/1226 | 0.0259 (0.00158-0.426) | 8.22e-07 | IN2015 |
| DRB1*07:01 | 100/874 | 85/1226 | 1.73 (1.28-2.35) | 0.000418 | IN2015 |

# Files created

- `project/results/p2_pillar1_forest_v2/paper2_6allele_or_fisher_primary.tsv`
- `project/results/p2_pillar1_forest_v2/paper2_6allele_or_fisher_primary.json`
- `project/results/p2_pillar1_forest_v2/paper2_6allele_exploratory_or_fisher_forest.png`
- `project/results/p2_pillar1_forest_v2/paper2_6allele_exploratory_or_fisher_forest.pdf`
- `project/papers_hub_2026_05_04/assets/paper2_hla/F13_exploratory_or_fisher_forest.png`
- `project/papers_hub_2026_05_04/assets/paper2_hla/F14_exploratory_stats_table.png`
- `project/papers_hub_2026_05_04/assets/paper2_hla/F15_exploratory_metric_inputs.png`

# Caveat

These are exploratory metric-mismatched statistics. They are useful as a decision/visualization layer, not as a clean population-genetic case-control estimate unless the metric rule is accepted in writing.
