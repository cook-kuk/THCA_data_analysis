# Paper 2 HLA Allele-vs-Baek 2021 NGS Baseline Validation

Date: 2026-05-06

## Executive Verdict

The carrier-vs-allele mismatch was addressed by recounting Korean PTC HLA calls as allele copies from the two allele columns per locus, then comparing those counts with Baek 2021 healthy Korean NGS allele counts. This creates a metric-matched allele-frequency validation table.

This does not make the result final clinical or population-genetic proof. It is still a sensitivity/validation layer because the PTC calls and Baek 2021 controls differ by cohort design, platform, callable-locus completeness, and control matching.

## Primary Results

| Allele | PTC allele count | Baek 2021 NGS control count | OR | 95% CI | Fisher p | Direction |
|---|---:|---:|---:|---:|---:|---|
| A*02:07 | 77/1748 (4.41%) | 18/346 (5.20%) | 0.84 | 0.496-1.42 | 0.4821 | no_clear_difference |
| B*46:01 | 93/1748 (5.32%) | 21/346 (6.07%) | 0.87 | 0.534-1.42 | 0.6035 | no_clear_difference |
| C*01:02 | 225/1746 (12.89%) | 72/346 (20.81%) | 0.563 | 0.419-0.756 | 2.00e-04 | depleted_in_ptc |
| DPB1*05:01 | 564/1558 (36.20%) | 118/346 (34.10%) | 1.1 | 0.858-1.4 | 0.4955 | no_clear_difference |
| DQB1*02:01 | 0/1262 (0.00%) | 8/346 (2.31%) | 0.0158 | 0.000908-0.274 | 4.31e-06 | depleted_in_ptc |
| DRB1*07:01 | 102/1722 (5.92%) | 27/346 (7.80%) | 0.744 | 0.479-1.16 | 0.1819 | no_clear_difference |

## Files Created

- `project/results/p2_pillar1_forest_v2/paper2_ptc_allele_dosage_counts.tsv`
- `project/results/p2_pillar1_forest_v2/paper2_allele_vs_baek2021_ngs_validation.tsv`
- `project/papers_hub_2026_05_04/assets/paper2_hla/F19_allele_vs_baek2021_ngs_forest.png`
- `project/papers_hub_2026_05_04/assets/paper2_hla/F20_allele_frequency_inputs_baek2021.png`
- `project/papers_hub_2026_05_04/assets/paper2_hla/F21_metric_fix_flow_allele_vs_allele.png`

## Method Boundary

- PTC denominator is the number of non-missing allele calls per locus.
- Baek 2021 denominator is 346 alleles from 173 healthy South Koreans.
- No carrier-to-allele conversion was performed.
- No allele-to-carrier conversion was performed.
- Fisher exact tests were run on allele-count tables.
- Haldane correction was used only for CI/OR display when a zero cell occurred.

## Remaining Caveats

- PTC HLA calls are derived from the Korean PTC pool, not from a matched case-control genotyping study.
- Baek 2021 is an open Korean healthy NGS baseline, but it is not matched by age, sex, region, platform, or technical pipeline to the PTC calls.
- DQB1 and DPB1 have lower PTC callable denominators than A/B/C.
- Final claim still needs matched Korean control carrier or allele-count validation and HLA QC review.
