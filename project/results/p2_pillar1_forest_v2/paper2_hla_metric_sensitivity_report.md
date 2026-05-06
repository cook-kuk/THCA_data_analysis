# Paper 2 HLA metric sensitivity report

Primary result is exploratory because Korean PTC values are carrier counts over individuals, while baseline values are published allele frequencies over 2n.

## Source/provenance audit
| allele_4digit | korean_ptc_metric_type | korean_ptc_value | baseline_primary_source_id | baseline_primary_metric_type | baseline_primary_value | caveat |
| --- | --- | --- | --- | --- | --- | --- |
| A*02:07 | carrier_frequency | 0.0812 | IN2015 | allele_frequency | 3.51% | PTC value is carrier frequency from existing Korean PTC pool; baseline is published allele frequency. Do not convert or mix without Yu-approved analysis rule. |
| B*46:01 | carrier_frequency | 0.103 | IN2015 | allele_frequency | 5.06% | PTC value is carrier frequency from existing Korean PTC pool; baseline is published allele frequency. Do not convert or mix without Yu-approved analysis rule. |
| C*01:02 | carrier_frequency | 0.2414 | IN2015 | allele_frequency | 17.81% | PTC value is carrier frequency from existing Korean PTC pool; baseline is published allele frequency. C*01:02 PTC value comes from v1 PTC pool file; v1 GD comparator is deprecated, but this row uses only the Korean PTC carrier-frequency provenance. |
| DPB1*05:01 | carrier_frequency | 0.532 | AFND_SOUTH_KOREA_EXISTING_V2 | allele_frequency | 0.3667 | Current v2 already has AFND South Korea primary baseline. PTC value is carrier frequency and baseline is allele frequency; Jung 2023 is cross-check only until full table is archived/approved. |
| DQB1*02:01 | carrier_frequency | 0 | IN2015 | allele_frequency | 2.12% | PTC value is carrier frequency from existing Korean PTC pool and is zero in Combined row; baseline is published allele frequency. Do not compute any continuity correction or OR in this table. |
| DRB1*07:01 | carrier_frequency | 0.1144 | IN2015 | allele_frequency | 6.93% | PTC value is carrier frequency from existing Korean PTC pool; baseline is published allele frequency. Harbin Korean fallback should move to sensitivity because direct Korean DRB1*07:01 values are available. |

## Metric sensitivity
| allele_4digit | scenario | or | ci95_low | ci95_high | fisher_p | haldane_ci_used | direction | direction_vs_current |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A*02:07 | A_current_carrier_vs_baseline_allele | 2.433 | 1.648 | 3.59 | 5.584e-06 | False | enrichment | same |
| A*02:07 | B_baseline_AF_to_HWE_expected_carrier | 1.202 | 0.8087 | 1.787 | 0.3731 | False | enrichment | same |
| A*02:07 | C_PTC_carrier_to_HWE_implied_AF | 1.199 | 0.8166 | 1.76 | 0.3871 | False | enrichment | same |
| B*46:01 | A_current_carrier_vs_baseline_allele | 2.155 | 1.54 | 3.016 | 7.377e-06 | False | enrichment | same |
| B*46:01 | B_baseline_AF_to_HWE_expected_carrier | 1.058 | 0.7498 | 1.493 | 0.7933 | False | enrichment | same |
| B*46:01 | C_PTC_carrier_to_HWE_implied_AF | 1.043 | 0.7493 | 1.452 | 0.8666 | False | enrichment | same |
| C*01:02 | A_current_carrier_vs_baseline_allele | 1.472 | 1.189 | 1.821 | 0.0004346 | False | enrichment | same |
| C*01:02 | B_baseline_AF_to_HWE_expected_carrier | 0.6621 | 0.5264 | 0.8327 | 0.0004976 | False | depletion | flips |
| C*01:02 | C_PTC_carrier_to_HWE_implied_AF | 0.6866 | 0.5608 | 0.8406 | 0.0003045 | False | depletion | flips |
| DPB1*05:01 | A_current_carrier_vs_baseline_allele | 1.962 | 1.651 | 2.331 | 1.655e-14 | False | enrichment | same |
| DPB1*05:01 | B_baseline_AF_to_HWE_expected_carrier | 0.7626 | 0.6226 | 0.9341 | 0.009969 | False | depletion | flips |
| DPB1*05:01 | C_PTC_carrier_to_HWE_implied_AF | 0.7964 | 0.6858 | 0.9247 | 0.002882 | False | depletion | flips |
| DQB1*02:01 | A_current_carrier_vs_baseline_allele | 0.0259 | 0.001576 | 0.4256 | 8.217e-07 | True | depletion | same |
| DQB1*02:01 | B_baseline_AF_to_HWE_expected_carrier | 0.01268 | 0.0007709 | 0.2084 | 7.175e-11 | True | depletion | same |
| DQB1*02:01 | C_PTC_carrier_to_HWE_implied_AF | 0.01295 | 0.0007887 | 0.2128 | 8.425e-11 | True | depletion | same |
| DRB1*07:01 | A_current_carrier_vs_baseline_allele | 1.734 | 1.281 | 2.349 | 0.0004184 | False | enrichment | same |
| DRB1*07:01 | B_baseline_AF_to_HWE_expected_carrier | 0.8366 | 0.6123 | 1.143 | 0.2619 | False | depletion | flips |
| DRB1*07:01 | C_PTC_carrier_to_HWE_implied_AF | 0.8405 | 0.6245 | 1.131 | 0.252 | False | depletion | flips |

## DPB1 source sensitivity
DPB1*05:01 is flagged because its baseline source is AFND_SOUTH_KOREA_EXISTING_V2, not IN2015. Excluding DPB1 leaves five IN2015-backed alleles for sensitivity discussion.

## DQB1 zero-cell stress
DQB1*02:01 has 0/874 PTC carriers versus 26/1226 baseline alleles. Exact Fisher p=8.22e-07; Haldane OR for CI=0.0259. This is a depletion signal with zero-cell caution.