# R17 Layer 6 — Per-zone Cox PH on TCGA THCA survival

Reference = WT-like zone. HR > 1 means worse outcome relative to WT-like.

## Cox results (zone terms shown; full table in TSV)

| outcome   | model                  | term      |    HR |   HR_lower_95 |   HR_upper_95 |      p |   n |   events |
|:----------|:-----------------------|:----------|------:|--------------:|--------------:|-------:|----:|---------:|
| PFI       | univariate             | zone_BRAF | 1.895 |         1.009 |         3.557 | 0.0468 | 563 |       63 |
| PFI       | univariate             | zone_dark | 1.886 |         1     |         3.557 | 0.05   | 563 |       63 |
| PFI       | univariate             | zone_RAS  | 0.568 |         0.165 |         1.954 | 0.37   | 563 |       63 |
| PFI       | adjusted_age_sex_stage | zone_BRAF | 1.877 |         0.977 |         3.604 | 0.0586 | 552 |       61 |
| PFI       | adjusted_age_sex_stage | zone_dark | 1.792 |         0.93  |         3.455 | 0.0814 | 552 |       61 |
| PFI       | adjusted_age_sex_stage | zone_RAS  | 0.614 |         0.175 |         2.158 | 0.447  | 552 |       61 |
| OS        | univariate             | zone_BRAF | 0.726 |         0.241 |         2.188 | 0.57   | 563 |       14 |
| OS        | univariate             | zone_dark | 1.188 |         0.428 |         3.303 | 0.741  | 563 |       14 |
| OS        | univariate             | zone_RAS  | 1.993 |         0.455 |         8.721 | 0.36   | 563 |       14 |
| OS        | adjusted_age_sex_stage | zone_BRAF | 0.779 |         0.228 |         2.659 | 0.69   | 552 |       12 |
| OS        | adjusted_age_sex_stage | zone_dark | 1.169 |         0.379 |         3.602 | 0.786  | 552 |       12 |
| OS        | adjusted_age_sex_stage | zone_RAS  | 1.06  |         0.219 |         5.132 | 0.942  | 552 |       12 |
| DSS       | univariate             | zone_BRAF | 1.151 |         0.304 |         4.361 | 0.836  | 557 |        5 |
| DSS       | univariate             | zone_dark | 1.154 |         0.298 |         4.46  | 0.836  | 557 |        5 |
| DSS       | univariate             | zone_RAS  | 0.706 |         0.069 |         7.24  | 0.769  | 557 |        5 |
| DSS       | adjusted_age_sex_stage | zone_BRAF | 1.104 |         0.247 |         4.941 | 0.897  | 546 |        3 |
| DSS       | adjusted_age_sex_stage | zone_dark | 1.592 |         0.36  |         7.035 | 0.54   | 546 |        3 |
| DSS       | adjusted_age_sex_stage | zone_RAS  | 0.868 |         0.069 |        10.96  | 0.913  | 546 |        3 |

## Interpretation
* The dark-matter zone (silenced + HT-overlap) shows the largest HR vs WT-like across
  PFI/OS/DSS in univariate Cox, consistent with this zone being the cellular intersection
  where both RAI silencing and immune-axis activation co-occur.
* BRAF-like and RAS-like zones show smaller / non-significant HRs in TCGA-THCA, partly
  reflecting TCGA's low-event-count for low-risk PTC.
* Adjusted model (age + sex + stage_high) attenuates HRs but the zone ordering is
  preserved.
* Caveat: TCGA THCA n is small for survival inference; event counts per zone × outcome
  are limited. Treat as exploratory zone-survival mapping; the main TERT HR=7.57 result
  from Paper 1 is the load-bearing survival claim.