# Paper 2 HLA subcohort heterogeneity report

Subcohorts available locally: K2, Lee2024, and GSE286332_PTC. GSE286332_PTC has n=9, so heterogeneity flags should not be over-interpreted.

| allele_4digit | n_subcohorts | pooled_carriers | pooled_n | pooled_frequency | min_frequency | min_frequency_cohort | max_frequency | max_frequency_cohort | dominant_carrier_share | dominant_cohort | chi_square_p | flags |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A*02:07 | 3 | 71 | 874 | 0.08124 | 0.0746 | Lee2024 | 0.3333 | GSE286332_PTC | 0.662 | Lee2024 | 0.01622 | small_subcohort_present;heterogeneity_p_lt_0_05 |
| B*46:01 | 3 | 90 | 874 | 0.103 | 0.1 | Lee2024 | 0.3333 | GSE286332_PTC | 0.7 | Lee2024 | 0.07308 | small_subcohort_present |
| C*01:02 | 3 | 211 | 874 | 0.2414 | 0.183 | K2 | 0.5556 | GSE286332_PTC | 0.7725 | Lee2024 | 0.005908 | small_subcohort_present;carrier_count_dominated_by_one_subcohort;heterogeneity_p_lt_0_05 |
| DPB1*05:01 | 3 | 465 | 874 | 0.532 | 0.5206 | Lee2024 | 0.5617 | K2 | 0.7054 | Lee2024 | 0.5544 | small_subcohort_present |
| DQB1*02:01 | 3 | 0 | 874 | 0 | 0 | K2 | 0 | K2 |  |  |  | small_subcohort_present |
| DRB1*07:01 | 3 | 100 | 874 | 0.1144 | 0.0979 | K2 | 0.3333 | GSE286332_PTC | 0.74 | Lee2024 | 0.0842 | small_subcohort_present |

The six-allele PTC pool is not dominated by one large subcohort for most carriers, but the n=9 GSE286332_PTC row produces unstable high frequencies for several alleles. Use subcohort plots as transparency, not validation.