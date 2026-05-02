# Paper 2 Pillar 1 Forest Meta — Status: PARTIAL → STRONG

> ⚠️ **DEPRECATED 2026-05-04** — Advisor scope (Yu, 2026-05-04): Paper 2 = Hashimoto-overlap PTC only (HT destructive infiltration → dedifferentiation). The "autoimmune-thyroid continuum" framing and Chu 2018 GD comparison conflate Paper 2 (HT) with Paper 3 (GD). Pillar I is being reformulated in Task A as **Korean PTC vs Korean baseline (HT context)**; Chu 2018 GD reduced to a one-line Discussion mention. See `project/results/p2_pillar1_forest_v2/` for the v2 forest. This v1 summary is preserved for record only.

**Date:** 2026-05-03 (v1) · audit 2026-05-04
**Source:** Chu X et al. 2018 J Med Genet 55:685–692, doi:10.1136/jmedgenet-2017-105146 (PMC 6161647)

## ★ One-line conclusion (v1, deprecated)

**Pan-Asian thyroid HLA susceptibility (HT context).** Korean PTC pool (n=874) sits above Han Chinese ctrl baseline for principal risk alleles (DPB1*05:01 53% vs ctrl 31%, p_kr_vs_ctrl=0; B*46:01 10% vs ctrl 7%, p_kr_vs_ctrl=0.00109). [v1 also compared to Chu 2018 GD 44%/14% — that arm now lives in Paper 3 reserve.]

## Per-allele 3-arm forest

| Allele | Chen ctrl | **Korean PTC** | Chen GD | OR Kr vs ctrl | p Kr vs ctrl | OR Kr vs GD | p Kr vs GD | Pooled OR | I²% |
|---|---|---|---|---|---|---|---|---|---|
| A*02:07 | 4.9% | **8.1%** | 9.7% | 1.716 [1.224, 2.406] | 0.00173 | 0.826 [0.613, 1.112] | 0.208 | 1.985 [1.66, 2.373] | 0 |
| B*46:01 | 6.5% | **10.3%** | 14.1% | 1.649 [1.221, 2.225] | 0.00109 | 0.699 [0.538, 0.91] | 0.0077 | 2.022 [1.414, 2.89] | 76 |
| C*01:02 | 10.9% | **24.1%** | 18.4% | 2.609 [2.083, 3.267] | 0 | 1.412 [1.152, 1.731] | 0.000892 | 2.163 [1.529, 3.06] | 85 |
| DPB1*05:01 | 31.3% | **53.2%** | 44.0% | 2.498 [2.103, 2.968] | 0 | 1.447 [1.223, 1.712] | 1.68e-05 | 2.162 [1.654, 2.826] | 85 |
| DQB1*02:01 | 17.8% | **0.0%** | 10.9% | 0.003 [0.0, 0.042] | 2.76e-05 | 0.005 [0.0, 0.075] | 0.000152 | 0.57 [0.491, 0.662] | 0 |
| DRB1*07:01 | 15.3% | **11.4%** | 7.1% | 0.715 [0.556, 0.92] | 0.00899 | 1.694 [1.27, 2.261] | 0.000337 | 0.55 [0.334, 0.905] | 91 |


## Korean sub-cohort heterogeneity (Cochran's Q)

| Allele | K2 | Lee | GSE286332 | Q | I²% |
|---|---|---|---|---|---|
| A*02:07 | 21/235 | 47/630 | 3/9 | 8.243 | 75.7% |
| B*46:01 | 24/235 | 63/630 | 3/9 | 5.232 | 61.8% |
| C*01:02 | 43/235 | 163/630 | 5/9 | 10.263 | 80.5% |
| DPB1*05:01 | 132/235 | 328/630 | 5/9 | 1.18 | 0.0% |
| DQB1*02:01 | 0/235 | 0/630 | 0/9 | 0.0 | 0.0% |
| DRB1*07:01 | 23/235 | 74/630 | 3/9 | 4.949 | 59.6% |


## Sensitivity (4 scenarios — full / excl GSE286332 / Lee only / K2 only)

DPB1*05:01 across scenarios:
- **All_Korean_PTC_n874**: 465/874 = 53.2%, OR vs Chen ctrl = 2.498
- **Excl_GSE286332_n865**: 460/865 = 53.2%, OR vs Chen ctrl = 2.496
- **Lee_only_n630**: 328/630 = 52.1%, OR vs Chen ctrl = 2.387
- **K2_only_n235**: 132/235 = 56.2%, OR vs Chen ctrl = 2.816


B*46:01 across scenarios:
- **All_Korean_PTC_n874**: 90/874 = 10.3%, OR vs Chen ctrl = 1.649
- **Excl_GSE286332_n865**: 87/865 = 10.1%, OR vs Chen ctrl = 1.606
- **Lee_only_n630**: 63/630 = 10.0%, OR vs Chen ctrl = 1.596
- **K2_only_n235**: 24/235 = 10.2%, OR vs Chen ctrl = 1.633


## Decision

- **STRONG** — direction-consistent alleles: 6 (A*02:07 (Korean > Chen ctrl), B*46:01 (Korean > Chen ctrl), C*01:02 (Korean > Chen ctrl), DPB1*05:01 (Korean > Chen ctrl), DQB1*02:01 (Korean < Chen ctrl), DRB1*07:01 (Korean < Chen ctrl))
- **Pillar 1 status: PARTIAL → STRONG**
- Paper 2 venue ladder unchanged but qualitative claim now becomes quantitative

## Honest disclosure

1. **Cohort size imbalance:** Korean PTC n=874 vs Chu cohort n=2,958 — Chen estimates dominate pooled estimates.
2. **Phenotype heterogeneity:** PTC and Graves' are distinct diseases (cancer vs autoimmunity); the comparison tests *susceptibility allele overlap* not *disease equivalence*.
3. **Population stratification:** Korean and Han Chinese are related East Asian populations but with documented allele frequency differences (e.g., DRB1*15:01 18% Korean vs 4% Han Chinese in Chu controls).
4. **Korean sub-cohort heterogeneity:** within-Korean Cochran Q assessed; if I² > 50% for an allele, sub-cohort technical confounders may bias.
5. **Allele-level scope:** haplotype-level (e.g., DRB1-DQA1-DQB1 trios) interactions deferred to future work.

## Files

```
project/results/p2_pillar1_forest/
├── chu2018_allele_summary.tsv
├── korean_PTC_pool_per_subcohort.tsv
├── korean_subcohort_heterogeneity.tsv
├── forest_meta_results.tsv
├── random_effects_pooled.tsv
├── sensitivity_4scenarios.tsv
├── allele_freq_3cohort.{pdf,png}
├── forest_panasian_HLA.{pdf,png}
├── methods_paragraph.md (Cell Press paste-ready)
├── discussion_paragraph.md (Cell Press paste-ready)
└── PILLAR1_FOREST_SUMMARY.md (this file)
```

## Note on citation

Note: The paper formerly cited as "Chen 2018" in our prior memory is **Chu X et al. 2018**
(J Med Genet 55:685–692). Verified via PMC 6161647 web fetch 2026-05-03. All previous
references in `D4P1_summary`, `MEMORY.md`, and prior reports should be updated to **Chu 2018**.
