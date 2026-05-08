# Track 8 — AFND Pan-Asian + South + Southeast Asian Healthy HLA Atlas

_Generated: 2026-05-08; Track 8 of HLA deep-dive sprint._

## 0 Boundary

Healthy populations only. **No** disease association is computed or implied. Boundary contract: `project/paper2_hla_boundary/HLA_CANCER_SEPARATION_RULES.md`. All allele frequencies below are AFND population baselines; cancer cohorts (TCGA, K2, Lee, GSE) are not included and **must not** be joined to these baselines except under Paper 2/4 explicit cohort-vs-AFND contrast pre-registered separately.

## 1 Data extraction methods

AFND `hla6006a.asp` country-locus query for 92 (country, locus, page) combinations. Budget 92 of 100 GETs; 8 remaining. 1.5s sleep between requests. Cache: `cache/{Country}__{LOCUS}__pN.html`. Parser extracts `tblNormal` rows → (allele, population, allele_freq, sample_size); cells [1]/[3]/[5]/[7] (header has 10 cells but body rows have 12 due to colspan). Allele names truncated to 2-field resolution (`A*02:01:01:01` → `A*02:01`). Within (population, allele_2f) duplicates are collapsed to the largest-sample-size record.

**Fetch failures / workarounds**: AFND throttled near the budget limit (slow response on the Indonesia/Malaysia DPB1 query around request ~92). The fetch was interrupted at 92/100 GETs and the manifest was reconstructed from the cache. Outgroup populations Germany / United Kingdom / USA were planned but not reached before the throttle. The atlas is therefore Pan-Asian-focused (16 Asian countries) with no European outgroup overlay. Adding outgroup later requires only ~18 more GETs (3 countries × 6 loci) and the cache+parse infrastructure handles re-runs idempotently.

Country-level fetch summary:

| country     | region         |   n_pops |   n_loci |   n_alleles |   total_records |
|:------------|:---------------|---------:|---------:|------------:|----------------:|
| China       | East Asia      |       78 |        6 |          58 |             537 |
| India       | South Asia     |       33 |        6 |         109 |             485 |
| Taiwan      | East Asia      |       33 |        6 |          45 |             513 |
| Japan       | East Asia      |       23 |        6 |         117 |             511 |
| Malaysia    | Southeast Asia |       18 |        6 |         164 |             419 |
| Mongolia    | East Asia      |       11 |        6 |         177 |             392 |
| Thailand    | Southeast Asia |       11 |        6 |         133 |             378 |
| Indonesia   | Southeast Asia |        8 |        5 |         118 |             355 |
| Pakistan    | South Asia     |        8 |        6 |         156 |             365 |
| Singapore   | Southeast Asia |        8 |        6 |         166 |             496 |
| South Korea | East Asia      |        8 |        6 |         283 |             519 |
| Vietnam     | Southeast Asia |        5 |        5 |         157 |             342 |
| Philippines | Southeast Asia |        2 |        4 |         280 |             320 |
| Sri Lanka   | South Asia     |        2 |        6 |         207 |             214 |

## 2 Master matrix overview

**245 populations × 6 loci × 695 unique 2-field alleles** across 14 countries. Per-locus wide pivots saved as `tables/master_matrix__{LOCUS}.tsv` (rows = populations, columns = top 100 alleles by sample-size-weighted pooled frequency). Long form with Clopper-Pearson 95% CIs at `tables/master_long_with_ci.tsv`.

By region:

| region         |   n_pops |   n_alleles |
|:---------------|---------:|------------:|
| East Asia      |      150 |         385 |
| South Asia     |       43 |         332 |
| Southeast Asia |       52 |         474 |

## 3 PCA + dendrogram

PCA on a population × allele-frequency vector built by **per-locus normalization** (each locus block sums to 1 per population). This avoids the artifact where populations with sparse locus coverage in AFND look anomalous because their missing loci are zero-filled. Populations are kept only if they have >=4 of 6 loci with at least one allele record. Dendrogram uses Jensen-Shannon divergence with average linkage. Outputs: `figures/F01_pca.{png,pdf}`, `figures/F02_umap.png`, `figures/F03_dendrogram.png`.

**Korean populations sit tightly with Japanese populations** (both at PC1 ≈ -0.4 to -0.7, PC2 near 0), distinct from Han Chinese populations (PC1 +0.6 to +0.9) and from Taiwan aboriginal populations (PC1 +0.5, PC2 +0.7). South Asian populations form a separate cluster on PC2 (negative PC2). This is consistent with the published East-Asia HLA literature (Korea-Japan ancestry overlap, Han structure, Taiwan aboriginal divergence).

Top 10 alleles distinguishing East Asia vs South Asia (East mean − South mean):

| allele_2f        |   east_minus_south |    abs |
|:-----------------|-------------------:|-------:|
| A::A*02:01       |             0.4606 | 0.4606 |
| A::A*01:01       |            -0.4282 | 0.4282 |
| C::C*01:02       |             0.3737 | 0.3737 |
| DRB1::DRB1*10:01 |            -0.3367 | 0.3367 |
| DQB1::DQB1*03:01 |            -0.2730 | 0.2730 |
| B::B*14:01       |            -0.1923 | 0.1923 |
| DQB1::DQB1*02:02 |            -0.1911 | 0.1911 |
| C::C*06:02       |            -0.1752 | 0.1752 |
| DRB1::DRB1*03:01 |            -0.1414 | 0.1414 |
| C::C*03:02       |            -0.1300 | 0.1300 |

## 4 Pairwise Reynolds distance (FST analog)

Per-locus Reynolds distance (Reynolds, Weir & Cockerham 1983) computed from allele frequencies, averaged across loci. Heatmap reordered by hierarchical clustering: `figures/F04_fst_heatmap.png`. Full distance matrix: `tables/fst_matrix.tsv`.

## 5 Per-locus diversity

Shannon entropy and expected heterozygosity (1−Σp²) per population × locus: `tables/locus_diversity.tsv`. Boxplot by region × locus: `figures/F05_diversity.png`.

Median diversity per region × locus:

| region         | locus   |   median_H |   median_het |   n_pops |
|:---------------|:--------|-----------:|-------------:|---------:|
| East Asia      | A       |     0.1956 |       0.0933 |       58 |
| East Asia      | B       |     0.4675 |       0.4421 |       53 |
| East Asia      | C       |     0.1113 |       0.0474 |       57 |
| East Asia      | DPB1    |     1.0808 |       0.6247 |       55 |
| East Asia      | DQB1    |     1.0438 |       0.6252 |       52 |
| East Asia      | DRB1    |     0.6484 |       0.5827 |      103 |
| South Asia     | A       |     0.8218 |       0.4880 |       29 |
| South Asia     | B       |     0.7486 |       0.5310 |       28 |
| South Asia     | C       |     0.6365 |       0.4420 |       26 |
| South Asia     | DPB1    |     2.2175 |       0.8079 |        5 |
| South Asia     | DQB1    |     0.8810 |       0.4991 |       34 |
| South Asia     | DRB1    |     0.6392 |       0.4032 |       35 |
| Southeast Asia | A       |     1.5353 |       0.7189 |       30 |
| Southeast Asia | B       |     1.5089 |       0.7142 |       31 |
| Southeast Asia | C       |     2.3057 |       0.8778 |       15 |
| Southeast Asia | DPB1    |     2.0490 |       0.8330 |       13 |
| Southeast Asia | DQB1    |     2.0157 |       0.8198 |       26 |
| Southeast Asia | DRB1    |     1.6293 |       0.7225 |       43 |

## 6 Korean baseline triangulation

Top-8 alleles per locus (by Korean weighted-mean frequency) compared across South Korea / Japan / Han-Beijing / East-Asian pool. Forest grid: `figures/F07_korean_triangulation_forest.png`. Table: `tables/korean_triangulation.tsv`.

Top-3 per locus (Korea weighted mean):

| locus   | allele     |   Korea_freq |   Korea_n |   Japan_freq |   Japan_n |   Han_Beijing_freq |   Han_Beijing_n |   EastAsia_pool_freq |
|:--------|:-----------|-------------:|----------:|-------------:|----------:|-------------------:|----------------:|---------------------:|
| A       | A*02:01    |       0.0956 |      7523 |       0.1155 |     20446 |             0.1395 |             790 |               0.1086 |
| A       | A*11:01    |       0.0852 |      5583 |       0.0912 |     20303 |           nan      |               0 |               0.0901 |
| A       | A*02:06    |       0.0850 |      4613 |       0.0912 |     20303 |           nan      |               0 |               0.0895 |
| B       | B*15:01    |       0.0797 |     14134 |       0.0765 |     20160 |           nan      |               0 |               0.0777 |
| B       | B*07:02    |       0.0316 |      6068 |       0.0564 |     20160 |             0.0198 |            1210 |               0.0351 |
| B       | B*13:02    |       0.0313 |      4613 |       0.0029 |     20110 |           nan      |               0 |               0.0083 |
| C       | C*01:02    |       0.1830 |       485 |       0.1761 |     20885 |             0.0916 |             844 |               0.1673 |
| C       | C*03:02    |       0.0908 |       809 |       0.0057 |     20553 |           nan      |               0 |               0.0092 |
| C       | C*08:01    |       0.0840 |       809 |       0.0746 |     20835 |           nan      |               0 |               0.0750 |
| DPB1    | DPB1*05:01 |       0.3667 |       680 |       0.3584 |      5632 |           nan      |               0 |               0.3747 |
| DPB1    | DPB1*02:01 |       0.2500 |       680 |       0.2399 |      5632 |             0.2310 |             180 |               0.1980 |
| DPB1    | DPB1*04:01 |       0.0873 |       680 |       0.0515 |      5632 |           nan      |               0 |               0.0665 |
| DQB1    | DQB1*03:02 |       0.1001 |      1227 |       0.1090 |      4723 |             0.0610 |             171 |               0.0776 |
| DQB1    | DQB1*04:01 |       0.0817 |      1227 |       0.1377 |      4723 |           nan      |               0 |               0.1072 |
| DQB1    | DQB1*03:01 |       0.0771 |      2197 |       0.1150 |      4723 |             0.2050 |             171 |               0.1706 |
| DRB1    | DRB1*04:05 |       0.0712 |      7274 |       0.1358 |     25498 |           nan      |               0 |               0.1182 |
| DRB1    | DRB1*01:01 |       0.0651 |      6304 |       0.0590 |     25448 |             0.0182 |             789 |               0.0438 |
| DRB1    | DRB1*04:06 |       0.0507 |      6304 |       0.0291 |      1053 |           nan      |               0 |               0.0418 |

## 7 Track 1 — 8-allele baseline check

Korean baseline frequencies for the 8 Track 1 focus alleles. Source: AFND v3 healthy extraction via Track 8 (`scripts/hla_deepdive_2026_05_08/track8/01_fetch_afnd.py`), with cross-fill from the Track 3 hand-curated focus-allele dump (`project/manuscript_p2_brief/lit_enrich_2026_05_02/data/afnd_alleles.json`) for alleles that were not on AFND page-1 of the Korean query. Forest plot: `figures/F08_track1_baseline.png`. Table: `tables/track1_korean_baseline.tsv`.

**Note**: 2 of 8 Track 1 alleles have no Korean record in either source: `DRB1*08:02`, `DRB1*16:02`. These are sub-1% Korean alleles that did not appear on AFND page-1 for any Korean study and were not in the Track 3 hand-curated dump. Re-running the fetcher with the freed 8-GET budget on Korea-specific page-2/3 queries would resolve this gap.

Korean baseline (sample-size weighted; Clopper-Pearson 95% CI on pooled chromosomes):

| allele     |   n_pops |   n_pooled_individuals |   weighted_mean_freq |   ci_lo |   ci_hi |
|:-----------|---------:|-----------------------:|---------------------:|--------:|--------:|
| DPB1*05:01 |        3 |                    680 |               0.3667 |  0.3412 |  0.3932 |
| A*02:07    |        2 |                   4613 |               0.0336 |  0.0300 |  0.0375 |
| C*03:02    |        2 |                    809 |               0.0908 |  0.0773 |  0.1059 |
| DQB1*03:02 |        4 |                   1227 |               0.1001 |  0.0886 |  0.1128 |
| B*46:01    |        3 |                  82197 |               0.0505 |  0.0494 |  0.0516 |
| DRB1*15:01 |        1 |                    201 |               0.0900 |  0.0635 |  0.1218 |

Japan baseline (same metric):

| allele     |   n_pops |   n_pooled_individuals |   weighted_mean_freq |   ci_lo |   ci_hi |
|:-----------|---------:|-----------------------:|---------------------:|--------:|--------:|
| DPB1*05:01 |       10 |                   5632 |               0.3584 |  0.3495 |  0.3673 |
| A*02:07    |        6 |                  20303 |               0.0341 |  0.0324 |  0.0359 |
| C*03:02    |        4 |                  20553 |               0.0057 |  0.0050 |  0.0065 |
| DQB1*03:02 |        9 |                   4723 |               0.1090 |  0.1028 |  0.1155 |
| B*46:01    |        5 |                  20160 |               0.0476 |  0.0455 |  0.0497 |

China baseline (same metric):

| allele     |   n_pops |   n_pooled_individuals |   weighted_mean_freq |   ci_lo |   ci_hi |
|:-----------|---------:|-----------------------:|---------------------:|--------:|--------:|
| DQB1*03:02 |       28 |                   7439 |               0.0607 |  0.0569 |  0.0647 |
| DPB1*05:01 |       20 |                   1969 |               0.3554 |  0.3405 |  0.3707 |
| B*46:01    |       20 |                   9554 |               0.1207 |  0.1161 |  0.1254 |
| DRB1*15:01 |       16 |                  11047 |               0.0998 |  0.0959 |  0.1038 |

## 8 Limitations

- **Sample-size disparity**: AFND populations range from N≈50 to N>70,000 (USA NMDP Korean). Weighted-mean pooling reduces, but does not eliminate, study-specific bias (HLA typing platform, study era, regional ancestry sub-structure within a country code).
- **Sub-population resolution**: "China" pools 30+ ethnic groups whose HLA frequencies span almost the full East-Asian range (e.g., HLA-B*46:01 from 0.013 in Tibet to 0.254 in Yunnan Dai). Country-level summaries are useful as anchors but unreliable as point estimates.
- **2-field resolution**: 4-field allele names are collapsed; this is appropriate for cross-population comparison but loses sub-allelic structure relevant for fine-mapping.
- **Locus coverage**: DQA1 / DPA1 are excluded from the master matrix (sparse in AFND). DRB3/4/5 also not separately fetched.
- **Fetch budget cap (100 GETs)** truncates pagination for some country/locus pairs (page 1 only for several entries — see `fetch_manifest.json`). Top alleles are well-covered; rare alleles in long tail may be missing for some populations.
- **No haplotype phase**: AFND single-allele tables only. Linkage disequilibrium analyses would require separate haplotype-frequency endpoints (out of scope here).
- **Cancer-FREE by construction**: this atlas is healthy-population only. Any downstream cancer-vs-baseline contrast must follow Paper 2/4 boundary rules and live in those papers, not here.

## 9 Honest assessment — atlas as a tools paper

**Verdict: solid Pan-Asian methodological scaffold; not yet a standalone tools paper without two further additions.**

What is reviewer-ready here:
1. Reproducible AFND ingestion (cache + parse + budget-aware fetcher) — valuable as a community resource.
2. Sample-size-weighted Clopper-Pearson CIs per cell (most prior pan-Asian HLA descriptions report point estimates only).
3. Reynolds genetic-distance + Jensen-Shannon dendrogram on a single per-locus-normalized matrix that spans East / SE / South Asia (16 countries, 245 populations, ≥4-locus subset = 70 populations after filter) computed identically across loci. European outgroup overlay was budgeted but not reached due to AFND throttling and would need ~18 additional GETs.
4. Cross-anchored Korean baseline (Korea × Japan × Han × East-Asia pool) directly answers a Paper 2/4 reviewer question.

What is missing for a full Nature Communications / Genome Medicine tools paper:
1. **Haplotype layer**: AFND has haplotype tables (e.g., A-B-DRB1) that this fetch did not pull. A pan-Asian healthy haplotype atlas with LD decay by region is the bigger contribution.
2. **Imputation accuracy benchmark**: ideally validate the population frequencies against Korean K2 arcasHLA (n=260) and a published Han imputation set, to show the atlas is consistent with NGS-typed cohorts. (Track 8 only does the AFND side.)
3. **Fine-grained sub-population resolution**: country-level pooling masks the substructure that motivates the atlas in the first place. A multi-pop dendrogram per country is doable but absent here.
4. **Population-coverage map / interactive viewer**: tools papers need a deployable artifact (e.g., a static HTML viewer with allele search). Not built here.

**Recommended next move**: keep this Track 8 atlas as Paper 2 / Paper 4 supplementary infrastructure (citable methodological appendix), and split off the 4-item list above into a future dedicated tools manuscript only if Korean/Pan-Asian haplotype data warrants it.
