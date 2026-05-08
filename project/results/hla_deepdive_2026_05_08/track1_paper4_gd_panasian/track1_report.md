# Track 1 — Paper 4: Pan-Asian Graves' disease HLA susceptibility deep-dive (v2)

**Date:** 2026-05-08
**Author:** Seungho Cook (Track 1 deep-dive)
**Scope:** Autoimmune Graves' disease (GD) HLA susceptibility, Pan-Asian (Korean / Han Chinese / Taiwanese / Japanese).
**Status:** Reviewer-grade upgrade of the prior Pan-Asian forest (`paper4_chu2018_gd_anchor_forest.tsv` and `panasian_forest_meta.tsv` superseded for forest analytics; raw count tables retained).

---

## 0. Boundary statement (read first)

This track is **autoimmune-only**. No claim about thyroid cancer susceptibility, prognosis, OS / DSS / PFI / RAI response, DM1 / DM2 status, BRAF / RAS / TERT status, or HLA-mediated tumour escape is made anywhere in this report. The analyses below treat HLA alleles solely as autoimmune Graves' disease risk loci, in line with the cancer-separation contract at `project/paper2_hla_boundary/HLA_CANCER_SEPARATION_RULES.md`. Any reader looking for HLA × thyroid cancer should consult Paper 1 / Paper 2 / Paper 3 — those analyses are deliberately walled off.

---

## 1. Cohort assembly

The Pan-Asian GD case-control evidence base used here is:

| # | Source | Country | Ancestry | n_GD | n_ctrl | Allele rows | Citation key |
|---|--------|---------|----------|-----:|-------:|------------:|--------------|
| 1 | Chu 2018 | China  | Han Chinese        | 1,468 | 1,490 | 6 | [Chu 2018 J Med Genet 55:685] |
| 2 | Shin 2019 | Korea | Korean pediatric    |    71 |   142 | 3 | [Shin 2019 Korean GD pediatric] |
| 3 | Chen 2011 | Taiwan | Taiwan ethnic Chinese | – |  – | 2 | [Chen 2011 Taiwan GD Bonferroni-corrected] |
| 4 | Inoue 1992 | Japan | Japanese            |   124 |   150 | 1 | [Inoue 1992 J Clin Endocrinol Metab 75:1444] |
| 5 | Dong 1992 | Japan | Japanese            |   121 |   140 | 1 | [Dong 1992 Tissue Antigens 39:185] |
| 6 | Park 2005 | Korea | Korean adult        |    88 |   104 | 1 | [Park 2005 Korean GD HLA-DPB1] |
| 7 | Cho 1987 | Korea | Korean adult         |    95 |   178 | 1 | [Cho 1987 Korean serology, B-locus only] |

Per-cohort source-rows table with inverse-variance weights:
`tables/T01_panasian_GD_source_rows_v2.tsv` (15 rows, 6 alleles).

Population-baseline anchor (no GD case-control): AFND South Korea, China, Japan, Taiwan pools, summarized in `tables/T08b_AFND_country_weighted_summary.tsv`. Korean baseline coverage at AFND v3.0.0 is sparse for this allele subset (only 5 of the 8 focus alleles have any AFND South-Korea entries), which is itself a documentation gap, not a finding.

K2 South-Korean NGS reference for LD context: `K2_arcasHLA_genotypes.tsv` (n=63 four-digit-typed RNA-seq germline calls; PRJEB11591 Yoo 2016). Treated as a healthy / non-GD ancestry-matched HLA reference, not as GD evidence.

---

## 2. Meta-analytic forest (random-effects, DerSimonian-Laird)

Random-effects pooled OR per allele (`tables/T02_panasian_GD_DL_random_effects_v2.tsv`). Forest plot:
`plots/forest/forest_panasian_GD_v2.{png,pdf,svg}`. Per-allele pages: `plots/per_allele/forest_<allele>.png` (6 alleles).

Headline numbers (Pan-Asian random-effects pooled OR for GD susceptibility):

| Allele | k | Pooled OR | 95% CI | I² (%) | τ² | 95% PI | Sources |
|--------|--:|----------:|--------|------:|---:|--------|---------|
| **DPB1\*05:01** | **5** | **2.01** | **1.75 – 2.32** | **22** | 0.006 | **1.50 – 2.71** | Chu / Shin / Chen / Inoue / Park |
| **B\*46:01**    | **4** | **2.08** | **1.34 – 3.23** | **83** | 0.137 | 0.52 – 8.26 | Chu / Shin / Chen / Cho |
| **A\*02:07**    | 2 | 2.06 | 1.70 – 2.50 | 0 | 0 | – | Chu / Dong |
| **C\*01:02**    | 2 | 1.85 | 1.59 – 2.14 | 0 | 0 | – | Chu / Shin |
| DQB1\*02:01 | 1 | 0.57 | 0.49 – 0.66 | – | – | – | Chu (anchor only) |
| DRB1\*07:01 | 1 | 0.43 | 0.36 – 0.51 | – | – | – | Chu (anchor only) |

**DPB1\*05:01** is the meta-analytically firmest Pan-Asian GD risk allele: pooled OR 2.01 (1.75 – 2.32), I² 22%, τ² 0.006, 95% prediction interval 1.50 – 2.71, all five contributing cohorts effect-size compatible. **B\*46:01** is also pooled-significant (OR 2.08, p=0.0012) but heterogeneous (I² 83%, τ² 0.14) because Chen 2011 Taiwan reports a much smaller effect (OR 1.33) than Chu 2018 Han (OR 2.38) and Shin 2019 Korean pediatric (OR 3.96). The wide 95% prediction interval (0.52 – 8.26) means the true B\*46:01 effect in a *new* Asian cohort is not yet pinned down.

The protective alleles (DRB1\*07:01, DQB1\*02:01, DQA1\*02:01) remain Chu-anchor-only; Korean / Japanese replication for these is the next gap.

The **previous "Chen 2018 Front Endocrinol 9:467" rows in `panasian_forest_meta.tsv` were a misattribution** of Chu 2018 J Med Genet 55:685 (DEPRECATION audit 2026-05-01). The v2 forest in this report uses Chu 2018 verified counts and the deprecated rows are not propagated.

---

## 3. Sub-population stratified meta

Per-ancestry random-effects pooled OR per allele (`tables/T04_subpopulation_stratified_meta.tsv`, `plots/stratified_panasian_GD.png`):

| Allele | Korean | Han Chinese | Taiwan | Japanese | POOLED |
|--------|-------:|------------:|-------:|---------:|-------:|
| DPB1\*05:01 | **2.70 (1.27 – 5.73)** k=2 | 1.90 (1.69 – 2.14) | 2.34 (1.80 – 3.05) | 1.78 (1.28 – 2.48) | 2.01 (1.75 – 2.32) |
| B\*46:01    | **2.78 (1.56 – 4.95)** k=2 | 2.38 (1.99 – 2.85) | 1.33 (1.07 – 1.66) | – | 2.08 (1.34 – 3.23) |
| A\*02:07    | – | 2.10 (1.70 – 2.59) | – | 1.86 (1.15 – 3.02) | 2.06 (1.70 – 2.50) |
| C\*01:02    | 2.51 (1.04 – 6.04) k=1 | 1.83 (1.57 – 2.13) | – | – | 1.85 (1.59 – 2.14) |

Notable: the **Korean-only sub-meta produces the highest point estimates** for both DPB1\*05:01 (OR 2.70) and B\*46:01 (OR 2.78), but those rest on Shin 2019 (n=71 pediatric cases) plus a single adult Korean source (Park 2005 / Cho 1987). The 95% CIs are wide. This is consistent with — but does not yet establish — a stronger Korean-specific DPB1\*05:01 / B\*46:01 effect than Han Chinese. Definitive resolution requires a properly powered Korean adult GD case-control NGS HLA cohort (Section 8 / Section 10).

---

## 4. Carrier-frequency vs allele-frequency reconciliation (HWE)

`tables/T05_carrier_vs_allele_HWE_reconciliation.tsv` reconciles Chu 2018's reported carrier frequencies (% of subjects carrying ≥1 copy) against allele frequencies under Hardy-Weinberg.

| Allele | Carrier f (GD) | Carrier f (ctrl) | q (GD) HWE | q (ctrl) HWE | Carrier OR | Allele OR (HWE) |
|--------|---------------:|-----------------:|-----------:|-------------:|-----------:|----------------:|
| DPB1\*05:01 | 0.440 | 0.313 | 0.252 | 0.171 | 1.90 | 1.63 |
| B\*46:01    | 0.141 | 0.065 | 0.073 | 0.033 | 2.38 | 2.31 |
| A\*02:07    | 0.097 | 0.049 | 0.050 | 0.025 | 2.10 | 2.06 |
| C\*01:02    | 0.184 | 0.109 | 0.097 | 0.056 | 1.83 | 1.80 |
| DQB1\*02:01 | 0.109 | 0.178 | 0.056 | 0.093 | 0.57 | 0.58 |
| DRB1\*07:01 | 0.071 | 0.153 | 0.036 | 0.080 | 0.43 | 0.43 |

**Reading**: for B\*46:01, A\*02:07, C\*01:02 and the protective DR/DQ pair, carrier-OR and allele-OR (HWE-back-converted) agree within rounding. For DPB1\*05:01 they diverge (1.90 vs 1.63) because DPB1\*05:01 has high carrier frequency (44 % GD / 31 % ctrl) where the HWE conversion is increasingly non-linear and homozygote contribution becomes large. **Operationally:** Pan-Asian GD HLA reviews should report DPB1\*05:01 as a carrier-frequency OR (1.90) and explicitly note that the underlying allele-frequency OR is 1.63 — the larger carrier-OR partly reflects allele frequency, not extra biology. This footnote was missing in the prior anchor table.

---

## 5. DPB1 sub-allele decomposition (K2 NGS Korean reference)

Is the GD signal specifically DPB1\***05:01** or any DPB1\***05:XX** allele? Decomposition in K2 (n=63 four-digit-typed Korean NGS): `tables/T06_DPB1_subAllele_decomposition_K2.tsv` and `tables/T06b_DPB1_top10_subAlleles_K2.tsv`.

In the K2 South-Korean NGS reference, **the entire DPB1\*05:XX family is DPB1\*05:01 — there are no DPB1\*05:02 / 05:03 / etc. detected** at four-digit resolution. n_chromosomes(DPB1\*05:01) = 44 / 116; allele frequency 0.379; carrier frequency 0.587. Top sub-alleles after DPB1\*05:01 are DPB1\*02:01 (0.224), DPB1\*04:01 (0.078), DPB1\*04:02 (0.069), DPB1\*03:01 (0.069), DPB1\*02:02 (0.043), DPB1\*13:01 (0.034). **Implication:** in Korean populations, "DPB1\*05" essentially equals "DPB1\*05:01"; the GD risk allele cannot be empirically separated from the family within current Korean NGS, so reporting at the 4-digit DPB1\*05:01 level is correct and not granularity-inflated.

(Caveat: the K2 NGS reference is small for sub-allele estimation. Definitive sub-allele decomposition would require a larger Korean NGS HLA panel, e.g. Korean Genome Project KOHLA + KoGES.)

---

## 6. Publication-bias diagnostics

`plots/galbraith_panasian_GD.png` (Galbraith) and `plots/funnel_DPB1_0501.png` (DPB1\*05:01-specific funnel).

- The Galbraith plot shows all 15 study-allele points stratified by precision; the GD-risk alleles (DPB1\*05:01, B\*46:01, A\*02:07, C\*01:02) sit in the upper half (z > +1.96) and the protective alleles (DRB1\*07:01, DQB1\*02:01) sit in the lower half (z < -1.96), with no obvious precision-asymmetry that would indicate small-study bias. Chu 2018 dominates precision for every allele — high precision, high z — which is expected given its sample size (n≈3,000).
- The DPB1\*05:01 funnel plot has all five points (Chu / Shin / Chen / Inoue / Park) within a triangular envelope around the RE pooled estimate, no obvious asymmetry — but with k=5 a formal Egger test is underpowered, so we restrict ourselves to a visual conclusion: **no flagrant small-study or selective-reporting pattern visible.**

A formal Begg / Egger asymmetry test was not performed because k≤5 per allele renders such tests unreliable.

---

## 7. LD context (K2 NGS Korean reference)

Carrier-level pairwise r² across eight focus alleles in K2 (n=63): `tables/T07_LD_matrix_K2_carrier_level.tsv` and `plots/LD_heatmap_K2.png`.

Top pairs in K2 Korean reference:

| allele 1 | allele 2 | r² | Fisher p | Joint carrier f |
|----------|----------|---:|--------:|----------------:|
| B\*46:01    | A\*02:07    | 0.337 | 0.0019 | 0.048 |
| A\*02:07    | DRB1\*08:03 | 0.217 | 0.018  | 0.032 |
| B\*46:01    | DRB1\*08:03 | 0.129 | 0.042  | 0.032 |
| C\*01:02    | DRB1\*15:01 | 0.127 | 0.017  | 0.063 |
| DPB1\*05:01 | B\*46:01   | 0.077 | 0.073  | 0.016 |

**Reading:** B\*46:01 and A\*02:07 carriers co-occur substantially in K2 Koreans (r²≈0.34). C\*01:02 is partly tagging the same B-locus haplotype block. **DPB1\*05:01 is not in carrier-level LD with B\*46:01 in K2 (r²≈0.08, p=0.07)** — the two GD risk alleles segregate independently in this Korean reference. That is consistent with their literature interpretation as separate independent GD susceptibility loci (DPB1 = class II antigen presentation, B-locus = class I) rather than a single haplotype effect.

Caveat: K2 n=63 four-digit-typed subjects is small; LD estimates have wide CIs. A larger Korean NGS reference (e.g. KOHLA full release) would tighten this.

---

## 8. Reviewer-grade open questions

1. **Is the DPB1\*05:01 effect Korean > Han Chinese?** The Korean sub-meta point estimate (OR 2.70, k=2) is higher than Chu 2018 Han Chinese (1.90), but the Korean 95% CI (1.27 – 5.73) is wide and overlaps the Han point. Resolution: a properly powered Korean adult GD NGS HLA cohort.
2. **B\*46:01 heterogeneity.** Why is Taiwan Chen 2011 OR 1.33 but Korean Shin 2019 OR 3.96? Possible answers: (a) genuinely different B\*46:01 LD blocks across Asian populations; (b) Bonferroni-corrected p in Chen 2011 attenuates the effect; (c) Shin 2019 small-sample inflation. A Pan-Asian B\*46:01 meta should not be reported as a single-effect estimate without flagging the I²=83 % spread.
3. **Are protective alleles (DRB1\*07:01, DQB1\*02:01, DQA1\*02:01) replicated in Korean GD?** Currently anchor-only. A Korean adult GD cohort with full DR/DQ typing would close this.
4. **DPA1 axis.** Chu 2018 reports DPA1\*02:02 as a top GD risk allele (OR 1.90), but no Korean replication source row exists in our table. DPA1 typing is non-routine in older studies.
5. **Functional bridge.** netMHCIIpan binding of TSHR peptides on HLA-DPA1\*02:02/DPB1\*05:01 (already computed, see `project/results/v17_korean/A2_netmhciipan/`) provides a candidate molecular hypothesis for *why* DPB1\*05:01 is a GD risk allele; but the deep-dive here makes no claim of mechanism, only association.

---

## 9. Limitations

- **k is small** (1–5 per allele); funnel / Egger tests are not powered.
- **SE for older studies** (Shin 2019, Chen 2011, Inoue 1992, Dong 1992, Cho 1987, Park 2005) is approximated from reported p-values via the normal-z relation when the original CI was not reported. Underlying counts are documented in `tables/T01_panasian_GD_source_rows_v2.tsv`.
- **No Japanese B\*46:01 row** is available in the curated source table at this time (Onuma 1994 mentions B\*46:01 but does not formally report the test). This is a documentation gap, not a negative result.
- **K2 NGS reference n=63** is small for LD and sub-allele decomposition; conclusions are direction-only.
- **AFND coverage gap**: only 5 of 8 focus alleles have a Korean South-Korea AFND row in the local snapshot. Three (DRB1\*08:02, DRB1\*16:02, DQB1\*03:02) lack any Korean AFND baseline in the cached JSON; for those, no Korean-specific population frequency is anchored.
- **No GD cases were typed by us in this track.** All GD-side numbers come from literature counts/ORs; we built only the meta-analytic synthesis and the Korean NGS reference context.

---

## 10. Recommended next study

**A Korean adult Graves' disease HLA case-control NGS cohort, n≥250 cases, n≥500 ancestry-matched controls, NGS HLA typing across all 8 classical loci.** Such a study would:

1. Replicate DPB1\*05:01 and B\*46:01 with adult-onset Korean cases (Shin 2019 is pediatric-only).
2. Resolve the Korean-vs-Han DPB1\*05:01 effect-size question.
3. Provide Korean replication of the protective DR/DQ axis.
4. Allow joint-locus haplotype analysis (B\*46:01–C\*01:02 block, DPB1\*05:01 vs DR/DQ extended haplotype).
5. Generate the substrate for a Korean-led Pan-Asian Graves' GWAS / immunochip integration.

This study is paper-blocking for upgrading the Paper 4 readiness scoreboard (Section 11 / `tables/T09_replication_readiness_scoreboard_v2.tsv`) from "A_full_panasian_meta_ready" to "A+_korean_adult_anchor_added".

---

## 11. Replication-readiness scoreboard v2

Full table: `tables/T09_replication_readiness_scoreboard_v2.tsv`. Summary:

| Allele | Grade | Bottleneck |
|--------|-------|------------|
| DPB1\*05:01 | A — full Pan-Asian meta ready | Korean adult NGS anchor (A→A+) |
| B\*46:01    | A — full Pan-Asian meta ready (with I² caveat) | Japanese row missing; Korean adult anchor for I² resolution |
| A\*02:07    | B — partial (Han + Japan) | Korean replication missing |
| C\*01:02    | B — partial (Han + Korea pediatric) | Taiwan / Japanese replication |
| DRB1\*07:01 | C — Chu anchor only | Korean + Taiwan + Japanese replication |
| DQB1\*02:01 | C — Chu anchor only | Korean + Taiwan + Japanese replication |
| DRB1\*08:02, DRB1\*15:01, DRB1\*16:02, C\*03:02, DQB1\*03:02 | D — underpowered (no anchor) | Need anchor cohort |

---

## 12. Output index

```
project/results/hla_deepdive_2026_05_08/track1_paper4_gd_panasian/
├── track1_report.md                                     ← this file
├── track1_summary.json
├── tables/
│   ├── T01_panasian_GD_source_rows_v2.tsv
│   ├── T02_panasian_GD_DL_random_effects_v2.tsv
│   ├── T03_leave_one_out_sensitivity.tsv
│   ├── T04_subpopulation_stratified_meta.tsv
│   ├── T05_carrier_vs_allele_HWE_reconciliation.tsv
│   ├── T06_DPB1_subAllele_decomposition_K2.tsv
│   ├── T06b_DPB1_top10_subAlleles_K2.tsv
│   ├── T07_LD_matrix_K2_carrier_level.tsv
│   ├── T08_AFND_panasian_baselines.tsv
│   ├── T08b_AFND_country_weighted_summary.tsv
│   └── T09_replication_readiness_scoreboard_v2.tsv
└── plots/
    ├── forest/forest_panasian_GD_v2.{png,pdf,svg}
    ├── per_allele/forest_<allele>.png    (6 alleles)
    ├── stratified_panasian_GD.png
    ├── galbraith_panasian_GD.png
    ├── funnel_DPB1_0501.png
    └── LD_heatmap_K2.png
```

---

*Track 1 closes here. No cancer claim is made anywhere in this report.*
