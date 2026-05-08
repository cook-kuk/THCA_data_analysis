# Track 10 — Korean published HLA literature extraction → Pan-Asian GD HLA forest v3

**Date:** 2026-05-08
**Scope:** AUTOIMMUNE-ONLY (Graves' disease, Hashimoto's disease, AITD).
**Boundary:** No thyroid cancer claim. No HLA × cancer outcome. Distinguish allele frequency from carrier frequency in every table. See `project/paper2_hla_boundary/HLA_CANCER_SEPARATION_RULES.md`.
**Built on:** Track 1 v2 forest (Pan-Asian GD HLA random-effects meta).

---

## 0. Boundary statement

This track curates Korean published HLA case-control literature for Graves' disease (and where reported, Hashimoto's disease / pooled AITD). All effect estimates pertain to autoimmune thyroid disease vs ancestry-matched controls. The five papers extracted here are autoimmune by design (the cohorts are GD or pediatric AITD); no cancer cohort is involved at any point. Allele frequency vs carrier frequency is recorded explicitly in every row of `T03_korean_lit_extraction_master.tsv`.

---

## 1. PMID / DOI resolution (5 / 5)

All five target papers were resolved on the first PubMed/Crossref pass.

| Paper | PMID | DOI | Title (truncated) | Journal · Year · Vol · Pages | Cohort |
|---|---|---|---|---|---|
| Shin 2019 | 31091281 | 10.1371/journal.pone.0216941 | HLA alleles, especially amino-acid signatures of HLA-DPB1, might contribute to the molecular pathogenesis of early-onset autoimmune thyroid disease | PLoS ONE · 2019 · 14 · e0216941 | Korean pediatric AITD (n_GD=71, n_HD=45, n_AITD=116, n_CTRL=142) |
| Cho 2011 | 21952423 | 10.1159/000331134 | Association of HLA alleles with autoimmune thyroid disease in Korean children | Horm Res Paediatr · 2011 · 76 · 328-334 | Korean pediatric AITD (n_GD=41, n_HD=32, n_CTRL=159) |
| Park 2005 | 15993720 | 10.1016/j.humimm.2005.03.001 | Association of HLA-DR and -DQ genes with Graves disease in Koreans | Hum Immunol · 2005 · 66 · 741-747 | Korean adult GD (n_GD=198, n_CTRL=200) |
| Jang 2011 | 21062236 | 10.3109/08820139.2010.525571 | Identification of HLA-DRB1 alleles associated with Graves' disease in Koreans by sequence-based typing | Immunol Invest · 2011 · 40 · 172-182 | Korean adult GD (n_GD=133, n_CTRL=200) |
| Baek 2021 | (not on PubMed) | 10.1111/tan.14134 | Distributions of HLA-A, -B, and -DRB1 alleles typed by amplicon-based next-generation sequencing in Korean volunteer donors for unrelated hematopoietic stem cell transplantation | HLA · 2021 · 97(2) · 112-126 | Korean unrelated HSCT volunteers (KONOS, n=26,202; class-I + DRB1 6-digit; population reference, NOT case-control) |

Full candidate scoreboard: `tables/T01_pmid_doi_candidates.tsv` (55 candidate hits across 10 query phrasings × 5 papers).
Authoritative citations: `tables/T02_resolved_citations.tsv`.
Raw resolution payload: `pmid_resolution_candidates.json`.

**Open-access status:** Shin 2019 is gold OA via PLoS (full PDF + S1 Table DOCX + S1 Fig TIFF retrieved). Cho 2011 / Park 2005 / Jang 2011 / Baek 2021 are paywalled (Karger, Elsevier, Taylor & Francis, Wiley); we extracted from PubMed-hosted full-abstract text only, which preserves the headline allele-level numbers each paper chose to publicize. Sci-Hub was NOT used.

---

## 2. Per-paper extractions

### 2.1 Shin 2019 — Korean pediatric AITD (PMID 31091281)

**Source:** S1 Table (`Shin_2019_S1_table.docx`) — "Allele-carrying frequencies of HLA-A, -B, -C, -DRB1, -DQB1 and -DPB1 alleles associated with GD or HD in Korean children with AITD (P < 0.05)."

**Typing:** Sanger PCR-SBT, high-resolution 4-digit, 6 loci.

**Frequency type:** carrier frequency (≥1 of 2 chromosomes carries the allele).

**Cohort layout:** Controls (n=142) | AITD (n=116, =GD+HD) | GD (n=71) | HD (n=45). 4 OR / 95% CI / p / Pc columns: Controls vs AITD; Controls vs GD; Controls vs HD; GD vs HD.

**Extracted:** **126 allele rows × 3 case-vs-control comparisons = 378 long-format rows**, saved to `extraction_worksheets/shin2019_long.tsv`. Wide format: `extraction_worksheets/shin2019_S1_table_full.tsv`.

**Headline numbers for Track-1 focus alleles (Controls vs GD):**

| Allele | GD carriers/N | CTRL carriers/N | OR | 95% CI | p | Pc |
|---|---|---|---|---|---|---|
| DPB1*05:01 | 63/71 | 88/142 | 4.60 | 2.07–10.21 | 0.0002 | 0.003 |
| B*46:01 | 24/71 | 16/142 | 3.96 | 1.94–8.09 | 0.0002 | 0.008 |
| A*02:07 | 14/71 | 8/142 | 3.99 | 1.59–10.01 | 0.003 | — |
| DRB1*08:02 | 5/71 | 10/142 | 1.04 | 0.35–3.15 | 0.94 | — |
| DRB1*15:01 | 13/71 | 21/142 | 1.30 | 0.61–2.78 | 0.49 | — |
| DRB1*16:02 | 0/71 | 5/142 | 0.18 | 0.01–4.22 | 0.28 | — |
| C*03:02 | 7/71 | 16/142 | 0.89 | 0.35–2.26 | 0.81 | — |
| DQB1*03:02 | 10/71 | 30/142 | 0.63 | 0.29–1.37 | 0.24 | — |

### 2.2 Cho 2011 — Korean pediatric AITD (PMID 21952423)

**Source:** PubMed abstract verbatim text (PDF paywalled).

**Typing:** PCR-SSP, 2-digit allele family resolution (HLA-A*02 etc).

**Frequency type:** allele family frequency (low resolution).

**Cohort:** n_HD=32, n_GD=41 (n_AITD=73), n_CTRL=159.

**Reported direction (no numeric ORs in abstract):**
- ↑ in AITD: HLA-A*02, B*46, Cw*01, DRB1*08
- ↓ in AITD: HLA-A*30, B*07, Cw*07, DRB1*01
- ↑ in HD: B*46, Cw*01
- ↓ in HD: DRB1*01, Cw*07
- ↑ in GD: A*02, B*46, Cw*01, DRB1*08
- ↓ in GD: DRB1*07, Cw*07
- Composite carrier statement: B*46 + Cw*01 co-carriage > either alone for AITD risk.

**Extracted:** 18 directional rows (no numeric meta inputs) → `extraction_worksheets/cho2011_long.tsv`. Cannot be added to a 4-digit forest as-is. Manual table extraction from the paywalled PDF is required for inclusion (audit A5).

### 2.3 Park 2005 — Korean adult GD (PMID 15993720)

**Source:** PubMed abstract verbatim text (PDF paywalled).

**Typing:** PCR-SSO; HLA-DRB1 + HLA-DQB1 ONLY (no DPB1, contra Track 1 v2).

**Frequency type:** allele frequency (per-chromosome). 2N_case=396, 2N_ctrl=400.

**Cohort:** n_GD=198, n_CTRL=200 (NOT 88/104 as Track 1 v2 mistakenly recorded; audit A2).

**Headline numbers from abstract:**
- DRB1*08:03 — GD 27.8% vs CTRL 14.5%, OR=2.27, Pc=0.03 (**susceptibility**)
- DRB1*16:02 — GD 5.1% vs CTRL 0%, OR=22.34, Pc=0.03 (**susceptibility, rare in controls**)
- DRB1*03:01 — male-only stratum: 12.5% vs 3.5%, OR=3.57, p<0.05 (excluded from main forest)
- DRB1*01:01, DRB1*07:01, DRB1*12:02, DRB1*13:02 — OR<0.5, p<0.05 (resistance; exact OR/CI not in abstract)

**Extracted:** 7 rows (3 with full numeric ORs; 4 weak-resistance rows OR/CI=NaN) → `extraction_worksheets/park2005_long.tsv`.

### 2.4 Jang 2011 — Korean adult GD (PMID 21062236)

**Source:** PubMed abstract verbatim text (PDF paywalled).

**Typing:** PCR-SBT high-resolution; HLA-DRB1 6-digit (e.g. DRB1*030101).

**Frequency type:** allele frequency (per-chromosome). 2N_case=266, 2N_ctrl=400.

**Cohort:** n_GD=133, n_CTRL=200.

**Headline numbers from abstract (raw p; corrected p NS for all):**
- DRB1*030101: 4.9% vs 1.8%, p=0.034 (↑ GD)
- DRB1*080201: 5.3% vs 2.3%, p=0.050 (↑ GD)
- DRB1*140301: 3.4% vs 1.0%, p=0.043 (↑ GD)
- DRB1*070101: 3.0% vs 7.3%, p=0.024 (↓ GD)
- DRB1*130201: 4.1% vs 9.0%, p=0.010 (↓ GD)

**Extracted:** 5 DRB1 rows (4-digit collapsed to allow alignment with Track 1 v2) → `extraction_worksheets/jang2011_long.tsv`. 6-digit retained in `allele_6d` column.

### 2.5 Baek 2021 — Korean class-I + DRB1 NGS reference (DOI 10.1111/tan.14134)

**Source:** Crossref full abstract (PDF paywalled, not on PubMed).

**Typing:** Amplicon-based MiSeqDx NGS; 6-digit; HLA-A, -B, -DRB1 ONLY.

**Cohort:** n=26,202 Korean unrelated HSCT volunteers (KONOS).

**Note (citation correction A4):** The data registry hint had described Baek 2021 as "Korean class-II NGS reference" suggesting DPB1/DPA1/DQA1/DQB1 coverage. The actual paper reports HLA-A, -B, -DRB1 only — class-I + DRB1, not class-II proper. **Baek 2021 is therefore NOT a DPB1*05:01 baseline source**; it does support HLA-A/B/DRB1 baselines.

**Extracted:** 17 well-known top alleles per locus from the paper's abstract / KONOS frequency atlas, saved to `tables/T08_Baek2021_class_I_DRB1_baseline.tsv` (caveat: requires PDF retrieval to populate the full 70 + 102 + 69 allele list with exact frequencies).

---

## 3. Pan-Asian GD HLA forest v3 — what changed vs v2

Source rows: **27 in v3** (vs 15 in v2). Source paper count: **8 in v3** (vs 6). Allele coverage: **15 alleles in v3** (vs 11). Detailed v2 vs v3 diff table: `tables/T06_v2_vs_v3_meta_diff.tsv`.

**Focus-allele meta v3 (DerSimonian-Laird random-effects):**

| Allele | k v2 → v3 | OR v2 | OR v3 | 95% CI v3 | I² v3 | Status v3 |
|---|---|---|---|---|---|---|
| **DPB1*05:01** | 5 → 4 | 2.01 | **2.10** | 1.70–2.60 | 55% | A (1 row dropped: Park 2005 citation error A1; Shin 2019 row replaced with full-CI) |
| **B*46:01** | 4 → 4 | 2.08 | **2.17** | 1.39–3.36 | 85% | A (Shin 2019 row replaced with full-CI; OR widened upward) |
| **A*02:07** | 2 → 3 | 2.06 | **2.12** | 1.73–2.61 | 5% | A (Shin 2019 added; promoted from B → A) |
| **C*01:02** | 2 → 2 | 1.85 | **1.88** | 1.58–2.24 | 6% | A (Shin 2019 row replaced) |
| **DRB1*07:01** | 1 → 2 | 0.43 | **0.43** | 0.36–0.51 | 0% | A (Jang 2011 added — replicates Chu 2018 protective signal) |
| **DRB1*08:02** | 0 → 2 | — | **1.71** | 0.76–3.84 | 29% | C (Shin 2019 + Jang 2011 — heterogeneous) |
| **DRB1*15:01** | 0 → 1 | — | **1.30** | 0.61–2.78 | 0% | C (Shin 2019; null effect) |
| **DRB1*16:02** | 0 → 2 | — | **2.07** | 0.02–234 | 81% | C (Park 2005 OR=22.3 pediatric Shin 2019 OR=0.18 — strong age-effect divergence) |
| **C*03:02** | 0 → 1 | — | **0.89** | 0.35–2.26 | (k=1) | C (Shin 2019; null effect) |
| **DQB1*03:02** | 0 → 1 | — | **0.63** | 0.29–1.37 | 0% | C (Shin 2019; trend toward protection) |
| DRB1*08:03 | new | — | **2.27** | 1.59–3.24 | 0% | C (Park 2005; bona fide Korean adult susceptibility) |
| DRB1*03:01 | new | — | **2.88** | 1.14–7.33 | 0% | C (Jang 2011; consistent with classic Caucasian DRB1*03:01 GD risk) |
| DRB1*14:03 | new | — | **3.47** | 1.06–11.4 | 0% | C (Jang 2011) |
| DRB1*13:02 | new | — | **0.44** | 0.22–0.87 | 0% | C (Jang 2011 — protective) |

**Headline shift:**
- 5 of 8 D-grade focus alleles now have ≥1 Korean source row (DRB1*08:02, DRB1*15:01, DRB1*16:02, C*03:02, DQB1*03:02 → all upgraded D → C).
- A*02:07 promoted B → A (k=3 ancestries: China, Japan, Korea).
- DPB1*05:01 / B*46:01 retain A-grade with cleaner per-cohort SE estimates.
- B*46:01 I² remains 85% — the 4 cohorts span China (OR 2.38), Korea (OR 3.96), Korea historic (OR 2.34), Taiwan (OR 1.33) — Taiwan is the outlier on the low side.

Output forest plot: `plots/forest_panasian_GD_v3.{png,pdf,svg}`. Per-allele forests for all 15 alleles: `plots/per_allele/forest_<allele>_v3.png`.

Source rows v3: `tables/T04_panasian_GD_source_rows_v3.tsv`.
Meta v3: `tables/T05_panasian_GD_DL_random_effects_v3.tsv`.
v2 vs v3 diff: `tables/T06_v2_vs_v3_meta_diff.tsv`.

---

## 4. Korean-only sub-meta v2

Restricting to Korean rows only (Cho 1987, Park 2005, Jang 2011, Shin 2019), 14 of the 15 alleles still have a Korean point estimate. Significant Korean-only susceptibility / protection signals (p<0.05):

| Allele | Korean RE OR | 95% CI | Source(s) | Direction |
|---|---|---|---|---|
| **B*46:01** | 3.03 | 1.81–5.08 | Cho 1987 + Shin 2019 | susceptibility |
| **DPB1*05:01** | 4.60 | 2.07–10.22 | Shin 2019 (single) | susceptibility |
| **A*02:07** | 3.99 | 1.59–10.01 | Shin 2019 (single) | susceptibility |
| **C*01:02** | 2.51 | 1.40–4.50 | Shin 2019 (single) | susceptibility |
| **DRB1*08:03** | 2.27 | 1.59–3.24 | Park 2005 | susceptibility |
| **DRB1*03:01** | 2.88 | 1.14–7.33 | Jang 2011 | susceptibility |
| **DRB1*14:03** | 3.47 | 1.06–11.4 | Jang 2011 | susceptibility |
| **DRB1*13:02** | 0.44 | 0.22–0.87 | Jang 2011 | resistance |
| **DRB1*07:01** | 0.40 | 0.18–0.88 | Jang 2011 | resistance |

Output: `tables/T07_korean_only_sub_meta_v2.tsv`, `plots/forest_korean_only_v2.{png,pdf,svg}`.

---

## 5. Baek 2021 baseline triangulation

Baek 2021 is HLA-A/B/DRB1 ONLY (n=26,202 Korean KONOS NGS). Its locus coverage does **not** overlap with the DPB1*05:01 / DQB1*03:02 / DQB1*02:01 forest core. Where it does overlap (HLA-A, B, DRB1 alleles), Baek 2021 frequencies broadly concord with the AFND South Korea pool — exact pairwise concordance scatter saved to `plots/baek_vs_afnd_concordance.png` (output uses the verbatim subset of well-known top alleles cited in the abstract; full per-allele Baek 2021 table requires the paywalled PDF for the remaining ~200 alleles).

Output: `tables/T08_Baek2021_class_I_DRB1_baseline.tsv`, `tables/T09_Baek_vs_AFND_concordance.tsv`, `plots/baek_vs_afnd_concordance.png`.

---

## 6. Pediatric vs adult Korean AITD HLA pattern

Three alleles have both pediatric (Shin 2019) and adult (Park 2005 or Jang 2011 or Cho 1987) Korean rows:

| Allele | Pediatric OR | Adult OR | Δ log OR | p (diff) | Same direction? |
|---|---|---|---|---|---|
| **B*46:01** | 3.96 (Shin 2019) | 2.34 (Cho 1987) | −0.53 | 0.30 | **YES** (consistent susceptibility, weaker in adult cohort) |
| **DRB1*08:02** | 1.04 (Shin 2019) | 2.41 (Jang 2011) | +0.84 | 0.24 | YES (both >1; significantly different not) |
| **DRB1*16:02** | 0.18 (Shin 2019) | 22.34 (Park 2005) | +4.82 | **0.022** | **NO — direction flip** |

The DRB1*16:02 flip is driven by 0/71 carriers in Shin's pediatric GD vs 5.1% chromosome frequency in Park's adult GD with 0% in adult controls. Both estimates are wide because of zero-cell instability; the pediatric "0.18" is essentially the Haldane-corrected lower bound of an unestimable cell. This is a KEY pediatric-vs-adult finding to flag for the GD HLA narrative: DRB1*16:02 may be an adult-onset GD signal that is silent in early-onset Korean AITD.

Output: `tables/T10_pediatric_vs_adult_korean.tsv`.

---

## 7. Replication-readiness scoreboard v3

| Allele | n_sources v3 | Grade v2 | Grade v3 | Changed? |
|---|---|---|---|---|
| DPB1*05:01 | 4 | A | A | no |
| B*46:01 | 4 | A | A | no |
| A*02:07 | 3 | B | **A** | **upgraded** |
| DRB1*08:02 | 2 | D | **C** | **upgraded** |
| DRB1*15:01 | 1 | D | **C** | **upgraded** |
| DRB1*16:02 | 2 | D | **C** | **upgraded** |
| C*03:02 | 1 | D | **C** | **upgraded** |
| DQB1*03:02 | 1 | D | **C** | **upgraded** |

Six of 8 focus alleles changed grade. No allele remains D. Heatmap: `plots/readiness_scoreboard_v3_heatmap.png`. Detail: `tables/T11_readiness_scoreboard_v3.tsv`.

**Definitions (v3):**
- A: Chu 2018 anchor + ≥3 ancestries (China + Japan + Taiwan or Korea) → meta-publishable.
- B: Chu 2018 anchor + 2 ancestries → close to A.
- C: ≥1 anchor (Chu / Korean pediatric / Korean adult) but <2 ancestries — needs +2 cohorts.
- D: no anchor — was true in v2 for 5 alleles; now zero in v3.

---

## 8. Limitations

1. **Paywalls bound resolution.** Cho 2011, Park 2005, Jang 2011, Baek 2021 — only PubMed/Crossref abstract text was usable. Cho 2011 yielded zero numeric meta inputs (directional only). Park 2005 / Jang 2011 yielded only the headline alleles (3 + 5 alleles). Full DRB1 / DQB1 frequency tables behind these paywalls are not yet machine-readable.
2. **Resolution heterogeneity.** Cho 2011 is 2-digit; Jang 2011 is 6-digit; Shin 2019 / Park 2005 are 4-digit; Baek 2021 is 6-digit. The forest harmonizes at 4-digit, which loses sub-allele information from Jang 2011 and Baek 2021 (sub-allele columns retained in the worksheets for future re-analysis).
3. **Carrier frequency vs allele frequency.** Shin 2019 reports carrier frequencies; Park 2005 / Jang 2011 report allele (per-chromosome) frequencies. The forest treats both as exchangeable on the OR scale, which is a standard assumption when the allele is rare (HWE → carrier OR ≈ allele OR for q<0.1) but may inflate effects for common alleles (DPB1*05:01, B*46:01). Track-1 v2 already addressed this via the HWE reconciliation table (T05) and the same caveat applies here. **Every row carries a `frequency_type` column.**
4. **DRB1*16:02 zero-cell instability.** Park 2005 reports 0% in controls; Shin 2019 reports 0/71 in pediatric GD. Both yield Haldane-corrected estimates with very wide CIs. The age-effect interpretation in Section 6 should be treated as hypothesis-generating, not confirmatory.
5. **Single-source rows do not benefit from random-effects shrinkage.** DRB1*15:01, DRB1*03:01, DRB1*14:03, DRB1*13:02, DQB1*03:02, C*03:02 are all k=1 in v3. Their CIs are study-specific Wald CIs; I² and τ² are undefined. They should not be over-interpreted as "Pan-Asian" signals.
6. **Baek 2021 incomplete.** Only the top ~17 alleles by frequency could be transcribed from the abstract; the full ~241-allele frequency table is behind the Wiley paywall. No DPB1 / DQB1 / DPA1 / DQA1 baseline can be sourced from Baek 2021 (locus coverage limitation, audit A4).
7. **No GD vs HD discrimination at the meta level.** Shin 2019 contributes a Controls_vs_GD row to the pan-Asian forest. The full Controls_vs_AITD and Controls_vs_HD tables (in `extraction_worksheets/shin2019_long.tsv`) are reserved for the HD-specific track and were not pooled into the GD forest to avoid mixing phenotypes.
8. **Shin 2019 cohort is pediatric only.** Adult Korean GD currently rests on Cho 1987 (B-locus serology), Park 2005 (DRB1/DQB1), and Jang 2011 (DRB1) — none of which fully replicate the DPB1*05:01 anchor. A Korean-adult-NGS GD cohort remains the highest-value next step (already named in Track 1 v2 next-steps).

---

## 9. Citation correction audit (A1–A6)

Six audit findings, recorded in `tables/T12_citation_audit.tsv`:

- **A1 (HIGH).** Track 1 v2 had a Park 2005 row for DPB1*05:01 (OR=2.05, n=88/104). Park 2005 reports DRB1+DQB1 only (n=198/200) — DPB1 was never typed. Row dropped from v3.
- **A2 (MEDIUM).** Track 1 v2 listed Park 2005 cohort as n=88/104; correct is n=198/200. v3 corrects.
- **A3 (MEDIUM).** Track 1 v2 derived Shin 2019 SEs from corrected p (Pc) — biased. v3 uses full 95% CIs from S1 Table.
- **A4 (LOW).** Manifest hint described Baek 2021 as "class-II NGS reference"; actual coverage is HLA-A/B/DRB1. Manifest update flagged.
- **A5 (INFO).** Cho 2011 is 2-digit PCR-SSP — cannot enter the 4-digit forest. Documented for future PDF-based re-curation.
- **A6 (INFO).** Jang 2011 6-digit collapsed to 4-digit for forest; 6-digit retained in worksheets.

Track 1's prior corrections (Chu 2018 ↔ "Chen 2018 Front Endocrinol" misattribution) were not re-flagged here, having already been resolved upstream.

---

## 10. Next steps

1. Acquire full PDFs for Cho 2011 / Park 2005 / Jang 2011 / Baek 2021 via institutional library; re-run the extractor on each full Table → expand from headline-allele only to full DRB1/DQB1 rows.
2. Korean-adult NGS GD cohort (planned Bundang SNUH outreach per memory `K2_vs_bundang_distinction`) — would push DRB1*08:02, DRB1*15:01, DRB1*16:02, C*03:02, DQB1*03:02 from C → A grade.
3. Manifest update: change Baek 2021 description from "class-II NGS reference" to "class-I + DRB1 NGS reference (HLA-A/B/DRB1, n=26,202)".
4. Update Track 1 v2 source-row file to mark the Park 2005 DPB1*05:01 row as DEPRECATED (or delete it) so downstream tools do not silently reuse it.
5. Pursue full DRB1 6-digit forest (use Jang 2011's `allele_6d` column + Baek 2021 6-digit baseline) — would let DPB1*05:01 narrative be accompanied by a parallel DRB1*030101 / DRB1*080201 sub-allele forest.

---

## Output index

```
project/results/hla_deepdive_2026_05_08/track10_korean_lit/
├── pmid_resolution_candidates.json       # raw PubMed/Crossref hits
├── supplement_fetch_log.json             # which fetches succeeded/failed
├── track10_report.md                     # this file
├── source_pdfs/
│   ├── Shin_2019.pdf                     # full OA PDF
│   ├── Shin_2019_S1_table.docx           # S1 carrier-frequency table (extracted)
│   └── Shin_2019_S1_fig.tif              # S1 ribbon-model figure
├── extraction_worksheets/
│   ├── shin2019_S1_table_full.tsv        # 126 allele-rows wide
│   ├── shin2019_long.tsv                 # 378 allele-comparison rows
│   ├── cho2011_long.tsv                  # 18 directional rows
│   ├── park2005_long.tsv                 # 7 numeric rows
│   └── jang2011_long.tsv                 # 5 numeric rows
├── tables/
│   ├── T01_pmid_doi_candidates.tsv
│   ├── T02_resolved_citations.tsv
│   ├── T03_korean_lit_extraction_master.tsv  # 390 unified rows
│   ├── T04_panasian_GD_source_rows_v3.tsv    # 27 rows
│   ├── T05_panasian_GD_DL_random_effects_v3.tsv  # 15 alleles × meta
│   ├── T06_v2_vs_v3_meta_diff.tsv
│   ├── T07_korean_only_sub_meta_v2.tsv
│   ├── T08_Baek2021_class_I_DRB1_baseline.tsv
│   ├── T09_Baek_vs_AFND_concordance.tsv
│   ├── T10_pediatric_vs_adult_korean.tsv
│   ├── T11_readiness_scoreboard_v3.tsv
│   └── T12_citation_audit.tsv
└── plots/
    ├── forest_panasian_GD_v3.{png,pdf,svg}
    ├── forest_korean_only_v2.{png,pdf,svg}
    ├── baek_vs_afnd_concordance.png
    ├── readiness_scoreboard_v3_heatmap.png
    └── per_allele/forest_<allele>_v3.png  (15 files)
```
