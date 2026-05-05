---
title: "Paper 2 Pillar I v2 — forest meta SPEC (Korean PTC vs AFND South Korea baseline)"
date: 2026-05-04
status: SPEC ONLY — no analysis execution. Awaiting advisor approval + Lee 2014 fill decision + Figure 1 re-render decision.
scope: Paper 2 (Hashimoto-overlap PTC) ONLY. GD/Graves = Paper 4 backlog (out of scope). BRAF/RAS/TERT/DM1 = Paper 1 (out of scope).
supersedes: v1 forest meta (Korean PTC vs Chu 2018 GD direct comparison) deprecated 2026-05-04
---

# Paper 2 Pillar I v2 — forest meta spec

---

## 1. Why v1 was deprecated

**v1 design (deprecated 2026-05-04)**:
- Forest meta directly compared Korean PTC pool (n=874) versus Chu et al. 2018 Han Chinese Graves' disease cohort (n=1,468)
- Top finding: "Korean PTC pool DPB1\*05:01 (53.2%) HIGHER than Chu Han Chinese GD itself (44.0%)"
- Framing: "Pan-Asian autoimmune-thyroid continuum confirmed"

**Why deprecated** (2 fundamental errors):

(a) **Cross-disease mismatch**. Korean PTC pool's autoimmune comorbidity is Hashimoto's thyroiditis (~20-30% within Korean PTC). Chu 2018 cohort is Graves' disease (n=1,468 GD vs n=1,490 controls). Direct frequency comparison places HT-substrate Korean PTC against GD-substrate Chinese cohort — two distinct autoimmune diseases with different mechanisms (HT = destructive lymphocytic infiltration → hypothyroidism; GD = TSH receptor stimulating antibody → hyperthyroidism). Reviewer would flag as cross-disease mismatch.

(b) **Ancestry baseline confound**. Korean DPB1\*05:01 baseline is reportedly ~38-42% (Lee et al. Tissue Antigens 2014, KOTRY donor cohort) — already higher than Chinese baseline ~38% (Chu 2018 controls). The "Korean PTC 53.2% > Chinese GD 44.0%" comparison conflates ancestry baseline differential with disease-specific elevation. The proper question is each cohort's within-population elevation (Korean PTC vs Korean baseline; Chinese GD vs Chinese baseline), not direct cross-population frequency comparison.

**Conclusion**: v1 forest meta and accompanying "Pan-Asian autoimmune-thyroid continuum" framing must be retracted from Paper 2. v2 redesign uses within-population comparator only.

---

## 2. v2 comparator definition

**Comparator**: Korean PTC pool (n=874, HT-overlap substrate) versus **AFND South Korea baseline pool** (Korean general population reference).

**Definition properties**:
- **Same ancestry** — both Korean populations
- **Disease-substrate-matched** — Korean PTC pool with HT-overlap (~20-30%) compared to general Korean population (mixed thyroid health states including normal thyroid). Elevation = enrichment for HT-overlap-associated alleles within Korean PTC.
- **AFND** = Allele Frequency Net Database (`http://www.allelefrequencies.net`), public reference. South Korea baseline = pooled Korean reference population entries (KOTRY donors / KBP / Lee 2014 / cookHLA Korean reference panel).
- **Single-population scope** — no Pan-Asian / Chinese / Japanese cross-population comparison in this spec. Cross-population framing is Paper 4 territory (GD reference background only, 1-line in Paper 2 Discussion).

**Comparator schema**:
```
Korean PTC pool (n=874, Paper 2 substrate)
  └─ vs ─
AFND South Korea baseline pool
  ├─ KOTRY donor cohort (Lee 2014 Tissue Antigens, n~5000+)
  ├─ KBP Korean Bone Marrow Donor Program (n~2000+)
  └─ cookHLA Korean reference panel (Cook et al. Nat Commun 2021)
```

Pooled AFND South Korea baseline n estimate: TBD (decision gate §8).

---

## 3. Cohort composition (Korean PTC pool n=874)

| Sub-cohort | n typeable | Source | HT comorbidity rate |
|---|---|---|---|
| K2 / PRJEB11591 (Yoo 2016) | **235** (subset of n=260 with valid HLA imputation) | ENA PRJEB11591, arcasHLA from RNA-seq | ~20-30% Korean PTC general (per published Korean PTC literature) |
| Lee 2024 / GSE213647 | **630** | GEO GSE213647, arcasHLA from RNA-seq | Hashimoto-like signature 22-28% (Otsu/GMM, R5-2 generalization) |
| GSE286332 PTC arm | **9** (PTC only, NOT PTC+HT) | GEO GSE286332, arcasHLA from RNA-seq | 0% (PTC-only arm by selection) |
| **Total** | **874** | — | — |

**Note**:
- GSE286332 PTC+HT arm (n=9) is excluded from Pillar I HLA pool — those go into Paper 2 Pillar II/III (PTC vs PTC+HT differential mechanism, mediation analysis), not Pillar I HLA susceptibility.
- K2 typeable n=235 reflects arcasHLA confidence cutoff; n=25 dropped due to insufficient read depth on HLA genes (per `v17_arcasHLA_korean_k2.md` memory).
- All Korean PTC HLA pool members have RNA-seq-imputed HLA alleles, NOT SNP-genotyped (cookHLA SNP imputation deferred per memory).

---

## 4. Alleles

**6 alleles in v2 forest** (down-scoped from v1's broader Pan-Asian list):

| # | Allele | East Asian relevance | Lee 2014 baseline status |
|---|---|---|---|
| 1 | **DPB1\*05:01** | East Asian thyroid autoimmunity primary risk allele (HT + GD shared); Korean PTC pool 53.2% per memory | ✅ available (~38-42%) |
| 2 | **A\*02:07** | East Asian-specific HLA-A allele | ✅ available |
| 3 | **B\*46:01** | East Asian thyroid autoimmunity allele | ✅ available |
| 4 | **DRB1\*07:01** | Candidate protective allele | ✅ available |
| 5 | C\*01:02 | East Asian HLA-C allele | ⏸ **PENDING Lee 2014 fill decision** (§8 decision gate) |
| 6 | DQB1\*02:01 | Class II HLA-DQ allele | ⏸ **PENDING Lee 2014 fill decision** (§8 decision gate) |

**Allele exclusion rule**: GD-specific alleles (TSH receptor antibody-associated) NOT included in Paper 2 — they belong to Paper 4 (Korean GD HLA / Pan-Asian).

---

## 5. Statistical plan

### 5.1 Per-allele frequency comparison

For each of the 6 alleles (or 4 if Lee 2014 fill not done):

- **Korean PTC pool**: count carriers / n=874 × 100% = frequency
- **AFND South Korea baseline pool**: count carriers / n_baseline × 100% = baseline frequency
- **Wilson 95% CI** per frequency (preferred over normal approximation for low/high frequencies)

### 5.2 Odds ratio

Per allele:
- 2×2 contingency: [carriers in PTC, non-carriers in PTC; carriers in baseline, non-carriers in baseline]
- **OR with 95% CI** (Wald or exact)
- **Fisher exact test** for two-sided p-value

### 5.3 Multiple testing correction

- **Benjamini-Hochberg FDR** at q=0.05 across the 6 (or 4) alleles
- Rationale: small allele panel, BH appropriate (Bonferroni would be overconservative)

### 5.4 Heterogeneity note

If AFND South Korea baseline pool aggregates multiple sub-populations (KOTRY + KBP + cookHLA + Lee 2014 KOTRY), report:
- Per sub-population baseline frequency separately
- Pooled baseline used for primary OR
- Cochran I² across sub-populations as heterogeneity metric

### 5.5 Sensitivity analysis

- **Sensitivity 1**: Excluding any "Harbin Korean fallback" entries from AFND (if such ethnic-Korean-from-Manchuria entries exist in AFND South Korea pool — verify in §8 decision gate)
- **Sensitivity 2**: K2 + Lee + GSE286332 separately (per-sub-cohort OR) to confirm no single sub-cohort drives the signal
- **Sensitivity 3**: Excluding GSE286332 PTC arm n=9 (small contribution) → robust to small-cohort dropout

### 5.6 Forest plot

- **Per-allele forest panel**: 6 (or 4) horizontal forest entries, each showing OR + 95% CI of (Korean PTC vs AFND South Korea baseline)
- Reference line at OR=1.0
- **NO Chu 2018 GD entries in this forest** (those move to Paper 4 spec)

---

## 6. Output schema

All outputs land under `project/results/paper2_pillar1_v2/` (new directory, created at execution).

| Output | Format | Schema |
|---|---|---|
| `korean_ptc_allele_freq.tsv` | TSV | `allele, n_carriers, n_typeable, frequency, wilson_lower, wilson_upper, sub_cohort` (long format with sub_cohort column for K2/Lee/GSE286332/pooled) |
| `korean_baseline_allele_freq.tsv` | TSV | `allele, n_carriers, n_total, frequency, ci_lower, ci_upper, source` (source = KOTRY/KBP/cookHLA/Lee2014/AFND_pooled) |
| `korean_PTC_vs_korean_baseline_forest.tsv` | TSV | `allele, ptc_freq, baseline_freq, OR, OR_lower_95, OR_upper_95, fisher_p, BH_FDR_q, I2_subpop` |
| `korean_PTC_vs_korean_baseline_forest.json` | JSON | metadata: `n_PTC, n_baseline, alleles_tested, multiple_testing, sensitivity_runs, generation_date, comparator_definition` |
| `forest_paper2_HT_only.pdf` | PDF | Forest plot, 6 alleles, Cell Press style, 1-page |
| `forest_paper2_HT_only.png` | PNG | 300 dpi raster version |
| `discussion_paragraph_v2.md` | MD | **Spec only — paragraph TEMPLATE for Paper 2 Methods + Results sections, with placeholder for actual numbers**. NOT prose generation — template skeleton with `{n_PTC}`, `{OR}`, `{p}` placeholders for fill at execution. Voice-protected prose (Paper 2 mechanism story interpretation) NOT included; only quantitative description. |

---

## 7. Interpretation guard

The following framings are **PROHIBITED** in v2 outputs and any downstream Paper 2 prose:

| Prohibited framing | Why |
|---|---|
| "Korean PTC HLA pattern shares Graves' disease susceptibility" | GD direct equivalence; conflation 1 (deprecated v1) |
| "Pan-Asian autoimmune-thyroid continuum" | Cross-disease (HT + GD) lumping |
| "Korean PTC > Chinese GD direct comparison" | Cross-disease + cross-ancestry conflation |
| "Korean PTC HLA = GD-like" | Disease equivalence claim |
| GD-specific allele effect citations as Paper 2 mechanism | Paper 4 territory |

**ALLOWED framings** (within Paper 2 v2):

| Allowed framing | Notes |
|---|---|
| "Korean PTC pool DPB1\*05:01 frequency exceeds Korean general population baseline (Lee 2014)" | Within-population elevation |
| "Korean PTC HLA enrichment is consistent with Hashimoto-overlap comorbidity (~20-30% in Korean PTC)" | HT-substrate framing |
| "DPB1\*05:01 is a known East-Asian thyroid autoimmunity susceptibility allele (cite Chu 2018 for GD context, 1-line only)" | Background reference, not mechanism |
| "Forest meta of 6 alleles (Korean PTC vs Korean baseline) shows X alleles direction-consistent with HT-associated risk" | Within-Korean comparison |

**1-line GD reference allowance** in Paper 2 Discussion:
> "DPB1\*05:01 is a published East-Asian thyroid autoimmunity allele (per Chu et al. 2018 *J Med Genet* in Han Chinese Graves' disease; OR=1.90 [1.69, 2.14]). Mechanistic interpretation in Paper 2 is restricted to Hashimoto's thyroiditis context."

This 1-line is the maximum GD reference; any further GD discussion belongs in Paper 4 backlog.

---

## 8. Decision gates (advisor + user before execution)

| Gate | Question | Required by |
|---|---|---|
| **G1** | Advisor (Yu professor) approval of v2 spec? | Before execution |
| **G2** | AFND South Korea baseline pool composition: KOTRY only / KOTRY+KBP / KOTRY+KBP+cookHLA / + Lee 2014? | Before allele frequency lookup |
| **G3** | Lee 2014 fill decision: include C\*01:02 + DQB1\*02:01 (6-allele forest) or skip (4-allele forest)? | Before forest plot generation |
| **G4** | "Harbin Korean fallback" entries in AFND South Korea: verify exist + decision to exclude in sensitivity? | Before sensitivity §5.5 |
| **G5** | Figure 1 (manuscript) re-render decision: replace v1 forest plot in Paper 2 Figure 1 (or Pillar I figure) with v2? Affects Paper 2 outline + figure caption v2. | After v2 execution complete |
| **G6** | Paper 2 manuscript paragraph (Methods + Results) v2 framing approval before Paper 2 prose draft? | After v2 execution + advisor review |

**Default decision pending user/advisor input**: All gates open until explicit confirmation. Spec frozen as of 2026-05-04, execution plan not started.

---

# END SPEC — execution requires user + advisor sign-off on G1-G6.
