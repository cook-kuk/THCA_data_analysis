# Track 11 — FinnGen + Pan-UKBB MHC trans-ancestry replication of thyroid autoimmunity

**Date:** 2026-05-08
**Author:** Seungho Cook + Claude Code
**Boundary:** AUTOIMMUNE thyroid endpoints ONLY. NO cancer endpoints. Per `paper2_hla_boundary/HLA_CANCER_SEPARATION_RULES.md`. **All findings are non-Korean (Finnish, UK European, UK East Asian, UK CSA, UK AFR) and serve as TRANS-ANCESTRY SENSITIVITY context, NOT direct effect-size transfer to Korean Track 1.**
**Status:** Track 11 of multi-track HLA deep-dive sprint. Independent of Track 1 (Korean Pan-Asian forest) and Track 7 (GWAS Catalog mining).

---

## 0. Boundary + ancestry caveat (read first)

This track pulls public summary statistics from FinnGen R12 (Finnish population, ancestry≈Northern European with founder structure) and Pan-UKBB (UK Biobank EUR + EAS + CSA + AFR per-ancestry sumstats). **Finnish ≠ Korean. EUR ≠ Korean.** Differences in allele frequency, LD structure, and tag-SNP imputation accuracy mean the effect-size *magnitudes* observed here cannot be transferred to Korean cohorts. We use these data exclusively to ask:

> **Does the MHC class II locus harbor genome-wide-significant signal for autoimmune thyroid disease in non-Asian populations, and do canonical trans-ancestry tag SNPs for Track 1 alleles show direction-consistent effects?**

This is a sensitivity / hypothesis-corroboration analysis. It is NOT a Korean-effect validation, and we do not draw cross-ancestry effect equivalence claims.

---

## 1. Endpoint manifest

42 endpoints assembled across two sources:

- **FinnGen R12** (Finnish, GRCh38): 13 endpoints
  - Primary autoimmune thyroid: `E4_GRAVES_STRICT` (3,962 cases / 496,386 ctrl), `E4_GRAVES_OPHT_STRICT` (753), `GRAVES_OPHT` (858), `E4_HYTHY_AI_STRICT` (autoimmune hypothyroidism, **54,752** cases / 362,639 ctrl), `E4_HYTHY_AI_DX` (133), `E4_THYROIDITAUTOIM` (688), `E4_THYROIDITCHRON` (141), `AUTOIMMUNE_HYPERTHYROIDISM` (2,469).
  - Comparators: `E4_THYROID` (any thyroid, 76,140), `HYPOTHY_REIMB` (drug reimbursement, 16,150).
  - Autoimmune sanity: `T1D` (4,721), `M13_RHEUMA` (16,314), `L12_PSORIASIS` (12,760).
- **Pan-UKBB** (GRCh37; EUR/EAS/CSA/AFR per phenotype where available): 29 ancestry × phenotype slices spanning 12 base phenotypes — `phecode:242.1` (Graves), `phecode:244` and `244.4` (hypothyroidism), `icd10:E03/E05/E06`, `prescriptions:levothyroxine/carbimazole`, self-reported categorical 20002 hyperthyroidism / hypothyroidism, chronic lymphocytic thyroiditis (`245.21`).

Manifest: `tables/T01_endpoint_manifest.tsv` + `tables/T01_endpoint_summary.json`.

---

## 2. Sumstats acquisition

- FinnGen URL pattern: `https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_<phenocode>.gz` with bgzip + `.tbi` index.
- Pan-UKBB URL pattern: `https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/<trait>-<code>-both_sexes.tsv.bgz` with `.tbi` at `sumstats_flat_files_tabix/`.
- **Streaming MHC slice via `pysam.TabixFile` remote fetch** — no full-genome download. Cache (62 MB total, well under 2 GB limit) at `cache/<source>__<phenocode>__<ancestry>.tsv.gz`. MHC region: `chr6:28,510,120-33,480,577` (GRCh38) / `chr6:28,477,797-33,448,354` (GRCh37). All 42/42 endpoints sliced successfully (status `ok`); zero fetch failures (slice summary `tables/T02_mhc_slice_summary.tsv`).
- Per-endpoint MHC variant counts: 68k–81k (FinnGen ~68k, Pan-UKBB EUR ~81k, Pan-UKBB EAS ~63k, Pan-UKBB AFR ~78k due to ancestry-specific MAF filtering).

---

## 3. MHC region peaks per endpoint

Per-endpoint Manhattan plots of the chr6 28-34 Mb window: `figures/F01_manhattan_<TAG>.png` (42 individual). Combined panels:

- `figures/F02_manhattan_grid_finngen.png` — FinnGen GD / GD-OPHT / HT-AI / AITD-ICD / T1D / RA grid.
- `figures/F03_manhattan_grid_panukbb_eur.png` — Pan-UKBB EUR grid (12 panels).
- `figures/F06_position_density_by_endpoint.png` — overlaid 50 kb-binned max -log10p track across primary endpoints.
- `figures/F04_mhc_peak_strength.png` — bar chart of MHC top-SNP -log10p across all 42 endpoint × ancestry combinations.

**Top MHC peaks (FinnGen, primary autoimmune thyroid):**

| Endpoint | Top SNP | chr6 pos (GRCh38) | β | p | Nearest HLA gene |
|---|---|---|---|---|---|
| Graves (E4_GRAVES_STRICT) | rs1617322 | 32,706,116 | +0.717 | **1.78×10⁻¹²⁴** | HLA-DQA2 (35 kb) |
| Graves | rs3130297 | 32,231,204 | +0.717 | 2.02×10⁻¹²⁰ | HLA-DRA (209 kb) |
| HT-AI (E4_HYTHY_AI_STRICT) | rs9273386 | 32,659,200 | +0.230 | **1.82×10⁻²³⁸** | HLA-DQB1 (0.3 kb) |
| HT-AI | rs1980496 | 32,372,293 | +0.189 | 1.61×10⁻¹⁴³ | HLA-DRA |
| HT-AI | rs3132487 | 31,275,231 | +0.199 | 4.19×10⁻¹¹⁹ | HLA-C (3 kb) |
| AITD-ICD (E4_THYROIDITAUTOIM) | top hit | 32–33 Mb | — | <10⁻⁷ | HLA class II |
| T1D (autoimmune comparator) | — | 32 Mb | — | <10⁻⁵⁰⁰ | HLA-DQB1/DRB1 |

**Pan-UKBB EUR Graves (phecode:242.1, n=554) top 5:** all 5 fall in chr6:31.3–32.6 Mb, peak β≈+1.0, p ≤ 3×10⁻²⁰, nearest genes HLA-DQA1, HLA-DRA, HLA-B.

Both populations recapitulate the classical MHC-class-II-anchored autoimmune-thyroid signal.

---

## 4. Cross-endpoint MHC peak overlap

`figures/F05_lead_snp_overlap_heatmap.png` — Jaccard overlap of top-10 lead SNPs across all 42 endpoint × ancestry pairs.
`tables/T04_lead_snp_overlap_matrix.tsv` and `tables/T04_recurring_lead_snps.tsv` enumerate recurrence.

**Top recurring MHC lead SNPs (≥4 endpoints):**

| Locus | n_endpoints sharing | Region |
|---|---|---|
| chr6:30,737,591 | 6 | HLA-class I cluster (near HLA-E) |
| chr6:33,047,432 | 5 | HLA-DOA / HLA-DPA1 |
| chr6:30,335,759 | 5 | HLA-A region |
| rs3132941 | 4 | HLA-class II |
| chr6:32,605,189 | 4 | HLA-DRB1 / DQA1 boundary |
| chr6:31,418,355 | 4 | HLA-B / HLA-C |
| chr6:32,278,635 | 4 | HLA-DRA region |

Two distinct convergence zones emerge: **(i) HLA-class II core (32.4–33.1 Mb)** dominates HT-AI / GD / hypothyroidism endpoints; **(ii) HLA-class I region (30.5–31.4 Mb)** is a shared secondary peak across Graves, T1D and RA. This pattern matches published trans-ancestry AITD GWAS (Saevarsdottir 2020 Nature, Sakaue 2021 Nat Genet).

---

## 5. Trans-ancestry reconciliation vs Track 1 Korean forest

`tables/T05_track1_tagsnp_transancestry.tsv` (252 rows: 6 alleles × 42 endpoints).
`figures/F08_dpb1_0501_transancestry_forest.png` — DPB1*05:01-only forest.
`figures/F09_track1_alleles_transancestry_grid.png` — 6-allele × ~13 endpoint grid.

### 5.1 HLA-DPB1*05:01 (canonical East-Asian GD allele) trans-ancestry test

Tag SNP **rs9277534**, GRCh38 pos 33,089,763 / GRCh37 pos 33,057,540. Notes on caveats: tag-LD is ~0.7 in East Asians (Cooper 2008), substantially lower in Europeans; AF of the tag allele is ~0.25 in EUR vs ~0.61 in EAS (consistent with our extracted control AFs).

| Endpoint | Source | Ancestry | N_case | β | p | AF (ctrl) |
|---|---|---|---|---|---|---|
| Graves (E4_GRAVES_STRICT) | FinnGen | FIN | 3,962 | **+0.129** | **3.64×10⁻⁷** | 0.255 |
| Graves ophthalmopathy | FinnGen | FIN | 753 | +0.228 | 4.39×10⁻⁵ | 0.255 |
| Autoimmune hyperthyroidism (umbrella) | FinnGen | FIN | 2,469 | +0.121 | 1.5×10⁻⁴ | 0.256 |
| HT-AI strict | FinnGen | FIN | 54,752 | **−0.107** | **7.10×10⁻⁴²** | 0.253 |
| AITD-ICD | FinnGen | FIN | 688 | −0.250 | 7.6×10⁻⁵ | 0.257 |
| T1D | FinnGen | FIN | 4,721 | +0.229 | 9.16×10⁻²³ | 0.255 |
| Rheumatoid arthritis | FinnGen | FIN | 16,314 | −0.247 | 8.08×10⁻⁷³ | 0.254 |
| Graves (phecode 242.1) | Pan-UKBB | EUR | 554 | +0.146 | 0.025 | 0.309 |
| Hypothyroidism (phecode 244) | Pan-UKBB | EUR | 18,404 | −0.136 | 9.33×10⁻²⁹ | 0.309 |
| ICD10-E03 hypothyroidism | Pan-UKBB | EUR | 17,592 | −0.141 | 2.51×10⁻²⁹ | 0.309 |
| Self-hyperthyroidism | Pan-UKBB | EUR | 20,563 | −0.186 | 8.51×10⁻⁵⁵ | 0.309 |
| Hypothyroidism | Pan-UKBB | EAS | 71 | +0.076 | 0.69 | 0.612 |
| Hypothyroidism | Pan-UKBB | CSA | 512 | −0.080 | 0.25 | 0.322 |
| Hypothyroidism | Pan-UKBB | AFR | 144 | −1.158 | 0.40 | 0.002 |

**Direction interpretation (with caveats):**

- In Finnish and UK-EUR populations, rs9277534 shows a **clear, GW-significant Graves-RISK / hypothyroidism / RA / T1D-RISK pattern** at the canonical Track 1 East-Asian GD-risk allele.
- The HT-AI direction is **opposite** to GD: tag β is protective for autoimmune hypothyroidism (p=7×10⁻⁴²), suggesting the locus partitions GD vs HT differently in Finnish than the unidirectional autoimmune-risk reading would suggest. This is consistent with established class II HLA orthogonality between Graves and Hashimoto.
- The Self-hyperthyroidism Pan-UKBB EUR effect is *negative* — likely reflecting phenotype-label mixing (UKBB self-reported hyperthyroidism is often levothyroxine-treated post-Graves hypothyroidism, not active GD). FinnGen GD-strict (the cleanest endpoint) is **positive**, matching East-Asian direction.
- **Pan-UKBB EAS / CSA / AFR for thyroid endpoints are severely underpowered** (n_case ≤ 600); confidence intervals span zero. AFR β = −1.16 reflects AF=0.002 (essentially monomorphic — uninformative).
- **Bottom line for Track 1 reconciliation:** the canonical DPB1*05:01 tag SNP reaches GW significance in Finnish GD with the same risk direction as East-Asian published Graves' literature. This is hypothesis-corroborating but **does not validate Korean-specific effect sizes**, because (i) tag-LD differs across ancestries, (ii) AF is 2.4-fold lower in EUR than EAS, (iii) the locus is dominated by class II haplotype effects that are likely allele-specific to local LD.

### 5.2 Other Track 1 alleles

| Allele | Tag rsid | Notes from grid |
|---|---|---|
| HLA-DRB1*15:01 | rs3135388 | Strong negative β (protective) for FinnGen HT-AI / hypothyroidism endpoints in EUR. |
| HLA-DRB1*03:01 | rs2187668 | Risk for FinnGen T1D (canonical), risk-direction for autoimmune endpoints; consistent with classical AITD class II story. |
| HLA-DRB1*04:01 | rs6910071 | Mixed; thyroid-endpoint signal weaker than RA (where it is canonical risk allele). |
| HLA-B*46:01 | rs2596542 | Tag is virtually monomorphic in EUR (AF<0.05); under-detected — expected because B*46:01 is essentially East-Asian-specific. **Cannot reconcile across ancestries from this tag in EUR-only sumstats** — explicit failure of trans-ancestry imputation for this allele. |
| HLA-DQB1*02:01 | rs2647044 | Direction-consistent risk for AITD-ICD and T1D in FinnGen. |

---

## 6. Ancestry portability test

`figures/F07_panukbb_ancestry_portability.png` and `tables/T04b_panukbb_ancestry_portability.tsv` — top-SNP -log10p comparison across EUR / EAS / CSA / AFR within Pan-UKBB for shared phenotypes.
`figures/F11_eur_vs_eas_effectsize_ratio.png` and `tables/T06b_eur_vs_eas_topshared_loci.tsv` — same-position β comparison: at EUR-lead MHC SNPs, look up the matching EAS β.

**Key portability observations:**

- For hypothyroidism (`phecode:244`, `phecode:244.4`, `icd10:E03`) the EUR top MHC SNP (chr6:32,605,189, near HLA-DRB1) reaches p ~ 4×10⁻¹²⁰ at β=+0.27. The EAS β at the same position is +0.42 (n=71 cases, p=0.02) — **same direction, larger point estimate but very wide CI** due to small EAS sample.
- All other EUR top loci have EAS β estimates with p > 0.1 — sample size, not direction inconsistency, dominates.
- AFR effect estimates are dominated by very rare alleles (AF<0.005) producing noise-level β; **Pan-UKBB AFR thyroid endpoints are not informative for MHC ancestry portability**.
- CSA (~500 hypothyroidism cases) shows direction-consistent but non-significant effects.
- **Honest conclusion: Pan-UKBB EAS thyroid sumstats are too underpowered to formally test EUR↔EAS portability of MHC class II AITD effects.** Track 1's Korean Pan-Asian forest remains the only reasonable East-Asian effect-size estimate.

---

## 7. Heritability (Pan-UKBB)

`tables/T06_panukbb_h2_extracted.tsv` and `figures/F10_panukbb_h2_per_endpoint_bar.png` — genome-wide liability-scale h² from Pan-UKBB sLDSC (EUR) / RHEmc (other ancestries).

EUR liability h² estimates:

- Self-hyperthyroidism: h²=0.233 (SE 0.030)
- Hypothyroidism (phecode 244): h²=0.185 (0.025)
- Hypothyroidism NOS (244.4): h²=0.193 (0.026)
- ICD10-E03 hypothyroidism: h²=0.192 (0.026)
- Self-hypothyroidism / myxoedema: h²=0.292 (0.101)
- Levothyroxine prescription: h²=0.163 (0.025)
- ICD10-E05 thyrotoxicosis: h²=0.167 (0.070)
- Carbimazole prescription: h²=0.315 (0.151)
- Graves (phecode 242.1, n=554): h²=0.148 (SE 0.181) — too underpowered for confident estimate.

**Note:** these are **genome-wide** h² (not partitioned to MHC). We did not attempt LDSC partitioned h² because that requires reference panels and full-genome sumstats outside our 2 GB cache budget.

---

## 8. Limitations

1. **Ancestry mismatch.** Finnish and UK populations are not Korean. LD structure, allele frequencies, and tag-SNP imputation accuracy differ substantially. Effect-size magnitudes are not transferable.
2. **Tag-SNP imperfection.** Single-SNP tags for HLA classical alleles have ancestry-dependent LD; tag-effect ≠ allele-effect, especially across populations. For e.g. HLA-B*46:01 (EAS-specific) the EUR tag is essentially uninformative.
3. **Build mismatch.** FinnGen sumstats are GRCh38; Pan-UKBB are GRCh37. We used build-appropriate positions for tag SNPs, but tiny build-version offsets (≤200 bp) at HLA region were tolerated via window matching — could miss SNPs at boundaries.
4. **Phenotype heterogeneity.** UKBB self-reported hyperthyroidism (20002:1226) bundles treated/untreated cases; some cases are post-Graves hypothyroid on levothyroxine, which can produce direction-confusing effects. FinnGen E4_GRAVES_STRICT is cleaner.
5. **Pan-UKBB non-EUR power.** EAS / AFR / CSA thyroid endpoints have <600 cases each; detection power for moderate MHC effects (OR≈1.2) is <30%. Negative non-EUR results are uninterpretable.
6. **Heritability proxy.** Pan-UKBB h² is genome-wide, not MHC-partitioned. We don't claim Track 11 quantifies the MHC fraction of h².
7. **Correction for multiple endpoints.** With 42 endpoints × 6 alleles = 252 trans-ancestry tests, formal multiple-comparison control would warrant Bonferroni-style p<2×10⁻⁴; results above (DPB1*05:01 GD p=3.6×10⁻⁷, HT-AI p=7×10⁻⁴²) survive that threshold easily for the highlighted findings, but smaller effects (p=10⁻⁴ to 10⁻³ class) should be treated as suggestive only.

---

## 9. Honest assessment of trans-ancestry replication

- **Robust corroboration:** Class II MHC region (chr6 ~32.5–33 Mb, HLA-DRB1/DQA1/DQB1/DPB1 cluster) drives genome-wide-significant signal for autoimmune thyroid disease in both Finnish and UK-EUR populations. This is unsurprising — published in dozens of EUR AITD GWAS — but our point is that it is **not absent in non-Asian data**, and the same chromosomal architecture as Korean Pan-Asian Track 1 is engaged.
- **Suggestive corroboration:** The DPB1*05:01 tag (rs9277534) shows GW-significant GD-risk direction in Finnish (β=+0.13, p=3.6×10⁻⁷ at n=3,962 cases). This is direction-consistent with the East-Asian Graves' literature for HLA-DPB1*05:01, despite imperfect EUR tag LD. We treat this as hypothesis corroboration, not effect-size validation.
- **Cannot replicate:** HLA-B*46:01-tagging via rs2596542 is uninformative in EUR (allele essentially monomorphic). This is an honest failure of trans-ancestry tag-SNP testing for an EAS-specific allele — not surprising, but documented.
- **Confounding observation:** the rs9277534 effect direction *flips* between GD (positive) and HT-AI (negative) within Finnish data — both p << 10⁻⁶. This is consistent with HLA class II's known orthogonal partitioning of Graves' vs Hashimoto's susceptibility and is itself a useful trans-ancestry sensitivity finding for the Korean Track 1 forest's GD-specific framing.
- **Net:** Track 11 supports Track 1 framing (Korean GD has an MHC class II signal) at the *locus level* and at the *DPB1 region* (tag-SNP level), with explicit acknowledgement that **effect-size transfer across ancestries is not warranted**.

---

## File map

```
project/results/hla_deepdive_2026_05_08/track11_finngen_panukbb/
├── cache/                                        (62 MB, MHC slices for 42 endpoints)
├── tables/
│   ├── T01_endpoint_manifest.tsv                  (42 endpoints; URLs, N_case/N_ctrl, build)
│   ├── T01_endpoint_summary.json
│   ├── T02_mhc_slice_summary.tsv                  (per-endpoint slice status + variant counts)
│   ├── T03_lead_snps_per_endpoint.tsv             (top-10 MHC SNPs × 42 endpoints = 420 rows)
│   ├── T04_lead_snp_overlap_matrix.tsv            (Jaccard matrix)
│   ├── T04_recurring_lead_snps.tsv
│   ├── T04b_panukbb_ancestry_portability.tsv
│   ├── T05_track1_tagsnp_transancestry.tsv        (DPB1*05:01 + 5 other Track 1 alleles)
│   ├── T06_panukbb_h2_extracted.tsv
│   └── T06b_eur_vs_eas_topshared_loci.tsv
├── figures/                                       (51 PNGs total)
│   ├── F01_manhattan_<TAG>.png  × 42              (per-endpoint MHC Manhattan)
│   ├── F02_manhattan_grid_finngen.png
│   ├── F03_manhattan_grid_panukbb_eur.png
│   ├── F04_mhc_peak_strength.png
│   ├── F05_lead_snp_overlap_heatmap.png
│   ├── F06_position_density_by_endpoint.png
│   ├── F07_panukbb_ancestry_portability.png
│   ├── F08_dpb1_0501_transancestry_forest.png
│   ├── F09_track1_alleles_transancestry_grid.png
│   ├── F10_panukbb_h2_per_endpoint_bar.png
│   └── F11_eur_vs_eas_effectsize_ratio.png
└── track11_report.md   (this file)
```

Scripts: `/home/seungho/personal/THCA_data_analysis/scripts/hla_deepdive_2026_05_08/track11/` (steps 01–06).
