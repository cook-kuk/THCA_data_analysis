# Track 3 — Pan-Asian Healthy HLA Atlas

**Date:** 2026-05-08
**Owner:** Seungho Cook (Track 3 of 6 in the HLA deep-dive sprint)
**Scope:** Healthy-only resource paper (tools / data-paper framing)
**Boundary contract:** `project/paper2_hla_boundary/HLA_CANCER_SEPARATION_RULES.md`

---

## 0. Hard boundary (healthy only — no disease association)

This deliverable is a **healthy-controls HLA frequency atlas**. No PTC, GD,
Hashimoto, or any other disease association is tested or implied. Every
cohort included contains only normal/healthy thyroid tissue or healthy
non-disease control populations. The atlas is a description of population HLA
allele frequencies and their cross-population structure — closer to a
"resource / methods" deliverable than to a hypothesis-testing paper.

The boundary doc forbids HLA × cancer-outcome joins; this Track 3 work does
not do any such join. It also does not import any AITD susceptibility claim.
Any downstream user of this resource who joins these tables to cancer
outcomes must satisfy their own boundary contract.

---

## 1. Data sources & filtering

| population | country | n_individuals | typing | source |
|---|---|---:|---|---|
| Baek2021_KoreanNGS | South Korea | 173 | NGS long-range PCR (TruSight HLA v2) | Baek IC et al. PLoS ONE 2021;16:e0253619 (S9 + S10 tables, parsed) |
| K2_normal_RNAseq | South Korea | 81 | RNA-seq imputed (arcasHLA) | PRJEB11591 (Yoo SK 2016) — adjacent-normal subset (`-N` suffix) |
| GSE213647_normal_RNAseq | South Korea | 263 | RNA-seq imputed (arcasHLA) | Lee SE 2024 Nat Commun (PMID 38331894) — `tissue_type == Normal` subset |
| AFND_Korea | South Korea | ≈82,197 (164,394 chr) | Sanger / PCR-SBT, weighted pool | Allele Frequency Net Database (AFND) — Korea + Korean diaspora pops |
| AFND_Japan | Japan | ≈20,110 | Sanger / PCR-SBT, weighted pool | AFND Japan (Ainu excluded) |
| AFND_HanChinese | China (Han) | ≈8,823 | Sanger / PCR-SBT, weighted pool | AFND Beijing/Hubei/Jiangsu/Canton/Shanghai/Yunnan-Han |
| AFND_Taiwan_Han | Taiwan (Han) | ≈1,371 | Sanger / PCR-SBT, weighted pool | AFND Taiwan Han + Hakka + Minnan + Tzu Chi Cord Blood Bank |
| AFND_Taiwan_Indigenous | Taiwan (Indigenous) | ≈212 | Sanger / PCR-SBT, weighted pool | AFND Atayal / Ami / Bunun / Pazeh / Siraya |
| AFND_EastAsianAggregate | East Asian aggregate | ≈113,282 | Sanger / PCR-SBT, weighted pool | AFND Korea + Japan + China + Taiwan-Han (excludes Indigenous) |

### 1.1 Excluded sources

* **GSE286332** — Korean PTC vs PTC+HT (n=18, all tumor) contains **zero
  healthy-normal samples**. Per the GEO sample sheet
  (`/data/thca/repo_results/d4p1_panasian_meta/GSE286332_arcasHLA_genotypes.tsv`)
  the 18 samples are split 9 PTC / 9 PTC+HT with no normal arm. Including
  GSE286332 would import tumor-derived HLA calls into a healthy atlas, a
  methodological boundary violation. It is therefore **dropped from this atlas**
  with the explicit note that the underlying pre-cancer Korean RNA-seq pool of
  this resource is already represented by GSE213647 normals (n=263) and K2
  adjacent-normals (n=81).

### 1.2 Filtering rules

* **K2 (PRJEB11591):** 261 thyroid samples; 81 with `-N` suffix on `sample_alias`
  (=adjacent-normal thyroid; SNU-GMI-PTC-N, FTC-N, FA-N, FV-N). All 81 had
  successful arcasHLA `genotype.json` files.
* **GSE213647:** sample sheet `tissue_type` field, `Normal` only — 263 of 632
  samples. arcasHLA result table uses run accession (SRR…) which we map back to
  the Lee clinical sheet by sequential GSM ordering when needed.
* **AFND:** sub-populations pooled with 2N-weighted means within each country
  bucket. Sample-size threshold not enforced because AFND already only retains
  large reference cohorts at the subnational level.

### 1.3 Important caveat on AFND coverage

The local AFND cache contains five 4-digit alleles per population
(`B*46:01`, `DPB1*05:01`, `DQB1*06:02`, `DRB1*04:05`, `DRB1*15:01`), the set
that was originally fetched for autoimmune-thyroid relevance. The Baek2021 NGS
table provides 146 4-digit alleles in Koreans, and the RNA-seq cohorts each
type ≈150 alleles per cohort. Therefore:

* The **master matrix** is rich on the Korean side (RNA-seq + Baek2021) but
  sparse on the AFND-Japan / AFND-Han / AFND-Taiwan side.
* For honest cross-population comparisons we provide a **dual-PCA**:
  (i) primary PCA on all 145 alleles (data-coverage-confounded), and
  (ii) bias-corrected PCA on the 5 AFND-shared alleles only.
* Pairwise Korean-vs-comparator deltas (`korean_vs_pairwise_signature_alleles.tsv`)
  use the inner-join on shared alleles only.

A future expansion of the AFND cache to all 4-digit class-I/II alleles would
unlock a full pan-Asian PCA. This is documented under `9. Future cohorts`.

---

## 2. Master matrix

* `master_population_x_allele_freq.tsv` — 145 alleles × 9 populations
  (alleles seen in ≥2 populations; sorted by total frequency; long-tail
  alleles seen in only one population are dropped to reduce sparsity).
* `master_population_metadata.tsv` — n, typing method, source per
  population.
* `per_locus_top_alleles_long.tsv` — long-format frequency table including
  all populations and loci.
* `Baek2021_korean_ngs_2field.tsv` — Baek 2021 NGS 2-field allele freqs.
* `K2_normal_genotypes_2field.tsv`,
  `GSE213647_normal_genotypes_2field.tsv` —
  per-sample 2-field genotype calls (raw arcasHLA results, healthy subset only).
* `afnd_subpopulation_long.tsv` — every AFND sub-population row with its
  weight before pooling.

---

## 3. PCA / UMAP

### 3.1 Primary PCA (`F1_pca_populations.png`, `pca_coords.tsv`)

PC1 (76.8 % var) splits the **Korean-typed populations** (Baek2021, K2, Lee
RNA-seq) from the **AFND-pooled populations** along the data-coverage axis.
This is largely a *technical* axis — AFND populations are zero-coverage on
the 140 alleles unique to the rich Korean datasets, so they cluster at low
PC1. PC2 (14.5 %) separates Taiwan Indigenous and Taiwan Han from the rest of
East Asia, which is biology consistent with prior AFND reports (Austronesian
substructure on the island).

### 3.2 Bias-corrected PCA (`F1_pca_populations_afnd_shared.png`)

Restricted to the 5 AFND-shared alleles. All Korean populations (AFND_Korea,
Baek2021, K2_normal, GSE213647_normal) cluster in a tight central region.
**Taiwan_Indigenous lies at the extreme of PC1** (driven by `DPB1*05:01`
frequency ≈0.79, ≈2× higher than mainland East Asia) and Taiwan_Han at the
opposite end of PC2. The Korean / Japanese / Han Chinese populations are
mutually close, consistent with expected mainland East Asian similarity.

### 3.3 Korean-only PCA (`F1_pca_populations_korean_only.png`)

Three Korean platforms (Baek2021 NGS / K2 RNA-seq / Lee RNA-seq) on Korean-
rich allele set. Baek2021 NGS sits between the two RNA-seq cohorts on PC1,
with PC2 cleanly separating PRJEB11591 (FFPE 2014) and Lee 2024 (Fresh
Frozen). PC scores differ but pairwise Pearson r=0.86–0.93 (Section 5.1).

### 3.4 UMAP (`F2_umap_populations.png`)

Top-100 most-variable alleles, n_neighbors=5, min_dist=0.3. UMAP recapitulates
the same population structure (Korean cluster vs Taiwan Indigenous vs AFND
non-Korean) without forcing linearity.

---

## 4. Hierarchical clustering / Jensen–Shannon

### 4.1 Euclidean Ward (`F3_dendrogram_euclidean.png`)

Two macro-clades: (a) the data-rich Korean cohorts (Baek2021 + RNA-seq), and
(b) the AFND-pooled populations. Inside the AFND clade, Korea / Japan / Han /
EA-aggregate cluster together; Taiwan Han + Taiwan Indigenous form a tight
sister clade.

### 4.2 Jensen–Shannon divergence (`F4_dendrogram_jsd.png`,
`jsd_distance_matrix.tsv`)

Within-locus normalized JS divergence averaged across the six loci. Reads
better than Euclidean for compositional data (HLA freq within a locus is
inherently a probability distribution).

| pair | JSD |
|---|---:|
| K2_normal_RNAseq ↔ GSE213647_normal_RNAseq | 0.206 |
| Baek2021_KoreanNGS ↔ GSE213647_normal_RNAseq | 0.215 |
| Baek2021_KoreanNGS ↔ K2_normal_RNAseq | 0.264 |
| AFND_Taiwan_Han ↔ AFND_Taiwan_Indigenous | 0.000 (only DPB1*05:01 drove this; tied) |
| AFND_HanChinese ↔ AFND_EastAsianAggregate | 0.001 |
| AFND_Korea ↔ AFND_Japan | 0.259 |

Three Korean platforms are mutually within 0.21–0.26 JSD — the closest of any
non-trivial cross-platform pair, confirming that NGS-typed Korean and
RNA-seq-imputed Korean atlas blocks describe the same population.

---

## 5. Korean signature alleles

### 5.1 Cross-Korean platform concordance (`F9_cross_korean_concordance.png`,
`cross_korean_concordance.tsv`)

| pair | n_alleles | Pearson r | Spearman ρ |
|---|---:|---:|---:|
| Baek2021_KoreanNGS vs K2_normal_RNAseq | 145 | 0.862 | 0.813 |
| Baek2021_KoreanNGS vs GSE213647_normal_RNAseq | 145 | 0.906 | 0.850 |
| K2_normal_RNAseq vs GSE213647_normal_RNAseq | 145 | 0.928 | 0.857 |

This is a tools-paper-grade concordance: NGS Sanger-equivalent vs arcasHLA
RNA-seq imputation aligns at r≈0.86–0.91 across 145 4-digit alleles in
healthy Koreans. Ties together the broader claim that arcasHLA RNA-seq is a
defensible substrate for population-level HLA frequency estimation when the
two NGS-typed and RNA-seq-typed datasets agree at this level.

### 5.2 Korean vs East Asian sub-populations
(`F6_korean_signature_bars.png`, `korean_vs_pairwise_signature_alleles.tsv`)

Restricted to the 5 AFND-shared alleles (the AFND cache limit). Using
Baek2021_KoreanNGS as the Korean reference:

| comparator | enriched_in_Korean | depleted_in_Korean |
|---|---|---|
| AFND_Japan | `B*46:01` (Δ=+0.013) | `DRB1*04:05` (Δ=−0.034), `DPB1*05:01` (−0.036) |
| AFND_HanChinese | `DRB1*04:05` (Δ=+0.032) | `B*46:01` (Δ=−0.059), `DRB1*15:01` (−0.033) |
| AFND_Taiwan_Han | (none in this 5-allele set) | `DPB1*05:01` (Δ=−0.107), `B*46:01` (−0.071) |
| AFND_Taiwan_Indigenous | `B*46:01` (+0.028) | `DPB1*05:01` (Δ=−0.452, Indigenous near-fixation) |

The Δ ≥ 0.05 absolute threshold is purely descriptive. For richer Korean-
specific alleles (those uniquely captured by NGS+RNA-seq but absent from the
AFND-cached comparator), see Section 9 future expansion.

### 5.3 Korean vs East Asian aggregate
(`korean_vs_eastasian_signature_alleles.tsv`)

Bottlenecked by the same 5-allele AFND aggregate. The single |Δ|≥0.05 hit is
`DQB1*06:02` enriched in K2_normal_RNAseq (Δ=+0.071). This is a single-cohort
signal at the 81-individual scale, not a meta-validated claim.

---

## 6. Locus diversity (`F8_locus_diversity_bars.png`,
`locus_diversity.tsv`)

Shannon entropy and expected heterozygosity per locus per population. The
RNA-seq + NGS Korean cohorts dominate diversity (because they observe the
full long tail of low-frequency alleles); AFND pools have intentionally
truncated diversity because only 5 alleles are cached.

| population | mean Shannon H (bits) across 6 loci |
|---|---:|
| GSE213647_normal_RNAseq | 4.05 |
| K2_normal_RNAseq | 4.01 |
| Baek2021_KoreanNGS | 3.77 |
| AFND_* (any) | 0–0.97 (capped by cache size, not biology) |

Within Korean cohorts, the most diverse locus is HLA-B (Shannon H ≈ 4.6–4.7
bits). The least diverse is DPB1 (≈ 2.7–2.8 bits) — consistent with prior
literature reporting `DPB1*05:01` as a dominant Korean / East Asian allele
at ≈34–45 % frequency.

---

## 7. Limitations

1. **AFND cache size.** Only 5 alleles per AFND population in the cached
   JSON. Pan-Asian PCA on all 145 alleles is therefore data-coverage-biased;
   we mitigate via (i) AFND-shared 5-allele PCA, and (ii) the pairwise Δ
   tables that use inner-join on shared alleles. A future re-pull of AFND
   for all 4-digit class-I and class-II alleles is needed to remove this.
2. **Typing-platform heterogeneity.** RNA-seq imputation (arcasHLA) and NGS
   Sanger typing (Baek2021, AFND) are not equivalent. We document this with
   the `typing_method` column in the metadata table and report cross-platform
   r≈0.86–0.93 within Koreans as an empirical validation. arcasHLA can miss
   alleles at low coverage and tends to call homozygosity in heterozygotes
   when one allele has lower expression.
3. **Sample-size disparity.** AFND populations carry 100–80,000 individuals
   (Sanger-typed pool); the Korean RNA-seq cohorts carry 81–263. CI widths
   are not directly comparable across rows.
4. **Adjacent-normal vs healthy.** K2 normals (n=81) are *adjacent-normal*
   thyroid tissue from cancer-bearing patients. Germline HLA is not affected
   by the adjacent tumor (HLA loss-of-heterozygosity is tumor-only), so the
   germline genotype is valid. Caveat documented.
5. **GSE213647 SRR↔GSM mapping.** The clinical sheet keys on GSM; we use
   sequential ordering of SRR runs onto sorted GSMs when no explicit run
   column is present. This is a heuristic and could mis-assign labels at the
   boundary if the GEO upload order differed from our SRR enumeration. A
   strict per-sample audit would resolve any rare residual mismatch (most of
   the 263 Normal labels can be cross-checked against the sample title
   `*-N`-style suffix).

---

## 8. Utility statement (downstream uses)

This atlas is meant to support, *as a healthy reference only*:

* **Pre-disease baseline frequency tables** for any East-Asian-focused
  HLA-association study that needs to draw "expected frequency under HW
  equilibrium" from a same-ancestry healthy panel.
* **Imputation-platform calibration.** The cross-Korean concordance (NGS vs
  RNA-seq, r≈0.86–0.93) lets future studies estimate how much of an observed
  cohort-vs-cohort delta is platform noise vs biology.
* **Sub-Asian ancestry stratification.** PC1/PC2 from the AFND-shared PCA
  (Section 3.2) can be used as nuisance covariates when combining East
  Asian sub-populations.
* **Public-data sample-size planning.** `master_population_metadata.tsv`
  documents per-source `n_individuals_approx` so a downstream analyst can
  pre-power a sub-allele test before cohort recruitment.

This atlas is **not** suitable for:

* Disease-association testing of any kind. (Out of scope by boundary.)
* Tumor / cancer HLA inference. (Boundary: cancer cohorts forbidden as
  HLA-typing substrate when used to make HLA-disease claims.)

---

## 9. Future cohorts to add

* **Full AFND pull** — all 4-digit class-I (A/B/C) and class-II
  (DRB1/DQA1/DQB1/DPA1/DPB1) alleles for Korea, Japan, Han, Taiwan, Mongolia,
  Vietnam, Thailand, Philippines, Indonesia. Lifts the 5-allele cap.
* **Japanese RNA-seq healthy reference** — GTEx-Asia, JPN-MoCA, or any
  Japanese normal-thyroid RNA-seq cohort would let us replicate the
  cross-platform concordance test in a non-Korean East Asian population.
* **Han Chinese RNA-seq healthy reference** — analogous (HCA-China).
* **Korean Genome Project (KOGES)** — 5,802 Koreans imputed from SNP arrays
  via Kim KOR_REF (memory: HLA imputation reference panel). Supplements
  Baek2021's 173 NGS individuals.
* **Vietnamese / Thai / Mongolian** — to triangulate East Asian / Southeast
  Asian boundary alleles.
* **HapMap3 CHB+JPT** — already on disk
  (`project/external_refs/kim_korean_hla_ref/HapMap3_CHB_JPT/`). A small
  imputed addition to round out the Pan-Asian frame.

---

## Files (relative to `project/results/hla_deepdive_2026_05_08/track3_pan_asian_atlas/`)

### Tables
* `master_population_x_allele_freq.tsv`
* `master_population_metadata.tsv`
* `per_locus_top_alleles_long.tsv`
* `Baek2021_korean_ngs_2field.tsv`
* `K2_normal_genotypes_2field.tsv`
* `GSE213647_normal_genotypes_2field.tsv`
* `afnd_subpopulation_long.tsv`
* `pca_coords.tsv`, `pca_coords_afnd_shared.tsv`, `pca_coords_korean_only.tsv`
* `umap_coords.tsv`
* `jsd_distance_matrix.tsv`
* `cross_korean_concordance.tsv`
* `korean_vs_eastasian_signature_alleles.tsv`
* `korean_vs_pairwise_signature_alleles.tsv`
* `locus_diversity.tsv`
* `track3_summary.json`

### Figures
* `F1_pca_populations.png` — primary PCA (all alleles, coverage-aware)
* `F1_pca_populations_afnd_shared.png` — bias-corrected PCA (5 AFND alleles)
* `F1_pca_populations_korean_only.png` — Korean-only platform PCA
* `F2_umap_populations.png`
* `F3_dendrogram_euclidean.png`
* `F4_dendrogram_jsd.png`
* `F5_per_locus_forest_top5.png` — 6-panel grid, top-5 alleles per locus
* `F6_korean_signature_bars.png` — Korean (Baek2021) Δ vs each comparator
* `F7_heatmap_top50_alleles.png` — top-50 alleles × 9 populations with
  row+col dendrograms
* `F8_locus_diversity_bars.png` — Shannon H + expected heterozygosity bars
* `F9_cross_korean_concordance.png` — pairwise scatter for 3 Korean
  platforms, Pearson + Spearman annotated

### Boundary statement
Healthy-only resource. No disease association tested. No cancer-outcome
join. Downstream users must reconcile their use against the boundary
contract at `project/paper2_hla_boundary/HLA_CANCER_SEPARATION_RULES.md`
before any HLA × cancer analysis.
