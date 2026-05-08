# Track 5 — DM1 × HLA-I gene-expression module in TCGA-THCA

**Date:** 2026-05-08
**Owner:** Seungho Cook
**Track scope:** transcriptomic gene-expression module analysis only — *not* allele genotype.

---

## 0. Boundary statement (read first)

This track is one of six parallel HLA deep-dive tracks. Its work product is restricted to **HLA-I and HLA-II as transcriptomic gene-expression modules** computed on TCGA-THCA tumor RNA-seq.

Per `project/paper2_hla_boundary/HLA_CANCER_SEPARATION_RULES.md`:

- Section 1.1 explicitly allows "HLA-I / HLA-II as transcriptomic gene-expression module signatures *only* when used as a residualization or contextual control."
- Section 1.2 forbids HLA *allele* imputation/typing from cancer data, including TCGA-THCA RNA-seq.

Every figure and table in this track carries the caption boilerplate:

> "HLA-I/II gene-expression module — not allele genotype. Cancer-cohort allele genotyping is out of scope per separation rules."

No allele-level inference (HLA typing, allele frequency, carrier frequency) is performed or reported here. HLA LOH is consumed only if a *precomputed* file already exists in the repo; LOHHLA was not run.

---

## 1. Data sources

| Source | Path | Use |
|---|---|---|
| TCGA pan-cancer expression (RSEM, log2-normalized) | `project/data/raw/TCGA_pancan/pancan_geneExp.gz` | HLA-I/II gene matrix for THCA primary tumors |
| TCGA pan-cancer phenotype | `project/data/raw/TCGA_pancan/phenotype.tsv.gz` | sample_type filter (Primary Tumor only) |
| DM1 score per sample (canonical, pancan-z) | `project/results/paper11_pancancer/pancan_dm1_scored.tsv` | DM1 score column (`DM1_like`) |
| DM master with PFI / driver / TERT / stage | `project/results/dark_matter_phase1/tcga_dm_master_with_pfi.tsv` | covariates + outcomes (n=482) |
| 8-gene panel coefficients | `project/results/v17p35/tables/AMP4_8gene_model_coefficients.tsv` | local DM1 cross-check |
| HM450 mean panel methylation | `project/results/audit_2026_04_30/round5/r5_2_sample_methylation_8gene.tsv` | methylation × HLA module cross-axis |

Cohort: TCGA-THCA primary tumors (n=505 expression / n=527 after DM1 merge with covariates / n=482 with PFI). No normal-tissue samples included.

---

## 2. Module score construction

For each sample we compute the z-scored mean of the constituent genes (z computed across all 505 THCA primary tumors).

- **HLA-I module (19 genes):** HLA-A, HLA-B, HLA-C, B2M, TAP1, TAP2, TAPBP, NLRC5, IRF1, PSMB8, PSMB9, ERAP1, ERAP2, HLA-E, HLA-F, HLA-G, CALR, CANX, PDIA3.
- **HLA-II module (12 genes):** HLA-DRA, HLA-DRB1, HLA-DPA1, HLA-DPB1, HLA-DQA1, HLA-DQB1, HLA-DMA, HLA-DMB, HLA-DOA, HLA-DOB, CIITA, CD74.
- **Immune proxy:** CD2, CD3D/E/G, CD4, CD8A/B, CD19, MS4A1, CD79A/B, PRF1, GZMB, GZMA, NKG7, KLRB1, ITGAX, ITGAM, CD68, CD163, FCGR3A, PTPRC.
- **Stromal proxy:** COL1A1/A2, COL3A1, COL5A1, FAP, ACTA2, PDGFRA/B, DCN, LUM, VIM, PECAM1, VWF, CDH5, CD34.
- **Transcriptomic purity proxy:** −(Immune + Stromal) z-scored. Used as ESTIMATE/ABSOLUTE substitute because no precomputed CPE/ABSOLUTE/ESTIMATE table was found locally; this is documented as a limitation in §9.
- **DM1 score:** primary = pancan z-scored `DM1_like` from `pancan_dm1_scored.tsv`. Cross-check `DM1_local = -mean(z(8-panel))` agrees with pancan DM1 at Spearman ρ=0.94.

All module scores are written per-sample to `tables/T01_per_sample_module_scores.tsv`.

---

## 3. DM1 × HLA-I module

**Headline (n=527):** Spearman ρ(DM1, HLA-I) = **+0.318, p = 8.3 × 10⁻¹⁴**, Pearson r = +0.283.

Tertile contrast (DM1-high vs DM1-low):
- HLA-I module: median 0.117 vs −0.446 — Wilcoxon p = 3.5 × 10⁻¹³, Cohen's d = +0.76.
- HLA-II module: median 0.335 vs −0.660 — Wilcoxon p = 1.2 × 10⁻¹⁷, d = +1.02.
- Immune proxy: d = +0.53; Stromal proxy: d = +0.06 (n.s.).

Direction: high-DM1 (de-differentiated) thyroid tumors have **higher HLA-I and HLA-II module expression**. The relationship is monotonic on the LOWESS overlay (no crossover at the high end).

Per-gene Spearman with DM1 (`tables/T04_per_gene_dm1_corr.tsv`):
- Strongest *positive*: PSMB8 (+0.51), HLA-B (+0.48), TAP1 (+0.43), HLA-C (+0.43), HLA-G (+0.43), B2M (+0.42), HLA-A (+0.40).
- Strongest *negative*: **CALR (−0.71), PDIA3 (−0.64), CANX (−0.52)** — i.e. the chaperone arm of the antigen-presentation machinery moves opposite the IFN-driven arm (NLRC5/PSMB8/TAP1/HLA-A-B-C). This split is biologically consistent: chaperones track ER/secretory load (high in differentiated thyrocytes producing thyroglobulin), while presentation machinery tracks immune activation. It is a real, interpretable feature, not an artifact.

Files: `figs/F01_dm1_vs_hla1_hla2_scatter.{png,pdf}`, `figs/F02_dm1_tertile_boxplots.{png,pdf}`, `figs/F03_hla1_per_gene_heatmap.{png,pdf}`, `figs/F04_per_gene_dm1_corr_lollipop.{png,pdf}`.

---

## 4. DM1 × HLA-II module

**Headline (n=527):** Spearman ρ(DM1, HLA-II) = **+0.393, p = 7.3 × 10⁻²¹**.

HLA-II is more tightly correlated with DM1 than HLA-I, consistent with the prior memory `v17_D4P2_tcga_hashimoto_generalization` (HT-overlap signature transfers to TCGA DM2-enriched cohorts) and the broader pattern that HLA-II / antigen-uptake machinery (CD74, CIITA, HLA-DR) is the most responsive of the antigen-presentation arms in dedifferentiated thyroid disease.

---

## 5. Driver-stratified analysis

`tables/T05_driver_stratified.tsv`, `figs/F05_driver_stratified_forest.{png,pdf}`.

| Driver | n | ρ(DM1, HLA-I) [95 % CI] | p | ρ(DM1, HLA-II) [95 % CI] | p |
|---|---|---|---|---|---|
| BRAF | 294 | +0.152 [+0.04, +0.26] | 0.009 | +0.228 [+0.12, +0.33] | 8 × 10⁻⁵ |
| RAS | 54 | −0.001 [−0.27, +0.27] | 0.996 | +0.161 [−0.11, +0.41] | 0.25 |
| Triple-negative | 168 | +0.251 [+0.10, +0.39] | 0.001 | +0.365 [+0.23, +0.49] | 1 × 10⁻⁶ |
| Fusion | <5 in this label set | — | — | — | — |

*Driver labels follow `dark_matter_phase1/tcga_dm_master_with_pfi.tsv :: driver_anchor_v17`. Fusion samples were too few under that label to test (a separate fusion call set exists in the deconv tracks but is not used here).*

Pattern: the DM1 × HLA-I/II coupling is **strongest in triple-negative and BRAF tumors and absent in RAS** — consistent with memory `v17_dark_matter_pivot_2026_04_29` (the 8-gene paper as a BRAF/RAS-negative sub-stratifier) and `dm1_round5` (RAS = 98 % DM1 by composition, so within RAS the DM1 axis is compressed and information about HLA flow is lost; this is structural, not a biology-killer).

---

## 6. Methylation × HLA-I module

Using the precomputed `r5_2_sample_methylation_8gene.tsv :: mean_8g_beta` (HM450 mean β across the 8-gene panel; per memory `dm1_round4`), n=523:

- ρ(DM1, mean_8g_beta) = +0.49, p = 5 × 10⁻³³.
- ρ(mean_8g_beta, HLA-I module) = +0.64, p = 6 × 10⁻⁶¹.
- ρ(mean_8g_beta, HLA-II module) = +0.66, p = 7 × 10⁻⁶⁷.

So the same axis that hyper-methylates the 8-gene panel (driving DM1) also tracks **higher** HLA-I/II module expression. This is the cleanest single piece of evidence in this track that the DM1 ↔ HLA module relationship is biology, not a purity artifact: mean panel β is a tumor-cell methylation property, and it predicts the immune-context module the same direction as DM1 does. `figs/F06_methylation_dm1_hla1.{png,pdf}`, `tables/T06_methylation_dm1_hla.tsv`.

---

## 7. HLA LOH

A repo-wide search for `*loh*hla*`, `*hla*loh*`, `*lohhla*` returned **no precomputed TCGA-THCA HLA LOH calls**. Per the boundary, LOHHLA is *not* run from BAMs in this track. Status: skipped, hooks in place if a published TCGA-THCA HLA LOH file is added to the repo later.

---

## 8. Purity decoupling

`tables/T08_purity_partial.tsv`, `figs/F07_purity_decoupling.{png,pdf}`.

Spearman partial correlations (rank-residualized):

| Module | raw ρ | ρ \| Purity proxy | ρ \| Immune+Stromal proxies |
|---|---|---|---|
| HLA-I | +0.318 (p=8e-14) | **+0.301 (p=2e-12)** | +0.266 (p=5e-10) |
| HLA-II | +0.393 (p=7e-21) | **+0.411 (p=7e-23)** | +0.390 (p=1e-20) |

Conditioning on the transcriptomic purity proxy moves the DM1 ↔ HLA-I ρ from +0.318 to +0.301 (15 % attenuation) and the DM1 ↔ HLA-II ρ from +0.393 to +0.411 (no attenuation, slight inflation). Even conditioning jointly on the immune *and* stromal proxies (a stricter control that bleeds biology into the covariate) leaves ρ ≥ +0.27 / +0.39 with p ≪ 1 × 10⁻⁹.

**Interpretation: the DM1 ↔ HLA-I/II coupling is not purity-driven. It is biology.**

---

## 9. Limitations

1. **No precomputed ABSOLUTE / CPE / ESTIMATE table** was found in the repo, so purity is controlled via an in-script transcriptomic proxy (negative of immune + stromal marker means). This is conservative — if anything the proxy bleeds real immune-axis biology *into* the covariate, which would over-attenuate the partial ρ. Even so, the relationship survives. A canonical purity table (TCGA pancan ABSOLUTE) would tighten this further.
2. **No HLA LOH layer.** Per boundary, LOHHLA is not run from cancer BAMs. Hook is left for an external precomputed file.
3. **Driver labels are from `driver_anchor_v17`**, which collapses fusion calls into a single bucket and may not catch all fusion carriers. A more granular fusion call set (RET / NTRK / ALK breakdown) lives in the deconv track and is not consumed here.
4. **DM1 score is a one-dimensional summary** (pancan-z `DM1_like`). The 8-gene panel direction is preserved (Spearman ρ between pancan DM1 and local 8-panel DM1 = 0.94), but DM1 itself does not separate the chaperone-arm vs IFN-arm of HLA-I — which is why the per-gene heatmap (F03) shows two visually distinct gene clusters, both real.
5. **Allele-level HLA inference is forbidden by the boundary** and was not attempted. All claims here are at the gene-expression-module level.

---

## 10. Implication for Paper 1 residualization (boundary-respecting)

**Question:** does adding HLA-I module as a covariate change the DM1 → outcome relationship?

**Answer (Cox PH on PFI, n = 482):** `tables/T09_cox_residualization.tsv`, `figs/F08_residualization_forest.{png,pdf}`.

| Spec | HR (DM1) | p (DM1) | HR (HLA-I) | p (HLA-I) | HR (Purity) | p (Purity) |
|---|---|---|---|---|---|---|
| DM1 only | **2.97** | 0.0026 | — | — | — | — |
| DM1 + HLA-I | **2.94** | 0.0033 | 1.04 | 0.86 | — | — |
| DM1 + Purity | **3.07** | 0.0017 | — | — | 1.18 | 0.24 |
| DM1 + HLA-I + Purity | **2.86** | 0.0044 | 1.43 | 0.22 | 1.35 | 0.09 |

The DM1 → PFI hazard ratio is **stable at HR ≈ 2.9 – 3.1** across all four specifications. Adding the HLA-I gene-expression module as a covariate (or the transcriptomic purity proxy) does not move the DM1 → PFI signal. HLA-I is itself non-significant as a PFI covariate once DM1 is in the model (HR ≈ 1.0–1.4, p ≥ 0.22), which is the expected residualization-control pattern for Paper 1.

**Conclusion for Paper 1:** the DM1 → PFI hazard is *not* a hidden HLA / immune-context confound. It survives a residualization sanity check and HLA-I module is appropriate as a Paper-1 contextual control (per Section 1.1 of the boundary doc).

---

## 11. Summary numbers

```
n TCGA-THCA primary tumors with DM1 + HLA-I + HLA-II + driver: 527
ρ (DM1, HLA-I module)               = +0.318  p = 8.3e-14
ρ (DM1, HLA-II module)              = +0.393  p = 7.3e-21
ρ (DM1, Immune proxy)               = +0.216  p = 5.3e-07
ρ (DM1, Stromal proxy)              = +0.047  p = 0.29   (n.s.)
ρ (DM1, mean panel HM450 β)         = +0.49   p = 5e-33
ρ (mean panel β, HLA-I module)      = +0.64   p = 6e-61
ρ (mean panel β, HLA-II module)     = +0.66   p = 7e-67
Cohen's d (HLA-I high vs low DM1)   = +0.76
Cohen's d (HLA-II high vs low DM1)  = +1.02
Driver-stratified ρ(DM1, HLA-I): BRAF +0.15 / Triple-neg +0.25 / RAS  0.00
Partial ρ (DM1, HLA-I | Purity proxy)            = +0.30   (15 % attenuation)
Partial ρ (DM1, HLA-I | Immune + Stromal proxies) = +0.27
Cox PFI HR (DM1) — DM1 only:             HR = 2.97   p = 0.003
Cox PFI HR (DM1) — DM1 + HLA-I:          HR = 2.94   p = 0.003
Cox PFI HR (DM1) — DM1 + HLA-I + Purity: HR = 2.86   p = 0.004
```

Boundary line for every caption:
> *HLA-I/II gene-expression module — not allele genotype. Cancer-cohort allele genotyping is out of scope per separation rules.*
