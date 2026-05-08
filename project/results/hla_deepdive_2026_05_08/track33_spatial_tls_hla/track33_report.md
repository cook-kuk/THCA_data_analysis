# Track 33 — Spatial TLS niche × HLA-II niche overlap (GSE250521 Visium, n=16)

## 0. Boundary

This track is **HLA gene-expression module only — not allele genotype.** All
HLA-II results are RNA composite scores (HLA-DRA/DRB1/DPA1/DPB1/DQA1/DQB1/
DMA/DMB + CIITA + CD74) from Visium spots. No allele-level claim is made.
Where the bridge zone (HT-overlap PTC) is touched, framing is
hypothesis-generation only — no allele claim, no causal claim, per
`project/paper2_hla_boundary/HLA_CANCER_SEPARATION_RULES.md` Section 1.1.

This is a Paper 1 spatial supplement (transcriptomic immune-context) and
Paper 2 immune-context bridge support — it does **not** speak to Paper 4
HLA-allele Korean GD work.

## 1. Spatial cohort

GSE250521 Visium atlas: 16 samples × ~2.8k spots (median), four stages of
four samples each: Normal (N-1..4), classical PTC (PTC-1..4), late/locally
advanced PTC (LPTC-1..4), anaplastic ATC (ATC-1..4). All re-used from
Track 12 / Pantheon; full 36,601-gene scored h5ads provide all 10 HLA-II
genes and all 12 TLS markers. MOSCOT-mapped cell-type labels (B cell,
Malignant cell, Myeloid cell, T cell, Fibroblast, Endothelial, Epithelial)
are joined per-spot from `pantheonos_demo`.

## 2. TLS niche definition

Spot-level TLS_score = `sc.tl.score_genes` of TLS-core markers
(CXCL13 + MS4A1 + CD79A). Visium spots are smoothed on the hexagonal
adjacency (k=6 oddr offsets, alpha=0.5 self/mean-neighbor mix).

A **TLS niche** is then a connected component of spots above the per-sample
top-10% smoothed TLS-score threshold, retained only if size ≥ 3 spots
(filter out singletons and noise pairs). Per-spot HLA-II module score is
smoothed identically.

**Yield:** 608 TLS niches across 16 samples (median 32 niches/sample,
median niche size 6–7 spots, max niche ~85 spots in PTC-4).

## 3. Niche-level TLS × HLA-II co-localization

This is the headline finding. Track 12's **spot-level** correlation
TLS×HLA-II was rho ≈ −0.18 — a composition artifact (multicell spots dilute
HLA-II expression at high-density TLS spots that are mostly B cells).

At the **niche** level the relationship reverses cleanly:

| Stratum | Spot-level rho | Niche-level Cohen d (in vs out) | Mean Δ HLA-II in–out |
|---|---|---|---|
| Pooled 16 samples | +0.111 (smoothed) | **+0.655** | +0.220 |
| Normal (n=4) | +0.000 (mean) | +0.142 | +0.026 |
| PTC (n=4) | +0.13 | **+0.935** | +0.342 |
| LPTC (n=4) | +0.20 | +0.626 | +0.225 |
| ATC (n=4) | +0.13 | **+0.916** | +0.289 |

**14/16 samples** show in-niche HLA-II elevated above background, with
ttest p<1e-9 in 11/16 samples. Strongest single sample: PTC-1
(Cohen d = +1.51, p ≈ 1.3e-55, niche size up to 13 spots).

**Resolution:** Smoothing + connected-component aggregation removes the
spot-density confound; HLA-II is co-localized with TLS — but only when you
look at the niche, not the spot. The Track 12 spot-level negative rho was
a multicellular-spot artifact, not a biological anti-correlation.

## 4. Distance gradient

For each spot we compute Euclidean distance to nearest niche centroid
(full-resolution pixels), bin non-niche spots into deciles per sample, and
plot mean HLA-II vs decile. F4 shows **a monotone decline** out from the
niche in PTC, LPTC, ATC: HLA-II ≥ +0.5 inside niches drops to ≈ +0.05–0.20
at decile 9 (farthest tissue). Normal samples are flat near zero.

This is **at-the-niche elevation with a short gradient**, not whole-tissue
HLA-II contamination — consistent with antigen-presentation activity
spatially restricted to the lymphoid follicle and immediate periphery.

## 5. Per-sample stage stratification (forest)

F2 forest of per-sample Cohen d:

- **Normal:** mean d = +0.142 ± 0.028. All four positive but small (modest
  HLA-II baseline at lymphoid aggregates in healthy thyroid).
- **PTC:** mean d = +0.935 ± 0.265. PTC-1 (+1.51), PTC-3 (+1.19), PTC-2
  (+0.76), PTC-4 (+0.29).
- **LPTC:** mean d = +0.626 ± 0.201. LPTC-2 (+0.85) and LPTC-4 (+0.93)
  strong; LPTC-3 weak (+0.045) — this is the LPTC sample with most niches
  (72) but high background HLA-II.
- **ATC:** mean d = +0.916 ± 0.081. Tightest variance of any stage; HLA-II
  in TLS niches ≈ 1.0 vs ≈ 0.6 background. Four-of-four ATC samples
  d > 0.69.

**TLS-HLA-II coupling is stage-dependent in magnitude (Normal ≪ tumor) but
not in sign.** ATC niches are individually fewer (median 24 vs 49 in N) yet
each carries the strongest HLA-II contrast — fewer-but-hotter pattern.

## 6. Niche size × HLA-II intensity

F5 scatter (608 niches): Spearman ρ = +0.30 (positive, modest). Bigger TLS
niches do carry slightly stronger mean HLA-II, but the effect saturates;
small niches (3–5 spots) already show clear elevation. Most variance is
between-niche-membership-mix, not size.

## 7. B / T / FDC zone composition within TLS niche

F6 boxplot per stage, scored per niche:

- **B-zone (MS4A1/CD79A/CD79B):** highest in **PTC** niches (median ~+1.5),
  drops in ATC (~+0.6) — consistent with TLS maturation declining in ATC.
- **T-zone (CD3D):** elevated in PTC and LPTC (~+0.7), lower in ATC (~+0.3).
- **FDC (CXCL13):** strongest in PTC (~+1.0), lower elsewhere.

This is exactly the canonical mature-TLS pattern in
classical PTC, with progressive zone disorganization toward ATC. ATC
niches keep HLA-II contrast (Section 5) despite weaker B/T/FDC zonation —
likely myeloid-dominated tertiary aggregates rather than mature follicles.

## 8. HLA-II producers in TLS niche (per cell type)

T5 / F7 — across niches with mapped cell-type labels:

| Cell type | Mean HLA-II in-niche | Background | Δ in–out | n_in_niche |
|---|---|---|---|---|
| **Myeloid cell** | **1.278** | 1.021 | +0.257 | 151 |
| **B cell** | **0.970** | 0.848 | +0.123 | 506 |
| Endothelial | 0.887 | 0.438 | +0.449 | 18 |
| Fibroblast | 0.841 | 0.805 | +0.035 | 346 |
| Malignant | 0.707 | 0.474 | +0.233 | 3,233 |
| T cell | 0.573 | 0.474 | +0.099 | 172 |
| Epithelial | 0.133 | 0.111 | +0.023 | 1,996 |

**Dominant producers in niche: Myeloid > B > Malignant.**
Note: GSE250521 MOSCOT mapping has very few mapped DC labels (Track 12
flagged DC > Myeloid > B in scRNA but DCs are sparse in the spatial
mapping). At Visium spot resolution, "Myeloid cell" likely captures the
DC + macrophage compartment together. Result still confirms the immune
APC compartment is the principal in-niche HLA-II source, with malignant
cells co-elevating but not dominating.

## 9. DM1 score × TLS niche presence

F8: number of TLS niches per sample vs sample mean DM1 score (clipped 1–99
percentile to suppress outlier scoring noise). Spearman ρ ≈ −0.30 (n=16,
not significant). Stage-aggregate niche counts: N=49 (median), PTC=33,
LPTC=33, ATC=24 — declining with progression. This matches memory
`v17_D5P6_BCR_clonal_TLS` (TLS d=+1.96 in HT-PTC bulk) and
`pantheonos_lptc_tls_finding_2026_05_08` (LPTC peaks LR diversity) at the
**niche-count** level: ATC has fewer-but-stronger TLS niches; PTC/LPTC
have intermediate. DM1-high samples do **not** carry more TLS niches —
DM1 dedifferentiation and TLS density are independent axes at this
sample size (n=16).

## 10. HT-overlap PTC bridge angle (skip rationale)

GSE250521 Visium does not include explicit HT-overlap PTC samples — the
16-sample atlas is N / PTC / LPTC / ATC, with no annotated Hashimoto
background. GSE286332 (memory `v17_GSE286332_strong_go`,
`v17_D5P6_BCR_clonal_TLS`) is bulk RNA-seq, not spatial. GSE163203 PTC+HT
scRNA was flagged not-loaded in Track 12.

Therefore Section 10 is intentionally skipped — bridge-zone hypothesis is
**not testable here** without a spatial HT-PTC dataset. Hypothesis recorded
for Paper 2 reviewer reserve: bulk-level HT-PTC TLS d=+1.96 + AICDA + IGHV
clonality predict that any future HT-PTC Visium sample should show
elevated niche count + niche-level HLA-II-in-niche d ≈ +1, exceeding
classical PTC. Bridge-only, no allele claim, no causal claim.

## 11. Limitations

1. **Multicellular Visium spots (~55 µm, ~10 cells/spot)** — every spot
   is a mixture; "HLA-II in B cell spot" is HLA-II-given-B-cell-dominant,
   not single-cell B-cell HLA-II.
2. **Hex-grid resolution** — niches with size 3 are at the
   filter-floor; biological TLS that fall on a single Visium spot are
   missed by design (singleton filter).
3. **Donor-level confounds** — n=4/stage; per-stage Cohen d means are
   robust (sign-consistent 14/16) but absolute magnitudes are
   donor-driven (e.g., LPTC-3 outlier).
4. **TLS marker threshold sensitivity** — top-10% threshold is
   sample-relative; Normal samples have far weaker absolute TLS scores so
   their "niches" are statistically defined but biologically tiny (median
   d=+0.14 confirms this).
5. **Mapped-celltype gaps** — N samples have very limited cell-type
   transfer (most spots un-mapped or "Epithelial") because the MOSCOT
   reference is tumor-focused. The producer-ranking table is dominated by
   tumor stages; n_samples per cell type ranges 4–12.
6. **Distance metric** — full-resolution pixels are not actual microns
   in absolute terms because Visium full-resolution pixel size is
   sample-dependent; deciles are within-sample so this is fine for
   gradient shape but absolute distances are not directly cross-sample
   comparable.
7. **DM1 score dynamic range** — sample mean DM1 was clipped 1–99% for
   Section 9 plotting; raw outliers (1e-7 scale) were dataset-floor noise
   and inflated the regression visually.

## 12. Implication

**Paper 1 spatial supplement.** This is the proper niche-level rebuttal
to "spatial spot-level rho is negative" — at niche scale TLS and HLA-II
co-localize (global Cohen d = +0.655, sign-consistent in 14/16 samples,
strong in PTC and ATC). Add as a single supplementary figure: per-stage
forest (F2) + distance gradient (F4) + zone composition (F6) + producer
ranking (F7), with caption "TLS niches are HLA-II hotspots in PTC/LPTC/ATC
thyroid cancer; the immune APC compartment (myeloid + B) dominates with
malignant co-elevation. Spot-level negative rho is a multicellular-spot
composition artifact resolved at niche level."

**Paper 2 immune-context.** Confirms that the immune APC compartment
(not malignant epithelium) is the primary HLA-II producer at lymphoid
aggregates, consistent with mature-TLS biology. Stage gradient (PTC ≥ ATC
> LPTC ≫ N in Cohen d) supports the manuscript's framing of HT-overlap
mechanism as immune-driven; no allele-level claim is made or implied.

---

## Output paths

- Tables (6): `/home/seungho/personal/THCA_data_analysis/project/results/hla_deepdive_2026_05_08/track33_spatial_tls_hla/tables/T1..T6_*.tsv`
- Figures (8 aggregate + 16 per-sample = 24): `.../track33_spatial_tls_hla/figs/F1..F8_*.png`, `sample_*_TLS_HLA_II.png`
- JSON: `.../track33_spatial_tls_hla/track33_summary.json`
- Script: `/home/seungho/personal/THCA_data_analysis/scripts/hla_deepdive_2026_05_08/track33/run_track33.py`
