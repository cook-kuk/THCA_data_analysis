# Paper 1 — Supp Fig SX (deconvolution v5) — manuscript drop-in

**Status:** v5 ABCE complete (D running). Author keyboard required for voice-protected
sections (Discussion 3.1, Limitations, Cover Para 1, Q9 in `09_reviewer_qa.md`).
This file is a factual drop-in: figure caption, STAR Methods, results bullets.
None of these touch the voice-protected list.

## Updated Supp Fig SX caption (panels A–K, ready to paste into figure_captions.md)

> **Supplementary Figure SX. Multi-method bulk cell-type deconvolution and
> compositional axes of the DM1↔DM2 phenotype.**
>
> **(A) Multi-method residualization grid.** Cohen's d (DM1−DM2) for the canonical
> 8-gene RAI score (`rai_score_recalc`) under five residualization stages
> (raw, +stromal [Endo+Fibro], +immune [T+Myel+B+NK], +Epithelial-only, +ALL 8 fractions)
> across four deconvolution methods (NNLS, Ridge-NNLS [L2 α=1], LR-clip, nu-SVR
> [3 ν, lowest-RMSE]); reference = Lu 2023 (GSE193581) `author_celltype` 67,678 cells × 8 cell types in 1,898 HVG ∩ TCGA gene symbols.
>
> **(B) Effect-retention ratio after full residualization.** NNLS / Ridge 24% (sparse-weight artifact); LR-clip 70%; nu-SVR 47% (matches canonical S4 56% retention).
>
> **(C) Per-cell-type DM1 vs DM2 fraction Cohen's d.** Direction-consistent enrichment of Malignant↑, Myeloid↑, Epithelial↓ in DM1 across all four methods.
>
> **(D) nu-SVR primary mean cell-type fractions** in TCGA-THCA (n=572).
>
> **(E) Methodology and caveats panel.**
>
> **(F) Cross-cohort direction consistency** TCGA nu-SVR Cohen's d (DM1−DM2) vs Lee/GSE213647 (n=632) Spearman ρ vs `panel_z` (sign-flipped to align with DM1-direction). Direction-consistent for 7/8 cell types; T-cell discordance reflects binary-vs-continuous score frame.
>
> **(G) Within-DM1 fusion+ vs fusion− cell-type fraction Cohen's d** (n=74 vs 391, kinase fusion: RET/NTRK/ALK/BRAF/PAX8/PPARG aggregated from cBioPortal SV table). All |d| ≤ 0.42 — within-DM1 fusion+/− tumors share near-identical compositions, direct support for the fusion-independent epigenetic silencing claim (Fig 8 / §2.4a).
>
> **(H) Per-driver-class TCGA-THCA cell-type composition** (BRAF V600E n=280 / RAS mutant n=54 / driver-negative n=179): RAS-mutant tumors are Epithelial-cluster–enriched (d_RAS−BRAF=+1.52) and Malignant-cluster–depleted (d=−1.66) versus BRAF V600E, consistent with Landa 2016 / Paper 1 §2.1 differentiation framing.
>
> **(I) HM450 8-gene mean β × cell-type fraction Spearman ρ heatmap** (TCGA n≈484 paired): mean β positively tracks Myeloid (ρ=+0.39), B cell (+0.27), Fibroblast (+0.23) and Malignant (+0.30); negatively tracks T cell (ρ=−0.42) and Epithelial (−0.31), connecting the Round-4 methylation layer (DM1 vs DM2 mean-β d=−1.75) to compartment composition.
>
> **(J) DM1 sub-A vs sub-B cell-type fraction Cohen's d** (sub-A n=84, sub-B n=56): sub-A is Malignant-cell rich (d=+1.22, p=2.5×10⁻⁹) and Epithelial-poor (d=−1.26, p=5.2×10⁻¹⁰); immune-compartment differences are non-significant. **The sub-A/B split is a tumor-purity-vs-thyrocyte split, not immune-hot vs immune-cold** — Paper-2 boundary marker (sub-B = fusion-/mutation-negative Hashimoto-overlap retains thyrocyte identity).
>
> **(K) Cell-type composition pseudotime along the canonical 8-gene score.** TCGA-THCA (n=513, score = `rai_score_recalc`) and Lee/GSE213647 (n=632, score = `panel_z`) samples were ranked by score, binned into 10 deciles; the mean per-decile cell-type fraction (nu-SVR against Lu 2023) defines a 10-step pseudotime trajectory. Four compartments — Malignant (↓), Epithelial (↑), Myeloid (↓), Endothelial (↑) — show monotonic decile-level Spearman ρ ≥ |0.95| in BOTH cohorts, defining a reproducible compositional pseudotime independent of cohort, scoring scheme, or sample size.

## STAR Methods — additions to "Bulk cell-type deconvolution"

> **Per-driver-class composition.** TCGA nu-SVR cell-type fractions were merged with the cBioPortal `v3_anchor_6class` driver call (BRAF V600E / RAS mutant / driver-negative; `fusion_calls_per_sample.tsv`); Cohen's d was computed per cell type for each pairwise contrast.
>
> **Methylation × composition.** Per-sample HM450 mean 8-gene β values (`r5_2_sample_methylation_8gene.tsv`, TCGA-THCA n≈484) were correlated (Spearman) with per-sample cell-type fractions for each cell type. Per-gene 8-gene heatmap reports gene × cell-type Spearman ρ.
>
> **DM1 sub-A vs sub-B teaser.** Sub-cluster labels from `d6p7_dm1_subcluster/dm1_subcluster_labels.tsv` were intersected with cell-type fractions; Cohen's d and Mann-Whitney U two-sided p were reported per cell type.
>
> **Pseudotime trajectory.** Samples were ranked by canonical 8-gene score (TCGA `rai_score_recalc`, Lee `panel_z`), binned into 10 deciles, and mean per-decile cell-type fraction was computed. Decile-level Spearman ρ between mean score and mean fraction quantifies monotonic trajectory.
>
> **Full-transcriptome reference robustness.** Pu 2021 raw counts (33,694 genes × 66,015 cells) were subsampled to 5,000 cells balanced across 7 patients (seed=42); each cell was assigned a Lu 2023 `author_celltype` label by maximum cosine similarity over the Lu HVG ∩ Pu intersection (1,898 genes after Ensembl→symbol conversion). A full-transcriptome pseudobulk per cell type was constructed (per-cell-type log-normalized mean, library-size-corrected, n_genes=33,694; `Pu_pseudobulk_full`). TCGA bulk was re-deconvolved by NNLS over the 21,369-gene Pu × TCGA intersection (`scipy.optimize.nnls`, sum-to-one normalization). Per-cell-type Cohen's d (DM1−DM2) was compared to the v2 Lu HVG nu-SVR primary result; all four informative axes (Malignant, Epithelial, Myeloid, Endothelial) sign-match between the two references with Pu full NNLS magnitudes ≥ Lu HVG. NNLS sparse-collapse zeros the B/Fibroblast/NK/T compartments in the full-transcriptome regime; these are interrogated by the v2 nu-SVR Lu HVG primary instead.

## Results bullets (paste into appropriate Results sub-section, NOT voice-protected)

- 8-gene methylation β positively tracks Myeloid (ρ=+0.39) and Malignant (+0.30); negatively tracks T cell (−0.42) and Epithelial (−0.31), connecting the methylation layer (Round-4: DM1 vs DM2 mean-β d=−1.75) to compartment composition (Supp Fig SX I).
- BRAF V600E and RAS-mutant tumors carry distinct compositional signatures: RAS retains Epithelial-cluster identity (d_RAS−BRAF=+1.52) while BRAF is Malignant-cell rich (d=−1.66), consistent with the Landa 2016 / §2.1 differentiation framing (Supp Fig SX H).
- DM1 sub-A vs sub-B is a *tumor-purity-vs-thyrocyte* split (Malignant d=+1.22 / Epithelial d=−1.26, p<10⁻⁹) — not immune-hot vs immune-cold; immune-compartment differences are NS. This is the Paper-2 boundary marker (Supp Fig SX J).
- Decile-binned pseudotime trajectory along the canonical 8-gene score reproduces 4 monotonic compositional axes in both TCGA (n=513) and Lee/GSE213647 (n=632): Malignant ↓, Epithelial ↑, Myeloid ↓, Endothelial ↑ (decile Spearman |ρ| ≥ 0.95 both cohorts; Supp Fig SX K).

## Files (v5)

| File | Purpose |
|---|---|
| `run_deconv_v5_abc.py` | A/B/C pipeline (driver / methylation / sub-A_B) |
| `run_deconv_v5_d_fast.py` | D pipeline (Pu full NNLS + Pu top-10K nu-SVR) |
| `run_deconv_v5_e_trajectory.py` | E pipeline (TCGA + Lee decile pseudotime) |
| `plot_deconv_v5.py` | A/B/C figure (Fig_SX_deconvolution_v5_abc) |
| `plot_deconv_v5_composite.py` | A–E composite (Fig_SX_deconvolution_v5_composite) |
| `Fig_SX_deconvolution_v5_abc.{png,pdf}` | A/B/C composite |
| `Fig_SX_deconvolution_v5_e_trajectory.{png,pdf}` | E pseudotime |
| `Fig_SX_deconvolution_v5_composite.{png,pdf}` | full v5 composite for manuscript |
| `v5A_per_driver_class.tsv` | A: cell-type × driver class means + d's |
| `v5B_methylation_celltype_corr.tsv` + `v5B_methylation_celltype_pivot.tsv` | B: gene × cell-type ρ |
| `v5C_dm1_subA_subB_celltype.tsv` | C: cell-type × subcluster d/p |
| `v5E_trajectory_TCGA.tsv`, `v5E_trajectory_Lee.tsv` | E: decile × cell-type means |
| `v5D_dm1_dm2_full_pu_NNLS.tsv` | D1: NNLS full-transcriptome DM d |
| `v5D_dm1_dm2_top10k_pu_nu_SVR.tsv` | D2: nu-SVR top-10K DM d |
| `v5D_three_way_concordance.tsv` | D3: Lu HVG / Pu full / Pu top-10K 3-way table |
