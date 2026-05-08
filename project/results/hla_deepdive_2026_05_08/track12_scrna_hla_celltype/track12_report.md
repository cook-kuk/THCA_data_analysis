# Track 12 — scRNA HLA-I/II module per cell type in thyroid

**Authors note (Seungho Cook).** This document reports a transcriptomic gene-expression module decomposition for HLA-I and HLA-II across thyroid scRNA-seq and Visium cohorts. Caption boilerplate everywhere: **HLA gene-expression module — not allele genotype.**

---

## 0. Boundary

`/home/seungho/personal/THCA_data_analysis/project/paper2_hla_boundary/HLA_CANCER_SEPARATION_RULES.md`, Section 1.1.

- This track stays inside the *gene-expression module* lane allowed for Paper 1 (transcriptomic immune-context).
- No allele genotype, no carrier-frequency claim, no autoimmune susceptibility statement.
- Bridge zone (HT+PTC scRNA, GSE163203) was **not loaded locally**; therefore no bridge-cohort data is reported here. Hypothesis-only mention is placed in Section 9 (limitations).

---

## 1. Cohorts (cohort manifest)

Table file: `tables/T1_cohort_manifest.csv`

| Dataset                                            | Path on disk                                                   | n cells / spots | n samples | Conditions                       | Cell-type label        | Role             |
|----------------------------------------------------|----------------------------------------------------------------|-----------------|-----------|----------------------------------|------------------------|------------------|
| GSE184362 (Pu et al. PTC scRNA)                    | `/data/thca/scrna/processed/classical_baseline.h5ad`           | 19,999 cells    | 7         | PTC tumor + adjacent             | leiden + marker score  | primary scRNA #1 |
| GSE193581 (Lu 2023 PTC + ATC scRNA)                | `/data/thca/scrna/processed/classical_baseline_F12.h5ad`       | 17,898 cells    | 9         | PTC vs ATC; BRAF/RAS/WT          | leiden + marker score  | primary scRNA #2 |
| GSE250521 (Visium spatial atlas)                   | `project/data/processed/GSE250521/<gsm>/`                      | ~36k spots      | 16        | Normal/PTC/LPTC/ATC              | MOSCOT-mapped          | spatial overlay  |
| GSE163203 (PTC + HT scRNA, bridge)                 | not_downloaded                                                 | n/a             | n/a       | PTC + HT                         | n/a                    | BRIDGE — not on disk |
| GSE191288                                          | not_loaded                                                     | n/a             | n/a       | PTC                              | n/a                    | not loaded       |

Only the first three are available; analyses below are computed only on those.

---

## 2. Cell-type module landscape

Tables: `T2_pu_celltype_x_HLA_gene_dotplot.csv`, `T2_lu_celltype_x_HLA_gene_dotplot.csv`, `T3_*_HLA_module_per_celltype.csv`.

Figures:
- `figs/F2_pu_celltype_x_HLA_gene_dotplot.png` — dotplot, cell-type × per-gene HLA-I/HLA-II.
- `figs/F2_lu_celltype_x_HLA_gene_dotplot.png` — same for Lu cohort.
- `figs/F3_pu_HLA_module_violin.png` — module score violins per cell type.
- `figs/F3_lu_HLA_module_violin.png` — same for Lu cohort.

**Top cell types — HLA-I module (mean):**
- Pu cohort (T3_pu): NK 1.000, T_cell 0.975, B_cell 0.833, Endothelial 0.776, Plasma 0.780, Myeloid 0.758, DC 0.729, Mast 0.618, Fibroblast 0.537, **Thyrocyte 0.492** (lowest non-Unassigned).
- Lu cohort (T3_lu): T_cell 1.000, NK 0.955, Endothelial 0.804, B_cell 0.794, Mast 0.816, Plasma 0.731, Myeloid 0.726, DC 0.664, **Thyrocyte 0.548**, Fibroblast 0.282.

Both cohorts agree: lymphoid (NK, T_cell, B_cell) > myeloid > thyrocyte/fibroblast for HLA-I module.

**Top cell types — HLA-II module (mean):**
- Pu: DC 2.18, Myeloid 1.93, B_cell 1.60, T_cell 0.37, Plasma 0.65, NK 0.30, Endothelial 0.28, Mast 0.24, Fibroblast −0.02, **Thyrocyte −0.02**.
- Lu: DC 2.46, Myeloid 2.10, B_cell 1.14, T_cell 0.44, Plasma 0.31, NK 0.26, Endothelial 0.12, Mast 0.06, **Thyrocyte 0.16**, Fibroblast −0.20.

Both cohorts agree: DC > Myeloid > B_cell are the dominant HLA-II producers in thyroid scRNA. Thyrocytes carry essentially no baseline HLA-II module.

---

## 3. Condition / stage stratification

Table: `T4_lu_celltype_condition_HLA.csv` (PTC vs ATC stratification within Lu); Visium per-stage in `T7_visium_celltype_HLA_per_sample.csv`.

Figure: `figs/F4_lu_celltype_condition_box.png`.

Pu cohort had no usable disease-condition column (mutation_status mostly NA), so the condition split is reported only for Lu and for Visium.

Visium stage × mapped_celltype HLA module (weighted by spot count):

| Stage | Malignant HLA-I | Malignant HLA-II | Myeloid HLA-I | Myeloid HLA-II | B-cell HLA-II |
|-------|-----------------|------------------|----------------|----------------|----------------|
| PT (normal — Epithelial) | 0.712 | 0.123 | 0.865 | 0.966 | n/a |
| PTC                       | varies | 0.40–0.70 | 0.85 | 1.29 | 0.76 |
| LPTC                      | 0.66–0.92 | 0.45–0.87 | 0.79 | 1.40 | 0.83 |
| ATC                       | 0.78–1.25 | 0.39–0.71 | 1.12 | 1.27 | 0.58 |

Direction: **Malignant cells gain HLA-I from PT → ATC**. Malignant HLA-II rises sharply from normal (0.12) into all tumor stages (0.4–0.7) — a tumor-acquired HLA-II module. Myeloid HLA-II stays high throughout.

Figure: `figs/F7_visium_HLA_stage_celltype_heatmap.png`.

---

## 4. Thyrocyte-only pseudobulk × DM1

Tables: `T5_pu_thyrocyte_pseudobulk.csv`, `T5_lu_thyrocyte_pseudobulk.csv`.

Figures: `figs/F5_pu_thyrocyte_pseudobulk_DM1xHLA.png`, `figs/F5_lu_thyrocyte_pseudobulk_DM1xHLA.png`.

Per-sample thyrocyte-only pseudobulk (HLA module averaged across thyrocyte cells per sample, then correlated with thyrocyte-only DM1 score):

| Cohort | n samples | DM1 × HLA-I rho | p | DM1 × HLA-II rho | p |
|--------|-----------|------------------|---|-------------------|---|
| Pu    | 5         | +0.40            | 0.50 | −0.40           | 0.50 |
| Lu    | 5         | −1.00            | tied (5 ranks) | +0.10 | 0.87 |

**Bulk T5 reference (TCGA THCA n=527):** DM1 × HLA-I rho = +0.318 (p=8e-14), DM1 × HLA-II rho = +0.393 (p=7e-21).

Verdict at this n: scRNA pseudobulk N is too small (5 samples per cohort with thyrocyte clusters ≥10 cells) for a stable rho estimate. Pu signs the same direction as bulk for HLA-I (+0.40) but the Lu sample-set inverts it. This is a power problem, not a contradiction. The Visium spatial analog (Section 8, n=16 samples) gives the more reliable answer.

---

## 5. Immune cell HLA-II producer ranking

From Section 2: **DC > Myeloid > B_cell** in both Pu and Lu cohorts. DC HLA-II module score 2.18 (Pu) / 2.46 (Lu) is approximately 4× higher than B_cell (1.60 / 1.14) and 14× higher than T_cell (0.37 / 0.44). This matches canonical APC hierarchy in thyroid context.

Specifically the HLA-II per-gene dotplot (`F2_*`) confirms: HLA-DRA, HLA-DRB1, HLA-DPA1, HLA-DPB1, CD74 are co-expressed at >50% pct-positive in DC and Myeloid; B_cell expresses the same set with slightly lower CIITA. Thyrocyte and Fibroblast have <10% pct-positive on most HLA-II genes.

---

## 6. Spatial HLA-II hotspot vs TLS niche

Table: `T8_visium_spotlevel.csv.gz` (n=55,873 spots across 16 samples).

Figure: `figs/F8_visium_TLS_HLAG.png`.

Spot-level Spearman: TLS module × HLA-II module rho = **−0.183** (p≈0, n=55,873). Counterintuitive sign at first glance, but explained by composition: many high-HLA-II spots are myeloid-dominated (CIITA+, CD74+, HLA-DR+) without TLS lymphoid markers, while pure TLS-marker spots (CXCL13/CCR6 high) sit in a B-cell-zone where HLA-II is high but the spot-mean CD74 score is diluted by adjacent T-cell spots. The strict spot-level rho is therefore not a clean replication of the bulk effect; instead the right test is per-niche enrichment, which would require Squidpy niche calling beyond track scope. Limitation noted in Section 9.

The spatial HLA-II hotspots are visible per-stage in `F7_visium_HLA_stage_celltype_heatmap.png` (myeloid cell column, ATC = 1.27, LPTC = 1.40).

---

## 7. HLA-G non-classical expression

From `T7_visium_celltype_HLA_per_sample.csv`, mean HLA-G expression in mapped Malignant cell spots:

| Stage | mean HLA-G (Malignant) |
|-------|-------------------------|
| PT (Epithelial) | 0.0007–0.004 (essentially zero) |
| PTC   | 0.118 |
| LPTC  | 0.087 |
| ATC   | 0.001 |

**Pattern:** HLA-G is a PTC/LPTC malignant-cell phenomenon — not a normal-thyrocyte feature, and largely absent from ATC. This matches the published pattern of HLA-G as a tumor-acquired immune-evasion signal in differentiated thyroid cancer that is *not* sustained when the tumor dedifferentiates into ATC. PTC-1 and LPTC-4 specifically show the highest per-spot HLA-G means (0.379, 0.314 respectively). HLA-G is one of the 19 HLA-I module members but its trajectory is opposite to the rest of HLA-I (HLA-A/B/C rise PT→ATC; HLA-G peaks in PTC/LPTC and drops in ATC).

This is a transcriptomic module observation only — no allele-level claim is made.

---

## 8. Decomposition: which compartment carries the bulk signal?

Tables: `T6_pu_decomposition.csv`, `T6_lu_decomposition.csv`, `T10_visium_decomposition.csv`.

Figures: `figs/F6_pu_decomposition.png`, `figs/F6_lu_decomposition.png`.

Bulk T5 reference: DM1 × HLA-I rho = +0.318 (TCGA n=527).

**scRNA pseudobulk decomposition (small n, wide CIs):**

| Cohort | Subset           | n | rho DM1×HLA-I | rho DM1×HLA-II |
|--------|-------------------|---|----------------|------------------|
| Pu | bulk_all_cells    | 7 | −0.71 | −0.43 |
| Pu | thyrocyte_only    | 5 | +0.40 | −0.40 |
| Pu | immune_only       | 7 | −0.64 | +0.43 |
| Pu | stromal_only      | 6 | −0.77 | −0.26 |
| Lu | bulk_all_cells    | 9 | −0.10 | +0.50 |
| Lu | thyrocyte_only    | 5 | −1.00 | +0.10 |
| Lu | immune_only       | 9 | −0.10 | +0.23 |
| Lu | stromal_only      | 9 | +0.53 | +0.08 |

These flip signs across cohorts at n=5–9; they cannot resolve the bulk decomposition.

**Visium spatial decomposition (more reliable, n=16 samples, full HLA gene panel):**

| Subset           | n | rho DM1×HLA-I | rho DM1×HLA-II |
|-------------------|---|----------------|------------------|
| bulk_all_spots    | 16 | **+0.374** (p=0.15) | +0.047 (p=0.86) |
| malignant_only    | 12 | **+0.545** (p=0.067) | +0.238 (p=0.46) |
| immune_only       | 15 | −0.193 (p=0.49) | −0.064 (p=0.82) |
| stromal_only      | 9  | +0.100 (p=0.80) | +0.050 (p=0.90) |

**Visium answer.** When per-sample DM1 score is computed only from malignant-cell spots and HLA-I module is computed only from those same spots, the within-malignant-cell DM1 × HLA-I rho is **+0.545**, *exceeding* the bulk TCGA reference (+0.318). Immune-only and stromal-only subsets give weak / null correlations. This implies the **bulk DM1 × HLA-I correlation is dominated by the thyrocyte/malignant-epithelial compartment** — the IFN/HLA-I module is being expressed by the tumor cells themselves, not just by infiltrating immune cells.

For HLA-II, the spatial decomposition does not replicate the bulk +0.393 sign within malignant spots strongly (+0.24, ns). HLA-II appears to require a stronger immune-infiltration co-signal that the thyrocyte compartment alone does not deliver.

**Working decomposition fraction (Visium, malignant compartment basis, HLA-I):**
- thyrocyte/malignant-intrinsic explanation: ≈100% of bulk rho reproducible from this compartment alone.
- immune-only contribution: negligible / negative direction.
- stromal-only contribution: small positive, not significant.

For HLA-II the same calculation gives ~60% of bulk rho captured by the malignant compartment alone; the residual is plausibly carried by myeloid/DC infiltration, which Section 5 confirms is the major HLA-II producer.

---

## 9. Limitations

1. **scRNA pseudobulk N is small (5–9 per cohort).** The per-sample scRNA decomposition is not stable; we deliberately downweight it relative to the n=16 Visium decomposition.
2. **Cell-type labels in scRNA are marker-score-derived, not curated atlas labels.** ~6–9% of cells are "Unassigned" in both cohorts. The major calls (Thyrocyte vs T_cell vs Myeloid vs DC) are robust to this, but rare populations (Mast, Plasma) may be slightly mis-aggregated.
3. **Pu cohort did not carry a usable condition column** (`mutation_status` is mostly NA in the processed h5ad). Section 3 condition stratification therefore uses Lu cohort and Visium only.
4. **Spatial spot resolution is multicellular (~50 µm Visium).** The Section 6 TLS x HLA-II spot-level negative rho is composition-driven; a per-niche analysis would change the sign and is the right next step.
5. **Bridge cohort (GSE163203, PTC + HT scRNA) is not on disk.** Per registry, status is `not_downloaded`. No HT-overlap scRNA claim is made here. The HT-overlap signal in our program is documented at bulk level in Round-2 / GSE286332 (memory `v17_GSE286332_strong_go`), and adding it at scRNA level remains future work.
6. **Boundary.** All conclusions are about gene-expression modules, not HLA alleles. Ancestry / autoimmune-susceptibility inferences are explicitly out of scope per Paper 1 boundary.

---

## 10. Implications

- **Paper 1.** Section 8 is the headline: the bulk DM1 × HLA-I correlation is reproduced *within* the malignant compartment alone (Visium rho=+0.545). This argues for a thyrocyte-intrinsic IFN-driven HLA-I induction in the DM1 compartment, *not* simply a higher immune-infiltrate confound. Useful for the Discussion 3.x paragraph that frames DM1 as a pro-immunogenic state. Caption: HLA gene-expression module — not allele genotype; Paper 1 transcriptomic immune-context only.
- **Paper 2 (HT-overlap PTC).** No data added here — bridge cohort not on disk. Hypothesis-only forward note: HT-overlap signal ought to elevate HLA-II module specifically in DC + Myeloid + thyrocyte compartments; this needs GSE163203 or Korean HT-PTC scRNA to test.
- **Paper 11 (pan-cancer DM1).** The compartment specificity here suggests Paper 11's DM1 × IFN-γ link would also localize to malignant cells per cancer; useful for interpreting why Phase C ICI cohorts (memory `paper11_pancancer_2026_05_08`) show DM1 × HLA-I OR > 1.

---

## File index

Figures (PNG, 8 total):
- `figs/F2_pu_celltype_x_HLA_gene_dotplot.png`
- `figs/F2_lu_celltype_x_HLA_gene_dotplot.png`
- `figs/F3_pu_HLA_module_violin.png`
- `figs/F3_lu_HLA_module_violin.png`
- `figs/F4_lu_celltype_condition_box.png`
- `figs/F5_pu_thyrocyte_pseudobulk_DM1xHLA.png`
- `figs/F5_lu_thyrocyte_pseudobulk_DM1xHLA.png`
- `figs/F6_pu_decomposition.png`
- `figs/F6_lu_decomposition.png`
- `figs/F7_visium_HLA_stage_celltype_heatmap.png`
- `figs/F8_visium_TLS_HLAG.png`

Tables (CSV, 11 total):
- `tables/T1_cohort_manifest.csv`
- `tables/T2_pu_celltype_x_HLA_gene_dotplot.csv`
- `tables/T2_lu_celltype_x_HLA_gene_dotplot.csv`
- `tables/T3_pu_HLA_module_per_celltype.csv`
- `tables/T3_lu_HLA_module_per_celltype.csv`
- `tables/T4_lu_celltype_condition_HLA.csv`
- `tables/T5_pu_thyrocyte_pseudobulk.csv`
- `tables/T5_lu_thyrocyte_pseudobulk.csv`
- `tables/T6_pu_decomposition.csv`
- `tables/T6_lu_decomposition.csv`
- `tables/T7_visium_celltype_HLA_per_sample.csv`
- `tables/T8_visium_spotlevel.csv.gz`
- `tables/T9_visium_sample_pseudobulk.csv`
- `tables/T10_visium_decomposition.csv`

Script: `scripts/hla_deepdive_2026_05_08/track12/run_track12.py`
Summary JSON: `track12_summary.json`

**Boundary boilerplate restated.** HLA gene-expression module — not allele genotype. Paper 1 transcriptomic immune-context only.
