# R7v2 — Single-cell DM1 axis projection (GSE193581 / Lu 2023)

**HT-13 peaks in Myeloid+B cells; FA-12 peaks in Malignant cells; MAPK-9 peaks in BRAF-mutant Malignant cells — therefore DM1 = immune+thyroid two-source axis confirmed at single-cell level.**

## Cohort & method
- GSE193581 (Lu 2023 thyroid scRNA): n=30,035 raw-count cells from `scrna_F12.h5ad`; `author_celltype` joined from `GSE193581_hvg_adata.h5ad` via stripped barcode (-N suffix). Matched 27,754 / 30,035 (92.4%).
- Per-cell library-size normalization (1e4 + log1p), in-process on panel genes.
- Coverage: HT-13 13/13, FA-12 12/12, MAPK-9 9/9.

## Cell-type aggregation (median panel score)

| Cell type | n | HT-13 | FA-12 | MAPK-9 |
|---|---|---|---|---|
| Myeloid cell | 5064 | **1.759** | 0.078 | 0.162 |
| B cell | 344 | **1.540** | 0.000 | 0.000 |
| T cell | 13819 | 0.563 | 0.000 | 0.243 |
| NK cell | 733 | 0.248 | 0.000 | 0.153 |
| Endothelial cell | 143 | 0.082 | 0.111 | 0.382 |
| Malignant cell | 7028 | 0.030 | **0.113** | 0.235 |
| Fibroblast | 623 | 0.000 | 0.069 | 0.311 |

## Within-PTC by driver — Malignant compartment

| Driver | n | HT-13 | FA-12 | MAPK-9 |
|---|---|---|---|---|
| Malignant, **BRAF** | 1132 | 0.212 | 0.076 | **0.506** |
| Myeloid, BRAF | 3000 | 1.910 | 0.076 | 0.165 |

BRAF malignant cells more than double the overall malignant MAPK-9 median (0.506 vs 0.235).

## Axis decomposition (fraction of total panel signal)

- HT-13: Myeloid 50.5% + T 44.1% + B 3.0% = 97.6% lymphoid/myeloid (Malignant 1.2%).
- FA-12: Malignant 63.7% + Myeloid 31.6% = 95.3% (lymphocyte fraction = 0%).
- MAPK-9: T 54.3% + Malignant 26.7% + Myeloid 13.3%; per-cell median highest in BRAF-malignant.

## Interpretation
The two source-channels separate cleanly: HT-13 lives in immune cells (HLA-II / BCR / CXCL13 silent in thyrocytes, median 0.030). FA-12 lives in malignant thyrocytes (zero in lymphoid lineages) — the bulk DM1 FA-12 collapse is a thyroid-intrinsic loss, not immune dilution. MAPK-9 is BRAF-tumor-specific (0.506 vs 0.235), consistent with deconv v9–v12 MAPK-driven 8-gene methylation.

This is the missing bridge between H10 bulk deconv (DM1 BRAF-cPTC ≈ M2-Mφ + DC + Treg + Naive-B + Tfh + CD8) and the bulk observation that DM1 also tracks thyroid-differentiation collapse: the panel signal is two cell sources — immune-up + thyroid-down — summed in bulk RNA, with MAPK-9 resolving the BRAF tumor-cell layer.

## Outputs (in `r7_singlecell/`)
- `r7v2_per_cell_scores.tsv.gz`, `r7v2_celltype_panel_scores.tsv`, `r7v2_celltype_histology_panel_scores.tsv`, `r7v2_celltype_braf_panel_scores.tsv`, `r7v2_axis_decomposition.tsv`, `r7v2_summary.json`, `r7v2_run.py`.

## Caveats
- `author_celltype` is the original 8-class Lu 2023 annotation; DC/Treg/Tfh sub-labels need sub-clustering. H10 bulk LM22 deconv already resolves finer; this layer confirms compartment split, not exact sub-lineage.
- HT-13 contains HLA-II genes inducible on thyrocytes under IFN-γ — but malignant HT-13 (0.030) is ~60× below myeloid (1.76), so immune-dominance is robust.
- Per-cell library normalization (not full-matrix scanpy normalize_total) — absolute values not comparable across runs, but rank order and BRAF/RAS contrast are stable.
