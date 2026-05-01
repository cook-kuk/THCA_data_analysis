# GSE241184 — single-cell DM-axis projection

**Author:** Seungho Cook  •  **Date:** 2026-04-28  •  **Pipeline:** v17p35
**Inputs:** `/data/thca/v17_gse241184/`  •  **Outputs:** `results/v17_gse241184/`, `reports/html/figs_interactive/v17/v17_gse241184_*.html`  •  **Log:** `logs/v17_gse241184.log`

## 0. Cohort mismatch (important)
The task brief described GSE241184 as the *"largest public Chinese thyroid scRNA cohort, 50 tumor + 14 normal"*. The GEO record (PubMed 38061122, Chen W & Zhong S, Affiliated Cancer Hospital of Nanjing Medical University, Aug 2023) actually contains **3 samples from a single 17-year-old female PTC patient**: thyroid tumor (TT), adjacent normal thyroid (NT), and lymph-node metastasis (LN). We proceeded with single-cell DM-axis projection using *per-tissue* statistics; per-patient dominance and immune-DM1 Spearman tests are reported descriptively only (n=1 patient). If a 50+14 Chinese atlas is desired, candidate datasets to investigate next: GSE184362, GSE193581 (Lu 2023), GSE253560, or HRA-CN equivalents.

## 1. Download & QC
- Tarball `GSE241184_RAW.tar` = **204 MB** (10x Genomics v3.1, Cell Ranger 5.0.0, hg38).
- Raw barcode totals: TT 9 864, NT 7 991, LN 13 175 = **31 030**.
- After QC (mito < 20%, n_genes 200–6 000) and gene filter (≥ 3 cells): **29 939 cells × 26 448 genes**, all stored as `csr_matrix`.

## 2. Integration & annotation
- normalize_total(1e4) → log1p → HVG 2 000 (per-batch) → PCA 50 → **Harmony** integration on `sample` (10 iters, stopped before convergence — acceptable for 3 batches) → neighbours (k=15) → UMAP → Leiden res=0.6.
- Cell types from marker-score `argmax` (Thyrocyte rule: `score_Thyrocyte > 0.25` overrides generic Epithelial). Counts: Thyrocyte **6 632**, T_cell 11 044, Myeloid 1 800, B_cell 1 097, NK_cell 622, Endothelial / Fibroblast / Epithelial fill the rest.

## 3. 8-gene panel projection
All 8 genes present in the matrix: **SLC5A5, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, DIO1** (zero missing — clean recovery on the Chinese cohort). Per-cell `DM_score` = mean Z(log-norm) across the 8 genes. DM1/DM2 = median split on the 6 632 Thyrocyte subset (median = **0.926**).

## 4. Per-sample (per-tissue) DM dominance — Thyrocytes only
| Sample | Tissue   | n_cells | n_thyrocyte | DM1 | DM2 | **DM1 %** | DM_score median |
|--------|----------|---------|-------------|-----|-----|-----------|------------------|
| NT     | Normal   | 7 769   | 1 134       | 807 | 327 | **71.2 %** | 1.309 |
| TT     | Tumor    | 9 637   | 2 100       | 1 185 | 915 | **56.4 %** | 1.008 |
| LN     | LN_Met   | 12 533  | 3 398       | 1 324 | 2 074 | **39.0 %** | 0.817 |

**Monotonic Normal → Tumor → LN_Met gradient**, exactly the direction predicted by the v17p35 differentiation hypothesis.

## 5. Statistics (Thyrocyte DM_score)
- Kruskal-Wallis across NT/TT/LN: **H = 571.7, p = 7.1 × 10⁻¹²⁵**.
- Mann-Whitney NT vs (TT+LN), one-sided greater: **p = 4.4 × 10⁻⁸⁴**.
- NT vs LN: **p = 2.5 × 10⁻¹⁰⁶**;  TT vs LN: **p = 1.8 × 10⁻⁴⁸**.
- Patient-level dominance (n=1): NT reaches the **70 % threshold (71.2 %)**; TT and LN do not. Reported descriptively.
- Spearman DM1 % vs immune %: rho = 0.50, p = 0.67 (n=3, **uninterpretable**).

## 6. Figures (Plotly, dark-bg `#0b0e12`)
1. `reports/html/figs_interactive/v17/v17_gse241184_umap_tissue.html` — Harmony-integrated UMAP coloured by tissue.
2. `reports/html/figs_interactive/v17/v17_gse241184_umap_dmscore.html` — UMAP coloured by per-cell DM_score (RdBu_r).
3. `reports/html/figs_interactive/v17/v17_gse241184_umap_celltype.html` — UMAP coloured by marker-score cell type.
4. `reports/html/figs_interactive/v17/v17_gse241184_dm1_per_sample.html` — DM1 % bar with 50 / 70 % reference lines.
5. `reports/html/figs_interactive/v17/v17_gse241184_immune_vs_dm1.html` — immune fraction vs DM1 % (n=3, descriptive).
6. `reports/html/figs_interactive/v17/v17_gse241184_dmscore_violin.html` — Thyrocyte DM_score violins per tissue with KW & MWU annotations.

## 7. Result files
- `results/v17_gse241184/summary.json` — machine-readable summary including all p-values & cohort caveat.
- `results/v17_gse241184/per_sample_summary.tsv` — table reproduced above.
- `results/v17_gse241184/per_cell_metadata.tsv.gz` — 29 939 cells × {sample, tissue, leiden, celltype, DM_score, DM_class, QC}.

## 8. Conclusions
- The 8-gene v17p35 DM panel projects cleanly to single-cell resolution on a Chinese PTC sample with **all 8 genes detected**.
- A monotonic decrease in mean DM_score and DM1 % across **Normal → primary Tumor → Lymph-node Metastasis** is observed in Thyrocytes (KW p = 7e-125), consistent with the dedifferentiation gradient seen in bulk Korean GSE213647 (Kruskal p = 8.7e-55).
- The **n=1 patient** scope means per-patient dominance and immune-DM1 correlation cannot be tested inferentially; a true Chinese multi-patient atlas (e.g. GSE184362 Pu 2021) is needed to fulfil the original brief.

## 9. Blockers
- **Cohort mismatch** (only 3 samples / 1 patient instead of 50 + 14). Not a runtime failure — flagged here so downstream v17p35 talk slides do not over-claim.
- Harmony stopped before convergence at iter 10 — expected with only 3 batches; not a problem.
- Disk: 204 MB tarball + ~1 GB working files; well within `/data` headroom (441 GB free). RAM peak < 4 GB of 31 GB available.

## 10. Reproducibility
- Script: `notebooks_or_scripts/v17_GSE241184_scrna.py`
- random_state = 42 throughout (numpy, scanpy, harmonypy, leiden, UMAP).
- Run: `cd /opt/thyroid-dash/project && . .venv/bin/activate && python notebooks_or_scripts/v17_GSE241184_scrna.py`
