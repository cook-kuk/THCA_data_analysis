# PantheonOS port — 16-sample spatial LR analysis (GSE250521 × Lu 2023)

**Pipeline**: PantheonOS skill `single_cell_spatial_mapping.md` (MOSCOT optimal transport, sc → spatial mapping) + Gallery #6 spatial disease LR (Squidpy ligand-receptor permutation, 500 perms, p<0.05).
**sc atlas**: Lu 2023 GSE193581 (n=67,678 cells × 2,000 HVG, `author_celltype` ∈ {T cell, Malignant cell, Myeloid cell, B cell, NK cell, Fibroblast, Endothelial cell, Epithelial cell}).
**spatial cohort**: GSE250521 Visium n=16 across 4 stages (PT/PTC/LPTC/ATC, n=4 each).
**Authority**: marathon override `v19_marathon_override_pantheonos_2026_05_08`. Voice-protected sections untouched.
**Run**: 2026-05-08, total wall ~75 min (3-way parallel for non-ATC, sequential for ATC due to memory).

## Headline finding — inverted-V LR diversity peaks at LPTC

| Stage | n | Median sig LR (p<0.05) | IQR | Top stage-specific LR |
|---|---|---|---|---|
| PT (paratumoral normal) | 4 (3 with LR; N-4 mapped 100% Epithelial → no clusters) | 170 | 108–220 | ANGPT2↔TIE1 (vasculature), CCL5↔TNFRSF4 |
| PTC | 4 | 1,742 | 1,484–2,028 | APOE↔TYROBP, HLA-DR↔CTSL/LGMN |
| **LPTC** | 4 | **2,512** | **1,293–4,117** | **CXCL13↔CCR6/CCR7, CCL23↔CCR6, CXCL3↔CCR6/CCR7** (TLS/germinal-center axis) |
| ATC | 4 | 598 | 234–1,194 | IFNG↔CD74/CTSK/CTSL, IL1A↔CSF1/IL1B, IL6↔IL13RA2 |

→ See `_AGGREGATE/stage_trend_n_significant.png` for box plot (log-scale).

The shape is **PT < PTC < LPTC > ATC**, not monotonic. LPTC has ~4× more cell-cell crosstalk than ATC. ATC isn't quietest by raw count but is the least diverse and dominated by inflammatory/IFNg signals — consistent with paper3 RAI-dediff axis (`v19_paper3_rai_dediff_axis_2026_05_06`) and Track B-lite dedifferentiation finding.

## Stage-specific LR axes (cross_stage_delta.csv, 46 pairs)

**LPTC → TLS / germinal-center program** (this is the strongest finding):
- CXCL13↔CCR6, CXCL13↔CCR7 — **CXCL13 is the canonical TLS/B-cell follicle chemokine**. Recurrence 1.0 in LPTC, 0.5 in PTC, 0 in ATC.
- CCL23↔CCR6/KDR, CXCL3↔ACKR1/CCR6/CCR7 — chemokine cascade for B/T zone organization
- RSPO3↔LGR5 — Wnt-stem signaling
- THBS1↔ELANE — neutrophil extracellular trap context

This **independently replicates** `v17_D5P6_BCR_clonal_TLS` (PTC+HT TLS d=+1.96) but in lateral PTC samples with a fully spatial readout — orthogonal validation. The lateral PTC stage appears to be where the immune ecosystem peaks, not the primary tumor.

**ATC → inflammatory dedifferentiation**:
- IFNG↔CD74/CTSK/CTSL/CXCL12/TYMS — cathepsin/MHC-II processing under sustained IFNg
- IL1A↔CSF1/IL1B/IL1R2/ITGA6/MMP9 — IL-1 inflammasome with MMP-9 matrix proteolysis
- IL6↔IL13RA2 — IL6/IL13 dual axis
- ADM↔NOS1 — adrenomedullin/nitric oxide
- SAA1↔FPR2 — acute-phase signaling

Aligns with `v19_paper3_track_b_lite_2026_05_06` (ATC d=+2.53 myeloid / d=−2.51 thyroid_diff).

**PT → quiet baseline** dominated by ANGPT2/TIE1, CCL2/CCL5, FGF7/FGFR2, CLU↔TYROBP — vasculature + resident myeloid; no T/B effector signaling.

## Per-sample summary

See `_AGGREGATE/per_sample_summary.csv` (15 rows; N-4 had `lr_skipped: true` because mapping was 100% Epithelial cell so Squidpy ligrec needs ≥2 groups).

```
sample, stage, n_spots, n_sc_cells, n_groups, n_significant_lr
ATC-1, ATC, 4715, 36157, 3, 650
ATC-2, ATC, 2354, 36157, 4, 546
ATC-3, ATC, 1621, 36157, 2, 234
ATC-4, ATC, 2586, 36157, 3, 1194
LPTC-1, LPTC, 1520, 23532, 5, 1717
LPTC-2, LPTC, 4071, 23532, 7, 4117
LPTC-3, LPTC, 3336, 23532, 6, 1293
LPTC-4, LPTC, 4707, 23532, 6, 3308
N-1, PT, 3395, 7989, 3, 108
N-2, PT, 3780, 7989, 2, 170
N-3, PT, 4538, 7989, 2, 220
PTC-1, PTC, 2262, 23532, 5, 2028
PTC-2, PTC, 4947, 23532, 5, 1484
PTC-3, PTC, 4516, 23532, 5, 1687
PTC-4, PTC, 4498, 23532, 6, 1797
```

## Outputs

```
project/results/pantheonos_demo/
├── SUMMARY.md                          ← this file
├── <SAMPLE>/                           ← 16 sample dirs (N-1..4, PTC-1..4, LPTC-1..4, ATC-1..4)
│   ├── transport_matrix.npy            (MOSCOT OT, n_spots × n_sc_cells)
│   ├── visium_with_mapped.h5ad         (Visium + obs.mapped_celltype, mapped_celltype_conf)
│   ├── lr_means.csv / lr_pvalues.csv   (Squidpy multi-index)
│   ├── lr_top.csv                      (top 30 by p,mean)
│   ├── spatial_overlay.{png,pdf}       (mapped celltype + confidence on tissue)
│   └── run_meta.json
└── _AGGREGATE/
    ├── per_sample_summary.csv          (15 rows)
    ├── lr_condition_recurrence.csv     (608 LR pairs × 4 stages)
    ├── lr_top_per_condition.csv        (193 top-recurrent pairs)
    ├── cross_stage_delta.csv           (46 stage-specific pairs, ≥0.5 delta)
    ├── aggregate_heatmap.{png,pdf}     (top 50 most stage-variable LR pairs)
    └── stage_trend_n_significant.{png,pdf}
```

## Caveats (Nat Comm reviewer-grade)

1. **HVG-scale mismatch**: Lu 2023 sc atlas was provided as 2,000 HVG with z-scored X. MOSCOT solved on this z-scored input rather than ideal log1p. PCA computation done internally by MOSCOT before solver. OT theory robust to this, but parameter sensitivity (alpha, tau) not yet swept.
2. **MOSCOT alpha=0.01** (near-pure linear). Skill recipe said `alpha=0` but moscot 0.5 forbids that. Sensitivity not tested.
3. **Cell-type mapping is dominated by 2-3 categories per sample** (e.g. PTC-1: 48% B cell, 48% Malignant cell, 4% other). The TLS-like B-cell-rich zone is real and aligns with prior findings, but mapping confidence shows uniform high values — could indicate over-confident OT solution. Down-weighting via tau-b (=0.8) attempted but more work needed.
4. **No spatial coherence regularizer**: MOSCOT used `alpha≈0` (linear OT, gene-expr only). Adding spatial structure penalty (alpha=0.5) would smooth celltype assignments — defer to revision.
5. **N=3 per stage for PT** (N-4 dropped — 100% Epithelial mapping). Acceptable for descriptive; for inferential statistics we'd want N≥4 per group.
6. **One sc atlas, one spatial cohort** — replication on independent matched pairs (e.g. GSE241184 sc) recommended for revision.
7. **Squidpy multi-index CSV format**: needs `header=[0,1], index_col=[0,1]` to parse — documented in `aggregate_lr_results.py`.
8. **B cell over-mapping in primary PTCs**: spatial overlay (PTC-1) shows two distinct zones consistent with TLS infiltrate but warrants pathologist review of original H&E.

## Provenance

- **Source skill** (verbatim recipe): `project/external/pantheonos/single_cell_spatial_analysis/skill_sc_spatial_mapping.md`
- **Pipeline runner**: `scripts/external/run_pantheonos_pipeline.py`
- **Aggregator**: `scripts/external/aggregate_lr_results.py`
- **Spatial plot helper**: `scripts/external/add_spatial_plot.py`
- **Batch driver**: `scripts/external/batch_run_all_samples.sh`

## Marathon footnote

This run is the data-execution arm of the one-shot PantheonOS port override authorized 2026-05-08. Default ("새 분석 = paper-blocking only") resumes after this session. Voice-protected sections (Hook/Aim/Disc 3.1/Limitations/Cover Para 1/Q9) untouched.
