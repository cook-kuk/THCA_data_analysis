# K2 full-transcriptome pilot verdict

Build: 2026-05-09 KST.

## Verdict

**Not a breakthrough result. Use as reviewer-reserve / negative-control only.**

The local kallisto run completed successfully for all available FASTQs, but the
available local data are only **n=9**: **5 tumors and 4 matched normals**. All
9 samples inherit the prior 8-gene mini-index call **DM2**, so the completed
full-transcriptome pilot cannot test DM1-vs-DM2 biology in K2.

## What succeeded

- Quantification completed: 9/9 local FASTQ runs.
- Gene-level matrix: 61,228 gene symbols.
- Panel coverage: HT-13 13/13,
  FA-12 12/12, MAPK-9
  9/9.
- kallisto index: `/data/thca/reference_kallisto/gencode.v44.kallisto.idx.NEW`.

## What it says

Tumor-vs-normal signal is weak and not deployable at n=9:

| Panel | AUC tumor vs normal | MW p |
|---|---:|---:|
| HT_13 | 0.400 | 0.730 |
| FA_12 | 0.750 | 0.286 |
| MAPK_9 | 0.750 | 0.286 |

Paired tumor-minus-normal medians across the 4 complete pairs:

| Panel | n pairs | median tumor-normal | Wilcoxon p |
|---|---:|---:|---:|
| HT_13 | 4 | -0.158 | 0.875 |
| FA_12 | 4 | 0.224 | 0.125 |
| MAPK_9 | 4 | 0.284 | 0.375 |

Spearman correlations across the 9 samples are descriptive only:
HT-13 vs MAPK-9 rho=-0.483; FA-12 vs MAPK-9 rho=0.817.
The signs are compatible with the broader two-axis story, but n=9 and all-DM2
labels make this non-deployable as evidence.

## Disposition

- Do **not** call this "K2 full cohort validation".
- Do **not** cite it as Korean DM1-vs-DM2 replication.
- Safe use: one reviewer-reserve sentence saying that the partial local K2
  re-quant was technically successful but underpowered/all-DM2 and therefore
  did not alter the Paper 1/2 decision.

## Files

- `k2_full_panel_scores.tsv`
- `k2_panel_aucs.tsv`
- `k2_paired_tumor_normal_deltas.tsv`
- `k2_paired_tumor_normal_summary.tsv`
- `k2_panel_spearman.tsv`
- `fig_k2_pilot_panel_scores.png`
