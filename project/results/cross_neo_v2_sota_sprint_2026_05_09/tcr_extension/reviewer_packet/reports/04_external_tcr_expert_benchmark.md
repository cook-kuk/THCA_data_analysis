# External TCR Expert Challenge Benchmark

Fast multi-panel benchmark using positive public TCR-pMHC rows and same-panel shuffled-TCR decoys.

Claim boundary: this is a model-selection diagnostic. Decoys are synthetic and public-source overlap is possible; do not use as external SOTA evidence.

## Pooled Result

| Model | n | AUPRC | AUROC | Positive mean | Negative mean |
|---|---:|---:|---:|---:|---:|
| mean_pMTnet_TEPCAM | 1191 | 0.7841 | 0.7652 | 0.6920 | 0.5108 |
| logreg_groupcv_pMTnet_TEPCAM | 1177 | 0.7803 | 0.7683 | 0.5970 | 0.4029 |
| pMTnet | 1177 | 0.7479 | 0.7255 | 0.7540 | 0.5413 |
| TEPCAM | 1191 | 0.6878 | 0.6970 | 0.6321 | 0.4812 |

## Panel Construction

| Panel | Eligible positives | Sampled positives | Kept decoys | Dropped decoys | Final rows |
|---|---:|---:|---:|---:|---:|
| all_public | 65815 | 100 | 100 | 0 | 200 |
| cancer_context | 2034 | 100 | 100 | 0 | 200 |
| paired_tcr | 42297 | 100 | 100 | 0 | 200 |
| vdjdb | 38989 | 100 | 100 | 0 | 200 |
| mcpas | 9334 | 100 | 100 | 0 | 200 |
| tenx | 17348 | 100 | 91 | 9 | 191 |

## Panel-Level Best Raw Expert

| Panel | Best raw expert | AUPRC | AUROC |
|---|---|---:|---:|
| all_public | TEPCAM | 0.6811 | 0.6584 |
| cancer_context | pMTnet | 0.7914 | 0.8073 |
| mcpas | pMTnet | 0.8402 | 0.8759 |
| paired_tcr | TEPCAM | 0.7077 | 0.7116 |
| pooled | pMTnet | 0.7479 | 0.7255 |
| tenx | pMTnet | 0.8488 | 0.8610 |
| vdjdb | TEPCAM | 0.7087 | 0.7220 |

## Source-Heldout Ensemble

| Heldout source | n | AUPRC | AUROC |
|---|---:|---:|---:|
| VDJdb | 538 | 0.6523 | 0.6040 |
| McPAS-TCR | 448 | 0.8025 | 0.8232 |
| VDJdb_10x_chunk | 191 | 0.8408 | 0.8340 |

## Decision

- Promote pMTnet as the first HLA-aware optional TCR expert.
- Promote TEPCAM as a second sequence-only optional TCR expert.
- Use the mean pMTnet/TEPCAM score as the first no-training diagnostic ensemble; only use logistic ensembles inside strict folds.
- Keep all outputs out of the main pMHC model unless paired TCR data are present.
