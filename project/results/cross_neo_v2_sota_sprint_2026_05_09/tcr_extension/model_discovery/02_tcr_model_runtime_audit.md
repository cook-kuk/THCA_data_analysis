# CROSS-Neo-TCR Runtime Model Audit

Date: 2026-05-09

## Fast Accurate Decision

The next TCR modeling layer should move from handcrafted sequence features to external TCR-specific predictors.

Immediate runnable models:

| Model | Runtime status | Input fit to current registry | Result |
|---|---|---|---|
| NetTCR-2.2 | Local repo + pretrained TFLite weights run successfully | Needs CDR1/CDR2 plus CDR3 alpha/beta; current registry mostly lacks CDR1/CDR2 | Repository sanity run AUPRC 0.932, AUROC 0.980; runtime check only |
| pMTnet | Local repo + trained library run successfully after Keras 3 `compile=False` load patch | Excellent fit: CDR3 beta + peptide + HLA | CROSS-Neo-TCR shuffled-decoy pilot AUPRC 0.749, AUROC 0.717 |
| TEPCAM | Local repo + checkpoint run successfully after HMMER/ANARCI runtime fix and PyTorch load patch | Good fit: CDR3 beta + peptide, no HLA | Same pilot AUPRC 0.753, AUROC 0.716 |
| ERGO-II | Local repo sanity and same-pilot adapter run successful after Lightning import/decorator and autoencoder-path patch | Good fit when CDR3 alpha/beta, peptide, MHC, V/J exist; beta-only/unknown-VJ adapter is weak | Same pilot beta-only/unknown-VJ AUPRC 0.559, AUROC 0.487 |
| PanPep | Local repo zero-shot run successful after CPU deserialization patch | Good fit: CDR3 beta + peptide, no HLA | Same pilot weak: AUPRC 0.525, AUROC 0.485 |
| TSpred | Local repo + pretrained example weights run successfully after CPU fallback patch | Needs CDR1/CDR2 plus paired alpha/beta | Runtime sanity successful; needs CDR1/CDR2 reconstruction before CROSS-Neo use |

Still blocked for broad CROSS-Neo use:

| Model | Blocker | Next action |
|---|---|---|
| NetTCR-2.2 | Requires CDR1/CDR2 plus CDR3 features | Reconstruct CDR1/CDR2 from V genes where possible; otherwise keep sanity-only |
| TSpred | Requires CDR1/CDR2 and paired alpha/beta context | Same CDR1/CDR2 reconstruction gate as NetTCR |
| Structure predictors | Full TCR chains and MHC sequences missing for most registry rows | Recover full chains/templates before making structure claims |

## pMTnet Pilot Result

Local pMTnet was run on positive public TCR-pMHC rows versus same peptide-HLA shuffled-TCR decoys.

| Metric | Value |
|---|---:|
| Evaluated pairs with valid pMTnet rank | 355 |
| Positives | 190 |
| Decoys | 165 |
| AUPRC | 0.7486 |
| AUROC | 0.7167 |
| Mean positive score | 0.6973 |
| Mean decoy score | 0.4608 |
| Mean positive rank | 0.3026 |
| Mean decoy rank | 0.5392 |

Interpretation: lower pMTnet rank is stronger binding, so the pilot used `1 - rank` as the model score.

## Same-Pilot External Model Comparison

| Model | Inputs | n | AUPRC | AUROC | Positive mean | Decoy mean | Read |
|---|---|---:|---:|---:|---:|---:|---|
| pMTnet | TCRbeta + peptide + HLA | 355 | 0.7486 | 0.7167 | 0.6973 | 0.4608 | P0 external TCR expert; HLA-aware |
| TEPCAM | TCRbeta + peptide | 374 | 0.7532 | 0.7156 | 0.6030 | 0.4339 | P1 sequence-only expert; surprisingly competitive on this pilot |
| ERGO-II-vdjdb | TCRbeta + peptide + coarse MHC | 374 | 0.5587 | 0.4870 | 0.3620 | 0.3821 | Runtime works, but beta-only/unknown-VJ adapter is weak |
| PanPep | TCRbeta + peptide | 374 | 0.5248 | 0.4855 | 0.4455 | 0.4545 | Diagnostic only for now |

Pairwise score correlations are low (`pMTnet` vs `TEPCAM` Spearman rho 0.241), so pMTnet and TEPCAM may provide complementary signals rather than redundant scores.

## Fast Multi-Panel Challenge Benchmark

We then reran pMTnet and TEPCAM on six 100-positive challenge panels with same-panel shuffled-TCR decoys: all public, cancer-context, paired-TCR, VDJdb, McPAS-TCR, and 10x-derived rows.

| Model | n | AUPRC | AUROC | Positive mean | Negative mean |
|---|---:|---:|---:|---:|---:|
| mean pMTnet + TEPCAM | 1,191 | 0.7841 | 0.7652 | 0.6920 | 0.5108 |
| peptide-HLA GroupKFold logistic pMTnet + TEPCAM | 1,177 | 0.7803 | 0.7683 | 0.5970 | 0.4029 |
| pMTnet | 1,177 | 0.7479 | 0.7255 | 0.7540 | 0.5413 |
| TEPCAM | 1,191 | 0.6878 | 0.6970 | 0.6321 | 0.4812 |

Source-heldout logistic ensemble was strongest for McPAS-TCR and 10x-derived heldouts (AUPRC 0.8025 and 0.8408), but weaker for VDJdb heldout (AUPRC 0.6523, AUROC 0.6040). This reinforces the source-shift boundary: useful diagnostic signal, not a universal TCR-recognition claim.

## Claim Boundary

This is a real improvement over the handcrafted shuffled-decoy pilot, which was approximately random at AUPRC/AUROC 0.50. However, the external-model pilot still uses synthetic shuffled decoys and public-source data, so it supports **model-priority escalation**, not a final SOTA claim.

## Next Order

1. Make pMTnet the first external TCR expert score for CROSS-Neo-TCR.
2. Add TEPCAM as a second sequence-only expert and use mean pMTnet/TEPCAM as the first no-training diagnostic ensemble.
3. Restrict ERGO-II to rows with richer V/J and alpha/beta fields; the beta-only/unknown-VJ adapter should not be promoted.
4. Keep PanPep as diagnostic/negative-control unless it improves under a better adapter or few-shot setting.
5. Build CDR1/CDR2 reconstruction from V genes to unlock NetTCR-2.2 and TSpred on CROSS-Neo rows.
6. Keep structure models behind full-chain TCR/MHC sequence recovery.
