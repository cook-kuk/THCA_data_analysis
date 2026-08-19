# External TCR Model Pilot Comparison

Same positive public TCR-pMHC rows versus same peptide-HLA shuffled-TCR decoys used for the pMTnet pilot.

| Model | Inputs | n | AUPRC | AUROC | Positive mean | Decoy mean | Claim boundary |
|---|---|---:|---:|---:|---:|---:|---|
| pMTnet | TCRbeta + peptide + HLA | 355 | 0.7486 | 0.7167 | 0.6973 | 0.4608 | P0 runnable; pilot claim-limited |
| PanPep | TCRbeta + peptide | 374 | 0.5248 | 0.4855 | 0.4455 | 0.4545 | P1 runnable; no HLA in score |
| TEPCAM | TCRbeta + peptide | 374 | 0.7532 | 0.7156 | 0.6030 | 0.4339 | P1 runnable; no HLA in score |
| ERGO-II-vdjdb | TCRbeta + peptide + coarse MHC | 374 | 0.5587 | 0.4870 | 0.3620 | 0.3821 | P1 runnable; adapter is beta-only/unknown-VJ |

Interpretation boundary: these are fast pilot numbers on synthetic shuffled decoys, not a clean external benchmark.
pMTnet has an input advantage because it uses HLA; TEPCAM and PanPep are TCRbeta-peptide models, and the ERGO-II pilot used only coarse MHC with unknown V/J and no alpha chain.
