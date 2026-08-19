# pMTnet CROSS-Neo-TCR Pilot

Fast pilot using local pMTnet on positive public TCR-pMHC rows versus same peptide-HLA shuffled-TCR decoys.

| Metric | Value |
|---|---:|
| model | pMTnet |
| pilot | positive_vs_same_pmhc_shuffled_tcr_decoy |
| n | 355 |
| n_pos | 190 |
| n_neg | 165 |
| auprc | 0.7486 |
| auroc | 0.7167 |
| positive_score_mean | 0.6973 |
| decoy_score_mean | 0.4608 |
| positive_rank_mean | 0.3026 |
| decoy_rank_mean | 0.5392 |

Interpretation: lower pMTnet rank means stronger predicted binding, so `pmtnet_score = 1 - rank` is used for AUPRC/AUROC.
This is still a pilot because shuffled decoys are synthetic and public-source leakage has not been fully eliminated.
