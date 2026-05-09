# Wave8B Strict ESMFold Structure Summary

Strict set: n=89, positives=21. Filters: no in-house overlap, no TCR exact hit, no self exact hit, HLA pseudo-sequence available.

## Top Same-Row Results

| method | category | n / pos | AUROC | AUPRC | note |
|---|---|---:|---:|---:|---|
| Structure_LR | existing_method | 89 / 21 | 0.679 [0.548, 0.817] | 0.435 [0.260, 0.656] |  |
| MHCflurry | existing_method | 89 / 21 | 0.657 [0.472, 0.813] | 0.488 [0.299, 0.705] |  |
| ESMFold_univariate__min_pLDDT_peptide | structure_univariate | 89 / 21 | 0.644 [0.527, 0.765] | 0.313 [0.219, 0.493] | direction=inverse; raw_auroc=0.356 |
| Wave8_TCR_SelfSim_no_exact+ESMFold_3D_LR_repeated5x5_CV | combo_cv | 89 / 21 | 0.637 [0.501, 0.775] | 0.321 [0.206, 0.525] | strict-set internal CV; includes existing score as one feature |
| Wave8_TCR_SelfSim_full+ESMFold_3D_LR_repeated5x5_CV | combo_cv | 89 / 21 | 0.636 [0.498, 0.772] | 0.318 [0.205, 0.514] | strict-set internal CV; includes existing score as one feature |
| Wave8_TCR_SelfSim_full | existing_method | 89 / 21 | 0.632 [0.470, 0.785] | 0.403 [0.231, 0.649] |  |
| Wave8_TCR_SelfSim_no_exact | existing_method | 89 / 21 | 0.629 [0.461, 0.777] | 0.472 [0.272, 0.662] |  |
| Wave8_TCR_motif_only | existing_method | 89 / 21 | 0.625 [0.482, 0.779] | 0.462 [0.283, 0.673] |  |
| ESMFold_univariate__mean_pLDDT_peptide | structure_univariate | 89 / 21 | 0.609 [0.475, 0.749] | 0.310 [0.197, 0.511] | direction=inverse; raw_auroc=0.391 |
| Structure_LR+ESMFold_3D_LR_repeated5x5_CV | combo_cv | 89 / 21 | 0.598 [0.442, 0.752] | 0.357 [0.218, 0.586] | strict-set internal CV; includes existing score as one feature |
| BigMHC_IM | existing_method | 89 / 21 | 0.596 [0.431, 0.726] | 0.288 [0.185, 0.474] |  |
| ESMFold_3D_LR_repeated5x5_CV | structure_cv | 89 / 21 | 0.569 [0.427, 0.716] | 0.279 [0.181, 0.438] | strict-set internal CV |

## Interpretation Boundary

- The ESMFold branch is a pseudo-complex proxy: peptide + `GGGGS` + 34-aa HLA pseudo-sequence folded as one chain.
- It is useful as a fast structure-derived feature generator, not as a physical pMHC structure claim.
- Internal CV on n=89 is hypothesis-generating. The paper-grade claim still requires a leakage-controlled external/time-split benchmark.
