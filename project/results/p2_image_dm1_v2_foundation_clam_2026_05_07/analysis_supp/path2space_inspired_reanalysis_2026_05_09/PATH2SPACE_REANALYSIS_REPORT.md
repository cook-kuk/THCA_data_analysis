# Path2Space-inspired spatial pathology reanalysis

## Strategy Ported From The Cell 2026 Paper

- Keep the useful pattern, not the breast-specific claim: H&E embeddings are used as a low-cost surrogate for local spatial RNA programs.
- Replace raw spot labels with KNN-8 spatial-smoothed labels before model evaluation, matching the paper's noise-reduction logic.
- Evaluate by leave-one-slide-out validation; no same-slide leakage.
- Move beyond direct MAPK x Panel spot correlation by adding SPAND-like heterogeneity and spatial-cluster composition.
- Treat this as GSE250521/Paper-2 support or caveat, not as causal Paper-1 mechanism prose.

## Headline Results

- Top smoothed target: **DM1_like_score**, LOSO Spearman rho = **0.392**.
- DM1/RAI axis: raw rho = 0.295; smoothed rho = **0.392**; gain = +0.097.
- MAPK output axis from H&E: smoothed rho = **0.298**.
- Predicted spatial clusters: best k = **3**, silhouette = 0.272, dominant-cluster/stage ARI = 0.000.

## Ranked Smoothed Targets

| rank | target | rho | p | median slide rho | slides rho>0.2 |
|---:|---|---:|---:|---:|---:|
| 1 | DM1_like_score | 0.392 | 4.73e-118 | 0.373 | 12/16 |
| 2 | TDS_like_score | 0.381 | 6.43e-111 | 0.338 | 12/16 |
| 3 | CAF_ECM_score | 0.318 | 2.63e-76 | 0.289 | 9/16 |
| 4 | MAPK_output_score | 0.298 | 9.00e-67 | 0.307 | 12/16 |
| 5 | Epithelial_score | 0.244 | 1.76e-44 | 0.274 | 12/16 |
| 6 | Hypoxia_score | 0.188 | 9.72e-27 | 0.289 | 10/16 |
| 7 | Proliferation_score | 0.108 | 8.64e-10 | 0.146 | 6/16 |
| 8 | EMT_score | 0.092 | 1.59e-07 | 0.101 | 5/16 |

## Predicted SPAND-Like Heterogeneity

Values are `1 - Moran's I` on rank-normalized predicted scores, so higher values mean more local intermixing / less spatial autocorrelation.

| axis              |    PT |   PTC |   LPTC |   ATC |
|:------------------|------:|------:|-------:|------:|
| DM1_like_score    | 0.45  | 0.313 |  0.278 | 0.506 |
| TDS_like_score    | 0.426 | 0.336 |  0.282 | 0.519 |
| MAPK_output_score | 0.44  | 0.326 |  0.269 | 0.331 |
| Epithelial_score  | 0.356 | 0.329 |  0.267 | 0.362 |
| CAF_ECM_score     | 0.603 | 0.305 |  0.302 | 0.541 |

## Boundary

This reanalysis improves the pathology-to-spatial-RNA framing but does not rescue the direct same-spot MAPK x Panel anti-correlation closed in v15-v18. The safe use is: spatial H&E carries a weak-to-moderate DM1/RAI signal and can support Paper 2 image-to-DM1 strategy; GSE250521 remains caveat-only for Paper 1 MAPK-silencing mechanism.
