# Paper 1 Deconvolution Rollup — 2026-05-09

**Purpose.** One-page operating map for the `p_deconv_2026_05_08` result directory. This is scaffolding/infra only: no voice-protected manuscript prose.

## Bottom Line

The deconvolution package is **paper-useful**, but it must be split by evidence tier:

| Tier | Components | Use |
|---|---|---|
| **Primary supplement** | v2, v3, v5A-E, v5D | Multi-method deconvolution, cross-cohort composition, driver/methylation/subcluster overlays, pseudotime, Pu full-reference robustness. |
| **Fig 8 mechanism support** | v9-v14 | MAPK/TDS-16/cross-cohort forest locks. This is the positive mechanism layer. |
| **Reviewer reserve** | v15A K2, v15C PRISM/DepMap | Generalizability/actionability reserve; not main mechanism proof. |
| **Caveat / stress test** | GSE250521 spatial v15B-v18 | Do not use as positive support. Use only to explain why Visium spatial MAPK × Panel anti-correlation is not a suitable validation. |
| **Rejected/deferred** | v1, scaden, Pu LinearSVR D2 | Retained for audit trail only. |

## Primary Deconvolution Story

| Version | Result | Disposition |
|---|---|---|
| v2 | Canonical labels + 4 methods. nu-SVR full residualization retains **47%** of the DM1 vs DM2 RAI effect; method bracket **24-70%**. | **Use.** Core purity/composition defense. |
| v3 | Lee/GSE213647 cross-cohort composition: **7/8** direction-consistent with TCGA; within-DM1 fusion+ vs fusion− all **|d| ≤ 0.42**. | **Use.** Cross-cohort + fusion-independence support. |
| v5A | Driver composition: RAS has more Epithelial and less Malignant than BRAF. | **Use as context.** |
| v5B | Mean 8-gene methylation beta tracks Myeloid **rho=+0.39**, Malignant **+0.30**, Epithelial **-0.31**, T cell **-0.42**. | **Use.** Methylation × compartment bridge. |
| v5C | DM1 sub-A vs sub-B: sub-A Malignant-rich **d=+1.22**, Epithelial-poor **d=-1.26**; immune compartments mostly NS. | **Use as boundary.** Not immune-hot/cold. |
| v5E | TCGA + Lee decile pseudotime: Malignant↓, Epithelial↑, Myeloid↓, Endothelial↑ from DM1→DM2 in both cohorts. | **Use.** Reproducible compositional trajectory. |
| v5D | Pu 2021 full-transcriptome NNLS reproduces Lu HVG informative axes **4/4**; Malignant/Epithelial magnitudes stronger. | **Use.** Refutes HVG-reference artifact. |

## Mechanism Lock

| Version | Result | Disposition |
|---|---|---|
| v9-v10 | MAPK output anti-correlates with thyroid differentiation in TCGA/Lee. | **Use.** Early mechanism layer. |
| v12 | Two-axis convergence: MAPK-active and HT-route can both reach panel silencing. | **Use.** Explains cohort heterogeneity. |
| v13 | Panel-8 tracks TDS-16: TCGA/Lee MAPK × Panel **-0.291/-0.395**, MAPK × TDS-16 **-0.306/-0.435**; ΔAUC TDS-16 vs Panel only **+0.007/+0.012**. | **Use.** Reviewer Q3 cherry-pick lock. |
| v14 | Five-cohort forest: pooled MAPK × Panel-8 **rho=-0.327 [−0.376, −0.278]**, n=1,287. Heterogeneity is predicted by HT-route and advanced saturation. | **Use.** Cross-cohort generalizability lock. |

## Reserve Evidence

| Component | Result | Boundary |
|---|---|---|
| v15A K2 | TCGA centered-panel OOF AUC **0.960**; K2 **246/260 DM2**, median p_DM2 **0.978**. | Score-distribution generalizability only; K2 lacks MAPK genes. |
| v15C PRISM | FDR<0.05 DM1-high selective hits: **7/11 MAPK-pathway**, hypergeom p **1.08e-14**. | Actionability reserve only; not thyroid-gene restoration. |
| v15C DepMap | Top dependencies: **MYC d=-0.499**, **NAMPT d=-0.445**. | Pan-cancer proxy; not thyroid-specific functional validation. |

## Spatial Caveat Closure

| Version | Result | Disposition |
|---|---|---|
| v15B | GSE250521 tumor epithelial MAPK × Panel raw **rho=+0.208**. Expected anti-correlation not observed. | **Do not use as support.** |
| v16 | Full covariate residualization attenuates raw positive signal by **83.1%**; ATC epithelial full residual **rho=-0.0007**. | Caveat/autopsy. |
| v17 | Full+spatial+detection residual pooled rho **+0.007**; observed raw median **+0.110** is generic vs random module median **+0.099**. | Strong reviewer-defense caveat. |
| v18 | KNN lag no rescue: raw KNN1-6 **+0.090**, residual KNN1-6 **+0.033**; anti-pockets not enriched (raw **0.88**, residual **1.01**). | Closure. Stop extending unless reviewer asks. |

## Recommended Integration

| Manuscript location | Include | Exclude |
|---|---|---|
| Main Results / Fig 8 support | v13/v14 factual cross-references only. | GSE250521 spatial positive/caveat details. |
| Supplementary deconvolution figure | v2 + v3 + v5D/E as the core figure set. | v1, scaden failed attempt, Pu LinearSVR D2. |
| Reviewer Q / reproducibility dossier | v15A K2, v15C PRISM/DepMap, v16-v18 caveat autopsy. | Any claim that GSE250521 spatial supports MAPK→Panel silencing. |
| Limitations / Q9 | Author keyboard only. | AI-generated prose. |

## Canonical Files

| File | Role |
|---|---|
| `SUMMARY.md` | Full chronological audit trail. |
| `DECONV_ROLLUP_2026_05_09.md` | This operating map. |
| `Fig_SX_deconvolution_v5_composite.{png,pdf}` | Core deconvolution composite candidate. |
| `Fig_SX_v13_TDS16_MAPK.{png,pdf}` | TDS-16/cherry-pick lock. |
| `Fig_SX_v14_cross_cohort_forest.{png,pdf}` | Cross-cohort MAPK × thyroid-score forest. |
| `Fig_SX_v18_spatial_lag_pockets.{png,pdf}` | Final spatial caveat closure. |
| `paper1_deconvolution_rollup_v18.html` | HTML dashboard version in papers hub. |

## Stop Rule

Do **not** extend the GSE250521 spatial branch again unless a reviewer specifically asks. The spatial branch is now closed as caveat evidence. Future deconvolution work should be limited to manuscript wiring, figure numbering, or methods reproducibility.
