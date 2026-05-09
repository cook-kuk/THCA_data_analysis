# Pathology finding scout - conclusion brief

## Bottom line

No new TCGA-level image-DM1 main finding should be claimed yet. The stronger
pathology lead is narrower and more honest:

**In spatial thyroid H&E, UNI local morphology moderately reconstructs
spot-level DM1/RAI/TDS-like transcriptomic programs under leave-one-slide-out
validation, but the transferable TCGA DM1 label signal is weak. The visual
atlas suggests that the image signal may be driven by follicular/colloid and
stromal/collagenized tissue states rather than a clean tumor-cell DM1 subtype.**

This is a pathologist-review hypothesis, not a main-text claim.

## What was tested

1. Spatial-to-TCGA DM1 morphology transfer
   - Train PCA(32)+ridge on GSE250521 H&E UNI spot embeddings to predict
     spot-level `DM1_like_score`.
   - Evaluate by leave-one-slide-out inside GSE250521.
   - Project the learned morphology score to TCGA-THCA UNI WSI tiles.

2. UNI image score as a thyroid differentiation readout
   - Test whether center-held-out UNI LOTO `prob_DM1` tracks RNA-derived
     `dedifferentiation_proxy_score`.
   - Residualize against DM label, TSS, histology, molecular subtype, sex,
     and age.

3. Multi-axis spatial pathology scan
   - Test which spatial RNA programs are most visible in local H&E morphology:
     DM1, RAI, TDS, CAF/ECM, EMT, hypoxia, proliferation, epithelial score.

## Key results

| Test | Result | Interpretation |
|---|---:|---|
| Spatial DM1 morphology, leave-slide-out | rho = 0.302, p = 2.5e-68 | Moderate spatial H&E/RNA coupling |
| Slides with per-slide rho > 0.2 | 10/16 | Reproducible but not strong enough for headline |
| TCGA DM1 vs DM2 from spatial-trained morphology | AUC = 0.550 | No TCGA label transfer |
| TCGA DM1 vs DM2 after TSS residualization | AUC = 0.503 | Transfer signal dies |
| UNI LOTO vs dediff proxy, raw | rho = 0.527, p = 4.2e-05 | Strong raw explanatory association |
| UNI LOTO vs dediff proxy, after DM label | rho = 0.184, p = 0.18 | Not independent of DM label |
| UNI LOTO vs dediff proxy, full residualization | rho = 0.107, p = 0.44 | Not a standalone finding |
| Top spatial RNA program visible in H&E | DM1_like / RAI_8 | DM1 axis is the strongest local morphology program, but only moderate |

## Visual pathology lead

From `fig_spatial_dm1_tile_atlas.png`:

- Actual spatial `DM1_like_score` high tiles are often follicular/colloid-rich
  or low-cellularity thyroid tissue regions.
- Actual spatial `DM1_like_score` low tiles include more cellular tumor-like
  or ATC-like regions.
- H&E-predicted DM1 morphology high tiles are enriched for pink collagenized /
  desmoplastic stromal areas, especially in LPTC/PT fields.
- H&E-predicted DM1 morphology low tiles show more compact epithelial
  follicular/papillary architecture.

## Safe conclusion

The pathology lead is not "H&E perfectly predicts DM1." The safer lead is:

**DM1-associated H&E signal appears to be a local tissue-state signal involving
follicular/colloid and stromal-collagenized morphology. It is visible in spatial
H&E/RNA data, but it does not yet transfer cleanly to TCGA slide-level DM1
labels.**

## Next action

Send the tile atlas to a thyroid pathologist for blinded review:

1. Are predicted-DM1-high tiles genuinely collagenized/desmoplastic stroma?
2. Are transcriptomic-DM1-high tiles follicular/colloid-rich or low-cellularity
   benign-like regions?
3. Is the axis biologically plausible for dedifferentiation / thyroid
   differentiation loss, or is it mostly tissue-composition confounding?

If the pathologist can annotate a consistent morphology state, this becomes a
real pathology hypothesis for K2 H&E validation.
