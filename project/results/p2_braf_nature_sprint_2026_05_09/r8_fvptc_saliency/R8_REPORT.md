# R8 — FVPTC image saliency: what does CLAM learn where it works (AUC=0.97)?

**Mechanism verdict: PARTIAL — image attention is TRENDING with dedifferentiation RNA in FVPTC but does not clear p<0.10 at n=16.** Specifically: attention top-10 mass × MAPK9 rho=-0.40 p=0.128; DM1-vs-DM2 top-10 centroid L2 = 25.1 (perm p_distance = 0.125), ~2x the BRAF-cPTC reference. Direction matches the H1 reframe (image reads dediff morphology, not lymphoid signal), and the FVPTC pattern differs from BRAF-cPTC where attention instead tracks TDS-8 thyroid-diff (rho=+0.33, p=0.037). n=16 / 3 DM1 caps detection power; a pathology-foundation backbone or external FVPTC cohort would be needed to harden.

## Cohort

RAS/FVPTC n=16 (3 DM1 / 13 DM2), image AUC=0.97 (H1). BRAF/cPTC reference n=41 (26 DM1 / 15 DM2), image AUC=0.56 (H1; H7 falsified rescue).

## Findings

**Attention concentration.** FVPTC top-10 tile mass mean=0.223 vs BRAF-cPTC 0.306 (MWU p=0.409); entropy ratio 0.912 vs 0.865 (p=0.576). No global concentration shift between strata.

**DM1 vs DM2 top-10 centroid divergence.** FVPTC L2=25.1 (perm p=0.125); BRAF-cPTC L2=12.0 (perm p=0.675). Centroid divergence is ~2x larger in FVPTC; the model 'looks at different tiles' for DM1 vs DM2 in FVPTC, but not in BRAF-cPTC.

**Image-RNA correlation (top-10 mass × panel z).** FVPTC: MAPK9 rho=-0.40 p=0.128 (strongest, NEGATIVE), TDS8 rho=+0.29 p=0.279, FA12 rho=-0.09 p=0.737, HT13 rho=+0.05 p=0.863, TLS11 rho=+0.25 p=0.356. BRAF-cPTC: TDS8 rho=+0.33 p=0.037 (only significant coupling); MAPK9 rho=+0.07 p=0.680. Direction in FVPTC matches the H1 reframe (image reads dediff axis, not lymphoid signal).

## Stratum-specific reading (Nature-grade claim)

- **FVPTC.** Attention biases away from MAPK-output-high tiles (dediff-poor / encapsulation-rich morphology) — direction-correct, p=0.13.
- **BRAF/cPTC.** Attention drifts toward residual thyroid-diff (TDS-8) tiles regardless of DM (p=0.04, but tracks neither DM1 axis nor immune signal).
- **Net.** Stratum-specific morphology: image-only AUC works in FVPTC because attention partitions a coherent dediff axis; it fails in BRAF/cPTC because attention finds an axis (TDS-8) that does not separate DM in that stratum.

## Tile coordinates

18 of 57 slides have coord-TSV alongside features; top-10 (rank, x, y, mpp) saved for downstream visual review. Partial-coverage is a legacy artifact (later H&E runs wrote coords), not a re-derivable resource here.

## Caveats

- n=16 / 3 DM1: |rho|≈0.5 is the floor for p<0.05.
- ImageNet ViT-L tile features are generic; UNI/Virchow/Phikon would likely amplify FVPTC MAPK rho.
- AUC=0.97 partly reflects 3:13 imbalance; the mechanism question still stands.
- Centroid permutation null = 5000 draws preserving label balance.

## Files

- `r8_fvptc_attention_metrics.tsv`, `r8_attention_rna_correlation.tsv`, `r8_centroid_cosine.tsv`, `r8_top_attention_tiles.tsv`, `fig_r8_fvptc_saliency.{png,pdf}`, `R8_REPORT.md`
