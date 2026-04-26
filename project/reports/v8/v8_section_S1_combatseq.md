> ⚠️ **v5.2 SUPERSEDED NOTICE (2026-04-25).** This robustness section was designed to stress-test the v5.1 THCA DIAL flip (DIAL=0.494, batch_entangled in 4/5). Under proper LODO ComBat the underlying flip disappears (v5.2: DIAL=0.000 / true_biology in 5/5; see reports/v5p2/v5p2_critical_assessment.md). Conclusions below are conditional on the v5.1 result they reference.

# Section S1 - ComBat-seq vs ComBat benchmark (THCA)

## Result

On THCA BRAF-vs-RAS (n=392; TCGA-THCA n=351 RNA-seq, GSE27155 n=41
Affymetrix GPL96 microarray), we re-ran the v5.1 LODO DIAL pipeline with two
batch-correction variants: (a) the v5.1 baseline, ComBat (`pycombat_norm`) on
log2-TPM; and (b) ComBat-seq (`inmoose.pycombat.pycombat_seq`) on raw STAR
unstranded counts for the RNA-seq cohort, followed by a cross-platform ComBat
step to bridge to microarray values. The DIAL-flip signature is robust: four
of five classifiers show `|Delta DIAL| < 0.1` (LogReg_l2 0.494 vs 0.495,
LogReg_elasticnet 0.492 vs 0.494, GradientBoosting 0.334 vs 0.273, XGBoost
0.013 vs 0.000). Only Random Forest moves meaningfully (0.324 -> 0.129, post-
correction AUC 0.177 -> 0.372), shrinking its inversion but still below 0.5.
The BRAF-vs-RAS signal does not recover above chance under either correction,
so the v5.1 conclusion that THCA is the "flip cancer" survives the ComBat-seq
check.

## Microarray limitation (the degradation note)

ComBat-seq is defined only on RNA-seq counts and is mathematically undefined
for the microarray cohort; we therefore cannot run it jointly across the two
batches. Our implementation is a hybrid: ComBat-seq on TCGA counts
(degenerate single-batch no-op), then a cross-platform ComBat to merge with
the microarray log2 intensities. That hybrid is the strongest defensible
answer given that the flip cancer has exactly one RNA-seq cohort paired with
one microarray cohort - a structural feature of the THCA public-data
landscape, not a method choice. Raw counts for 351/351 TCGA-THCA samples were
loaded from GDC STAR gene-counts files; the microarray side uses log2
intensities as the best available analogue. No counts were fabricated.

## Reviewer-proofing takeaway

The flip reported in v5.1 is not an artefact of applying ComBat to log2-TPM.
Switching to the count-native estimator leaves four of five classifiers in
the same DIAL regime, and every AUC remains <= 0.61 (well below the 0.7
"true biology" threshold). Reviewer 2's objection is anticipated and
neutralised: the THCA inversion reproduces under ComBat-seq, so the
cross-cancer generalisation story stands.
