# Manuscript Results Scaffold KR

## Result 1. Candidate universe and track separation
We standardized candidate records from existing CLEAN-Neobench/cross-neo/neoantigen hub artifacts into a canonical schema. The key design is a hard split between a clean-science track and a production-stack track.

Claim boundary: this is data/model integration, not clinical validation.

## Result 2. Public predictor scores are useful but cannot define clean novelty
MHCflurry, BigMHC, PRIME, and NetMHCpan-style artifacts can be collected as frozen comparators/features. They are excluded from clean-science training.

Claim boundary: binding/presentation predictors are not immunogenicity proof.

## Result 3. Local algorithms provide the clean scientific core
Structure_LR, Wave8/TCR-self-similarity, ESM2_Bayesian, W7A/W7B, and quantum-kernel families are registered as local branches. In strict no-existing-overlap views, TCR/self-similarity branches become a key biological story.

Claim boundary: strict public subset signal is promising but must be confirmed under true source-heldout/patient-heldout splits.

## Result 4. Leakage audit changes the interpretation of performance
The pipeline reports exact peptide, peptide-HLA, mutant-WT, study, patient, HLA, source-window, and public predictor training-contamination risk where available.

Claim boundary: random split and existing-artifact performance are smoke-test evidence only.

## Result 5. Production stack is strong as an operating layer
The production stack combines local branches with frozen external predictors and patient-context fields when available.

Claim boundary: production stacking can support prioritization, not a clean algorithm novelty claim.

## Result 6. Patient-level top-N is the right endpoint but needs hospital data
The code writes patient top-20 outputs and patient-level metrics, but the current public integrated table lacks real patient IDs.

Claim boundary: patient-level Recall@20/hit rate is not claimable until real patient-level candidate sets are obtained.

## Result 7. Wet-lab handoff is now explicit
The package defines what to test, what to hold, and what metadata blocks a candidate. This converts algorithm output into a collaborator-facing assay queue.

Claim boundary: candidates are hypotheses, not validated vaccine products.
