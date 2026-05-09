
# CLEAN-NeoBench Method Card

## What This Framework Does

CLEAN-NeoBench is a leakage-aware AI neoantigen predictor benchmarking framework. It normalizes candidate identity, HLA metadata, labels, public/internal predictor outputs, overlap flags, split contracts, and reviewer-facing metrics into one rerunnable benchmark surface.

## What It Does Not Claim

- No clinical vaccine selection.
- No new SOTA predictor claim.
- No external validation proven.
- No standalone QK/quantum superiority claim.
- No statement that public pretrained tools are clean baselines without row-level training-corpus overlap audit.
- No unified Class I/Class II predictor.

## Benchmark Contracts

- Exact peptide-HLA holdout.
- Near peptide cluster holdout using existing clusters or k-mer/identity fallback.
- Source protein/window holdout when metadata exist.
- Source/study heldout, including CEDAR, NEPdb, and TESLA-like stress groups when present.
- HLA allele and HLA supertype heldout.
- Patient/study heldout when metadata exist.
- Low-prevalence heldout for TESLA-like settings.
- Korean-HLA focus board for A*24:02, A*11:01, A*02:01, B*15:01, B*40:01, C*01:02, C*03:03, C*07:02.

## Method Zoo

- `internal_candidate`: 57054
- `bounded_fallback`: 20198
- `anchor`: 6066
- `caveated_public_comparator`: 3224
- `uncertainty_only`: 400

## Leakage Audit Design

Overlap flags include exact peptide, exact peptide-HLA, near peptide cluster, source protein/window, study, patient, and public-tool training overlap status. Public pretrained methods are caveated unless official training corpora are row-audited.

## Metrics

CLEAN-NeoBench emphasizes AUPRC, top-k precision/recall/hit, positive rank, calibration Brier/ECE, confidence coverage, risk at high confidence, and abstention rate. AUROC is reported but is not the sole decision metric.

## Limitations

Missing WT peptide, source window, patient, expression, clonality, HLA LOH, and B2M metadata reduce the strength of source-window and patient-gated claims. Missing metrics are reported as unavailable rather than fabricated.

## Claim Boundary

Allowed claims: leakage-aware benchmark, benchmark-adaptive reliability ranking, calibrated candidate prioritization with abstention, reviewer-safe comparison of internal and public AI predictors, and research triage framework.
