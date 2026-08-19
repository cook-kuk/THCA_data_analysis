# Bridge From Review Paper to CROSS-Neo Original Paper

## Review Paper Role
Define the field problem:
- public predictors are powerful but overlap-sensitive;
- vaccine deployment is a top-k ranking problem;
- negative labels are often incomplete;
- source/HLA/study shifts can dominate performance;
- evaluation must be contamination-controlled and OOD-aware.

## Original CROSS-Neo Paper Role
Demonstrate one implementation:
- pan-allele mutant-WT counterfactual prioritizer;
- hard-decoy and PU-aware pressure for ambiguous negatives;
- explicit public-comparator separation;
- retrieval evidence marked as safe/unsafe;
- source-aware and OOD-aware triage;
- business/product mode separated from clean scientific claims.

## Paper Pairing
1. Review title: From neoantigen discovery to vaccine-ready prioritization.
2. Original title: CROSS-Neo: contamination-controlled pan-allele neoantigen prioritization for vaccine candidate triage.

## Clean Claim Boundary for Original Paper
- Internal locked/split-safe validation only unless an independent prospective set is obtained.
- No quantum advantage claim.
- Public predictor scores are comparators or product-assist signals, not clean model features.
- Main outcome should be AUPRC/top-k/enrichment and failure-mode audit.
