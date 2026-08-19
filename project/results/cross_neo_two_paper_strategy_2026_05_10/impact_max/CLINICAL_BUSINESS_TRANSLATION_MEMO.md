# Clinical and Business Translation Memo

## Product Wedge
CROSS-Neo should be packaged as an auditable vaccine-candidate triage system, not as a standalone biological truth engine.

## Buyer/User Pain
Personalized vaccine teams face too many candidates and too few assay/manufacturing slots. They need ranked shortlists, uncertainty flags, and a way to justify which candidates move forward.

## Practical Product Modules
1. Candidate queue: top-k ranking with enrichment and prevalence baseline.
2. Evidence tags: clean-no-reference, near-retrieval-supported, exact-retrieval-supported, OOD.
3. Abstention mode: flag unsafe source/HLA/length/structure/retrieval cases.
4. Comparator panel: optional public predictor-assisted view, clearly labeled.
5. Audit report: NEO-PRIOR checklist, data provenance, overlap and claim-boundary table.

## Demo Claim
Acceptable:

> In internal retrospective testing, the product-assisted queue concentrated positives in the top candidate list better than any single available score in this dataset.

Not acceptable:

> Externally validated, universally superior, or quantum-advantaged cancer vaccine predictor.

## Business-Safe Differentiation
The differentiator is not a secret score. It is the audit layer: leakage controls, source-risk labeling, top-k decision metrics, and abstention.
