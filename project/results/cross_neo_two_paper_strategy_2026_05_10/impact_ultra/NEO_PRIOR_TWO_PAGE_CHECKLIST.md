# NEO-PRIOR Two-Page Checklist

## Page 1: Minimum Reporting Items
1. Define the prediction task: binding, presentation, immunogenicity, triage or benchmark ranking.
2. Report data provenance: study, patient, assay, HLA, candidate-generation route and date/year when available.
3. Disclose public predictor use as feature, comparator, filter or product assist.
4. Report exact peptide and exact peptide-HLA overlap.
5. Report near-peptide, source-window, mutation-pair and TCR motif overlap when available.
6. Keep retrieval, scaling, imputation, calibration, feature selection and thresholding inside train folds.
7. Use exact, near, source/study, HLA allele/supertype, time/assay and retrieval-clean splits where feasible.
8. Treat negatives as ambiguous unless experimentally established as true non-immunogenic cases.
9. Report AUPRC, top-k precision, enrichment over prevalence and recall@k as primary metrics.
10. Report AUROC as secondary, not as the primary evidence of shortlist utility.

## Page 2: Deployment and Claim Boundary
11. Report Brier/ECE, calibration curves and score reliability.
12. Report OOD flags and abstention coverage versus precision.
13. Report per-source and source-heldout results.
14. Report per-HLA and HLA-heldout/supertype-heldout results.
15. Publish split manifests, row IDs, seeds and scripts where possible.
16. Publish model cards with allowed inputs, forbidden inputs and known failures.
17. Separate clean scientific validation from product-assisted triage.
18. State whether public predictor scores are forbidden, allowed or disclosed.
19. State whether the evidence is internal, locked retrospective, external retrospective or prospective.
20. Include a one-paragraph claim boundary in the abstract or discussion.
