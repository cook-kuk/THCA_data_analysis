# Nature Biomedical Engineering Presubmission Synopsis

## Proposed Format
Article

## Proposed Title
Benchmarking vaccine-ready AI systems for neoantigen prioritization

## 175-Word Abstract Draft
Personalized neoantigen vaccines require computational systems that choose a small number of candidate peptides for manufacture or immune testing. Yet most neoantigen predictors are benchmarked as retrospective classifiers, often under public-corpus overlap, source shift, HLA shortcuts and incomplete negative labels. Here, we present CROSS-Neo, a contamination-controlled benchmarking and prioritization system for vaccine-ready neoantigen ranking. CROSS-Neo integrates mutant-WT counterfactual sequence features, hard-decoy pressure, fixed sequence/QK fallback signals, explicit retrieval-risk flags, source-aware calibration and OOD-aware abstention. We evaluate the system using locked split families that remove exact peptide-HLA overlap, near-peptide overlap, HLA/supertype shortcuts and source/study leakage. Primary outcomes are AUPRC, top-k precision, enrichment over prevalence and abstention behavior rather than AUROC alone. In internal retrospective benchmarks, the pan-allele product ranker improves top-k prioritization over public predictor-only baselines, while source-heldout stress tests reveal where abstention remains necessary. CROSS-Neo is therefore presented as a benchmarked validation framework, not an external-valid clinical product.

## Editorial Pitch
The article fits Nature Biomedical Engineering because it concerns the validation and deployment requirements of a computational system that may facilitate personalized therapy design. The central contribution is not a black-box model, but a leakage-controlled evaluation architecture and reporting standard for AI-guided vaccine candidate selection.

## Display Items
1. CROSS-Neo validation-system overview.
2. Leakage-safe split contract.
3. Internal locked competitor gauntlet.
4. Source-heldout stress and abstention boundary.
5. Claim-boundary decision tree.
6. NEO-PRIOR reporting checklist.
