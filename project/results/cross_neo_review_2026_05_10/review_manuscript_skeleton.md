# From Neoantigen Discovery to Vaccine-Ready Prioritization

## Abstract Draft
Personalized neoantigen vaccines have re-emerged as a clinically plausible immunotherapy strategy, supported by durable T cell responses in pancreatic cancer, renal cell carcinoma, and melanoma studies. Yet the translational bottleneck is no longer simply whether candidate peptides bind HLA. In clinical and product settings, only a small number of candidates can be manufactured, assayed, or administered, making prioritization a top-k ranking problem under severe class imbalance, incomplete negative labels, source shift, and public benchmark overlap. This review synthesizes recent progress in neoantigen immunogenicity prediction, including binding-centric predictors, transfer-learned immunogenicity models, feature-based patient-context models, and structure/TCR-aware approaches. We argue that the next generation of neoantigen prioritizers should be evaluated by contamination-controlled splits, AUPRC, top-k precision, enrichment over prevalence, calibration, and OOD/abstention behavior rather than AUROC alone. We conclude with a practical evaluation contract for vaccine-ready computational prioritization.

## 1. Clinical Re-Emergence of Personalized Neoantigen Vaccines
- Melanoma: individualized mRNA therapy plus pembrolizumab.
- Pancreatic cancer: autogene cevumeran and durable CD8 T cell responses.
- Renal cell carcinoma: feasibility and immunogenicity in a lower mutation-burden tumor.
- Practical implication: better candidate ranking matters because manufacturing slots are scarce.

## 2. Why Binding Prediction Was Necessary but Is Not Sufficient
- MHC binding and presentation are required filters.
- Immunogenicity also depends on mutant-WT contrast, TCR-facing features, clonality, expression, tumor context, immune state, and tolerance.
- Binding-centric success can inflate apparent utility when evaluated on overlap-sensitive public corpora.

## 3. Current Predictor Families
- Binding/presentation predictors: NetMHCpan, MHCflurry, MixMHCpred, BigMHC-EL.
- Immunogenicity predictors: PRIME, DeepImmuno, BigMHC-IM, TransPHLA, T-SCAPE, IMPROVE.
- Dataset-harmonization and feature models: Muller et al., IMPROVE.
- Structure/TCR-aware and emerging PLM approaches.

## 4. Benchmarking Failure Modes
- Exact peptide-HLA overlap.
- Near-peptide and source-protein overlap.
- HLA allele and supertype shortcuts.
- Source/study/assay shift.
- Positive-label source bias and ambiguous negatives.
- AUROC-only reporting under low prevalence.

## 5. Metrics That Match Vaccine Decisions
- AUPRC over AUROC for rare positives.
- top5/top10 precision and recall@k.
- enrichment over prevalence.
- Brier/ECE calibration.
- abstention coverage versus precision.
- per-HLA, per-study, source-heldout, and time-heldout reporting.

## 6. Toward Vaccine-Ready Prioritization
- Mutant-WT counterfactual modeling.
- Retrieval evidence with explicit contamination flags.
- Structure/geometry and TCR-facing features.
- PU-aware treatment of negatives.
- Pan-allele generalization with HLA as context, not a memorized boundary.
- OOD-aware abstention and wet-lab triage.

## 7. Practical Evaluation Contract
- Train-only retrieval indices.
- No public predictor scores as clean training features if used as comparators.
- Exact/near/source/HLA/study/time split families.
- Separate clean manuscript claims from product-assisted triage claims.

## 8. Bridge to Original Work
- CROSS-Neo should be introduced later as an implementation of this review's contract.
- The review creates the rationale: not another MHC-binding predictor, but a contamination-controlled prioritization stack.

## Figure Plan
See `review_figure_plan.tsv`.

## Table Plan
See `review_table_plan.tsv`.
