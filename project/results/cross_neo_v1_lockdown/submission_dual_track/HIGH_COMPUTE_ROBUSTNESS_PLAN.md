# High-Compute Robustness Plan

Status: execution plan.

## Principle

Training time is not important here. Reviewer confidence is more important than runtime. Any expensive run is acceptable if it is split-safe, reproducible, and answers a reviewer-relevant question.

## Do Not Optimize

- wall-clock training time
- number of models tried
- AUROC-only gains
- leaderboard-style internal repeated CV

## Optimize

- AUPRC
- top5/top10/top20 precision
- enrichment over prevalence
- source/HLA/near-neighbor robustness
- calibration and abstention
- leakage audit completeness
- case-level biological interpretability

## High-Compute Experiments To Add

1. Full nested weight selection
   - Replace quick inner OOF selection with repeated inner CV per outer fold.
   - Record selected weights, variance, and failure cases.
   - No outer test labels.

2. Large frozen PLM embeddings
   - ESM2/ProtT5 embeddings for mutant, WT, delta, source window, and HLA pseudo/full sequence.
   - No fine-tuning in first pass.
   - Compare to current AA/k-mer fallback.

3. Structure-token branch
   - Use all completed ESMFold pseudo-complexes.
   - If available, add Boltz/AF3 structures.
   - Extract geometry plus structure-token embeddings.
   - Keep structure_missing/confidence as explicit flags.

4. Source-robust training
   - Leave-source-out training with source-balanced batches.
   - GroupDRO-style source loss weighting.
   - Label-shift calibration using train-only source priors.
   - Evaluate only on heldout source labels.

5. Public overlap completion
   - Download/parse official training corpora.
   - Exact peptide, exact peptide-HLA, near peptide, HLA-normalized, source-protein overlap.
   - Comparator scores remain forbidden as features.

6. Case-level biological audit
   - Top 20 positives rescued by fusion/QK.
   - Top 20 negatives harmed by fusion/QK.
   - Annotate HLA, peptide length, source, retrieval flag, structure confidence, WT availability.

## GPU/RunPod Priority

1. Frozen PLM embeddings.
2. Structure-token extraction.
3. Repeated nested source-robust training.
4. Bootstrap CIs and permutation tests.

## Promotion Criteria

Promote from HOLD to KEEP only if:

- fold-safe fusion improves AUPRC/top10 on at least 3 of 4 primary locked splits;
- source-heldout NEPdb top10 no longer collapses;
- TESLA top20 or abstained top-k shows nonzero positive recovery;
- public overlap audit is complete enough to state comparator cleanliness boundaries;
- QK rescue/harm supports bounded fallback role.

