# Pre-Submission Synopsis: Nature Cancer Perspective

## Proposed Article Type
Perspective

## Proposed Title
Neoantigen prediction needs vaccine-ready benchmarks, not another AUROC leaderboard

## Alternative Titles
1. The benchmark problem in personalized neoantigen vaccines
2. Why top-k utility should replace AUROC leaderboards in neoantigen prioritization
3. Toward contamination-controlled evaluation of neoantigen vaccine rankers

## Core Argument
The field has built increasingly powerful neoantigen predictors, but evaluation has not kept pace with deployment. In vaccine design, the clinically relevant question is not whether a model separates positives from negatives across an entire retrospective corpus, but whether it can prioritize the handful of candidates worth manufacturing or testing for a new patient without relying on public-corpus overlap, HLA shortcuts or source-specific label bias.

## Perspective Thesis
Neoantigen prediction should be reframed as vaccine-ready prioritization: a top-k, source-aware, contamination-controlled ranking problem under label uncertainty.

## Why This Belongs in Nature Cancer
Nature Cancer Perspectives are intended to discuss models and ideas from a personal viewpoint, be forward-looking, and stimulate new research efforts. This Perspective would argue for a concrete shift in how the field evaluates neoantigen predictors:
- from AUROC to AUPRC/top-k/enrichment;
- from random splits to exact/near/source/HLA/study/time holdouts;
- from hidden retrieval gains to explicit contamination flags;
- from hard negatives to PU-aware and assay-aware interpretation;
- from universal scores to OOD-aware abstention.

## Proposed Structure
1. The clinical success of personalized vaccines makes prioritization urgent.
2. Binding prediction solved one layer but not vaccine selection.
3. Public benchmarks can reward memorization and source shortcuts.
4. Negatives are often unlabeled, assay-limited or patient-context-dependent.
5. Top-k precision and abstention should become primary endpoints.
6. A vaccine-ready evaluation contract.
7. Implications for trials, regulators and computational method papers.

## Proposed Display Items
1. A conceptual figure contrasting AUROC leaderboard thinking with vaccine-ready top-k selection.
2. A failure-mode map: exact overlap, near overlap, source shift, HLA shortcut, assay bias and incomplete negatives.
3. A practical evaluation checklist for future neoantigen predictor papers.

## Strong Opinion, Balanced Wording
The manuscript should not argue that public predictors are flawed or unusable. It should argue that they are powerful but often evaluated under conditions that do not match vaccine deployment. The strongest line is:

> The next advance in neoantigen prediction may come less from another model architecture than from asking a more clinically faithful benchmarking question.

## CROSS-Neo Boundary
CROSS-Neo should not be the centerpiece of this Perspective. It can be mentioned, if at all, as an example of the type of contamination-controlled, OOD-aware, top-k evaluation stack motivated by the Perspective.

## Expected Length
3,000-4,000 words, 2-4 display items, up to 100 references.
