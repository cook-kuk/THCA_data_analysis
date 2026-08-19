# NeoBench-Vax Public Launch Spec

## Minimum Public Website Sections
1. Benchmark overview.
2. Dataset registry and provenance cards.
3. Locked split downloads.
4. Public overlap audit reports.
5. Clean comparator leaderboard.
6. Product-assisted leaderboard.
7. Source-stress and OOD/abstention leaderboard.
8. Model-card template.
9. Submission instructions.
10. Claim-boundary policy.

## Leaderboard Columns
- model name;
- allowed inputs;
- public predictor score use;
- split;
- n;
- prevalence;
- AUPRC;
- AUPRC / prevalence;
- top5 precision;
- top10 precision;
- enrichment@10;
- recall@10;
- Brier;
- ECE;
- abstention coverage;
- source-heldout collapse penalty;
- claim eligibility.

## Governance Rule
A model may be useful in the product-assisted track while being ineligible for clean scientific superiority claims. This rule is the core trust mechanism.

## First Public Release
Release as a static website and GitHub-style benchmark repository before or alongside Paper 2 submission. The first version can include internal locked examples, missing-data manifests and public-overlap status, provided the claim boundary is explicit.
