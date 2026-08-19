# Nature Machine Intelligence Analysis Synopsis

## Proposed Format
Analysis

## Proposed Title
Contamination-controlled benchmarking of neoantigen prioritization models

## Core Claim
Neoantigen prediction leaderboards can change substantially when models are evaluated under exact/near overlap controls, source/HLA holdouts and top-k metrics that reflect vaccine candidate selection.

## Why Analysis
Nature Machine Intelligence defines Analysis as a new analysis of existing data or comparative analysis leading to novel conclusions for a broad audience. This manuscript would be framed as a comparative AI evaluation study, not a wet-lab validation paper.

## What Must Be Strengthened Before Submission
- Expand public competitor coverage where possible: IMPROVE, MixMHCpred standalone, DeepHLApan, GraphMHC if runnable locally.
- Add bootstrap/permutation tests for the final candidate.
- Add model-card and data-card package.
- Make the benchmark registry reproducible from raw or intermediate files.
- Keep biological claims modest and emphasize AI evaluation.
