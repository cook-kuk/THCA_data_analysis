# Nature Medicine / Nature Cancer High-Impact Synopsis

## Proposed Title
Vaccine-ready AI for personalized cancer immunotherapy

## Standfirst
Personalized neoantigen vaccines have re-entered clinical oncology, but the AI layer that selects candidate antigens is still evaluated like a prediction leaderboard rather than a therapeutic decision system.

## 150-Word Abstract
Personalized neoantigen vaccines are becoming clinically credible across melanoma, pancreatic cancer and renal cell carcinoma, renewing the need to select patient-specific vaccine targets. Yet the computational layer that prioritizes neoantigens is often evaluated with retrospective AUROC leaderboards, random splits and public benchmark corpora that may not reflect the clinical decision: choosing a small number of manufacturable candidates for a new patient. We argue that neoantigen prediction should be reframed as vaccine-ready AI: contamination-controlled, source-aware top-k ranking under incomplete-label uncertainty. This Perspective proposes NEO-PRIOR, a practical reporting and evaluation checklist for neoantigen prioritization, covering public-corpus overlap, exact and near holdouts, HLA and study shift, AUPRC, top-k precision, enrichment, calibration and OOD-aware abstention. Establishing such standards is essential before AI-driven vaccine candidate selection can be interpreted as clinically meaningful rather than benchmark-specific.

## Killer Sentence
The next advance in personalized cancer vaccines may come less from predicting more HLA binders than from trusting the shortlist that decides what gets manufactured.

## Core Claims
1. Personalized vaccine trials have made neoantigen selection clinically consequential.
2. The clinically relevant output is a short ranked list, not a genome-wide AUROC.
3. Benchmark contamination and source shift can make weak deployment models look strong.
4. Public predictors are valuable, but their scores must be separated from clean model claims.
5. NEO-PRIOR can become the minimum reporting standard for vaccine-ready neoantigen AI.

## Display Items
1. Grand challenge overview: clinical reality -> evaluation gap -> field standard.
2. NEO-PRIOR reporting checklist.
3. Benchmark failure modes and split controls.
4. Vaccine-ready metric panel: AUPRC, top-k, enrichment, calibration, OOD.

## Why It Is High Impact
This is not a methods review. It is a standards paper for a clinically emerging AI decision layer.
