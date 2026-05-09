# CROSS-Neo Next 72h Action Board

Status: execution checklist.

## Priority 1: Public Overlap Closure

- Download or locate official training sets:
  - MHCflurry
  - NetMHCpan
  - NetMHCstabpan
  - BigMHC
  - PRIME
  - MixMHCpred
  - IEDB export
- Parse into a common schema:
  - peptide
  - HLA
  - label or measurement
  - source/study
  - publication/date if available
- Re-run exact/near overlap.

## Priority 2: High-Compute Embeddings

- Generate ESM2 or ProtT5 frozen embeddings for:
  - mutant peptide
  - WT peptide
  - mutant-WT delta
  - source window
  - HLA pseudo/full sequence
- Save fold-independent features only.
- Re-run v1 lockdown with same splits.

## Priority 3: Source-Robust Training

- Add source-balanced and GroupDRO-style training.
- Leave each source out.
- Tune only on train sources.
- Evaluate top5/top10/top20 on heldout source.

## Priority 4: Reviewer Figures

- Rebuild seven v1 figures with final high-compute branches.
- Add prevalence baselines everywhere.
- Mark small-n splits as descriptive.

## Priority 5: Dual Submission Packs

- ED/evaluation paper:
  - emphasize evaluation design and failure modes.
- Bio/workshop paper:
  - emphasize CROSS-Neo as stress-tested candidate prioritizer.

## Stop Conditions

- If public overlap shows hidden contamination explaining gains, downgrade predictor claim and convert paper fully to evaluation/audit paper.
- If source-heldout remains collapsed after high-compute features, do not claim predictor readiness.
- If QK harm exceeds rescue in top-k, keep QK only as diagnostic/fallback, not headline.

