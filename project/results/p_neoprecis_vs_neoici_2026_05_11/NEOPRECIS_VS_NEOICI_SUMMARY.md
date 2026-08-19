# NeoPrecis vs NeoICI reanalysis

## What is directly comparable

- NeoPrecis is a published immunogenicity + landscape framework with CEDAR/NCI/melanoma/NSCLC benchmarks.
- NeoICI is our local combo-therapy prioritization layer tested on the current gauntlet tables.

## What is not yet directly comparable

- I do not have a same-dataset score table for NeoPrecis and NeoICI in this workspace.
- So this is a performance-contract comparison, not a strict paired re-run.

## Current local result

- `KG_GA_evolved` remains the strongest local scorer on the validation-like split: AUPRC 0.978.
- `NeoICI_GA_RL_combo_v1` is strong but not the top local scorer: AUPRC 0.863.
- `BioDarwin_mode_bank_v4` still anchors the locked industrial mini-set: AUPRC 0.628.

## Decision

- No, we cannot say the new package is universally better than NeoPrecis yet.
- Yes, we can say the package is broader in scope: it adds combo-therapy readiness, GA/RL search, and explicit claim boundaries.
