# NeoICI-combo algorithm gauntlet

Generated: 2026-05-11T00:47:16

## What was tested

Compared BigMHC, CROSS-Neo variants, KG-GA, BioDarwin v2/v3/v4, and new NeoICI combo scores on the labeled academic validation-like split and the locked public industrial mini-set.

## Academic validation-like top 5

- KG_GA_evolved: AUPRC 0.978, AUROC 0.999, top10 9/10.
- BioDarwin_PanVax_v2: AUPRC 0.886, AUROC 0.972, top10 9/10.
- BioDarwin_mode_bank_v4: AUPRC 0.886, AUROC 0.972, top10 9/10.
- NeoICI_GA_RL_combo_v1: AUPRC 0.863, AUROC 0.983, top10 8/10.
- CROSS_claimsafe: AUPRC 0.846, AUROC 0.992, top10 8/10.

## Industrial locked top 5

- BioDarwin_mode_bank_v4: AUPRC 0.628, AUROC 0.639, top10 5/10.
- BioDarwin_public_anchor_v3: AUPRC 0.628, AUROC 0.639, top10 5/10.
- BigMHC_IM: AUPRC 0.578, AUROC 0.611, top10 4/10.
- NeoICI_GA_RL_combo_v1: AUPRC 0.563, AUROC 0.597, top10 4/10.
- NeoICI_public_anchor_v1: AUPRC 0.535, AUROC 0.618, top10 5/10.

## Decision

Use NeoICI_GA_RL_combo_v1 as the production/experiment-priority score and NeoICI_claimsafe_combo_v1 as the reviewer-safe fallback. GA/RL remains valuable, but not as clinical proof.

## Files

- `neoici_scored_academic_candidates.tsv`
- `neoici_scored_industrial_candidates.tsv`
- `neoici_algorithm_gauntlet_metrics.tsv`
- `neoici_algorithm_winloss.tsv`
- `neoici_all_methods_context.tsv`
- `figures/Fig_NeoICI_algorithm_gauntlet_AUPRC.png`
- `project/papers_hub_2026_05_04/neo_ici_combo_algorithm_gauntlet.html`
