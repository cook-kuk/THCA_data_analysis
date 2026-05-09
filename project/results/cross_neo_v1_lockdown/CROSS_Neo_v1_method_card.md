# CROSS-Neo v1 Method Card

## Intended Use

Internal locked-split neoantigen prioritization research. Not clinical vaccine selection and not external validation.

## Locked Components

- Clean anchors: `anchor_rf`, `anchor_lr` from mutant-WT counterfactual features.
- Fallback: fixed QK branches, used only through prespecified or fold-safe fusion.
- Forbidden: public predictor scores as model features.

## Best Locked Evidence

- Best primary fusion row: `prespecified_rf_qk_no_anchor_w0.5` on `exact_peptide_hla_holdout`; AUPRC 0.587, AUROC 0.757, top10 0.700.
- Fold-safe fusion improves the clean anchor on all four primary locked splits by AUPRC or top10.

## Known Failure

- Source-heldout top-k collapse remains in: NEPdb, TESLA_mmc4, TESLA_mmc7_validation.
- Public overlap for several pretrained comparators remains unresolved.

## QK Boundary

- Total audited rescued positives: 119.
- Total audited harmed negatives: 419.
- Interpretation: fallback/fusion candidate only; no quantum advantage claim.
