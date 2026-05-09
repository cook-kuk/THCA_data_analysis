# CROSS-Neo 2.0 Initial Diagnosis

## What Is Already Strong

- v1 has reviewer-safe internal locked evidence: fold-safe fusion improves AUPRC or top-k over `anchor_rf` on the four primary locked splits.
- The clean anchor is not retrieval-dependent: `C_counterfactual / rf_secondary` and `C_counterfactual / elastic_net_lr` remain the defensible base.
- QK has a bounded rescue role in some locked splits, especially exact peptide-HLA and near peptide cluster holdouts, but it is not a headline claim.
- The v1 package already separates public predictor scores from model features.

## What Blocks SOTA Claims

- Source-heldout performance still collapses, especially NEPdb/TESLA top-k recovery.
- Public overlap with MHCflurry, NetMHCpan, BigMHC, PRIME, MixMHCpred, NetMHCstabpan is unresolved.
- Fusion weights are useful but not fully stable, so any weight larger than prespecified/fold-safe choices must remain exploratory.
- Current evidence is internal/locked only. It is not external validation.

## Promotion Targets From HOLD To KEEP

- Beat both `anchor_rf` and the v1 best fold-safe fusion on at least 3 of 4 primary locked splits by AUPRC or top10 precision.
- NEPdb source-heldout must no longer have top10 precision equal to 0.
- TESLA_mmc4 or TESLA_mmc7_validation must show nonzero top20 recovery, or a justified abstained-top-k recovery.
- Public overlap audit must be completed, or all public-comparator/SOTA claims must be downgraded.
- QK must remain fallback/diagnostic unless rescue/harm is clearly favorable and harms are bounded.

## Forbidden Claims

- No external validation claim.
- No clinical utility or vaccine-selection claim.
- No quantum advantage claim.
- No public predictor cleanliness claim until public overlap is resolved.
- No SOTA claim based only on internal locked splits.
