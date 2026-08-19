# BAR-Neo-NG Algorithm Card

## What It Is

BAR-Neo-NG is a claim-gated BAR-Neo layer. It combines stress-guarded benchmark reliability, reviewer kill-audit status, public/fallback dependency penalties, HLA support, assay readiness, and translational metadata gates.

## What It Does Not Claim

It is not a new SOTA neoantigen predictor, not clinical vaccine selection, not external validation, and not a QK headline claim.

## Decision Logic

- High leakage or exact train overlap caps/blocks clean claims.
- Near-neighbor caveats become `T2_resolve_before_claim`.
- Kill-audit pass plus pHLA assay readiness becomes `T1_phla_assay_design_candidate`.
- Missing translational metadata caps the T1 score and blocks patient/translational claims.
- Public pretrained support is penalized unless row-level overlap audit is clean.
