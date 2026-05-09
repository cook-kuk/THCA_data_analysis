# CROSS-Neo-TCR Extension Decision Report

Decision: **PROMOTE_TO_DIAGNOSTIC_CASE_STUDY** and **PROMOTE_TO_WETLAB_PRIORITIZATION_TOOL** for the current sprint. Hold main-method promotion until paired TCR benchmarks, stronger TCR sequence models, and structure validation are complete.

## Required Answers

1. Paired TCR alpha/beta rows: **64,037** in the local TCR registry.
2. Peptide-HLA-TCR labeled rows: **98,718** with any TCR sequence; **56,329** with paired alpha/beta.
3. Cancer neoantigen overlap: **2,016/2,715** CROSS-Neo rows have some TCR-resource overlap; **101** exact paired TCR-pMHC rows; **1,876** rows have cancer-context evidence.
4. Top-k/ranking improvement: internal/exact/near pilots improve with TCR evidence. Conservative exact holdout: AUPRC 0.588 -> 0.666 (delta +0.078), AUROC delta +0.094; near holdout: AUPRC 0.393 -> 0.558 (delta +0.165), AUROC delta +0.133.
5. Structure over sequence-only: **not yet supported**. Structure files currently contain missingness/QC features unless PDB/template evidence exists.
6. Mutant-WT interface delta: **not yet supported**. Needs parsed mutant/WT TCR-pMHC structures.
7. Source-heldout rescue: partial. Conservative NEPdb source-heldout readout: AUPRC 0.233 -> 0.259 (delta +0.026), AUROC delta +0.050. Raw TCR evidence is stronger but treated as leakage-prone.
8. False positives: case audit files identify main-high/TCR-low and TCR-high-only negatives, but explanation remains diagnostic until structures or external TCR assays support it.
9. Main claim or supplement: **supplement/diagnostic branch now**, not the main CROSS-Neo ranking claim.
10. Wetlab candidates: prioritize rows in `tcr_wetlab_candidate_prioritization.tsv` and the de-duplicated `tcr_wetlab_candidate_prioritization_unique_pmhc.tsv` with high conservative TCR-augmented score, positive delta, exact paired or cancer-context TCR evidence, and low pathogen-only dependence.

## Decoy Recognition Pilot

- Paired positive versus shuffled-TCR decoy pilot best result: pmhc_only AUPRC 0.500, AUROC 0.500. This indicates the simple handcrafted TCR sequence features are not sufficient for de novo cognate-recognition prediction.

## Claim Boundary

- Do not replace the main pMHC neoantigen model with the TCR branch.
- Do not transfer pathogen epitope labels to cancer neoantigens.
- Do not claim clinical utility or universal TCR-aware prediction.
- Safe current claim: optional TCR evidence/diagnostic layer plus wetlab prioritization scaffold.
