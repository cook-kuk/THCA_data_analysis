# Non-Confidential Teaser: CROSS-Neo KG-GA Selection Layer

## One-Line Positioning

CROSS-Neo KG-GA is an algorithmic neoantigen-selection layer designed to improve which candidates enter a fixed-slot individualized mRNA vaccine payload, before manufacturing.

## Why This Fits A V940-Like Workflow

- Public Moderna/Merck materials describe V940/mRNA-4157 as an individualized neoantigen therapy encoding up to 34 algorithmically derived neoantigens.
- The value gap is therefore not mRNA manufacturing; it is better candidate selection under a constrained payload budget.
- CROSS-Neo outputs ranked candidates, uncertainty, reason codes, and wetlab-control-ready exports.

## Masked Evidence Snapshot

- Same-board 2,715-candidate result: KG-GA AUPRC 0.933 / AUROC 0.943.
- Fixed integrated comparator: AUPRC 0.878.
- Public BigMHC context on the same board: IM AUPRC 0.672; EL AUPRC 0.648.
- Frozen validation-like subset: KG-GA AUPRC 0.978.
- Strict no-overlap reliability remains a separate lane: BAR-Neo_confidence AUPRC 0.857; KG-GA-controller AUPRC 0.624.

## Validation Ask

Under CDA/NDA, provide a blinded candidate pool with sponsor labels withheld. We return frozen scores, top-34 ranking, uncertainty, and assay-control recommendations. Labels are unblinded only after score lock.

## Pre-NDA Boundary

No code, full feature list, candidate tables, training corpus details, or proprietary scoring weights should be shared before IP filing and CDA/NDA.

## Official Context Sources

- [Merck/Moderna melanoma Phase 3 announcement](https://www.merck.com/news/merck-and-moderna-initiate-phase-3-study-evaluating-v940-mrna-4157-in-combination-with-keytruda-pembrolizumab-for-adjuvant-treatment-of-patients-with-resected-high-riskstage-iib-iv-melanom/) — V940/mRNA-4157 Phase 3 adjuvant melanoma program; approximately 1,089 patients; RFS primary endpoint.
- [Merck/Moderna NSCLC Phase 3 announcement](https://www.merck.com/news/merck-and-moderna-initiate-phase-3-trial-evaluating-adjuvant-v940-mrna-4157-in-combination-with-keytruda-pembrolizumab-after-neoadjuvant-keytruda-and-chemotherapy-in-patients-with-certain-ty/) — INTerpath-009 expansion; V940 described as mRNA coding for up to 34 algorithmically derived neoantigens.

## Contact

- Seungho Cook
- Email: kukshomr@snu.ac.kr
- X/Twitter: @Ho15421542

