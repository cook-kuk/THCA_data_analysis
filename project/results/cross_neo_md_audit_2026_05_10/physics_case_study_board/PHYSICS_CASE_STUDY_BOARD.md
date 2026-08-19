# CROSS-Neo physics case study board

## Executive verdict

The strongest story now is no longer just a ranked list. It is a two-candidate decision board with a physics escalation path and explicit hold rules for PMF/FEP.

## Summary board

| lead                 | peptide    | hla_4digit   | md_label       |   md_score |   preclinical_readiness |   cross_neo_i_readiness | physics_tier        |   physics_priority_score | launch_sequence                                                                                          | endpoint_step                                                        | next_decision                                                         | blocked_claim                                             | claim_boundary                                                   |
|:---------------------|:-----------|:-------------|:---------------|-----------:|------------------------:|------------------------:|:--------------------|-------------------------:|:---------------------------------------------------------------------------------------------------------|:---------------------------------------------------------------------|:----------------------------------------------------------------------|:----------------------------------------------------------|:-----------------------------------------------------------------|
| TP53_R175H_HLA_A0201 | HMTEVVRHC  | HLA-A*02:01  | MD_VERY_STRONG |   0.821789 |                0.935447 |                0.605838 | P0_FLAGSHIP_PHYSICS |                 0.839191 | endpoint energy now -> PMF if endpoint supports mutant specificity -> FEP only after WT mapping is clean | compare mutant vs WT vs decoy endpoint energies over ensemble frames | launch WT/decoy-controlled HLA stability + multimer/activation assays | immunogenicity until activation exceeds WT/decoy controls | physics and assay support ranking only, not immunogenicity proof |
| KRAS_G12D_HLA_C0802  | GADGVGKSAL | HLA-C*08:02  | MD_MODERATE    |   0.577335 |                0.796652 |                0.61065  | P0_FLAGSHIP_PHYSICS |                 0.656129 | endpoint energy now -> PMF only if TCR-pMHC model stays stable                                           | compare mutant vs WT vs decoy endpoint energies over ensemble frames | launch WT/decoy-controlled HLA stability + multimer/activation assays | immunogenicity until activation exceeds WT/decoy controls | physics and assay support ranking only, not immunogenicity proof |

## Claim boundary

This board can prioritize experiments and reject weak candidates. It cannot prove immunogenicity without WT/decoy-controlled activation and killing labels.

