# CROSS-Neo-TCR Case/Error Analysis

This is a diagnostic pilot over the current subset-comparison predictions. The conservative comparison uses source-excluded, no-public-label TCR evidence features.

## Registry And Linkage Counts

- TCR registry rows: 537,618
- Paired alpha/beta TCR rows: 64,037
- Any-TCR peptide-HLA labeled rows: 98,718
- Paired-TCR peptide-HLA labeled rows: 56,329
- CROSS-Neo rows with any TCR-resource overlap: 2,016 / 2,715
- CROSS-Neo rows with exact paired TCR-pMHC overlap: 101 / 2,715
- CROSS-Neo rows with cancer-context TCR evidence: 1,876 / 2,715

## Conservative Case Counts

- TCR rescue cases: 62
- TCR harm cases: 85
- False-positive audit rows: 252
- False-negative audit rows: 86

## Split-Level Mean Score Shift

| Split | n | positives | mean pMHC | mean TCR-aug | mean delta |
|---|---:|---:|---:|---:|---:|
| exact_peptide_hla_holdout | 89 | 21 | 0.379 | 0.374 | -0.005 |
| near_peptide_cluster_holdout | 89 | 21 | 0.390 | 0.373 | -0.017 |
| repeated_stratified_5x5_internal | 445 | 105 | 0.384 | 0.379 | -0.005 |
| source_heldout_NEPdb | 572 | 151 | 0.471 | 0.218 | -0.253 |

## Interpretation

- Raw TCR evidence remains useful as a diagnostic annotation, but can be inflated by source self-evidence.
- The conservative source-excluded/no-label branch is the safer readout for claim discussions.
- Handcrafted sequence/TCR-evidence features do not establish de novo TCR recognition. A true TCR expert still needs ERGO-II/NetTCR/pMTnet/PLM or structure-backed validation.
