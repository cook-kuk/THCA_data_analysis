# CROSS-Neo-TCR Linkage Report

## Claim Boundary

TCR evidence is added as a diagnostic annotation layer only. The linker does not transfer labels from public TCR resources into the main CROSS-Neo immunogenicity task.

## Linkage Counts

- Main CROSS-Neo rows: 2,715
- `exact_tcr_pmhc_match`: 101
- `peptide_hla_match`: 1,731
- `peptide_only_match`: 14
- `near_peptide_match`: 170
- `recurrent_cancer_gene_context_match`: 0
- `no_tcr_match`: 699

## Summary Table

| tcr_link_category    |   n_neo_rows |   mean_tcr_evidence_count |   total_tcr_evidence_count |   paired_tcr_evidence_count |   cancer_context_evidence_count |   pathogen_context_evidence_count |
|:---------------------|-------------:|--------------------------:|---------------------------:|----------------------------:|--------------------------------:|----------------------------------:|
| peptide_hla_match    |         1731 |                   1.79318 |                       3104 |                           0 |                            2626 |                                32 |
| no_tcr_match         |          699 |                   0       |                          0 |                           0 |                               0 |                                 0 |
| near_peptide_match   |          170 |                  19.9529  |                       3392 |                         790 |                             498 |                              1249 |
| exact_tcr_pmhc_match |          101 |                  31.8317  |                       3215 |                         810 |                            2535 |                                93 |
| peptide_only_match   |           14 |                   1.92857 |                         27 |                           0 |                              20 |                                 0 |

## Interpretation Rules

- Exact peptide+HLA with paired TCR is the strongest local evidence, but still diagnostic until source, assay, and split compatibility are checked.
- Peptide-only matches are not definitive because MHC restriction and antigen source can differ.
- Pathogen-derived matches are not transferred to cancer neoantigens without clear cancer annotation.
- Recurrent-gene matches are context only; they do not imply TCR recognition of the same neoepitope.

## Outputs

- `tcr_neo_linked_registry.parquet` / `tcr_neo_linked_registry.tsv`
- `tcr_neo_linkage_summary.tsv`
- `tcr_neo_linkage_evidence_examples.tsv`
