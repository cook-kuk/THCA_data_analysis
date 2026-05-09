# CROSS-Neo 2.0 Canonical Registry Schema

The canonical registry normalizes local neoantigen immunogenicity rows without silently dropping missing fields.

Required identifiers:
- `row_id`
- `source_dataset`
- `study_id`

Core biology fields:
- `peptide`
- `mutant_peptide`
- `wildtype_peptide`
- `hla_raw`
- `hla_4digit`
- `hla_gene`
- `hla_supertype`
- `peptide_length`
- `mutation_position`
- `source_protein`
- `source_window`

Labels:
- `label_raw`
- `label_binary`
- `label_strength`

Evidence/context:
- `expression`, `TPM`, `RNA_evidence`
- `clonality`, `VAF`
- `train_test_original`
- `public_overlap_flags`

Deduplication/leakage keys:
- `exact_peptide_hla_key`
- `peptide_only_key`
- `mutant_wt_hla_key`
- `near_peptide_cluster_key`

Missingness is explicit via `missing_*` columns and must be treated as signal.
