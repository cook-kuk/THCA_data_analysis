# CROSS-Neo-TCR Feature Report

Generated deterministic TCR sequence, fallback k-mer embedding, and structure-missingness feature tables from the local TCR registry.

## Counts

| Field | Count |
|---|---:|
| registry_rows | 537618 |
| paired_alpha_beta | 64037 |
| beta_only | 49573 |
| alpha_only | 9425 |
| any_cdr3 | 123035 |
| peptide_hla_available | 406977 |
| known_structure_or_pdb | 899 |
| modelable_paired_tcr_pmhc | 56329 |

## Local PLM / Tool Status

| Component | Status |
|---|---|
| TCR_BERT_local_module | False |
| TCRpeg_local_module | False |
| transformers_local_module | True |
| torch_local_module | True |
| esm_local_module | True |
| prot_t5_status | not_run; fallback hashed k-mer embeddings emitted |
| embedding_dim_per_family | 64 |
| embedding_rows | 123035 |
| sequence_rows | 537618 |
| structure_rows | 537618 |

## Claim Boundary

- `tcr_plm_embeddings.parquet` currently contains deterministic hashed k-mer fallback embeddings, not trained TCR-BERT/TCRpeg/ESM2/ProtT5 outputs.
- Structure features are explicit missingness/QC placeholders unless a PDB/template exists; they should not be interpreted as interface evidence without parsed structures.
- These features are suitable for diagnostic pilots and wetlab prioritization triage, not as a standalone claim of full TCR-aware neoantigen SOTA.
