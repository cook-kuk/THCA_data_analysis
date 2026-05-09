# CROSS-Neo-TCR Registry Report

## Scope

This registry is an optional TCR-aware evidence layer. It is not a replacement for the main CROSS-Neo pMHC ranking registry because most neoantigen rows do not contain paired TCR alpha/beta chains.

## Headline Counts

- Registry rows retained: 537,618
- Rows with paired TCR alpha/beta CDR3: 64,037
- Rows with beta-only TCR CDR3: 49,573
- Rows with peptide + normalized HLA: 406,977
- Rows with cancer-context annotation: 72,291
- Rows with pathogen-context annotation: 298,369
- Rows with PDB/structure identifier: 899

## Source Balance

| source_dataset         |      n |   paired_tcr |   beta_only_tcr |   alpha_only_tcr |   any_tcr |   peptide_hla |   positives |   prevalence |   cancer_context |   pathogen_context |   structures |
|:-----------------------|-------:|-------------:|----------------:|-----------------:|----------:|--------------:|------------:|-------------:|-----------------:|-------------------:|-------------:|
| IEDB_tcell_full_export | 397627 |            0 |               0 |                0 |         0 |        291583 |      137415 |     0.345588 |            52112 |             193413 |          308 |
| VDJdb                  |  62177 |        30594 |           24128 |             7455 |     62177 |         62177 |       62177 |     1        |              343 |              52765 |          278 |
| McPAS-TCR              |  40731 |        12717 |           25430 |             1966 |     40113 |         16304 |       40731 |     1        |             3880 |              31855 |            0 |
| VDJdb_10x_chunk        |  20358 |        20358 |               0 |                0 |     20358 |         20358 |       20358 |     1        |                0 |              19863 |            0 |
| NEPdb                  |  15912 |           58 |              15 |                1 |        74 |         15912 |         298 |     0.018728 |            15912 |                  0 |            0 |
| IEDB_tcell_api         |    500 |            0 |               0 |                0 |         0 |           330 |         500 |     1        |               43 |                377 |            0 |
| VDJdb_PDB_chunk        |    313 |          310 |               0 |                3 |       313 |           313 |         313 |     1        |                1 |                 96 |          313 |

## Missingness Policy

- No parsed source rows are silently discarded; non-peptidic IEDB rows are retained with `peptide_raw` populated and normalized `peptide` left empty.
- Paired alpha/beta, beta-only, alpha-only, and no-TCR rows are explicitly flagged.
- Public TCR evidence remains diagnostic unless labels and leakage boundaries are compatible with the downstream task.

## Parsed Primary Local Sources

- `/data/neoantigen_vaccine_hub/data_raw/tcr/vdjdb.txt`
- `/data/neoantigen_vaccine_hub/data_raw/tcr/mcpas/McPAS-TCR.csv`
- `/data/neoantigen_vaccine_hub/data_raw/nepdb/nepdb_all.csv`
- `/data/neoantigen_vaccine_hub/data_raw/iedb/iedb_tcell_human_positive.tsv`
- `/data/neoantigen_vaccine_hub/data_raw/iedb/iedb_tcell_human_cancer.tsv`
- VDJdb 10x and PDB chunks are parsed as paired/structure supplements when present.

## Outputs

- `tcr_registry.parquet` / `tcr_registry.tsv`
- `tcr_registry_source_counts.tsv`
- `tcr_registry_missingness.tsv`
- `tcr_registry_label_balance.tsv`
- `tcr_registry_source_inventory.tsv`
- `tcr_registry_partition_counts.tsv`
