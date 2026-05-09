
# CLEAN-NeoBench Split Audit

## Candidate Counts

- Candidates: 2715
- Labeled candidates: 2715
- Method score rows: 90292
- Metric rows: 7144

## Source Distribution

- `CEDAR`: 909
- `TESLA_mmc4`: 605
- `NEPdb`: 572
- `TESLA_mmc7_validation`: 310
- `ITSNdb_main`: 199
- `ITSNdb_Val`: 120

## Positive Counts by Source

| source_name           |   count |   sum |      mean |
|:----------------------|--------:|------:|----------:|
| CEDAR                 |     909 |   851 | 0.936194  |
| TESLA_mmc4            |     605 |    37 | 0.061157  |
| NEPdb                 |     572 |   151 | 0.263986  |
| TESLA_mmc7_validation |     310 |     4 | 0.0129032 |
| ITSNdb_main           |     199 |   129 | 0.648241  |
| ITSNdb_Val            |     120 |     7 | 0.0583333 |

## Overlap Flag Summary

| flag                                | value                     |   count |
|:------------------------------------|:--------------------------|--------:|
| exact_peptide_train_overlap         | True                      |    2275 |
| exact_peptide_train_overlap         | False                     |     440 |
| exact_peptide_hla_train_overlap     | True                      |    2274 |
| exact_peptide_hla_train_overlap     | False                     |     441 |
| near_peptide_train_overlap          | True                      |    2300 |
| near_peptide_train_overlap          | False                     |     415 |
| source_protein_window_train_overlap | False                     |    2715 |
| study_train_overlap                 | True                      |    2086 |
| study_train_overlap                 | False                     |     629 |
| patient_train_overlap               | False                     |    2715 |
| public_tool_training_overlap_any    | False                     |    2502 |
| public_tool_training_overlap_any    | True                      |     213 |
| public_tool_training_overlap_detail | no_inhouse_master_overlap |    2502 |
| public_tool_training_overlap_detail | in_master                 |     213 |
| leakage_risk_level                  | high                      |    2304 |
| leakage_risk_level                  | low                       |     388 |
| leakage_risk_level                  | medium                    |      23 |

## Leakage Risk Distribution

- `high`: 2304
- `low`: 388
- `medium`: 23

## Missing Metadata Summary

| column                    |   missing_or_empty |
|:--------------------------|-------------------:|
| patient_id                |               2715 |
| mutation_id               |               2715 |
| protein_id                |               2715 |
| wt_peptide                |               2715 |
| disease_context           |               2715 |
| cancer_type               |               2715 |
| gene                      |               2715 |
| source_protein_window     |               2715 |
| hla_loh                   |               2715 |
| cytolytic_score           |               2715 |
| tls_score                 |               2715 |
| immune_context_score      |               2715 |
| treatment_context         |               2715 |
| ifng_score                |               2715 |
| tumor_stage               |               2715 |
| antigen_processing_status |               2715 |
| mutant_expression         |               2715 |
| b2m_status                |               2715 |
| clonality                 |               2715 |
| vaf                       |               2715 |

## Split Contract Rows

- `hla_heldout`: 3039
- `supertype_heldout`: 2599
- `existing_prediction_context`: 416
- `source_heldout`: 295
- `study_heldout`: 295
- `overall_labeled`: 100
- `exact_peptide_hla_holdout`: 100
- `near_peptide_cluster_holdout`: 100
- `low_prevalence_heldout`: 100
- `korean_hla_focus`: 100
