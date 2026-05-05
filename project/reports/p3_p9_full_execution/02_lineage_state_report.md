# Lineage State Report

- Output: `/home/seungho/personal/THCA_data_analysis/project/results/p3_p9_full_execution/shared/lineage_state_master.tsv`
- Rows: 756
- Datasets: 5
- Platform pooling rule: no absolute cross-platform pooling; available score columns were retained or z-scored within dataset.
- Claim boundary: these are research-use lineage-state variables, not clinical treatment variables.

## Dataset Counts
- TCGA-THCA: 513
- GSE33630: 105
- GSE29265: 49
- GSE53157: 27
- GSE65144: 25

## Notes
- External processed score files were prioritized where present.
- TCGA lineage state was derived from already-local processed TCGA table columns, not from new protected data.
- Missing score columns are left as unavailable rather than imputed.
