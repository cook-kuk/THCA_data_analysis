# TCGA-BRCA pathology-AI audit input status

**Proposed task:** candidate binary pathology task to lock after metadata assembly (e.g. basal-like vs luminal or ER-negative vs ER-positive)

**Ready for 12-test audit:** `False`

## Required inputs

| | |   | I | n | p | u | t |   | | |   | E | x | p | e | c | t | e | d |   | p | a | t | h |   | | |   | S | t | a | t | u | s |   | | |
| | | - | - | - | | | - | - | - | | | - | - | - | | |
| Feature dir | /home/seungho/personal/THCA_data_analysis/project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/phase2_tcga_brca_UNI/features | 0 .pt files |
| Manifest | /home/seungho/personal/THCA_data_analysis/project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/phase2_tcga_brca_UNI/slide_manifest.tsv | missing |
| Metadata | /home/seungho/personal/THCA_data_analysis/project/metadata/tcga_brca_audit_metadata.tsv | missing |
| OOF predictions | /home/seungho/personal/THCA_data_analysis/project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/phase2_tcga_brca_UNI/clam_per_slide_predictions.tsv | missing |

## Blockers

- missing per-slide UNI .pt feature files
- missing slide_manifest.tsv
- missing cohort audit metadata TSV
- missing OOF predictions; needed after first CLAM run, not for input staging

## Expected manifest schema

`file_id`, `submitter_id`, `label`, `tss`; optional but recommended: `histology`, `sex`, `subtype`, `age`, `stage`, `site`.

## Next command after feature extraction

```bash
python3 /home/seungho/personal/THCA_data_analysis/project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/scripts/audit_tcga_brca.py
```
