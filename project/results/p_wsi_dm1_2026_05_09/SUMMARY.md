# WSI DM1/DM2 RunPod Recovery Summary

Date: 2026-05-09

## Inputs

- Manifest: `project/data/manifests/tcga_thca_wsi_dm_balanced.tsv`
- RunPod split artifacts: `project/data/processed/TCGA-THCA-WSI-DM/runpod_split_2026_05_09_v2/`
- Backbone: DINOv2 ViT-L 518 fallback; UNI unavailable due gated Hugging Face access.

## RunPod Outcome

- `uvp9i2r9s6l85y` / `thca-spark-dm-a6000-v4`: completed first 25 WSI split, outputs pulled, pod stopped.
- `r8v801csxxscw0` / `thca-neo-bayesian-aux`: completed second 25 WSI split, outputs pulled, pod stopped.
- `ytnsyachant0a4` / `thca-img-dm1-a6000`: remained EXITED; host GPU unavailable for resume.

## Key Fix

The first Phase 3 pass used level dimensions as if they were level-0 coordinates, which sampled only an upper-left region and produced 8 zero-tile slides. The v2 pass converts coordinates by `slide.level_downsamples[level]`, retiled all slides, and reached 50/50 embedding coverage.

## V2 Metrics

- Manifest slides: 50
- Embedded slides: 50
- Missing embeddings: 0
- Class counts: DM1=25, DM2=25
- Merged case-grouped LOSO AUC using DM2 probability: 0.7456
- Accuracy at 0.5: 0.72

## Outputs

- Dossier: `project/papers_hub_2026_05_04/wsi_dm1_recovery_dossier_2026_05_09.html`
- Live page: `/var/www/papers/papers_hub_2026_05_04/wsi_dm1_recovery_dossier_2026_05_09.html`
- Figures: `fig_wsi_dm1_v2_roc.{png,pdf}`, `fig_wsi_dm1_v2_score_distribution.{png,pdf}`
- Tables: `wsi_dm1_v2_predictions.tsv`, `wsi_dm1_v2_per_bucket_summary.tsv`, `wsi_dm1_v2_top_error_slides.tsv`
- `project/data/processed/TCGA-THCA-WSI-DM/runpod_split_2026_05_09_v2/merged/merged_loso_metrics.json`
- `project/data/processed/TCGA-THCA-WSI-DM/runpod_split_2026_05_09_v2/merged/merged_dm1_vs_dm2_loso_predictions.tsv`
- `project/data/processed/TCGA-THCA-WSI-DM/runpod_split_2026_05_09_v2/merged/embedding_coverage.tsv`
- `project/data/processed/TCGA-THCA-WSI-DM/runpod_split_2026_05_09_v2/merged/missing_embedding_slides.tsv`
- Per-pod archives: `pod_20779/dm_wsi_artifacts.tgz`, `pod_20527/dm_wsi_artifacts.tgz`

## Caveat

This is a recovered 50-WSI foundation-embedding pilot, not a manuscript-ready pathology claim. The corrected v2 result supports a moderate image-DM1 signal, but it should be treated as reserve or methods/audit evidence unless a larger WSI run is explicitly authorized.
