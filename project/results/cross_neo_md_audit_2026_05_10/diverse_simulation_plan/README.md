# Diverse MD Simulation Launch Pack

## Boundary

This is a GPU launch scaffold for structural robustness testing. It does not prove immunogenicity, cancer specificity, vaccine efficacy, or SOTA superiority.

## Summary

- total planned jobs: 17
- immediate ready jobs: 5
- runnable total jobs: 17
- immediate requested GPU trajectory length: 50.0 ns

## Immediate Jobs

| job_id | peptide | hla_4digit | control_type | complex_scope | ns | random_seed | chain_quality_flag |
| --- | --- | --- | --- | --- | --- | --- | --- |
| tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_mutant_TCR-pMHC_6UON_GADGVGKSAL_HLA-C0802_chainF_rep01_10p0ns | GADGVGKSAL | HLA-C*08:02 | mutant | TCR-pMHC | 10 | 21011 | ok |
| tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_mutant_TCR-pMHC_6UON_GADGVGKSAL_HLA-C0802_chainF_rep02_10p0ns | GADGVGKSAL | HLA-C*08:02 | mutant | TCR-pMHC | 10 | 21012 | ok |
| tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_same_or_similar_hla_positive_control_pMHC_6ULR_GADGVGKSA_HLA-C0802_chainC_rep01_10p0ns | GADGVGKSAL | HLA-C*08:02 | same_or_similar_hla_positive_control | pMHC | 10 | 21021 | tcr_chain_count_lt_2 |
| tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_same_or_similar_hla_positive_control_pMHC_6ULR_GADGVGKSA_HLA-C0802_chainC_rep02_10p0ns | GADGVGKSAL | HLA-C*08:02 | same_or_similar_hla_positive_control | pMHC | 10 | 21022 | tcr_chain_count_lt_2 |
| tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_same_or_similar_hla_positive_control_pMHC_6ULR_GADGVGKSA_HLA-C0802_chainC_rep03_10p0ns | GADGVGKSAL | HLA-C*08:02 | same_or_similar_hla_positive_control | pMHC | 10 | 21023 | tcr_chain_count_lt_2 |

## Hold / After-Sync Jobs

| job_id | peptide | hla_4digit | control_type | complex_scope | ns | hold_reason | chain_quality_flag |
| --- | --- | --- | --- | --- | --- | --- | --- |
| tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_6VRM_HMTEVVRHC_HLA-A0201_chainP_rep01_0p5ns | HMTEVVRHC | HLA-A*02:01 | mutant_alt_template | TCR-pMHC | 0.5 | finish_or_sync_current_HMTEVVRHC_10ns_before_more_HMTEVVRHC_GPU_work | ok |
| tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_6VQO_HMTEVVRHC_HLA-A0201_chainP_rep01_0p5ns | HMTEVVRHC | HLA-A*02:01 | mutant_alt_template | TCR-pMHC | 0.5 | finish_or_sync_current_HMTEVVRHC_10ns_before_more_HMTEVVRHC_GPU_work | ok |
| tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_6VQO_HMTEVVRHC_HLA-A0201_chainQ_rep01_0p5ns | HMTEVVRHC | HLA-A*02:01 | mutant_alt_template | TCR-pMHC | 0.5 | finish_or_sync_current_HMTEVVRHC_10ns_before_more_HMTEVVRHC_GPU_work | ok |
| tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_7RM4_HMTEVVRHC_HLA-A0201_chainH_rep01_0p5ns | HMTEVVRHC | HLA-A*02:01 | mutant_alt_template | TCR-pMHC | 0.5 | finish_or_sync_current_HMTEVVRHC_10ns_before_more_HMTEVVRHC_GPU_work | ok |
| tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_7RM4_HMTEVVRHC_HLA-A0201_chainC_rep01_0p5ns | HMTEVVRHC | HLA-A*02:01 | mutant_alt_template | TCR-pMHC | 0.5 | finish_or_sync_current_HMTEVVRHC_10ns_before_more_HMTEVVRHC_GPU_work | ok |
| tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_7RM4_HMTEVVRHC_HLA-A0201_chainR_rep01_0p5ns | HMTEVVRHC | HLA-A*02:01 | mutant_alt_template | TCR-pMHC | 0.5 | finish_or_sync_current_HMTEVVRHC_10ns_before_more_HMTEVVRHC_GPU_work | ok |
| tier0_short_template_diversity_screen_HMTEVVRHC_HLA-A_02_01_mutant_alt_template_TCR-pMHC_7RM4_HMTEVVRHC_HLA-A0201_chainM_rep01_0p5ns | HMTEVVRHC | HLA-A*02:01 | mutant_alt_template | TCR-pMHC | 0.5 | finish_or_sync_current_HMTEVVRHC_10ns_before_more_HMTEVVRHC_GPU_work | ok |
| tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_mutant_TCR-pMHC_6VRN_HMTEVVRHC_HLA-A0201_chainP_rep01_10p0ns | HMTEVVRHC | HLA-A*02:01 | mutant | TCR-pMHC | 10 | finish_or_sync_current_HMTEVVRHC_10ns_before_duplicate_or_control_batch | ok |
| tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_mutant_TCR-pMHC_6VRN_HMTEVVRHC_HLA-A0201_chainP_rep02_10p0ns | HMTEVVRHC | HLA-A*02:01 | mutant | TCR-pMHC | 10 | finish_or_sync_current_HMTEVVRHC_10ns_before_duplicate_or_control_batch | ok |
| tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_same_or_similar_hla_positive_control_TCR-pMHC_5NQK_ELAGIGILTV_HLA-A0201_chainP_rep01_10p0ns | HMTEVVRHC | HLA-A*02:01 | same_or_similar_hla_positive_control | TCR-pMHC | 10 | finish_or_sync_current_HMTEVVRHC_10ns_before_duplicate_or_control_batch | ok |
| tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_same_or_similar_hla_positive_control_TCR-pMHC_5NQK_ELAGIGILTV_HLA-A0201_chainP_rep02_10p0ns | HMTEVVRHC | HLA-A*02:01 | same_or_similar_hla_positive_control | TCR-pMHC | 10 | finish_or_sync_current_HMTEVVRHC_10ns_before_duplicate_or_control_batch | ok |
| tier1_full_10ns_replicate_or_control_HMTEVVRHC_HLA-A_02_01_same_or_similar_hla_positive_control_TCR-pMHC_5NQK_ELAGIGILTV_HLA-A0201_chainP_rep03_10p0ns | HMTEVVRHC | HLA-A*02:01 | same_or_similar_hla_positive_control | TCR-pMHC | 10 | finish_or_sync_current_HMTEVVRHC_10ns_before_duplicate_or_control_batch | ok |

## Blocked Work

WT and anchor-preserved decoy simulations remain structure/sequence-curation work. See `blocked_wt_decoy_structure_work.tsv`.

## Run Commands

On a CUDA/OpenMM RunPod after upload:

```bash
bash run_immediate_ready_batch.sh
bash run_after_hmtevv_sync.sh
```

For upload/start, edit environment variables if the pod endpoint changes, then run:

```bash
bash runpod_upload_and_start_template.sh
```
