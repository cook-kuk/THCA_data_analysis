# MD Live Status Report

This report parses OpenMM `state.tsv` files. It does not infer immunogenicity.

## Latest Status

| run_id | candidate | condition | time_ps | target_ns | completion_fraction | temperature_k | potential_energy_kj_mol | speed_ns_per_day | complete_by_metadata |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_mutant_TCR-pMHC_6UON_GADGVGKSAL_HLA-C0802_chainF_rep01_10p0ns | GADGVGKSAL/HLA-C*08:02 | explicit_cuda | 1e+04 | 10 | 1 | 301 | -5.361e+06 | 44.1 | True |
| prod_10ns_6VRN_1fs300K | HMTEVVRHC/HLA-A*02:01 | explicit_cuda_10ns | 1e+04 | 10 | 1 | 299.9 | -5.875e+06 | 21 | True |
| 6VRM_HMTEVVRHC_HLA-A0201_chainP_0p5ns | HMTEVVRHC/HLA-A*02:01 | explicit_cuda_replicate_screen | 500 | 0.5 | 1 | 300.4 | -5.388e+06 | 22.4 | True |
| prod_10ns_6VRN_1fs300K | HMTEVVRHC/HLA-A*02:01 | explicit_cuda_10ns | 1e+04 | 10 | 1 | 299.9 | -5.875e+06 | 21 | True |
| prod_10ns_6UON_2fs300K | GADGVGKSAL/HLA-C*08:02 | explicit_cuda_10ns | 1e+04 | 10 | 1 | 300 | -5.332e+06 | 41.9 | True |
| 6UON_GADGVGKSAL_HLA-C0802_chainC_1p0ns | GADGVGKSAL/HLA-C*08:02 | explicit_cuda_replicate_screen | 1000 | 1 | 1 | 301.1 | -5.154e+06 | 43.1 | True |
| 6UON_GADGVGKSAL | GADGVGKSAL/HLA-C*08:02 | vacuum_smoke | 0.1 | 10 | 1e-05 | 28.07 | -8.517e+04 | 0.417 | True |
| 6VRN_HMTEVVRHC | HMTEVVRHC/HLA-A*02:01 | vacuum_smoke | 0.1 | 10 | 1e-05 | 30.96 | -8.609e+04 | 0.267 | True |

## Sanity Checks

- Temperature is expected to stay near 300 K for production/pilot runs.
- Energy should not show abrupt unbounded growth.
- Missing trajectory files are treated as partial or state-only runs.
