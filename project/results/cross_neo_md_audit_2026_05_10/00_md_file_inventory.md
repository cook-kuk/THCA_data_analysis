# MD File Inventory

Inventory of locally available and RunPod-synced OpenMM outputs. Remote runs may still be partial if the trajectory was still running when synced.

| run_id | candidate | condition | trajectory_bytes | latest_time_ps | target_ns | completion_fraction | run_complete | missing_files |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| prod_10ns_6VRN_1fs300K | HMTEVVRHC/HLA-A*02:01 | explicit_cuda_10ns | 0 | 6600 | 10 | 0.66 | False | trajectory,metadata |
| prod_10ns_6UON_2fs300K | GADGVGKSAL/HLA-C*08:02 | explicit_cuda_10ns | 763257076 | 1e+04 | 10 | 1 | True |  |
| 6UON_GADGVGKSAL_HLA-C0802_chainC_1p0ns | GADGVGKSAL/HLA-C*08:02 | explicit_cuda_replicate_screen | 73828756 | 1000 | 1 | 1 | True |  |
| 6UON_GADGVGKSAL | GADGVGKSAL/HLA-C*08:02 | vacuum_smoke | 1533876 | 0.1 | 10 | 1e-05 | True |  |
| 6VRN_HMTEVVRHC | HMTEVVRHC/HLA-A*02:01 | vacuum_smoke | 1533516 | 0.1 | 10 | 1e-05 | True |  |
