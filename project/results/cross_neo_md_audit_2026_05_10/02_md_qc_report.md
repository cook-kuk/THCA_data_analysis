# MD QC Report

QC is descriptive. Short smoke/screen trajectories are startup evidence only; 10 ns pilots support qualitative stability screening, not immunogenicity proof.

| run_id | candidate | condition | frames | total_time_ps | trajectory_complete_fraction | nan_coordinates | peptide_rmsd_final_nm | peptide_com_drift_final_nm | qc_status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_mutant_TCR-pMHC_6UON_GADGVGKSAL_HLA-C0802_chainF_rep01_10p0ns | GADGVGKSAL/HLA-C*08:02 | explicit_cuda | 100 | 1e+04 | 1 | False | 0.1949 | 0.5873 | ok |
| prod_10ns_6VRN_1fs300K | HMTEVVRHC/HLA-A*02:01 | explicit_cuda_10ns | 200 | 1e+04 | 1 | False | 0.06511 | 0.4046 | ok |
| 6VRM_HMTEVVRHC_HLA-A0201_chainP_0p5ns | HMTEVVRHC/HLA-A*02:01 | explicit_cuda_replicate_screen | 10 | 500 | 1 | False | 0.03782 | 0.04359 | ok |
| prod_10ns_6VRN_1fs300K | HMTEVVRHC/HLA-A*02:01 | explicit_cuda_10ns | 200 | 1e+04 | 1 | False | 0.06511 | 0.4046 | ok |
| prod_10ns_6UON_2fs300K | GADGVGKSAL/HLA-C*08:02 | explicit_cuda_10ns | 200 | 1e+04 | 1 | False | 0.1388 | 0.6598 | ok |
| 6UON_GADGVGKSAL_HLA-C0802_chainC_1p0ns | GADGVGKSAL/HLA-C*08:02 | explicit_cuda_replicate_screen | 20 | 1000 | 1 | False | 0.1119 | 0.2434 | ok |
| 6UON_GADGVGKSAL | GADGVGKSAL/HLA-C*08:02 | vacuum_smoke | 10 | 0.1 | 1e-05 | False | 0.01004 | 0.004788 | ok |
| 6VRN_HMTEVVRHC | HMTEVVRHC/HLA-A*02:01 | vacuum_smoke | 10 | 0.1 | 1e-05 | False | 0.004762 | 0.002221 | ok |
