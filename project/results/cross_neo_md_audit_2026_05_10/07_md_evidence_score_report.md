# MD Evidence Score Report

Weights are transparent heuristics: pMHC stability, TCR contacts, runtime/QC, replicate consistency, and counterfactual specificity. Counterfactual specificity is zero until WT/decoy trajectories exist.

| run_id | candidate | MD_evidence_score | MD_evidence_label | pMHC_stability_score | TCR_recognition_score | simulation_qc_score |
| --- | --- | --- | --- | --- | --- | --- |
| tier1_full_10ns_replicate_or_control_GADGVGKSAL_HLA-C_08_02_mutant_TCR-pMHC_6UON_GADGVGKSAL_HLA-C0802_chainF_rep01_10p0ns | GADGVGKSAL/HLA-C*08:02 | 0.5773 | MD_MODERATE | 0.4666 | 0.3704 | 1 |
| prod_10ns_6VRN_1fs300K | HMTEVVRHC/HLA-A*02:01 | 0.8218 | MD_VERY_STRONG | 0.6626 | 1 | 1 |
| 6VRM_HMTEVVRHC_HLA-A0201_chainP_0p5ns | HMTEVVRHC/HLA-A*02:01 | 0.9158 | MD_VERY_STRONG | 0.9101 | 1 | 1 |
| prod_10ns_6VRN_1fs300K | HMTEVVRHC/HLA-A*02:01 | 0.8218 | MD_VERY_STRONG | 0.6626 | 1 | 1 |
| prod_10ns_6UON_2fs300K | GADGVGKSAL/HLA-C*08:02 | 0.6166 | MD_MODERATE | 0.5203 | 0.4404 | 1 |
| 6UON_GADGVGKSAL_HLA-C0802_chainC_1p0ns | GADGVGKSAL/HLA-C*08:02 | 0.7101 | MD_STRONG | 0.7212 | 0.504 | 1 |
| 6UON_GADGVGKSAL | GADGVGKSAL/HLA-C*08:02 | 0.6215 | MD_INSUFFICIENT_RUNTIME | 0.9886 | 0.54 | 1e-05 |
| 6VRN_HMTEVVRHC | HMTEVVRHC/HLA-A*02:01 | 0.7455 | MD_INSUFFICIENT_RUNTIME | 0.9882 | 1 | 1e-05 |
