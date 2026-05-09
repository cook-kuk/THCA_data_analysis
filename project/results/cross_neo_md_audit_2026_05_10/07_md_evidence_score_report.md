# MD Evidence Score Report

Weights are transparent heuristics: pMHC stability, TCR contacts, runtime/QC, replicate consistency, and counterfactual specificity. Counterfactual specificity is zero until WT/decoy trajectories exist.

| run_id | candidate | MD_evidence_score | MD_evidence_label | pMHC_stability_score | TCR_recognition_score | simulation_qc_score |
| --- | --- | --- | --- | --- | --- | --- |
| prod_10ns_6UON_2fs300K | GADGVGKSAL/HLA-C*08:02 | 0.5566 | MD_MODERATE | 0.5203 | 0.4404 | 1 |
| 6UON_GADGVGKSAL_HLA-C0802_chainC_1p0ns | GADGVGKSAL/HLA-C*08:02 | 0.6501 | MD_MODERATE | 0.7212 | 0.504 | 1 |
| 6UON_GADGVGKSAL | GADGVGKSAL/HLA-C*08:02 | 0.5615 | MD_INSUFFICIENT_RUNTIME | 0.9886 | 0.54 | 1e-05 |
| 6VRN_HMTEVVRHC | HMTEVVRHC/HLA-A*02:01 | 0.6655 | MD_INSUFFICIENT_RUNTIME | 0.9882 | 1 | 1e-05 |
