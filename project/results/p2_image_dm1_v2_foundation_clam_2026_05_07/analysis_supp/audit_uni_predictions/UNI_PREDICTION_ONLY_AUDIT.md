# UNI-final prediction-only audit

This is an OOF-table audit only. It does not replace UNI LOTO retraining.

| Metric | Value |
|---|---:|
| n slides | 54 |
| overall AUC | 0.874 |
| clinical-only OOF AUC | 0.735 |
| TSS-only within-RAS-like AUC | 1.000 |
| molecular_subtype=RAS_like AUC | 0.795 |
| histology_subtype=FVPTC AUC | 0.795 |
| molecular_subtype=BRAF_like AUC | 0.891 |
| histology_subtype=cPTC AUC | 0.891 |
| sex=Female AUC | 0.893 |
| sex=Male AUC | 0.812 |

## Interpretation

UNI-final OOF predictions improve overall performance relative to the ViT-L audit table and do not reproduce the RAS-like AUC=1.000. However, the RAS-like TSS-label structure is unchanged: DJ/FK contain the 3 RAS-like DM1 cases and EM contains 10 RAS-like DM2 cases. Therefore K2 H&E or UNI-specific LOTO retraining is still required before any center-generalized image-DM1 subgroup claim.
