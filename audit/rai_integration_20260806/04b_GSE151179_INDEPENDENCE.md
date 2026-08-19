# AUDIT 04b — GSE151179 patient structure and independence

Generated 2026-08-06 by `scripts/audit/05_gse151179_patient_structure.py`

The uptake label is recorded once per patient, while the series contributes several specimens per patient. Sample-level rows therefore replicate the same exposure and are not independent. The first pass ran the primary test at patient level, which was correct, but reported sample-level sensitivity rows beside it without flagging the violation. Those rows are re-derived here.

| check                                   | value           |
|:----------------------------------------|:----------------|
| specimens (all tissue)                  | 52              |
| patients (all tissue)                   | 33              |
| tumour specimens                        | 39              |
| patients contributing tumour specimens  | 32              |
| patients with >1 tumour specimen        | 6               |
| max tumour specimens per patient        | 3               |
| patients with discordant uptake label   | 0               |
| mean tumour specimens per patient       | 1.219           |
| ICC of panel z within patient           | 0.0             |
| design effect                           | 1.0             |
| effective sample size (tumour)          | 39.0            |
| patient-level d (uptake Yes - No)       | 0.373           |
| patient-level 95% CI                    | -0.338 to 0.980 |
| patient-level Mann-Whitney P            | 0.534           |
| patient-level permutation P             | 0.3107          |
| sample-level naive SE (uptake)          | 0.1751          |
| sample-level cluster-robust SE (uptake) | 0.156           |
| SE inflation from clustering            | 0.891           |
| negative control patient-level d        | 0.653           |
| negative control patient-level P        | 0.1807          |

No patient carries discordant uptake labels, confirming the field is patient-level and that a within-patient (paired) comparison is impossible in this dataset. Only the between-patient contrast exists.
