# CLEAN-NeoBench Public Tool Overlap Audit

## Summary

This audit scaffold separates documentation-level provenance from row-level training-corpus overlap auditing.

- Public pretrained methods detected: 9
- Internal/local methods detected: 91
- Public methods promoted to clean comparators: 0
- Row-level clean-pass status: unresolved unless a future method-specific candidate/peptide-HLA audit table proves otherwise.

## Public Method Disposition

| method_name     | method_role                | training_overlap_audit_status           | clean_comparator_allowed_after_audit   | reviewer_disposition            | caveat                                                                                   |
|:----------------|:---------------------------|:----------------------------------------|:---------------------------------------|:--------------------------------|:-----------------------------------------------------------------------------------------|
| BigMHC_IM       | caveated_public_comparator | unresolved_no_row_level_training_corpus | False                                  | caveated_public_comparator_only | public deep presentation/immunogenicity model; row-level training overlap unresolved     |
| DeepImmuno      | caveated_public_comparator | unresolved_no_row_level_training_corpus | False                                  | caveated_public_comparator_only | public immunogenicity model; row-level training overlap unresolved                       |
| MHCflurry       | caveated_public_comparator | unresolved_no_row_level_training_corpus | False                                  | caveated_public_comparator_only | public pretrained MHC-I presentation/binding tool; row-level training overlap unresolved |
| MHCnuggets_2    | caveated_public_comparator | unresolved_no_row_level_training_corpus | False                                  | caveated_public_comparator_only | public MHC binding predictor; row-level training overlap unresolved                      |
| NetMHCpan_4.1   | caveated_public_comparator | unresolved_no_row_level_training_corpus | False                                  | caveated_public_comparator_only | public pretrained NetMHC-family tool; row-level training overlap unresolved              |
| NetMHCstabpan   | caveated_public_comparator | unresolved_no_row_level_training_corpus | False                                  | caveated_public_comparator_only | public pretrained NetMHC-family tool; row-level training overlap unresolved              |
| PRIME           | caveated_public_comparator | unresolved_no_row_level_training_corpus | False                                  | caveated_public_comparator_only | public immunogenicity predictor; row-level training overlap unresolved                   |
| TSCAPE_TITANiAN | caveated_public_comparator | unresolved_no_row_level_training_corpus | False                                  | caveated_public_comparator_only | public multidomain immunogenicity/TCR model; row-level training overlap unresolved       |
| TransPHLA       | caveated_public_comparator | unresolved_no_row_level_training_corpus | False                                  | caveated_public_comparator_only | public pHLA binding transformer; row-level training overlap unresolved                   |

## Audit Artifacts Discovered

- `/home/seungho/personal/THCA_data_analysis/project/results/cross_neo_v0/public_overlap_audit.tsv`
- `/home/seungho/personal/THCA_data_analysis/project/results/cross_neo_v1/public_overlap_audit_v1.tsv`
- `/home/seungho/personal/THCA_data_analysis/project/results/cross_neo_v1_lockdown/public_overlap_audit_if_available.tsv`
- `/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09/method_training_audit_2026_05_09/method_training_audit.tsv`
- `/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09/wave11/wave11_overlap_audit.tsv`

## Rule

Public pretrained outputs may be reported as caveated comparators. They must not be used as clean training features or described as clean external baselines until row-level training-corpus overlap is audited.

## Forbidden Claims

- Public tools are clean baselines without overlap audit.
- External validation is proven by public pretrained tool agreement.
- A public pretrained score can rescue a high-leakage candidate into a clean benchmark claim.
