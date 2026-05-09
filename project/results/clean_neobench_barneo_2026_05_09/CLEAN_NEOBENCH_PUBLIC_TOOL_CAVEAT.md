
# CLEAN-NeoBench Public Tool Caveat

## Included Public or Public-Pretrained Tools

| method_name     | method_family       | method_role                | uses_public_pretraining   | training_overlap_audited   | clean_comparator_allowed   | caveat                                                              |
|:----------------|:--------------------|:---------------------------|:--------------------------|:---------------------------|:---------------------------|:--------------------------------------------------------------------|
| BigMHC_IM       | deep_immunogenicity | caveated_public_comparator | True                      | False                      | False                      | public pretrained comparator; row-level training overlap unresolved |
| DeepImmuno      | deep_immunogenicity | caveated_public_comparator | True                      | False                      | False                      | public pretrained comparator; row-level training overlap unresolved |
| MHCflurry       | deep_immunogenicity | caveated_public_comparator | True                      | False                      | False                      | public pretrained comparator; row-level training overlap unresolved |
| MHCnuggets_2    | deep_immunogenicity | caveated_public_comparator | True                      | False                      | False                      | public pretrained comparator; row-level training overlap unresolved |
| NetMHCpan_4.1   | deep_immunogenicity | caveated_public_comparator | True                      | False                      | False                      | public pretrained comparator; row-level training overlap unresolved |
| NetMHCstabpan   | deep_immunogenicity | caveated_public_comparator | True                      | False                      | False                      | public pretrained comparator; row-level training overlap unresolved |
| PRIME           | deep_immunogenicity | caveated_public_comparator | True                      | False                      | False                      | public pretrained comparator; row-level training overlap unresolved |
| TSCAPE_TITANiAN | deep_immunogenicity | caveated_public_comparator | True                      | False                      | False                      | public pretrained comparator; row-level training overlap unresolved |
| TransPHLA       | deep_immunogenicity | caveated_public_comparator | True                      | False                      | False                      | public pretrained comparator; row-level training overlap unresolved |

## Caveat Rule

Public pretrained tools are caveated comparators until their official training corpora are row-audited against benchmark candidates. A no-overlap flag against our local master table is not the same as a no-overlap guarantee against MHCflurry, NetMHCpan, BigMHC, PRIME, MixMHCpred, NetMHCstabpan, or other public training data.

## Allowed Comparisons

- Public tools may be shown as caveated upper-bound or field-context comparators.
- Internal models may be compared against them only with the training-overlap caveat visible.
- Clean internal anchor status is reserved for locally audited methods such as Structure_LR.

## Forbidden Claims

- Public pretrained tools are clean external baselines without row-level training-corpus audit.
- A public tool win proves external validation.
- A public score can be used as a clean training feature before overlap audit.
