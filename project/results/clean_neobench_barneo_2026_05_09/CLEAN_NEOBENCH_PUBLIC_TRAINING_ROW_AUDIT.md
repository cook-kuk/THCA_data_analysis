# CLEAN-NeoBench Public Training Row-Level Audit

## Purpose

This report audits whether benchmark candidates overlap public pretrained tool training corpora at row level. Documentation-level provenance is not enough. A method becomes a cleaner comparator only when the relevant training corpus is present and exact peptide-HLA/peptide overlap is audited.

## Corpus Files Discovered

- No public training corpus files found in configured search directories.

## Row-Level Summary

No row-level audit summary available.

## Example Overlap Rows

No overlap rows found or no corpora available.

## Required Corpus Inputs

| public_tool     | filename_tokens         | required_key                                                          | minimum_columns                                                              | clean_pass_rule                                                                               | drop_location                                            |
|:----------------|:------------------------|:----------------------------------------------------------------------|:-----------------------------------------------------------------------------|:----------------------------------------------------------------------------------------------|:---------------------------------------------------------|
| MHCflurry       | mhcflurry               | peptide + HLA allele preferred; peptide-only accepted as weaker audit | peptide plus HLA when available; source/split/label optional but recommended | no exact peptide-HLA and no exact peptide overlap; near-peptide overlap summarized separately | project/data/public_training_corpora/<tool>_training.tsv |
| NetMHCpan_4.1   | netmhcpan, netmhc       | peptide + HLA allele                                                  | peptide plus HLA when available; source/split/label optional but recommended | no exact peptide-HLA overlap against BA/EL training rows                                      | project/data/public_training_corpora/<tool>_training.tsv |
| BigMHC_IM       | bigmhc                  | peptide + HLA allele + split/source if available                      | peptide plus HLA when available; source/split/label optional but recommended | no exact peptide-HLA or exact peptide overlap against im_train/el_train rows                  | project/data/public_training_corpora/<tool>_training.tsv |
| PRIME           | prime                   | peptide + HLA allele preferred                                        | peptide plus HLA when available; source/split/label optional but recommended | no exact peptide-HLA and no exact peptide overlap against immunogenicity training rows        | project/data/public_training_corpora/<tool>_training.tsv |
| DeepImmuno      | deepimmuno, deephimmuno | peptide + HLA allele                                                  | peptide plus HLA when available; source/split/label optional but recommended | no exact peptide-HLA overlap against IEDB-derived training rows                               | project/data/public_training_corpora/<tool>_training.tsv |
| MHCnuggets_2    | mhcnuggets              | peptide + HLA allele                                                  | peptide plus HLA when available; source/split/label optional but recommended | no exact peptide-HLA overlap against IEDB-derived binding rows                                | project/data/public_training_corpora/<tool>_training.tsv |
| NetMHCstabpan   | netmhcstab              | peptide + HLA allele                                                  | peptide plus HLA when available; source/split/label optional but recommended | no exact peptide-HLA overlap against stability training rows                                  | project/data/public_training_corpora/<tool>_training.tsv |
| TransPHLA       | transphla               | peptide + HLA allele                                                  | peptide plus HLA when available; source/split/label optional but recommended | no exact peptide-HLA overlap against pHLA training rows                                       | project/data/public_training_corpora/<tool>_training.tsv |
| TSCAPE_TITANiAN | tscape, titanian        | peptide + HLA allele + TCR fields when available                      | peptide plus HLA when available; source/split/label optional but recommended | no exact peptide-HLA overlap; TCR/pMHC train rows audited separately when available           | project/data/public_training_corpora/<tool>_training.tsv |

## Warnings

- None.

## Claim Rule

If no corpus file is available for a public pretrained tool, its training overlap status remains unresolved. If overlap is found, the method remains caveated. If a method-specific corpus is present and no overlap is found, it can be marked row-audited for this benchmark, subject to the corpus being complete and version-matched.

