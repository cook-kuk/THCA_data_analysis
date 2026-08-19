# CLEAN-NeoBench Public Overlap Audit Intake KR

## 한 줄 결론

public tool을 clean comparator로 올리려면 row-level training corpus 파일이 필요하다. 현재는 파일이 0개라 public clean comparator allowed는 0이 맞다.

## 어디에 넣나

`/home/seungho/personal/THCA_data_analysis/project/data/public_training_corpora`

## 무엇을 넣나

| public_tool     | expected_drop_file                                                                                          | template_file                                                          | required_key                                                          | clean_pass_rule                                                                               | current_disposition             |
|:----------------|:------------------------------------------------------------------------------------------------------------|:-----------------------------------------------------------------------|:----------------------------------------------------------------------|:----------------------------------------------------------------------------------------------|:--------------------------------|
| MHCflurry       | /home/seungho/personal/THCA_data_analysis/project/data/public_training_corpora/MHCflurry_training.tsv       | public_tool_overlap_audit_intake/MHCflurry_training_template.tsv       | peptide + HLA allele preferred; peptide-only accepted as weaker audit | no exact peptide-HLA and no exact peptide overlap; near-peptide overlap summarized separately | caveated_public_comparator_only |
| NetMHCpan_4.1   | /home/seungho/personal/THCA_data_analysis/project/data/public_training_corpora/NetMHCpan_4.1_training.tsv   | public_tool_overlap_audit_intake/NetMHCpan_4.1_training_template.tsv   | peptide + HLA allele                                                  | no exact peptide-HLA overlap against BA/EL training rows                                      | caveated_public_comparator_only |
| BigMHC_IM       | /home/seungho/personal/THCA_data_analysis/project/data/public_training_corpora/BigMHC_IM_training.tsv       | public_tool_overlap_audit_intake/BigMHC_IM_training_template.tsv       | peptide + HLA allele + split/source if available                      | no exact peptide-HLA or exact peptide overlap against im_train/el_train rows                  | caveated_public_comparator_only |
| PRIME           | /home/seungho/personal/THCA_data_analysis/project/data/public_training_corpora/PRIME_training.tsv           | public_tool_overlap_audit_intake/PRIME_training_template.tsv           | peptide + HLA allele preferred                                        | no exact peptide-HLA and no exact peptide overlap against immunogenicity training rows        | caveated_public_comparator_only |
| DeepImmuno      | /home/seungho/personal/THCA_data_analysis/project/data/public_training_corpora/DeepImmuno_training.tsv      | public_tool_overlap_audit_intake/DeepImmuno_training_template.tsv      | peptide + HLA allele                                                  | no exact peptide-HLA overlap against IEDB-derived training rows                               | caveated_public_comparator_only |
| MHCnuggets_2    | /home/seungho/personal/THCA_data_analysis/project/data/public_training_corpora/MHCnuggets_2_training.tsv    | public_tool_overlap_audit_intake/MHCnuggets_2_training_template.tsv    | peptide + HLA allele                                                  | no exact peptide-HLA overlap against IEDB-derived binding rows                                | caveated_public_comparator_only |
| NetMHCstabpan   | /home/seungho/personal/THCA_data_analysis/project/data/public_training_corpora/NetMHCstabpan_training.tsv   | public_tool_overlap_audit_intake/NetMHCstabpan_training_template.tsv   | peptide + HLA allele                                                  | no exact peptide-HLA overlap against stability training rows                                  | caveated_public_comparator_only |
| TransPHLA       | /home/seungho/personal/THCA_data_analysis/project/data/public_training_corpora/TransPHLA_training.tsv       | public_tool_overlap_audit_intake/TransPHLA_training_template.tsv       | peptide + HLA allele                                                  | no exact peptide-HLA overlap against pHLA training rows                                       | caveated_public_comparator_only |
| TSCAPE_TITANiAN | /home/seungho/personal/THCA_data_analysis/project/data/public_training_corpora/TSCAPE_TITANiAN_training.tsv | public_tool_overlap_audit_intake/TSCAPE_TITANiAN_training_template.tsv | peptide + HLA allele + TCR fields when available                      | no exact peptide-HLA overlap; TCR/pMHC train rows audited separately when available           | caveated_public_comparator_only |

## 최소 컬럼

`peptide`, `hla` 또는 `hla_allele_4digit`, `public_tool`이 핵심이다. 가능하면 `source_dataset`, `train_split`, `label`, `assay_type`, `publication_or_url`도 넣는다.

## 다시 실행

```bash
python scripts/clean_neobench/run_clean_neobench_pipeline.py \
  --repo-root . \
  --output-root project/results/clean_neobench_barneo_2026_05_09
```

## 정책

- exact peptide-HLA overlap 있으면 clean comparator 불가.
- exact peptide overlap이 unresolved이면 clean comparator 불가.
- near-peptide overlap은 별도 caveat로 보고한다.
- public tool이 점수상 이겨도 audit 전에는 caveated comparator다.
