# RAW_DUMP_for_external_review

## SEC 0. 실행 메타

- orchestrator.log absolute path: `/home/seungho/personal/THCA_data_analysis/project/results/v17p3/orchestrator.log`
```text
[2026-04-27 01:52:27] start v17p3_F2_msigdb_proper.py
[2026-04-27 01:52:27] start v17p3_F1_external_recovery.py
[2026-04-27 01:52:27] start v17p3_F3_clinical_reframe.py
[2026-04-27 01:52:28] v17p3 F3 start
[2026-04-27 01:52:28] v17p3 F1 start
[2026-04-27 01:52:28] v17p3 F2 start
[2026-04-27 01:52:32] v17p3 F3 done {'n_endpoints_p_lt_0_05': 4, 'best_endpoint': 'rai_score_recalc'}
[2026-04-27 01:52:32] v17p3_F3_clinical_reframe.py stdout [2026-04-27 01:52:28] v17p3 F3 start
[2026-04-27 01:52:32] v17p3 F3 done {'n_endpoints_p_lt_0_05': 4, 'best_endpoint': 'rai_score_recalc'}

[2026-04-27 01:52:32] end v17p3_F3_clinical_reframe.py code=0 sec=5.6
[2026-04-27 01:54:04] v17p3 F1 done {'gse213647_gene_overlap_tiera': 55, 'robust_cohort_count_dia_auc_gt_0_85': 2}
[2026-04-27 01:54:05] v17p3_F1_external_recovery.py stdout [2026-04-27 01:52:28] v17p3 F1 start
[2026-04-27 01:54:04] v17p3 F1 done {'gse213647_gene_overlap_tiera': 55, 'robust_cohort_count_dia_auc_gt_0_85': 2}

[2026-04-27 01:54:05] end v17p3_F1_external_recovery.py code=0 sec=98.0
[2026-04-27 01:55:32] v17p3_F2_msigdb_proper.py stdout [2026-04-27 01:52:28] v17p3 F2 start

[2026-04-27 01:55:32] v17p3_F2_msigdb_proper.py stderr ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/seungho/personal/THCA_data_analysis/project/.venv/lib/python3.12/site-packages/pandas/core/frame.py", line 10859, in merge
    return merge(
           ^^^^^^
  File "/home/seungho/personal/THCA_data_analysis/project/.venv/lib/python3.12/site-packages/pandas/core/reshape/merge.py", line 170, in merge
    op = _MergeOperation(
         ^^^^^^^^^^^^^^^^
  File "/home/seungho/personal/THCA_data_analysis/project/.venv/lib/python3.12/site-packages/pandas/core/reshape/merge.py", line 794, in __init__
    ) = self._get_merge_keys()
        ^^^^^^^^^^^^^^^^^^^^^^
  File "/home/seungho/personal/THCA_data_analysis/project/.venv/lib/python3.12/site-packages/pandas/core/reshape/merge.py", line 1311, in _get_merge_keys
    left_keys.append(left._get_label_or_level_values(lk))
                     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/seungho/personal/THCA_data_analysis/project/.venv/lib/python3.12/site-packages/pandas/core/generic.py", line 1928, in _get_label_or_level_values
    raise ValueError(
ValueError: The column label 'pathway' is not unique.

[2026-04-27 01:55:32] end v17p3_F2_msigdb_proper.py code=1 sec=185.8
[2026-04-27 01:55:32] start v17p3_A1_scrna_heterogeneity.py
[2026-04-27 01:55:32] start v17p3_A3_drug_response.py
[2026-04-27 01:55:32] start v17p3_A2_dm_in_braf_ras.py
[2026-04-27 01:55:34] v17p3 A3 start
[2026-04-27 01:55:34] v17p3 A2 start
[2026-04-27 01:55:35] v17p3 A1 start
[2026-04-27 01:55:38] v17p3_A1_scrna_heterogeneity.py stdout [2026-04-27 01:55:35] v17p3 A1 start

[2026-04-27 01:55:38] v17p3_A1_scrna_heterogeneity.py stderr  "/home/seungho/personal/THCA_data_analysis/project/.venv/lib/python3.12/site-packages/pandas/core/frame.py", line 10784, in join
    return merge(
           ^^^^^^
  File "/home/seungho/personal/THCA_data_analysis/project/.venv/lib/python3.12/site-packages/pandas/core/reshape/merge.py", line 184, in merge
    return op.get_result(copy=copy)
           ^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/seungho/personal/THCA_data_analysis/project/.venv/lib/python3.12/site-packages/pandas/core/reshape/merge.py", line 888, in get_result
    result = self._reindex_and_concat(
             ^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/seungho/personal/THCA_data_analysis/project/.venv/lib/python3.12/site-packages/pandas/core/reshape/merge.py", line 840, in _reindex_and_concat
    llabels, rlabels = _items_overlap_with_suffix(
                       ^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/seungho/personal/THCA_data_analysis/project/.venv/lib/python3.12/site-packages/pandas/core/reshape/merge.py", line 2721, in _items_overlap_with_suffix
    raise ValueError(f"columns overlap but no suffix specified: {to_rename}")
ValueError: columns overlap but no suffix specified: Index(['patient_id'], dtype='object')

[2026-04-27 01:55:38] end v17p3_A1_scrna_heterogeneity.py code=1 sec=5.9
[2026-04-27 01:55:40] v17p3_A3_drug_response.py stdout [2026-04-27 01:55:34] v17p3 A3 start

[2026-04-27 01:55:40] v17p3_A3_drug_response.py stderr ow="left").dropna(subset=["symbol"]).drop_duplicates("symbol").set_index("symbol")
         ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/seungho/personal/THCA_data_analysis/project/.venv/lib/python3.12/site-packages/pandas/core/frame.py", line 10859, in merge
    return merge(
           ^^^^^^
  File "/home/seungho/personal/THCA_data_analysis/project/.venv/lib/python3.12/site-packages/pandas/core/reshape/merge.py", line 170, in merge
    op = _MergeOperation(
         ^^^^^^^^^^^^^^^^
  File "/home/seungho/personal/THCA_data_analysis/project/.venv/lib/python3.12/site-packages/pandas/core/reshape/merge.py", line 794, in __init__
    ) = self._get_merge_keys()
        ^^^^^^^^^^^^^^^^^^^^^^
  File "/home/seungho/personal/THCA_data_analysis/project/.venv/lib/python3.12/site-packages/pandas/core/reshape/merge.py", line 1298, in _get_merge_keys
    right_keys.append(right._get_label_or_level_values(rk))
                      ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/seungho/personal/THCA_data_analysis/project/.venv/lib/python3.12/site-packages/pandas/core/generic.py", line 1914, in _get_label_or_level_values
    raise KeyError(key)
KeyError: 'entrezGeneId'

[2026-04-27 01:55:40] end v17p3_A3_drug_response.py code=1 sec=7.5
[2026-04-27 01:55:41] v17p3 A2 done {'driver_dm_correlation_braf_vs_probdm2': -0.5779189851110272, 'braf_dm2_like_n': 2, 'ras_dm1_like_n': 15}
[2026-04-27 01:55:41] v17p3_A2_dm_in_braf_ras.py stdout [2026-04-27 01:55:34] v17p3 A2 start
[2026-04-27 01:55:41] v17p3 A2 done {'driver_dm_correlation_braf_vs_probdm2': -0.5779189851110272, 'braf_dm2_like_n': 2, 'ras_dm1_like_n': 15}

[2026-04-27 01:55:41] end v17p3_A2_dm_in_braf_ras.py code=0 sec=8.6
[2026-04-27 01:55:41] start v17p3_A4_pancancer.py
[2026-04-27 01:55:41] start v17p3_A5_genomic_immune.py
[2026-04-27 01:55:41] start v17p3_A6_tf_circuit.py
[2026-04-27 01:55:42] v17p3 A5 start
[2026-04-27 01:55:43] v17p3 A4 start
[2026-04-27 01:55:43] v17p3 A6 start
[2026-04-27 01:55:43] v17p3 A4 done {'applicable_cancers': 5, 'max_cluster_delta': 3.3378477096557617}
[2026-04-27 01:55:44] v17p3_A4_pancancer.py stdout [2026-04-27 01:55:43] v17p3 A4 start
[2026-04-27 01:55:43] v17p3 A4 done {'applicable_cancers': 5, 'max_cluster_delta': 3.3378477096557617}

[2026-04-27 01:55:44] end v17p3_A4_pancancer.py code=0 sec=2.5
[2026-04-27 01:55:46] v17p3_A5_genomic_immune.py stdout [2026-04-27 01:55:42] v17p3 A5 start

[2026-04-27 01:55:46] v17p3_A5_genomic_immune.py stderr roject/.venv/lib/python3.12/site-packages/pandas/core/apply.py", line 190, in agg
    return self.agg_dict_like()
           ^^^^^^^^^^^^^^^^^^^^
  File "/home/seungho/personal/THCA_data_analysis/project/.venv/lib/python3.12/site-packages/pandas/core/apply.py", line 423, in agg_dict_like
    return self.agg_or_apply_dict_like(op_name="agg")
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/seungho/personal/THCA_data_analysis/project/.venv/lib/python3.12/site-packages/pandas/core/apply.py", line 1603, in agg_or_apply_dict_like
    result_index, result_data = self.compute_dict_like(
                                ^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/seungho/personal/THCA_data_analysis/project/.venv/lib/python3.12/site-packages/pandas/core/apply.py", line 462, in compute_dict_like
    func = self.normalize_dictlike_arg(op_name, selected_obj, func)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/seungho/personal/THCA_data_analysis/project/.venv/lib/python3.12/site-packages/pandas/core/apply.py", line 663, in normalize_dictlike_arg
    raise KeyError(f"Column(s) {list(cols)} do not exist")
KeyError: "Column(s) ['gene'] do not exist"

[2026-04-27 01:55:46] end v17p3_A5_genomic_immune.py code=1 sec=5.1
[2026-04-27 01:55:46] v17p3 A6 done {'pdtc_auc_using_low_rai': 0.011764705882352951}
[2026-04-27 01:55:47] v17p3_A6_tf_circuit.py stdout [2026-04-27 01:55:43] v17p3 A6 start
[2026-04-27 01:55:46] v17p3 A6 done {'pdtc_auc_using_low_rai': 0.011764705882352951}

[2026-04-27 01:55:47] end v17p3_A6_tf_circuit.py code=0 sec=5.6
[2026-04-27 01:55:47] start v17p3_build_report.py
[2026-04-27 01:55:48] end v17p3_build_report.py code=0 sec=1.3
[2026-04-27 01:55:48] start v17p3_make_docs.py
[2026-04-27 01:55:49] end v17p3_make_docs.py code=0 sec=1.3
[2026-04-27 01:55:55] v17p3 F2 start
[2026-04-27 01:59:12] v17p3 A6 start
[2026-04-27 01:59:13] v17p3 A3 start
[2026-04-27 01:59:13] v17p3 A5 start
[2026-04-27 01:59:13] v17p3 A4 start
[2026-04-27 01:59:13] v17p3 A2 start
[2026-04-27 01:59:13] v17p3 A1 start
[2026-04-27 01:59:14] v17p3 A4 done {'applicable_cancers': 5, 'max_cluster_delta': 3.3378477096557617}
[2026-04-27 01:59:16] v17p3 A6 done {'pdtc_auc_using_low_rai': 0.011764705882352951}
[2026-04-27 01:59:20] v17p3 A2 done {'driver_dm_correlation_braf_vs_probdm2': -0.5779189851110272, 'braf_dm2_like_n': 2, 'ras_dm1_like_n': 15}
[2026-04-27 02:00:10] v17p3 A3 start
[2026-04-27 02:00:10] v17p3 A5 start
[2026-04-27 02:00:10] v17p3 A1 start
[2026-04-27 02:00:13] v17p3 A5 done {'tmb_dm2_minus_dm1': 0.09215167548500883}
[2026-04-27 02:00:14] v17p3 A1 done {'patients': 7, 'mixed_patients': 6, 'dominance_ratio_range': [0.7566824644549763, 0.9956058813587967]}
[2026-04-27 02:00:30] v17p3 A3 start
[2026-04-27 02:02:51] v17p3 A3 start
[2026-04-27 02:02:51] v17p3 A3 done {'n_drugs_fdr_lt_0_1': 0, 'best_dm1_drug': None}
```

### task exit/raw lines
| timestamp           | event   | task   | code   | sec   |
|:--------------------|:--------|:-------|:-------|:------|
| 2026-04-27 01:52:27 | start   | v      |        |       |
| 2026-04-27 01:52:27 | start   | v      |        |       |
| 2026-04-27 01:52:27 | start   | v      |        |       |
| 2026-04-27 01:52:32 | end     | v      |        |       |
| 2026-04-27 01:54:05 | end     | v      |        |       |
| 2026-04-27 01:55:32 | end     | v      |        |       |
| 2026-04-27 01:55:32 | start   | v      |        |       |
| 2026-04-27 01:55:32 | start   | v      |        |       |
| 2026-04-27 01:55:32 | start   | v      |        |       |
| 2026-04-27 01:55:38 | end     | v      |        |       |
| 2026-04-27 01:55:40 | end     | v      |        |       |
| 2026-04-27 01:55:41 | end     | v      |        |       |
| 2026-04-27 01:55:41 | start   | v      |        |       |
| 2026-04-27 01:55:41 | start   | v      |        |       |
| 2026-04-27 01:55:41 | start   | v      |        |       |
| 2026-04-27 01:55:44 | end     | v      |        |       |
| 2026-04-27 01:55:46 | end     | v      |        |       |
| 2026-04-27 01:55:47 | end     | v      |        |       |
| 2026-04-27 01:55:47 | start   | v      |        |       |
| 2026-04-27 01:55:48 | end     | v      |        |       |
| 2026-04-27 01:55:48 | start   | v      |        |       |
| 2026-04-27 01:55:49 | end     | v      |        |       |

### Python env
| package      | version   |
|:-------------|:----------|
| python       | 3.12.3    |
| pandas       | 2.3.3     |
| numpy        | 2.4.4     |
| scikit-learn | 1.8.0     |
| scanpy       | 1.12.1    |
| plotly       | 6.7.0     |
| gseapy       | 1.1.13    |
| mygene       | MISSING   |
| lifelines    | 0.30.3    |
| hdbscan      | MISSING   |
| umap-learn   | 0.5.12    |
| scipy        | 1.17.1    |
| requests     | 2.33.1    |
| pyscenic     | MISSING   |

### 가장 큰 raw 데이터 fetch 성공/실패 list
| source                                                                               | script                        | status          | evidence                                                                                            |
|:-------------------------------------------------------------------------------------|:------------------------------|:----------------|:----------------------------------------------------------------------------------------------------|
| https://mygene.info/v3/query                                                         | v17p3_F1_external_recovery.py | observed_output | /home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/F1_gene_recovery_mapping.tsv |
| https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2021_Human | v17p3_F2_msigdb_proper.py     | observed_output | /home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/F2_gsea_kegg_proper.tsv      |
| local MSigDB Hallmark cache                                                          | v17p3_F2_msigdb_proper.py     | observed_input  | /data/thca/data_raw/msigdb_cache/MSigDB_Hallmark_2020.json                                          |
| local Reactome cache                                                                 | v17p3_F2_msigdb_proper.py     | observed_input  | /data/thca/data_raw/msigdb_cache/Reactome_2022.txt                                                  |
| /data/thca/data_processed/bulk_rnaseq_v3/GSE213647_v3_log2.tsv                       | v17p3_F1_external_recovery.py | observed_input  | /data/thca/data_processed/bulk_rnaseq_v3/GSE213647_v3_log2.tsv                                      |

## SEC 1. F1 GSE213647 recovery — RAW

### 1.1 진단 결과
- file: `/home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/F1_gene_id_diagnosis.tsv`
| dataset   | first_column   |   ensembl_like_fraction_top500 | failure_mode                                        |   raw_rows |   raw_samples |
|:----------|:---------------|-------------------------------:|:----------------------------------------------------|-----------:|--------------:|
| GSE213647 | gene_symbol    |                              1 | column named gene_symbol but values are Ensembl IDs |      51389 |           632 |

MISSING: TCGA expression head read error: [Errno 2] No such file or directory: '/data/thca/data_processed/bulk_rnaseq_v3/TCGA_THCA_tumor_log2tpm.tsv'
- GSE213647 expression file: `/data/thca/data_processed/bulk_rnaseq_v3/GSE213647_v3_log2.tsv`
- GSE213647 first column: `gene_symbol`
| gse213647_identifier_head10   |
|:------------------------------|
| ENSG00000223972               |
| ENSG00000227232               |
| ENSG00000278267               |
| ENSG00000243485               |
| ENSG00000237613               |
| ENSG00000240361               |
| ENSG00000186092               |
| ENSG00000238009               |
| ENSG00000233750               |
| ENSG00000268903               |

- identifier convention raw strings:
| dataset | convention |
|---|---|
| TCGA expression | `MISSING` |
| GSE213647 expression | `gene_symbol` |

### 1.2 매핑 결과
- file: `/home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/F1_gene_recovery_mapping.tsv`
| ensembl         | symbol   |
|:----------------|:---------|
| ENSG00000000003 | TSPAN6   |
| ENSG00000000005 | TNMD     |
| ENSG00000000419 | DPM1     |
| ENSG00000000457 | SCYL3    |
| ENSG00000000460 | FIRRM    |
| ENSG00000000938 | FGR      |
| ENSG00000000971 | CFH      |
| ENSG00000001036 | FUCA2    |
| ENSG00000001084 | GCLC     |
| ENSG00000001167 | NFYA     |
| ENSG00000001460 | STPG1    |
| ENSG00000001461 | NIPAL3   |
| ENSG00000001497 | LAS1L    |
| ENSG00000001561 | ENPP4    |
| ENSG00000001617 | SEMA3F   |
| ENSG00000001626 | CFTR     |
| ENSG00000001629 | ANKIB1   |
| ENSG00000001630 | CYP51A1  |
| ENSG00000001631 | KRIT1    |
| ENSG00000002016 | RAD52    |
| ENSG00000002330 | BAD      |
| ENSG00000002549 | LAP3     |
| ENSG00000002586 | CD99     |
| ENSG00000002587 | HS3ST1   |
| ENSG00000002726 | AOC1     |
| ENSG00000002745 | WNT16    |
| ENSG00000002746 | HECW1    |
| ENSG00000002822 | MAD1L1   |
| ENSG00000002834 | LASP1    |
| ENSG00000002919 | SNX11    |
| ENSG00000002933 | TMEM176A |
| ENSG00000003056 | M6PR     |
| ENSG00000003096 | KLHL13   |
| ENSG00000003137 | CYP26B1  |
| ENSG00000003147 | ICA1     |
| ENSG00000003249 | DBNDD1   |
| ENSG00000003393 | ALS2     |
| ENSG00000003400 | CASP10   |
| ENSG00000003402 | CFLAR    |
| ENSG00000003436 | TFPI     |
| ENSG00000003509 | NDUFAF7  |
| ENSG00000003756 | RBM5     |
| ENSG00000003987 | MTMR7    |
| ENSG00000003989 | SLC7A2   |
| ENSG00000004059 | ARF5     |
| ENSG00000004139 | SARM1    |
| ENSG00000004142 | POLDIP2  |
| ENSG00000004399 | PLXND1   |
| ENSG00000004455 | AK2      |
| ENSG00000004468 | CD38     |

...

| ensembl         | symbol       |
|:----------------|:-------------|
| ENSG00000284505 | LYNX1-SLURP2 |
| ENSG00000284508 | MIR133A2     |
| ENSG00000284516 | VAMP9P       |
| ENSG00000284523 | LOC105375421 |
| ENSG00000284526 | LOC122455342 |
| ENSG00000284543 | LINC01226    |
| ENSG00000284546 | SSU72L3      |
| ENSG00000284575 | MIR4793      |
| ENSG00000284584 | MIR5006      |
| ENSG00000284594 | MIR7847      |
| ENSG00000284601 | LOC105378736 |
| ENSG00000284605 | ENPP7P15     |
| ENSG00000284607 | LOC105374988 |
| ENSG00000284609 | OR8B3        |
| ENSG00000284623 | LINC02786    |
| ENSG00000284624 | LOC105374312 |
| ENSG00000284629 | SCYGR1       |
| ENSG00000284631 | SCYGR4       |
| ENSG00000284634 | LOC374443    |
| ENSG00000284635 | SCYGR8       |
| ENSG00000284638 | SMIM44       |
| ENSG00000284640 | LOC105378644 |
| ENSG00000284641 | LOC124903871 |
| ENSG00000284643 | SCYGR2       |
| ENSG00000284645 | LOC105378710 |
| ENSG00000284650 | LOC124904027 |
| ENSG00000284668 | LINC02780    |
| ENSG00000284670 | RPSAP56      |
| ENSG00000284674 | LINC02781    |
| ENSG00000284677 | ZFP69B-DT    |
| ENSG00000284678 | LINC02811    |
| ENSG00000284680 | OR8B2        |
| ENSG00000284688 | BTN1A1P1     |
| ENSG00000284689 | OR7E84P      |
| ENSG00000284690 | CD300H       |
| ENSG00000284691 | LOC728743    |
| ENSG00000284693 | LINC02606    |
| ENSG00000284696 | LINC02808    |
| ENSG00000284701 | TMEM247      |
| ENSG00000284704 | SCYGR3       |
| ENSG00000284706 | LOC100420006 |
| ENSG00000284713 | SMIM38       |
| ENSG00000284718 | SCYGR7       |
| ENSG00000284719 | LOC130932201 |
| ENSG00000284720 | LOC124904023 |
| ENSG00000284723 | OR8S1        |
| ENSG00000284730 | TMDD1        |
| ENSG00000284731 | LOC124906284 |
| ENSG00000284736 | RNA5SP343    |
| ENSG00000284738 | LOC101928059 |

- mapped TierA67_clean count: `66` / `76`
- mapped TierA67_clean gene list:
| mapped_tiera67_gene   |
|:----------------------|
| AKT1                  |
| ALK                   |
| APC                   |
| ATM                   |
| BRAF                  |
| CD274                 |
| CD8A                  |
| CDH1                  |
| CDH2                  |
| CDKN2A                |
| CDKN2B                |
| CTNNB1                |
| DIO1                  |
| DIO2                  |
| DUOX1                 |
| DUOX2                 |
| DUSP4                 |
| DUSP5                 |
| DUSP6                 |
| EIF1AX                |
| ETV4                  |
| ETV5                  |
| FOSL1                 |
| FOXE1                 |
| FOXP3                 |
| GLIS3                 |
| HLA-DRA               |
| HRAS                  |
| IDO1                  |
| IYD                   |
| KLK10                 |
| KRAS                  |
| LOX                   |
| MET                   |
| MMP9                  |
| MSH2                  |
| NKX2-1                |
| NRAS                  |
| NTRK1                 |
| NTRK3                 |
| PAX8                  |
| PHLDA1                |
| PIK3CA                |
| PPARG                 |
| PTEN                  |
| RET                   |
| SLC26A4               |
| SLC5A5                |
| SLC5A8                |
| SNAI1                 |
| SNAI2                 |
| SPRY1                 |
| SPRY2                 |
| SPRY4                 |
| TERT                  |
| TG                    |
| THADA                 |
| THRA                  |
| THRB                  |
| TP53                  |
| TPO                   |
| TSHR                  |
| TWIST1                |
| VIM                   |
| ZEB1                  |
| ZEB2                  |

- mapping failed TierA67_clean gene list:
| missing_tiera67_gene                             |
|:-------------------------------------------------|
| # HARDCODED COMPOSITE TIERA67 ENTRIES            |
| # NOTE: PAX8 APPEARS IN TWO CATEGORIES BY DESIGN |
| # [AGGRESSIVE_MARKER]                            |
| # [DEDIFF_INVASION]                              |
| # [DRIVER_ANCHOR]                                |
| # [IMMUNE_STROMAL_LIGHT]                         |
| # [MAPK_OUTPUT_ERK]                              |
| # [TDS_CORE]                                     |
| # [THYROID_LINEAGE_EXTRA]                        |

### 1.3 5-cohort transfer
- file: `/home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/F1_external_5cohort_recovery.tsv`
| cohort    | classifier   |   n |   gene_overlap |      AUC |   DIA_AUC |
|:----------|:-------------|----:|---------------:|---------:|----------:|
| GSE27155  | LogReg       |  54 |             50 | 0.698113 |  0.698113 |
| GSE27155  | RandomForest |  54 |             50 | 0.339623 |  0.660377 |
| GSE27155  | SVM_RBF      |  54 |             50 | 0.150943 |  0.849057 |
| GSE27155  | kNN          |  54 |             50 | 0.556604 |  0.556604 |
| GSE76039  | LogReg       |  37 |             54 | 1        |  1        |
| GSE76039  | RandomForest |  37 |             54 | 0.972222 |  0.972222 |
| GSE76039  | SVM_RBF      |  37 |             54 | 0.277778 |  0.722222 |
| GSE76039  | kNN          |  37 |             54 | 0.916667 |  0.916667 |
| GSE213647 | LogReg       | 381 |             55 | 0.985185 |  0.985185 |
| GSE213647 | RandomForest | 381 |             55 | 0.981699 |  0.981699 |
| GSE213647 | SVM_RBF      | 381 |             55 | 0.963181 |  0.963181 |
| GSE213647 | kNN          | 381 |             55 | 0.975447 |  0.975447 |

- expected columns requested by prompt: `cohort, n_dark_candidates, n_gene_overlap, classifier, AUC, DIA-AUC, identifiability_no_combat, identifiability_combat`
- MISSING columns in file: `['n_dark_candidates', 'n_gene_overlap', 'identifiability_no_combat', 'identifiability_combat']`

### 1.4 GSE126698, GSE97466 forced classification
- cohort: `GSE126698`
MISSING: no rows in `/home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/F1_external_predictions.tsv` for `GSE126698`
- cohort: `GSE97466`
MISSING: no rows in `/home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/F1_external_predictions.tsv` for `GSE97466`
## SEC 2. F2 MSigDB GSEA — RAW

### 2.1 패키지/library status
| item             | value   |
|:-----------------|:--------|
| gseapy version   | 1.1.13  |
| Hallmark cache   | present |
| Reactome cache   | present |
| KEGG output file | present |
| library         |   pathway_count |
|:----------------|----------------:|
| Hallmark        |              50 |
| Reactome        |            1818 |
| KEGG_2021_Human |             319 |
### 2.2 Hallmark 결과
- file: `/home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/F2_gsea_hallmark_proper.tsv`
| pathway.1                         |       NES |   NOM p-val |   FDR q-val | Lead_genes                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
|:----------------------------------|----------:|------------:|------------:|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Oxidative Phosphorylation         |  1.90388  |  0          | 0           | ALDH6A1;VDAC3;OGDH;ACADSB;NNT;NDUFA5;AFG3L2;COX10;ATP5F1B;MPC1;ACADM;SUCLG1;CYCS;ETFDH;MRPS11;PHYH;GOT2;DLD;ISCA1;SLC25A20;NDUFB2;SLC25A4;ACAT1;ATP5F1A;FDX1;AIFM1;HSPA9;UQCRC2;ATP6V1H;ATP5MC1;NDUFC2;FH;COX8A;MDH2;CYC1;FXN;COX15;ECI1;SUCLA2;IDH3B;SDHD;ACO2;ATP5PD;POR;SDHC;MRPL15;TOMM70;UQCRFS1;HSD17B10;NDUFS1;RETSAT;ATP5PB;PDHA1;HADHB;SLC25A5;NDUFB5;ATP5MC3;MRPL34;HCCS;DLAT;OAT;DLST;ATP6V1D;COX7B;UQCRH;PRDX3;NDUFA8;ECHS1;CS;NDUFA1;NDUFS2;IMMT;UQCR10;COX11;ATP5F1C;NDUFV1;OPA1;NDUFA6;SLC25A3;HADHA;PDHX;TIMM17A;NDUFB3;MFN2;ATP5PF;GRPEL1;MTX2;CPT1A;COX5A;ETFA;RHOT1;SLC25A11;COX6C;NDUFB6;ATP6AP1;MRPS30;ATP6V0E1;PDHB;NDUFS4;GPI;COX6A1;SDHB;COX7C;VDAC2;ALAS1;ATP6V0B;SDHA;COX7A2;IDH3A;LRPPRC;NDUFA2;ETFB;OXA1L;NDUFA9;MRPS22;UQCRQ;DECR1;NDUFS3;CYB5A;UQCRC1;SURF1                                                              |
| Allograft Rejection               | -1.56194  |  0          | 0           | TIMP1;ICAM1;GBP2;HLA-DQA1;STAT4;CCL13;ITGB2;MAP4K1;LTB;HLA-DOA;IGSF6;INHBB;CCL22;SPI1;HLA-DRA;NCF4;CD2;CD74;IL2RG;CCL5;CD3E;CD86;FGR;CD40;FCGR2B;CCR1;WAS;HLA-DOB;CCR2;HLA-G;STAT1;CCR5;LCK;CD3D;PTPRC;CD7;IL2RA;CD4;CTSS;IL2RB;HLA-DMA;LY86;ST8SIA4;C2;FLNA;ZAP70;CXCL9;FAS;IL18;MMP9;CCL19;ITK;GZMA;HCLS1;NLRP3;IL4R;HLA-A;FYB1;BCAT1;LCP2;TLR2;CXCL13;CD47;CD3G;CCL2;IL12RB1;CD247;LYN;TGFB1;IRF8;LIF;CDKN2A;SIT1;ITGAL;CD28;INHBA;CSF1;BCL3;CD79A;TAP1;CRTAM;CCL4;IL7;WARS1;GZMB;CD40LG;TAP2;CD8B;GPR65;TLR6;EREG;CAPG;HLA-E;B2M;IL16;CD8A;CSK;PSMB10;TRAT1;CXCR3;RIPK2;CD96;IFNAR2;PRF1;DEGS1;SRGN;FASLG;IFNGR1;IL6;KLRD1;PTPN6;IL27RA;IRF7;IL1B;ELF4;TLR1;CCND3;SOCS1;CCND2;HLA-DMB;RPS9;IL18RAP;LY75                                                                                                                                            |
| Interferon Gamma Response         | -1.55656  |  0          | 0           | ICAM1;HLA-DQA1;STAT4;NOD1;PELI1;MVP;HLA-DRB1;LY6E;CFH;UPP1;CD74;CCL5;IRF5;CD86;PTGS2;CD40;CIITA;METTL7B;HLA-G;STAT1;PSMB9;HLA-B;IL10RA;CSF2RB;PARP14;CD69;IL2RB;HLA-DMA;PSMB8;ST8SIA4;C1S;CXCL9;FAS;CFB;TRIM14;SERPING1;EPSTI1;SECTM1;SP110;C1R;GZMA;SLAMF7;ITGB7;SAMHD1;IL4R;CXCL10;HLA-A;LCP2;CCL2;MARCHF1;FPR1;IRF8;PNP;CASP1;CASP4;ZBP1;IFIT2;NLRC5;FGL2;IRF1;OASL;TAP1;PLA2G4A;SELP;IFIT3;ST3GAL5;CXCL11;DDX60;SAMD9L;IL7;APOL6;WARS1;CMKLR1;AUTS2;UBE2L6;MX2;HERC6;PARP12;B2M;RSAD2;FCGR1A;ISG15;CASP8;TNFAIP2;IL15RA;PSMB10;RIPK2;XCL1;VCAM1;PLSCR1;IFNAR2;BST2;TRIM21;IRF9;TNFAIP6;ISG20;VAMP8;MTHFD2;MYD88;ARID5B;IL6;PTPN6;IRF7;SOCS3;GPR18;TNFAIP3;OAS3;RTP4;CD274;PSME2;BATF2;SOCS1;NMI;IL18BP;PML;PSME1;VAMP5;IFITM2;LYSMD2;IFI44;NAMPT;BTG1;IFITM3;DDX58;GBP6;IFIH1;PTPN2;KLRK1;IFIT1;TAPBP;CASP7;GCH1;CASP3;STAT2;XAF1;GBP4;HIF1A;IFI35 |
| TNF-alpha Signaling via NF-kB     | -1.54312  |  0          | 0           | LAMB3;PLAU;PTPRE;PLAUR;DUSP5;ICAM1;TNC;SDC4;G0S2;CXCL2;LDLR;DUSP4;BIRC3;DRAM1;TSC22D1;CCL5;IL7R;PTGS2;OLR1;TIPARP;NFKBIE;PTGER4;PLEK;BCL2A1;SERPINB8;CSF2;CD69;SLC2A3;IL18;NFKB2;MAFF;PHLDA2;IER3;DUSP2;BHLHE40;CXCL10;TLR2;CCL2;RELB;TNFSF9;ZC3H12A;CCL20;B4GALT5;FOSL1;LIF;TNFRSF9;INHBA;CSF1;CD83;PFKFB3;BCL3;IFIT2;IRF1;TAP1;CLCF1;AREG;MYC;STAT5A;PMEPA1;CXCL11;CCL4;KYNU;MSC;TRAF1;SLC2A6;GPR183;TNFAIP2;IL15RA;CEBPD;RIPK2;F3;CXCL3;TNFAIP6;HBEGF;SERPINE1;IER5;RNF19B;KLF4;IL6;MCL1;KLF6;GADD45A;CFLAR;SOCS3;IL1B;SPHK1;TNFAIP3;TGIF1;CD44;CXCL1;GEM;PPP1R15A;IL23A;SAT1;TRIP10;CCND1;REL;TNF;CCRL2;NAMPT;CD80;EGR3;BTG1;TANK;DDX58;ATF3;BTG3;IFIH1;JUNB;ZFP36;FOSL2;IFNGR2;RELA;GCH1                                                                                                                                                          |
| Interferon Alpha Response         | -1.4886   |  0          | 0.000241774 | GBP2;SELL;LY6E;CD74;PSMB9;PARP14;PSMB8;C1S;TRIM14;EPSTI1;SP110;IL4R;CXCL10;HLA-C;CD47;LAMP3;CASP1;IFITM1;CSF1;PARP9;IFIT2;IRF1;OASL;TAP1;IFIT3;CXCL11;DDX60;SAMD9L;IL7;WARS1;UBE2L6;HERC6;PARP12;B2M;SAMD9;RSAD2;ISG15;CASP8;PROCR;NCOA7;RIPK2;PLSCR1;BST2;TRIM21;IRF9;ISG20;OAS1;IRF7;TENT5A;RTP4;PSME2;BATF2;NMI;PSME1;IFITM2;MOV10;UBA7;IFI44;CCRL2;IFITM3;IFIH1;STAT2;GBP4;IFI35;LPAR6                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| Coagulation                       | -1.47606  |  0          | 0.000276313 | SERPINA1;TIMP1;PLAU;FN1;DPP4;PROS1;CTSH;APOC1;C3;ANXA1;PRSS23;DUSP6;CLU;TMPRSS6;CFH;COMP;CFI;CTSE;OLR1;MMP7;PLEK;C2;C1S;CTSK;CFB;MAFF;SERPING1;MMP9;C1R;C1QA;MMP11;HPN;BMP1;GSN;ITGA2;THBS1;MMP2;PROC;F3;FYN;SERPINE1;ANG                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| Epithelial Mesenchymal Transition | -1.47825  |  0          | 0.000322365 | PDLIM4;TIMP1;FN1;COL8A2;PLAUR;TNC;SDC4;EMP3;LGALS1;CDH11;TNFRSF12A;COL7A1;PRSS2;SLIT2;COL1A1;NNMT;LUM;COMP;FSTL3;DPYSL3;CDH6;ANPEP;CXCL8;ECM1;QSOX1;TGM2;VCAN;CRLF1;FLNA;DCN;LOXL1;FAS;DST;NT5E;SFRP4;BASP1;FAP;WIPF1;CTHRC1;TGFB1;BMP1;SPP1;INHBA;LAMA3;RGS4;SPOCK1;AREG;ITGA2;PCOLCE;IL32;PMEPA1;ADAM12;FBLN2;COL11A1;COL1A2;FBN2;CAPG;COL3A1;LOX;THBS2;NTM;FBLN1;THBS1;FGF2;SERPINE2;MMP2;VCAM1;MFAP5;TGFBI;COL6A3;SERPINE1;IL6;COL5A1;GADD45A;TPM4;TNFAIP3;GLIPR1;CD44;CXCL1;GEM;COLGALT1;COL16A1;LRP1;MMP14;CD59;LRRC15;COL6A2;WNT5A;SAT1;LAMA1;PRRX1;FOXC2;MYL9;MGP;SCG2;NOTCH2                                                                                                                                                                                                                                                                  |
| Complement                        | -1.48094  |  0          | 0.000386838 | SERPINA1;TIMP1;FN1;DPP4;LGALS3;PLAUR;DUSP5;SPOCK2;CTSH;APOC1;CTSC;C3;DUSP6;CLU;TMPRSS6;CFH;FCER1G;S100A9;CD55;ITGAM;CCL5;OLR1;WAS;PSMB9;LCK;PLEK;HPCAL4;PHEX;PIK3R5;CTSS;GNB4;C2;F5;C1S;CFB;MAFF;SERPING1;C1R;GZMA;C1QC;C1QA;RASGRP1;LCP2;LYN;APOBEC3G;CASP1;CASP4;IRF1;PLA2G4A;MMP12;RHOG;KYNU;CR1;GZMB;CD40LG;TIMP2;SRC;FCN1;GZMK;CP;PIK3CG;GNAI2;CR2;PLSCR1;DOCK9;F3;CDK5R1;FYN;SERPINE1;ANG;IL6;DOCK10;IRF7;ANXA5;TNFAIP3;CASP10;PFN1;CXCL1;LRP1;S100A13;APOBEC3F;MMP14;GNGT2;DGKH;CD59;MMP13;LIPA;CTSD                                                                                                                                                                                                                                                                                                                                            |
| Inflammatory Response             | -1.47076  |  0          | 0.000483548 | TIMP1;PTPRE;PLAUR;ICAM1;KCNJ2;EMP3;MARCO;HAS2;SELL;LY6E;SCN1B;LDLR;EBI3;CCL22;AHR;MET;CD55;CCL5;IL7R;GNA15;RGS1;CXCL8;ITGB8;CD40;OLR1;TNFRSF1B;PTGER4;LCK;CLEC5A;CYBB;IL10RA;PIK3R5;CD69;IL2RB;CCL17;OSMR;NOD2;CXCL9;IL18;SLAMF1;PTAFR;IRAK2;HPN;NLRP3;CHST2;IL4R;CXCR6;CXCL10;RASGRP1;LCP2;TLR2;C3AR1;MSR1;HRH1;CCL2;CD48;TNFSF9;CD70;CCR7;LYN;LAMP3;CCL20;FPR1;LIF;AQP9;TNFRSF9;PDPN;IFITM1;INHBA;CSF1;IRF1;MYC;SGMS2;RHOG;AXL;CXCL11;C5AR1;CMKLR1;CD14;P2RX7;OSM;EREG;CX3CL1;ADM;LTA;GPR183;IL15RA;RIPK2;TNFSF15;PTGER2;BST2;F3;TNFAIP6;HBEGF;SERPINE1;IL6;CSF3R;KLF6;IRF7;IL1B;SPHK1;TLR1;RTP4;ROS1;IL18RAP;MMP14;NMI;PTGIR;CCL24                                                                                                                                                                                                                  |
| KRAS Signaling Up                 | -1.4538   |  0          | 0.00139692  | KCNN4;PLAU;SLPI;PLAUR;G0S2;DUSP6;TMEM100;ITGB2;MAP4K1;CFH;SCN1B;FCER1G;BIRC3;LAT2;LCP1;ADAM8;LAPTM5;IL2RG;CD37;CSF2RA;IL7R;PTGS2;MPZL2;TMEM176A;TMEM176B;TNFRSF1B;IL10RA;CTSS;CSF2;PSMB8;IKZF1;CFB;MMP9;DOCK2;LY96;MMP11;APOD;CXCL10;GPNMB;C3AR1;HSD11B1;ETV4;SPON1;CCL20;PCSK1N;IRF8;SPP1;LIF;SCG5;INHBA;CXCR4;ITGA2;TLR8;CMKLR1;ADAMDEC1;TRAF1;PRDM1;ETV5;EREG;VWA5A;ALDH1A3;PDCD1LG2;TRIB2;SNAP25;ETV1;HBEGF;ANGPTL4;SPRY2;KLF4;IL1B;TNFAIP3;RETN;EMP1;CCND2;GLRX;PPP1R15A;MALL;PIGR;WNT7A;TMEM158;MAP3K1;ITGBL1;RBP4;PRRX1;ANKH;DCBLD2                                                                                                                                                                                                                                                                                                             |
| Apoptosis                         | -1.44761  |  0          | 0.00145064  | IGFBP6;TIMP1;LGALS3;ANXA1;CLU;TNFRSF12A;BIRC3;BID;LUM;CD2;CDC25B;GNA15;CCNA1;PLCB2;PMAIP1;CD69;DCN;FAS;IL18;ERBB3;IER3;LMNA;GSN;CASP1;DPYD;CASP4;SLC20A1;IRF1;TAP1;PEA15;GPX1;CD14;TIMP2;EREG;HMOX1;BCL2L1;CASP8;MMP2;PRF1;ISG20;FASLG;IFNGR1;PLPPR4;IL6;MCL1;GADD45A;CFLAR;IL1B;RARA;EMP1;TSPO;PPT1;CD44;CCND2;SAT1;CCND1;GUCY2D;ANKH;TNF;EGR3;BAX;IFITM3;ATF3;BTG3;GPX4;BGN;CASP7;TGFB2;RELA;GCH1;CASP3;IGF2R;HMGB2                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| Apical Junction                   | -1.44769  |  0.00110254 | 0.00158252  | LAMB3;ALOX15B;RAC2;ICAM1;ICAM5;CDH3;CDH11;ITGB4;AMIGO2;NECTIN4;MDK;SDC3;SLIT2;CDH6;CD86;MPZL2;SIRPA;DSC3;PTPRC;ITGA3;GRB7;VCAN;MMP9;RRAS;ACTN1;ZYX;MYH10;FYB1;CD276;VASP;BMP1;LAMA3;ITGA2;CDH4;CADM3;MSN;NLGN2;TRAF1;CDH15;SRC;CX3CL1;YWHAH;ARPC2;MAPK13;GNAI2;MMP2;VCAM1;CLDN4;NECTIN1;TGFBI;ADAM15;BAIAP2;CD274;CALB2;PFN1;COL16A1;CLDN11;CNN2;MYL12B;ADRA1B;SGCE;RASA1;CLDN9;MYL9;SKAP2;CAP1;CLDN7;ITGA9;ADAM23                                                                                                                                                                                                                                                                                                                                                                                                                                     |
| p53 Pathway                       | -1.44841  |  0          | 0.00174077  | PTPRE;ALOX15B;SFN;S100A4;ITGB4;S100A10;INHBB;KRT17;TM4SF1;DRAM1;UPP1;PLK3;TSC22D1;TP63;TGFA;DGKA;FAS;IER3;ZFP36L1;TNFSF9;DEF6;TGFB1;EPHA2;LIF;CDKN2A;RAP2B;CASP1;RRAD;TAP1;HMOX1;VWA5A;PHLDA3;PROCR;ZMAT3;PERP;DDB2;CDK5R1;CLCA2;VAMP8;HBEGF;IER5;RNF19B;TCHH;RALGDS;KLF4;H2AW;STOM;GADD45A;RPL36;SPHK1;BAIAP2;CCND3;SOCS1;MAPKAPK3;CCND2;CDKN2B;VDR;PPP1R15A;AK1;PLXNB2;SAT1;DDIT4;CTSD;TNNI1;TPD52L1;PRMT2;BTG1;BAX;ATF3;SERPINB5                                                                                                                                                                                                                                                                                                                                                                                                                    |
| IL-6/JAK/STAT3 Signaling          | -1.42369  |  0.0025     | 0.00407562  | TNFRSF12A;LTB;EBI3;IL2RG;CSF2RA;CCR1;TNFRSF1B;STAT1;IL2RA;TNFRSF21;PIK3R5;CSF2RB;CSF2;OSMR;CXCL9;FAS;TNFRSF1A;IL4R;CXCL10;TLR2;CXCL13;IL12RB1;TGFB1;CSF1;IL1R2;IRF1;ITGA4;CXCL11;IL7;CD14;CRLF2;LEPR;HMOX1;IL15RA;IRF9;CXCL3;IFNGR1;MYD88;IL6;CSF3R;SOCS3;IL1B;SOCS1;CD44;CXCL1;PLA2G2A;TNF;IL13RA1;IL18R1;PTPN2;IL17RA;IFNGR2;STAT2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| IL-2/STAT5 Signaling              | -1.42464  |  0          | 0.00424034  | MUC1;SELL;LTB;SLC1A5;AHR;TNFRSF18;CD86;BATF3;TIAM1;ECM1;TNFRSF1B;DHRS3;BATF;TGM2;IL10RA;IL2RA;TNFRSF21;CSF2;IL2RB;SLC2A3;MAFF;CST7;NT5E;CTLA4;RHOH;BHLHE40;IL4R;IL1RL1;CXCL10;TNFSF11;CCR4;ETV4;CD48;NFKBIZ;IRF8;PNP;SPP1;LIF;CDCP1;TNFRSF9;CSF1;CD83;IL1R2;FGL2;ICOS;SELP;ADAM19;MYC;ODC1;TRAF1;CAPN3;GPR65;CAPG;BCL2L1;PLEC;IKZF2;EOMES;SYNGR2;AHNAK;PLSCR1;PTGER2;CD79B;IFNGR1;PLIN2;KLF6;TLR7;SYT11;CCND3;EMP1;WLS;GLIPR2;SOCS1;CD44;CCND2;TNFRSF8;AGER;MAPKAPK2;CTSZ;IFITM3;ANXA4;IL18R1;TTC39B;GPX4;GUCY1B1;MYO1E;S100A1;CASP3;GBP4;IGF2R;IL10;IRF4                                                                                                                                                                                                                                                                                              |
| Angiogenesis                      | -1.37298  |  0.0292826  | 0.0243708   | TIMP1;S100A4;LUM;OLR1;TNFRSF21;VCAN;SPP1;COL3A1;LPL;CCND2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| Cholesterol Homeostasis           | -1.20017  |  0.185606   | 0.370056    | LGALS3;PLAUR;CLU;TNFRSF12A;LDLR;SCD;S100A11;CXCL16;ETHE1;NIBAN1;ALDOC;PDK3;PLSCR1;LPL;ANXA5;TP53INP1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| Notch Signaling                   | -1.2036   |  0.192982   | 0.379887    | DTX4;LFNG;DTX1;DTX2;WNT5A;CCND1;NOTCH2;MAML2;WNT2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| Estrogen Response Late            | -1.19022  |  0.137625   | 0.383185    | SERPINA1;CXCL14;GJB3;SFN;KRT19;KLK10;PRSS23;GALE;MDK;ALDH3B1;S100A9;RET;CCNA1;TIAM1;NAB2;BATF;KLK11;IGSF1;NBL1;DUSP2;AGR2;EMP2;SNX10;LSR;TPSAB1;PTGES;AREG;MICB;RAB31;PGR;TFAP2C;SCNN1A;NRIP1;MAPK13;KRT13;PERP;SLC7A5;ISG20;NMU;KLF4;PTPN6;CD44;RBBP8;CCN5;SLC27A2;KCNK5;CCND1;PLXNB1;SIAH2;TPD52L1;EGR3;TPBG;SULT2B1;BTG3;ZFP36;PLK4;RAPGEFL1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
| UV Response Up                    | -1.17993  |  0.173034   | 0.394677    | ICAM1;CXCL2;BID;RET;RAB27A;AQP3;HLA-F;PTPRD;RASGRP1;LYN;ARRB2;RRAD;IRF1;TAP1;NPTXR;TUBA4A;NTRK3;CDO1;H2AX;HMOX1;HYAL2;TCHH;IL6;CCND3;SLC6A12;PPT1;DLG4;FGF18;CDKN2B;GLS;MMP14;TYRO3;BTG1;ATF3;SHOX2;BTG3;JUNB;NXF1;SELENOW;GCH1;CASP3;HTR7                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| Hypoxia                           | -1.12193  |  0.266667   | 0.543323    | PLAUR;S100A4;SDC4;SDC3;ANXA2;TIPARP;TGM2;SLC2A3;DCN;MAFF;IER3;CHST2;BHLHE40;HS3ST1;PFKFB3;SLC6A6;CXCR4;PAM;CAVIN1;KDELR3;PKP1;ALDOC;SRPX;PDK3;HMOX1;CP;ADM;LOX;PLAC8;TGFB3;F3;ISG20;TGFBI;SERPINE1;PLIN2;ANGPTL4;IL6;COL5A1;KLF6;TNFAIP3;WSB1;GLRX;PPP1R15A;AK4;CCN5;MAP3K1;DDIT4;SIAH2;SLC2A5;ERRFI1;BTG1;TPBG;SULT2B1;ATF3;DTNA;PPFIA4;ZFP36;FOSL2;STBD1;BGN;P4HA2;IGFBP3;ETS1                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| G2-M Checkpoint                   | -1.10909  |  0.268935   | 0.561575    | CDC25B;HMGA1;TGFB1;BCL3;JPT1;MYC;ODC1;H2AX;RAD23B;DMD;MKI67;AURKB;MYBL2;SLC7A5;BUB1;UBE2C;MCM5;KIF23;NDC80;BRCA2;EFNA5;PML;CCND1;KNL1;MCM6;CCNB2;KIF4A;NOTCH2;CCNA2;CENPF;BIRC5;PLK4;CDKN2C;TMPO;POLQ;HIF1A;PLK1;CDC45;ESPL1;LMNB1;TOP2A;ORC6;EZH2;CENPE;STMN1;TPX2;CKS2;SNRPD1;RASAL2;CKS1B;E2F2;CENPA;TROAP;KIF11;CDC7;KPNA2;NEK2;KIF15;BARD1;DR1;SFPQ;CDKN3;HMMR;FOXN3                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| Mitotic Spindle                   | -1.1222   |  0.277162   | 0.569426    | TIAM1;ARHGAP4;FLNA;DST;CDC42EP1;DOCK2;MYH10;SYNPO;WASF1;GSN;TUBA4A;PLEKHG2;CDK5RAP2;MYO9B;ARHGAP10;DLGAP5;ARF6;ABR;MAP1S;BUB1;CLIP2;KIF23;NDC80;CTTN;BRCA2;RASA2;ARHGEF2;RASA1;LLGL1;SEPTIN9;PIF1;KNTC1;NET1;CCNB2;KIF4A;NOTCH2;CNTRL;CENPF;HOOK3;BIRC5;MYO1E;ARHGAP27;MAP3K11;SHROOM1;ANLN;PLK1;ESPL1;LMNB1;TOP2A;NEDD9;CENPE;TPX2;AKAP13;RASAL2;PALLD;SPTAN1;RAPGEF5;VCL;ABI1;ACTN4;KIF11                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
| Hedgehog Signaling                | -1.08465  |  0.386724   | 0.591983    | NRP2;SLIT1;AMOT;HEY2;NRCAM;CDK5R1;PML;UNC5C;RASA1;GLI1;SCG2;CRMP1;L1CAM;CELSR1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| Estrogen Response Early           | -1.08798  |  0.326689   | 0.606957    | SYT12;SFN;KRT19;MUC1;KRT15;KLK10;PRSS23;ALDH3B1;INHBB;RET;TIAM1;TIPARP;DHRS3;AQP3;TGM2;PMAIP1;UGCG;ELF3;NBL1;ENDOD1;BHLHE40;RASGRP1;LAD1;PTGES;FCMR;AREG;MICB;RAB31;MYC;MREG;PGR;TFAP2C;SCNN1A;NRIP1;KRT13;SLC7A5;GREB1;OLFML3;KLF4;RARA;KCNK15;CALB2;CD44;NAV2;RBBP8;CCN5;SLC27A2;KCNK5;CCND1;SIAH2;DLC1;SVIL;TPD52L1;EGR3;CLDN7;TPBG;SULT2B1;GAB2;THSD4;RAPGEFL1;CELSR1                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| E2F Targets                       | -1.032    |  0.429688   | 0.609141    | CDC25B;HMGA1;CDKN2A;JPT1;MYC;H2AX;RRM2;MKI67;AURKB;DLGAP5;MYBL2;MTHFD2;CHEK2;SPC24;MCM5;BRCA2;ASF1B;DNMT1;TK1;MCM6;MELK;CCNB2;KIF4A;DONSON;NAP1L1;BIRC5;PLK4;HELLS;CDKN2C;TMPO;PLK1;HMGB2;ESPL1;LMNB1;KIF18B;TOP2A;ORC6;TRIP13;EZH2;CENPE;STMN1;ATAD2;CKS2;CKS1B;CDCA8;POLE4;CENPM;ANP32E;MCM4;MCM7;BUB1B;EXOSC8;SPC25;KPNA2;LYAR;DIAPH3;CDKN1A;E2F8;MMS22L;BARD1;CDKN3;HMMR;MXD3;UBE2T;XRCC6;TACC3;POLD1;RAD1;CDCA3;TUBB;EED;DDX39A;SNRPB                                                                                                                                                                                                                                                                                                                                                                                                             |
| Apical Surface                    | -1.03226  |  0.45839    | 0.629546    | PLAUR;IL2RG;IL2RB;LYN;CRYBG1;SULF2;SRPX;CX3CL1;EFNA5;LYPD3;DCBLD2;GATA3                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |
| PI3K/AKT/mTOR  Signaling          | -1.03913  |  0.435644   | 0.633136    | SFN;NOD1;IL2RG;TIAM1;LCK;DAPP1;TNFRSF1A;CXCR4;GRK2;FASLG;CAMK4;MYD88;PFN1;ITPR2;PRKCB;CFL1;STAT2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| Bile Acid Metabolism              |  1.1344   |  0.194286   | 0.637081    | DIO1;TFCP2L1;DIO2;CYP7B1;ACSL1;CYP46A1;PFKM;MLYCD;ALDH9A1;PHYH;HSD17B4;ALDH1A1;PEX11A;GSTK1;RETSAT;IDI1;APOA1;PEX12;HACL1;PEX1;SCP2;CROT;PNPLA8;PEX11G;EPHX2;SLC27A5;SLC23A2;PXMP2;NUDT12;RXRA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| Myogenesis                        | -1.04453  |  0.391403   | 0.640432    | EPHB3;CLU;RYR1;ITGB4;SCD;COL1A1;LSP1;MYBPH;DAPK2;PDLIM7;ERBB3;APOD;IGF1;BHLHE40;TGFB1;GSN;PDE4DIP;ADAM12;CTF1;COL3A1;AEBP1;KLF5;FGF2;DMD;FLII;SYNGR2;PTP4A3;FST;COL6A3;HBEGF;TNNI2;CNN3;SPHK1;MYL4;TPM3;NAV2;AK1;PLXNB2;COL6A2;TNNI1;SVIL;TPD52L1;DTNA;PPFIA4;MYLPF;AGRN;ITGA7;PTGIS                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| Pancreas Beta Cells               | -1.0515   |  0.445714   | 0.644495    | DPP4;PCSK2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| UV Response Dn                    | -1.05722  |  0.387471   | 0.652712    | RUNX1;HAS2;LDLR;MET;COL1A1;ANXA2;LTBP1;MMP16;PIK3CD;RND3;BHLHE40;IGFBP5;RGS4;MYC;COL11A1;COL1A2;COL3A1;CELF2;F3;FYN;SERPINE1;KCNMA1;SIPA1L1;RASA2;AMPH;DLC1;ANXA4;NOTCH2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| Xenobiotic Metabolism             | -0.996906 |  0.507151   | 0.655781    | PROS1;CYP2S1;APOE;SLC1A5;UPP1;TMEM176B;ABCC3;FAS;CFB;TNFRSF1A;IGF1;BCAT1;HSD11B1;EPHA2;IRF8;AQP9;SLC6A6;PTGES;CDO1;KYNU;HES6;HMOX1;FBLN1;VNN1;PSMB10;TDO2;GCKR;SERPINE1;ESR1;SLC6A12;MAN1A1;XDH;CYP27A1;ABCD2;RBP4;NPC1;CYP2E1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| Glycolysis                        | -0.986708 |  0.539313   | 0.6595      | B3GNT3;GALE;PKP2;MET;SDC3;QSOX1;TGFA;HS6ST2;VCAN;ELF3;DCN;NT5E;IER3;CHST2;SPAG4;CXCR4;PAM;KDELR3;DSC2;PDK3;ISG20;TGFBI;ANG;ANGPTL4;COL5A1;CD44;GLRX;PKM;AK4;CHST1;PYGL;DDIT4;CLDN9;TPBG;IL13RA1;SLC16A3;PPFIA4;AGRN;P4HA2;FUT8;B3GAT1;IGFBP3;EGLN3;STMN1;PC;EGFR;LDHA;TXN;NDST3;ERO1A;HS2ST1;MERTK;CENPA;CHST6;EFNA3                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| KRAS Signaling Dn                 | -1.00055  |  0.487805   | 0.668668    | CLDN16;KRT15;RYR1;KLK7;CCNA1;CD207;TGM1;SPTBN2;PLAG1;KLHDC8A;CHST2;RGS11;SYNPO;PDCD1;SLC6A14;PKP1;KRT5;CD40LG;NUDT11;LFNG;RSAD2;KCNMB1;KRT13;MAST3;PNMT;TCL1A;KCND1;LYPD3;IFNG;SIDT1;SKIL;CD80;HTR1D;SHOX2;TGFB2;GPR19;CCR8;PRODH;IL12B;NPHS1;EDN1;MTHFR;MEFV;YPEL1;SLC12A3;CACNA1F                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| Myc Targets V1                    |  1.02716  |  0.336634   | 0.679527    | VDAC3;AIMP2;POLD2;GOT2;CANX;CYC1;TOMM70;PGK1;STARD7;PRDX3;C1QBP;NOLC1;TUFM;SNRPD3;SLC25A3;EIF4E;MRPS18B;COX5A;GSPT1;XPOT;ACP1;EIF1AX;PHB;EIF4H;TXNL4A;EIF2S1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
| Fatty Acid Metabolism             |  1.06278  |  0.270677   | 0.70355     | ADIPOR2;TP53INP2;ACSS1;ACSL1;CA4;MLYCD;CD36;ACADM;CPT2;HADH;SUCLG1;HSDL2;ETFDH;ALDH9A1;DLD;HSD17B4;ALDH1A1;BCKDHB;ACADS;EHHADH;FH;CRAT;MDH2;BPHL;ECI1;SUCLA2;IDH3B;SDHD;ACO2;HMGCS1;SDHC;KMT5A;HSD17B10;RETSAT;IDI1;PDHA1;HADHB;HCCS;DLST;SUCLG2;ME1;RDH11;GCDH;ECHS1;HMGCL;ENO2;GSTZ1;HMGCS2;ALAD;ENO3;AQP7;NTHL1;GABARAPL1;CPT1A;ERP29;CA2;AUH;ECI2;AADAT;PDHB;PPARA;ALDOA                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
| mTORC1 Signaling                  | -0.917268 |  0.664062   | 0.793758    | CTSC;CORO1A;ITGB2;LDLR;SCD;SLC1A5;DAPP1;SLC2A3;BHLHE40;BCAT1;IGFBP5;PNP;SLC6A6;FGL2;CXCR4;TUBA4A;WARS1;NIBAN1;EDEM1;RRM2;SLC7A5;MTHFD2;BUB1;GLRX;PPP1R15A;AK4;SYTL2;DDIT4;SKAP2;NAMPT;ELOVL6;ARPC5L;DHCR24;SLC1A4;PLK1;DHCR7;EGLN3;MLLT11;SLC2A1;IFI30;LDHA;ERO1A;ACTR3;SERPINH1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| TGF-beta Signaling                | -0.889181 |  0.668027   | 0.827695    | TGFBR1;TGFB1;SLC20A1;RAB31;PMEPA1;THBS1;LTBP2;SERPINE1;TGIF1;PPP1R15A;CDK9;SMURF2;ACVR1;SKIL;JUNB;IFNGR2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| Adipogenesis                      |  1.16025  |  0.0934579  | 0.836095    | MTARC2;LIFR;FABP4;ADIPOR2;GPAM;NDUFA5;ADIG;COQ9;CD36;ACADM;DBT;CPT2;HADH;SUCLG1;CHCHD10;SCARB1;ALDH2;PHYH;DLD;QDPR;ACADS;AIFM1;RTN3;COQ5;COX8A;CRAT;MDH2;CYC1;HSPB8;ESRRA;DHRS7;ACO2;POR;SDHC;MRPL15;REEP5;TOB1;SPARCL1;RETSAT;GHITM;DLAT;FZD4;COX7B;CMBL;ME1;PRDX3;ABCB8;VEGFB;GPX3;ECHS1;CS;LPCAT3;IMMT;UQCR10;MCCC1;SCP2;NKIRAS1;GRPEL1;SLC19A1;MTCH2;PREB;COQ3;EPHX2;DHRS7B;SAMM50;CAVIN2;COX6A1;ALDOA;SDHB;DNAJB9;SLC25A1;IDH3A;ETFB;AGPAT3;UQCRQ;COL4A1;DECR1;NDUFS3                                                                                                                                                                                                                                                                                                                                                                             |
| heme Metabolism                   | -0.87098  |  0.757976   | 0.837773    | PDZK1IP1;C3;ACP5;CTSE;AQP3;ACKR1;ENDOD1;MBOAT2;RBM38;SLC25A37;OPTN;SLC22A4;TRAK2;CCND3;UCP2;MYL4                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| Wnt-beta Catenin Signaling        | -0.836373 |  0.740896   | 0.871823    | HEY2;TCF7;MYC;CSNK1E;CCND2;WNT5B;NKD1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| Protein Secretion                 | -0.780342 |  0.842893   | 0.905821    | CTSC;DST;PAM;PPT1;STAM;CLTA;VAMP4;SGMS1;IGF2R;SCRN1;CD63;EGFR;KRT18;ANP32E                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| Androgen Response                 | -0.795243 |  0.829152   | 0.910037    | KRT19;SCD;TSC22D1;PTK2B;ACTN1;PMEPA1;B2M;ALDH1A3;UBE2I;ARID5B;CCND3;STEAP4;MYL12A;SAT1;CCND1;ANKH;DHCR24;STK39;ADAMTS1;MAF                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| Myc Targets V2                    |  0.834593 |  0.71374    | 0.950592    | SORD;UNG;PRMT3;AIMP2;TFB2M;TMEM97;NOLC1;SLC19A1;NDUFAF4;DCTPP1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| Reactive Oxygen Species Pathway   |  0.74474  |  0.883803   | 0.957659    | IPCEF1;GCLC;HHEX;ERCC2;GPX3;NDUFS2;ATOX1;NDUFA6;TXNRD1;PRDX1;G6PD;TXNRD2;GLRX2;PRDX2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| Unfolded Protein Response         |  0.873436 |  0.683908   | 0.968385    | WFS1;PSAT1;HYOU1;HSPA9;SRPRB;HSP90B1;CALR;ATP6V0D1;PDIA6;EIF2AK3;VEGFA;HSPA5;NOLC1;ALDH18A1;GOSR2;CEBPG;EIF4E;SRPRA;DNAJC3;XPOT;PREB;CNOT4;STC2;WIPI1;SPCS3;DNAJB9;HERPUD1;POP4;EIF2S1;TTC37;SLC30A5;SEC31A                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
| Spermatogenesis                   | -0.698173 |  0.950704   | 0.977661    | CCNA1;PCSK1N;SCG5;GFI1;CAMK4;NCAPH;BUB1;RPL39L;OAZ3;SLC2A5;TLE4;CCNB2;EZH2;CNIH2;ADCYAP1;DCC;SHE;NEK2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| DNA Repair                        | -0.623306 |  0.998833   | 0.994716    | HCLS1;PNP;DDB2;AK1;RAD51;GPX4;ELL;POLR3C;STX3;POLE4;GTF2B;GUK1;ADA;NPR2;IMPDH2;NME3;PDE6G;POLD1;MPG;POLR2D;RFC2;SUPT5H;RBX1;CANT1;POLD4;LIG1;TYMS;CSTF3;CMPK2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| Pperoxisome                       |  0.895584 |  0.695402   | 1           | DIO1;CRABP1;ACSL1;HSD11B2;MLYCD;ALDH9A1;ABCC8;HSD17B4;SLC25A4;ALDH1A1;EHHADH;CRAT;PEX11A;GSTK1;RETSAT;IDI1;RDH11;HMGCL;SCP2;PRDX1;SLC25A17;ECI2;EPHX2;SLC23A2;CLN6                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |

### 2.3 Reactome top 30
- file: `/home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/F2_gsea_reactome_proper.tsv`
- DM2/positive NES head 15
| pathway.1                                                                                                                   |     NES |   NOM p-val |   FDR q-val | Lead_genes                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
|:----------------------------------------------------------------------------------------------------------------------------|--------:|------------:|------------:|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Respiratory Electron Transport R-HSA-611105                                                                                 | 1.70788 |  0          |  0.00114921 | NDUFA5;SCO1;NUBPL;CYCS;ETFDH;NDUFB2;UQCRC2;NDUFC2;COX8A;ECSIT;SDHD;UQCRFS1;NDUFS1;NDUFB5;COX7B;UQCRH;NDUFA8;NDUFA10;NDUFA1;NDUFS2;TRAP1;UQCR10;COX18;NDUFV1;NDUFA6;NDUFV3;NDUFB3;NDUFA12;COX5A;ETFA;NDUFB11;NDUFAF4;TACO1;COX6C;NDUFB6;COA1;NDUFS4;COX6A1;SDHB;COX7C;SDHA;TMEM186;LRPPRC;TIMMDC1;NDUFA2;ETFB;NDUFA9;UQCRQ;NDUFS3;UQCRC1;COX14;SURF1;NDUFAF1;COX5B;NDUFS6;NDUFB1;NDUFB10;NDUFB9;NDUFS7;COX19;NDUFC1;NDUFAF5;NDUFAF2;COX4I1                                                                                                                                                                             |
| Respiratory Electron Transport, ATP Synthesis By Chemiosmotic Coupling, Heat Production By Uncoupling Proteins R-HSA-163200 | 1.71622 |  0          |  0.00229843 | NDUFA5;SCO1;NUBPL;ATP5F1B;CYCS;ETFDH;NDUFB2;ATP5F1A;UQCRC2;ATP5MC1;NDUFC2;COX8A;ECSIT;DMAC2L;SDHD;ATP5PD;UQCRFS1;NDUFS1;ATP5PB;NDUFB5;ATP5MC3;COX7B;UQCRH;NDUFA8;NDUFA10;NDUFA1;NDUFS2;TRAP1;UQCR10;COX18;ATP5F1C;NDUFV1;NDUFA6;NDUFV3;NDUFB3;ATP5PF;NDUFA12;COX5A;ETFA;NDUFB11;NDUFAF4;TACO1;COX6C;NDUFB6;COA1;NDUFS4;COX6A1;SDHB;COX7C;SDHA;TMEM186;LRPPRC;TIMMDC1;NDUFA2;ETFB;NDUFA9;UQCRQ;NDUFS3;UQCRC1;COX14;SURF1;ATP5MC2;NDUFAF1;COX5B;NDUFS6;NDUFB1;NDUFB10;NDUFB9;ATP5PO;NDUFS7;COX19                                                                                                                        |
| Citric Acid (TCA) Cycle And Respiratory Electron Transport R-HSA-1428517                                                    | 1.6506  |  0          |  0.00689528 | LDHD;OGDH;NNT;NDUFA5;SCO1;L2HGDH;NUBPL;ATP5F1B;MPC1;SUCLG1;CYCS;ETFDH;NDUFB2;ATP5F1A;UQCRC2;ATP5MC1;NDUFC2;COX8A;MDH2;ECSIT;DMAC2L;SUCLA2;IDH3B;SDHD;ACO2;ATP5PD;UQCRFS1;NDUFS1;ATP5PB;PDHA1;NDUFB5;ATP5MC3;DLAT;DLST;COX7B;SUCLG2;UQCRH;NDUFA8;NDUFA10;CS;NDUFA1;PDPK1;NDUFS2;TRAP1;GSTZ1;UQCR10;COX18;ATP5F1C;NDUFV1;NDUFA6;NDUFV3;PDHX;NDUFB3;ATP5PF;NDUFA12;COX5A;ETFA;NDUFB11;NDUFAF4;TACO1;COX6C;NDUFB6;COA1;PDHB;NDUFS4;RXRA;COX6A1;SDHB;COX7C;SDHA;IDH3A;TMEM186;LRPPRC;TIMMDC1;NDUFA2;PDPR;ETFB;NDUFA9;UQCRQ;NDUFS3;PLPP6;UQCRC1;PDK2;COX14;SURF1;PDP2;LDHB;HAGH;ATP5MC2;NDUFAF1;COX5B;NDUFS6;NDUFB1;NDUFB10 |
| Protein Localization R-HSA-9609507                                                                                          | 1.63511 |  0          |  0.00900217 | LDHD;MLYCD;ATP5F1B;TIMM17B;CHCHD10;PHYH;ACBD5;HSD17B4;SLC25A4;ATP5F1A;EHHADH;HSPA9;ATP5MC1;CRAT;FXN;ACO2;TOMM70;GSTK1;ACOX2;TIMM10B;PEX12;TIMM21;PXMP4;SLC25A13;HMGCL;CS;CHCHD3;MTX1;CROT;TIMM17A;NUDT7;GRPEL1;MTX2;SLC25A17;ECI2;EPHX2;PXMP2;SAMM50;TYSND1;EMD;GRPEL2;ACOT2;TIMM22;PMPCB;ACOX3;SGTA;NUDT19;CYB5A;PEX10;DHRS4;CMC2;PITRM1;GET3;PEX5;DNAJC19;CHCHD2;GDAP1;TIMM44;ZNF260;COX19;PEX11B;TIMM23;ACOT4;UBE2D2;PEX2;ACOT8;SEC61G;TIMM50;PEX3;UBB;TIMM8A;TOMM40;HSPD1;SLC25A12;PMPCA;USP9X;TOMM20;GNPAT;COA6;COX17                                                                                            |
| Mitochondrial Protein Import R-HSA-1268020                                                                                  | 1.60793 |  0.00411523 |  0.0176213  | LDHD;ATP5F1B;TIMM17B;CHCHD10;SLC25A4;ATP5F1A;HSPA9;ATP5MC1;FXN;TOMM70;TIMM10B;TIMM21;SLC25A13;CS;CHCHD3;MTX1;TIMM17A;GRPEL1;MTX2;SAMM50;GRPEL2;TIMM22;PMPCB;CMC2;PITRM1;DNAJC19;CHCHD2;TIMM44;COX19;TIMM23;TIMM50;TIMM8A;TOMM40;HSPD1;SLC25A12;PMPCA;TOMM20;COA6;COX17                                                                                                                                                                                                                                                                                                                                                |
| Mitochondrial Biogenesis R-HSA-1592230                                                                                      | 1.57791 |  0.00515464 |  0.0357533  | PPARGC1A;ATP5F1B;CYCS;SIRT4;ATP5F1A;PRKAA2;APOO;HSPA9;ATP5MC1;ESRRA;DMAC2L;ATP5PD;ATP5PB;TFB2M;ATP5MC3;NCOA6;IMMT;ATP5F1C;APOOL;CHCHD3;MTX1;SIRT5;ATP5PF;MTX2;SAMM50;RXRA;PPARA;ALAS1;PRKAG2;DNAJC11;PRKAB2;PRKAG1;TBL1X;ATP5MC2;HDAC3;TMEM11                                                                                                                                                                                                                                                                                                                                                                         |
| Complex I Biogenesis R-HSA-6799198                                                                                          | 1.55886 |  0.00369004 |  0.0453118  | NDUFA5;NUBPL;NDUFB2;NDUFC2;ECSIT;NDUFS1;NDUFB5;NDUFA8;NDUFA10;NDUFA1;NDUFS2;NDUFV1;NDUFA6;NDUFV3;NDUFB3;NDUFA12;NDUFB11;NDUFAF4;NDUFB6;COA1;NDUFS4;TMEM186;TIMMDC1;NDUFA2;NDUFA9;NDUFS3;NDUFAF1;NDUFS6;NDUFB1;NDUFB10;NDUFB9;NDUFS7;NDUFC1;NDUFAF5;NDUFAF2                                                                                                                                                                                                                                                                                                                                                            |
| Mitochondrial Translation R-HSA-5368287                                                                                     | 1.48175 |  0.0157895  |  0.207369   | MRPS33;MRPS11;MRPS10;MRPS35;GFM1;MRPL15;GFM2;MRPS7;MRPL34;MRPS27;MRPS26;MTIF2;MRPS17;MRPS14;MRPL19;MRPL42;TUFM;ERAL1;MRPS18B;MRPL30;MRPL49;MRPS30;MRPS18A;OXA1L;MRPS22;MRPS36;MRPS16;MRPL45;MRRF;MRPS31;MTIF3;MRPS2;MRPS34;MTRF1L;MRPS9;MRPL44;MRPL24;MRPL16;MRPL13;MRPL32;MRPL36;MRPS15;MRPL12;PTCD3                                                                                                                                                                                                                                                                                                                 |
| Glyoxylate Metabolism And Glycine Degradation R-HSA-389661                                                                  | 1.48299 |  0          |  0.226586   | LDHD;ALDH4A1;GLDC;DHTKD1;OGDH;GCSH;DBT;GOT2;BCKDHB;HOGA1;PDHA1;DLAT;DLST;LIAS;LIPT2;PDHX;PDHB;PXMP2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| Mitochondrial Translation Initiation R-HSA-5368286                                                                          | 1.43437 |  0.0251256  |  0.335217   | MRPS33;MRPS11;MRPS10;MRPS35;MRPL15;MRPS7;MRPL34;MRPS27;MRPS26;MTIF2;MRPS17;MRPS14;MRPL19;MRPL42;ERAL1;MRPS18B;MRPL30;MRPL49;MRPS30;MRPS18A;OXA1L;MRPS22;MRPS36;MRPS16;MRPL45;MRPS31;MTIF3;MRPS2;MRPS34;MRPS9;MRPL44;MRPL24;MRPL16;MRPL13;MRPL32;MRPL36;MRPS15;MRPL12;PTCD3                                                                                                                                                                                                                                                                                                                                            |
| Intraflagellar Transport R-HSA-5620924                                                                                      | 1.44373 |  0.0132013  |  0.33564    | TTC30A;DYNLL2;TTC30B;HSPB11;IFT22;KIF3B;IFT140;DYNC2I1;DYNLRB2;IFT20;TNPO1;KIF17;DYNC2H1;TTC26;IFT88                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| Mitochondrial Translation Termination R-HSA-5419276                                                                         | 1.43765 |  0.0193237  |  0.343359   | MRPS33;MRPS11;MRPS10;MRPS35;MRPL15;GFM2;MRPS7;MRPL34;MRPS27;MRPS26;MRPS17;MRPS14;MRPL19;MRPL42;ERAL1;MRPS18B;MRPL30;MRPL49;MRPS30;MRPS18A;OXA1L;MRPS22;MRPS36;MRPS16;MRPL45;MRRF;MRPS31;MRPS2;MRPS34;MTRF1L;MRPS9;MRPL44;MRPL24;MRPL16;MRPL13;MRPL32;MRPL36;MRPS15;MRPL12;PTCD3                                                                                                                                                                                                                                                                                                                                       |
| Mitochondrial Translation Elongation R-HSA-5389840                                                                          | 1.44551 |  0.029703   |  0.357635   | MRPS33;MRPS11;MRPS10;MRPS35;GFM1;MRPL15;MRPS7;MRPL34;MRPS27;MRPS26;MRPS17;MRPS14;MRPL19;MRPL42;TUFM;ERAL1;MRPS18B;MRPL30;MRPL49;MRPS30;MRPS18A;OXA1L;MRPS22;MRPS36;MRPS16;MRPL45;MRPS31;MRPS2;MRPS34;MRPS9;MRPL44;MRPL24;MRPL16;MRPL13;MRPL32;MRPL36;MRPS15;MRPL12;PTCD3                                                                                                                                                                                                                                                                                                                                              |
| Metabolism Of Amine-Derived Hormones R-HSA-209776                                                                           | 1.4167  |  0.00266667 |  0.394257   | TPO;DIO2;IYD                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| Cristae Formation R-HSA-8949613                                                                                             | 1.41883 |  0.0309598  |  0.406657   | ATP5F1B;ATP5F1A;HSPA9;ATP5MC1;DMAC2L;ATP5PD;ATP5PB;ATP5MC3;IMMT;ATP5F1C;APOOL;CHCHD3;MTX1;ATP5PF;MTX2;SAMM50;DNAJC11;ATP5MC2;TMEM11;ATP5PO                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |

- DM1/negative NES head 15
| pathway.1                                                                             |      NES |   NOM p-val |   FDR q-val | Lead_genes                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
|:--------------------------------------------------------------------------------------|---------:|------------:|------------:|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Immunoregulatory Interactions Between A Lymphoid And A non-Lymphoid Cell R-HSA-198933 | -1.50539 |  0          |  0.00111748 | ICAM1;C3;ICAM5;JAML;ITGB2;SELL;TREM1;OSCAR;TYROBP;CD3E;CD1C;COLEC12;CD40;FCGR2B;LAIR1;TREM2;CD1E;HLA-G;CD3D;CD300LF;LILRB4;CLEC2D;HLA-B;CD33;HLA-F;PIANP;SLAMF7;ITGB7;HLA-A;SIGLEC10;SLAMF6;KLRB1;HLA-C;CD247;SIGLEC6;CD22;PILRA;SIGLEC9;LILRB2;ITGAL;IFITM1;CD300A;SH2D1A;ITGA4;MICB;LILRB1;CRTAM;CD40LG;ICAM3;CD8B;NPDC1;CD19;HLA-E;CD200R1;CD8A;KLRC1;FCGR1A;CD1B;LILRA6;CD300LB;CD96;VCAM1;RAET1E;LILRA2;TREML2;IGLL5;SIGLEC1;KLRD1;NCR3;SIGLEC7;CD226;HCST;CD200;LAIR2                                                                                                                                                                                                                                                                                                                                                                                                                           |
| Neutrophil Degranulation R-HSA-6798695                                                | -1.49355 |  0          |  0.00167622 | CHI3L1;SERPINA1;PLAU;ALOX5;SLPI;LRG1;LGALS3;PLAUR;CTSH;CTSC;C3;TMC6;STING1;ITGB2;LCN2;ALDH3B1;SELL;ITGAX;FCER1G;OSCAR;S100A11;S100A9;CEACAM6;CD55;PRSS2;ITGAM;QPCT;ADAM8;CD53;LYZ;ANXA2;TYROBP;UNC13D;ANPEP;FGR;SIRPA;ARHGAP9;OLR1;CHIT1;MNDA;ARHGAP45;HPSE;TNFRSF1B;RAB27A;LAIR1;ADGRE5;QSOX1;FCGR2A;NCKAP1L;PYCARD;KCNAB2;PTPRC;CLEC5A;HLA-B;BIN2;CYBB;S100A8;CTSS;SLC2A3;CD33;CD58;MMP9;PTAFR;CRACR2A;NPC2;DOCK2;HLA-A;TCIRG1;TLR2;C3AR1;TRPM2;HLA-C;CD47;NFAM1;DOK3;CLEC12A;GSN;SIGLEC9;FPR1;RAP2B;ATP8B4;LILRB2;ITGAL;CD300A;GPR84;FGL2;MMP12;RAB31;RHOG;MCEMP1;C5AR1;PKP1;CD177;ATP11A;CR1;SFTPA2;ALDOC;TIMP2;SLC15A4;FCN1;CEACAM21;B2M;SIGLEC14;SIRPB1;SLC11A1;HSPA6;VNN1;LILRA6;PLAC8;NBEAL2;ADA2;BST2;TCN1;DEGS1;SNAP25;TNFAIP6;VAMP8;TMEM63A;IQGAP1;STOM;PTPN6;CYBA;CRISPLD2;SERPINB1;MAN2B1;RETN;GLIPR1;P2RX1;RNASE2;CD44;CXCL1;CNN2;PKM;PIGR;CD59;SLC27A2;PYGL;CTSD;DLC1;GMFG;SLC2A5;ANO6 |
| Cell-Cell Communication R-HSA-1500931                                                 | -1.48067 |  0          |  0.00391119 | LAMB3;CLDN1;CLDN16;CLDN10;SDK1;CDH11;ITGB4;NECTIN4;CDH6;TYROBP;CLDN2;SIRPA;PTK2B;FLNA;DST;FYB1;SIRPG;CD47;CD151;VASP;CDH4;CADM3;CDH15;SRC;SIRPB1;PLEC;CDH24;CLDN4;NECTIN1;FYN;ANG;IQGAP1;PTPN6                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
| Cell Surface Interactions At Vascular Wall R-HSA-202733                               | -1.48073 |  0          |  0.00521491 | FN1;PROS1;SDC4;JAML;ITGB2;SELL;TREM1;ITGAX;FCER1G;S100A9;CEACAM6;SELPLG;SDC3;ITGAM;DOK2;CD2;CD74;SIRPA;OLR1;LCK;ITGA3;INPP5D;CD84;GRB7;CD58;TNFRSF10A;JCHAIN;SIRPG;CD47;CD48;LYN;TGFB1;ITGAL;ITGA4;CD177;SRC;SLC7A7;CEACAM21;PROCR;PROC;SLC7A5;FYN;IGLL5;PTPN6;CD244;CD44;TNFRSF10B                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| TCR Signaling R-HSA-202403                                                            | -1.46154 |  0.00121359 |  0.0136333  | HLA-DQB2;CARD11;HLA-DQA1;HLA-DRB1;HLA-DRA;HLA-DPA1;HLA-DRB5;CD3E;WAS;HLA-DQA2;PSMB9;LCK;CD3D;PTPRC;CD4;INPP5D;PSMB8;PTPN22;ZAP70;ITK;FYB1;LCP2;VASP;CD247;PAG1;CD101;GRAP2;CSK;PSMB10;TRAT1;RIPK2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |
| Cell Junction Organization R-HSA-446728                                               | -1.45732 |  0          |  0.014341   | LAMB3;CLDN1;CLDN16;CLDN10;SDK1;CDH11;ITGB4;NECTIN4;CDH6;CLDN2;FLNA;DST;CD151;VASP;CDH4;CADM3;CDH15;PLEC;CDH24;CLDN4;NECTIN1;ANG                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| Signaling By Interleukins R-HSA-449147                                                | -1.44674 |  0          |  0.0258617  | TIMP1;ALOX5;FN1;JAK3;ICAM1;IL1RAP;MUC1;ANXA1;DUSP6;STAT4;NOD1;PELI1;CXCL2;ITGB2;LCN2;S100B;EBI3;ITGAX;DUSP4;CCL22;S100A9;ITGAM;IL2RG;HCK;ANXA2;PTPN7;CCL5;CSF2RA;IL7R;CD86;PTGS2;CXCL8;CCR1;PTK2B;HAVCR2;TNFRSF1B;CCR2;STAT1;PSMB9;CCR5;LCK;VAV1;BATF;IL10RA;IL2RA;CD4;CSF2RB;INPP5D;CSF2;IL2RB;OSMR;PSMB8;NOD2;CRLF1;CSF1R;IL18;LGALS9;NFKB2;MMP9;PTAFR;CCL19;TNFRSF1A;IRAK2;IL4R;IL1RL1;CXCL10;IL1RN;IL6R;CCL2;IL12RB1;LYN;TGFB1;CCL20;CCM2;FPR1;IL21R;LIF;CASP1;CSF1;IL1R2;NLRC5;CLCF1;MYC;IL32;STAT5A;MSN;CCL4;IL7;CRLF2;COL1A2;CTF1;HMOX1;IL16;BCL2L1;CASP8;IL15RA;PSMB10;MMP2;CEBPD;RIPK2;VCAM1;RHOU;SNAP25;FASLG;FYN;MYD88;IL6;MCL1;CSF3R;PTPN6;IL27RA;SOCS3;IL1B;FCER2;PSME2;SOCS1;IL36RN;MAPKAPK3;CXCL1;CNN2;IL18RAP;IL18BP;IL23A;TWIST1;PSME1;STAT6;TEC;CAPZA1;IL22RA2;CCND1;IFNG;MAPKAPK2;TNF;GATA3;CD80;TANK;STXBP2;IL13RA1;IL18R1;IL24;PTPN5;PTPN2;BIRC5;GAB2;F13A1;MAP2K7               |
| Keratinization R-HSA-6805567                                                          | -1.44139 |  0          |  0.0307307  | KRT19;KRT15;PKP2;KRT17;IVL;DSC3;TGM1;CSTA;KRT80;KRT6A;PKP1;KRT5;PPL;KRT14;KRT13;PERP;KRT16;TCHH;KLK12;PKP4;KLK13;SPRR1B;KRTAP2-3;KRT6B                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |
| Antimicrobial Peptides R-HSA-6803157                                                  | -1.43878 |  0.00129199 |  0.0322828  | CLU;LCN2;S100A9;PRSS2;LYZ;GNLY;CCR2;CD4;S100A8;RNASE6;TLR2;SLC11A1;TRIM21;TLR1;BPIFB1;PLA2G2A                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| Interferon Signaling R-HSA-913531                                                     | -1.43586 |  0          |  0.0337479  | ICAM1;GBP2;HLA-DQB2;HLA-DQA1;HLA-DRB1;HLA-DRA;HLA-DPA1;HLA-DRB5;IRF5;CIITA;HLA-DQA2;HLA-G;STAT1;HLA-B;PSMB8;FLNA;TRIM14;HLA-F;PTAFR;GBP5;IFIT5;SAMHD1;HLA-A;HLA-C;TRIM22;IRF8;IFITM1;IFIT2;IRF1;OASL;TRIM46;IFIT3;UBE2L6;MX2;SP100;HLA-E;B2M;RSAD2;FCGR1A;ISG15;HERC5;VCAM1;IFNAR2;BST2;TRIM21;IRF9;ISG20;IFNGR1;IRF3;OAS1;PTPN6;IRF7;SOCS3;GBP1;OAS3;SOCS1;CD44;PML;IFITM2;TRIM38;UBA7;IFNG;IFITM3;GBP6;PTPN2;IFIT1;IFNGR2;TRIM29;STAT2;XAF1;GBP4;IFI35;OAS2;IRF4                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| TNFR2 Non-Canonical NF-kB Pathway R-HSA-5668541                                       | -1.43344 |  0.00121655 |  0.0360642  | TNFRSF12A;LTB;BIRC3;TNFRSF18;CD40;TNFRSF1B;PSMB9;PSMB8;NFKB2;TNFRSF1A;TNFSF14;TNFRSF25;TNFSF11;TNFSF13B;RELB;TNFSF9;CD70;TNFSF4;TNFRSF9;MAP3K14;CD40LG;LTA;PSMB10;TNFSF8;TNFSF15;FASLG;TNFRSF13C;TNFSF18;PSME2;TNFRSF13B;TNFRSF8;TNFRSF11A;TNFRSF17;PSME1;TNFSF12;TNF;TANK                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
| Regulation Of IGF Transport And Uptake By IGFBPs R-HSA-381426                         | -1.42853 |  0.00352113 |  0.0428368  | SERPINA1;IGFBP6;TIMP1;FN1;FAM20A;C3;PRSS23;LGALS1;EVA1A;APOE;MFGE8;MXRA8;FSTL3;APOL1;FAM20C;QSOX1;LTBP1;VCAN;F5;IGF1;IGFBP5;CSF1;GZMH;C4A;CP;CHRDL1;MELTF;MMP2;PROC;SHISA5;IL6                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
| Interferon Alpha/Beta Signaling R-HSA-909733                                          | -1.42364 |  0          |  0.0504463  | GBP2;IRF5;HLA-G;STAT1;HLA-B;PSMB8;HLA-F;IFIT5;SAMHD1;HLA-A;HLA-C;IRF8;IFITM1;IFIT2;IRF1;OASL;IFIT3;MX2;HLA-E;RSAD2;ISG15;IFNAR2;BST2;IRF9;ISG20;IRF3;OAS1;PTPN6;IRF7;SOCS3;OAS3;SOCS1;IFITM2;IFITM3;IFIT1;STAT2;XAF1;IFI35;OAS2;IRF4                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| Antigen processing-Cross Presentation R-HSA-1236975                                   | -1.42226 |  0          |  0.0507337  | MRC2;S100A9;NCF2;MRC1;NCF4;CD207;HLA-G;PSMB9;HLA-B;CYBB;S100A8;CTSS;BTK;PSMB8;HLA-F;LY96;HLA-A;NCF1;TLR2;HLA-C;TLR6;HLA-E;FCGR1A;PSMB10;VAMP8;MYD88;CYBA;TLR1;PSME2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| Class A/1 (Rhodopsin-like Receptors) R-HSA-373076                                     | -1.42397 |  0          |  0.0530374  | CYSLTR2;C3;ANXA1;LPAR5;CXCL2;CCL13;CCL22;P2RY6;CCL5;CXCL8;SSTR3;CXCL16;CCR1;CCR2;PTGER4;CCR5;FPR3;CCL17;CXCL9;PTAFR;CCL19;S1PR2;NLRP3;ADORA1;CXCR6;CXCL10;GPNMB;CXCL13;C3AR1;P2RY10;CCL2;CCR7;GPR68;CCL20;FPR1;ECE1;P2RY13;NPFFR1;FFAR4;CXCR4;GPR55;ADRA2A;INSL3;CXCL11;C5AR1;CCL4;CMKLR1;SFTPA2;S1PR4;GPR65;CX3CL1;CD200R1;GPR37;GPR183;CXCR3;RIPK2;XCL1;PTGER2;ADRB2;XCR1;CXCL3;NMU;PTGER1;PLPPR4;P2RY12;GPR18;XCL2;CCL21;CXCL5;ADORA3;CCL23;CXCL1;KISS1R;ADRA1B;EDNRA;PTGIR;PNOC;CCRL2;HTR1D;LTB4R;ACKR4;UTS2;HTR7;C5AR2;HCAR2;TACR1;GPBAR1;LPAR6;CCR8;SAA2;PLPPR1;OXTR;PTGDR;F2RL2;CCR3;FPR2;CXCR2;ADRB3;EDN1;CCL7;LHB;NPY1R;PPBP;S1PR3                                                                                                                                                                                                                                                           |

### 2.4 KEGG top 20
- file: `/home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/F2_gsea_kegg_proper.tsv`
| pathway.1                                                     |      NES |   NOM p-val |   FDR q-val | Lead_genes                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
|:--------------------------------------------------------------|---------:|------------:|------------:|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Staphylococcus aureus infection                               | -1.5271  |  0          |  0          | ICAM1;C3;KRT19;KRT15;HLA-DQA1;ITGB2;HLA-DRB1;CFH;HLA-DOA;KRT17;SELPLG;ITGAM;HLA-DPB1;HLA-DRA;HLA-DPA1;HLA-DQB1;HLA-DRB5;CFI;FCGR2B;HLA-DQA2;HLA-DOB;FCGR2A;FPR3;HLA-DMA;C2;C1S;CFB;C1QB;PTAFR;C1R;C1QC;C1QA;C3AR1;FCGR3A;FPR1;ITGAL;C4B;FCGR2C;SELP;C5AR1;C4A;KRT14;FCGR1A;KRT13;KRT16                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| Cell adhesion molecules                                       | -1.49814 |  0          |  0.00190547 | CLDN1;ICAM1;CLDN16;CLDN10;SDC4;CDH3;HLA-DQA1;ITGB2;HLA-DRB1;SELL;HLA-DOA;SELPLG;SDC3;ITGAM;HLA-DPB1;SPN;HLA-DRA;HLA-DPA1;HLA-DQB1;CD2;HLA-DRB5;CLDN2;CD86;ITGB8;CD40;HLA-DQA2;HLA-DOB;HLA-G;CD6;PTPRC;HLA-B;CD4;VCAN;HLA-DMA;TIGIT;CD58;HLA-F;CTLA4;ITGB7;HLA-A;HLA-C;CD276;CD22;CNTNAP1;ITGAL;CD28;ICOS;PDCD1;ITGA4;SELP;CDH4;CADM3;NLGN2;NRCAM;CD40LG;ICAM3;CD8B;CDH15;HLA-E;PDCD1LG2;VTCN1;CD8A;VCAM1;CLDN4;NECTIN1;SIGLEC1;CD274;CD226;HLA-DMB;CLDN11;PTPRF;CLDN9;VSIR;CD80;NCAM2;CLDN7;NTNG2;ITGA9                                                                                                                                                                                                                                                                                     |
| Cytokine-cytokine receptor interaction                        | -1.46889 |  0          |  0.00647861 | CXCL14;IL1RAP;CCL18;CXCL2;CCL13;TNFRSF12A;LTB;EBI3;INHBB;CXCL17;CCL22;TNFRSF18;IL2RG;CCL5;CSF2RA;IL7R;NGFR;CXCL8;CXCL16;CD40;CCR1;TNFRSF1B;CCR2;CCR5;IL10RA;IL2RA;TNFRSF21;CD4;CSF2RB;CSF2;IL2RB;CCL17;OSMR;CSF1R;CXCL9;FAS;IL18;IFNE;CCL19;TGFBR1;TNFRSF1A;TNFSF14;TNFRSF10A;IL4R;CXCR6;IL1RL1;CXCL10;TNFRSF25;IL1RN;CXCL13;TNFSF11;IL6R;CCR4;CCL2;TNFSF13B;TNFSF9;CD70;CCR7;IL12RB1;TGFB1;CCL20;TNFSF4;IL21R;LIF;RELT;TNFRSF9;INHBA;CSF1;IL1R2;GDF15;CXCR4;CLCF1;IL32;CXCL11;CCL4;IL7;CRLF2;CD40LG;LEPR;OSM;CTF1;CX3CL1;LTA;IL16;IL15RA;TNFSF8;CXCR3;XCL1;TNFSF15;TGFB3;IFNAR2;XCR1;CXCL3;FASLG;IFNGR1;IL6;CSF3R;TNFRSF13C;IL27RA;CCL4L2;IL1B;TNFSF18;XCL2;CCL21;CXCL5;CCL23;IL36RN;CXCL1;TNFRSF13B;TNFRSF10B;TNFRSF8;IL18RAP;TNFRSF11A;TNFRSF17;IL23A;TNFSF12;CCL24;ACVR1;TNFRSF10C;IFNG |
| Complement and coagulation cascades                           | -1.47575 |  0          |  0.00666915 | SERPINA1;PLAU;PLAUR;PROS1;C3;CLU;ITGB2;CFH;ITGAX;CD55;ITGAM;CFI;C7;C2;F5;C1S;CFB;SERPING1;C1QB;C1R;C1QC;C1QA;VSIG4;C3AR1;C4B;C5AR1;CR1;C4A;PROCR;CR2;PROC;SERPINF2;F3;SERPINE1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| Systemic lupus erythematosus                                  | -1.47136 |  0          |  0.00666915 | C3;HLA-DQA1;HLA-DRB1;HLA-DOA;HLA-DPB1;HLA-DRA;HLA-DPA1;HLA-DQB1;HLA-DRB5;CD86;CD40;HLA-DQA2;HLA-DOB;FCGR2A;C7;HLA-DMA;C2;C1S;C1QB;C1R;C1QC;C1QA;ACTN1;FCGR3A;CD28;C4B;CD40LG;C4A;H2AX;FCGR1A;TRIM21;H2AW                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| Natural killer cell mediated cytotoxicity                     | -1.46151 |  0          |  0.00762189 | RAC2;ICAM1;ITGB2;FCER1G;BID;TYROBP;PTK2B;HLA-G;LCK;VAV1;PIK3CD;HLA-B;CSF2;ZAP70;FAS;SH3BP2;TNFRSF10A;HLA-A;LCP2;HLA-C;CD48;CD247;FCGR3A;ITGAL;SH2D1A;MICB;GZMB;HLA-E;SHC3;KLRC1;IFNAR2;PRF1;FASLG;FYN;RAET1E;IFNGR1;KLRD1;PTPN6;CD244;NCR3;ULBP2;HCST;TNFRSF10B;IFNG;RAET1G;TNF;PRKCB;ULBP1;KLRK1;IFNGR2;CASP3                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| Phagosome                                                     | -1.46467 |  0          |  0.00825705 | MRC2;C3;HLA-DQA1;CORO1A;ITGB2;MARCO;HLA-DRB1;HLA-DOA;NCF2;ITGAM;HLA-DPB1;HLA-DRA;MRC1;HLA-DPA1;NCF4;HLA-DQB1;COMP;HLA-DRB5;CLEC7A;COLEC12;RAB7B;OLR1;FCGR2B;HLA-DQA2;HLA-DOB;HLA-G;FCGR2A;HLA-B;CYBB;TUBA1A;CTSS;HLA-DMA;HLA-F;C1R;SFTPA1;HLA-A;TCIRG1;NCF1;TLR2;HLA-C;MSR1;FCGR3A;TAP1;FCGR2C;TUBA4A;ITGA2;TUBB6;SFTPA2;CD14;ATP6V1B2;TAP2;TLR6;HLA-E;THBS2;THBS1;FCGR1A                                                                                                                                                                                                                                                                                                                                                                                                                   |
| Transcriptional misregulation in cancer                       | -1.45096 |  0          |  0.00831479 | RUNX2;PLAU;HMGA2;RUNX1;DUSP6;RXRG;SPI1;BIRC3;MET;ITGAM;NGFR;CD86;CXCL8;CCNA1;BAIAP3;ETV7;CD40;BCL2A1;GRIA3;CSF2;IL2RB;CSF1R;MMP9;ITGB7;IGF1;ETV4;NFKBIZ;IL1R2;MYC;GZMB;TRAF1;CD14;ETV5;BCL2L1;FCGR1A;SPINT1;DDB2;PAX5;ETV1;FLT3;IL6;GADD45A;RARA;CCND2;SLC45A3;PML;SSX1;CDK9;REL;BAX;CCNA2;LYL1;RELA;CDKN2C;FUT8;IGFBP3;TCF3;MAF;PTCRA;BCL6                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| Th17 cell differentiation                                     | -1.45733 |  0.00240674 |  0.00869372 | JAK3;IL1RAP;RUNX1;RXRG;HLA-DQA1;HLA-DRB1;HLA-DOA;AHR;HLA-DPB1;HLA-DRA;HLA-DPA1;HLA-DQB1;IL2RG;HLA-DRB5;CD3E;HLA-DQA2;HLA-DOB;NFKBIE;STAT1;LCK;CD3D;IL2RA;CD4;IL2RB;HLA-DMA;ZAP70;FOXP3;TGFBR1;IL4R;IL6R;CD3G;IL12RB1;CD247;TGFB1;IL21R;STAT5A;MAPK13;IFNGR1;IL6;IL27RA;IL1B;RARA;HLA-DMB;IL23A;STAT6;IFNG;TBX21;GATA3                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| Influenza A                                                   | -1.45124 |  0          |  0.00895572 | TMPRSS4;ICAM1;HLA-DQA1;HLA-DRB1;HLA-DOA;PRSS2;HLA-DPB1;HLA-DRA;BID;HLA-DPA1;HLA-DQB1;HLA-DRB5;CCL5;CXCL8;CIITA;HLA-DQA2;HLA-DOB;STAT1;PYCARD;PIK3CD;HLA-DMA;FAS;IL18;TNFRSF1A;TNFRSF10A;NLRP3;CXCL10;PRSS1;CCL2;TPSB2;CASP1;TPSAB1;IKBKE;MX2;RSAD2;CASP8;IFNAR2;IRF9;FASLG;IFNGR1;IRF3;OAS1;MYD88;IL6;IRF7;TLR7;SOCS3;IL1B;OAS3;CCND3;HLA-DMB;TNFRSF10B;PML;IFNG;TICAM1;TNF;BAX;PRKCB;DDX58;IFIH1;NXF1;ACTB;IFNGR2;RELA;CASP3;STAT2                                                                                                                                                                                                                                                                                                                                                         |
| Chemokine signaling pathway                                   | -1.45247 |  0          |  0.00952736 | CXCL14;JAK3;RAC2;CCL18;CXCL2;CCL13;CCL22;HCK;CCL5;FGR;CXCL8;TIAM1;CXCL16;ADCY8;CCR1;WAS;PTK2B;CCR2;PLCB2;STAT1;ADCY7;CCR5;VAV1;PIK3CD;PIK3R5;CCL17;GNB4;CXCL9;CCL19;ITK;DOCK2;CXCR6;CXCL10;NCF1;CXCL13;CCR4;CCL2;CCR7;LYN;CCL20;ARRB2;CXCR4;CXCL11;CCL4;SRC;GRK2;CX3CL1;PIK3R6;SHC3;PIK3CG;GNAI2;CXCR3;XCL1;XCR1;CXCL3;GRK6;CCL4L2;XCL2;CCL21;CXCL5;CCL23;CXCL1;GNGT2                                                                                                                                                                                                                                                                                                                                                                                                                       |
| Rheumatoid arthritis                                          | -1.43982 |  0.00124688 |  0.0122268  | ICAM1;ACP5;HLA-DQA1;CXCL2;ITGB2;HLA-DRB1;LTB;HLA-DOA;HLA-DPB1;HLA-DRA;HLA-DPA1;HLA-DQB1;HLA-DRB5;CCL5;CD86;CXCL8;HLA-DQA2;HLA-DOB;CSF2;HLA-DMA;IL18;CTSK;CTLA4;TCIRG1;TLR2;TNFSF11;CCL2;TNFSF13B;TGFB1;CCL20;ITGAL;CD28;CSF1;ATP6V1B2;TGFB3;CXCL3;IL6;IL1B;CXCL5;CXCL1;HLA-DMB;TNFRSF11A;IL23A;IFNG;TNF;CD80                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |
| Leukocyte transendothelial migration                          | -1.43806 |  0          |  0.012752   | CLDN1;RAC2;ICAM1;CLDN16;CLDN10;ITGB2;NCF2;ITGAM;NCF4;CLDN2;PTK2B;VAV1;PIK3CD;CYBB;RASSF5;MMP9;ITK;RHOH;SIPA1;ACTN1;NCF1;VASP;ITGAL;ITGA4;CXCR4;MSN;MAPK13;GNAI2;MMP2;VCAM1;CLDN4;TXK;CYBA;CLDN11;MYL12B;MYL12A;CLDN9;MYL9;PRKCB;CLDN7                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| Leishmaniasis                                                 | -1.43047 |  0          |  0.0147674  | C3;HLA-DQA1;ITGB2;HLA-DRB1;HLA-DOA;NCF2;ITGAM;HLA-DPB1;HLA-DRA;HLA-DPA1;NCF4;HLA-DQB1;HLA-DRB5;PTGS2;HLA-DQA2;HLA-DOB;FCGR2A;STAT1;CYBB;HLA-DMA;NCF1;TLR2;TGFB1;FCGR3A;ITGA4;FCGR2C;CR1;MAPK13;FCGR1A;TGFB3;IFNGR1;MYD88;PTPN6;CYBA;IL1B;HLA-DMB;IFNG;TNF;ELK1;PRKCB;TGFB2;IFNGR2;RELA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| Hematopoietic cell lineage                                    | -1.43068 |  0.00123762 |  0.0154978  | HLA-DQA1;HLA-DRB1;HLA-DOA;CD55;ITGAM;HLA-DPB1;HLA-DRA;HLA-DPA1;HLA-DQB1;CD2;HLA-DRB5;CD37;CD3E;CD1C;CSF2RA;IL7R;ANPEP;HLA-DQA2;HLA-DOB;CD1E;CD5;CD3D;CD1A;CD7;IL2RA;CD4;ITGA3;CSF2;HLA-DMA;CSF1R;CD33;IL4R;MS4A1;IL6R;CD3G;CD22;CSF1;IL1R2;ITGA4;ITGA2;IL7;CR1;CD14;CD8B;CD19;CD8A;FCGR1A;CD1B;CR2;FLT3;IL6;CSF3R;IL1B;FCER2;CD44;HLA-DMB;MME;CD59                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| Epstein-Barr virus infection                                  | -1.43205 |  0          |  0.0156521  | JAK3;ICAM1;HLA-DQA1;HLA-DRB1;HLA-DOA;HLA-DPB1;HLA-DRA;BID;HLA-DPA1;HLA-DQB1;HLA-DRB5;CD3E;CCNA1;TRAF5;CD40;HLA-DQA2;HLA-DOB;NFKBIE;HLA-G;STAT1;PIK3CD;CD3D;HLA-B;BTK;HLA-DMA;FAS;NFKB2;CD58;HLA-F;CXCL10;HLA-A;TLR2;HLA-C;CD3G;RELB;CD247;LYN;RUNX3;ITGAL;ENTPD8;TAP1;IKBKE;MYC;MAP3K14;TAP2;CD19;HLA-E;B2M;MAPK13;ENTPD1;ISG15;CASP8;CR2;DDB2;IFNAR2;IRF9;IRF3;OAS1;MYD88;IL6;GADD45A;IRF7;FCER2;TNFAIP3;OAS3;CCND3;CD44;CCND2;HLA-DMB                                                                                                                                                                                                                                                                                                                                                     |
| Coronavirus disease                                           | -1.42788 |  0.0010917  |  0.0161965  | C3;STING1;CXCL8;FCGR2A;STAT1;PIK3CD;CYBB;C7;CSF2;C2;C1S;CFB;C1QB;C1R;C1QC;C1QA;TNFRSF1A;NLRP3;CXCL10;TLR2;C3AR1;IL6R;CCL2;CASP1;C4B;IKBKE;SELP;TLR8;C5AR1;MX2;C4A;CGAS;MAPK13;ISG15;IFNAR2;RPLP1;IRF9;HBEGF;IRF3;OAS1;MYD88;IL6;TLR7;IL1B;RPL36;OAS3;RPS9;RPL22L1;RPL12;RPS2;RPL28;RPS4X;TNF;RPLP2;PRKCB;DDX58;IFIH1;RPS8;F13A1;RPL14;RPS20;RPSA;RELA;RPS19;STAT2;RPS3;OAS2;RPL24;STAT3;IL12B;RPS27;RPL35;RPL7A;CFD;MAPK11;RPL32;EGFR;PIK3R1;MMP1;C6;RPL13;RPL13A;RPLP0;RPL4;RPS15;MX1;RPL3;RPL8;IFNA1;RPL18;CSF3;RPS29;RPS18;RPL37A;MMP3;RPL34;FAU;TYK2;RPS17;IL6ST;PLCG2;RPL10A;RPS28                                                                                                                                                                                                     |
| Viral protein interaction with cytokine and cytokine receptor | -1.42643 |  0          |  0.0163553  | CXCL14;CCL18;CXCL2;CCL13;CCL22;IL2RG;CCL5;CXCL8;CCR1;TNFRSF1B;CCR2;CCR5;IL10RA;IL2RA;IL2RB;CCL17;CSF1R;CXCL9;IL18;CCL19;TNFRSF1A;TNFSF14;TNFRSF10A;CXCL10;CXCL13;IL6R;CCR4;CCL2;CCR7;CCL20;CSF1;CXCR4;CXCL11;CCL4;CX3CL1;LTA;CXCR3;XCL1;XCR1;CXCL3;IL6;CCL4L2;XCL2;CCL21;CXCL5;CCL23;CXCL1;TNFRSF10B;IL18RAP;CCL24;TNFRSF10C;TNF;IL18R1;IL24;CCL26;ACKR4                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| Inflammatory bowel disease                                    | -1.42414 |  0          |  0.0176005  | HLA-DQA1;STAT4;HLA-DRB1;HLA-DOA;HLA-DPB1;HLA-DRA;HLA-DPA1;HLA-DQB1;IL2RG;HLA-DRB5;HLA-DQA2;HLA-DOB;STAT1;HLA-DMA;NOD2;IL18;FOXP3;IL4R;TLR2;IL12RB1;TGFB1;IL21R;TGFB3;IFNGR1;TLR5;IL6;IL1B;HLA-DMB;IL18RAP;IL23A;STAT6;IFNG;TNF;TBX21;GATA3;IL18R1;TGFB2;IFNGR2;RELA;MAF;IL10;STAT3;IL12B                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| Th1 and Th2 cell differentiation                              | -1.4107  |  0.00251256 |  0.0286297  | JAK3;HLA-DQA1;STAT4;HLA-DRB1;HLA-DOA;HLA-DPB1;HLA-DRA;HLA-DPA1;HLA-DQB1;IL2RG;HLA-DRB5;CD3E;HLA-DQA2;HLA-DOB;NFKBIE;STAT1;LCK;CD3D;IL2RA;CD4;IL2RB;HLA-DMA;ZAP70;IL4R;CD3G;IL12RB1;CD247;RUNX3;STAT5A;MAPK13;IFNGR1;HLA-DMB;STAT6;IFNG;TBX21;GATA3;NOTCH2;MAML2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |

### 2.5 Proxy vs proper agreement
- file: `/home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/F2_proxy_vs_proper_comparison.tsv`
| pathway_name   | NES_proper   | FDR q-val_proper   | Term   | NES_proxy   | FDR q-val_proxy   | sign_agree   |
|----------------|--------------|--------------------|--------|-------------|-------------------|--------------|

## SEC 3. F3 Clinical reframe — RAW

### 3.1 시도된 endpoint 전체
- file: `/home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/F3_alternative_endpoints.tsv`
| endpoint                |   effect_dm2_minus_dm1 |      pvalue |   n_valid |
|:------------------------|-----------------------:|------------:|----------:|
| age                     |             12.8901    | 6.1852e-07  |       174 |
| stage_ord               |             -0.0779703 | 0.592739    |       165 |
| rai_score_recalc        |              1.59578   | 6.26563e-24 |       178 |
| tds16_score             |              1.3134    | 1.83254e-23 |       178 |
| histology_cPTC_vs_FVPTC |             86         | 4.16716e-08 |       178 |

- expected columns requested by prompt: `endpoint, n_total, n_events_or_positive, test_method, statistic, p_value, effect_size, effect_size_CI`
- MISSING columns in file: `['n_total', 'n_events_or_positive', 'test_method', 'statistic', 'p_value', 'effect_size', 'effect_size_CI']`
### 3.2 Logistic regression 전체 결과
- file: `/home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/F3_logistic_regression_cluster.tsv`
| feature                        |       coef |
|:-------------------------------|-----------:|
| num__rai_score_recalc          |  1.34629   |
| num__age                       |  0.977069  |
| cat__histology_subtype_FVPTC   |  0.7581    |
| cat__histology_subtype_cPTC    | -0.742405  |
| num__tds16_score               |  0.730597  |
| num__stage_ord                 | -0.440194  |
| cat__sex_Male                  |  0.0217735 |
| cat__sex_Female                | -0.0206043 |
| cat__histology_subtype_unknown | -0.0145254 |

- requested extra fields: `SE`, `p`, `AIC`, `BIC`, `pseudo-R2`
- MISSING columns in file: `['SE', 'p', 'AIC', 'BIC', 'pseudo-R2']`
### 3.3 Multivariable Cox (가능한 경우)
MISSING: no Phase 3 Cox output file found under `/home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables`
### 3.4 Effect size summary
- file: `/home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/F3_effect_size_summary.tsv`
| endpoint                |   effect_dm2_minus_dm1 |      pvalue |   n_valid |   abs_effect |
|:------------------------|-----------------------:|------------:|----------:|-------------:|
| age                     |             12.8901    | 6.1852e-07  |       174 |   12.8901    |
| stage_ord               |             -0.0779703 | 0.592739    |       165 |    0.0779703 |
| rai_score_recalc        |              1.59578   | 6.26563e-24 |       178 |    1.59578   |
| tds16_score             |              1.3134    | 1.83254e-23 |       178 |    1.3134    |
| histology_cPTC_vs_FVPTC |             86         | 4.16716e-08 |       178 |   86         |

## SEC 4. A1 scRNA heterogeneity — RAW

### 4.1 환자별 dominance
- file: `/home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/A1_per_patient_dominance.tsv`
| patient_id   |   n_cells |   dm1_dominance |   dm2_dominance |   dominance_ratio |
|:-------------|----------:|----------------:|----------------:|------------------:|
| Patient1     |      5917 |        0.995606 |      0.00439412 |          0.995606 |
| Patient10    |     22424 |        0.931502 |      0.068498   |          0.931502 |
| Patient2     |      5102 |        0.986476 |      0.0135241  |          0.986476 |
| Patient3     |      7436 |        0.954142 |      0.045858   |          0.954142 |
| Patient5     |      4477 |        0.903954 |      0.0960465  |          0.903954 |
| Patient8     |     10550 |        0.756682 |      0.243318   |          0.756682 |
| Patient9     |     10109 |        0.773766 |      0.226234   |          0.773766 |

### 4.2 Patient feature correlation
- file: `/home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/A1_dominance_vs_patient_features.tsv`
| patient_id   |   n_cells |   dm1_dominance |   dm2_dominance |   dominance_ratio | dominant_label   |   rank |
|:-------------|----------:|----------------:|----------------:|------------------:|:-----------------|-------:|
| Patient1     |      5917 |        0.995606 |      0.00439412 |          0.995606 | DM1              |      1 |
| Patient10    |     22424 |        0.931502 |      0.068498   |          0.931502 | DM1              |      2 |
| Patient2     |      5102 |        0.986476 |      0.0135241  |          0.986476 | DM1              |      3 |
| Patient3     |      7436 |        0.954142 |      0.045858   |          0.954142 | DM1              |      4 |
| Patient5     |      4477 |        0.903954 |      0.0960465  |          0.903954 | DM1              |      5 |
| Patient8     |     10550 |        0.756682 |      0.243318   |          0.756682 | DM1              |      6 |
| Patient9     |     10109 |        0.773766 |      0.226234   |          0.773766 | DM1              |      7 |

### 4.3 Cell type per cluster
- file: `/home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/A1_cell_type_per_cluster.tsv`
| cell_type   | dominant_dm   |     n |
|:------------|:--------------|------:|
| epithelial  | DM1           | 24518 |
| epithelial  | DM2           |  6240 |
| immune      | DM1           | 31089 |
| immune      | DM2           |   543 |
| stromal     | DM1           |  3152 |
| stromal     | DM2           |   473 |

## SEC 5. A2 BRAF/RAS-positive 안의 DM stratification — RAW

### 5.1 전체 cohort DM score
- file: `/home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/A2_dm_score_full_cohort.tsv`
| sample_id        | dataset   | platform                         | modality     | tissue_type   | normal_vs_tumor   | histology_subtype   | driver_anchor   | molecular_subtype   | tds_group   | dediff_flag   | aggressive_flag   | ajcc_stage_group   |   ata_risk_proxy |   bethesda_category | clinical_subtype_tag                    | label_confidence   | sex    |     age | outcome_available   | survival_available   | recurrence_available   | source_note                                  |   tds_score |   dedifferentiation_proxy_score |   brs_like_surrogate_score | v3_anchor         |    prob_dm2 |    prob_dm1 | dm_like   |   rai_score_recalc |
|:-----------------|:----------|:---------------------------------|:-------------|:--------------|:------------------|:--------------------|:----------------|:--------------------|:------------|:--------------|:------------------|:-------------------|-----------------:|--------------------:|:----------------------------------------|:-------------------|:-------|--------:|:--------------------|:---------------------|:-----------------------|:---------------------------------------------|------------:|--------------------------------:|---------------------------:|:------------------|------------:|------------:|:----------|-------------------:|
| TCGA-DJ-A2Q6-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | low         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 38.8583 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.679781  |                       0.679781  |                        nan | mutation_verified | 1.43348e-05 | 0.999986    | DM1_like  |            6.5982  |
| TCGA-FK-A3SE-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 31.4086 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   0.0449715 |                      -0.0449715 |                        nan | histology_proxy   | 6.10355e-05 | 0.999939    | DM1_like  |            7.919   |
| TCGA-DJ-A2QA-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | BRAF            | BRAF_like           | mid         | no            | no                | Stage III          |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Female | 56.9035 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   0.262867  |                      -0.262867  |                        nan | mutation_verified | 0.00022033  | 0.99978     | DM1_like  |            7.92517 |
| TCGA-FY-A2QD-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | high        | no            | no                | nan                |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 61.577  | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   1.70812   |                      -1.70812   |                        nan | histology_proxy   | 0.99996     | 3.95746e-05 | DM2_like  |            9.55679 |
| TCGA-EL-A3GR-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 31.7536 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   0.180948  |                      -0.180948  |                        nan | mutation_verified | 0.000111032 | 0.999889    | DM1_like  |            7.7278  |
| TCGA-EL-A3T7-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | low         | no            | no                | Stage II           |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 47.3703 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -1.0941    |                       1.0941    |                        nan | mutation_verified | 4.30828e-06 | 0.999996    | DM1_like  |            6.02518 |
| TCGA-FE-A230-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | low         | no            | no                | nan                |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 31.5537 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.726078  |                       0.726078  |                        nan | mutation_verified | 0.00267896  | 0.997321    | DM1_like  |            6.89569 |
| TCGA-EM-A3FO-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | low         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Male   | 40.6954 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.896989  |                       0.896989  |                        nan | mutation_verified | 6.77084e-06 | 0.999993    | DM1_like  |            6.20445 |
| TCGA-DJ-A3VB-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | low         | no            | no                | Stage III          |              nan |                 nan | Papillary carcinoma, columnar cell      | high               | Male   | 52.6982 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -1.1937    |                       1.1937    |                        nan | mutation_verified | 4.41358e-07 | 1           | DM1_like  |            6.19359 |
| TCGA-BJ-A28X-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | low         | no            | no                | Stage I            |              nan |                 nan | Papillary carcinoma, columnar cell      | high               | Female | 32.1287 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.630202  |                       0.630202  |                        nan | mutation_verified | 9.39099e-06 | 0.999991    | DM1_like  |            7.00976 |
| TCGA-DJ-A2QB-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | RAS             | RAS_like            | high        | no            | no                | Stage IVA          |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Female | 53.6646 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   0.973591  |                      -0.973591  |                        nan | mutation_verified | 0.991557    | 0.00844256  | DM2_like  |            8.40678 |
| TCGA-DO-A1K0-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 30.7077 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.513224  |                       0.513224  |                        nan | mutation_verified | 0.00576734  | 0.994233    | DM1_like  |            7.00289 |
| TCGA-EL-A4KD-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | high        | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Male   | 41.8645 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   0.923545  |                      -0.923545  |                        nan | histology_proxy   | 0.000390315 | 0.99961     | DM1_like  |            8.51266 |
| TCGA-BJ-A192-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | unknown             | unknown         | unknown             | high        | no            | no                | Stage III          |              nan |                 nan | Oxyphilic adenocarcinoma                | medium             | Female | 54.423  | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   0.765738  |                      -0.765738  |                        nan | ambiguous         | 0.938247    | 0.0617532   | DM2_like  |            8.52521 |
| TCGA-L6-A4EU-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | low         | no            | no                | Stage III          |              nan |                 nan | Papillary carcinoma, columnar cell      | high               | Female | 58.0205 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -1.05144   |                       1.05144   |                        nan | mutation_verified | 1.14056e-06 | 0.999999    | DM1_like  |            6.42539 |
| TCGA-DJ-A3V0-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | low         | no            | no                | Stage IVA          |              nan |                 nan | Papillary carcinoma, columnar cell      | high               | Male   | 56.115  | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.73031   |                       0.73031   |                        nan | mutation_verified | 3.80476e-05 | 0.999962    | DM1_like  |            6.71761 |
| TCGA-EM-A2P2-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | BRAF            | BRAF_like           | mid         | no            | no                | Stage III          |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Male   | 46.256  | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.116025  |                       0.116025  |                        nan | mutation_verified | 0.0416284   | 0.958372    | DM1_like  |            7.23888 |
| TCGA-BJ-A45G-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | RAS             | RAS_like            | high        | no            | no                | Stage III          |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Female | 48.8186 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   1.52976   |                      -1.52976   |                        nan | mutation_verified | 0.898523    | 0.101477    | DM2_like  |            9.27759 |
| TCGA-FY-A3R7-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | high        | no            | no                | Stage III          |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 50.3847 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   0.993874  |                      -0.993874  |                        nan | mutation_verified | 0.0145066   | 0.985493    | DM1_like  |            8.63756 |
| TCGA-DJ-A3VG-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | RAS             | RAS_like            | high        | no            | no                | Stage I            |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Male   | 40.9856 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   1.25342   |                      -1.25342   |                        nan | mutation_verified | 0.711085    | 0.288915    | DM2_like  |            8.99578 |
| TCGA-ET-A25O-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 35.384  | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.516503  |                       0.516503  |                        nan | mutation_verified | 0.0643165   | 0.935683    | DM1_like  |            6.67671 |
| TCGA-DJ-A3UR-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | BRAF            | BRAF_like           | mid         | no            | no                | Stage III          |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Female | 51.4935 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.244623  |                       0.244623  |                        nan | mutation_verified | 0.00610315  | 0.993897    | DM1_like  |            7.0817  |
| TCGA-E8-A415-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 39.1184 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   0.534554  |                      -0.534554  |                        nan | mutation_verified | 0.0416558   | 0.958344    | DM1_like  |            8.25852 |
| TCGA-DJ-A4UR-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | RAS             | RAS_like            | mid         | no            | no                | Stage II           |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Female | 36.8597 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   0.621185  |                      -0.621185  |                        nan | mutation_verified | 0.0518471   | 0.948153    | DM1_like  |            8.45109 |
| TCGA-DJ-A13L-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | unknown             | BRAF            | BRAF_like           | low         | no            | no                | nan                |              nan |                 nan | Adenocarcinoma, NOS                     | medium             | Male   | 85.6482 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -1.13349   |                       1.13349   |                        nan | mutation_verified | 6.79447e-07 | 0.999999    | DM1_like  |            6.5678  |
| TCGA-BJ-A45I-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | low         | no            | no                | Stage III          |              nan |                 nan | Papillary carcinoma, columnar cell      | high               | Female | 51.3812 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -1.54428   |                       1.54428   |                        nan | mutation_verified | 3.23445e-07 | 1           | DM1_like  |            5.95478 |
| TCGA-ET-A40S-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | mid         | no            | no                | Stage II           |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Male   | 62.2396 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.0829078 |                       0.0829078 |                        nan | histology_proxy   | 7.4676e-05  | 0.999925    | DM1_like  |            7.3831  |
| TCGA-EL-A3CZ-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | RAS             | RAS_like            | high        | no            | no                | nan                |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 41.462  | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   1.68504   |                      -1.68504   |                        nan | mutation_verified | 0.966366    | 0.0336345   | DM2_like  |            9.45622 |
| TCGA-EM-A2P3-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 47.8987 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.521843  |                       0.521843  |                        nan | mutation_verified | 0.000471716 | 0.999528    | DM1_like  |            7.21193 |
| TCGA-BJ-A2N7-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | RAS             | RAS_like            | high        | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 30.6311 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   0.893842  |                      -0.893842  |                        nan | mutation_verified | 0.839354    | 0.160646    | DM2_like  |            8.56438 |
| TCGA-EL-A3GW-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | RAS             | RAS_like            | high        | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 37.4264 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   0.65131   |                      -0.65131   |                        nan | mutation_verified | 0.790059    | 0.209941    | DM2_like  |            7.99915 |
| TCGA-DJ-A13P-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | mid         | no            | no                | Stage II           |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 52.8323 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.124384  |                       0.124384  |                        nan | mutation_verified | 0.00413669  | 0.995863    | DM1_like  |            7.50493 |
| TCGA-ET-A39O-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | low         | no            | no                | Stage I            |              nan |                 nan | Papillary carcinoma, columnar cell      | high               | Male   | 39.1376 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -1.14894   |                       1.14894   |                        nan | mutation_verified | 1.33841e-05 | 0.999987    | DM1_like  |            6.41195 |
| TCGA-IM-A3ED-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 58.256  | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.0741486 |                       0.0741486 |                        nan | mutation_verified | 0.00103178  | 0.998968    | DM1_like  |            7.63893 |
| TCGA-ET-A39L-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 20.9062 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   0.649935  |                      -0.649935  |                        nan | histology_proxy   | 0.0171387   | 0.982861    | DM1_like  |            8.3034  |
| TCGA-FY-A3YR-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | low         | no            | no                | Stage III          |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 64.052  | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -1.24124   |                       1.24124   |                        nan | mutation_verified | 1.10242e-05 | 0.999989    | DM1_like  |            5.97503 |
| TCGA-KS-A4IB-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | low         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 41.4976 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -1.02635   |                       1.02635   |                        nan | histology_proxy   | 0.012991    | 0.987009    | DM1_like  |            6.67275 |
| TCGA-EM-A3ST-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | RAS             | RAS_like            | mid         | no            | no                | Stage III          |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Female | 62.9432 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.368785  |                       0.368785  |                        nan | mutation_verified | 0.000254689 | 0.999745    | DM1_like  |            6.70307 |
| TCGA-E8-A413-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 38.0808 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   0.0300379 |                      -0.0300379 |                        nan | mutation_verified | 0.000496309 | 0.999504    | DM1_like  |            7.68533 |
| TCGA-EL-A3TA-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Male   | 42.0698 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   0.56188   |                      -0.56188   |                        nan | mutation_verified | 0.000134255 | 0.999866    | DM1_like  |            8.1372  |
| TCGA-DJ-A2PW-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | low         | no            | no                | Stage III          |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Male   | 65.6728 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.843437  |                       0.843437  |                        nan | mutation_verified | 7.01573e-06 | 0.999993    | DM1_like  |            6.1883  |
| TCGA-DJ-A2PO-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Male   | 54.9158 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.401619  |                       0.401619  |                        nan | histology_proxy   | 0.00170659  | 0.998293    | DM1_like  |            7.24729 |
| TCGA-ET-A3BT-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | mid         | no            | no                | nan                |              nan |                 nan | Papillary carcinoma, columnar cell      | high               | Female | 61.2594 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.43588   |                       0.43588   |                        nan | mutation_verified | 0.000209988 | 0.99979     | DM1_like  |            6.78953 |
| TCGA-EM-A2CJ-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | unknown         | RAS_like            | high        | no            | no                | Stage I            |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Female | 49.7413 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   1.41073   |                      -1.41073   |                        nan | histology_proxy   | 0.993745    | 0.00625466  | DM2_like  |            9.5395  |
| TCGA-DJ-A3VD-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 32.9062 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.29762   |                       0.29762   |                        nan | mutation_verified | 3.37192e-06 | 0.999997    | DM1_like  |            7.39856 |
| TCGA-DJ-A1QE-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | low         | no            | no                | Stage III          |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 61.9986 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.735947  |                       0.735947  |                        nan | mutation_verified | 0.00047234  | 0.999528    | DM1_like  |            6.81355 |
| TCGA-BJ-A0ZB-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | low         | no            | no                | Stage IVA          |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Male   | 66.8446 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.593441  |                       0.593441  |                        nan | mutation_verified | 0.00360169  | 0.996398    | DM1_like  |            6.67423 |
| TCGA-IM-A3U3-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | low         | no            | no                | Stage I            |              nan |                 nan | Papillary carcinoma, columnar cell      | high               | Female | 55.0773 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.792471  |                       0.792471  |                        nan | mutation_verified | 0.000454713 | 0.999545    | DM1_like  |            6.50312 |
| TCGA-BJ-A45F-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | RAS             | RAS_like            | high        | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 59.447  | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   0.741333  |                      -0.741333  |                        nan | mutation_verified | 0.965246    | 0.0347537   | DM2_like  |            8.23863 |
| TCGA-EL-A4KG-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | low         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 35.5318 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -1.43156   |                       1.43156   |                        nan | mutation_verified | 1.20151e-06 | 0.999999    | DM1_like  |            5.8709  |

...

| sample_id        | dataset   | platform                         | modality     | tissue_type   | normal_vs_tumor   | histology_subtype   | driver_anchor   | molecular_subtype   | tds_group   | dediff_flag   | aggressive_flag   | ajcc_stage_group   |   ata_risk_proxy |   bethesda_category | clinical_subtype_tag                    | label_confidence   | sex    |      age | outcome_available   | survival_available   | recurrence_available   | source_note                                  |   tds_score |   dedifferentiation_proxy_score |   brs_like_surrogate_score | v3_anchor         |    prob_dm2 |    prob_dm1 | dm_like   |   rai_score_recalc |
|:-----------------|:----------|:---------------------------------|:-------------|:--------------|:------------------|:--------------------|:----------------|:--------------------|:------------|:--------------|:------------------|:-------------------|-----------------:|--------------------:|:----------------------------------------|:-------------------|:-------|---------:|:--------------------|:---------------------|:-----------------------|:---------------------------------------------|------------:|--------------------------------:|---------------------------:|:------------------|------------:|------------:|:----------|-------------------:|
| TCGA-EL-A3GQ-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | unknown             | BRAF            | BRAF_like           | high        | no            | no                | nan                |              nan |                 nan | Basal cell carcinoma, NOS               | medium             | Female | nan      | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   0.733398  |                      -0.733398  |                        nan | mutation_verified | 0.00526831  | 0.994732    | DM1_like  |            8.44395 |
| TCGA-J8-A3YD-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | high        | no            | no                | Stage III          |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female |  47.6824 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   1.65061   |                      -1.65061   |                        nan | mutation_verified | 0.714097    | 0.285903    | DM2_like  |            9.62603 |
| TCGA-EM-A4FQ-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | BRAF            | BRAF_like           | low         | no            | no                | Stage I            |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Female |  42.5298 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.948848  |                       0.948848  |                        nan | mutation_verified | 4.93572e-05 | 0.999951    | DM1_like  |            6.64931 |
| TCGA-FY-A3WA-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | unknown         | RAS_like            | high        | no            | no                | Stage I            |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Female |  52      | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   1.71443   |                      -1.71443   |                        nan | histology_proxy   | 0.999965    | 3.52391e-05 | DM2_like  |            9.53063 |
| TCGA-ET-A2MX-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | low         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Male   |  27.499  | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -1.15993   |                       1.15993   |                        nan | histology_proxy   | 0.00273062  | 0.997269    | DM1_like  |            6.38221 |
| TCGA-FY-A3I5-01B | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | unknown         | RAS_like            | high        | no            | no                | Stage III          |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Female |  64.7611 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   1.47654   |                      -1.47654   |                        nan | histology_proxy   | 0.983477    | 0.0165228   | DM2_like  |            9.42581 |
| TCGA-IM-A41Y-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | unknown             | BRAF            | BRAF_like           | low         | no            | no                | Stage I            |              nan |                 nan | Nonencapsulated sclerosing carcinoma    | medium             | Female |  42.9678 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.885173  |                       0.885173  |                        nan | mutation_verified | 0.00222072  | 0.997779    | DM1_like  |            6.55222 |
| TCGA-EM-A2OV-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | unknown         | RAS_like            | high        | no            | no                | Stage II           |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Female |  64.7666 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   1.34848   |                      -1.34848   |                        nan | histology_proxy   | 0.999923    | 7.68104e-05 | DM2_like  |            9.16685 |
| TCGA-E3-A3E3-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | mid         | no            | no                | nan                |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female |  50.883  | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.458454  |                       0.458454  |                        nan | mutation_verified | 0.0015896   | 0.99841     | DM1_like  |            7.07539 |
| TCGA-ET-A3BQ-01B | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female |  27.4716 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.106844  |                       0.106844  |                        nan | mutation_verified | 0.000206642 | 0.999793    | DM1_like  |            7.68235 |
| TCGA-J8-A3O0-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | RAS             | RAS_like            | high        | no            | no                | Stage I            |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Male   |  38.9158 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   0.719564  |                      -0.719564  |                        nan | mutation_verified | 0.11268     | 0.88732     | DM1_like  |            8.27157 |
| TCGA-EL-A3CL-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | mid         | no            | no                | Stage II           |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female |  70.5298 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.0104837 |                       0.0104837 |                        nan | mutation_verified | 0.000622655 | 0.999377    | DM1_like  |            7.12519 |
| TCGA-EL-A3T3-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | low         | no            | no                | Stage II           |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Male   |  63.0883 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.937253  |                       0.937253  |                        nan | mutation_verified | 1.44894e-05 | 0.999986    | DM1_like  |            6.58186 |
| TCGA-ET-A25J-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female |  40.6379 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.34043   |                       0.34043   |                        nan | mutation_verified | 0.010486    | 0.989514    | DM1_like  |            7.42785 |
| TCGA-EL-A3ZT-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Male   |  35.8741 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.486777  |                       0.486777  |                        nan | mutation_verified | 0.00011179  | 0.999888    | DM1_like  |            7.14244 |
| TCGA-EM-A3AI-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | unknown         | RAS_like            | high        | no            | no                | Stage II           |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Female |  67.9535 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   1.68016   |                      -1.68016   |                        nan | histology_proxy   | 0.999979    | 2.12176e-05 | DM2_like  |            9.55198 |
| TCGA-EM-A22O-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | low         | no            | no                | Stage IVA          |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Male   |  75.9973 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.865168  |                       0.865168  |                        nan | mutation_verified | 0.0245477   | 0.975452    | DM1_like  |            6.30636 |
| TCGA-ET-A40Q-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | low         | no            | no                | nan                |              nan |                 nan | Papillary carcinoma, columnar cell      | high               | Male   |  39.2991 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -1.35605   |                       1.35605   |                        nan | mutation_verified | 4.01111e-06 | 0.999996    | DM1_like  |            6.08607 |
| TCGA-FY-A4B0-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | unknown         | RAS_like            | high        | no            | no                | Stage I            |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Male   |  76.512  | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   1.22771   |                      -1.22771   |                        nan | histology_proxy   | 0.989291    | 0.0107086   | DM2_like  |            9.14519 |
| TCGA-ET-A39K-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | low         | no            | no                | Stage III          |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female |  48.7036 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.961011  |                       0.961011  |                        nan | mutation_verified | 1.51098e-07 | 1           | DM1_like  |            6.43268 |
| TCGA-EM-A2CQ-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | RAS             | RAS_like            | high        | no            | no                | Stage I            |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Female |  40.4545 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   1.20616   |                      -1.20616   |                        nan | mutation_verified | 0.994166    | 0.00583413  | DM2_like  |            8.8546  |
| TCGA-EM-A2CS-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | high        | no            | no                | Stage IVA          |              nan |                 nan | Papillary carcinoma, oxyphilic cell     | high               | Female |  51.9918 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   0.687231  |                      -0.687231  |                        nan | mutation_verified | 0.00823889  | 0.991761    | DM1_like  |            8.22848 |
| TCGA-EM-A2CS-06A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Metastatic    | tumor             | cPTC                | BRAF            | BRAF_like           | mid         | no            | no                | Stage IVA          |              nan |                 nan | Papillary carcinoma, oxyphilic cell     | high               | Female |  51.9918 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   0.422968  |                      -0.422968  |                        nan | mutation_verified | 0.000216877 | 0.999783    | DM1_like  |            7.94268 |
| TCGA-FY-A40M-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | low         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female |  51.258  | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.752259  |                       0.752259  |                        nan | histology_proxy   | 0.117857    | 0.882143    | DM1_like  |            6.21784 |
| TCGA-ET-A4KN-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | low         | no            | no                | Stage III          |              nan |                 nan | Papillary carcinoma, columnar cell      | high               | Female |  51.7153 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -2.84402   |                       2.84402   |                        nan | histology_proxy   | 1.21224e-08 | 1           | DM1_like  |            4.55186 |
| TCGA-KS-A41I-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | RAS             | RAS_like            | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Female |  47.0801 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.412044  |                       0.412044  |                        nan | mutation_verified | 0.00115629  | 0.998844    | DM1_like  |            7.11252 |
| TCGA-DJ-A4UL-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | low         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female |  68.2793 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.85556   |                       0.85556   |                        nan | mutation_verified | 7.05423e-05 | 0.999929    | DM1_like  |            6.6537  |
| TCGA-ET-A40P-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | unknown         | RAS_like            | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Female |  32.4244 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   0.520114  |                      -0.520114  |                        nan | histology_proxy   | 0.000245517 | 0.999754    | DM1_like  |            8.12036 |
| TCGA-DJ-A3V5-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | mid         | no            | no                | Stage III          |              nan |                 nan | Papillary carcinoma, columnar cell      | high               | Female |  76.4353 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.509806  |                       0.509806  |                        nan | mutation_verified | 0.000177172 | 0.999823    | DM1_like  |            7.06297 |
| TCGA-FE-A22Z-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | unknown         | RAS_like            | low         | no            | no                | nan                |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Female |  66.6886 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -1.46452   |                       1.46452   |                        nan | histology_proxy   | 2.50798e-05 | 0.999975    | DM1_like  |            6.00023 |
| TCGA-FE-A237-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | mid         | no            | no                | nan                |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female |  21.629  | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.524478  |                       0.524478  |                        nan | mutation_verified | 0.000112758 | 0.999887    | DM1_like  |            6.88321 |
| TCGA-FK-A3SG-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female |  21.1088 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.145206  |                       0.145206  |                        nan | histology_proxy   | 3.16845e-05 | 0.999968    | DM1_like  |            7.38313 |
| TCGA-GE-A2C6-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | low         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female |  33.0404 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.634083  |                       0.634083  |                        nan | mutation_verified | 2.8026e-06  | 0.999997    | DM1_like  |            6.70911 |
| TCGA-EM-A22Q-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | unknown         | RAS_like            | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Male   |  36.5257 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   0.568189  |                      -0.568189  |                        nan | histology_proxy   | 0.140435    | 0.859565    | DM1_like  |            8.13347 |
| TCGA-FY-A3W9-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | RAS             | RAS_like            | high        | no            | no                | Stage II           |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Female |  66.0808 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   1.42703   |                      -1.42703   |                        nan | mutation_verified | 0.978507    | 0.0214925   | DM2_like  |            9.4773  |
| TCGA-DE-A4MC-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | low         | no            | no                | nan                |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female |  44.7036 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -1.16468   |                       1.16468   |                        nan | histology_proxy   | 2.35653e-06 | 0.999998    | DM1_like  |            6.31178 |
| TCGA-EL-A4K9-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | mid         | no            | no                | Stage II           |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Male   |  68.3915 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.276513  |                       0.276513  |                        nan | mutation_verified | 0.000366013 | 0.999634    | DM1_like  |            7.06056 |
| TCGA-FK-A3SB-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female |  28.1807 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.235783  |                       0.235783  |                        nan | mutation_verified | 0.00120187  | 0.998798    | DM1_like  |            7.38088 |
| TCGA-H2-A422-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female |  40.3039 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.139735  |                       0.139735  |                        nan | mutation_verified | 0.118188    | 0.881812    | DM1_like  |            7.62761 |
| TCGA-EM-A3FM-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | low         | no            | no                | Stage IVA          |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Male   |  57.3251 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -1.46987   |                       1.46987   |                        nan | mutation_verified | 5.16788e-06 | 0.999995    | DM1_like  |            5.96717 |
| TCGA-DJ-A13S-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | unknown         | RAS_like            | high        | no            | no                | Stage I            |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Female |  19.436  | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   0.842244  |                      -0.842244  |                        nan | histology_proxy   | 0.99358     | 0.0064197   | DM2_like  |            8.40564 |
| TCGA-EM-A3SY-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | unknown         | RAS_like            | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Male   |  36.8542 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.194991  |                       0.194991  |                        nan | histology_proxy   | 0.998978    | 0.00102188  | DM2_like  |            7.30696 |
| TCGA-EM-A2P1-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | low         | no            | no                | nan                |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Male   |  34.7242 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -1.02382   |                       1.02382   |                        nan | mutation_verified | 1.18107e-06 | 0.999999    | DM1_like  |            6.73862 |
| TCGA-EL-A3T1-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | low         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female |  38.3984 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -1.5159    |                       1.5159    |                        nan | mutation_verified | 3.75446e-05 | 0.999962    | DM1_like  |            5.72277 |
| TCGA-DJ-A2PS-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female |  40.8378 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   0.331201  |                      -0.331201  |                        nan | mutation_verified | 0.00215273  | 0.997847    | DM1_like  |            8.0854  |
| TCGA-EL-A3TB-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | mid         | no            | no                | Stage III          |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female |  47.0773 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   0.486016  |                      -0.486016  |                        nan | histology_proxy   | 0.00365558  | 0.996344    | DM1_like  |            7.87649 |
| TCGA-BJ-A45C-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | high        | no            | no                | Stage III          |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Male   |  78.0014 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   1.00152   |                      -1.00152   |                        nan | histology_proxy   | 0.995313    | 0.0046866   | DM2_like  |            8.82944 |
| TCGA-DJ-A3UX-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | BRAF            | BRAF_like           | high        | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female |  46.883  | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   0.892034  |                      -0.892034  |                        nan | mutation_verified | 0.0113022   | 0.988698    | DM1_like  |            8.71396 |
| TCGA-DJ-A2PX-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | unknown         | RAS_like            | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Female |  54.9979 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   0.160796  |                      -0.160796  |                        nan | histology_proxy   | 0.00062748  | 0.999373    | DM1_like  |            7.59361 |
| TCGA-EL-A3ZS-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female |  22.1109 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   0.371598  |                      -0.371598  |                        nan | histology_proxy   | 2.00388e-05 | 0.99998     | DM1_like  |            8.0443  |

### 5.2 Driver별 DM score 분포 (numerical)
| driver_anchor   |   n |      mean |      median |       std |         q25 |        q75 |
|:----------------|----:|----------:|------------:|----------:|------------:|-----------:|
| BRAF            | 281 | 0.0093832 | 0.000119084 | 0.0578957 | 7.81818e-06 | 0.00112512 |
| RAS             |  54 | 0.695648  | 0.876811    | 0.360537  | 0.458747    | 0.983369   |
| unknown         | 176 | 0.386463  | 0.0165871   | 0.472326  | 0.000237726 | 0.991994   |
| driver_anchor   | dm_like   |   n |
|:----------------|:----------|----:|
| BRAF            | DM1_like  | 279 |
| BRAF            | DM2_like  |   2 |
| NTRK            | DM1_like  |   1 |
| RAS             | DM1_like  |  15 |
| RAS             | DM2_like  |  39 |
| TP53            | DM2_like  |   1 |
| unknown         | DM1_like  | 108 |
| unknown         | DM2_like  |  68 |
### 5.3 4-group 정의 후 임상
- file: `/home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/A2_4group_clinical.tsv`
| driver_dm_group   |   n |   age_mean |   tds_mean |   rai_mean |
|:------------------|----:|-----------:|-----------:|-----------:|
| BRAF/DM1_like     | 279 |    48.3201 |  -0.609535 |    6.90339 |
| BRAF/DM2_like     |   2 |    44.8036 |   1.5307   |    9.46536 |
| RAS/DM1_like      |  15 |    45.243  |   0.590267 |    8.14408 |
| RAS/DM2_like      |  39 |    46.3422 |   1.11078  |    8.76678 |

MISSING: pairwise comparison p-value table under project/results/v17p3/tables
### 5.4 DM score vs driver 독립성
|             |   prob_dm1 |   prob_dm2 |   BRAF_binary |   RAS_binary |
|:------------|-----------:|-----------:|--------------:|-------------:|
| prob_dm1    |   1        |  -1        |      0.577919 |    -0.427334 |
| prob_dm2    |  -1        |   1        |     -0.577919 |     0.427334 |
| BRAF_binary |   0.577919 |  -0.577919 |      1        |    -0.377485 |
| RAS_binary  |  -0.427334 |   0.427334 |     -0.377485 |     1        |
MISSING: partial correlation controlling for age, sex not found as standalone output file
## SEC 6. A3 Drug response — RAW

### 6.1 Cell line DM score
- file: `/home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/A3_celline_dm_scores.tsv`
| sample            |   DM1_score |   DM2_score | dm_like   | v14_prediction   | mutation_label   |
|:------------------|------------:|------------:|:----------|:-----------------|:-----------------|
| 8305C_THYROID     |    16.7269  |   0.0435833 | DM1_like  | BRAF_like        | BRAF             |
| 8505C_THYROID     |    24.7592  |   0.02075   | DM1_like  | BRAF_like        | BRAF             |
| BCPAP_THYROID     |    15.2622  |   0.0495    | DM1_like  | BRAF_like        | BRAF             |
| BHT101_THYROID    |    30.0718  |   0.02965   | DM1_like  | RAS_like         | BRAF             |
| CAL62_THYROID     |    19.8638  |   0.147667  | DM1_like  | RAS_like         | RAS              |
| FTC133_THYROID    |    20.2357  |   0.450867  | DM1_like  | RAS_like         | other            |
| FTC238_THYROID    |    12.9766  |   0.08185   | DM1_like  | BRAF_like        | other            |
| MB1_THYROID       |    12.539   |   0.141733  | DM1_like  | BRAF_like        | other            |
| ML1_THYROID       |    21.6877  |  28.5213    | DM2_like  | RAS_like         | other            |
| S117_SOFT_TISSUE  |    11.6939  |   0.04765   | DM1_like  | BRAF_like        | other            |
| SW579_THYROID     |     6.22314 |   0.0576    | DM1_like  | RAS_like         | other            |
| TT2609C02_THYROID |    18.7048  |   0.0411    | DM1_like  | BRAF_like        | RAS              |
| TT_THYROID        |    11.9353  |   6.08668   | DM1_like  | RAS_like         | other            |

### 6.2 DM1-selective drug top 20
- file: `/home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/A3_top_drugs_dm1_selective.tsv`
| compound   | dm1_mean_lfc   | dm2_mean_lfc   | delta_dm1_minus_dm2   | pvalue   | fdr   |
|------------|----------------|----------------|-----------------------|----------|-------|

### 6.3 DM2-selective drug top 20
- file: `/home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/A3_top_drugs_dm2_selective.tsv`
| compound   | dm1_mean_lfc   | dm2_mean_lfc   | delta_dm1_minus_dm2   | pvalue   | fdr   |
|------------|----------------|----------------|-----------------------|----------|-------|

### 6.4 Pathway-drug coherence
- file: `/home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/A3_pathway_drug_coherence.tsv`
| compound   | dm1_mean_lfc   | dm2_mean_lfc   | delta_dm1_minus_dm2   | pvalue   | fdr   |
|------------|----------------|----------------|-----------------------|----------|-------|

### 6.5 Tumor predicted response (PERCEPTION-style)
- file: `/home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/A3_tumor_predicted_response.tsv`
MISSING: empty file `/home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/A3_tumor_predicted_response.tsv`
## SEC 7. A4 Pan-cancer — RAW

### 7.1 적용된 cancer
- file: `/home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/A4_pancancer_dm_signature_transfer.tsv`
| cancer   |   n |   dm1_overlap |   dm2_overlap |   dm_score_corr |   cluster_delta |
|:---------|----:|--------------:|--------------:|----------------:|----------------:|
| THCA     | 392 |            18 |            18 |               1 |         2.7727  |
| LUAD     | 236 |            16 |            16 |               1 |         3.33785 |
| COAD     | 372 |            14 |            14 |               1 |         1.86393 |
| LGG      | 142 |            14 |            14 |               1 |         1.07802 |
| SKCM     | 212 |             6 |             6 |               1 |         1.69904 |

### 7.2 Conserved genes
- file: `/home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/A4_conserved_genes_across_cancers.tsv`
| gene    |   cancer_count |
|:--------|---------------:|
| CDKN2A  |              5 |
| FOSL1   |              5 |
| ETV4    |              5 |
| DUSP6   |              5 |
| THRA    |              4 |
| DUSP4   |              4 |
| DIO2    |              4 |
| DUSP5   |              4 |
| ETV5    |              4 |
| FOXP3   |              4 |
| HLA-DRA |              4 |
| MMP9    |              4 |
| KLK10   |              3 |
| MET     |              3 |
| DIO1    |              2 |
| IYD     |              2 |
| TPO     |              2 |
| SLC26A4 |              2 |
| FOXE1   |              1 |
| SLC5A8  |              1 |

## SEC 8. A5 Genomic + immune — RAW

### 8.1 Genomic instability
- file: `/home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/A5_cnv_aneuploidy_per_cluster.tsv`
| metric     | status               |
|:-----------|:---------------------|
| CNV_burden | missing_local_source |
| aneuploidy | missing_local_source |

### 8.2 TMB / MSI
- file: `/home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/A5_tmb_msi_per_cluster.tsv`
| sample_id        | dataset   | platform                         | modality     | tissue_type   | normal_vs_tumor   | histology_subtype   | driver_anchor   | molecular_subtype   | tds_group   | dediff_flag   | aggressive_flag   | ajcc_stage_group   |   ata_risk_proxy |   bethesda_category | clinical_subtype_tag                    | label_confidence   | sex    |     age | outcome_available   | survival_available   | recurrence_available   | source_note                                  |   tds_score |   dedifferentiation_proxy_score |   brs_like_surrogate_score | v3_anchor       | v17_dark_cluster   |   tds16_score_v17 |   brs71_score_v17 |   rai_score_v17 | tcga12       |   total_mutations |   unique_driver_strings |   cytolytic_score |    CD274 |     PDCD1 |     CTLA4 |     IDO1 |    HLA-A |    HLA-B |    HLA-C |      B2M |
|:-----------------|:----------|:---------------------------------|:-------------|:--------------|:------------------|:--------------------|:----------------|:--------------------|:------------|:--------------|:------------------|:-------------------|-----------------:|--------------------:|:----------------------------------------|:-------------------|:-------|--------:|:--------------------|:---------------------|:-----------------------|:---------------------------------------------|------------:|--------------------------------:|---------------------------:|:----------------|:-------------------|------------------:|------------------:|----------------:|:-------------|------------------:|------------------------:|------------------:|---------:|----------:|----------:|---------:|---------:|---------:|---------:|---------:|
| TCGA-FK-A3SE-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 31.4086 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   0.0449715 |                      -0.0449715 |                        nan | histology_proxy | DM1                |           6.96256 |               nan |         8.02363 | TCGA-FK-A3SE |                 0 |                       0 |          0.838522 | 1.89636  | 0.557719  | 0.420972  | 0.557719 |  8.94436 |  8.786   |  8.77273 | 10.7395  |
| TCGA-FY-A2QD-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | high        | no            | no                | nan                |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 61.577  | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   1.70812   |                      -1.70812   |                        nan | histology_proxy | DM2                |           8.62571 |               nan |         9.6369  | TCGA-FY-A2QD |                 0 |                       0 |          0.553785 | 1.5418   | 0.365738  | 0.0257854 | 2.28607  |  8.6665  |  7.44263 |  7.0048  | 10.3523  |
| TCGA-EL-A4KD-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | high        | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Male   | 41.8645 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   0.923545  |                      -0.923545  |                        nan | histology_proxy | DM1                |           7.84114 |               nan |         8.95363 | TCGA-EL-A4KD |                 0 |                       0 |          2.63608  | 3.02656  | 1.3217    | 1.13347   | 1.56118  |  9.54072 |  9.91064 |  9.70173 | 10.8807  |
| TCGA-BJ-A192-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | unknown             | unknown         | unknown             | high        | no            | no                | Stage III          |              nan |                 nan | Oxyphilic adenocarcinoma                | medium             | Female | 54.423  | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   0.765738  |                      -0.765738  |                        nan | ambiguous       | DM2                |           7.68333 |               nan |         9.01191 | TCGA-BJ-A192 |                 0 |                       0 |          1.32534  | 2.82566  | 0.0240007 | 0.0476087 | 1.32755  |  9.31237 |  9.16514 |  8.84897 |  9.5853  |
| TCGA-ET-A40S-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | mid         | no            | no                | Stage II           |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Male   | 62.2396 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.0829078 |                       0.0829078 |                        nan | histology_proxy | DM1                |           6.83468 |               nan |         7.48305 | TCGA-ET-A40S |                 0 |                       0 |          3.42978  | 4.02443  | 2.61006   | 2.39935   | 2.31553  | 10.789   | 11.8027  | 11.457   | 12.5629  |
| TCGA-ET-A39L-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 20.9062 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   0.649935  |                      -0.649935  |                        nan | histology_proxy | DM1                |           7.56753 |               nan |         8.43548 | TCGA-ET-A39L |                 0 |                       0 |          2.37656  | 1.60287  | 0.966218  | 0.100643  | 0.644832 |  9.13171 |  9.20913 |  9.35162 | 10.3672  |
| TCGA-KS-A4IB-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | low         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 41.4976 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -1.02635   |                       1.02635   |                        nan | histology_proxy | DM1                |           5.89124 |               nan |         6.32345 | TCGA-KS-A4IB |                 1 |                       1 |          1.33339  | 3.47426  | 0.683118  | 0.190688  | 0.489273 |  9.864   | 10.721   |  9.7035  | 11.5164  |
| TCGA-DJ-A2PO-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Male   | 54.9158 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.401619  |                       0.401619  |                        nan | histology_proxy | DM1                |           6.51597 |               nan |         7.00901 | TCGA-DJ-A2PO |                 0 |                       0 |          3.72882  | 4.1104   | 2.80782   | 1.43315   | 1.62828  | 11.1348  | 11.8827  | 10.9039  | 12.4521  |
| TCGA-EM-A2CJ-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | unknown         | RAS_like            | high        | no            | no                | Stage I            |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Female | 49.7413 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   1.41073   |                      -1.41073   |                        nan | histology_proxy | DM2                |           8.32832 |               nan |         9.76971 | TCGA-EM-A2CJ |                 0 |                       0 |          0.710587 | 1.58     | 0.763954  | 0.0249628 | 0.864186 |  8.19517 |  8.43376 |  8.84677 |  9.85105 |
| TCGA-DJ-A1QL-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | unknown         | RAS_like            | high        | no            | no                | Stage II           |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Male   | 70.3984 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   1.5134    |                      -1.5134    |                        nan | histology_proxy | DM2                |           8.43099 |               nan |         9.51776 | TCGA-DJ-A1QL |                 1 |                       1 |          2.39803  | 1.36615  | 0.832221  | 0.0242711 | 0.950928 |  9.34672 |  9.89435 |  9.01656 | 10.0544  |
| TCGA-EM-A1CS-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 55.5346 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.0232746 |                       0.0232746 |                        nan | histology_proxy | DM1                |           6.89432 |               nan |         7.36269 | TCGA-EM-A1CS |                 0 |                       0 |          1.53293  | 3.35319  | 0.539891  | 0.136562  | 1.16872  |  9.16148 |  9.85622 |  9.35173 | 10.6516  |
| TCGA-CE-A484-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 37.8234 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   0.202154  |                      -0.202154  |                        nan | histology_proxy | DM1                |           7.11974 |               nan |         7.97868 | TCGA-CE-A484 |                 1 |                       1 |          1.67212  | 1.26336  | 2.44381   | 0.367771  | 0.423952 |  9.4519  |  9.49423 |  9.06358 | 10.3684  |
| TCGA-EL-A3ZP-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Male   | 19.4415 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   0.0380733 |                      -0.0380733 |                        nan | histology_proxy | DM1                |           6.95566 |               nan |         8.0507  | TCGA-EL-A3ZP |                 0 |                       0 |          4.08105  | 3.89232  | 3.58949   | 2.1773    | 4.02435  | 10.8179  | 11.5279  | 10.7554  | 12.177   |
| TCGA-BJ-A0YZ-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | mid         | no            | no                | Stage II           |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Male   | 65.9877 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   0.491545  |                      -0.491545  |                        nan | histology_proxy | DM2                |           7.40914 |               nan |         8.21595 | TCGA-BJ-A0YZ |                 0 |                       0 |          1.31121  | 0.944485 | 0.407402  | 0.0259182 | 4.27118  |  8.98814 |  8.56389 |  9.13352 | 10.3622  |
| TCGA-FY-A3NM-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | high        | no            | no                | Stage III          |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 48.553  | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   2.16569   |                      -2.16569   |                        nan | histology_proxy | DM2                |           9.08328 |               nan |        10.4552  | TCGA-FY-A3NM |                 0 |                       0 |          1.05508  | 2.01097  | 0.175578  | 0.229657  | 1.46358  |  7.71764 |  8.42761 |  8.18974 |  9.83521 |
| TCGA-BJ-A3PR-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 69.9932 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.254757  |                       0.254757  |                        nan | histology_proxy | DM1                |           6.66283 |               nan |         7.84422 | TCGA-BJ-A3PR |                 0 |                       0 |          5.28118  | 4.45473  | 5.55695   | 3.91943   | 5.40976  | 11.3847  | 12.3287  | 10.9473  | 13.2891  |
| TCGA-EL-A3CO-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | high        | no            | no                | Stage III          |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Male   | 88.011  | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   0.914544  |                      -0.914544  |                        nan | histology_proxy | DM1                |           7.83213 |               nan |         9.23233 | TCGA-EL-A3CO |                 1 |                       1 |          3.61151  | 4.87503  | 1.97346   | 1.77862   | 3.98344  | 10.3312  | 10.5465  |  9.64063 | 11.556   |
| TCGA-CE-A482-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 27.9042 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.489247  |                       0.489247  |                        nan | histology_proxy | DM1                |           6.42834 |               nan |         6.77373 | TCGA-CE-A482 |                 0 |                       0 |          1.68446  | 3.42153  | 1.05377   | 0.713111  | 0.893474 | 10.0325  |  9.90853 | 10.0734  | 10.9514  |
| TCGA-EM-A3FQ-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 19.6468 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   0.613279  |                      -0.613279  |                        nan | histology_proxy | DM1                |           7.53087 |               nan |         8.68297 | TCGA-EM-A3FQ |                 0 |                       0 |          3.11261  | 3.6919   | 2.4086    | 2.31915   | 3.51788  | 10.7221  | 11.0455  | 10.7051  | 11.6551  |
| TCGA-EL-A3GO-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | high        | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 31.77   | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   0.667869  |                      -0.667869  |                        nan | histology_proxy | DM2                |           7.58546 |               nan |         8.6053  | TCGA-EL-A3GO |                 1 |                       1 |          2.42594  | 1.6536   | 0.187088  | 0.332988  | 1.02749  |  8.49518 |  9.26351 |  8.54644 | 10.4997  |
| TCGA-DJ-A3UZ-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | high        | no            | no                | Stage III          |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 70.7625 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   0.935377  |                      -0.935377  |                        nan | histology_proxy | DM1                |           7.85297 |               nan |         8.85158 | TCGA-DJ-A3UZ |                 0 |                       0 |          3.14528  | 3.49142  | 2.26853   | 2.02758   | 3.35269  | 10.3334  | 10.8213  |  9.89604 | 12.292   |
| TCGA-E8-A432-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | mid         | no            | no                | nan                |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 59.1211 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.506581  |                       0.506581  |                        nan | histology_proxy | DM1                |           6.41101 |               nan |         6.48077 | TCGA-E8-A432 |                 0 |                       0 |          2.01167  | 3.67658  | 1.82633   | 0.755854  | 2.11472  | 10.086   |  9.9902  |  9.9968  | 11.9653  |
| TCGA-DJ-A4UT-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 32.3422 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.25731   |                       0.25731   |                        nan | histology_proxy | DM1                |           6.66028 |               nan |         7.38698 | TCGA-DJ-A4UT |                 0 |                       0 |          4.1742   | 4.21709  | 3.84097   | 2.95011   | 4.72238  | 11.3085  | 11.4732  | 10.8588  | 12.6204  |
| TCGA-ET-A40R-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | unknown         | RAS_like            | high        | no            | no                | Stage I            |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Female | 24.5394 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   1.33359   |                      -1.33359   |                        nan | histology_proxy | DM1                |           8.25118 |               nan |         9.6007  | TCGA-ET-A40R |                 0 |                       0 |          2.31913  | 2.3605   | 2.31361   | 1.74155   | 2.02168  | 10.2623  | 10.6706  |  9.4312  | 11.7106  |
| TCGA-EM-A3FR-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | high        | no            | no                | Stage III          |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 55.7125 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   0.91294   |                      -0.91294   |                        nan | histology_proxy | DM2                |           7.83053 |               nan |         8.91875 | TCGA-EM-A3FR |                 0 |                       0 |          1.59883  | 2.05461  | 0.45039   | 0.594345  | 1.71986  | 10.4909  |  8.71298 |  8.10718 | 10.891   |
| TCGA-DJ-A4UQ-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | unknown             | unknown         | unknown             | mid         | no            | no                | Stage IVA          |              nan |                 nan | Nonencapsulated sclerosing carcinoma    | medium             | Male   | 60.4983 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   0.293168  |                      -0.293168  |                        nan | ambiguous       | DM1                |           7.21076 |               nan |         7.99091 | TCGA-DJ-A4UQ |                 0 |                       0 |          2.30964  | 2.95887  | 1.51783   | 1.09227   | 1.1836   |  9.7953  |  9.90735 |  9.50239 | 11.113   |
| TCGA-ET-A4KQ-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | unknown         | RAS_like            | high        | no            | no                | Stage II           |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Female | 45.681  | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   1.58759   |                      -1.58759   |                        nan | histology_proxy | DM2                |           8.50518 |               nan |         9.59345 | TCGA-ET-A4KQ |                 1 |                       1 |          1.50554  | 1.58223  | 0.773004  | 0.0677156 | 0.405314 |  8.49104 |  8.8669  |  8.44995 | 10.2597  |
| TCGA-CE-A481-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 41.8015 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   0.488726  |                      -0.488726  |                        nan | histology_proxy | DM1                |           7.40632 |               nan |         8.15195 | TCGA-CE-A481 |                 0 |                       0 |          3.21979  | 3.42083  | 2.75479   | 1.60313   | 3.41412  | 10.2826  | 10.415   |  9.88526 | 11.7685  |
| TCGA-EM-A2CP-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | unknown         | RAS_like            | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Female | 26.7187 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   0.455291  |                      -0.455291  |                        nan | histology_proxy | DM1                |           7.37288 |               nan |         8.6882  | TCGA-EM-A2CP |                 0 |                       0 |          5.05209  | 4.66712  | 4.66712   | 3.03829   | 5.55109  | 11.2894  | 11.8604  | 11.8164  | 12.9118  |
| TCGA-EL-A3H3-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | low         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 19.165  | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.711258  |                       0.711258  |                        nan | histology_proxy | DM1                |           6.20633 |               nan |         6.83075 | TCGA-EL-A3H3 |                 0 |                       0 |          3.28696  | 4.23398  | 2.42324   | 1.58135   | 2.74363  | 11.0029  | 11.848   | 11.002   | 13.0102  |
| TCGA-EL-A3CS-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | low         | no            | no                | Stage IVA          |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 62.9514 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -1.13313   |                       1.13313   |                        nan | histology_proxy | DM1                |           5.78447 |               nan |         6.48032 | TCGA-EL-A3CS |                 0 |                       0 |          0.81823  | 4.01225  | 0.448196  | 0.674625  | 1.01917  | 10.0589  | 10.2652  |  9.99991 | 11.0548  |
| TCGA-FK-A3SD-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | high        | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 61.7659 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   1.62862   |                      -1.62862   |                        nan | histology_proxy | DM2                |           8.54621 |               nan |         9.9413  | TCGA-FK-A3SD |                 0 |                       0 |          0.392661 | 2.03151  | 1.37869   | 0.0201115 | 5.03482  |  9.71723 |  9.07666 |  8.8891  | 10.5986  |
| TCGA-DJ-A13W-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | unknown         | RAS_like            | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Female | 45.3333 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.465718  |                       0.465718  |                        nan | histology_proxy | DM1                |           6.45187 |               nan |         7.31335 | TCGA-DJ-A13W |                 0 |                       0 |          1.42938  | 0.611683 | 0.576884  | 0.285106  | 2.83883  |  8.90942 |  8.42171 |  6.75186 |  7.8071  |
| TCGA-DJ-A13R-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | unknown         | RAS_like            | mid         | no            | no                | Stage III          |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Male   | 50.2505 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.415482  |                       0.415482  |                        nan | histology_proxy | DM1                |           6.50211 |               nan |         7.5934  | TCGA-DJ-A13R |                 0 |                       0 |          5.10047  | 3.24778  | 4.49662   | 1.79485   | 6.04866  | 11.5168  | 11.1044  | 10.3054  | 11.9406  |
| TCGA-EL-A3ZK-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 41.399  | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.397993  |                       0.397993  |                        nan | histology_proxy | DM1                |           6.5196  |               nan |         7.04622 | TCGA-EL-A3ZK |               nan |                     nan |          1.95303  | 5.8787   | 1.12687   | 1.96561   | 2.11133  |  9.96024 | 10.4567  |  9.94812 | 11.8822  |
| TCGA-ET-A3DS-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | high        | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 33.5441 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   2.23116   |                      -2.23116   |                        nan | histology_proxy | DM2                |           9.14875 |               nan |        10.5403  | TCGA-ET-A3DS |                 0 |                       0 |          1.64834  | 2.02254  | 1.37098   | 0.587501  | 1.59754  |  8.3752  |  8.73162 |  9.22236 | 10.3235  |
| TCGA-EL-A3D4-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | mid         | no            | no                | Stage IVA          |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Male   | 62.7598 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   0.617167  |                      -0.617167  |                        nan | histology_proxy | DM1                |           7.53476 |               nan |         8.50464 | TCGA-EL-A3D4 |                 0 |                       0 |          1.89577  | 2.60715  | 1.022     | 0.112155  | 0.822112 | 10.34    |  9.58341 |  9.52904 | 10.9795  |
| TCGA-EM-A3AN-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 36.1068 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   0.272879  |                      -0.272879  |                        nan | histology_proxy | DM1                |           7.19047 |               nan |         7.72235 | TCGA-EM-A3AN |                 0 |                       0 |          1.89258  | 2.80769  | 0.98378   | 0.450607  | 2.2306   |  9.85776 |  9.85459 |  9.98242 | 11.5034  |
| TCGA-DE-A0Y2-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | low         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 30.7379 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.722768  |                       0.722768  |                        nan | histology_proxy | DM1                |           6.19482 |               nan |         6.61181 | TCGA-DE-A0Y2 |                 0 |                       0 |          0.929393 | 3.71768  | 0.640706  | 0.399971  | 0.451439 | 10.1925  | 10.3802  |  9.96885 | 11.741   |
| TCGA-CE-A27D-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | high        | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 28.846  | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   1.16773   |                      -1.16773   |                        nan | histology_proxy | DM1                |           8.08533 |               nan |         9.30818 | TCGA-CE-A27D |                 0 |                       0 |          4.41596  | 4.02803  | 4.00495   | 2.24097   | 6.20814  | 11.3559  | 11.4857  | 10.8933  | 12.7246  |
| TCGA-BJ-A0Z0-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | high        | no            | no                | Stage II           |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Male   | 55.7454 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   1.78021   |                      -1.78021   |                        nan | histology_proxy | DM2                |           8.6978  |               nan |         9.54684 | TCGA-BJ-A0Z0 |                 0 |                       0 |          0.801932 | 1.82566  | 0.264691  | 0         | 0.714356 |  8.77174 |  7.97964 |  7.48515 | 10.1068  |
| TCGA-FK-A4UB-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | mid         | no            | no                | Stage III          |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Male   | 51.4497 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.0620745 |                       0.0620745 |                        nan | histology_proxy | DM1                |           6.85552 |               nan |         7.43138 | TCGA-FK-A4UB |                 0 |                       0 |          1.48553  | 2.31998  | 1.46108   | 0.572387  | 1.28005  | 10.1503  | 10.2387  |  8.60565 | 10.5602  |
| TCGA-H2-A421-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | low         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 34.653  | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.685849  |                       0.685849  |                        nan | histology_proxy | DM1                |           6.23174 |               nan |         7.09984 | TCGA-H2-A421 |               nan |                     nan |          3.35964  | 4.1397   | 3.3995    | 2.67348   | 2.49112  | 11.2538  | 11.9274  | 11.0112  | 12.7889  |
| TCGA-CE-A483-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | low         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 34.9569 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -0.553398  |                       0.553398  |                        nan | histology_proxy | DM1                |           6.36419 |               nan |         6.72499 | TCGA-CE-A483 |               nan |                     nan |          1.55166  | 3.48745  | 0.367617  | 0.334795  | 0.399709 |  9.70301 | 10.1817  |  9.1897  | 11.8955  |
| TCGA-FY-A76V-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | high        | no            | no                | Stage II           |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Male   | 54.4066 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   0.656027  |                      -0.656027  |                        nan | histology_proxy | DM1                |           7.57362 |               nan |         8.30574 | TCGA-FY-A76V |               nan |                     nan |          1.9471   | 3.24788  | 0.907667  | 0.535896  | 1.24684  |  8.63893 |  8.8831  |  8.67799 | 10.5831  |
| TCGA-DO-A2HM-01B | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | low         | no            | no                | Stage IVA          |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Male   | 49.4675 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -1.11342   |                       1.11342   |                        nan | histology_proxy | DM1                |           5.80417 |               nan |         6.40334 | TCGA-DO-A2HM |                 0 |                       0 |          4.27281  | 5.01979  | 2.91461   | 3.36734   | 5.64062  | 10.5497  | 11.3927  | 10.5274  | 12.4981  |
| TCGA-EL-A3T0-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | high        | no            | no                | Stage III          |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 45.9904 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   0.968623  |                      -0.968623  |                        nan | histology_proxy | DM1                |           7.88621 |               nan |         8.92073 | TCGA-EL-A3T0 |               nan |                     nan |          2.29402  | 3.28772  | 1.55712   | 0.914156  | 1.51788  |  9.693   |  9.99665 |  9.31306 | 11.3485  |
| TCGA-BJ-A0ZJ-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | high        | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Male   | 36      | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   0.678984  |                      -0.678984  |                        nan | histology_proxy | DM1                |           7.59657 |               nan |         8.4326  | TCGA-BJ-A0ZJ |               nan |                     nan |          3.8632   | 3.3286   | 3.77705   | 2.05736   | 3.02592  | 11.0077  | 11.4858  | 10.6053  | 12.422   |
| TCGA-EL-A3T9-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | NTRK            | BRAF_like           | low         | no            | no                | Stage IV           |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female | 69.7714 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  -1.75854   |                       1.75854   |                        nan | histology_proxy | DM1                |           5.15906 |               nan |         5.9317  | TCGA-EL-A3T9 |                 1 |                       1 |          0.81239  | 2.66079  | 0.490295  | 2.27184   | 1.67577  |  8.9158  |  9.1876  |  9.05455 | 10.0959  |
| TCGA-EM-A1YD-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | unknown         | RAS_like            | high        | no            | no                | Stage II           |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Female | 53.0212 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |   1.86334   |                      -1.86334   |                        nan | histology_proxy | DM2                |           8.78093 |               nan |        10.1612  | TCGA-EM-A1YD |               nan |                     nan |          1.96123  | 2.13086  | 0.667018  | 0.526911  | 2.92491  |  9.02557 |  9.37094 |  9.365   | 10.4999  |

...

| sample_id        | dataset   | platform                         | modality     | tissue_type   | normal_vs_tumor   | histology_subtype   | driver_anchor   | molecular_subtype   | tds_group   | dediff_flag   | aggressive_flag   | ajcc_stage_group   |   ata_risk_proxy |   bethesda_category | clinical_subtype_tag                    | label_confidence   | sex    |      age | outcome_available   | survival_available   | recurrence_available   | source_note                                  |    tds_score |   dedifferentiation_proxy_score |   brs_like_surrogate_score | v3_anchor       | v17_dark_cluster   |   tds16_score_v17 |   brs71_score_v17 |   rai_score_v17 | tcga12       |   total_mutations |   unique_driver_strings |   cytolytic_score |    CD274 |     PDCD1 |     CTLA4 |      IDO1 |    HLA-A |    HLA-B |    HLA-C |      B2M |
|:-----------------|:----------|:---------------------------------|:-------------|:--------------|:------------------|:--------------------|:----------------|:--------------------|:------------|:--------------|:------------------|:-------------------|-----------------:|--------------------:|:----------------------------------------|:-------------------|:-------|---------:|:--------------------|:---------------------|:-----------------------|:---------------------------------------------|-------------:|--------------------------------:|---------------------------:|:----------------|:-------------------|------------------:|------------------:|----------------:|:-------------|------------------:|------------------------:|------------------:|---------:|----------:|----------:|----------:|---------:|---------:|---------:|---------:|
| TCGA-ET-A39T-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | unknown         | RAS_like            | mid         | no            | no                | Stage III          |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Female |  55.206  | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs | -0.232089    |                     0.232089    |                        nan | histology_proxy | DM1                |           6.6855  |               nan |         7.12165 | TCGA-ET-A39T |               nan |                     nan |          1.81449  | 4.43952  | 0.548777  | 0.28944   | 0.584841  | 10.3141  | 10.9952  | 10.0005  | 12.2804  |
| TCGA-DE-A69K-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | mid         | no            | no                | Stage III          |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female |  58.5982 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs | -0.0754255   |                     0.0754255   |                        nan | histology_proxy | DM1                |           6.84217 |               nan |         7.14371 | TCGA-DE-A69K |                 0 |                       0 |          1.88734  | 3.03546  | 0.702847  | 0.575701  | 1.05581   | 10.3214  | 10.366   | 10.0911  | 11.432   |
| TCGA-ET-A25M-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Male   |  33.0869 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs | -0.517999    |                     0.517999    |                        nan | histology_proxy | DM1                |           6.39959 |               nan |         6.91913 | TCGA-ET-A25M |               nan |                     nan |          1.45799  | 3.4845   | 0.536337  | 0.726797  | 0.957176  | 10.6822  | 10.9951  | 10.2333  | 11.4777  |
| TCGA-DJ-A4V5-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | mid         | no            | no                | Stage IVA          |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Male   |  55.0527 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  0.158631    |                    -0.158631    |                        nan | histology_proxy | DM1                |           7.07622 |               nan |         7.85746 | TCGA-DJ-A4V5 |                 0 |                       0 |          3.07046  | 3.80054  | 2.73639   | 2.33859   | 3.39673   | 10.4453  | 10.7648  | 10.4563  | 11.8103  |
| TCGA-DE-A4M9-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | high        | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Male   |  28.4298 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  1.03978     |                    -1.03978     |                        nan | histology_proxy | DM2                |           7.95737 |               nan |         8.83342 | TCGA-DE-A4M9 |                 0 |                       0 |          0.374329 | 1.25107  | 0.276033  | 0.0273936 | 3.59909   |  7.90309 |  7.3696  |  8.33327 |  9.62239 |
| TCGA-ET-A39N-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female |  35.3292 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  0.484382    |                    -0.484382    |                        nan | histology_proxy | DM2                |           7.40197 |               nan |         8.33539 | TCGA-ET-A39N |                 0 |                       0 |          2.31254  | 2.05003  | 0.519905  | 0.120049  | 2.76108   | 10.0741  |  9.06469 |  8.71653 | 10.5462  |
| TCGA-E8-A416-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | high        | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female |  51.833  | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  1.38522     |                    -1.38522     |                        nan | histology_proxy | DM2                |           8.30281 |               nan |         9.20746 | TCGA-E8-A416 |                 0 |                       0 |          2.37402  | 3.20263  | 1.69456   | 1.01201   | 3.16808   |  8.96605 |  9.18986 |  8.6278  | 10.6151  |
| TCGA-KS-A4ID-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | unknown         | RAS_like            | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Female |  55.8877 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  0.545188    |                    -0.545188    |                        nan | histology_proxy | DM1                |           7.46278 |               nan |         8.62262 | TCGA-KS-A4ID |                 0 |                       0 |          3.55695  | 4.18752  | 3.22754   | 1.96566   | 4.1244    | 11.6266  | 11.867   | 11.3372  | 12.924   |
| TCGA-DJ-A13U-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | mid         | no            | no                | Stage III          |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Male   |  60.5886 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs | -0.0483846   |                     0.0483846   |                        nan | histology_proxy | DM1                |           6.86921 |               nan |         7.42806 | TCGA-DJ-A13U |               nan |                     nan |          1.50784  | 4.09683  | 0.5745    | 1.10561   | 1.731     | 10.4111  | 11.0796  | 10.5805  | 12.0631  |
| TCGA-EL-A3T2-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | high        | no            | no                | nan                |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female |  58.4504 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  1.29095     |                    -1.29095     |                        nan | histology_proxy | DM2                |           8.20854 |               nan |         9.03334 | TCGA-EL-A3T2 |               nan |                     nan |          1.11392  | 2.19088  | 1.07223   | 0.0997706 | 0.558346  |  9.27003 |  8.4039  |  7.86335 | 10.3602  |
| TCGA-DE-A4MA-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | low         | no            | no                | Stage III          |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female |  52.293  | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs | -1.64479     |                     1.64479     |                        nan | histology_proxy | DM1                |           5.2728  |               nan |         6.10887 | TCGA-DE-A4MA |               nan |                     nan |          1.48861  | 2.77727  | 1.26094   | 2.09036   | 1.36065   | 11.066   | 11.9338  | 11.1282  | 11.8597  |
| TCGA-DJ-A2Q0-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | mid         | no            | no                | Stage II           |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female |  57.4209 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs | -0.241805    |                     0.241805    |                        nan | histology_proxy | DM2                |           6.67579 |               nan |         7.38092 | TCGA-DJ-A2Q0 |                 0 |                       0 |          1.29118  | 5.64954  | 0.0499676 | 0.144993  | 0.508834  | 10.9489  | 11.439   | 11.0801  | 12.225   |
| TCGA-EL-A3T6-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | low         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female |  34.9815 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs | -0.611075    |                     0.611075    |                        nan | histology_proxy | DM1                |           6.30652 |               nan |         6.8982  | TCGA-EL-A3T6 |               nan |                     nan |          3.84689  | 4.74243  | 2.68743   | 2.6542    | 4.16177   | 11.406   | 11.4432  | 11.2644  | 13.4979  |
| TCGA-EM-A3O8-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | unknown         | RAS_like            | high        | no            | no                | Stage I            |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Female |  33.963  | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  1.35824     |                    -1.35824     |                        nan | histology_proxy | DM2                |           8.27583 |               nan |         9.02142 | TCGA-EM-A3O8 |               nan |                     nan |          0.947813 | 2.08706  | 1.54226   | 0.198011  | 0.330504  |  8.68282 |  8.47862 |  7.84663 | 10.644   |
| TCGA-DJ-A2Q2-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | unknown         | RAS_like            | high        | no            | no                | Stage III          |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Female |  57.3607 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  1.64468     |                    -1.64468     |                        nan | histology_proxy | DM2                |           8.56227 |               nan |         9.69552 | TCGA-DJ-A2Q2 |               nan |                     nan |          2.30305  | 2.37084  | 1.11694   | 1.015     | 1.39371   | 10.4143  | 10.1307  |  9.72643 | 12.7487  |
| TCGA-BJ-A0ZG-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | unknown         | RAS_like            | high        | no            | no                | Stage I            |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Male   |  80.4244 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  1.70407     |                    -1.70407     |                        nan | histology_proxy | DM2                |           8.62166 |               nan |         9.53485 | TCGA-BJ-A0ZG |               nan |                     nan |          2.38312  | 2.54286  | 0.879965  | 0.421919  | 2.77161   |  9.57331 |  9.52805 |  8.9026  | 10.9132  |
| TCGA-EM-A2OW-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | unknown         | RAS_like            | high        | no            | no                | Stage I            |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Female |  37.6756 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  1.59551     |                    -1.59551     |                        nan | histology_proxy | DM2                |           8.5131  |               nan |        10.0407  | TCGA-EM-A2OW |                 0 |                       0 |          1.52219  | 3.16084  | 0.696972  | 0.277963  | 2.21055   |  9.29244 |  9.69774 |  9.1976  | 11.5701  |
| TCGA-BJ-A3F0-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | high        | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female |  64.4627 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  1.93187     |                    -1.93187     |                        nan | histology_proxy | DM2                |           8.84946 |               nan |        10.8516  | TCGA-BJ-A3F0 |                 0 |                       0 |          2.26825  | 1.19085  | 0.325616  | 0         | 4.25387   |  9.41142 |  9.9289  |  9.80266 | 11.0015  |
| TCGA-EL-A3CY-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | low         | no            | no                | nan                |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Male   | nan      | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs | -0.728505    |                     0.728505    |                        nan | histology_proxy | DM1                |           6.18909 |               nan |         7.04248 | TCGA-EL-A3CY |                 0 |                       0 |          4.30411  | 4.06207  | 3.41655   | 4.37119   | 5.9273    | 10.6931  | 11.0637  | 10.3786  | 12.2385  |
| TCGA-EM-A3AQ-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female |  85.4155 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs | -0.427982    |                     0.427982    |                        nan | histology_proxy | DM1                |           6.48961 |               nan |         7.03765 | TCGA-EM-A3AQ |                 0 |                       0 |          2.33784  | 2.45413  | 0.507546  | 0.507546  | 0.643555  | 10.4523  | 10.159   | 10.3399  | 11.2277  |
| TCGA-J8-A4HW-06A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Metastatic    | tumor             | unknown             | unknown         | unknown             | mid         | no            | no                | nan                |              nan |                 nan | Adenocarcinoma, NOS                     | medium             | Female |  44.2327 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs | -0.140721    |                     0.140721    |                        nan | ambiguous       | DM1                |           6.77687 |               nan |         7.52017 | TCGA-J8-A4HW |                 0 |                       0 |          3.14192  | 3.5326   | 2.83774   | 2.93823   | 3.03567   |  9.95081 | 10.24    |  9.93827 | 11.416   |
| TCGA-J8-A4HW-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | unknown             | unknown         | unknown             | mid         | no            | no                | nan                |              nan |                 nan | Adenocarcinoma, NOS                     | medium             | Female |  44.2327 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  0.28295     |                    -0.28295     |                        nan | ambiguous       | DM1                |           7.20054 |               nan |         8.31276 | TCGA-J8-A4HW |                 0 |                       0 |          1.00841  | 1.94306  | 0.665653  | 0.072766  | 0.072766  |  8.94781 |  8.84439 |  8.72318 | 10.379   |
| TCGA-BJ-A291-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | unknown             | unknown         | unknown             | high        | no            | no                | Stage 0a           |              nan |                 nan | Not Reported                            | medium             | Female | nan      | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  0.737047    |                    -0.737047    |                        nan | ambiguous       | DM2                |           7.65464 |               nan |         8.64342 | TCGA-BJ-A291 |               nan |                     nan |          1.50487  | 1.98915  | 0.367567  | 0.0484291 | 4.70363   |  9.72741 | 10.2594  | 10.2287  | 11.8091  |
| TCGA-CE-A3MD-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Male   |  31.2088 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs | -0.000831004 |                     0.000831004 |                        nan | histology_proxy | DM1                |           6.91676 |               nan |         7.60137 | TCGA-CE-A3MD |                 0 |                       0 |          2.77592  | 4.11513  | 1.08739   | 0.931438  | 1.63443   | 11.6585  | 12.056   | 10.6968  | 12.6417  |
| TCGA-EL-A3ZN-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female |  28.1068 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  0.137906    |                    -0.137906    |                        nan | histology_proxy | DM1                |           7.0555  |               nan |         7.52852 | TCGA-EL-A3ZN |                 0 |                       0 |          2.03035  | 2.65297  | 0.843455  | 0.641223  | 1.2894    |  9.90845 |  9.34216 |  8.86132 | 11.1911  |
| TCGA-BJ-A2N9-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | high        | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female |  42.193  | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  1.27548     |                    -1.27548     |                        nan | histology_proxy | DM2                |           8.19307 |               nan |         9.36015 | TCGA-BJ-A2N9 |                 0 |                       0 |          3.26485  | 2.86733  | 2.7452    | 1.07497   | 4.32422   | 10.5725  | 11.0634  | 10.5267  | 11.7405  |
| TCGA-DJ-A4UP-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | high        | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female |  15.4825 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  0.682614    |                    -0.682614    |                        nan | histology_proxy | DM1                |           7.6002  |               nan |         8.48931 | TCGA-DJ-A4UP |                 0 |                       0 |          3.52344  | 3.38774  | 3.58265   | 3.1539    | 3.46025   | 10.6449  | 11.0328  | 10.9081  | 11.6986  |
| TCGA-E8-A44M-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female |  23.3347 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  0.482463    |                    -0.482463    |                        nan | histology_proxy | DM1                |           7.40005 |               nan |         8.3633  | TCGA-E8-A44M |                 0 |                       0 |          4.90855  | 4.98962  | 4.93611   | 3.50786   | 5.66814   | 11.7933  | 12.8073  | 12.0333  | 13.3879  |
| TCGA-EM-A3O6-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | unknown         | RAS_like            | high        | no            | no                | Stage II           |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Female |  54.3628 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  1.48252     |                    -1.48252     |                        nan | histology_proxy | DM2                |           8.40011 |               nan |         9.75449 | TCGA-EM-A3O6 |                 0 |                       0 |          1.79882  | 1.92361  | 0         | 0.100643  | 1.1637    |  9.2427  |  9.35026 |  9.13379 | 11.255   |
| TCGA-EM-A4FH-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | unknown         | RAS_like            | high        | no            | no                | Stage III          |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Male   |  69.8152 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  1.61709     |                    -1.61709     |                        nan | histology_proxy | DM2                |           8.53469 |               nan |         9.52686 | TCGA-EM-A4FH |                 0 |                       0 |          1.5767   | 1.99965  | 0.243551  | 0.0402511 | 1.94769   |  8.76247 |  8.74231 |  9.02927 | 10.89    |
| TCGA-FY-A40N-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | unknown         | RAS_like            | high        | no            | no                | Stage I            |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Female |  55.6331 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  1.11602     |                    -1.11602     |                        nan | histology_proxy | DM2                |           8.03361 |               nan |         9.31664 | TCGA-FY-A40N |                 0 |                       0 |          2.6512   | 3.77594  | 2.67262   | 1.43119   | 2.59371   | 10.2865  | 10.5696  | 10.2531  | 12.2893  |
| TCGA-FY-A3WA-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | unknown         | RAS_like            | high        | no            | no                | Stage I            |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Female |  52      | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  1.71443     |                    -1.71443     |                        nan | histology_proxy | DM2                |           8.63202 |               nan |         9.52753 | TCGA-FY-A3WA |                 0 |                       0 |          1.66395  | 2.08467  | 1.24737   | 0.075002  | 0.32593   |  7.95351 |  7.57327 |  7.28075 |  9.8444  |
| TCGA-ET-A2MX-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | low         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Male   |  27.499  | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs | -1.15993     |                     1.15993     |                        nan | histology_proxy | DM1                |           5.75766 |               nan |         6.22166 | TCGA-ET-A2MX |                 0 |                       0 |          0.780114 | 4.07524  | 0.236095  | 0.112124  | 0.50615   | 10.2598  | 10.9956  |  9.80295 | 11.8258  |
| TCGA-FY-A3I5-01B | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | unknown         | RAS_like            | high        | no            | no                | Stage III          |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Female |  64.7611 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  1.47654     |                    -1.47654     |                        nan | histology_proxy | DM2                |           8.39413 |               nan |         9.5674  | TCGA-FY-A3I5 |                 0 |                       0 |          1.79016  | 4.90811  | 1.4369    | 2.10465   | 2.91279   |  8.88371 |  9.5269  |  8.5177  | 11.3859  |
| TCGA-EM-A2OV-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | unknown         | RAS_like            | high        | no            | no                | Stage II           |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Female |  64.7666 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  1.34848     |                    -1.34848     |                        nan | histology_proxy | DM2                |           8.26607 |               nan |         9.38485 | TCGA-EM-A2OV |                 1 |                       1 |          1.86344  | 2.31324  | 0.222407  | 0.0216966 | 1.04312   |  8.63969 |  9.13736 |  8.91299 | 10.4868  |
| TCGA-EM-A3AI-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | unknown         | RAS_like            | high        | no            | no                | Stage II           |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Female |  67.9535 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  1.68016     |                    -1.68016     |                        nan | histology_proxy | DM2                |           8.59775 |               nan |         9.79107 | TCGA-EM-A3AI |                 0 |                       0 |          0.792159 | 1.33603  | 0.133975  | 0         | 1.14822   |  8.1537  |  8.70069 |  8.4584  | 10.0087  |
| TCGA-FY-A4B0-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | unknown         | RAS_like            | high        | no            | no                | Stage I            |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Male   |  76.512  | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  1.22771     |                    -1.22771     |                        nan | histology_proxy | DM2                |           8.1453  |               nan |         9.35685 | TCGA-FY-A4B0 |               nan |                     nan |          1.63314  | 2.98861  | 1.25789   | 0.467962  | 1.95745   |  8.7535  |  9.02307 |  8.24455 | 10.1586  |
| TCGA-FY-A40M-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | low         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female |  51.258  | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs | -0.752259    |                     0.752259    |                        nan | histology_proxy | DM1                |           6.16533 |               nan |         6.99326 | TCGA-FY-A40M |                 0 |                       0 |          1.05638  | 2.20033  | 2.12544   | 0.0867803 | 0.535537  |  8.49742 |  9.95372 |  9.73302 | 10.0683  |
| TCGA-ET-A4KN-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | low         | no            | no                | Stage III          |              nan |                 nan | Papillary carcinoma, columnar cell      | high               | Female |  51.7153 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs | -2.84402     |                     2.84402     |                        nan | histology_proxy | DM1                |           4.07357 |               nan |         4.96019 | TCGA-ET-A4KN |                 0 |                       0 |          3.94966  | 2.56264  | 3.38779   | 2.48183   | 3.04285   | 12.6946  | 12.7836  | 11.6991  | 11.1088  |
| TCGA-ET-A40P-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | unknown         | RAS_like            | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Female |  32.4244 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  0.520114    |                    -0.520114    |                        nan | histology_proxy | DM1                |           7.4377  |               nan |         8.28714 | TCGA-ET-A40P |               nan |                     nan |          2.72343  | 2.51246  | 3.1774    | 0.750638  | 0.825659  |  9.99533 | 10.4813  |  9.71131 | 11.1298  |
| TCGA-FE-A22Z-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | unknown         | RAS_like            | low         | no            | no                | nan                |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Female |  66.6886 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs | -1.46452     |                     1.46452     |                        nan | histology_proxy | DM1                |           5.45307 |               nan |         5.91866 | TCGA-FE-A22Z |               nan |                     nan |          1.00351  | 3.635    | 0.118239  | 0.596892  | 0.513032  |  9.75841 |  9.45729 |  8.91194 | 10.892   |
| TCGA-FK-A3SG-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female |  21.1088 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs | -0.145206    |                     0.145206    |                        nan | histology_proxy | DM1                |           6.77239 |               nan |         7.62733 | TCGA-FK-A3SG |                 0 |                       0 |          1.31869  | 3.33566  | 1.00208   | 0.531795  | 0.400938  | 10.4821  | 10.144   |  9.88231 | 11.8442  |
| TCGA-EM-A22Q-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | unknown         | RAS_like            | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Male   |  36.5257 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  0.568189    |                    -0.568189    |                        nan | histology_proxy | DM1                |           7.48578 |               nan |         8.81481 | TCGA-EM-A22Q |               nan |                     nan |          1.80933  | 0.680788 | 0.680788  | 0.15009   | 2.1007    | 10.0904  |  9.63303 |  8.22585 |  8.44138 |
| TCGA-DE-A4MC-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | low         | no            | no                | nan                |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female |  44.7036 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs | -1.16468     |                     1.16468     |                        nan | histology_proxy | DM1                |           5.75291 |               nan |         6.05471 | TCGA-DE-A4MC |               nan |                     nan |          2.6435   | 4.1209   | 0.988352  | 3.02069   | 2.42699   | 11.2968  | 11.5036  | 10.7122  | 13.0263  |
| TCGA-DJ-A13S-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | unknown         | RAS_like            | high        | no            | no                | Stage I            |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Female |  19.436  | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  0.842244    |                    -0.842244    |                        nan | histology_proxy | DM2                |           7.75983 |               nan |         8.79407 | TCGA-DJ-A13S |               nan |                     nan |          1.08683  | 3.28956  | 1.78609   | 0.0381027 | 0.0934339 |  8.63326 |  9.80262 |  9.47191 | 11.2573  |
| TCGA-EM-A3SY-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | unknown         | RAS_like            | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Male   |  36.8542 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs | -0.194991    |                     0.194991    |                        nan | histology_proxy | DM2                |           6.7226  |               nan |         7.64784 | TCGA-EM-A3SY |                 0 |                       0 |          1.56291  | 4.47498  | 0.779811  | 0.219452  | 0.876873  | 11.0291  | 10.9153  | 10.6943  | 13.6045  |
| TCGA-EL-A3TB-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | mid         | no            | no                | Stage III          |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female |  47.0773 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  0.486016    |                    -0.486016    |                        nan | histology_proxy | DM1                |           7.40361 |               nan |         8.17911 | TCGA-EL-A3TB |               nan |                     nan |          2.02936  | 3.23903  | 1.58683   | 0.889119  | 2.16069   | 10.7454  | 10.9455  | 10.2707  | 12.1126  |
| TCGA-BJ-A45C-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | high        | no            | no                | Stage III          |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Male   |  78.0014 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  1.00152     |                    -1.00152     |                        nan | histology_proxy | DM2                |           7.91911 |               nan |         8.8994  | TCGA-BJ-A45C |               nan |                     nan |          1.29729  | 2.82482  | 0.587478  | 0         | 0.216058  |  9.32897 |  8.45422 |  8.42642 |  9.75937 |
| TCGA-DJ-A2PX-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | FVPTC               | unknown         | RAS_like            | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary carcinoma, follicular variant | high               | Female |  54.9979 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  0.160796    |                    -0.160796    |                        nan | histology_proxy | DM1                |           7.07839 |               nan |         7.66879 | TCGA-DJ-A2PX |                 0 |                       0 |          4.26099  | 3.77181  | 3.25107   | 2.16807   | 3.71631   | 10.6272  | 10.9515  | 10.7004  | 12.3857  |
| TCGA-EL-A3ZS-01A | TCGA-THCA | Illumina HiSeq / GDC STAR Counts | bulk RNA-seq | Primary Tumor | tumor             | cPTC                | unknown         | BRAF_like           | mid         | no            | no                | Stage I            |              nan |                 nan | Papillary adenocarcinoma, NOS           | high               | Female |  22.1109 | no                  | no                   | no                     | GDC API cases + masked somatic mutation MAFs |  0.371598    |                    -0.371598    |                        nan | histology_proxy | DM1                |           7.28919 |               nan |         8.45883 | TCGA-EL-A3ZS |               nan |                     nan |          3.33537  | 3.32306  | 3.48878   | 3.24166   | 4.20417   | 10.5952  | 11.4584  | 10.97    | 12.5984  |

### 8.3 Immune evasion genes
- file: `/home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/A5_immune_evasion_genes.tsv`
| sample_id        | v17_dark_cluster   |    CD274 |     PDCD1 |     CTLA4 |     IDO1 |    HLA-A |    HLA-B |    HLA-C |      B2M |   cytolytic_score |
|:-----------------|:-------------------|---------:|----------:|----------:|---------:|---------:|---------:|---------:|---------:|------------------:|
| TCGA-FK-A3SE-01A | DM1                | 1.89636  | 0.557719  | 0.420972  | 0.557719 |  8.94436 |  8.786   |  8.77273 | 10.7395  |          0.838522 |
| TCGA-FY-A2QD-01A | DM2                | 1.5418   | 0.365738  | 0.0257854 | 2.28607  |  8.6665  |  7.44263 |  7.0048  | 10.3523  |          0.553785 |
| TCGA-EL-A4KD-01A | DM1                | 3.02656  | 1.3217    | 1.13347   | 1.56118  |  9.54072 |  9.91064 |  9.70173 | 10.8807  |          2.63608  |
| TCGA-BJ-A192-01A | DM2                | 2.82566  | 0.0240007 | 0.0476087 | 1.32755  |  9.31237 |  9.16514 |  8.84897 |  9.5853  |          1.32534  |
| TCGA-ET-A40S-01A | DM1                | 4.02443  | 2.61006   | 2.39935   | 2.31553  | 10.789   | 11.8027  | 11.457   | 12.5629  |          3.42978  |
| TCGA-ET-A39L-01A | DM1                | 1.60287  | 0.966218  | 0.100643  | 0.644832 |  9.13171 |  9.20913 |  9.35162 | 10.3672  |          2.37656  |
| TCGA-KS-A4IB-01A | DM1                | 3.47426  | 0.683118  | 0.190688  | 0.489273 |  9.864   | 10.721   |  9.7035  | 11.5164  |          1.33339  |
| TCGA-DJ-A2PO-01A | DM1                | 4.1104   | 2.80782   | 1.43315   | 1.62828  | 11.1348  | 11.8827  | 10.9039  | 12.4521  |          3.72882  |
| TCGA-EM-A2CJ-01A | DM2                | 1.58     | 0.763954  | 0.0249628 | 0.864186 |  8.19517 |  8.43376 |  8.84677 |  9.85105 |          0.710587 |
| TCGA-DJ-A1QL-01A | DM2                | 1.36615  | 0.832221  | 0.0242711 | 0.950928 |  9.34672 |  9.89435 |  9.01656 | 10.0544  |          2.39803  |
| TCGA-EM-A1CS-01A | DM1                | 3.35319  | 0.539891  | 0.136562  | 1.16872  |  9.16148 |  9.85622 |  9.35173 | 10.6516  |          1.53293  |
| TCGA-CE-A484-01A | DM1                | 1.26336  | 2.44381   | 0.367771  | 0.423952 |  9.4519  |  9.49423 |  9.06358 | 10.3684  |          1.67212  |
| TCGA-EL-A3ZP-01A | DM1                | 3.89232  | 3.58949   | 2.1773    | 4.02435  | 10.8179  | 11.5279  | 10.7554  | 12.177   |          4.08105  |
| TCGA-BJ-A0YZ-01A | DM2                | 0.944485 | 0.407402  | 0.0259182 | 4.27118  |  8.98814 |  8.56389 |  9.13352 | 10.3622  |          1.31121  |
| TCGA-FY-A3NM-01A | DM2                | 2.01097  | 0.175578  | 0.229657  | 1.46358  |  7.71764 |  8.42761 |  8.18974 |  9.83521 |          1.05508  |
| TCGA-BJ-A3PR-01A | DM1                | 4.45473  | 5.55695   | 3.91943   | 5.40976  | 11.3847  | 12.3287  | 10.9473  | 13.2891  |          5.28118  |
| TCGA-EL-A3CO-01A | DM1                | 4.87503  | 1.97346   | 1.77862   | 3.98344  | 10.3312  | 10.5465  |  9.64063 | 11.556   |          3.61151  |
| TCGA-CE-A482-01A | DM1                | 3.42153  | 1.05377   | 0.713111  | 0.893474 | 10.0325  |  9.90853 | 10.0734  | 10.9514  |          1.68446  |
| TCGA-EM-A3FQ-01A | DM1                | 3.6919   | 2.4086    | 2.31915   | 3.51788  | 10.7221  | 11.0455  | 10.7051  | 11.6551  |          3.11261  |
| TCGA-EL-A3GO-01A | DM2                | 1.6536   | 0.187088  | 0.332988  | 1.02749  |  8.49518 |  9.26351 |  8.54644 | 10.4997  |          2.42594  |
| TCGA-DJ-A3UZ-01A | DM1                | 3.49142  | 2.26853   | 2.02758   | 3.35269  | 10.3334  | 10.8213  |  9.89604 | 12.292   |          3.14528  |
| TCGA-E8-A432-01A | DM1                | 3.67658  | 1.82633   | 0.755854  | 2.11472  | 10.086   |  9.9902  |  9.9968  | 11.9653  |          2.01167  |
| TCGA-DJ-A4UT-01A | DM1                | 4.21709  | 3.84097   | 2.95011   | 4.72238  | 11.3085  | 11.4732  | 10.8588  | 12.6204  |          4.1742   |
| TCGA-ET-A40R-01A | DM1                | 2.3605   | 2.31361   | 1.74155   | 2.02168  | 10.2623  | 10.6706  |  9.4312  | 11.7106  |          2.31913  |
| TCGA-EM-A3FR-01A | DM2                | 2.05461  | 0.45039   | 0.594345  | 1.71986  | 10.4909  |  8.71298 |  8.10718 | 10.891   |          1.59883  |
| TCGA-DJ-A4UQ-01A | DM1                | 2.95887  | 1.51783   | 1.09227   | 1.1836   |  9.7953  |  9.90735 |  9.50239 | 11.113   |          2.30964  |
| TCGA-ET-A4KQ-01A | DM2                | 1.58223  | 0.773004  | 0.0677156 | 0.405314 |  8.49104 |  8.8669  |  8.44995 | 10.2597  |          1.50554  |
| TCGA-CE-A481-01A | DM1                | 3.42083  | 2.75479   | 1.60313   | 3.41412  | 10.2826  | 10.415   |  9.88526 | 11.7685  |          3.21979  |
| TCGA-EM-A2CP-01A | DM1                | 4.66712  | 4.66712   | 3.03829   | 5.55109  | 11.2894  | 11.8604  | 11.8164  | 12.9118  |          5.05209  |
| TCGA-EL-A3H3-01A | DM1                | 4.23398  | 2.42324   | 1.58135   | 2.74363  | 11.0029  | 11.848   | 11.002   | 13.0102  |          3.28696  |
| TCGA-EL-A3CS-01A | DM1                | 4.01225  | 0.448196  | 0.674625  | 1.01917  | 10.0589  | 10.2652  |  9.99991 | 11.0548  |          0.81823  |
| TCGA-FK-A3SD-01A | DM2                | 2.03151  | 1.37869   | 0.0201115 | 5.03482  |  9.71723 |  9.07666 |  8.8891  | 10.5986  |          0.392661 |
| TCGA-DJ-A13W-01A | DM1                | 0.611683 | 0.576884  | 0.285106  | 2.83883  |  8.90942 |  8.42171 |  6.75186 |  7.8071  |          1.42938  |
| TCGA-DJ-A13R-01A | DM1                | 3.24778  | 4.49662   | 1.79485   | 6.04866  | 11.5168  | 11.1044  | 10.3054  | 11.9406  |          5.10047  |
| TCGA-EL-A3ZK-01A | DM1                | 5.8787   | 1.12687   | 1.96561   | 2.11133  |  9.96024 | 10.4567  |  9.94812 | 11.8822  |          1.95303  |
| TCGA-ET-A3DS-01A | DM2                | 2.02254  | 1.37098   | 0.587501  | 1.59754  |  8.3752  |  8.73162 |  9.22236 | 10.3235  |          1.64834  |
| TCGA-EL-A3D4-01A | DM1                | 2.60715  | 1.022     | 0.112155  | 0.822112 | 10.34    |  9.58341 |  9.52904 | 10.9795  |          1.89577  |
| TCGA-EM-A3AN-01A | DM1                | 2.80769  | 0.98378   | 0.450607  | 2.2306   |  9.85776 |  9.85459 |  9.98242 | 11.5034  |          1.89258  |
| TCGA-DE-A0Y2-01A | DM1                | 3.71768  | 0.640706  | 0.399971  | 0.451439 | 10.1925  | 10.3802  |  9.96885 | 11.741   |          0.929393 |
| TCGA-CE-A27D-01A | DM1                | 4.02803  | 4.00495   | 2.24097   | 6.20814  | 11.3559  | 11.4857  | 10.8933  | 12.7246  |          4.41596  |
| TCGA-BJ-A0Z0-01A | DM2                | 1.82566  | 0.264691  | 0         | 0.714356 |  8.77174 |  7.97964 |  7.48515 | 10.1068  |          0.801932 |
| TCGA-FK-A4UB-01A | DM1                | 2.31998  | 1.46108   | 0.572387  | 1.28005  | 10.1503  | 10.2387  |  8.60565 | 10.5602  |          1.48553  |
| TCGA-H2-A421-01A | DM1                | 4.1397   | 3.3995    | 2.67348   | 2.49112  | 11.2538  | 11.9274  | 11.0112  | 12.7889  |          3.35964  |
| TCGA-CE-A483-01A | DM1                | 3.48745  | 0.367617  | 0.334795  | 0.399709 |  9.70301 | 10.1817  |  9.1897  | 11.8955  |          1.55166  |
| TCGA-FY-A76V-01A | DM1                | 3.24788  | 0.907667  | 0.535896  | 1.24684  |  8.63893 |  8.8831  |  8.67799 | 10.5831  |          1.9471   |
| TCGA-DO-A2HM-01B | DM1                | 5.01979  | 2.91461   | 3.36734   | 5.64062  | 10.5497  | 11.3927  | 10.5274  | 12.4981  |          4.27281  |
| TCGA-EL-A3T0-01A | DM1                | 3.28772  | 1.55712   | 0.914156  | 1.51788  |  9.693   |  9.99665 |  9.31306 | 11.3485  |          2.29402  |
| TCGA-BJ-A0ZJ-01A | DM1                | 3.3286   | 3.77705   | 2.05736   | 3.02592  | 11.0077  | 11.4858  | 10.6053  | 12.422   |          3.8632   |
| TCGA-EL-A3T9-01A | DM1                | 2.66079  | 0.490295  | 2.27184   | 1.67577  |  8.9158  |  9.1876  |  9.05455 | 10.0959  |          0.81239  |
| TCGA-EM-A1YD-01A | DM2                | 2.13086  | 0.667018  | 0.526911  | 2.92491  |  9.02557 |  9.37094 |  9.365   | 10.4999  |          1.96123  |

...

| sample_id        | v17_dark_cluster   |    CD274 |     PDCD1 |     CTLA4 |      IDO1 |    HLA-A |    HLA-B |    HLA-C |      B2M |   cytolytic_score |
|:-----------------|:-------------------|---------:|----------:|----------:|----------:|---------:|---------:|---------:|---------:|------------------:|
| TCGA-ET-A39T-01A | DM1                | 4.43952  | 0.548777  | 0.28944   | 0.584841  | 10.3141  | 10.9952  | 10.0005  | 12.2804  |          1.81449  |
| TCGA-DE-A69K-01A | DM1                | 3.03546  | 0.702847  | 0.575701  | 1.05581   | 10.3214  | 10.366   | 10.0911  | 11.432   |          1.88734  |
| TCGA-ET-A25M-01A | DM1                | 3.4845   | 0.536337  | 0.726797  | 0.957176  | 10.6822  | 10.9951  | 10.2333  | 11.4777  |          1.45799  |
| TCGA-DJ-A4V5-01A | DM1                | 3.80054  | 2.73639   | 2.33859   | 3.39673   | 10.4453  | 10.7648  | 10.4563  | 11.8103  |          3.07046  |
| TCGA-DE-A4M9-01A | DM2                | 1.25107  | 0.276033  | 0.0273936 | 3.59909   |  7.90309 |  7.3696  |  8.33327 |  9.62239 |          0.374329 |
| TCGA-ET-A39N-01A | DM2                | 2.05003  | 0.519905  | 0.120049  | 2.76108   | 10.0741  |  9.06469 |  8.71653 | 10.5462  |          2.31254  |
| TCGA-E8-A416-01A | DM2                | 3.20263  | 1.69456   | 1.01201   | 3.16808   |  8.96605 |  9.18986 |  8.6278  | 10.6151  |          2.37402  |
| TCGA-KS-A4ID-01A | DM1                | 4.18752  | 3.22754   | 1.96566   | 4.1244    | 11.6266  | 11.867   | 11.3372  | 12.924   |          3.55695  |
| TCGA-DJ-A13U-01A | DM1                | 4.09683  | 0.5745    | 1.10561   | 1.731     | 10.4111  | 11.0796  | 10.5805  | 12.0631  |          1.50784  |
| TCGA-EL-A3T2-01A | DM2                | 2.19088  | 1.07223   | 0.0997706 | 0.558346  |  9.27003 |  8.4039  |  7.86335 | 10.3602  |          1.11392  |
| TCGA-DE-A4MA-01A | DM1                | 2.77727  | 1.26094   | 2.09036   | 1.36065   | 11.066   | 11.9338  | 11.1282  | 11.8597  |          1.48861  |
| TCGA-DJ-A2Q0-01A | DM2                | 5.64954  | 0.0499676 | 0.144993  | 0.508834  | 10.9489  | 11.439   | 11.0801  | 12.225   |          1.29118  |
| TCGA-EL-A3T6-01A | DM1                | 4.74243  | 2.68743   | 2.6542    | 4.16177   | 11.406   | 11.4432  | 11.2644  | 13.4979  |          3.84689  |
| TCGA-EM-A3O8-01A | DM2                | 2.08706  | 1.54226   | 0.198011  | 0.330504  |  8.68282 |  8.47862 |  7.84663 | 10.644   |          0.947813 |
| TCGA-DJ-A2Q2-01A | DM2                | 2.37084  | 1.11694   | 1.015     | 1.39371   | 10.4143  | 10.1307  |  9.72643 | 12.7487  |          2.30305  |
| TCGA-BJ-A0ZG-01A | DM2                | 2.54286  | 0.879965  | 0.421919  | 2.77161   |  9.57331 |  9.52805 |  8.9026  | 10.9132  |          2.38312  |
| TCGA-EM-A2OW-01A | DM2                | 3.16084  | 0.696972  | 0.277963  | 2.21055   |  9.29244 |  9.69774 |  9.1976  | 11.5701  |          1.52219  |
| TCGA-BJ-A3F0-01A | DM2                | 1.19085  | 0.325616  | 0         | 4.25387   |  9.41142 |  9.9289  |  9.80266 | 11.0015  |          2.26825  |
| TCGA-EL-A3CY-01A | DM1                | 4.06207  | 3.41655   | 4.37119   | 5.9273    | 10.6931  | 11.0637  | 10.3786  | 12.2385  |          4.30411  |
| TCGA-EM-A3AQ-01A | DM1                | 2.45413  | 0.507546  | 0.507546  | 0.643555  | 10.4523  | 10.159   | 10.3399  | 11.2277  |          2.33784  |
| TCGA-J8-A4HW-06A | DM1                | 3.5326   | 2.83774   | 2.93823   | 3.03567   |  9.95081 | 10.24    |  9.93827 | 11.416   |          3.14192  |
| TCGA-J8-A4HW-01A | DM1                | 1.94306  | 0.665653  | 0.072766  | 0.072766  |  8.94781 |  8.84439 |  8.72318 | 10.379   |          1.00841  |
| TCGA-BJ-A291-01A | DM2                | 1.98915  | 0.367567  | 0.0484291 | 4.70363   |  9.72741 | 10.2594  | 10.2287  | 11.8091  |          1.50487  |
| TCGA-CE-A3MD-01A | DM1                | 4.11513  | 1.08739   | 0.931438  | 1.63443   | 11.6585  | 12.056   | 10.6968  | 12.6417  |          2.77592  |
| TCGA-EL-A3ZN-01A | DM1                | 2.65297  | 0.843455  | 0.641223  | 1.2894    |  9.90845 |  9.34216 |  8.86132 | 11.1911  |          2.03035  |
| TCGA-BJ-A2N9-01A | DM2                | 2.86733  | 2.7452    | 1.07497   | 4.32422   | 10.5725  | 11.0634  | 10.5267  | 11.7405  |          3.26485  |
| TCGA-DJ-A4UP-01A | DM1                | 3.38774  | 3.58265   | 3.1539    | 3.46025   | 10.6449  | 11.0328  | 10.9081  | 11.6986  |          3.52344  |
| TCGA-E8-A44M-01A | DM1                | 4.98962  | 4.93611   | 3.50786   | 5.66814   | 11.7933  | 12.8073  | 12.0333  | 13.3879  |          4.90855  |
| TCGA-EM-A3O6-01A | DM2                | 1.92361  | 0         | 0.100643  | 1.1637    |  9.2427  |  9.35026 |  9.13379 | 11.255   |          1.79882  |
| TCGA-EM-A4FH-01A | DM2                | 1.99965  | 0.243551  | 0.0402511 | 1.94769   |  8.76247 |  8.74231 |  9.02927 | 10.89    |          1.5767   |
| TCGA-FY-A40N-01A | DM2                | 3.77594  | 2.67262   | 1.43119   | 2.59371   | 10.2865  | 10.5696  | 10.2531  | 12.2893  |          2.6512   |
| TCGA-FY-A3WA-01A | DM2                | 2.08467  | 1.24737   | 0.075002  | 0.32593   |  7.95351 |  7.57327 |  7.28075 |  9.8444  |          1.66395  |
| TCGA-ET-A2MX-01A | DM1                | 4.07524  | 0.236095  | 0.112124  | 0.50615   | 10.2598  | 10.9956  |  9.80295 | 11.8258  |          0.780114 |
| TCGA-FY-A3I5-01B | DM2                | 4.90811  | 1.4369    | 2.10465   | 2.91279   |  8.88371 |  9.5269  |  8.5177  | 11.3859  |          1.79016  |
| TCGA-EM-A2OV-01A | DM2                | 2.31324  | 0.222407  | 0.0216966 | 1.04312   |  8.63969 |  9.13736 |  8.91299 | 10.4868  |          1.86344  |
| TCGA-EM-A3AI-01A | DM2                | 1.33603  | 0.133975  | 0         | 1.14822   |  8.1537  |  8.70069 |  8.4584  | 10.0087  |          0.792159 |
| TCGA-FY-A4B0-01A | DM2                | 2.98861  | 1.25789   | 0.467962  | 1.95745   |  8.7535  |  9.02307 |  8.24455 | 10.1586  |          1.63314  |
| TCGA-FY-A40M-01A | DM1                | 2.20033  | 2.12544   | 0.0867803 | 0.535537  |  8.49742 |  9.95372 |  9.73302 | 10.0683  |          1.05638  |
| TCGA-ET-A4KN-01A | DM1                | 2.56264  | 3.38779   | 2.48183   | 3.04285   | 12.6946  | 12.7836  | 11.6991  | 11.1088  |          3.94966  |
| TCGA-ET-A40P-01A | DM1                | 2.51246  | 3.1774    | 0.750638  | 0.825659  |  9.99533 | 10.4813  |  9.71131 | 11.1298  |          2.72343  |
| TCGA-FE-A22Z-01A | DM1                | 3.635    | 0.118239  | 0.596892  | 0.513032  |  9.75841 |  9.45729 |  8.91194 | 10.892   |          1.00351  |
| TCGA-FK-A3SG-01A | DM1                | 3.33566  | 1.00208   | 0.531795  | 0.400938  | 10.4821  | 10.144   |  9.88231 | 11.8442  |          1.31869  |
| TCGA-EM-A22Q-01A | DM1                | 0.680788 | 0.680788  | 0.15009   | 2.1007    | 10.0904  |  9.63303 |  8.22585 |  8.44138 |          1.80933  |
| TCGA-DE-A4MC-01A | DM1                | 4.1209   | 0.988352  | 3.02069   | 2.42699   | 11.2968  | 11.5036  | 10.7122  | 13.0263  |          2.6435   |
| TCGA-DJ-A13S-01A | DM2                | 3.28956  | 1.78609   | 0.0381027 | 0.0934339 |  8.63326 |  9.80262 |  9.47191 | 11.2573  |          1.08683  |
| TCGA-EM-A3SY-01A | DM2                | 4.47498  | 0.779811  | 0.219452  | 0.876873  | 11.0291  | 10.9153  | 10.6943  | 13.6045  |          1.56291  |
| TCGA-EL-A3TB-01A | DM1                | 3.23903  | 1.58683   | 0.889119  | 2.16069   | 10.7454  | 10.9455  | 10.2707  | 12.1126  |          2.02936  |
| TCGA-BJ-A45C-01A | DM2                | 2.82482  | 0.587478  | 0         | 0.216058  |  9.32897 |  8.45422 |  8.42642 |  9.75937 |          1.29729  |
| TCGA-DJ-A2PX-01A | DM1                | 3.77181  | 3.25107   | 2.16807   | 3.71631   | 10.6272  | 10.9515  | 10.7004  | 12.3857  |          4.26099  |
| TCGA-EL-A3ZS-01A | DM1                | 3.32306  | 3.48878   | 3.24166   | 4.20417   | 10.5952  | 11.4584  | 10.97    | 12.5984  |          3.33537  |

### 8.4 Immune subtype distribution
- file: `/home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/A5_immune_subtype_distribution.tsv`
| immune_subtype   | status                           |
|:-----------------|:---------------------------------|
| missing          | Thorsson subtype table not local |

## SEC 9. A6 TF circuit + RAI prediction — RAW

### 9.1 TF activity
- file: `/home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/A6_tf_activity_per_sample.tsv`
| sample_id        | cluster   | tf   |   activity |
|:-----------------|:----------|:-----|-----------:|
| TCGA-FK-A3SE-01A | DM1       | PAX8 |    8.06643 |
| TCGA-FY-A2QD-01A | DM2       | PAX8 |   10.1899  |
| TCGA-EL-A4KD-01A | DM1       | PAX8 |    8.66641 |
| TCGA-BJ-A192-01A | DM2       | PAX8 |    8.94602 |
| TCGA-ET-A40S-01A | DM1       | PAX8 |    7.25262 |
| TCGA-ET-A39L-01A | DM1       | PAX8 |    8.45643 |
| TCGA-KS-A4IB-01A | DM1       | PAX8 |    6.29045 |
| TCGA-DJ-A2PO-01A | DM1       | PAX8 |    7.1391  |
| TCGA-EM-A2CJ-01A | DM2       | PAX8 |   10.0188  |
| TCGA-DJ-A1QL-01A | DM2       | PAX8 |    9.65814 |
| TCGA-EM-A1CS-01A | DM1       | PAX8 |    6.83195 |
| TCGA-CE-A484-01A | DM1       | PAX8 |    7.68702 |
| TCGA-EL-A3ZP-01A | DM1       | PAX8 |    7.51392 |
| TCGA-BJ-A0YZ-01A | DM2       | PAX8 |    7.7592  |
| TCGA-FY-A3NM-01A | DM2       | PAX8 |   10.8455  |
| TCGA-BJ-A3PR-01A | DM1       | PAX8 |    7.91498 |
| TCGA-EL-A3CO-01A | DM1       | PAX8 |    9.62434 |
| TCGA-CE-A482-01A | DM1       | PAX8 |    6.43279 |
| TCGA-EM-A3FQ-01A | DM1       | PAX8 |    8.51871 |
| TCGA-EL-A3GO-01A | DM2       | PAX8 |    8.87632 |
| TCGA-DJ-A3UZ-01A | DM1       | PAX8 |    8.77912 |
| TCGA-E8-A432-01A | DM1       | PAX8 |    5.87749 |
| TCGA-DJ-A4UT-01A | DM1       | PAX8 |    7.06489 |
| TCGA-ET-A40R-01A | DM1       | PAX8 |    9.26958 |
| TCGA-EM-A3FR-01A | DM2       | PAX8 |    9.12143 |
| TCGA-DJ-A4UQ-01A | DM1       | PAX8 |    7.92801 |
| TCGA-ET-A4KQ-01A | DM2       | PAX8 |    9.68026 |
| TCGA-CE-A481-01A | DM1       | PAX8 |    7.94368 |
| TCGA-EM-A2CP-01A | DM1       | PAX8 |    9.13855 |
| TCGA-EL-A3H3-01A | DM1       | PAX8 |    6.32722 |
| TCGA-EL-A3CS-01A | DM1       | PAX8 |    5.70783 |
| TCGA-FK-A3SD-01A | DM2       | PAX8 |   10.2028  |
| TCGA-DJ-A13W-01A | DM1       | PAX8 |    7.12219 |
| TCGA-DJ-A13R-01A | DM1       | PAX8 |    7.55411 |
| TCGA-EL-A3ZK-01A | DM1       | PAX8 |    7.00292 |
| TCGA-ET-A3DS-01A | DM2       | PAX8 |   11.0515  |
| TCGA-EL-A3D4-01A | DM1       | PAX8 |    8.0657  |
| TCGA-EM-A3AN-01A | DM1       | PAX8 |    7.48818 |
| TCGA-DE-A0Y2-01A | DM1       | PAX8 |    6.07517 |
| TCGA-CE-A27D-01A | DM1       | PAX8 |    9.27925 |
| TCGA-BJ-A0Z0-01A | DM2       | PAX8 |    9.86506 |
| TCGA-FK-A4UB-01A | DM1       | PAX8 |    7.02439 |
| TCGA-H2-A421-01A | DM1       | PAX8 |    6.81497 |
| TCGA-CE-A483-01A | DM1       | PAX8 |    6.87992 |
| TCGA-FY-A76V-01A | DM1       | PAX8 |    8.06675 |
| TCGA-DO-A2HM-01B | DM1       | PAX8 |    6.52478 |
| TCGA-EL-A3T0-01A | DM1       | PAX8 |    8.82307 |
| TCGA-BJ-A0ZJ-01A | DM1       | PAX8 |    8.62698 |
| TCGA-EL-A3T9-01A | DM1       | PAX8 |    4.78754 |
| TCGA-EM-A1YD-01A | DM2       | PAX8 |   10.4595  |

...

| sample_id        | cluster   | tf   |   activity |
|:-----------------|:----------|:-----|-----------:|
| TCGA-ET-A39T-01A | DM1       | HHEX |    9.78168 |
| TCGA-DE-A69K-01A | DM1       | HHEX |    9.77683 |
| TCGA-ET-A25M-01A | DM1       | HHEX |    9.74662 |
| TCGA-DJ-A4V5-01A | DM1       | HHEX |    9.88319 |
| TCGA-DE-A4M9-01A | DM2       | HHEX |   10.7291  |
| TCGA-ET-A39N-01A | DM2       | HHEX |    9.70663 |
| TCGA-E8-A416-01A | DM2       | HHEX |   10.7623  |
| TCGA-KS-A4ID-01A | DM1       | HHEX |    9.86756 |
| TCGA-DJ-A13U-01A | DM1       | HHEX |    9.92733 |
| TCGA-EL-A3T2-01A | DM2       | HHEX |   10.4996  |
| TCGA-DE-A4MA-01A | DM1       | HHEX |    8.95319 |
| TCGA-DJ-A2Q0-01A | DM2       | HHEX |    7.33406 |
| TCGA-EL-A3T6-01A | DM1       | HHEX |    9.66586 |
| TCGA-EM-A3O8-01A | DM2       | HHEX |   10.6361  |
| TCGA-DJ-A2Q2-01A | DM2       | HHEX |   10.462   |
| TCGA-BJ-A0ZG-01A | DM2       | HHEX |   10.7309  |
| TCGA-EM-A2OW-01A | DM2       | HHEX |    9.82088 |
| TCGA-BJ-A3F0-01A | DM2       | HHEX |   10.3183  |
| TCGA-EL-A3CY-01A | DM1       | HHEX |    8.673   |
| TCGA-EM-A3AQ-01A | DM1       | HHEX |    9.7309  |
| TCGA-J8-A4HW-06A | DM1       | HHEX |    9.8674  |
| TCGA-J8-A4HW-01A | DM1       | HHEX |   10.2931  |
| TCGA-BJ-A291-01A | DM2       | HHEX |    9.57979 |
| TCGA-CE-A3MD-01A | DM1       | HHEX |    9.98241 |
| TCGA-EL-A3ZN-01A | DM1       | HHEX |    9.79607 |
| TCGA-BJ-A2N9-01A | DM2       | HHEX |   10.1504  |
| TCGA-DJ-A4UP-01A | DM1       | HHEX |    9.95012 |
| TCGA-E8-A44M-01A | DM1       | HHEX |    9.84993 |
| TCGA-EM-A3O6-01A | DM2       | HHEX |   10.4144  |
| TCGA-EM-A4FH-01A | DM2       | HHEX |   10.8314  |
| TCGA-FY-A40N-01A | DM2       | HHEX |   10.5819  |
| TCGA-FY-A3WA-01A | DM2       | HHEX |   10.7303  |
| TCGA-ET-A2MX-01A | DM1       | HHEX |    8.8001  |
| TCGA-FY-A3I5-01B | DM2       | HHEX |   10.4098  |
| TCGA-EM-A2OV-01A | DM2       | HHEX |    9.96614 |
| TCGA-EM-A3AI-01A | DM2       | HHEX |   10.8707  |
| TCGA-FY-A4B0-01A | DM2       | HHEX |   10.3047  |
| TCGA-FY-A40M-01A | DM1       | HHEX |    8.91026 |
| TCGA-ET-A4KN-01A | DM1       | HHEX |    7.66665 |
| TCGA-ET-A40P-01A | DM1       | HHEX |    9.93503 |
| TCGA-FE-A22Z-01A | DM1       | HHEX |    9.45779 |
| TCGA-FK-A3SG-01A | DM1       | HHEX |    9.55418 |
| TCGA-EM-A22Q-01A | DM1       | HHEX |    9.94692 |
| TCGA-DE-A4MC-01A | DM1       | HHEX |    9.26422 |
| TCGA-DJ-A13S-01A | DM2       | HHEX |   10.178   |
| TCGA-EM-A3SY-01A | DM2       | HHEX |    7.96841 |
| TCGA-EL-A3TB-01A | DM1       | HHEX |   10.2835  |
| TCGA-BJ-A45C-01A | DM2       | HHEX |   10.7223  |
| TCGA-DJ-A2PX-01A | DM1       | HHEX |    9.83557 |
| TCGA-EL-A3ZS-01A | DM1       | HHEX |    9.84097 |

| cluster   | tf     |     mean |      std |
|:----------|:-------|---------:|---------:|
| DM1       | FOXE1  |  9.07092 | 1.26362  |
| DM1       | HHEX   |  9.62571 | 0.801288 |
| DM1       | NKX2-1 |  8.23412 | 1.08333  |
| DM1       | PAX8   |  7.30025 | 1.36705  |
| DM2       | FOXE1  | 10.7461  | 0.912325 |
| DM2       | HHEX   | 10.3138  | 0.765134 |
| DM2       | NKX2-1 |  9.36129 | 0.860003 |
| DM2       | PAX8   |  9.38622 | 0.916146 |
### 9.2 TF target overlap
- file: `/home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/A6_tf_target_overlap.tsv`
| tf     |   n_targets_present | targets                                   |
|:-------|--------------------:|:------------------------------------------|
| PAX8   |                   8 | TG,TPO,TSHR,SLC5A5,SLC26A4,DIO1,DIO2,IYD  |
| NKX2-1 |                   8 | TG,TPO,TSHR,SLC5A5,PAX8,FOXE1,DUOX1,DUOX2 |
| FOXE1  |                   6 | TG,TPO,PAX8,NKX2-1,SLC26A4,IYD            |
| HHEX   |                   5 | FOXE1,PAX8,NKX2-1,TG,TSHR                 |

### 9.3 RAI uptake score
- file: `/home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/A6_rai_uptake_score.tsv`
| sample_id        | cluster   |   rai_score |
|:-----------------|:----------|------------:|
| TCGA-FK-A3SE-01A | DM1       |     8.02363 |
| TCGA-FY-A2QD-01A | DM2       |     9.6369  |
| TCGA-EL-A4KD-01A | DM1       |     8.95363 |
| TCGA-BJ-A192-01A | DM2       |     9.01191 |
| TCGA-ET-A40S-01A | DM1       |     7.48305 |
| TCGA-ET-A39L-01A | DM1       |     8.43548 |
| TCGA-KS-A4IB-01A | DM1       |     6.32345 |
| TCGA-DJ-A2PO-01A | DM1       |     7.00901 |
| TCGA-EM-A2CJ-01A | DM2       |     9.76971 |
| TCGA-DJ-A1QL-01A | DM2       |     9.51776 |
| TCGA-EM-A1CS-01A | DM1       |     7.36269 |
| TCGA-CE-A484-01A | DM1       |     7.97868 |
| TCGA-EL-A3ZP-01A | DM1       |     8.0507  |
| TCGA-BJ-A0YZ-01A | DM2       |     8.21595 |
| TCGA-FY-A3NM-01A | DM2       |    10.4552  |
| TCGA-BJ-A3PR-01A | DM1       |     7.84422 |
| TCGA-EL-A3CO-01A | DM1       |     9.23233 |
| TCGA-CE-A482-01A | DM1       |     6.77373 |
| TCGA-EM-A3FQ-01A | DM1       |     8.68297 |
| TCGA-EL-A3GO-01A | DM2       |     8.6053  |
| TCGA-DJ-A3UZ-01A | DM1       |     8.85158 |
| TCGA-E8-A432-01A | DM1       |     6.48077 |
| TCGA-DJ-A4UT-01A | DM1       |     7.38698 |
| TCGA-ET-A40R-01A | DM1       |     9.6007  |
| TCGA-EM-A3FR-01A | DM2       |     8.91875 |
| TCGA-DJ-A4UQ-01A | DM1       |     7.99091 |
| TCGA-ET-A4KQ-01A | DM2       |     9.59345 |
| TCGA-CE-A481-01A | DM1       |     8.15195 |
| TCGA-EM-A2CP-01A | DM1       |     8.6882  |
| TCGA-EL-A3H3-01A | DM1       |     6.83075 |
| TCGA-EL-A3CS-01A | DM1       |     6.48032 |
| TCGA-FK-A3SD-01A | DM2       |     9.9413  |
| TCGA-DJ-A13W-01A | DM1       |     7.31335 |
| TCGA-DJ-A13R-01A | DM1       |     7.5934  |
| TCGA-EL-A3ZK-01A | DM1       |     7.04622 |
| TCGA-ET-A3DS-01A | DM2       |    10.5403  |
| TCGA-EL-A3D4-01A | DM1       |     8.50464 |
| TCGA-EM-A3AN-01A | DM1       |     7.72235 |
| TCGA-DE-A0Y2-01A | DM1       |     6.61181 |
| TCGA-CE-A27D-01A | DM1       |     9.30818 |
| TCGA-BJ-A0Z0-01A | DM2       |     9.54684 |
| TCGA-FK-A4UB-01A | DM1       |     7.43138 |
| TCGA-H2-A421-01A | DM1       |     7.09984 |
| TCGA-CE-A483-01A | DM1       |     6.72499 |
| TCGA-FY-A76V-01A | DM1       |     8.30574 |
| TCGA-DO-A2HM-01B | DM1       |     6.40334 |
| TCGA-EL-A3T0-01A | DM1       |     8.92073 |
| TCGA-BJ-A0ZJ-01A | DM1       |     8.4326  |
| TCGA-EL-A3T9-01A | DM1       |     5.9317  |
| TCGA-EM-A1YD-01A | DM2       |    10.1612  |

...

| sample_id        | cluster   |   rai_score |
|:-----------------|:----------|------------:|
| TCGA-ET-A39T-01A | DM1       |     7.12165 |
| TCGA-DE-A69K-01A | DM1       |     7.14371 |
| TCGA-ET-A25M-01A | DM1       |     6.91913 |
| TCGA-DJ-A4V5-01A | DM1       |     7.85746 |
| TCGA-DE-A4M9-01A | DM2       |     8.83342 |
| TCGA-ET-A39N-01A | DM2       |     8.33539 |
| TCGA-E8-A416-01A | DM2       |     9.20746 |
| TCGA-KS-A4ID-01A | DM1       |     8.62262 |
| TCGA-DJ-A13U-01A | DM1       |     7.42806 |
| TCGA-EL-A3T2-01A | DM2       |     9.03334 |
| TCGA-DE-A4MA-01A | DM1       |     6.10887 |
| TCGA-DJ-A2Q0-01A | DM2       |     7.38092 |
| TCGA-EL-A3T6-01A | DM1       |     6.8982  |
| TCGA-EM-A3O8-01A | DM2       |     9.02142 |
| TCGA-DJ-A2Q2-01A | DM2       |     9.69552 |
| TCGA-BJ-A0ZG-01A | DM2       |     9.53485 |
| TCGA-EM-A2OW-01A | DM2       |    10.0407  |
| TCGA-BJ-A3F0-01A | DM2       |    10.8516  |
| TCGA-EL-A3CY-01A | DM1       |     7.04248 |
| TCGA-EM-A3AQ-01A | DM1       |     7.03765 |
| TCGA-J8-A4HW-06A | DM1       |     7.52017 |
| TCGA-J8-A4HW-01A | DM1       |     8.31276 |
| TCGA-BJ-A291-01A | DM2       |     8.64342 |
| TCGA-CE-A3MD-01A | DM1       |     7.60137 |
| TCGA-EL-A3ZN-01A | DM1       |     7.52852 |
| TCGA-BJ-A2N9-01A | DM2       |     9.36015 |
| TCGA-DJ-A4UP-01A | DM1       |     8.48931 |
| TCGA-E8-A44M-01A | DM1       |     8.3633  |
| TCGA-EM-A3O6-01A | DM2       |     9.75449 |
| TCGA-EM-A4FH-01A | DM2       |     9.52686 |
| TCGA-FY-A40N-01A | DM2       |     9.31664 |
| TCGA-FY-A3WA-01A | DM2       |     9.52753 |
| TCGA-ET-A2MX-01A | DM1       |     6.22166 |
| TCGA-FY-A3I5-01B | DM2       |     9.5674  |
| TCGA-EM-A2OV-01A | DM2       |     9.38485 |
| TCGA-EM-A3AI-01A | DM2       |     9.79107 |
| TCGA-FY-A4B0-01A | DM2       |     9.35685 |
| TCGA-FY-A40M-01A | DM1       |     6.99326 |
| TCGA-ET-A4KN-01A | DM1       |     4.96019 |
| TCGA-ET-A40P-01A | DM1       |     8.28714 |
| TCGA-FE-A22Z-01A | DM1       |     5.91866 |
| TCGA-FK-A3SG-01A | DM1       |     7.62733 |
| TCGA-EM-A22Q-01A | DM1       |     8.81481 |
| TCGA-DE-A4MC-01A | DM1       |     6.05471 |
| TCGA-DJ-A13S-01A | DM2       |     8.79407 |
| TCGA-EM-A3SY-01A | DM2       |     7.64784 |
| TCGA-EL-A3TB-01A | DM1       |     8.17911 |
| TCGA-BJ-A45C-01A | DM2       |     8.8994  |
| TCGA-DJ-A2PX-01A | DM1       |     7.66879 |
| TCGA-EL-A3ZS-01A | DM1       |     8.45883 |

MISSING: requested columns `predicted_dm, NIS_expr, TPO_expr, TG_expr, TSHR_expr, PAX8_expr, NKX2-1_expr, FOXE1_expr, rai_score_normalized` not present in file
### 9.4 PDTC validation
- file: `/home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/A6_pdtc_rai_validation.tsv`
| sample_id   | dataset   | platform   | modality   | tissue_type                         | normal_vs_tumor   | histology_subtype   | driver_anchor   | molecular_subtype   | tds_group   | dediff_flag   | aggressive_flag   |   ajcc_stage_group |   ata_risk_proxy |   bethesda_category | clinical_subtype_tag   | label_confidence   | sex    |   age | outcome_available   | survival_available   | recurrence_available   | source_note                                                                                                    |   tds_score |   dedifferentiation_proxy_score |   brs_like_surrogate_score | v3_anchor       |   rai_score |   is_pdtc |
|:------------|:----------|:-----------|:-----------|:------------------------------------|:------------------|:--------------------|:----------------|:--------------------|:------------|:--------------|:------------------|-------------------:|-----------------:|--------------------:|:-----------------------|:-------------------|:-------|------:|:--------------------|:---------------------|:-----------------------|:---------------------------------------------------------------------------------------------------------------|------------:|--------------------------------:|---------------------------:|:----------------|------------:|----------:|
| GSM2024824  | GSE76039  | GPL570     | microarray | Anaplastic thyroid tumor            | tumor             | ATC                 | unknown         | dedifferentiated    | mid         | yes           | yes               |                nan |              nan |                 nan | s_JF_thy_001_P         | high               | Female |   nan | no                  | no                   | no                     | {"gender": "female", "tissue": "Thyroid", "tumor type": "Primary"}                                             |  -1.6492    |                       1.6492    |                        nan | histology_proxy |     5.55064 |         0 |
| GSM2024825  | GSE76039  | GPL570     | microarray | Anaplastic thyroid tumor            | tumor             | ATC                 | unknown         | dedifferentiated    | low         | yes           | yes               |                nan |              nan |                 nan | s_JF_thy_002_P         | high               | Female |   nan | no                  | no                   | no                     | {"gender": "female", "tissue": "Thyroid", "tumor type": "Primary"}                                             |  -2.30175   |                       2.30175   |                        nan | histology_proxy |     4.86282 |         0 |
| GSM2024826  | GSE76039  | GPL570     | microarray | Anaplastic thyroid tumor            | tumor             | ATC                 | unknown         | dedifferentiated    | low         | yes           | yes               |                nan |              nan |                 nan | s_JF_thy_003_P         | high               | Male   |   nan | no                  | no                   | no                     | {"gender": "male", "tissue": "Thyroid", "tumor type": "Primary"}                                               |  -2.44153   |                       2.44153   |                        nan | histology_proxy |     5.28801 |         0 |
| GSM2024827  | GSE76039  | GPL570     | microarray | Poorly-differentiated thyroid tumor | tumor             | PDTC                | unknown         | dedifferentiated    | high        | yes           | yes               |                nan |              nan |                 nan | s_JF_thy_004_P         | high               | Female |   nan | no                  | no                   | no                     | {"gender": "female", "tissue": "Thyroid", "tumor type": "Primary"}                                             |   2.61909   |                      -2.61909   |                        nan | histology_proxy |    11.1485  |         1 |
| GSM2024828  | GSE76039  | GPL570     | microarray | Poorly-differentiated thyroid tumor | tumor             | PDTC                | unknown         | dedifferentiated    | high        | yes           | yes               |                nan |              nan |                 nan | s_JF_thy_005_P         | high               | Male   |   nan | no                  | no                   | no                     | {"gender": "male", "tissue": "Thyroid", "tumor type": "Primary"}                                               |   2.67015   |                      -2.67015   |                        nan | histology_proxy |    11.0833  |         1 |
| GSM2024829  | GSE76039  | GPL570     | microarray | Poorly-differentiated thyroid tumor | tumor             | PDTC                | unknown         | dedifferentiated    | high        | yes           | yes               |                nan |              nan |                 nan | s_JF_thy_006_P         | high               | Female |   nan | no                  | no                   | no                     | {"gender": "female", "tissue": "Thyroid", "tumor type": "Primary"}                                             |   2.25557   |                      -2.25557   |                        nan | histology_proxy |    11.0072  |         1 |
| GSM2024830  | GSE76039  | GPL570     | microarray | Poorly-differentiated thyroid tumor | tumor             | PDTC                | unknown         | dedifferentiated    | mid         | yes           | yes               |                nan |              nan |                 nan | s_JF_thy_007_P         | high               | Female |   nan | no                  | no                   | no                     | {"gender": "female", "tissue": "Thyroid", "tumor type": "Primary"}                                             |   1.75655   |                      -1.75655   |                        nan | histology_proxy |     9.58108 |         1 |
| GSM2024831  | GSE76039  | GPL570     | microarray | Anaplastic thyroid tumor            | tumor             | ATC                 | unknown         | dedifferentiated    | low         | yes           | yes               |                nan |              nan |                 nan | s_JF_thy_008_P         | high               | Female |   nan | no                  | no                   | no                     | {"gender": "female", "tissue": "Thyroid", "tumor type": "Primary"}                                             |  -3.22468   |                       3.22468   |                        nan | histology_proxy |     4.06814 |         0 |
| GSM2024832  | GSE76039  | GPL570     | microarray | Poorly-differentiated thyroid tumor | tumor             | PDTC                | unknown         | dedifferentiated    | high        | yes           | yes               |                nan |              nan |                 nan | s_JF_thy_009_P         | high               | Female |   nan | no                  | no                   | no                     | {"gender": "female", "tissue": "Thyroid", "tumor type": "Primary"}                                             |   2.93594   |                      -2.93594   |                        nan | histology_proxy |    11.5096  |         1 |
| GSM2024833  | GSE76039  | GPL570     | microarray | Poorly-differentiated thyroid tumor | tumor             | PDTC                | unknown         | dedifferentiated    | mid         | yes           | yes               |                nan |              nan |                 nan | s_JF_thy_010_P         | high               | Female |   nan | no                  | no                   | no                     | {"gender": "female", "tissue": "Thyroid", "tumor type": "Recurrent tumor in neck"}                             |   0.711825  |                      -0.711825  |                        nan | histology_proxy |     8.3416  |         1 |
| GSM2024834  | GSE76039  | GPL570     | microarray | Anaplastic thyroid tumor            | tumor             | ATC                 | unknown         | dedifferentiated    | low         | yes           | yes               |                nan |              nan |                 nan | s_JF_thy_011_P         | high               | Male   |   nan | no                  | no                   | no                     | {"gender": "male", "tissue": "Thyroid", "tumor type": "Primary"}                                               |  -2.13233   |                       2.13233   |                        nan | histology_proxy |     5.41818 |         0 |
| GSM2024835  | GSE76039  | GPL570     | microarray | Poorly-differentiated thyroid tumor | tumor             | PDTC                | unknown         | dedifferentiated    | mid         | yes           | yes               |                nan |              nan |                 nan | s_JF_thy_012_P         | high               | Female |   nan | no                  | no                   | no                     | {"gender": "female", "tissue": "Thyroid", "tumor type": "Recurrent/persistent metastasic tumor to lymph node"} |   0.327024  |                      -0.327024  |                        nan | histology_proxy |     7.98736 |         1 |
| GSM2024836  | GSE76039  | GPL570     | microarray | Anaplastic thyroid tumor            | tumor             | ATC                 | unknown         | dedifferentiated    | low         | yes           | yes               |                nan |              nan |                 nan | s_JF_thy_013_P         | high               | Female |   nan | no                  | no                   | no                     | {"gender": "female", "tissue": "Thyroid", "tumor type": "Primary"}                                             |  -2.80435   |                       2.80435   |                        nan | histology_proxy |     4.53462 |         0 |
| GSM2024837  | GSE76039  | GPL570     | microarray | Anaplastic thyroid tumor            | tumor             | ATC                 | unknown         | dedifferentiated    | mid         | yes           | yes               |                nan |              nan |                 nan | s_JF_thy_014_P         | high               | Female |   nan | no                  | no                   | no                     | {"gender": "female", "tissue": "Thyroid", "tumor type": "Metastasis"}                                          |   0.0661645 |                      -0.0661645 |                        nan | histology_proxy |     7.55865 |         0 |
| GSM2024838  | GSE76039  | GPL570     | microarray | Anaplastic thyroid tumor            | tumor             | ATC                 | unknown         | dedifferentiated    | low         | yes           | yes               |                nan |              nan |                 nan | s_JF_thy_015_P         | high               | Female |   nan | no                  | no                   | no                     | {"gender": "female", "tissue": "Thyroid", "tumor type": "Primary"}                                             |  -3.29458   |                       3.29458   |                        nan | histology_proxy |     3.9969  |         0 |
| GSM2024839  | GSE76039  | GPL570     | microarray | Anaplastic thyroid tumor            | tumor             | ATC                 | unknown         | dedifferentiated    | mid         | yes           | yes               |                nan |              nan |                 nan | s_JF_thy_016_P         | high               | Male   |   nan | no                  | no                   | no                     | {"gender": "male", "tissue": "Thyroid", "tumor type": "Primary"}                                               |  -1.5201    |                       1.5201    |                        nan | histology_proxy |     6.57181 |         0 |
| GSM2024840  | GSE76039  | GPL570     | microarray | Anaplastic thyroid tumor            | tumor             | ATC                 | unknown         | dedifferentiated    | low         | yes           | yes               |                nan |              nan |                 nan | s_JF_thy_017_P         | high               | Male   |   nan | no                  | no                   | no                     | {"gender": "male", "tissue": "Thyroid", "tumor type": "Primary"}                                               |  -3.23456   |                       3.23456   |                        nan | histology_proxy |     4.2285  |         0 |
| GSM2024841  | GSE76039  | GPL570     | microarray | Poorly-differentiated thyroid tumor | tumor             | PDTC                | unknown         | dedifferentiated    | high        | yes           | yes               |                nan |              nan |                 nan | s_JF_thy_018_P         | high               | Female |   nan | no                  | no                   | no                     | {"gender": "female", "tissue": "Thyroid", "tumor type": "Recurrent tumor"}                                     |   2.1894    |                      -2.1894    |                        nan | histology_proxy |    10.7246  |         1 |
| GSM2024842  | GSE76039  | GPL570     | microarray | Poorly-differentiated thyroid tumor | tumor             | PDTC                | unknown         | dedifferentiated    | high        | yes           | yes               |                nan |              nan |                 nan | s_JF_thy_019_P         | high               | Female |   nan | no                  | no                   | no                     | {"gender": "female", "tissue": "Thyroid", "tumor type": "Primary"}                                             |   2.65314   |                      -2.65314   |                        nan | histology_proxy |    11.2276  |         1 |
| GSM2024843  | GSE76039  | GPL570     | microarray | Poorly-differentiated thyroid tumor | tumor             | PDTC                | unknown         | dedifferentiated    | high        | yes           | yes               |                nan |              nan |                 nan | s_JF_thy_020_P         | high               | Female |   nan | no                  | no                   | no                     | {"gender": "female", "tissue": "Thyroid", "tumor type": "Recurrent tumor in neck"}                             |   2.54389   |                      -2.54389   |                        nan | histology_proxy |    11.1946  |         1 |
| GSM2024844  | GSE76039  | GPL570     | microarray | Poorly-differentiated thyroid tumor | tumor             | PDTC                | unknown         | dedifferentiated    | high        | yes           | yes               |                nan |              nan |                 nan | s_JF_thy_021_P         | high               | Female |   nan | no                  | no                   | no                     | {"gender": "female", "tissue": "Thyroid", "tumor type": "Primary"}                                             |   2.19793   |                      -2.19793   |                        nan | histology_proxy |    10.6374  |         1 |
| GSM2024845  | GSE76039  | GPL570     | microarray | Anaplastic thyroid tumor            | tumor             | ATC                 | unknown         | dedifferentiated    | mid         | yes           | yes               |                nan |              nan |                 nan | s_JF_thy_022_P         | high               | Female |   nan | no                  | no                   | no                     | {"gender": "female", "tissue": "Thyroid", "tumor type": "Primary"}                                             |  -1.81797   |                       1.81797   |                        nan | histology_proxy |     6.40236 |         0 |
| GSM2024846  | GSE76039  | GPL570     | microarray | Anaplastic thyroid tumor            | tumor             | ATC                 | unknown         | dedifferentiated    | mid         | yes           | yes               |                nan |              nan |                 nan | s_JF_thy_023_P         | high               | Female |   nan | no                  | no                   | no                     | {"gender": "female", "tissue": "Thyroid", "tumor type": "Primary"}                                             |   1.62615   |                      -1.62615   |                        nan | histology_proxy |    10.3169  |         0 |
| GSM2024847  | GSE76039  | GPL570     | microarray | Anaplastic thyroid tumor            | tumor             | ATC                 | unknown         | dedifferentiated    | mid         | yes           | yes               |                nan |              nan |                 nan | s_JF_thy_024_P         | high               | Male   |   nan | no                  | no                   | no                     | {"gender": "male", "tissue": "Thyroid", "tumor type": "Primary"}                                               |  -1.11901   |                       1.11901   |                        nan | histology_proxy |     6.91571 |         0 |
| GSM2024848  | GSE76039  | GPL570     | microarray | Anaplastic thyroid tumor            | tumor             | ATC                 | unknown         | dedifferentiated    | low         | yes           | yes               |                nan |              nan |                 nan | s_JF_thy_025_P         | high               | Female |   nan | no                  | no                   | no                     | {"gender": "female", "tissue": "Thyroid", "tumor type": "Primary"}                                             |  -2.31208   |                       2.31208   |                        nan | histology_proxy |     4.82584 |         0 |
| GSM2024849  | GSE76039  | GPL570     | microarray | Poorly-differentiated thyroid tumor | tumor             | PDTC                | unknown         | dedifferentiated    | mid         | yes           | yes               |                nan |              nan |                 nan | s_JF_thy_026_P         | high               | Female |   nan | no                  | no                   | no                     | {"gender": "female", "tissue": "Thyroid", "tumor type": "Primary"}                                             |   0.530719  |                      -0.530719  |                        nan | histology_proxy |     8.5666  |         1 |
| GSM2024850  | GSE76039  | GPL570     | microarray | Poorly-differentiated thyroid tumor | tumor             | PDTC                | unknown         | dedifferentiated    | high        | yes           | yes               |                nan |              nan |                 nan | s_JF_thy_027_P         | high               | Female |   nan | no                  | no                   | no                     | {"gender": "female", "tissue": "Thyroid", "tumor type": "Primary"}                                             |   2.59107   |                      -2.59107   |                        nan | histology_proxy |    10.9828  |         1 |
| GSM2024851  | GSE76039  | GPL570     | microarray | Anaplastic thyroid tumor            | tumor             | ATC                 | unknown         | dedifferentiated    | low         | yes           | yes               |                nan |              nan |                 nan | s_JF_thy_028_P         | high               | Male   |   nan | no                  | no                   | no                     | {"gender": "male", "tissue": "Thyroid", "tumor type": "Primary"}                                               |  -3.48775   |                       3.48775   |                        nan | histology_proxy |     3.58334 |         0 |
| GSM2024852  | GSE76039  | GPL570     | microarray | Anaplastic thyroid tumor            | tumor             | ATC                 | unknown         | dedifferentiated    | low         | yes           | yes               |                nan |              nan |                 nan | s_JF_thy_029_P         | high               | Female |   nan | no                  | no                   | no                     | {"gender": "female", "tissue": "Thyroid", "tumor type": "Primary"}                                             |  -2.36262   |                       2.36262   |                        nan | histology_proxy |     4.73931 |         0 |
| GSM2024853  | GSE76039  | GPL570     | microarray | Anaplastic thyroid tumor            | tumor             | ATC                 | unknown         | dedifferentiated    | low         | yes           | yes               |                nan |              nan |                 nan | s_JF_thy_030_P         | high               | Female |   nan | no                  | no                   | no                     | {"gender": "female", "tissue": "Thyroid", "tumor type": "Primary"}                                             |  -2.10583   |                       2.10583   |                        nan | histology_proxy |     4.94083 |         0 |
| GSM2024854  | GSE76039  | GPL570     | microarray | Poorly-differentiated thyroid tumor | tumor             | PDTC                | unknown         | dedifferentiated    | high        | yes           | yes               |                nan |              nan |                 nan | s_JF_thy_031_P         | high               | Female |   nan | no                  | no                   | no                     | {"gender": "female", "tissue": "Thyroid", "tumor type": "Primary"}                                             |   2.82367   |                      -2.82367   |                        nan | histology_proxy |    11.3616  |         1 |
| GSM2024855  | GSE76039  | GPL570     | microarray | Poorly-differentiated thyroid tumor | tumor             | PDTC                | unknown         | dedifferentiated    | mid         | yes           | yes               |                nan |              nan |                 nan | s_JF_thy_032_P         | high               | Male   |   nan | no                  | no                   | no                     | {"gender": "male", "tissue": "Thyroid", "tumor type": "Primary"}                                               |   2.09792   |                      -2.09792   |                        nan | histology_proxy |    10.4408  |         1 |
| GSM2024856  | GSE76039  | GPL570     | microarray | Poorly-differentiated thyroid tumor | tumor             | PDTC                | unknown         | dedifferentiated    | high        | yes           | yes               |                nan |              nan |                 nan | s_JF_thy_033_P         | high               | Female |   nan | no                  | no                   | no                     | {"gender": "female", "tissue": "Thyroid", "tumor type": "Primary (residual)"}                                  |   2.88438   |                      -2.88438   |                        nan | histology_proxy |    11.6519  |         1 |
| GSM2024857  | GSE76039  | GPL570     | microarray | Anaplastic thyroid tumor            | tumor             | ATC                 | unknown         | dedifferentiated    | low         | yes           | yes               |                nan |              nan |                 nan | s_JF_thy_034_P         | high               | Female |   nan | no                  | no                   | no                     | {"gender": "female", "tissue": "Thyroid", "tumor type": "Primary"}                                             |  -3.92721   |                       3.92721   |                        nan | histology_proxy |     2.89789 |         0 |
| GSM2024858  | GSE76039  | GPL570     | microarray | Anaplastic thyroid tumor            | tumor             | ATC                 | unknown         | dedifferentiated    | low         | yes           | yes               |                nan |              nan |                 nan | s_JF_thy_035_P         | high               | Male   |   nan | no                  | no                   | no                     | {"gender": "male", "tissue": "Thyroid", "tumor type": "Recurrent tumor in neck"}                               |  -2.19222   |                       2.19222   |                        nan | histology_proxy |     5.0224  |         0 |
| GSM2024859  | GSE76039  | GPL570     | microarray | Poorly-differentiated thyroid tumor | tumor             | PDTC                | unknown         | dedifferentiated    | high        | yes           | yes               |                nan |              nan |                 nan | s_JF_thy_036_P         | high               | Male   |   nan | no                  | no                   | no                     | {"gender": "male", "tissue": "Thyroid", "tumor type": "Primary"}                                               |   2.56737   |                      -2.56737   |                        nan | histology_proxy |    10.8574  |         1 |
| GSM2024860  | GSE76039  | GPL570     | microarray | Anaplastic thyroid tumor            | tumor             | ATC                 | unknown         | dedifferentiated    | mid         | yes           | yes               |                nan |              nan |                 nan | s_JF_thy_037_P         | high               | Female |   nan | no                  | no                   | no                     | {"gender": "female", "tissue": "Thyroid", "tumor type": "Recurrent tumor in neck"}                             |  -0.388749  |                       0.388749  |                        nan | histology_proxy |     7.6173  |         0 |

MISSING: standalone ROC threshold table with sensitivity/specificity not found under project/results/v17p3/tables
## SEC 10. Cross-Phase 일관성 raw check

### 10.1 Sample membership 일관성
- Phase 1 sample list file: `/home/seungho/personal/THCA_data_analysis/project/results/v17/tables/dark_matter_cohort.tsv`
| sample_id        |
|:-----------------|
| TCGA-FK-A3SE-01A |
| TCGA-FY-A2QD-01A |
| TCGA-EL-A4KD-01A |
| TCGA-BJ-A192-01A |
| TCGA-ET-A40S-01A |
| TCGA-ET-A39L-01A |
| TCGA-KS-A4IB-01A |
| TCGA-DJ-A2PO-01A |
| TCGA-EM-A2CJ-01A |
| TCGA-DJ-A1QL-01A |
| TCGA-EM-A1CS-01A |
| TCGA-CE-A484-01A |
| TCGA-EL-A3ZP-01A |
| TCGA-BJ-A0YZ-01A |
| TCGA-FY-A3NM-01A |
| TCGA-BJ-A3PR-01A |
| TCGA-EL-A3CO-01A |
| TCGA-CE-A482-01A |
| TCGA-EM-A3FQ-01A |
| TCGA-EL-A3GO-01A |
| TCGA-DJ-A3UZ-01A |
| TCGA-E8-A432-01A |
| TCGA-DJ-A4UT-01A |
| TCGA-ET-A40R-01A |
| TCGA-EM-A3FR-01A |
| TCGA-DJ-A4UQ-01A |
| TCGA-ET-A4KQ-01A |
| TCGA-CE-A481-01A |
| TCGA-EM-A2CP-01A |
| TCGA-EL-A3H3-01A |
| TCGA-EL-A3CS-01A |
| TCGA-FK-A3SD-01A |
| TCGA-DJ-A13W-01A |
| TCGA-DJ-A13R-01A |
| TCGA-EL-A3ZK-01A |
| TCGA-ET-A3DS-01A |
| TCGA-EL-A3D4-01A |
| TCGA-EM-A3AN-01A |
| TCGA-DE-A0Y2-01A |
| TCGA-CE-A27D-01A |
| TCGA-BJ-A0Z0-01A |
| TCGA-FK-A4UB-01A |
| TCGA-H2-A421-01A |
| TCGA-CE-A483-01A |
| TCGA-FY-A76V-01A |
| TCGA-DO-A2HM-01B |
| TCGA-EL-A3T0-01A |
| TCGA-BJ-A0ZJ-01A |
| TCGA-EL-A3T9-01A |
| TCGA-EM-A1YD-01A |

...

| sample_id        |
|:-----------------|
| TCGA-ET-A39T-01A |
| TCGA-DE-A69K-01A |
| TCGA-ET-A25M-01A |
| TCGA-DJ-A4V5-01A |
| TCGA-DE-A4M9-01A |
| TCGA-ET-A39N-01A |
| TCGA-E8-A416-01A |
| TCGA-KS-A4ID-01A |
| TCGA-DJ-A13U-01A |
| TCGA-EL-A3T2-01A |
| TCGA-DE-A4MA-01A |
| TCGA-DJ-A2Q0-01A |
| TCGA-EL-A3T6-01A |
| TCGA-EM-A3O8-01A |
| TCGA-DJ-A2Q2-01A |
| TCGA-BJ-A0ZG-01A |
| TCGA-EM-A2OW-01A |
| TCGA-BJ-A3F0-01A |
| TCGA-EL-A3CY-01A |
| TCGA-EM-A3AQ-01A |
| TCGA-J8-A4HW-06A |
| TCGA-J8-A4HW-01A |
| TCGA-BJ-A291-01A |
| TCGA-CE-A3MD-01A |
| TCGA-EL-A3ZN-01A |
| TCGA-BJ-A2N9-01A |
| TCGA-DJ-A4UP-01A |
| TCGA-E8-A44M-01A |
| TCGA-EM-A3O6-01A |
| TCGA-EM-A4FH-01A |
| TCGA-FY-A40N-01A |
| TCGA-FY-A3WA-01A |
| TCGA-ET-A2MX-01A |
| TCGA-FY-A3I5-01B |
| TCGA-EM-A2OV-01A |
| TCGA-EM-A3AI-01A |
| TCGA-FY-A4B0-01A |
| TCGA-FY-A40M-01A |
| TCGA-ET-A4KN-01A |
| TCGA-ET-A40P-01A |
| TCGA-FE-A22Z-01A |
| TCGA-FK-A3SG-01A |
| TCGA-EM-A22Q-01A |
| TCGA-DE-A4MC-01A |
| TCGA-DJ-A13S-01A |
| TCGA-EM-A3SY-01A |
| TCGA-EL-A3TB-01A |
| TCGA-BJ-A45C-01A |
| TCGA-DJ-A2PX-01A |
| TCGA-EL-A3ZS-01A |

- Phase 2 robustness input sample_id list
MISSING: no standalone Phase 2 input sample list file found under `/home/seungho/personal/THCA_data_analysis/project/results/v17p2/tables`
- Phase 3 input sample_id lists
A2 full cohort sample_id head/tail
| sample_id        |
|:-----------------|
| TCGA-DJ-A2Q6-01A |
| TCGA-FK-A3SE-01A |
| TCGA-DJ-A2QA-01A |
| TCGA-FY-A2QD-01A |
| TCGA-EL-A3GR-01A |
| TCGA-EL-A3T7-01A |
| TCGA-FE-A230-01A |
| TCGA-EM-A3FO-01A |
| TCGA-DJ-A3VB-01A |
| TCGA-BJ-A28X-01A |
| TCGA-DJ-A2QB-01A |
| TCGA-DO-A1K0-01A |
| TCGA-EL-A4KD-01A |
| TCGA-BJ-A192-01A |
| TCGA-L6-A4EU-01A |
| TCGA-DJ-A3V0-01A |
| TCGA-EM-A2P2-01A |
| TCGA-BJ-A45G-01A |
| TCGA-FY-A3R7-01A |
| TCGA-DJ-A3VG-01A |
| TCGA-ET-A25O-01A |
| TCGA-DJ-A3UR-01A |
| TCGA-E8-A415-01A |
| TCGA-DJ-A4UR-01A |
| TCGA-DJ-A13L-01A |
| TCGA-BJ-A45I-01A |
| TCGA-ET-A40S-01A |
| TCGA-EL-A3CZ-01A |
| TCGA-EM-A2P3-01A |
| TCGA-BJ-A2N7-01A |
| TCGA-EL-A3GW-01A |
| TCGA-DJ-A13P-01A |
| TCGA-ET-A39O-01A |
| TCGA-IM-A3ED-01A |
| TCGA-ET-A39L-01A |
| TCGA-FY-A3YR-01A |
| TCGA-KS-A4IB-01A |
| TCGA-EM-A3ST-01A |
| TCGA-E8-A413-01A |
| TCGA-EL-A3TA-01A |
| TCGA-DJ-A2PW-01A |
| TCGA-DJ-A2PO-01A |
| TCGA-ET-A3BT-01A |
| TCGA-EM-A2CJ-01A |
| TCGA-DJ-A3VD-01A |
| TCGA-DJ-A1QE-01A |
| TCGA-BJ-A0ZB-01A |
| TCGA-IM-A3U3-01A |
| TCGA-BJ-A45F-01A |
| TCGA-EL-A4KG-01A |

...

| sample_id        |
|:-----------------|
| TCGA-EL-A3GQ-01A |
| TCGA-J8-A3YD-01A |
| TCGA-EM-A4FQ-01A |
| TCGA-FY-A3WA-01A |
| TCGA-ET-A2MX-01A |
| TCGA-FY-A3I5-01B |
| TCGA-IM-A41Y-01A |
| TCGA-EM-A2OV-01A |
| TCGA-E3-A3E3-01A |
| TCGA-ET-A3BQ-01B |
| TCGA-J8-A3O0-01A |
| TCGA-EL-A3CL-01A |
| TCGA-EL-A3T3-01A |
| TCGA-ET-A25J-01A |
| TCGA-EL-A3ZT-01A |
| TCGA-EM-A3AI-01A |
| TCGA-EM-A22O-01A |
| TCGA-ET-A40Q-01A |
| TCGA-FY-A4B0-01A |
| TCGA-ET-A39K-01A |
| TCGA-EM-A2CQ-01A |
| TCGA-EM-A2CS-01A |
| TCGA-EM-A2CS-06A |
| TCGA-FY-A40M-01A |
| TCGA-ET-A4KN-01A |
| TCGA-KS-A41I-01A |
| TCGA-DJ-A4UL-01A |
| TCGA-ET-A40P-01A |
| TCGA-DJ-A3V5-01A |
| TCGA-FE-A22Z-01A |
| TCGA-FE-A237-01A |
| TCGA-FK-A3SG-01A |
| TCGA-GE-A2C6-01A |
| TCGA-EM-A22Q-01A |
| TCGA-FY-A3W9-01A |
| TCGA-DE-A4MC-01A |
| TCGA-EL-A4K9-01A |
| TCGA-FK-A3SB-01A |
| TCGA-H2-A422-01A |
| TCGA-EM-A3FM-01A |
| TCGA-DJ-A13S-01A |
| TCGA-EM-A3SY-01A |
| TCGA-EM-A2P1-01A |
| TCGA-EL-A3T1-01A |
| TCGA-DJ-A2PS-01A |
| TCGA-EL-A3TB-01A |
| TCGA-BJ-A45C-01A |
| TCGA-DJ-A3UX-01A |
| TCGA-DJ-A2PX-01A |
| TCGA-EL-A3ZS-01A |

A5 dark cohort sample_id head/tail
| sample_id        |
|:-----------------|
| TCGA-FK-A3SE-01A |
| TCGA-FY-A2QD-01A |
| TCGA-EL-A4KD-01A |
| TCGA-BJ-A192-01A |
| TCGA-ET-A40S-01A |
| TCGA-ET-A39L-01A |
| TCGA-KS-A4IB-01A |
| TCGA-DJ-A2PO-01A |
| TCGA-EM-A2CJ-01A |
| TCGA-DJ-A1QL-01A |
| TCGA-EM-A1CS-01A |
| TCGA-CE-A484-01A |
| TCGA-EL-A3ZP-01A |
| TCGA-BJ-A0YZ-01A |
| TCGA-FY-A3NM-01A |
| TCGA-BJ-A3PR-01A |
| TCGA-EL-A3CO-01A |
| TCGA-CE-A482-01A |
| TCGA-EM-A3FQ-01A |
| TCGA-EL-A3GO-01A |
| TCGA-DJ-A3UZ-01A |
| TCGA-E8-A432-01A |
| TCGA-DJ-A4UT-01A |
| TCGA-ET-A40R-01A |
| TCGA-EM-A3FR-01A |
| TCGA-DJ-A4UQ-01A |
| TCGA-ET-A4KQ-01A |
| TCGA-CE-A481-01A |
| TCGA-EM-A2CP-01A |
| TCGA-EL-A3H3-01A |
| TCGA-EL-A3CS-01A |
| TCGA-FK-A3SD-01A |
| TCGA-DJ-A13W-01A |
| TCGA-DJ-A13R-01A |
| TCGA-EL-A3ZK-01A |
| TCGA-ET-A3DS-01A |
| TCGA-EL-A3D4-01A |
| TCGA-EM-A3AN-01A |
| TCGA-DE-A0Y2-01A |
| TCGA-CE-A27D-01A |
| TCGA-BJ-A0Z0-01A |
| TCGA-FK-A4UB-01A |
| TCGA-H2-A421-01A |
| TCGA-CE-A483-01A |
| TCGA-FY-A76V-01A |
| TCGA-DO-A2HM-01B |
| TCGA-EL-A3T0-01A |
| TCGA-BJ-A0ZJ-01A |
| TCGA-EL-A3T9-01A |
| TCGA-EM-A1YD-01A |

...

| sample_id        |
|:-----------------|
| TCGA-ET-A39T-01A |
| TCGA-DE-A69K-01A |
| TCGA-ET-A25M-01A |
| TCGA-DJ-A4V5-01A |
| TCGA-DE-A4M9-01A |
| TCGA-ET-A39N-01A |
| TCGA-E8-A416-01A |
| TCGA-KS-A4ID-01A |
| TCGA-DJ-A13U-01A |
| TCGA-EL-A3T2-01A |
| TCGA-DE-A4MA-01A |
| TCGA-DJ-A2Q0-01A |
| TCGA-EL-A3T6-01A |
| TCGA-EM-A3O8-01A |
| TCGA-DJ-A2Q2-01A |
| TCGA-BJ-A0ZG-01A |
| TCGA-EM-A2OW-01A |
| TCGA-BJ-A3F0-01A |
| TCGA-EL-A3CY-01A |
| TCGA-EM-A3AQ-01A |
| TCGA-J8-A4HW-06A |
| TCGA-J8-A4HW-01A |
| TCGA-BJ-A291-01A |
| TCGA-CE-A3MD-01A |
| TCGA-EL-A3ZN-01A |
| TCGA-BJ-A2N9-01A |
| TCGA-DJ-A4UP-01A |
| TCGA-E8-A44M-01A |
| TCGA-EM-A3O6-01A |
| TCGA-EM-A4FH-01A |
| TCGA-FY-A40N-01A |
| TCGA-FY-A3WA-01A |
| TCGA-ET-A2MX-01A |
| TCGA-FY-A3I5-01B |
| TCGA-EM-A2OV-01A |
| TCGA-EM-A3AI-01A |
| TCGA-FY-A4B0-01A |
| TCGA-FY-A40M-01A |
| TCGA-ET-A4KN-01A |
| TCGA-ET-A40P-01A |
| TCGA-FE-A22Z-01A |
| TCGA-FK-A3SG-01A |
| TCGA-EM-A22Q-01A |
| TCGA-DE-A4MC-01A |
| TCGA-DJ-A13S-01A |
| TCGA-EM-A3SY-01A |
| TCGA-EL-A3TB-01A |
| TCGA-BJ-A45C-01A |
| TCGA-DJ-A2PX-01A |
| TCGA-EL-A3ZS-01A |

| set_name            |   n |
|:--------------------|----:|
| phase1_dark         | 178 |
| phase3_A5_tmb       | 178 |
| intersection        | 178 |
| phase1_minus_phase3 |   0 |
| phase3_minus_phase1 |   0 |
- phase1_minus_phase3 sample_id list
| sample_id   |
|-------------|

- phase3_minus_phase1 sample_id list
| sample_id   |
|-------------|

### 10.2 DM1/DM2 라벨 일관성
| sample_id        | phase1_cluster   | sample_master_v17_full_cluster   |
|:-----------------|:-----------------|:---------------------------------|
| TCGA-FK-A3SE-01A | DM1              | DM1                              |
| TCGA-FY-A2QD-01A | DM2              | DM2                              |
| TCGA-EL-A4KD-01A | DM1              | DM1                              |
| TCGA-BJ-A192-01A | DM2              | DM2                              |
| TCGA-ET-A40S-01A | DM1              | DM1                              |
| TCGA-ET-A39L-01A | DM1              | DM1                              |
| TCGA-KS-A4IB-01A | DM1              | DM1                              |
| TCGA-DJ-A2PO-01A | DM1              | DM1                              |
| TCGA-EM-A2CJ-01A | DM2              | DM2                              |
| TCGA-DJ-A1QL-01A | DM2              | DM2                              |
| TCGA-EM-A1CS-01A | DM1              | DM1                              |
| TCGA-CE-A484-01A | DM1              | DM1                              |
| TCGA-EL-A3ZP-01A | DM1              | DM1                              |
| TCGA-BJ-A0YZ-01A | DM2              | DM2                              |
| TCGA-FY-A3NM-01A | DM2              | DM2                              |
| TCGA-BJ-A3PR-01A | DM1              | DM1                              |
| TCGA-EL-A3CO-01A | DM1              | DM1                              |
| TCGA-CE-A482-01A | DM1              | DM1                              |
| TCGA-EM-A3FQ-01A | DM1              | DM1                              |
| TCGA-EL-A3GO-01A | DM2              | DM2                              |
| TCGA-DJ-A3UZ-01A | DM1              | DM1                              |
| TCGA-E8-A432-01A | DM1              | DM1                              |
| TCGA-DJ-A4UT-01A | DM1              | DM1                              |
| TCGA-ET-A40R-01A | DM1              | DM1                              |
| TCGA-EM-A3FR-01A | DM2              | DM2                              |
| TCGA-DJ-A4UQ-01A | DM1              | DM1                              |
| TCGA-ET-A4KQ-01A | DM2              | DM2                              |
| TCGA-CE-A481-01A | DM1              | DM1                              |
| TCGA-EM-A2CP-01A | DM1              | DM1                              |
| TCGA-EL-A3H3-01A | DM1              | DM1                              |
| TCGA-EL-A3CS-01A | DM1              | DM1                              |
| TCGA-FK-A3SD-01A | DM2              | DM2                              |
| TCGA-DJ-A13W-01A | DM1              | DM1                              |
| TCGA-DJ-A13R-01A | DM1              | DM1                              |
| TCGA-EL-A3ZK-01A | DM1              | DM1                              |
| TCGA-ET-A3DS-01A | DM2              | DM2                              |
| TCGA-EL-A3D4-01A | DM1              | DM1                              |
| TCGA-EM-A3AN-01A | DM1              | DM1                              |
| TCGA-DE-A0Y2-01A | DM1              | DM1                              |
| TCGA-CE-A27D-01A | DM1              | DM1                              |
| TCGA-BJ-A0Z0-01A | DM2              | DM2                              |
| TCGA-FK-A4UB-01A | DM1              | DM1                              |
| TCGA-H2-A421-01A | DM1              | DM1                              |
| TCGA-CE-A483-01A | DM1              | DM1                              |
| TCGA-FY-A76V-01A | DM1              | DM1                              |
| TCGA-DO-A2HM-01B | DM1              | DM1                              |
| TCGA-EL-A3T0-01A | DM1              | DM1                              |
| TCGA-BJ-A0ZJ-01A | DM1              | DM1                              |
| TCGA-EL-A3T9-01A | DM1              | DM1                              |
| TCGA-EM-A1YD-01A | DM2              | DM2                              |
| TCGA-EL-A3H1-01A | DM2              | DM2                              |
| TCGA-DJ-A3VK-01A | DM1              | DM1                              |
| TCGA-EM-A4FR-01A | DM1              | DM1                              |
| TCGA-EL-A3ZR-01A | DM2              | DM2                              |
| TCGA-EM-A2P1-06A | DM1              | DM1                              |
| TCGA-ET-A39R-01A | DM1              | DM1                              |
| TCGA-ET-A3DR-01A | DM1              | DM1                              |
| TCGA-ET-A2N3-01B | DM2              | DM2                              |
| TCGA-CE-A485-01A | DM1              | DM1                              |
| TCGA-ET-A3DT-01A | DM1              | DM1                              |
| TCGA-CE-A13K-01A | DM1              | DM1                              |
| TCGA-DJ-A2Q9-01A | DM1              | DM1                              |
| TCGA-BJ-A190-01A | DM2              | DM2                              |
| TCGA-DE-A3KN-01A | DM1              | DM1                              |
| TCGA-EM-A3AL-01A | DM2              | DM2                              |
| TCGA-EM-A1YC-01A | DM2              | DM2                              |
| TCGA-BJ-A0Z5-01A | DM1              | DM1                              |
| TCGA-DO-A1JZ-01A | DM1              | DM1                              |
| TCGA-EM-A3AO-01A | DM1              | DM1                              |
| TCGA-EM-A1CW-01A | DM2              | DM2                              |
| TCGA-DJ-A3VM-01A | DM2              | DM2                              |
| TCGA-EM-A2CL-01A | DM2              | DM2                              |
| TCGA-ET-A40T-01A | DM1              | DM1                              |
| TCGA-BJ-A0ZF-01A | DM2              | DM2                              |
| TCGA-ET-A25L-01A | DM1              | DM1                              |
| TCGA-EM-A3FN-01A | DM2              | DM2                              |
| TCGA-EM-A1YA-01A | DM2              | DM2                              |
| TCGA-E3-A3E0-01A | DM1              | DM1                              |
| TCGA-BJ-A0ZC-01A | DM2              | DM2                              |
| TCGA-EM-A2CU-01A | DM1              | DM1                              |
| TCGA-DJ-A3V3-01A | DM1              | DM1                              |
| TCGA-BJ-A45K-01A | DM2              | DM2                              |
| TCGA-FY-A3R6-01A | DM1              | DM1                              |
| TCGA-IM-A41Z-01A | DM2              | DM2                              |
| TCGA-EL-A3CX-01A | DM1              | DM1                              |
| TCGA-EM-A1YB-01A | DM2              | DM2                              |
| TCGA-ET-A3DV-01A | DM2              | DM2                              |
| TCGA-BJ-A28V-01A | DM2              | DM2                              |
| TCGA-EM-A3FP-01A | DM2              | DM2                              |
| TCGA-EM-A22N-01A | DM2              | DM2                              |
| TCGA-DJ-A3VL-01A | DM1              | DM1                              |
| TCGA-EM-A3FL-01A | DM2              | DM2                              |
| TCGA-EM-A2OY-01A | DM2              | DM2                              |
| TCGA-EM-A3AJ-01A | DM1              | DM1                              |
| TCGA-FE-A238-01A | DM1              | DM1                              |
| TCGA-DJ-A2Q1-01A | DM1              | DM1                              |
| TCGA-EM-A1YE-01A | DM2              | DM2                              |
| TCGA-ET-A25N-01A | DM1              | DM1                              |
| TCGA-E8-A438-01A | DM1              | DM1                              |
| TCGA-ET-A3DQ-01A | DM1              | DM1                              |
| TCGA-J8-A3O1-01A | DM1              | DM1                              |
| TCGA-FE-A239-01A | DM2              | DM2                              |
| TCGA-DE-A4M8-01A | DM1              | DM1                              |
| TCGA-BJ-A28W-01A | DM1              | DM1                              |
| TCGA-BJ-A28T-01A | DM1              | DM1                              |
| TCGA-EM-A3FQ-06A | DM1              | DM1                              |
| TCGA-DJ-A4V0-01A | DM1              | DM1                              |
| TCGA-EL-A3H2-01A | DM2              | DM2                              |
| TCGA-EM-A2CO-01A | DM2              | DM2                              |
| TCGA-EL-A3ZM-01A | DM1              | DM1                              |
| TCGA-FE-A3PA-01A | DM1              | DM1                              |
| TCGA-ET-A25P-01A | DM1              | DM1                              |
| TCGA-BJ-A191-01A | DM2              | DM2                              |
| TCGA-EM-A3O9-01A | DM2              | DM2                              |
| TCGA-DJ-A3UV-01A | DM1              | DM1                              |
| TCGA-FE-A3PD-01A | DM1              | DM1                              |
| TCGA-ET-A3BN-01A | DM1              | DM1                              |
| TCGA-BJ-A45E-01A | DM1              | DM1                              |
| TCGA-EL-A3D5-01A | DM2              | DM2                              |
| TCGA-EL-A3CW-01A | DM1              | DM1                              |
| TCGA-KS-A41J-01A | DM1              | DM1                              |
| TCGA-BJ-A0ZA-01A | DM2              | DM2                              |
| TCGA-KS-A41L-01A | DM2              | DM2                              |
| TCGA-BJ-A28S-01A | DM2              | DM2                              |
| TCGA-BJ-A28Z-01A | DM1              | DM1                              |
| TCGA-DJ-A3US-01A | DM1              | DM1                              |
| TCGA-DE-A2OL-01A | DM2              | DM2                              |
| TCGA-FK-A3S3-01A | DM1              | DM1                              |
| TCGA-ET-A39T-01A | DM1              | DM1                              |
| TCGA-DE-A69K-01A | DM1              | DM1                              |
| TCGA-ET-A25M-01A | DM1              | DM1                              |
| TCGA-DJ-A4V5-01A | DM1              | DM1                              |
| TCGA-DE-A4M9-01A | DM2              | DM2                              |
| TCGA-ET-A39N-01A | DM2              | DM2                              |
| TCGA-E8-A416-01A | DM2              | DM2                              |
| TCGA-KS-A4ID-01A | DM1              | DM1                              |
| TCGA-DJ-A13U-01A | DM1              | DM1                              |
| TCGA-EL-A3T2-01A | DM2              | DM2                              |
| TCGA-DE-A4MA-01A | DM1              | DM1                              |
| TCGA-DJ-A2Q0-01A | DM2              | DM2                              |
| TCGA-EL-A3T6-01A | DM1              | DM1                              |
| TCGA-EM-A3O8-01A | DM2              | DM2                              |
| TCGA-DJ-A2Q2-01A | DM2              | DM2                              |
| TCGA-BJ-A0ZG-01A | DM2              | DM2                              |
| TCGA-EM-A2OW-01A | DM2              | DM2                              |
| TCGA-BJ-A3F0-01A | DM2              | DM2                              |
| TCGA-EL-A3CY-01A | DM1              | DM1                              |
| TCGA-EM-A3AQ-01A | DM1              | DM1                              |
| TCGA-J8-A4HW-06A | DM1              | DM1                              |
| TCGA-J8-A4HW-01A | DM1              | DM1                              |
| TCGA-BJ-A291-01A | DM2              | DM2                              |
| TCGA-CE-A3MD-01A | DM1              | DM1                              |
| TCGA-EL-A3ZN-01A | DM1              | DM1                              |
| TCGA-BJ-A2N9-01A | DM2              | DM2                              |
| TCGA-DJ-A4UP-01A | DM1              | DM1                              |
| TCGA-E8-A44M-01A | DM1              | DM1                              |
| TCGA-EM-A3O6-01A | DM2              | DM2                              |
| TCGA-EM-A4FH-01A | DM2              | DM2                              |
| TCGA-FY-A40N-01A | DM2              | DM2                              |
| TCGA-FY-A3WA-01A | DM2              | DM2                              |
| TCGA-ET-A2MX-01A | DM1              | DM1                              |
| TCGA-FY-A3I5-01B | DM2              | DM2                              |
| TCGA-EM-A2OV-01A | DM2              | DM2                              |
| TCGA-EM-A3AI-01A | DM2              | DM2                              |
| TCGA-FY-A4B0-01A | DM2              | DM2                              |
| TCGA-FY-A40M-01A | DM1              | DM1                              |
| TCGA-ET-A4KN-01A | DM1              | DM1                              |
| TCGA-ET-A40P-01A | DM1              | DM1                              |
| TCGA-FE-A22Z-01A | DM1              | DM1                              |
| TCGA-FK-A3SG-01A | DM1              | DM1                              |
| TCGA-EM-A22Q-01A | DM1              | DM1                              |
| TCGA-DE-A4MC-01A | DM1              | DM1                              |
| TCGA-DJ-A13S-01A | DM2              | DM2                              |
| TCGA-EM-A3SY-01A | DM2              | DM2                              |
| TCGA-EL-A3TB-01A | DM1              | DM1                              |
| TCGA-BJ-A45C-01A | DM2              | DM2                              |
| TCGA-DJ-A2PX-01A | DM1              | DM1                              |
| TCGA-EL-A3ZS-01A | DM1              | DM1                              |
### 10.3 Expression matrix 일관성
| file                                                                 | md5                                                                           |
|:---------------------------------------------------------------------|:------------------------------------------------------------------------------|
| /data/thca/data_processed/bulk_rnaseq_v3/TCGA_THCA_tumor_log2tpm.tsv | MISSING: /data/thca/data_processed/bulk_rnaseq_v3/TCGA_THCA_tumor_log2tpm.tsv |
| /data/thca/data_processed/bulk_rnaseq_v3/GSE213647_v3_log2.tsv       | f14994d312fef08645f0208c57a70820                                              |
| script                                                                                               | reference                                                      |
|:-----------------------------------------------------------------------------------------------------|:---------------------------------------------------------------|
| /home/seungho/personal/THCA_data_analysis/project/notebooks_or_scripts/v17p3_F1_external_recovery.py | /data/thca/data_processed/bulk_rnaseq_v3/GSE213647_v3_log2.tsv |
| /home/seungho/personal/THCA_data_analysis/project/notebooks_or_scripts/v17p3_common.py               | dark_matter_cohort.tsv                                         |

## SEC 11. Failure / Skip / MISSING list

### 11.1 실패한 fetch
| source                     | status       | reason                                              |
|:---------------------------|:-------------|:----------------------------------------------------|
| F2 proper comparison merge | error_in_log | ValueError: The column label pathway is not unique  |
| A1 first run               | error_in_log | columns overlap but no suffix specified: patient_id |
| A3 first run               | error_in_log | KeyError: entrezGeneId                              |
| A5 first run               | error_in_log | KeyError: Column(s) [gene] do not exist             |
### 11.2 Skip된 분석
| skip_item                                                                 |
|:--------------------------------------------------------------------------|
| F1 identifiability_no_combat / identifiability_combat columns not emitted |
| F3 multivariable Cox standalone output not emitted                        |
| A3 tumor_predicted_response content not emitted                           |
| A5 CNV burden / aneuploidy local source not present                       |
| A5 Thorsson immune subtype local table not present                        |
| A6 ROC threshold table not emitted                                        |
| Phase 2 robustness input sample_id standalone file not found              |

### 11.3 MISSING 컬럼/메트릭
| file                                   | missing_columns                                                                                                        |
|:---------------------------------------|:-----------------------------------------------------------------------------------------------------------------------|
| F1_external_5cohort_recovery.tsv       | n_dark_candidates | n_gene_overlap | DIA-AUC | identifiability_no_combat | identifiability_combat                      |
| F3_alternative_endpoints.tsv           | n_total | n_events_or_positive | test_method | statistic | p_value | effect_size | effect_size_CI                      |
| F3_logistic_regression_cluster.tsv     | coefficient | SE | p | AIC | BIC | pseudo-R2                                                                           |
| A1_per_patient_dominance.tsv           | dm1_dominant_cells | dm2_dominant_cells | dominant_phenotype                                                           |
| A2_dm_score_full_cohort.tsv            | dm1_score | dm2_score | predicted_dm                                                                                   |
| A3_top_drugs_dm1_selective.tsv         | drug_name | MOA | target | AUC_dm1_lines_mean | AUC_dm2_lines_mean | log2FC | p_value | FDR                            |
| A4_pancancer_dm_signature_transfer.tsv | n_samples | dm_signature_AUC_self_supervised | conserved_genes_n | NES_top_pathway | FDR                               |
| A6_rai_uptake_score.tsv                | predicted_dm | NIS_expr | TPO_expr | TG_expr | TSHR_expr | PAX8_expr | NKX2-1_expr | FOXE1_expr | rai_score_normalized |
## SEC 12. 모든 figure list (catalog)
| file                                                                                                             |   size_bytes | mtime                         |
|:-----------------------------------------------------------------------------------------------------------------|-------------:|:------------------------------|
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/figs/A1_dominance_clinical_correlation.html      |         8455 | 2026-04-26 17:00:14.278794765 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/figs/A1_dominance_per_patient_bar.html           |         8351 | 2026-04-26 17:00:14.241406202 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/figs/A1_scrna_celltype_dm_overlap.html           |         8560 | 2026-04-26 17:00:14.316019058 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/figs/A1_single_cell_trajectory.html              |       756362 | 2026-04-26 17:00:14.395924568 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/figs/A2_4group_outcomes.html                     |         9487 | 2026-04-26 16:59:20.658708096 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/figs/A2_braf_internal_umap.html                  |        18342 | 2026-04-26 16:59:20.618307352 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/figs/A2_dm_score_by_driver.html                  |        15259 | 2026-04-26 16:59:20.578813076 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/figs/A4_conservation_radar.html                  |         8149 | 2026-04-26 16:59:14.200968742 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/figs/A4_pancancer_heatmap.html                   |         8751 | 2026-04-26 16:59:14.145699024 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/figs/A5_genomic_instability_panel.html           |        11696 | 2026-04-26 17:00:13.697471380 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/figs/A5_immune_evasion_heatmap.html              |        25852 | 2026-04-26 17:00:13.767874956 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/figs/A5_immune_landscape.html                    |        41848 | 2026-04-26 17:00:13.736979723 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/figs/A6_pdtc_rai_validation.html                 |         9324 | 2026-04-26 16:59:16.874364376 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/figs/A6_rai_score_distribution.html              |        11680 | 2026-04-26 16:59:16.837241173 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/figs/A6_tf_circuit_diagram.html                  |         8709 | 2026-04-26 16:59:16.722294331 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/figs/F1_5cohort_dia_auc_heatmap.html             |         8581 | 2026-04-26 16:54:04.610234976 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/figs/F1_5cohort_forest.html                      |         8810 | 2026-04-26 16:54:04.568074703 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/figs/F1_recovered_gse213647_dm_distribution.html |        13507 | 2026-04-26 16:54:04.652387381 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/figs/F2_gsea_top_pathways_dual_dir.html          |         9700 | 2026-04-26 16:58:51.582866669 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/figs/F2_hallmark_heatmap_per_sample.html         |        49344 | 2026-04-26 16:58:51.622041464 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/figs/F2_pathway_volcano.html                     |        10560 | 2026-04-26 16:58:51.661794662 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/figs/F3_cluster_independence.html                |         8400 | 2026-04-26 16:52:32.462523222 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/figs/F3_decision_tree.html                       |         9196 | 2026-04-26 16:52:32.515458107 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/figs/F3_endpoint_forest.html                     |         8902 | 2026-04-26 16:52:32.426957607 |
## SEC 13. 모든 table list (catalog)
| file                                                                                                          | rows   | cols   |   size_bytes |
|:--------------------------------------------------------------------------------------------------------------|:-------|:-------|-------------:|
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/A1_cell_type_per_cluster.tsv           | 6      | 3      |          130 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/A1_dominance_vs_patient_features.tsv   | 7      | 7      |          633 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/A1_per_patient_dominance.tsv           | 7      | 5      |          571 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/A2_4group_clinical.tsv                 | 4      | 5      |          330 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/A2_braf_internal_stratification.tsv    | 4      | 6      |          329 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/A2_dm_score_full_cohort.tsv            | 513    | 31     |       192288 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/A3_celline_dm_scores.tsv               | 13     | 6      |          942 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/A3_pathway_drug_coherence.tsv          | 0      | 6      |           66 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/A3_top_drugs_dm1_selective.tsv         | 0      | 6      |           66 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/A3_top_drugs_dm2_selective.tsv         | 0      | 6      |           66 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/A3_tumor_predicted_response.tsv        | EMPTY  | EMPTY  |            1 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/A4_conserved_genes_across_cancers.tsv  | 20     | 2      |          172 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/A4_pancancer_dm_signature_transfer.tsv | 5      | 6      |          261 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/A5_cnv_aneuploidy_per_cluster.tsv      | 2      | 2      |           78 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/A5_immune_evasion_genes.tsv            | 178    | 11     |        33237 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/A5_immune_subtype_distribution.tsv     | 1      | 2      |           63 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/A5_tmb_msi_per_cluster.tsv             | 178    | 43     |        95557 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/A6_pdtc_rai_validation.tsv             | 37     | 29     |        11840 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/A6_rai_uptake_score.tsv                | 178    | 3      |         6951 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/A6_tf_activity_per_sample.tsv          | 712    | 4      |        31993 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/A6_tf_target_overlap.tsv               | 4      | 3      |          200 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/F1_external_5cohort_recovery.tsv       | 12     | 6      |          752 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/F1_external_predictions.tsv            | 1888   | 5      |        94980 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/F1_gene_id_diagnosis.tsv               | 1      | 6      |          172 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/F1_gene_recovery_mapping.tsv           | 40053  | 2      |       956137 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/F2_gsea_hallmark_proper.tsv            | 50     | 12     |        21183 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/F2_gsea_kegg_proper.tsv                | 319    | 12     |        93855 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/F2_gsea_reactome_proper.tsv            | 1698   | 12     |       386599 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/F2_proxy_vs_proper_comparison.tsv      | 0      | 7      |           83 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/F2_ssgsea_per_sample.tsv               | 178    | 1781   |      5875652 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/F3_alternative_endpoints.tsv           | 5      | 4      |          312 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/F3_effect_size_summary.tsv             | 5      | 5      |          405 |
| /home/seungho/personal/THCA_data_analysis/project/results/v17p3/tables/F3_logistic_regression_cluster.tsv     | 9      | 2      |          373 |
## SEC 14. 모든 script + line count
| file                                                                                                   |   line_count |
|:-------------------------------------------------------------------------------------------------------|-------------:|
| /home/seungho/personal/THCA_data_analysis/project/notebooks_or_scripts/v17p3_A1_scrna_heterogeneity.py |           71 |
| /home/seungho/personal/THCA_data_analysis/project/notebooks_or_scripts/v17p3_A2_dm_in_braf_ras.py      |           59 |
| /home/seungho/personal/THCA_data_analysis/project/notebooks_or_scripts/v17p3_A3_drug_response.py       |           95 |
| /home/seungho/personal/THCA_data_analysis/project/notebooks_or_scripts/v17p3_A4_pancancer.py           |           58 |
| /home/seungho/personal/THCA_data_analysis/project/notebooks_or_scripts/v17p3_A5_genomic_immune.py      |           52 |
| /home/seungho/personal/THCA_data_analysis/project/notebooks_or_scripts/v17p3_A6_tf_circuit.py          |           54 |
| /home/seungho/personal/THCA_data_analysis/project/notebooks_or_scripts/v17p3_F1_external_recovery.py   |          173 |
| /home/seungho/personal/THCA_data_analysis/project/notebooks_or_scripts/v17p3_F2_msigdb_proper.py       |          103 |
| /home/seungho/personal/THCA_data_analysis/project/notebooks_or_scripts/v17p3_F3_clinical_reframe.py    |           86 |
| /home/seungho/personal/THCA_data_analysis/project/notebooks_or_scripts/v17p3_build_report.py           |           53 |
| /home/seungho/personal/THCA_data_analysis/project/notebooks_or_scripts/v17p3_common.py                 |          147 |
| /home/seungho/personal/THCA_data_analysis/project/notebooks_or_scripts/v17p3_make_docs.py              |          120 |
| /home/seungho/personal/THCA_data_analysis/project/notebooks_or_scripts/v17p3_orchestrator.py           |           54 |
## SEC 15. Console output 핵심 메트릭 (자동 출력)
| metric                                         |      value |
|:-----------------------------------------------|-----------:|
| 1. F1 GSE213647 gene overlap recovered count   | 55         |
| 2. F1 5-cohort 중 DIA-AUC > 0.85 cohort 수       |  2         |
| 3. F2 Hallmark FDR < 0.01 pathway 수            | 15         |
| 4. F3 alternative endpoint 중 p<0.05 endpoint 수 |  4         |
| 5. A2 BRAF/DM2-like 환자 수                       |  2         |
| 6. A3 DM1-selective drug FDR<0.1 수             |  0         |
| 7. A4 DM signature applicable cancer 수         |  5         |
| 8. A6 RAI prediction AUC in PDTC               |  0.0117647 |