# THCA Project — Consolidated Status Report

_Generated: 2026-04-24 16:25:39_
- 32 HTML pages live on public URL, 100+ interactive Plotly figures, 35MB downloadable zip.
- v2 baseline TCGA internal AUC=1.00 on TierA67 / variance_top50 (LogReg+GB); external GSE27155 AUC 0.82–0.98 but bACC_youden collapses to 0.5 on linear models at default threshold.
- v3 honesty audit complete: strict-CV AUC 0.992 at k=0, drops to 0.895 at k=30; MAPK ablation nearly zero effect (0.923 → 0.924); permutation p=0.001.
- Dataset-identifiability OvR AUC = 1.00 across 6 cohorts → strong batch leakage signal; external AUC=1.00 likely reflects platform not biology.
- Top 3 gaps: (1) PRJEB11591 mutation-anchored external not obtained (supplementary URLs 403/404), (2) TCGA fusion callset empty, (3) BRS71 original vs proxy overlap 1.4% → proxy must be replaced.
## TL;DR
<!-- placeholder --> 


## 1. DATA INVENTORY

**Total samples:** 1509

**Samples per dataset** (n = 1509)

| Value | n |
|---|---:|
| GSE213647 | 632 |
| TCGA-THCA | 572 |
| GSE97466 | 141 |
| GSE27155 | 99 |
| GSE76039 | 37 |
| GSE126698 | 28 |

**Samples per platform** (n = 1509)

| Value | n |
|---|---:|
| GPL18573 | 632 |
| Illumina HiSeq / GDC STAR Counts | 572 |
| GPL13534 | 141 |
| GPL96 | 99 |
| GPL570 | 37 |
| GPL15456 | 28 |

**Samples per molecular_subtype** (n = 1509)

| Value | n |
|---|---:|
| BRAF_like | 1075 |
| unknown | 193 |
| RAS_like | 159 |
| dedifferentiated | 82 |

**Samples per histology_subtype** (n = 1509)

| Value | n |
|---|---:|
| cPTC | 1066 |
| unknown | 162 |
| FVPTC | 109 |
| normal | 58 |
| ATC | 57 |
| FTC | 30 |
| PDTC | 25 |
| MTC | 2 |

**Samples per label_confidence** (n = 1509)

| Value | n |
|---|---:|
| high | 1336 |
| low | 98 |
| medium | 64 |
| medium→excluded_relabeling | 11 |


### Dataset master

| dataset_name    | accession       | source   | modality                           | platform                       |   sample_count_reported |   sample_count_downloaded |   normal_count |   tumor_count | histology_labels_available   | molecular_labels_available   | clinical_info_available   | raw_files_available   | processed_matrix_available   | download_status   | notes                                                                                                                                                                       | download_url                                                                              |
|:----------------|:----------------|:---------|:-----------------------------------|:-------------------------------|------------------------:|--------------------------:|---------------:|--------------:|:-----------------------------|:-----------------------------|:--------------------------|:----------------------|:-----------------------------|:------------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:------------------------------------------------------------------------------------------|
| GSE213647       | GSE213647       | GEO      | bulk RNA-seq                       | GPL18573                       |                     632 |                       632 |            262 |           370 | yes                          | partial                      | partial                   | yes                   | yes                          | ok                | {"family_soft": "exists", "series_matrix": "exists", "GSE213647_RAW.tar": "exists", "filelist.txt": "exists"}                                                               | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213647/soft/GSE213647_family.soft.gz |
| GSE126698       | GSE126698       | GEO      | bulk RNA-seq                       | GPL15456                       |                      28 |                        28 |              6 |            22 | yes                          | partial                      | partial                   | no                    | yes                          | ok                | {"family_soft": "exists", "series_matrix": "exists", "GSE126698_DE_Thyroid_totalRNA_all.csv.gz": "exists", "GSE126698_DE_Thyroid_totalRNA_ATC_others_all.csv.gz": "exists"} | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126698/soft/GSE126698_family.soft.gz |
| GSE27155        | GSE27155        | GEO      | microarray                         | GPL96                          |                      99 |                        99 |              4 |            66 | yes                          | partial                      | partial                   | yes                   | yes                          | ok                | {"family_soft": "exists", "series_matrix": "exists", "GSE27155_Gior_logs.xls.gz": "exists", "GSE27155_RAW.tar": "exists", "filelist.txt": "exists"}                         | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE27nnn/GSE27155/soft/GSE27155_family.soft.gz    |
| GSE76039        | GSE76039        | GEO      | microarray                         | GPL570                         |                      37 |                        37 |              0 |            37 | yes                          | partial                      | partial                   | yes                   | yes                          | ok                | {"family_soft": "exists", "series_matrix": "exists", "GSE76039_RAW.tar": "exists", "filelist.txt": "exists"}                                                                | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE76nnn/GSE76039/soft/GSE76039_family.soft.gz    |
| GSE97466        | GSE97466        | GEO      | methylation                        | GPL13534                       |                     141 |                       141 |             67 |            74 | yes                          | no                           | partial                   | yes                   | yes                          | ok                | {"family_soft": "exists", "series_matrix": "exists", "GSE97466_RAW.tar": "exists", "filelist.txt": "exists"}                                                                | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE97nnn/GSE97466/soft/GSE97466_family.soft.gz    |
| TCGA-THCA       | TCGA-THCA       | GDC      | bulk RNA-seq                       | GDC STAR Counts                |                     572 |                       572 |             59 |           513 | yes                          | partial                      | yes                       | yes                   | yes                          | ok                | Counts via GDC API; clinical via GDC cases API; driver anchors via open masked somatic mutation MAFs.                                                                       | https://api.gdc.cancer.gov/files                                                          |
| Yoo_SK_2019_EGA | EGAD00001004845 | EGA      | controlled access / reference only | WGS/WES/RNA-seq/targeted mixed |                         |                         0 |                |               | yes                          | yes                          | partial                   | controlled            | unknown                      | held              | EGA controlled access; DAC approval and account required. Public reference only in this run.                                                                                | https://ega-archive.org/datasets/EGAD00001004845                                          |

### Gene panel coverage per dataset

| dataset   |   total_gene_count |   TDS16_coverage_n |   TDS16_coverage_frac |   BRS71_proxy_coverage_n |   BRS71_proxy_coverage_frac |   TierA67_entries_coverage_n |   TierA67_entries_coverage_frac |
|:----------|-------------------:|-------------------:|----------------------:|-------------------------:|----------------------------:|-----------------------------:|--------------------------------:|
| TCGA-THCA |              51711 |                 16 |                 1     |                       71 |                       1     |                           67 |                           1     |
| GSE213647 |              51389 |                  0 |                 0     |                        0 |                       0     |                            0 |                           0     |
| GSE126698 |              57769 |                 14 |                 0.875 |                       55 |                       0.775 |                           64 |                           0.955 |
| GSE27155  |              13236 |                 14 |                 0.875 |                       40 |                       0.563 |                           62 |                           0.925 |
| GSE76039  |              22878 |                 16 |                 1     |                       56 |                       0.789 |                           66 |                           0.985 |

### External raw-data availability

- ❌ **PRJEB11591 (ENA Korean SNU cohort)** — NOT FOUND
- ❌ **TCGA fusion callset (cBioPortal)** — NOT FOUND
- ✅ **TCGA fusion callset (alt path)** — `/opt/thyroid-dash/project/metadata/v3_fusion_anchor_tcga.tsv`
- ❌ **TCGA SCNA / GISTIC2** — NOT FOUND
- ✅ **GSE33630 SOFT family** — `/opt/thyroid-dash/project/data_raw/v3_ext/GSE33630/GSE33630_family.soft.gz`
- ✅ **GSE29265 SOFT family** — `/opt/thyroid-dash/project/data_raw/v3_ext/GSE29265/GSE29265_family.soft.gz`

## 2. ML PERFORMANCE (v2 baseline)


### ml_report_v2.md (verbatim)

```markdown
# ml_report_v2

## Setup

- Training cohort: TCGA-THCA
- Task: BRAF_like vs RAS_like
- Label source: open GDC masked somatic MAF anchor
- BRAF_like definition: BRAF V600E only
- RAS_like definition: KRAS/HRAS/NRAS mutant
- Feature sets: TDS16, TierA67, variance_top50
- BRS71 handling: proxy list derived internally from TCGA-THCA (71 genes); not used as a direct model feature in v2 because the required feature sets were fixed by task
- GSE213647 relabel status: excluded

## Internal CV

```tsv
task	dataset	feature_set	model	n_samples	cv_auc	cv_pr_auc	cv_auc_ci_lo	cv_auc_ci_hi	cv_brier	cv_balanced_accuracy	cv_balanced_accuracy_youden	cv_f1	cv_mcc	cv_youden_threshold	median_feature_count	status
BRAF_like_vs_RAS_like	TCGA-THCA	TierA67_clean	LogReg_elasticnet	333	1.0	0.9999999999999999	0.9999999999999999	1.0	0.0022207079372489847	0.9982078853046594	1.0	0.9982046678635548	0.989090046575439	0.495244215394235	55	ok
BRAF_like_vs_RAS_like	TCGA-THCA	TierA67_clean	LogReg_l2	333	1.0	0.9999999999999999	0.9999999999999999	1.0	0.001707639644917168	1.0	1.0	1.0	1.0	0.6315386087705629	55	ok
BRAF_like_vs_RAS_like	TCGA-THCA	TierA67	LogReg_elasticnet	333	1.0	0.9999999999999999	0.9999999999999999	1.0	0.00250750152371717	1.0	1.0	1.0	1.0	0.5300059365748957	66	ok
BRAF_like_vs_RAS_like	TCGA-THCA	TierA67	LogReg_l2	333	1.0	0.9999999999999999	0.9999999999999999	1.0	0.0015468648981861904	1.0	1.0	1.0	1.0	0.6360760907926487	66	ok
BRAF_like_vs_RAS_like	TCGA-THCA	variance_top50	LogReg_elasticnet	333	1.0	0.9999999999999999	0.9999999999999999	1.0	0.001863255492212747	0.9982078853046594	1.0	0.9982046678635548	0.989090046575439	0.4249419278760764	50	ok
BRAF_like_vs_RAS_like	TCGA-THCA	variance_top50	RandomForest	333	1.0	1.0	0.9999999999999999	1.0	0.003951567567567569	0.9814814814814814	1.0	0.9964285714285714	0.9778083319572793	0.722	50	ok
BRAF_like_vs_RAS_like	TCGA-THCA	variance_top50	GradientBoosting	333	1.0	1.0	0.9999999999999999	1.0	2.366541183619004e-10	1.0	1.0	1.0	1.0	0.9999925867600599	50	ok
BRAF_like_vs_RAS_like	TCGA-THCA	variance_top50	LogReg_l2	333	1.0	0.9999999999999999	0.9999999999999999	1.0	0.001218017359424287	1.0	1.0	1.0	1.0	0.5893460197847719	50	ok
BRAF_like_vs_RAS_like	TCGA-THCA	variance_top50	XGBoost	333	0.9999336253816541	0.9999871991807476	0.9996854356715948	1.0	0.002408626085184269	0.9889486260454002	0.9982078853046594	0.996415770609319	0.9778972520908005	0.8583704233169556	50	ok
BRAF_like_vs_RAS_like	TCGA-THCA	TierA67	RandomForest	333	0.9998672507633082	0.9999743524804225	0.9993553378029912	1.0	0.01204386786786787	0.9907407407407407	0.9964157706093191	0.998211091234347	0.9889267872174312	0.698	66	ok
```

## External validation

```tsv
task	dataset	feature_set	model	n_samples	auc	pr_auc	auc_ci_lo	auc_ci_hi	balanced_accuracy	balanced_accuracy_youden	f1	f1_youden	mcc	threshold_applied_youden	status
BRAF_like_vs_RAS_like	GSE126698	TierA67_clean	RandomForest	12	1.0	1.0	1.0	1.0	1.0	0.75	1.0	0.6666666666666666	1.0	0.702	ok
BRAF_like_vs_RAS_like	GSE126698	TierA67_clean	GradientBoosting	12	1.0	1.0	1.0	1.0	1.0	0.9166666666666667	1.0	0.9090909090909091	1.0	0.9984426182476762	ok
BRAF_like_vs_RAS_like	GSE126698	TierA67_clean	XGBoost	12	1.0	1.0	0.9999999999999999	1.0	0.9166666666666667	0.9166666666666667	0.9090909090909091	0.9090909090909091	0.8451542547285166	0.5997940301895142	ok
BRAF_like_vs_RAS_like	GSE126698	TierA67_clean	LogReg_elasticnet	12	1.0	1.0	0.9999999999999999	1.0	0.8333333333333333	0.8333333333333333	0.8	0.8	0.7071067811865476	0.495244215394235	ok
BRAF_like_vs_RAS_like	GSE126698	TierA67	GradientBoosting	12	1.0	1.0	1.0	1.0	1.0	0.9166666666666667	1.0	0.9090909090909091	1.0	0.9999925867600599	ok
BRAF_like_vs_RAS_like	GSE126698	TierA67	XGBoost	12	1.0	1.0	0.9999999999999999	1.0	1.0	1.0	1.0	1.0	1.0	0.47693729400634766	ok
BRAF_like_vs_RAS_like	GSE126698	TierA67	LogReg_elasticnet	12	1.0	1.0	0.9999999999999999	1.0	0.8333333333333333	0.8333333333333333	0.8	0.8	0.7071067811865476	0.5300059365748957	ok
BRAF_like_vs_RAS_like	GSE126698	TierA67	RandomForest	12	1.0	1.0	0.9999999999999999	1.0	1.0	0.75	1.0	0.6666666666666666	1.0	0.698	ok
BRAF_like_vs_RAS_like	GSE27155	variance_top50	LogReg_elasticnet	72	0.9799382716049383	0.9742637136524049	0.9388351393188854	1.0	0.5	0.5	0.0	0.0	0.0	0.4249419278760764	ok
BRAF_like_vs_RAS_like	GSE27155	TierA67_clean	LogReg_elasticnet	72	0.9783950617283951	0.9774733236461688	0.9456000982042648	0.9984459984459985	0.5	0.5	0.0	0.0	0.0	0.495244215394235	ok
```

## Interpretation

- External validation based on GSE27155 is still label-inferred, not mutation-verified.
- GSE213647 is no longer treated as a valid BRAF_like vs RAS_like external set because the supplementary sample table resolves it as PTC/PDTC/ATC/Normal without a comparable RAS-like differentiated tumor class.
- GSE126698 can be used only as a very small exploratory external set after restricting to PTC vs FTC.
```

### baseline_ml_results.tsv — internal CV

| task                  | dataset   | feature_set    | model             |   cv_auc |   cv_auc_ci_lo |   cv_auc_ci_hi |   cv_pr_auc |   cv_brier |   cv_balanced_accuracy |   cv_balanced_accuracy_youden |   cv_f1 |   cv_mcc |   cv_youden_threshold |
|:----------------------|:----------|:---------------|:------------------|---------:|---------------:|---------------:|------------:|-----------:|-----------------------:|------------------------------:|--------:|---------:|----------------------:|
| BRAF_like_vs_RAS_like | TCGA-THCA | TDS16          | LogReg_l2         |    0.996 |          0.99  |          0.999 |       0.999 |   0.0225   |                  0.971 |                         0.971 |   0.986 |    0.915 |                 0.514 |
| BRAF_like_vs_RAS_like | TCGA-THCA | TDS16          | LogReg_elasticnet |    0.996 |          0.99  |          0.999 |       0.999 |   0.0212   |                  0.971 |                         0.971 |   0.986 |    0.915 |                 0.524 |
| BRAF_like_vs_RAS_like | TCGA-THCA | TDS16          | RandomForest      |    0.988 |          0.972 |          0.997 |       0.997 |   0.0322   |                  0.893 |                         0.964 |   0.975 |    0.839 |                 0.744 |
| BRAF_like_vs_RAS_like | TCGA-THCA | TDS16          | GradientBoosting  |    0.978 |          0.939 |          0.998 |       0.99  |   0.0283   |                  0.904 |                         0.962 |   0.979 |    0.863 |                 0.995 |
| BRAF_like_vs_RAS_like | TCGA-THCA | TDS16          | XGBoost           |    0.991 |          0.981 |          0.998 |       0.998 |   0.0266   |                  0.906 |                         0.962 |   0.973 |    0.831 |                 0.947 |
| BRAF_like_vs_RAS_like | TCGA-THCA | TierA67        | LogReg_l2         |    1     |          1     |          1     |       1     |   0.00155  |                  1     |                         1     |   1     |    1     |                 0.636 |
| BRAF_like_vs_RAS_like | TCGA-THCA | TierA67        | LogReg_elasticnet |    1     |          1     |          1     |       1     |   0.00251  |                  1     |                         1     |   1     |    1     |                 0.53  |
| BRAF_like_vs_RAS_like | TCGA-THCA | TierA67        | RandomForest      |    1     |          0.999 |          1     |       1     |   0.012    |                  0.991 |                         0.996 |   0.998 |    0.989 |                 0.698 |
| BRAF_like_vs_RAS_like | TCGA-THCA | TierA67        | GradientBoosting  |    0.979 |          0.937 |          0.999 |       0.992 |   0.0153   |                  0.976 |                         0.978 |   0.991 |    0.945 |                 1     |
| BRAF_like_vs_RAS_like | TCGA-THCA | TierA67        | XGBoost           |    0.998 |          0.996 |          1     |       1     |   0.0117   |                  0.984 |                         0.985 |   0.991 |    0.946 |                 0.477 |
| BRAF_like_vs_RAS_like | TCGA-THCA | TierA67_clean  | LogReg_l2         |    1     |          1     |          1     |       1     |   0.00171  |                  1     |                         1     |   1     |    1     |                 0.632 |
| BRAF_like_vs_RAS_like | TCGA-THCA | TierA67_clean  | LogReg_elasticnet |    1     |          1     |          1     |       1     |   0.00222  |                  0.998 |                         1     |   0.998 |    0.989 |                 0.495 |
| BRAF_like_vs_RAS_like | TCGA-THCA | TierA67_clean  | RandomForest      |    1     |          0.999 |          1     |       1     |   0.0127   |                  0.991 |                         0.995 |   0.998 |    0.989 |                 0.702 |
| BRAF_like_vs_RAS_like | TCGA-THCA | TierA67_clean  | GradientBoosting  |    0.977 |          0.935 |          0.998 |       0.991 |   0.0234   |                  0.947 |                         0.98  |   0.984 |    0.9   |                 0.998 |
| BRAF_like_vs_RAS_like | TCGA-THCA | TierA67_clean  | XGBoost           |    0.998 |          0.995 |          1     |       1     |   0.0132   |                  0.965 |                         0.984 |   0.987 |    0.923 |                 0.6   |
| BRAF_like_vs_RAS_like | TCGA-THCA | variance_top50 | LogReg_l2         |    1     |          1     |          1     |       1     |   0.00122  |                  1     |                         1     |   1     |    1     |                 0.589 |
| BRAF_like_vs_RAS_like | TCGA-THCA | variance_top50 | LogReg_elasticnet |    1     |          1     |          1     |       1     |   0.00186  |                  0.998 |                         1     |   0.998 |    0.989 |                 0.425 |
| BRAF_like_vs_RAS_like | TCGA-THCA | variance_top50 | RandomForest      |    1     |          1     |          1     |       1     |   0.00395  |                  0.981 |                         1     |   0.996 |    0.978 |                 0.722 |
| BRAF_like_vs_RAS_like | TCGA-THCA | variance_top50 | GradientBoosting  |    1     |          1     |          1     |       1     |   2.37e-10 |                  1     |                         1     |   1     |    1     |                 1     |
| BRAF_like_vs_RAS_like | TCGA-THCA | variance_top50 | XGBoost           |    1     |          1     |          1     |       1     |   0.00241  |                  0.989 |                         0.998 |   0.996 |    0.978 |                 0.858 |

### baseline_ml_external.tsv — external validation

| task                  | dataset   | feature_set    | model             |   n_samples |   auc |   auc_ci_lo |   auc_ci_hi |   pr_auc |   balanced_accuracy |   balanced_accuracy_youden |    f1 |   f1_youden |   mcc | status   |
|:----------------------|:----------|:---------------|:------------------|------------:|------:|------------:|------------:|---------:|--------------------:|---------------------------:|------:|------------:|------:|:---------|
| BRAF_like_vs_RAS_like | GSE27155  | TDS16          | LogReg_l2         |          72 | 0.854 |       0.75  |       0.939 |    0.779 |               0.5   |                      0.5   | 0.667 |       0.667 | 0     | ok       |
| BRAF_like_vs_RAS_like | GSE126698 | TDS16          | LogReg_l2         |          12 | 0.75  |       0.407 |       1     |    0.85  |               0.583 |                      0.583 | 0.667 |       0.667 | 0.192 | ok       |
| BRAF_like_vs_RAS_like | GSE27155  | TDS16          | LogReg_elasticnet |          72 | 0.866 |       0.765 |       0.946 |    0.803 |               0.5   |                      0.5   | 0.667 |       0.667 | 0     | ok       |
| BRAF_like_vs_RAS_like | GSE126698 | TDS16          | LogReg_elasticnet |          12 | 0.778 |       0.429 |       1     |    0.859 |               0.583 |                      0.583 | 0.667 |       0.667 | 0.192 | ok       |
| BRAF_like_vs_RAS_like | GSE27155  | TDS16          | RandomForest      |          72 | 0.821 |       0.723 |       0.912 |    0.744 |               0.5   |                      0.5   | 0.667 |       0.667 | 0     | ok       |
| BRAF_like_vs_RAS_like | GSE126698 | TDS16          | RandomForest      |          12 | 0.778 |       0.444 |       1     |    0.859 |               0.5   |                      0.667 | 0.667 |       0.667 | 0     | ok       |
| BRAF_like_vs_RAS_like | GSE27155  | TDS16          | GradientBoosting  |          72 | 0.5   |       0.5   |       0.5   |    0.5   |               0.5   |                      0.5   | 0.667 |       0.667 | 0     | ok       |
| BRAF_like_vs_RAS_like | GSE126698 | TDS16          | GradientBoosting  |          12 | 0.75  |       0.428 |       1     |    0.708 |               0.667 |                      0.75  | 0.714 |       0.727 | 0.354 | ok       |
| BRAF_like_vs_RAS_like | GSE27155  | TDS16          | XGBoost           |          72 | 0.797 |       0.702 |       0.889 |    0.742 |               0.5   |                      0.5   | 0.667 |       0.667 | 0     | ok       |
| BRAF_like_vs_RAS_like | GSE126698 | TDS16          | XGBoost           |          12 | 0.778 |       0.469 |       1     |    0.837 |               0.667 |                      0.75  | 0.714 |       0.727 | 0.354 | ok       |
| BRAF_like_vs_RAS_like | GSE27155  | TierA67        | LogReg_l2         |          72 | 0.973 |       0.937 |       0.997 |    0.971 |               0.5   |                      0.5   | 0     |       0     | 0     | ok       |
| BRAF_like_vs_RAS_like | GSE126698 | TierA67        | LogReg_l2         |          12 | 0.861 |       0.571 |       1     |    0.897 |               0.833 |                      0.833 | 0.8   |       0.8   | 0.707 | ok       |
| BRAF_like_vs_RAS_like | GSE27155  | TierA67        | LogReg_elasticnet |          72 | 0.977 |       0.94  |       0.999 |    0.973 |               0.5   |                      0.5   | 0     |       0     | 0     | ok       |
| BRAF_like_vs_RAS_like | GSE126698 | TierA67        | LogReg_elasticnet |          12 | 1     |       1     |       1     |    1     |               0.833 |                      0.833 | 0.8   |       0.8   | 0.707 | ok       |
| BRAF_like_vs_RAS_like | GSE27155  | TierA67        | RandomForest      |          72 | 0.893 |       0.808 |       0.963 |    0.85  |               0.5   |                      0.5   | 0.667 |       0     | 0     | ok       |
| BRAF_like_vs_RAS_like | GSE126698 | TierA67        | RandomForest      |          12 | 1     |       1     |       1     |    1     |               1     |                      0.75  | 1     |       0.667 | 1     | ok       |
| BRAF_like_vs_RAS_like | GSE27155  | TierA67        | GradientBoosting  |          72 | 0.744 |       0.635 |       0.84  |    0.705 |               0.736 |                      0.5   | 0.771 |       0     | 0.496 | ok       |
| BRAF_like_vs_RAS_like | GSE126698 | TierA67        | GradientBoosting  |          12 | 1     |       1     |       1     |    1     |               1     |                      0.917 | 1     |       0.909 | 1     | ok       |
| BRAF_like_vs_RAS_like | GSE27155  | TierA67        | XGBoost           |          72 | 0.769 |       0.642 |       0.868 |    0.761 |               0.5   |                      0.5   | 0     |       0     | 0     | ok       |
| BRAF_like_vs_RAS_like | GSE126698 | TierA67        | XGBoost           |          12 | 1     |       1     |       1     |    1     |               1     |                      1     | 1     |       1     | 1     | ok       |
| BRAF_like_vs_RAS_like | GSE27155  | TierA67_clean  | LogReg_l2         |          72 | 0.968 |       0.929 |       0.995 |    0.97  |               0.5   |                      0.5   | 0     |       0     | 0     | ok       |
| BRAF_like_vs_RAS_like | GSE126698 | TierA67_clean  | LogReg_l2         |          12 | 0.889 |       0.656 |       1     |    0.911 |               0.75  |                      0.75  | 0.667 |       0.667 | 0.577 | ok       |
| BRAF_like_vs_RAS_like | GSE27155  | TierA67_clean  | LogReg_elasticnet |          72 | 0.978 |       0.946 |       0.998 |    0.977 |               0.5   |                      0.5   | 0     |       0     | 0     | ok       |
| BRAF_like_vs_RAS_like | GSE126698 | TierA67_clean  | LogReg_elasticnet |          12 | 1     |       1     |       1     |    1     |               0.833 |                      0.833 | 0.8   |       0.8   | 0.707 | ok       |
| BRAF_like_vs_RAS_like | GSE27155  | TierA67_clean  | RandomForest      |          72 | 0.912 |       0.828 |       0.973 |    0.86  |               0.5   |                      0.5   | 0.667 |       0     | 0     | ok       |
| BRAF_like_vs_RAS_like | GSE126698 | TierA67_clean  | RandomForest      |          12 | 1     |       1     |       1     |    1     |               1     |                      0.75  | 1     |       0.667 | 1     | ok       |
| BRAF_like_vs_RAS_like | GSE27155  | TierA67_clean  | GradientBoosting  |          72 | 0.744 |       0.635 |       0.84  |    0.705 |               0.736 |                      0.542 | 0.771 |       0.154 | 0.496 | ok       |
| BRAF_like_vs_RAS_like | GSE126698 | TierA67_clean  | GradientBoosting  |          12 | 1     |       1     |       1     |    1     |               1     |                      0.917 | 1     |       0.909 | 1     | ok       |
| BRAF_like_vs_RAS_like | GSE27155  | TierA67_clean  | XGBoost           |          72 | 0.851 |       0.747 |       0.931 |    0.824 |               0.5   |                      0.5   | 0     |       0     | 0     | ok       |
| BRAF_like_vs_RAS_like | GSE126698 | TierA67_clean  | XGBoost           |          12 | 1     |       1     |       1     |    1     |               0.917 |                      0.917 | 0.909 |       0.909 | 0.845 | ok       |
| BRAF_like_vs_RAS_like | GSE27155  | variance_top50 | LogReg_l2         |          72 | 0.978 |       0.939 |       1     |    0.973 |               0.5   |                      0.5   | 0     |       0     | 0     | ok       |
| BRAF_like_vs_RAS_like | GSE126698 | variance_top50 | LogReg_l2         |          12 | 0.944 |       0.75  |       1     |    0.958 |               0.75  |                      0.667 | 0.667 |       0.5   | 0.577 | ok       |
| BRAF_like_vs_RAS_like | GSE27155  | variance_top50 | LogReg_elasticnet |          72 | 0.98  |       0.939 |       1     |    0.974 |               0.5   |                      0.5   | 0     |       0     | 0     | ok       |
| BRAF_like_vs_RAS_like | GSE126698 | variance_top50 | LogReg_elasticnet |          12 | 0.944 |       0.75  |       1     |    0.958 |               0.667 |                      0.75  | 0.5   |       0.667 | 0.447 | ok       |
| BRAF_like_vs_RAS_like | GSE27155  | variance_top50 | RandomForest      |          72 | 0.975 |       0.931 |       1     |    0.965 |               0.5   |                      0.5   | 0     |       0     | 0     | ok       |
| BRAF_like_vs_RAS_like | GSE126698 | variance_top50 | RandomForest      |          12 | 0.861 |       0.562 |       1     |    0.924 |               0.75  |                      0.667 | 0.667 |       0.5   | 0.577 | ok       |
| BRAF_like_vs_RAS_like | GSE27155  | variance_top50 | GradientBoosting  |          72 | 0.5   |       0.5   |       0.5   |    0.5   |               0.5   |                      0.5   | 0     |       0     | 0     | ok       |
| BRAF_like_vs_RAS_like | GSE126698 | variance_top50 | GradientBoosting  |          12 | 0.75  |       0.5   |       1     |    0.679 |               0.75  |                      0.75  | 0.769 |       0.769 | 0.507 | ok       |
| BRAF_like_vs_RAS_like | GSE27155  | variance_top50 | XGBoost           |          72 | 0.723 |       0.635 |       0.808 |    0.711 |               0.5   |                      0.5   | 0     |       0     | 0     | ok       |
| BRAF_like_vs_RAS_like | GSE126698 | variance_top50 | XGBoost           |          12 | 0.792 |       0.5   |       1     |    0.774 |               0.75  |                      0.75  | 0.769 |       0.769 | 0.507 | ok       |

### Quantum algorithm comparison

```markdown
# Quantum ML full benchmark (publication-style)

Classical LogReg_l2 reference AUC = **1.0000**.

| Algorithm | Kind | CV AUC [95% CI] | PR-AUC | Wall (s) | Verdict vs classical |
|---|---|---|---|---|---|
| LogReg_l2 (classical baseline) | classical | 1.000 [1.000, 1.000] | 1.000 | 0.0 | TIE or WIN (CI above classical baseline) |
| MLP (classical baseline) | classical | 0.981 [0.955, 0.997] | 0.995 | 0.2 | LOSE |
| VQC (ZZ feature map) | quantum | 0.488 [0.410, 0.571] | 0.841 | 88.0 | LOSE |
| VQC (PAULI feature map) | quantum | 0.782 [0.729, 0.833] | 0.955 | 44.6 | LOSE |
| VQC (AMP feature map) | quantum | 0.844 [0.782, 0.902] | 0.964 | 34.2 | LOSE |
| QSVM (fidelity kernel) | quantum | 0.754 [0.685, 0.818] | 0.935 | 76.6 | LOSE |
| QNN (EstimatorQNN + 3 ansatz layers) | quantum | 0.594 [0.516, 0.672] | 0.881 | 25.6 | LOSE |
| Quantum Kitchen Sinks | quantum-inspired | 0.491 [0.420, 0.565] | 0.835 | 23.3 | LOSE |
| QAOA feature selection -> LogReg | quantum | 1.000 [1.000, 1.000] | 1.000 | 0.0 | TIE or WIN (CI above classical baseline) |
| dwave-neal SA feature selection -> LogReg | classical-SA | 1.000 [1.000, 1.000] | 1.000 | 0.0 | TIE or WIN (CI above classical baseline) |
| Grover-inspired AA (demo, sqrt(N) asymptotic) | quantum-demo | 1.000 [1.000, 1.000] | 1.000 | 0.0 | TIE or WIN (CI above classical baseline) |
| Quantum k-means (swap-test, k=2) | quantum-unsupervised | — | nan | 1.6 | N/A (unsupervised) |
| QBoost (16 weak stumps, QUBO w/ neal) | quantum-inspired-ensemble | 0.976 [0.962, 0.989] | 0.993 | 0.2 | LOSE |

## Honesty clause
- `dwave-neal` is a **classical** simulated annealer; it is NOT quantum hardware.
- QAOA + Aer runs on a noiseless classical simulator (`qiskit_aer.primitives.Sampler`).
- PennyLane `lightning.qubit` is also a classical statevector simulator.
- All AUC CIs are 95% percentile bootstraps with 800 resamples on held-out CV predictions.
- Grover-inspired demo has O(sqrt(N)) asymptotic speedup but at n=55 genes that is irrelevant.
- Quantum k-means is unsupervised, so AUC is not reported; ARI vs label is the metric.

## qPCA vs classical PCA (first 8 eigenvalues)

| Rank | Classical PCA | Quantum-style PCA (rho eig) |
|---|---|---|
| 1 | 13.1760 | 0.9075 |
| 2 | 7.8659 | 0.0242 |
| 3 | 5.2459 | 0.0109 |
| 4 | 4.2053 | 0.0090 |
| 5 | 3.0079 | 0.0059 |
| 6 | 1.9111 | 0.0049 |
| 7 | 1.6845 | 0.0033 |
| 8 | 1.4480 | 0.0030 |

## Grover demo
- qubits: 3, iterations: 2
- target index: 4
- target probability after AA: 0.9453124999999959
```

## 3. BIOMARKER + DRUG (v2)


### biomarker_analysis.md TL;DR

## TL;DR

Tested 51,711 genes on TCGA-THCA BRAF_like (n=392) vs RAS_like (n=111). 148 of the 221 hardcoded + panel reference genes were present; 70 replicated (same direction, external p<0.05) in ≥1 external cohort. 2773 genes passed the novel-validated filter (FDR<0.05, |d|>0.5, replication≥1, median log2>1, present in ≥2 cohorts). Top novel candidates: **TACSTD2** (d=+2.52, replicated 2/2), **TMPRSS4** (d=+2.09, replicated 2/2), **PLEKHA6** (d=+2.22, replicated 2/2). Genes with large TCGA effects that fail to replicate in BOTH externals are **not** called biomarkers.


### biomarker_validated.tsv — top 20 by score

**Total rows:** 2843

| gene       |   log2FC_tcga |   cohens_d_tcga |   pval_tcga |   braf_mean_tcga |   ras_mean_tcga |   n_braf_tcga |   n_ras_tcga |   fdr_tcga |   log2FC_27155 |
|:-----------|--------------:|----------------:|------------:|-----------------:|----------------:|--------------:|-------------:|-----------:|---------------:|
| DCSTAMP    |          5.35 |            2.24 |    1.35e-86 |             6.25 |           0.906 |           392 |          111 |   1.75e-82 |          0.915 |
| KCNN4      |          3.78 |            2.37 |    2.85e-62 |             5.43 |           1.65  |           392 |          111 |   4.47e-59 |          0.971 |
| TACSTD2    |          4.87 |            2.52 |    6.63e-55 |             8.42 |           3.56  |           392 |          111 |   6.23e-52 |          0.913 |
| TMPRSS4    |          4.13 |            2.09 |    2.01e-62 |             5.62 |           1.49  |           392 |          111 |   3.25e-59 |          0.815 |
| PLEKHA6    |          2.35 |            2.22 |    2.53e-57 |             4.11 |           1.76  |           392 |          111 |   3.04e-54 |          0.135 |
| CYP1B1     |          3.47 |            2.22 |    9.04e-56 |             6.28 |           2.81  |           392 |          111 |   8.99e-53 |          0.863 |
| BNC1       |          1.65 |            1.63 |    1.54e-69 |             2.03 |           0.381 |           392 |          111 |   4.99e-66 |          0.111 |
| LDLR       |          2.24 |            2.02 |    2.56e-56 |             4.8  |           2.56  |           392 |          111 |   2.76e-53 |          0.337 |
| GABRB2     |          3.7  |            2.24 |    5.11e-50 |             6.49 |           2.79  |           392 |          111 |   2.9e-47  |          0.139 |
| ST6GALNAC5 |          2.91 |            1.71 |    7.06e-65 |             3.55 |           0.645 |           392 |          111 |   1.49e-61 |          0.52  |
| B3GNT3     |          3.78 |            2.02 |    8.83e-55 |             5.02 |           1.24  |           392 |          111 |   8.01e-52 |          0.333 |
| PTPRE      |          2.38 |            2.26 |    8.16e-48 |             6.25 |           3.88  |           392 |          111 |   3.83e-45 |          0.533 |
| KCNQ3      |          2.48 |            1.96 |    6.13e-54 |             5    |           2.52  |           392 |          111 |   4.95e-51 |          0.119 |
| CREB5      |          2.09 |            2.03 |    4.92e-52 |             3.99 |           1.89  |           392 |          111 |   3.59e-49 |          0.104 |
| LY6E       |          2.13 |            2.16 |    1.61e-45 |             7.66 |           5.54  |           392 |          111 |   5.84e-43 |          0.527 |
| TAGLN2     |          1.22 |            2.18 |    6.8e-45  |             8.68 |           7.45  |           392 |          111 |   2.35e-42 |          0.277 |
| PDLIM4     |          3.44 |            2.25 |    1.14e-42 |             6.93 |           3.49  |           392 |          111 |   3.07e-40 |          0.997 |
| ITGA3      |          1.52 |            2.14 |    1.55e-44 |             9.82 |           8.3   |           392 |          111 |   5.01e-42 |          0.323 |
| SERPINA1   |          4.19 |            2.38 |    7.89e-40 |            10.7  |           6.5   |           392 |          111 |   1.68e-37 |          0.847 |
| BID        |          1.35 |            2.13 |    5.78e-44 |             5.36 |           4.01  |           392 |          111 |   1.74e-41 |          0.325 |

### biomarker_to_drug_report.md highlights

## Top 8 targets

| Gene | Novel | Classification | #Cmpds | Best pChEMBL | Thyroid papers | Top PMID | Top Compound |
|---|---|---|---|---|---|---|---|
| TACSTD2 | True | novel_target | 0 | — | 4 | 41413112 | — |
| TMPRSS4 | True | novel_target | 0 | — | 5 | 41656803 | — |
| PLEKHA6 | True | novel_target | 0 | — | 0 | — | — |
| CYP1B1 | True | validated_target | 10 | 8.75 | 5 | 40900788 | CHEMBL3132932 |
| LDLR | True | emerging_target | 4 | 6.22 | 5 | 41340871 | CHEMBL115992 |
| GABRB2 | True | emerging_target | 1 | 6.89 | 5 | 40900788 | CANNABIDIOL |
| B3GNT3 | True | novel_target | 0 | — | 0 | — | — |
| PTPRE | True | novel_target | 0 | — | 3 | 36654463 | — |

## Pharma-quality highlights (3 novel targets)
### TACSTD2 — novel_target
**TACSTD2** surfaces from biomarker discovery with high novelty_score and no known high-confidence chemical starting points in ChEMBL. For a pharma portfolio this is a **first-in-class opportunity** — the target hypothesis is de-risked by replicated differential expression in thyroid cancer cohorts, while the chemical matter is greenfield (fragment screening / DNA-encoded library campaigns would be the logical next step). Thyroid-specific literature count = 4.

### TMPRSS4 — novel_target
**TMPRSS4** surfaces from biomarker discovery with high novelty_score and no known high-confidence chemical starting points in ChEMBL. For a pharma portfolio this is a **first-in-class opportunity** — the target hypothesis is de-risked by replicated differential expression in thyroid cancer cohorts, while the chemical matter is greenfield (fragment screening / DNA-encoded library campaigns would be the logical next step). Thyroid-specific literature count = 5.

### PLEKHA6 — novel_target
**PLEKHA6** surfaces from biomarker discovery with high novelty_score and no known high-confidence chemical starting points in ChEMBL. For a pharma portfolio this is a **first-in-class opportunity** — the target hypothesis is de-risked by replicated differential expression in thyroid cancer cohorts, while the chemical matter is greenfield (fragment screening / DNA-encoded library campaigns would be the logical next step). Thyroid-specific literature count = 0.


## Caveat

### drug_discovery_compounds.tsv

| target   | chembl_id     | target_relation_confidence   |
|:---------|:--------------|:-----------------------------|
| CYP1B1   | CHEMBL3132932 | DIRECT                       |
| CYP1B1   | CHEMBL46909   | DIRECT                       |
| CYP1B1   | CHEMBL3132924 | DIRECT                       |
| CYP1B1   | CHEMBL243664  | DIRECT                       |
| CYP1B1   | CHEMBL379064  | DIRECT                       |
| CYP1B1   | CHEMBL226034  | DIRECT                       |
| CYP1B1   | CHEMBL214321  | DIRECT                       |
| CYP1B1   | CHEMBL117     | DIRECT                       |
| CYP1B1   | CHEMBL40919   | DIRECT                       |
| CYP1B1   | CHEMBL309490  | DIRECT                       |
| LDLR     | CHEMBL115992  | DIRECT                       |
| LDLR     | CHEMBL114890  | DIRECT                       |
| LDLR     | CHEMBL324616  | DIRECT                       |
| LDLR     | CHEMBL112449  | DIRECT                       |
| GABRB2   | CHEMBL190461  | DIRECT                       |

## 4. v3 PROGRESS

### ✅ `results/ml/v3_ext_validation_results.tsv`
_Rows: 76_

| task   | external   |   n_ext | feature_set           |   n_features | model             |   auc_pre |   auc_lo |   auc_hi |   pr_auc_pre |   brier_pre |   ece_pre |   bACC_Youden_pre |   F1_pre |   NPV_at_Se95_pre |   auc_post |   brier_post |   ece_post |   bACC_Youden_post |   n_classes |   n_samples |   auc_ovr_cv | classes                     | class_counts                                        | pool_tag   |   auc |   brier |   ece |   bACC_Youden | per_cohort_bacc                  | per_cohort_auc                   |
|:-------|:-----------|--------:|:----------------------|-------------:|:------------------|----------:|---------:|---------:|-------------:|------------:|----------:|------------------:|---------:|------------------:|-----------:|-------------:|-----------:|-------------------:|------------:|------------:|-------------:|:----------------------------|:----------------------------------------------------|:-----------|------:|--------:|------:|--------------:|:---------------------------------|:---------------------------------|
| V1     | GSE27155   |      41 | TDS16                 |           14 | LogReg_l2         |     0.997 |    0.985 |    1     |        0.999 |      0.317  |     0.317 |             0.982 |    0.982 |             0.929 |      0.5   |       0.317  |     0.317  |              0.5   |             |             |              | nan                         | nan                                                 | nan        |       |         |       |               | nan                              | nan                              |
| V1     | GSE27155   |      41 | TDS16                 |           14 | LogReg_elasticnet |     0.997 |    0.985 |    1     |        0.999 |      0.317  |     0.317 |             0.982 |    0.982 |             0.929 |      0.5   |       0.317  |     0.317  |              0.5   |             |             |              | nan                         | nan                                                 | nan        |       |         |       |               | nan                              | nan                              |
| V1     | GSE27155   |      41 | TDS16                 |           14 | RandomForest      |     0.828 |    0.621 |    0.989 |        0.873 |      0.271  |     0.243 |             0.828 |    0.915 |             0.9   |      0.5   |       0.317  |     0.317  |              0.5   |             |             |              | nan                         | nan                                                 | nan        |       |         |       |               | nan                              | nan                              |
| V1     | GSE27155   |      41 | TDS16                 |           14 | GradientBoosting  |     0.5   |    0.5   |    0.5   |        0.683 |      0.317  |     0.317 |             0.5   |    0     |             0     |      0.5   |       0.317  |     0.317  |              0.5   |             |             |              | nan                         | nan                                                 | nan        |       |         |       |               | nan                              | nan                              |
| V1     | GSE27155   |      41 | TDS16                 |           14 | XGBoost           |     0.861 |    0.701 |    1     |        0.897 |      0.316  |     0.316 |             0.885 |    0.949 |             0.909 |      0.5   |       0.317  |     0.317  |              0.5   |             |             |              | nan                         | nan                                                 | nan        |       |         |       |               | nan                              | nan                              |
| V1     | GSE27155   |      41 | TierA67_clean         |           50 | LogReg_l2         |     0.997 |    0.973 |    1     |        0.999 |      0.682  |     0.682 |             0.982 |    0.982 |             0.929 |      0.5   |       0.683  |     0      |              0.5   |             |             |              | nan                         | nan                                                 | nan        |       |         |       |               | nan                              | nan                              |
| V1     | GSE27155   |      41 | TierA67_clean         |           50 | LogReg_elasticnet |     0.964 |    0.892 |    1     |        0.988 |      0.683  |     0.683 |             0.964 |    0.963 |             0.923 |      0.5   |       0.683  |     0      |              0.5   |             |             |              | nan                         | nan                                                 | nan        |       |         |       |               | nan                              | nan                              |
| V1     | GSE27155   |      41 | TierA67_clean         |           50 | RandomForest      |     0.902 |    0.812 |    0.967 |        0.955 |      0.27   |     0.31  |             0.857 |    0.833 |             1     |      0.902 |       0.345  |     0.399  |              0.857 |             |             |              | nan                         | nan                                                 | nan        |       |         |       |               | nan                              | nan                              |
| V1     | GSE27155   |      41 | TierA67_clean         |           50 | GradientBoosting  |     0.592 |    0.424 |    0.734 |        0.8   |      0.448  |     0.491 |             0.727 |    0.723 |             0     |      0.592 |       0.448  |     0.491  |              0.727 |             |             |              | nan                         | nan                                                 | nan        |       |         |       |               | nan                              | nan                              |
| V1     | GSE27155   |      41 | TierA67_clean         |           50 | XGBoost           |     0.949 |    0.859 |    1     |        0.965 |      0.473  |     0.538 |             0.908 |    0.926 |             0.917 |      0.643 |       0.597  |     0.143  |              0.643 |             |             |              | nan                         | nan                                                 | nan        |       |         |       |               | nan                              | nan                              |
| V1     | GSE27155   |      41 | TierA67_clean_no_MAPK |           40 | LogReg_l2         |     0.995 |    0.966 |    1     |        0.997 |      0.127  |     0.266 |             0.964 |    0.963 |             1     |      0.857 |       0.2    |     0.011  |              0.857 |             |             |              | nan                         | nan                                                 | nan        |       |         |       |               | nan                              | nan                              |
| V1     | GSE27155   |      41 | TierA67_clean_no_MAPK |           40 | LogReg_elasticnet |     1     |    1     |    1     |        1     |      0.68   |     0.681 |             1     |    1     |             1     |      0.5   |       0.683  |     0      |              0.5   |             |             |              | nan                         | nan                                                 | nan        |       |         |       |               | nan                              | nan                              |
| V1     | GSE27155   |      41 | TierA67_clean_no_MAPK |           40 | RandomForest      |     0.933 |    0.803 |    0.997 |        0.959 |      0.24   |     0.201 |             0.905 |    0.947 |             0.917 |      0.933 |       0.242  |     0.292  |              0.905 |             |             |              | nan                         | nan                                                 | nan        |       |         |       |               | nan                              | nan                              |
| V1     | GSE27155   |      41 | TierA67_clean_no_MAPK |           40 | GradientBoosting  |     0.87  |    0.717 |    0.982 |        0.927 |      0.288  |     0.379 |             0.87  |    0.909 |             1     |      0.87  |       0.288  |     0.379  |              0.87  |             |             |              | nan                         | nan                                                 | nan        |       |         |       |               | nan                              | nan                              |
| V1     | GSE27155   |      41 | TierA67_clean_no_MAPK |           40 | XGBoost           |     0.995 |    0.978 |    1     |        0.995 |      0.45   |     0.548 |             0.964 |    0.963 |             1     |      0.679 |       0.557  |     0.161  |              0.679 |             |             |              | nan                         | nan                                                 | nan        |       |         |       |               | nan                              | nan                              |
| V1     | GSE27155   |      41 | BRS71_original        |           59 | LogReg_l2         |     1     |    1     |    1     |        1     |      0.213  |     0.256 |             1     |    1     |             1     |      0.654 |       0.303  |     0.309  |              0.654 |             |             |              | nan                         | nan                                                 | nan        |       |         |       |               | nan                              | nan                              |
| V1     | GSE27155   |      41 | BRS71_original        |           59 | LogReg_elasticnet |     1     |    1     |    1     |        1     |      0.0725 |     0.166 |             1     |    1     |             1     |      0.962 |       0.0279 |     0.0375 |              0.962 |             |             |              | nan                         | nan                                                 | nan        |       |         |       |               | nan                              | nan                              |
| V1     | GSE27155   |      41 | BRS71_original        |           59 | RandomForest      |     1     |    1     |    1     |        1     |      0.254  |     0.259 |             1     |    1     |             0.929 |      1     |       0.285  |     0.507  |              1     |             |             |              | nan                         | nan                                                 | nan        |       |         |       |               | nan                              | nan                              |
| V1     | GSE27155   |      41 | BRS71_original        |           59 | GradientBoosting  |     0.931 |    0.841 |    0.993 |        0.973 |      0.653  |     0.666 |             0.908 |    0.926 |             0.9   |      0.931 |       0.653  |     0.666  |              0.908 |             |             |              | nan                         | nan                                                 | nan        |       |         |       |               | nan                              | nan                              |
| V1     | GSE27155   |      41 | BRS71_original        |           59 | XGBoost           |     0.837 |    0.695 |    0.973 |        0.928 |      0.196  |     0.333 |             0.908 |    0.926 |             0     |      0.837 |       0.196  |     0.293  |              0.908 |             |             |              | nan                         | nan                                                 | nan        |       |         |       |               | nan                              | nan                              |
| V1     | GSE27155   |      41 | v3_panel_compact      |           22 | LogReg_l2         |     0.997 |    0.986 |    1     |        0.999 |      0.311  |     0.313 |             0.982 |    0.982 |             0.929 |      0.5   |       0.317  |     0.317  |              0.5   |             |             |              | nan                         | nan                                                 | nan        |       |         |       |               | nan                              | nan                              |
| V1     | GSE27155   |      41 | v3_panel_compact      |           22 | LogReg_elasticnet |     1     |    1     |    1     |        1     |      0.314  |     0.315 |             1     |    1     |             1     |      0.5   |       0.317  |     0.317  |              0.5   |             |             |              | nan                         | nan                                                 | nan        |       |         |       |               | nan                              | nan                              |
| V1     | GSE27155   |      41 | v3_panel_compact      |           22 | RandomForest      |     1     |    1     |    1     |        1     |      0.174  |     0.404 |             1     |    1     |             1     |      1     |       0.174  |     0.231  |              1     |             |             |              | nan                         | nan                                                 | nan        |       |         |       |               | nan                              | nan                              |
| V1     | GSE27155   |      41 | v3_panel_compact      |           22 | GradientBoosting  |     0.5   |    0.5   |    0.5   |        0.683 |      0.683  |     0.683 |             0.5   |    0     |             0     |      0.5   |       0.683  |     0.683  |              0.5   |             |             |              | nan                         | nan                                                 | nan        |       |         |       |               | nan                              | nan                              |
| V1     | GSE27155   |      41 | v3_panel_compact      |           22 | XGBoost           |     0.5   |    0.5   |    0.5   |        0.683 |      0.448  |     0.481 |             0.5   |    0     |             0     |      0.5   |       0.552  |     0.579  |              0.5   |             |             |              | nan                         | nan                                                 | nan        |       |         |       |               | nan                              | nan                              |
| V3     | nan        |         | nan                   |              | nan               |           |          |          |              |             |           |                   |          |                   |            |              |            |                    |           3 |         572 |        0.939 | BRAF_V600E,RAS_mutant,other | {"BRAF_V600E": 280, "other": 238, "RAS_mutant": 54} | nan        |       |         |       |               | nan                              | nan                              |
| V4     | nan        |         | TDS16                 |              | LogReg_l2         |           |    0.985 |    1     |              |             |           |                   |          |                   |            |              |            |                    |             |             |              | nan                         | nan                                                 | pre_combat | 0.997 |  0.317  | 0.317 |         0.982 | {"GSE27155": 0.9821428571428572} | {"GSE27155": 0.9972527472527473} |
| V4     | nan        |         | TDS16                 |              | LogReg_elasticnet |           |    0.985 |    1     |              |             |           |                   |          |                   |            |              |            |                    |             |             |              | nan                         | nan                                                 | pre_combat | 0.997 |  0.317  | 0.317 |         0.982 | {"GSE27155": 0.9821428571428572} | {"GSE27155": 0.9972527472527473} |
| V4     | nan        |         | TDS16                 |              | RandomForest      |           |    0.621 |    0.989 |              |             |           |                   |          |                   |            |              |            |                    |             |             |              | nan                         | nan                                                 | pre_combat | 0.828 |  0.271  | 0.243 |         0.828 | {"GSE27155": 0.8282967032967032} | {"GSE27155": 0.8282967032967034} |
| V4     | nan        |         | TDS16                 |              | GradientBoosting  |           |    0.5   |    0.5   |              |             |           |                   |          |                   |            |              |            |                    |             |             |              | nan                         | nan                                                 | pre_combat | 0.5   |  0.317  | 0.317 |         0.5   | {"GSE27155": 0.5}                | {"GSE27155": 0.5}                |
| V4     | nan        |         | TDS16                 |              | XGBoost           |           |    0.701 |    1     |              |             |           |                   |          |                   |            |              |            |                    |             |             |              | nan                         | nan                                                 | pre_combat | 0.861 |  0.316  | 0.316 |         0.885 | {"GSE27155": 0.8846153846153846} | {"GSE27155": 0.8612637362637363} |
| V4     | nan        |         | TierA67_clean         |              | LogReg_l2         |           |    0.973 |    1     |              |             |           |                   |          |                   |            |              |            |                    |             |             |              | nan                         | nan                                                 | pre_combat | 0.997 |  0.682  | 0.682 |         0.982 | {"GSE27155": 0.9821428571428572} | {"GSE27155": 0.9972527472527473} |
| V4     | nan        |         | TierA67_clean         |              | LogReg_elasticnet |           |    0.892 |    1     |              |             |           |                   |          |                   |            |              |            |                    |             |             |              | nan                         | nan                                                 | pre_combat | 0.964 |  0.683  | 0.683 |         0.964 | {"GSE27155": 0.9642857142857143} | {"GSE27155": 0.9642857142857143} |
| V4     | nan        |         | TierA67_clean         |              | RandomForest      |           |    0.812 |    0.967 |              |             |           |                   |          |                   |            |              |            |                    |             |             |              | nan                         | nan                                                 | pre_combat | 0.902 |  0.27   | 0.31  |         0.857 | {"GSE27155": 0.8571428571428572} | {"GSE27155": 0.9024725274725276} |
| V4     | nan        |         | TierA67_clean         |              | GradientBoosting  |           |    0.424 |    0.734 |              |             |           |                   |          |                   |            |              |            |                    |             |             |              | nan                         | nan                                                 | pre_combat | 0.592 |  0.448  | 0.491 |         0.727 | {"GSE27155": 0.7266483516483516} | {"GSE27155": 0.5920329670329669} |
| V4     | nan        |         | TierA67_clean         |              | XGBoost           |           |    0.859 |    1     |              |             |           |                   |          |                   |            |              |            |                    |             |             |              | nan                         | nan                                                 | pre_combat | 0.949 |  0.473  | 0.538 |         0.908 | {"GSE27155": 0.9079670329670331} | {"GSE27155": 0.9491758241758241} |
| V4     | nan        |         | TierA67_clean_no_MAPK |              | LogReg_l2         |           |    0.966 |    1     |              |             |           |                   |          |                   |            |              |            |                    |             |             |              | nan                         | nan                                                 | pre_combat | 0.995 |  0.127  | 0.266 |         0.964 | {"GSE27155": 0.9642857142857143} | {"GSE27155": 0.9945054945054945} |
| V4     | nan        |         | TierA67_clean_no_MAPK |              | LogReg_elasticnet |           |    1     |    1     |              |             |           |                   |          |                   |            |              |            |                    |             |             |              | nan                         | nan                                                 | pre_combat | 1     |  0.68   | 0.681 |         1     | {"GSE27155": 1.0}                | {"GSE27155": 1.0}                |
| V4     | nan        |         | TierA67_clean_no_MAPK |              | RandomForest      |           |    0.803 |    0.997 |              |             |           |                   |          |                   |            |              |            |                    |             |             |              | nan                         | nan                                                 | pre_combat | 0.933 |  0.24   | 0.201 |         0.905 | {"GSE27155": 0.9052197802197802} | {"GSE27155": 0.9326923076923077} |
| V4     | nan        |         | TierA67_clean_no_MAPK |              | GradientBoosting  |           |    0.717 |    0.982 |              |             |           |                   |          |                   |            |              |            |                    |             |             |              | nan                         | nan                                                 | pre_combat | 0.87  |  0.288  | 0.379 |         0.87  | {"GSE27155": 0.8695054945054945} | {"GSE27155": 0.8695054945054945} |
| V4     | nan        |         | TierA67_clean_no_MAPK |              | XGBoost           |           |    0.978 |    1     |              |             |           |                   |          |                   |            |              |            |                    |             |             |              | nan                         | nan                                                 | pre_combat | 0.995 |  0.45   | 0.548 |         0.964 | {"GSE27155": 0.9642857142857143} | {"GSE27155": 0.9945054945054945} |
| V4     | nan        |         | BRS71_original        |              | LogReg_l2         |           |    1     |    1     |              |             |           |                   |          |                   |            |              |            |                    |             |             |              | nan                         | nan                                                 | pre_combat | 1     |  0.213  | 0.256 |         1     | {"GSE27155": 1.0}                | {"GSE27155": 1.0}                |
| V4     | nan        |         | BRS71_original        |              | LogReg_elasticnet |           |    1     |    1     |              |             |           |                   |          |                   |            |              |            |                    |             |             |              | nan                         | nan                                                 | pre_combat | 1     |  0.0725 | 0.166 |         1     | {"GSE27155": 1.0}                | {"GSE27155": 1.0}                |
| V4     | nan        |         | BRS71_original        |              | RandomForest      |           |    1     |    1     |              |             |           |                   |          |                   |            |              |            |                    |             |             |              | nan                         | nan                                                 | pre_combat | 1     |  0.254  | 0.259 |         1     | {"GSE27155": 1.0}                | {"GSE27155": 1.0}                |
| V4     | nan        |         | BRS71_original        |              | GradientBoosting  |           |    0.841 |    0.993 |              |             |           |                   |          |                   |            |              |            |                    |             |             |              | nan                         | nan                                                 | pre_combat | 0.931 |  0.653  | 0.666 |         0.908 | {"GSE27155": 0.9079670329670331} | {"GSE27155": 0.9313186813186813} |
| V4     | nan        |         | BRS71_original        |              | XGBoost           |           |    0.695 |    0.973 |              |             |           |                   |          |                   |            |              |            |                    |             |             |              | nan                         | nan                                                 | pre_combat | 0.837 |  0.196  | 0.333 |         0.908 | {"GSE27155": 0.9079670329670331} | {"GSE27155": 0.8365384615384615} |
| V4     | nan        |         | v3_panel_compact      |              | LogReg_l2         |           |    0.986 |    1     |              |             |           |                   |          |                   |            |              |            |                    |             |             |              | nan                         | nan                                                 | pre_combat | 0.997 |  0.311  | 0.313 |         0.982 | {"GSE27155": 0.9821428571428572} | {"GSE27155": 0.9972527472527473} |
| V4     | nan        |         | v3_panel_compact      |              | LogReg_elasticnet |           |    1     |    1     |              |             |           |                   |          |                   |            |              |            |                    |             |             |              | nan                         | nan                                                 | pre_combat | 1     |  0.314  | 0.315 |         1     | {"GSE27155": 1.0}                | {"GSE27155": 1.0}                |
| V4     | nan        |         | v3_panel_compact      |              | RandomForest      |           |    1     |    1     |              |             |           |                   |          |                   |            |              |            |                    |             |             |              | nan                         | nan                                                 | pre_combat | 1     |  0.174  | 0.404 |         1     | {"GSE27155": 1.0}                | {"GSE27155": 1.0}                |
| V4     | nan        |         | v3_panel_compact      |              | GradientBoosting  |           |    0.5   |    0.5   |              |             |           |                   |          |                   |            |              |            |                    |             |             |              | nan                         | nan                                                 | pre_combat | 0.5   |  0.683  | 0.683 |         0.5   | {"GSE27155": 0.5}                | {"GSE27155": 0.5}                |

_(76 rows total — showing first 50)_

### ✅ `results/tables/v3_leakage_curve.tsv`
_Rows: 8_

|   k_removed |   mean_auc |   auc_lo |   auc_hi |   n_features | removed_genes                                                                                                                                                              |
|------------:|-----------:|---------:|---------:|-------------:|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|           0 |      0.992 |    0.984 |    0.999 |           55 | nan                                                                                                                                                                        |
|           3 |      0.991 |    0.981 |    0.999 |           52 | MET,KLK10,DUSP5                                                                                                                                                            |
|           5 |      0.99  |    0.98  |    0.998 |           50 | MET,KLK10,DUSP5,SLC5A8,DIO1                                                                                                                                                |
|          10 |      0.984 |    0.971 |    0.995 |           45 | MET,KLK10,DUSP5,SLC5A8,DIO1,TPO,DUSP6,LOX,DIO2,HLA-DRA                                                                                                                     |
|          15 |      0.981 |    0.966 |    0.993 |           40 | MET,KLK10,DUSP5,SLC5A8,DIO1,TPO,DUSP6,LOX,DIO2,HLA-DRA,CDKN2B,SLC26A4,FOXP3,DUSP4,TG                                                                                       |
|          20 |      0.979 |    0.964 |    0.99  |           35 | MET,KLK10,DUSP5,SLC5A8,DIO1,TPO,DUSP6,LOX,DIO2,HLA-DRA,CDKN2B,SLC26A4,FOXP3,DUSP4,TG,DUOX2,CD274,IYD,MMP9,PAX8                                                             |
|          25 |      0.96  |    0.937 |    0.98  |           30 | MET,KLK10,DUSP5,SLC5A8,DIO1,TPO,DUSP6,LOX,DIO2,HLA-DRA,CDKN2B,SLC26A4,FOXP3,DUSP4,TG,DUOX2,CD274,IYD,MMP9,PAX8,FOSL1,PIK3CA,DUOX1,CDKN2A,CTNNB1                            |
|          30 |      0.895 |    0.855 |    0.933 |           25 | MET,KLK10,DUSP5,SLC5A8,DIO1,TPO,DUSP6,LOX,DIO2,HLA-DRA,CDKN2B,SLC26A4,FOXP3,DUSP4,TG,DUOX2,CD274,IYD,MMP9,PAX8,FOSL1,PIK3CA,DUOX1,CDKN2A,CTNNB1,ETV4,SPRY4,FOXE1,MSH2,CDH1 |

### ❌ NOT YET RUN — `results/ml/v3_lodo.tsv`

### ✅ `results/tables/v3_brs71_proxy_vs_original.tsv`
_Rows: 142_

| gene       | in_original   | in_proxy   | shared   |
|:-----------|:--------------|:-----------|:---------|
| AC002401.4 | False         | True       | False    |
| AC004847.1 | False         | True       | False    |
| AC009549.1 | False         | True       | False    |
| ACTBL2     | False         | True       | False    |
| ADAMTS14   | False         | True       | False    |
| AL096865.1 | False         | True       | False    |
| AL137026.1 | False         | True       | False    |
| ALDH1A3    | True          | False      | False    |
| ANKRD18B   | True          | False      | False    |
| AP002358.1 | False         | True       | False    |
| APOE       | True          | False      | False    |
| ARSI       | False         | True       | False    |
| BCHE       | True          | False      | False    |
| BEND6      | False         | True       | False    |
| BNC1       | False         | True       | False    |
| BRINP2     | False         | True       | False    |
| C16ORF89   | False         | True       | False    |
| CA2        | True          | False      | False    |
| CDH1       | True          | False      | False    |
| CDH11      | True          | False      | False    |
| CDK5RAP2   | False         | True       | False    |
| CEACAM6    | False         | True       | False    |
| CFB        | False         | True       | False    |
| CITED1     | True          | False      | False    |
| CLDN1      | True          | False      | False    |
| CLDN10     | False         | True       | False    |
| CLDN16     | True          | False      | False    |
| COL9A3     | True          | False      | False    |
| CPE        | True          | False      | False    |
| CREB5      | False         | True       | False    |
| CRLF2      | False         | True       | False    |
| CRYBG2     | False         | True       | False    |
| CST2       | False         | True       | False    |
| CST5       | False         | True       | False    |
| CYP1B1-AS1 | False         | True       | False    |
| DCSTAMP    | False         | True       | False    |
| DIO1       | True          | False      | False    |
| DIO2       | True          | False      | False    |
| DLK1       | True          | False      | False    |
| DMBX1      | False         | True       | False    |
| DRAXIN     | False         | True       | False    |
| DSC3       | False         | True       | False    |
| DUOX1      | True          | False      | False    |
| DUOX2      | True          | False      | False    |
| DUOXA2     | True          | False      | False    |
| ELFN2      | False         | True       | False    |
| ERBB4      | True          | False      | False    |
| EREG       | False         | True       | False    |
| FAM111A-DT | False         | True       | False    |
| FAM155B    | False         | True       | False    |

_(142 rows total — showing first 50)_

### ✅ `results/tables/v3_bethesda_sim.tsv`
_Rows: 455_

|   threshold |   prev |   sens |   spec |   npv |   ppv |   unnecessary_surgery |   missed_cancer |   nb_model |   nb_treat_all |   nb_treat_none |   tp |   fp |   tn |   fn |
|------------:|-------:|-------:|-------:|------:|------:|----------------------:|----------------:|-----------:|---------------:|----------------:|-----:|-----:|-----:|-----:|
|        0.05 |    0.1 |  0.994 |  0.39  | 0.998 | 0.153 |                 0.549 |          0.0006 |    0.0705  |         0.0526 |               0 |  994 | 5494 | 3506 |    6 |
|        0.06 |    0.1 |  0.994 |  0.39  | 0.998 | 0.153 |                 0.549 |          0.0006 |    0.0644  |         0.0426 |               0 |  994 | 5488 | 3512 |    6 |
|        0.07 |    0.1 |  0.994 |  0.391 | 0.998 | 0.153 |                 0.548 |          0.0006 |    0.0581  |         0.0323 |               0 |  994 | 5484 | 3516 |    6 |
|        0.08 |    0.1 |  0.994 |  0.392 | 0.998 | 0.154 |                 0.547 |          0.0006 |    0.0518  |         0.0217 |               0 |  994 | 5471 | 3529 |    6 |
|        0.09 |    0.1 |  0.994 |  0.394 | 0.998 | 0.154 |                 0.546 |          0.0006 |    0.0454  |         0.011  |               0 |  994 | 5458 | 3542 |    6 |
|        0.1  |    0.1 |  0.994 |  0.395 | 0.998 | 0.154 |                 0.545 |          0.0006 |    0.0389  |         0      |               0 |  994 | 5448 | 3552 |    6 |
|        0.11 |    0.1 |  0.994 |  0.396 | 0.998 | 0.155 |                 0.543 |          0.0006 |    0.0322  |        -0.0112 |               0 |  994 | 5435 | 3565 |    6 |
|        0.12 |    0.1 |  0.994 |  0.397 | 0.998 | 0.155 |                 0.543 |          0.0006 |    0.0253  |        -0.0227 |               0 |  994 | 5431 | 3569 |    6 |
|        0.13 |    0.1 |  0.993 |  0.398 | 0.998 | 0.155 |                 0.542 |          0.0007 |    0.0183  |        -0.0345 |               0 |  993 | 5422 | 3578 |    7 |
|        0.14 |    0.1 |  0.993 |  0.399 | 0.998 | 0.155 |                 0.541 |          0.0007 |    0.0112  |        -0.0465 |               0 |  993 | 5411 | 3589 |    7 |
|        0.15 |    0.1 |  0.993 |  0.4   | 0.998 | 0.155 |                 0.54  |          0.0007 |    0.00399 |        -0.0588 |               0 |  993 | 5401 | 3599 |    7 |
|        0.16 |    0.1 |  0.993 |  0.401 | 0.998 | 0.156 |                 0.539 |          0.0007 |   -0.00339 |        -0.0714 |               0 |  993 | 5391 | 3609 |    7 |
|        0.17 |    0.1 |  0.993 |  0.403 | 0.998 | 0.156 |                 0.538 |          0.0007 |   -0.0108  |        -0.0843 |               0 |  993 | 5377 | 3623 |    7 |
|        0.18 |    0.1 |  0.993 |  0.403 | 0.998 | 0.156 |                 0.537 |          0.0007 |   -0.0186  |        -0.0976 |               0 |  993 | 5372 | 3628 |    7 |
|        0.19 |    0.1 |  0.993 |  0.404 | 0.998 | 0.156 |                 0.536 |          0.0007 |   -0.0265  |        -0.111  |               0 |  993 | 5362 | 3638 |    7 |
|        0.2  |    0.1 |  0.993 |  0.405 | 0.998 | 0.157 |                 0.535 |          0.0007 |   -0.0345  |        -0.125  |               0 |  993 | 5352 | 3648 |    7 |
|        0.21 |    0.1 |  0.993 |  0.406 | 0.998 | 0.157 |                 0.534 |          0.0007 |   -0.0427  |        -0.139  |               0 |  993 | 5343 | 3657 |    7 |
|        0.22 |    0.1 |  0.993 |  0.407 | 0.998 | 0.157 |                 0.533 |          0.0007 |   -0.0512  |        -0.154  |               0 |  993 | 5335 | 3665 |    7 |
|        0.23 |    0.1 |  0.993 |  0.408 | 0.998 | 0.157 |                 0.532 |          0.0007 |   -0.0598  |        -0.169  |               0 |  993 | 5325 | 3675 |    7 |
|        0.24 |    0.1 |  0.993 |  0.41  | 0.998 | 0.157 |                 0.531 |          0.0007 |   -0.0685  |        -0.184  |               0 |  993 | 5314 | 3686 |    7 |
|        0.25 |    0.1 |  0.992 |  0.41  | 0.998 | 0.157 |                 0.531 |          0.0008 |   -0.0777  |        -0.2    |               0 |  992 | 5307 | 3693 |    8 |
|        0.26 |    0.1 |  0.991 |  0.411 | 0.998 | 0.158 |                 0.53  |          0.0009 |   -0.087   |        -0.216  |               0 |  991 | 5298 | 3702 |    9 |
|        0.27 |    0.1 |  0.932 |  0.843 | 0.991 | 0.398 |                 0.141 |          0.0068 |    0.0411  |        -0.233  |               0 |  932 | 1409 | 7591 |   68 |
|        0.28 |    0.1 |  0.932 |  0.844 | 0.991 | 0.399 |                 0.14  |          0.0068 |    0.0386  |        -0.25   |               0 |  932 | 1403 | 7597 |   68 |
|        0.29 |    0.1 |  0.932 |  0.845 | 0.991 | 0.4   |                 0.14  |          0.0068 |    0.0361  |        -0.268  |               0 |  932 | 1397 | 7603 |   68 |
|        0.3  |    0.1 |  0.932 |  0.846 | 0.991 | 0.402 |                 0.139 |          0.0068 |    0.0337  |        -0.286  |               0 |  932 | 1389 | 7611 |   68 |
|        0.31 |    0.1 |  0.932 |  0.847 | 0.991 | 0.403 |                 0.138 |          0.0068 |    0.0312  |        -0.304  |               0 |  932 | 1381 | 7619 |   68 |
|        0.32 |    0.1 |  0.931 |  0.847 | 0.991 | 0.404 |                 0.137 |          0.0069 |    0.0285  |        -0.324  |               0 |  931 | 1373 | 7627 |   69 |
|        0.33 |    0.1 |  0.931 |  0.848 | 0.991 | 0.405 |                 0.137 |          0.0069 |    0.0258  |        -0.343  |               0 |  931 | 1367 | 7633 |   69 |
|        0.34 |    0.1 |  0.931 |  0.849 | 0.991 | 0.406 |                 0.136 |          0.0069 |    0.0229  |        -0.364  |               0 |  931 | 1363 | 7637 |   69 |
|        0.35 |    0.1 |  0.931 |  0.849 | 0.991 | 0.407 |                 0.136 |          0.0069 |    0.0199  |        -0.385  |               0 |  931 | 1359 | 7641 |   69 |
|        0.36 |    0.1 |  0.931 |  0.85  | 0.991 | 0.408 |                 0.135 |          0.0069 |    0.0172  |        -0.406  |               0 |  931 | 1349 | 7651 |   69 |
|        0.37 |    0.1 |  0.931 |  0.851 | 0.991 | 0.409 |                 0.134 |          0.0069 |    0.0142  |        -0.429  |               0 |  931 | 1343 | 7657 |   69 |
|        0.38 |    0.1 |  0.931 |  0.852 | 0.991 | 0.411 |                 0.134 |          0.0069 |    0.0113  |        -0.452  |               0 |  931 | 1335 | 7665 |   69 |
|        0.39 |    0.1 |  0.931 |  0.853 | 0.991 | 0.413 |                 0.132 |          0.0069 |    0.00851 |        -0.475  |               0 |  931 | 1323 | 7677 |   69 |
|        0.4  |    0.1 |  0.93  |  0.854 | 0.991 | 0.414 |                 0.132 |          0.007  |    0.00527 |        -0.5    |               0 |  930 | 1316 | 7684 |   70 |
|        0.41 |    0.1 |  0.929 |  0.855 | 0.991 | 0.415 |                 0.131 |          0.0071 |    0.00207 |        -0.525  |               0 |  929 | 1307 | 7693 |   71 |
|        0.42 |    0.1 |  0.928 |  0.855 | 0.991 | 0.416 |                 0.13  |          0.0072 |   -0.00141 |        -0.552  |               0 |  928 | 1301 | 7699 |   72 |
|        0.43 |    0.1 |  0.928 |  0.856 | 0.991 | 0.418 |                 0.129 |          0.0072 |   -0.00482 |        -0.579  |               0 |  928 | 1294 | 7706 |   72 |
|        0.44 |    0.1 |  0.927 |  0.857 | 0.991 | 0.418 |                 0.129 |          0.0073 |   -0.00858 |        -0.607  |               0 |  927 | 1289 | 7711 |   73 |
|        0.45 |    0.1 |  0.927 |  0.858 | 0.991 | 0.42  |                 0.128 |          0.0073 |   -0.0119  |        -0.636  |               0 |  927 | 1279 | 7721 |   73 |
|        0.46 |    0.1 |  0.927 |  0.859 | 0.991 | 0.421 |                 0.127 |          0.0073 |   -0.0157  |        -0.667  |               0 |  927 | 1273 | 7727 |   73 |
|        0.47 |    0.1 |  0.926 |  0.859 | 0.991 | 0.422 |                 0.127 |          0.0074 |   -0.0198  |        -0.698  |               0 |  926 | 1267 | 7733 |   74 |
|        0.48 |    0.1 |  0.926 |  0.86  | 0.991 | 0.424 |                 0.126 |          0.0074 |   -0.0237  |        -0.731  |               0 |  926 | 1260 | 7740 |   74 |
|        0.49 |    0.1 |  0.926 |  0.86  | 0.991 | 0.424 |                 0.126 |          0.0074 |   -0.0281  |        -0.765  |               0 |  926 | 1256 | 7744 |   74 |
|        0.5  |    0.1 |  0.926 |  0.861 | 0.991 | 0.426 |                 0.125 |          0.0074 |   -0.0324  |        -0.8    |               0 |  926 | 1250 | 7750 |   74 |
|        0.51 |    0.1 |  0.926 |  0.861 | 0.991 | 0.426 |                 0.125 |          0.0074 |   -0.0372  |        -0.837  |               0 |  926 | 1247 | 7753 |   74 |
|        0.52 |    0.1 |  0.925 |  0.862 | 0.99  | 0.427 |                 0.124 |          0.0075 |   -0.042   |        -0.875  |               0 |  925 | 1242 | 7758 |   75 |
|        0.53 |    0.1 |  0.925 |  0.862 | 0.99  | 0.428 |                 0.124 |          0.0075 |   -0.0471  |        -0.915  |               0 |  925 | 1238 | 7762 |   75 |
|        0.54 |    0.1 |  0.925 |  0.863 | 0.99  | 0.429 |                 0.123 |          0.0075 |   -0.0522  |        -0.957  |               0 |  925 | 1233 | 7767 |   75 |

_(455 rows total — showing first 50)_

### ✅ `results/tables/v3_survival_cox.tsv`
_Rows: 13_

| covariate                  | mode         |   n |    hr |     lo |     hi |        p | note   |
|:---------------------------|:-------------|----:|------:|-------:|-------:|---------:|:-------|
| tds_tertile_low            | univariate   | 504 | 0.856 | 0.308  |  2.38  | 0.766    |        |
| tds_tertile_mid            | univariate   | 504 | 0.12  | 0.0148 |  0.976 | 0.0474   |        |
| dediff_tertile_mid         | univariate   | 504 | 0.12  | 0.0148 |  0.976 | 0.0474   |        |
| dediff_tertile_high        | univariate   | 504 | 0.856 | 0.308  |  2.38  | 0.766    |        |
| molecular_subtype_RAS_like | univariate   | 504 | 0.296 | 0.0388 |  2.25  | 0.24     |        |
| molecular_subtype_unknown  | univariate   | 504 | 3.57  | 0.462  | 27.5   | 0.223    |        |
| stage_num                  | univariate   | 460 | 2.53  | 1.51   |  4.22  | 0.000413 |        |
| sex_num                    | univariate   | 504 | 1.95  | 0.704  |  5.38  | 0.2      |        |
| age_at_diagnosis           | univariate   | 492 | 1.18  | 1.11   |  1.25  | 1.12e-07 |        |
| subtype_RAS_like           | multivariate | 456 | 0.691 | 0.0791 |  6.03  | 0.738    |        |
| stage_num                  | multivariate | 456 | 1.28  | 0.636  |  2.58  | 0.488    |        |
| age_at_diagnosis           | multivariate | 456 | 1.17  | 1.1    |  1.25  | 1.89e-06 |        |
| sex_num                    | multivariate | 456 | 1.52  | 0.439  |  5.25  | 0.51     |        |

### ✅ `results/tables/v3_ext_batch_diagnostics.json`
```json
{
  "d1": {
    "status": "ok",
    "cohorts": [
      "GSE126698",
      "GSE27155",
      "GSE29265",
      "GSE33630",
      "GSE76039",
      "TCGA-THCA"
    ],
    "macro_auc": 1.0,
    "flag_identifiable": true,
    "n_shared_genes": 11731
  },
  "d2": {
    "status": "ok",
    "n_samples": 890,
    "n_shared_genes": 11731
  },
  "d3": {
    "status": "ok",
    "tds16_std_pre": 3.586282902473211,
    "tds16_std_post": 3.063863884889448,
    "delta_std": 0.5224190175837631
  },
  "d4": {
    "status": "ok",
    "n_genes": 13
  }
}
```

### ✅ `results/tables/v3_three_class_metrics.tsv`
_Rows: 8_

| model   | class_name       |   n |   ovr_auc |   precision |   recall |    f1 | small_n_flag   |
|:--------|:-----------------|----:|----------:|------------:|---------:|------:|:---------------|
| logreg  | BRAF_like        | 747 |     0.917 |       0.956 |    0.929 | 0.942 |                |
| logreg  | RAS_like         | 117 |     0.965 |       0.673 |    0.915 | 0.775 |                |
| logreg  | dedifferentiated |  68 |     0.93  |       0.979 |    0.676 | 0.8   |                |
| logreg  | MACRO            | 932 |     0.937 |             |          |       |                |
| gb      | BRAF_like        | 747 |     0.911 |       0.938 |    0.971 | 0.954 |                |
| gb      | RAS_like         | 117 |     0.971 |       0.814 |    0.786 | 0.8   |                |
| gb      | dedifferentiated |  68 |     0.931 |       0.978 |    0.662 | 0.789 |                |
| gb      | MACRO            | 932 |     0.938 |             |          |       |                |

### ✅ `results/tables/v3_three_class_shap.tsv`
_Rows: 45_

| class_name       |   rank | gene    |   mean_abs_shap |
|:-----------------|-------:|:--------|----------------:|
| BRAF_like        |      1 | TPO     |         0.0313  |
| BRAF_like        |      2 | DIO1    |         0.0303  |
| BRAF_like        |      3 | DIO2    |         0.0286  |
| BRAF_like        |      4 | MET     |         0.0202  |
| BRAF_like        |      5 | SLC5A8  |         0.016   |
| BRAF_like        |      6 | SLC26A4 |         0.0125  |
| BRAF_like        |      7 | TWIST1  |         0.0102  |
| BRAF_like        |      8 | IYD     |         0.0101  |
| BRAF_like        |      9 | DUOX2   |         0.00969 |
| BRAF_like        |     10 | DUOX1   |         0.00934 |
| BRAF_like        |     11 | KLK10   |         0.00904 |
| BRAF_like        |     12 | VIM     |         0.00665 |
| BRAF_like        |     13 | TG      |         0.00605 |
| BRAF_like        |     14 | MSH2    |         0.00596 |
| BRAF_like        |     15 | AKT1    |         0.00596 |
| RAS_like         |      1 | TPO     |         0.0334  |
| RAS_like         |      2 | DIO1    |         0.0305  |
| RAS_like         |      3 | DIO2    |         0.0276  |
| RAS_like         |      4 | SLC5A8  |         0.0182  |
| RAS_like         |      5 | MET     |         0.0168  |
| RAS_like         |      6 | SLC26A4 |         0.0129  |
| RAS_like         |      7 | DUOX2   |         0.0112  |
| RAS_like         |      8 | IYD     |         0.011   |
| RAS_like         |      9 | DUOX1   |         0.00874 |
| RAS_like         |     10 | KLK10   |         0.00742 |
| RAS_like         |     11 | TG      |         0.00676 |
| RAS_like         |     12 | PAX8    |         0.0062  |
| RAS_like         |     13 | LOX     |         0.00437 |
| RAS_like         |     14 | CDH1    |         0.00384 |
| RAS_like         |     15 | THADA   |         0.00383 |
| dedifferentiated |      1 | TWIST1  |         0.0114  |
| dedifferentiated |      2 | CTNNB1  |         0.00704 |
| dedifferentiated |      3 | VIM     |         0.00671 |
| dedifferentiated |      4 | MSH2    |         0.00517 |
| dedifferentiated |      5 | SNAI2   |         0.00498 |
| dedifferentiated |      6 | PTEN    |         0.00495 |
| dedifferentiated |      7 | APC     |         0.00474 |
| dedifferentiated |      8 | MET     |         0.00423 |
| dedifferentiated |      9 | AKT1    |         0.00359 |
| dedifferentiated |     10 | CDKN2A  |         0.0035  |
| dedifferentiated |     11 | THADA   |         0.00345 |
| dedifferentiated |     12 | CDH1    |         0.00303 |
| dedifferentiated |     13 | FOXE1   |         0.00252 |
| dedifferentiated |     14 | SLC5A8  |         0.00226 |
| dedifferentiated |     15 | TSHR    |         0.00223 |

### ❌ NOT YET RUN — `results/tables/v3_survival_results.tsv`

### ✅ `results/ml/v3_ext_permutation_null.tsv`
_Rows: 1_

|   observed_auc |   n_iter |   null_mean |   null_std |   p_value |   null_q97.5 |
|---------------:|---------:|------------:|-----------:|----------:|-------------:|
|              1 |     1000 |       0.469 |      0.168 |         0 |         0.78 |

### ✅ `results/tables/v3_bethesda_surgery_reduction.tsv`
_Rows: 5_

|   prev |   chosen_threshold |   sens |   spec |   npv |   ppv |   model_surgery_rate |   reduction_vs_treat_all |
|-------:|-------------------:|-------:|-------:|------:|------:|---------------------:|-------------------------:|
|   0.1  |               0.26 |  0.991 |  0.411 | 0.998 | 0.158 |                0.629 |                    0.371 |
|   0.15 |               0.26 |  0.987 |  0.406 | 0.995 | 0.227 |                0.653 |                    0.347 |
|   0.2  |               0.26 |  0.995 |  0.408 | 0.997 | 0.296 |                0.673 |                    0.327 |
|   0.25 |               0.26 |  0.992 |  0.397 | 0.993 | 0.354 |                0.7   |                    0.3   |
|   0.3  |               0.26 |  0.992 |  0.415 | 0.991 | 0.421 |                0.707 |                    0.293 |

### ✅ `results/tables/v3_three_class_confusion.tsv`
_Rows: 18_

| model   | true             | pred             |   count |
|:--------|:-----------------|:-----------------|--------:|
| logreg  | BRAF_like        | BRAF_like        |     694 |
| logreg  | BRAF_like        | RAS_like         |      52 |
| logreg  | BRAF_like        | dedifferentiated |       1 |
| logreg  | RAS_like         | BRAF_like        |      10 |
| logreg  | RAS_like         | RAS_like         |     107 |
| logreg  | RAS_like         | dedifferentiated |       0 |
| logreg  | dedifferentiated | BRAF_like        |      22 |
| logreg  | dedifferentiated | RAS_like         |       0 |
| logreg  | dedifferentiated | dedifferentiated |      46 |
| gb      | BRAF_like        | BRAF_like        |     725 |
| gb      | BRAF_like        | RAS_like         |      21 |
| gb      | BRAF_like        | dedifferentiated |       1 |
| gb      | RAS_like         | BRAF_like        |      25 |
| gb      | RAS_like         | RAS_like         |      92 |
| gb      | RAS_like         | dedifferentiated |       0 |
| gb      | dedifferentiated | BRAF_like        |      23 |
| gb      | dedifferentiated | RAS_like         |       0 |
| gb      | dedifferentiated | dedifferentiated |      45 |

### ✅ `results/tables/v3_decision_curve_gse27155.tsv`
_Rows: 46_

|   threshold |   tp |   fp |   tn |   fn |   sens |   spec |   nb_model |   nb_treat_all |   nb_treat_none |
|------------:|-----:|-----:|-----:|-----:|-------:|-------:|-----------:|---------------:|----------------:|
|        0.05 |   45 |   18 |    2 |    1 |  0.978 |   0.1  |      0.667 |          0.681 |               0 |
|        0.06 |   45 |   18 |    2 |    1 |  0.978 |   0.1  |      0.664 |          0.678 |               0 |
|        0.07 |   45 |   18 |    2 |    1 |  0.978 |   0.1  |      0.661 |          0.674 |               0 |
|        0.08 |   45 |   18 |    2 |    1 |  0.978 |   0.1  |      0.658 |          0.671 |               0 |
|        0.09 |   45 |   18 |    2 |    1 |  0.978 |   0.1  |      0.655 |          0.667 |               0 |
|        0.1  |   45 |   18 |    2 |    1 |  0.978 |   0.1  |      0.652 |          0.663 |               0 |
|        0.11 |   45 |   18 |    2 |    1 |  0.978 |   0.1  |      0.648 |          0.66  |               0 |
|        0.12 |   45 |   18 |    2 |    1 |  0.978 |   0.1  |      0.645 |          0.656 |               0 |
|        0.13 |   45 |   18 |    2 |    1 |  0.978 |   0.1  |      0.641 |          0.652 |               0 |
|        0.14 |   45 |   18 |    2 |    1 |  0.978 |   0.1  |      0.637 |          0.648 |               0 |
|        0.15 |   45 |   18 |    2 |    1 |  0.978 |   0.1  |      0.634 |          0.643 |               0 |
|        0.16 |   45 |   18 |    2 |    1 |  0.978 |   0.1  |      0.63  |          0.639 |               0 |
|        0.17 |   45 |   15 |    5 |    1 |  0.978 |   0.25 |      0.635 |          0.635 |               0 |
|        0.18 |   45 |   15 |    5 |    1 |  0.978 |   0.25 |      0.632 |          0.63  |               0 |
|        0.19 |   45 |   15 |    5 |    1 |  0.978 |   0.25 |      0.629 |          0.626 |               0 |
|        0.2  |   45 |   15 |    5 |    1 |  0.978 |   0.25 |      0.625 |          0.621 |               0 |
|        0.21 |   45 |   13 |    7 |    1 |  0.978 |   0.35 |      0.629 |          0.616 |               0 |
|        0.22 |   45 |   13 |    7 |    1 |  0.978 |   0.35 |      0.626 |          0.611 |               0 |
|        0.23 |   45 |   12 |    8 |    1 |  0.978 |   0.4  |      0.628 |          0.606 |               0 |
|        0.24 |   45 |   12 |    8 |    1 |  0.978 |   0.4  |      0.624 |          0.601 |               0 |
|        0.25 |   45 |   12 |    8 |    1 |  0.978 |   0.4  |      0.621 |          0.596 |               0 |
|        0.26 |   44 |   12 |    8 |    2 |  0.957 |   0.4  |      0.603 |          0.59  |               0 |
|        0.27 |   44 |   12 |    8 |    2 |  0.957 |   0.4  |      0.599 |          0.585 |               0 |
|        0.28 |   44 |   12 |    8 |    2 |  0.957 |   0.4  |      0.596 |          0.579 |               0 |
|        0.29 |   44 |   12 |    8 |    2 |  0.957 |   0.4  |      0.592 |          0.573 |               0 |
|        0.3  |   44 |   12 |    8 |    2 |  0.957 |   0.4  |      0.589 |          0.567 |               0 |
|        0.31 |   44 |   12 |    8 |    2 |  0.957 |   0.4  |      0.585 |          0.561 |               0 |
|        0.32 |   44 |   12 |    8 |    2 |  0.957 |   0.4  |      0.581 |          0.554 |               0 |
|        0.33 |   44 |   12 |    8 |    2 |  0.957 |   0.4  |      0.577 |          0.548 |               0 |
|        0.34 |   43 |   12 |    8 |    3 |  0.935 |   0.4  |      0.558 |          0.541 |               0 |
|        0.35 |   43 |   12 |    8 |    3 |  0.935 |   0.4  |      0.554 |          0.534 |               0 |
|        0.36 |   43 |   12 |    8 |    3 |  0.935 |   0.4  |      0.549 |          0.527 |               0 |
|        0.37 |   43 |   12 |    8 |    3 |  0.935 |   0.4  |      0.545 |          0.519 |               0 |
|        0.38 |   43 |   12 |    8 |    3 |  0.935 |   0.4  |      0.54  |          0.511 |               0 |
|        0.39 |   43 |   12 |    8 |    3 |  0.935 |   0.4  |      0.535 |          0.503 |               0 |
|        0.4  |   43 |   12 |    8 |    3 |  0.935 |   0.4  |      0.53  |          0.495 |               0 |
|        0.41 |   43 |   12 |    8 |    3 |  0.935 |   0.4  |      0.525 |          0.486 |               0 |
|        0.42 |   43 |   12 |    8 |    3 |  0.935 |   0.4  |      0.52  |          0.478 |               0 |
|        0.43 |   43 |   12 |    8 |    3 |  0.935 |   0.4  |      0.514 |          0.468 |               0 |
|        0.44 |   43 |   12 |    8 |    3 |  0.935 |   0.4  |      0.509 |          0.459 |               0 |
|        0.45 |   42 |   12 |    8 |    4 |  0.913 |   0.4  |      0.488 |          0.449 |               0 |
|        0.46 |   42 |   12 |    8 |    4 |  0.913 |   0.4  |      0.481 |          0.439 |               0 |
|        0.47 |   42 |   12 |    8 |    4 |  0.913 |   0.4  |      0.475 |          0.428 |               0 |
|        0.48 |   42 |   12 |    8 |    4 |  0.913 |   0.4  |      0.469 |          0.417 |               0 |
|        0.49 |   42 |   12 |    8 |    4 |  0.913 |   0.4  |      0.462 |          0.406 |               0 |
|        0.5  |   36 |   10 |   10 |   10 |  0.783 |   0.5  |      0.394 |          0.394 |               0 |

### ❌ NOT YET RUN — `results/ml/v3_three_class_results.tsv`

### ✅ `results/tables/v3_calibration_gse27155.tsv`
_Rows: 26_

| mode     |   bin_mean_pred |   bin_frac_pos |   brier |   ece |   n_samples |
|:---------|----------------:|---------------:|--------:|------:|------------:|
| uncal    |           0.941 |          0.143 |   0.284 | 0.285 |          66 |
| uncal    |           0.969 |          0.429 |   0.284 | 0.285 |          66 |
| uncal    |           0.976 |          0.5   |   0.284 | 0.285 |          66 |
| uncal    |           0.98  |          0.714 |   0.284 | 0.285 |          66 |
| uncal    |           0.984 |          0.5   |   0.284 | 0.285 |          66 |
| uncal    |           0.987 |          0.857 |   0.284 | 0.285 |          66 |
| uncal    |           0.99  |          1     |   0.284 | 0.285 |          66 |
| uncal    |           0.991 |          0.857 |   0.284 | 0.285 |          66 |
| uncal    |           0.994 |          1     |   0.284 | 0.285 |          66 |
| uncal    |           0.996 |          1     |   0.284 | 0.285 |          66 |
| isotonic |           0.112 |          0.125 |   0.161 | 0.147 |          66 |
| isotonic |           0.439 |          0.75  |   0.161 | 0.147 |          66 |
| isotonic |           0.587 |          0.333 |   0.161 | 0.147 |          66 |
| isotonic |           0.707 |          0.5   |   0.161 | 0.147 |          66 |
| isotonic |           0.881 |          0.929 |   0.161 | 0.147 |          66 |
| isotonic |           0.994 |          0.947 |   0.161 | 0.147 |          66 |
| platt    |           0.692 |          0.714 |   0.211 | 0.195 |          66 |
| platt    |           0.693 |          0.714 |   0.211 | 0.195 |          66 |
| platt    |           0.697 |          0.167 |   0.211 | 0.195 |          66 |
| platt    |           0.698 |          0.429 |   0.211 | 0.195 |          66 |
| platt    |           0.698 |          0.5   |   0.211 | 0.195 |          66 |
| platt    |           0.698 |          0.714 |   0.211 | 0.195 |          66 |
| platt    |           0.698 |          0.833 |   0.211 | 0.195 |          66 |
| platt    |           0.699 |          0.857 |   0.211 | 0.195 |          66 |
| platt    |           0.699 |          1     |   0.211 | 0.195 |          66 |
| platt    |           0.699 |          1     |   0.211 | 0.195 |          66 |

### ✅ `results/tables/v3_panel_size_curve.tsv`
_Rows: 18_

| strategy     |   k |   tcga_mean_auc |   tcga_auc_lo |   tcga_auc_hi |   gse27155_mean_auc |   gse126698_mean_auc |   gse27155_mean_npv_se95 |   gse126698_mean_npv_se95 |   n_tcga |   n_gse27155 |   n_gse126698 |
|:-------------|----:|----------------:|--------------:|--------------:|--------------------:|---------------------:|-------------------------:|--------------------------:|---------:|-------------:|--------------:|
| univariate_d |   4 |           0.918 |         0.879 |         0.955 |               0.909 |                0.875 |                        1 |                         1 |      503 |           66 |            12 |
| univariate_d |   8 |           0.927 |         0.901 |         0.953 |               0.899 |                0.487 |                        1 |                         1 |      503 |           66 |            12 |
| univariate_d |  16 |           0.932 |         0.898 |         0.961 |               0.929 |                0.487 |                        1 |                         1 |      503 |           66 |            12 |
| univariate_d |  24 |           0.936 |         0.911 |         0.961 |               0.902 |                0.688 |                        1 |                         1 |      503 |           66 |            12 |
| univariate_d |  32 |           0.928 |         0.903 |         0.96  |               0.863 |                0.625 |                        1 |                         1 |      503 |           66 |            12 |
| univariate_d |  40 |           0.923 |         0.896 |         0.955 |               0.839 |                0.662 |                        1 |                           |      503 |           66 |            12 |
| tds16_union  |   4 |           0.903 |         0.871 |         0.935 |               0.857 |                0.55  |                        1 |                         1 |      503 |           66 |            12 |
| tds16_union  |   8 |           0.91  |         0.89  |         0.933 |               0.877 |                0.512 |                          |                         1 |      503 |           66 |            12 |
| tds16_union  |  16 |           0.928 |         0.908 |         0.949 |               0.818 |                0.688 |                        1 |                         1 |      503 |           66 |            12 |
| tds16_union  |  24 |           0.939 |         0.909 |         0.964 |               0.899 |                0.562 |                        1 |                         1 |      503 |           66 |            12 |
| tds16_union  |  32 |           0.934 |         0.91  |         0.96  |               0.844 |                0.725 |                        1 |                         1 |      503 |           66 |            12 |
| tds16_union  |  40 |           0.922 |         0.895 |         0.955 |               0.844 |                0.662 |                        1 |                           |      503 |           66 |            12 |
| qubo_neal    |   4 |           0.935 |         0.908 |         0.957 |               0.874 |                0.575 |                        1 |                         1 |      503 |           66 |            12 |
| qubo_neal    |   8 |           0.93  |         0.894 |         0.965 |               0.826 |                0.738 |                        1 |                         1 |      503 |           66 |            12 |
| qubo_neal    |  16 |           0.921 |         0.888 |         0.961 |               0.866 |                0.662 |                        1 |                         1 |      503 |           66 |            12 |
| qubo_neal    |  24 |           0.922 |         0.889 |         0.958 |               0.848 |                0.787 |                        1 |                         1 |      503 |           66 |            12 |
| qubo_neal    |  32 |           0.928 |         0.903 |         0.958 |               0.843 |                0.812 |                          |                         1 |      503 |           66 |            12 |
| qubo_neal    |  40 |           0.928 |         0.9   |         0.96  |               0.801 |                0.525 |                        1 |                           |      503 |           66 |            12 |

### ✅ `results/tables/v3_mapk_ablation.tsv`
_Rows: 2_

| variant            |   n_features |   mean_cv_auc |   cv_auc_lo |   cv_auc_hi |
|:-------------------|-------------:|--------------:|------------:|------------:|
| full_TierA67_clean |           55 |         0.923 |       0.902 |       0.957 |
| minus_MAPK_OUTPUT  |           45 |         0.924 |       0.899 |       0.954 |

### ✅ `results/tables/v3_methylation_results.tsv`
_Rows: 1_

| task                        | label           |   n_samples |   n_positive |   mean_cv_auc | note                                   |
|:----------------------------|:----------------|------------:|-------------:|--------------:|:---------------------------------------|
| methylation_within_gse97466 | tumor_vs_normal |         141 |           74 |         0.986 | ; LODO infeasible: no TCGA methylation |

### ❌ NOT YET RUN — `results/ml/v3_dataset_identifiability.tsv`

### ✅ `results/tables/v3_lodo.tsv`
_Rows: 6_

| external   | mode                   |   n_test |   auc |   pr_auc |   bacc |   npv_at_se95 |   brier |   youden | note                |
|:-----------|:-----------------------|---------:|------:|---------:|-------:|--------------:|--------:|---------:|:--------------------|
| GSE27155   | train_tcga_only        |       66 | 0.843 |    0.929 |  0.798 |         0.889 |   0.284 |    0.985 | nan                 |
| GSE27155   | train_tcga_plus_others |       66 | 0.883 |    0.951 |  0.815 |         0.846 |   0.263 |    0.975 | nan                 |
| GSE126698  | train_tcga_only        |       12 | 0.694 |    0.809 |  0.75  |               |   0.45  |    0.986 | nan                 |
| GSE126698  | train_tcga_plus_others |       12 | 0.778 |    0.862 |  0.833 |         1     |   0.281 |    0.879 | ⚠ small n           |
| GSE213647  | train_tcga_only        |      349 |       |          |        |               |         |          | skipped: <2 classes |
| GSE213647  | train_tcga_plus_others |      349 |       |          |        |               |         |          | skipped: <2 classes |

### ❌ NOT YET RUN — `results/ml/v3_calibration_gse27155.tsv`

### ✅ `results/tables/v3_permutation_null_samples.tsv`
_Rows: 1000_

|   iteration |   null_auc |
|------------:|-----------:|
|           0 |      0.742 |
|           1 |      0.727 |
|           2 |      0.713 |
|           3 |      0.718 |
|           4 |      0.727 |
|           5 |      0.752 |
|           6 |      0.732 |
|           7 |      0.714 |
|           8 |      0.748 |
|           9 |      0.694 |
|          10 |      0.718 |
|          11 |      0.702 |
|          12 |      0.715 |
|          13 |      0.703 |
|          14 |      0.679 |
|          15 |      0.733 |
|          16 |      0.701 |
|          17 |      0.704 |
|          18 |      0.727 |
|          19 |      0.719 |
|          20 |      0.747 |
|          21 |      0.72  |
|          22 |      0.745 |
|          23 |      0.734 |
|          24 |      0.748 |
|          25 |      0.726 |
|          26 |      0.727 |
|          27 |      0.735 |
|          28 |      0.735 |
|          29 |      0.755 |
|          30 |      0.725 |
|          31 |      0.725 |
|          32 |      0.715 |
|          33 |      0.697 |
|          34 |      0.727 |
|          35 |      0.697 |
|          36 |      0.703 |
|          37 |      0.74  |
|          38 |      0.722 |
|          39 |      0.7   |
|          40 |      0.736 |
|          41 |      0.721 |
|          42 |      0.716 |
|          43 |      0.734 |
|          44 |      0.711 |
|          45 |      0.708 |
|          46 |      0.73  |
|          47 |      0.703 |
|          48 |      0.72  |
|          49 |      0.715 |

_(1000 rows total — showing first 50)_

### ✅ `results/ml/v3_ext_lodo_v3.tsv`
_Rows: 50_

| held_out   | feature_set           | model             |   n_train |   n_test |   auc |   bACC_Youden |   brier |    ece |
|:-----------|:----------------------|:------------------|----------:|---------:|------:|--------------:|--------:|-------:|
| GSE27155   | TDS16                 | LogReg_l2         |       351 |       41 | 0.997 |         0.982 |  0.317  | 0.317  |
| GSE27155   | TDS16                 | LogReg_elasticnet |       351 |       41 | 0.997 |         0.982 |  0.317  | 0.317  |
| GSE27155   | TDS16                 | RandomForest      |       351 |       41 | 0.828 |         0.828 |  0.271  | 0.243  |
| GSE27155   | TDS16                 | GradientBoosting  |       351 |       41 | 0.5   |         0.5   |  0.317  | 0.317  |
| GSE27155   | TDS16                 | XGBoost           |       351 |       41 | 0.861 |         0.885 |  0.316  | 0.316  |
| GSE27155   | TierA67_clean         | LogReg_l2         |       351 |       41 | 0.997 |         0.982 |  0.682  | 0.682  |
| GSE27155   | TierA67_clean         | LogReg_elasticnet |       351 |       41 | 0.964 |         0.964 |  0.683  | 0.683  |
| GSE27155   | TierA67_clean         | RandomForest      |       351 |       41 | 0.902 |         0.857 |  0.27   | 0.31   |
| GSE27155   | TierA67_clean         | GradientBoosting  |       351 |       41 | 0.592 |         0.727 |  0.448  | 0.491  |
| GSE27155   | TierA67_clean         | XGBoost           |       351 |       41 | 0.949 |         0.908 |  0.473  | 0.538  |
| GSE27155   | TierA67_clean_no_MAPK | LogReg_l2         |       351 |       41 | 0.995 |         0.964 |  0.127  | 0.266  |
| GSE27155   | TierA67_clean_no_MAPK | LogReg_elasticnet |       351 |       41 | 1     |         1     |  0.68   | 0.681  |
| GSE27155   | TierA67_clean_no_MAPK | RandomForest      |       351 |       41 | 0.933 |         0.905 |  0.24   | 0.201  |
| GSE27155   | TierA67_clean_no_MAPK | GradientBoosting  |       351 |       41 | 0.87  |         0.87  |  0.288  | 0.379  |
| GSE27155   | TierA67_clean_no_MAPK | XGBoost           |       351 |       41 | 0.995 |         0.964 |  0.45   | 0.548  |
| GSE27155   | BRS71_original        | LogReg_l2         |       351 |       41 | 1     |         1     |  0.213  | 0.256  |
| GSE27155   | BRS71_original        | LogReg_elasticnet |       351 |       41 | 1     |         1     |  0.0725 | 0.166  |
| GSE27155   | BRS71_original        | RandomForest      |       351 |       41 | 1     |         1     |  0.254  | 0.259  |
| GSE27155   | BRS71_original        | GradientBoosting  |       351 |       41 | 0.931 |         0.908 |  0.653  | 0.666  |
| GSE27155   | BRS71_original        | XGBoost           |       351 |       41 | 0.837 |         0.908 |  0.196  | 0.333  |
| GSE27155   | v3_panel_compact      | LogReg_l2         |       351 |       41 | 0.997 |         0.982 |  0.311  | 0.313  |
| GSE27155   | v3_panel_compact      | LogReg_elasticnet |       351 |       41 | 1     |         1     |  0.314  | 0.315  |
| GSE27155   | v3_panel_compact      | RandomForest      |       351 |       41 | 1     |         1     |  0.174  | 0.404  |
| GSE27155   | v3_panel_compact      | GradientBoosting  |       351 |       41 | 0.5   |         0.5   |  0.683  | 0.683  |
| GSE27155   | v3_panel_compact      | XGBoost           |       351 |       41 | 0.5   |         0.5   |  0.448  | 0.481  |
| TCGA-THCA  | TDS16                 | LogReg_l2         |        41 |      351 | 0.897 |         0.865 |  0.835  | 0.835  |
| TCGA-THCA  | TDS16                 | LogReg_elasticnet |        41 |      351 | 0.907 |         0.877 |  0.834  | 0.834  |
| TCGA-THCA  | TDS16                 | RandomForest      |        41 |      351 | 0.82  |         0.822 |  0.743  | 0.784  |
| TCGA-THCA  | TDS16                 | GradientBoosting  |        41 |      351 | 0.821 |         0.821 |  0.835  | 0.835  |
| TCGA-THCA  | TDS16                 | XGBoost           |        41 |      351 | 0.841 |         0.831 |  0.721  | 0.769  |
| TCGA-THCA  | TierA67_clean         | LogReg_l2         |        41 |      351 | 0.967 |         0.964 |  0.832  | 0.833  |
| TCGA-THCA  | TierA67_clean         | LogReg_elasticnet |        41 |      351 | 0.962 |         0.954 |  0.392  | 0.435  |
| TCGA-THCA  | TierA67_clean         | RandomForest      |        41 |      351 | 0.975 |         0.923 |  0.296  | 0.418  |
| TCGA-THCA  | TierA67_clean         | GradientBoosting  |        41 |      351 | 0.5   |         0.5   |  0.165  | 0.165  |
| TCGA-THCA  | TierA67_clean         | XGBoost           |        41 |      351 | 0.788 |         0.79  |  0.169  | 0.223  |
| TCGA-THCA  | TierA67_clean_no_MAPK | LogReg_l2         |        41 |      351 | 0.964 |         0.955 |  0.835  | 0.835  |
| TCGA-THCA  | TierA67_clean_no_MAPK | LogReg_elasticnet |        41 |      351 | 0.956 |         0.942 |  0.831  | 0.832  |
| TCGA-THCA  | TierA67_clean_no_MAPK | RandomForest      |        41 |      351 | 0.976 |         0.93  |  0.382  | 0.508  |
| TCGA-THCA  | TierA67_clean_no_MAPK | GradientBoosting  |        41 |      351 | 0.5   |         0.5   |  0.165  | 0.165  |
| TCGA-THCA  | TierA67_clean_no_MAPK | XGBoost           |        41 |      351 | 0.791 |         0.79  |  0.19   | 0.265  |
| TCGA-THCA  | BRS71_original        | LogReg_l2         |        41 |      351 | 0.973 |         0.974 |  0.0458 | 0.0483 |
| TCGA-THCA  | BRS71_original        | LogReg_elasticnet |        41 |      351 | 0.967 |         0.974 |  0.0596 | 0.0629 |
| TCGA-THCA  | BRS71_original        | RandomForest      |        41 |      351 | 0.975 |         0.954 |  0.218  | 0.395  |
| TCGA-THCA  | BRS71_original        | GradientBoosting  |        41 |      351 | 0.517 |         0.517 |  0.16   | 0.16   |
| TCGA-THCA  | BRS71_original        | XGBoost           |        41 |      351 | 0.517 |         0.517 |  0.149  | 0.123  |
| TCGA-THCA  | v3_panel_compact      | LogReg_l2         |        41 |      351 | 0.964 |         0.971 |  0.177  | 0.194  |
| TCGA-THCA  | v3_panel_compact      | LogReg_elasticnet |        41 |      351 | 0.974 |         0.98  |  0.0619 | 0.0744 |
| TCGA-THCA  | v3_panel_compact      | RandomForest      |        41 |      351 | 0.971 |         0.976 |  0.169  | 0.339  |
| TCGA-THCA  | v3_panel_compact      | GradientBoosting  |        41 |      351 | 0.909 |         0.909 |  0.0598 | 0.0598 |
| TCGA-THCA  | v3_panel_compact      | XGBoost           |        41 |      351 | 0.933 |         0.933 |  0.049  | 0.0404 |

### ✅ `results/tables/v3_dataset_identifiability.tsv`
_Rows: 4_

| dataset   |   n_samples |   ovr_auc_mean |   ovr_auc_std | ovr_auc_folds                 |
|:----------|------------:|---------------:|--------------:|:------------------------------|
| TCGA-THCA |         513 |              1 |             0 | 1.000,1.000,1.000,1.000,1.000 |
| GSE27155  |          66 |              1 |             0 | 1.000,1.000,1.000,1.000,1.000 |
| GSE76039  |          37 |              1 |             0 | 1.000,1.000,1.000,1.000,1.000 |
| GSE126698 |          22 |              1 |             0 | 1.000,1.000,1.000,1.000,1.000 |

### ❌ NOT YET RUN — `results/tables/v3_bethesda_triage.tsv`

### ✅ `results/ml/v3_ext_calibration.tsv`
_Rows: 25_

| task   | feature_set           | model             | y                                                                                 | p_pre                                                                                                                                                                                                                                                                                          | p_post                                                                                                                                                                                                                                                                                         |
|:-------|:----------------------|:------------------|:----------------------------------------------------------------------------------|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| V1     | TDS16                 | LogReg_l2         | 0;0;0;0;0;0;0;0;0;0;0;0;0;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1 | 1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000 | 1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000 |
| V1     | TDS16                 | LogReg_elasticnet | 0;0;0;0;0;0;0;0;0;0;0;0;0;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1 | 1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000 | 1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000 |
| V1     | TDS16                 | RandomForest      | 0;0;0;0;0;0;0;0;0;0;0;0;0;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1 | 0.9125;0.9250;0.9300;0.9200;0.9350;0.9125;0.9200;0.9050;0.9100;0.9075;0.9125;0.9300;0.9125;0.9200;0.9300;0.9325;0.9300;0.9250;0.9300;0.9475;0.9250;0.9250;0.9275;0.9325;0.9300;0.9300;0.9250;0.9300;0.9300;0.9300;0.9300;0.9300;0.9300;0.9300;0.9275;0.9250;0.9300;0.9300;0.9300;0.9500;0.9300 | 1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000 |
| V1     | TDS16                 | GradientBoosting  | 0;0;0;0;0;0;0;0;0;0;0;0;0;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1 | 0.9997;0.9997;0.9997;0.9997;0.9997;0.9997;0.9997;0.9997;0.9997;0.9997;0.9997;0.9997;0.9997;0.9997;0.9997;0.9997;0.9997;0.9997;0.9997;0.9997;0.9997;0.9997;0.9997;0.9997;0.9997;0.9997;0.9997;0.9997;0.9997;0.9997;0.9997;0.9997;0.9997;0.9997;0.9997;0.9997;0.9997;0.9997;0.9997;0.9997;0.9997 | 1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000 |
| V1     | TDS16                 | XGBoost           | 0;0;0;0;0;0;0;0;0;0;0;0;0;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1 | 0.9985;0.9991;0.9991;0.9985;0.9985;0.9985;0.9985;0.9982;0.9985;0.9985;0.9982;0.9990;0.9985;0.9989;0.9991;0.9990;0.9991;0.9991;0.9989;0.9991;0.9991;0.9991;0.9989;0.9989;0.9989;0.9991;0.9990;0.9991;0.9991;0.9991;0.9991;0.9989;0.9991;0.9989;0.9991;0.9989;0.9989;0.9991;0.9989;0.9991;0.9989 | 1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000 |
| V1     | TierA67_clean         | LogReg_l2         | 0;0;0;0;0;0;0;0;0;0;0;0;0;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1 | 0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0001;0.0000;0.0000;0.0000;0.0001;0.0000;0.0002;0.0002;0.0013;0.0005;0.0005;0.0011;0.0009;0.0010;0.0040;0.0020;0.0033;0.0009;0.0002;0.0008;0.0006;0.0021;0.0009;0.0003;0.0007;0.0001;0.0010;0.0001;0.0001;0.0009;0.0003;0.0014;0.0003;0.0024 | 0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000 |
| V1     | TierA67_clean         | LogReg_elasticnet | 0;0;0;0;0;0;0;0;0;0;0;0;0;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1 | 0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000 | 0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000 |
| V1     | TierA67_clean         | RandomForest      | 0;0;0;0;0;0;0;0;0;0;0;0;0;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1 | 0.4075;0.4200;0.4225;0.3925;0.4275;0.4225;0.4000;0.4000;0.4275;0.4075;0.4300;0.4075;0.4300;0.4425;0.4250;0.4400;0.4250;0.4375;0.4475;0.4175;0.4175;0.4375;0.4500;0.4475;0.4525;0.4275;0.4450;0.4350;0.4625;0.4300;0.4375;0.4550;0.4325;0.4400;0.4200;0.4325;0.4450;0.4275;0.4500;0.4625;0.4350 | 0.2364;0.2667;0.2727;0.2000;0.2848;0.2727;0.2182;0.2182;0.2848;0.2364;0.2909;0.2364;0.2909;0.3212;0.2788;0.3152;0.2788;0.3091;0.3333;0.2606;0.2606;0.3091;0.3394;0.3333;0.3455;0.2848;0.3273;0.3030;0.3697;0.2909;0.3091;0.3515;0.2970;0.3152;0.2667;0.2970;0.3273;0.2848;0.3394;0.3697;0.3030 |
| V1     | TierA67_clean         | GradientBoosting  | 0;0;0;0;0;0;0;0;0;0;0;0;0;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1 | 0.2267;0.2267;0.2203;0.0804;0.2203;0.2267;0.1060;0.0904;0.2203;0.2203;0.2847;0.1174;0.2846;0.2847;0.2500;0.2435;0.2500;0.2500;0.0706;0.2847;0.2847;0.2847;0.0742;0.0881;0.0967;0.2570;0.0724;0.2847;0.0718;0.2847;0.2847;0.0894;0.2144;0.2776;0.2847;0.2144;0.2776;0.2847;0.0851;0.0914;0.2755 | 0.2267;0.2267;0.2203;0.0803;0.2203;0.2267;0.1060;0.0904;0.2203;0.2203;0.2846;0.1174;0.2846;0.2846;0.2500;0.2435;0.2500;0.2500;0.0706;0.2846;0.2846;0.2846;0.0742;0.0881;0.0967;0.2570;0.0724;0.2846;0.0718;0.2846;0.2846;0.0894;0.2143;0.2776;0.2846;0.2143;0.2776;0.2846;0.0851;0.0914;0.2755 |
| V1     | TierA67_clean         | XGBoost           | 0;0;0;0;0;0;0;0;0;0;0;0;0;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1 | 0.0557;0.0715;0.0715;0.0557;0.0715;0.0715;0.0624;0.0585;0.0715;0.0485;0.1230;0.0585;0.0879;0.1273;0.1230;0.1230;0.0921;0.1230;0.1071;0.1230;0.0921;0.1230;0.4388;0.2622;0.1793;0.0850;0.4475;0.1230;0.4663;0.0954;0.1230;0.4388;0.0741;0.1273;0.1230;0.0715;0.1273;0.0921;0.1793;0.4113;0.1153 | 0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.3829;0.1447;0.0329;0.0000;0.3946;0.0000;0.4200;0.0000;0.0000;0.3829;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0329;0.3458;0.0000 |
| V1     | TierA67_clean_no_MAPK | LogReg_l2         | 0;0;0;0;0;0;0;0;0;0;0;0;0;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1 | 0.0193;0.0739;0.0366;0.0318;0.0048;0.0590;0.0723;0.1372;0.0933;0.0919;0.1305;0.2659;0.0729;0.4758;0.3694;0.8051;0.6148;0.5327;0.7540;0.6571;0.7051;0.8734;0.8206;0.9382;0.7242;0.2899;0.6860;0.5982;0.8323;0.7365;0.3830;0.6581;0.2098;0.7344;0.2044;0.4713;0.6648;0.3677;0.7593;0.8169;0.9495 | 0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;1.0000;1.0000;0.5484;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;0.0000;1.0000;1.0000;1.0000;1.0000;0.0000;1.0000;0.0000;1.0000;0.0000;0.0000;1.0000;0.0000;1.0000;1.0000;1.0000 |
| V1     | TierA67_clean_no_MAPK | LogReg_elasticnet | 0;0;0;0;0;0;0;0;0;0;0;0;0;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1 | 0.0000;0.0001;0.0001;0.0000;0.0000;0.0000;0.0001;0.0001;0.0001;0.0001;0.0001;0.0003;0.0001;0.0008;0.0006;0.0026;0.0012;0.0010;0.0031;0.0015;0.0015;0.0048;0.0037;0.0110;0.0021;0.0004;0.0023;0.0012;0.0037;0.0022;0.0007;0.0022;0.0006;0.0022;0.0003;0.0009;0.0015;0.0006;0.0032;0.0003;0.0063 | 0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000 |
| V1     | TierA67_clean_no_MAPK | RandomForest      | 0;0;0;0;0;0;0;0;0;0;0;0;0;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1 | 0.4475;0.4575;0.4575;0.4275;0.4575;0.4525;0.4250;0.4300;0.4550;0.4475;0.5000;0.4550;0.4675;0.5100;0.4900;0.4875;0.4850;0.4975;0.4875;0.4625;0.4575;0.4725;0.5300;0.5100;0.5225;0.5050;0.5275;0.4850;0.5275;0.4975;0.4800;0.5350;0.4700;0.4975;0.4625;0.4650;0.4850;0.4625;0.5225;0.5350;0.4900 | 0.3591;0.3812;0.3812;0.3149;0.3812;0.3702;0.3094;0.3204;0.3757;0.3591;0.4751;0.3757;0.4033;0.4972;0.4530;0.4475;0.4420;0.4696;0.4475;0.3923;0.3812;0.4144;0.5414;0.4972;0.5249;0.4862;0.5359;0.4420;0.5359;0.4696;0.4309;0.5525;0.4088;0.4696;0.3923;0.3978;0.4420;0.3923;0.5249;0.5525;0.4530 |
| V1     | TierA67_clean_no_MAPK | GradientBoosting  | 0;0;0;0;0;0;0;0;0;0;0;0;0;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1 | 0.1034;0.0958;0.1034;0.0038;0.0958;0.0958;0.0114;0.0114;0.1034;0.1034;0.3546;0.0134;0.3355;0.3355;0.3355;0.3355;0.3355;0.3355;0.1366;0.3355;0.3546;0.3355;0.7621;0.7619;0.4298;0.3355;0.8990;0.3355;0.7621;0.3355;0.3355;0.7621;0.0958;0.3355;0.3355;0.0958;0.3355;0.3355;0.4298;0.4021;0.0958 | 0.1033;0.0957;0.1033;0.0038;0.0957;0.0957;0.0114;0.0114;0.1033;0.1033;0.3545;0.0134;0.3355;0.3355;0.3355;0.3355;0.3355;0.3355;0.1365;0.3355;0.3545;0.3355;0.7621;0.7620;0.4297;0.3355;0.8990;0.3355;0.7621;0.3355;0.3355;0.7621;0.0957;0.3355;0.3355;0.0957;0.3355;0.3355;0.4297;0.4020;0.0957 |
| V1     | TierA67_clean_no_MAPK | XGBoost           | 0;0;0;0;0;0;0;0;0;0;0;0;0;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1 | 0.0605;0.1046;0.1046;0.0605;0.0996;0.0996;0.0622;0.0622;0.0996;0.0996;0.1022;0.0655;0.1022;0.1075;0.1075;0.1075;0.1075;0.1075;0.1894;0.1075;0.1075;0.1075;0.4869;0.4869;0.1967;0.1075;0.6125;0.1075;0.5445;0.1075;0.1075;0.4869;0.1046;0.1075;0.1075;0.1046;0.1075;0.1075;0.1910;0.5747;0.1464 | 0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0876;0.0000;0.0000;0.0000;0.4612;0.4612;0.0968;0.0000;0.6189;0.0000;0.5335;0.0000;0.0000;0.4612;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0896;0.5715;0.0336 |
| V1     | BRS71_original        | LogReg_l2         | 0;0;0;0;0;0;0;0;0;0;0;0;0;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1 | 0.8483;0.8787;0.8418;0.7110;0.7478;0.7933;0.8004;0.8391;0.6781;0.8212;0.9754;0.9306;0.7436;0.9918;0.9950;0.9912;0.9947;0.9962;0.9972;0.9925;0.9972;0.9955;0.9980;0.9981;0.9978;0.9930;0.9988;0.9963;0.9992;0.9938;0.9975;0.9981;0.9970;0.9991;0.9948;0.9958;0.9984;0.9963;0.9983;0.9995;0.9977 | 1.0000;1.0000;1.0000;0.8959;0.9997;1.0000;1.0000;1.0000;0.8034;1.0000;1.0000;1.0000;0.9879;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000 |
| V1     | BRS71_original        | LogReg_elasticnet | 0;0;0;0;0;0;0;0;0;0;0;0;0;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1 | 0.4336;0.5614;0.3807;0.2967;0.2563;0.3207;0.3869;0.5652;0.4332;0.4319;0.8269;0.5071;0.4980;0.9326;0.9550;0.9685;0.9523;0.9622;0.9677;0.9865;0.9711;0.9797;0.9746;0.9808;0.9685;0.9675;0.9707;0.9599;0.9895;0.9687;0.9687;0.9741;0.9596;0.9893;0.9518;0.9519;0.9869;0.9603;0.9756;0.9658;0.9489 | 0.0000;0.2508;0.0000;0.0000;0.0000;0.0000;0.0000;0.2865;0.0000;0.0000;1.0000;0.0000;0.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000 |
| V1     | BRS71_original        | RandomForest      | 0;0;0;0;0;0;0;0;0;0;0;0;0;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1 | 0.4075;0.4100;0.4050;0.4075;0.4075;0.4150;0.4175;0.4225;0.4125;0.4150;0.3975;0.4150;0.4200;0.4450;0.4500;0.4575;0.4600;0.4550;0.4625;0.4600;0.4350;0.4525;0.4725;0.4600;0.4675;0.4575;0.4650;0.4575;0.4775;0.4675;0.4600;0.4725;0.4875;0.4775;0.4475;0.4350;0.4700;0.4450;0.4725;0.4275;0.4350 | 0.2547;0.2609;0.2484;0.2547;0.2547;0.2733;0.2795;0.2919;0.2671;0.2733;0.2298;0.2733;0.2857;0.3478;0.3602;0.3789;0.3851;0.3727;0.3913;0.3851;0.3230;0.3665;0.4161;0.3851;0.4037;0.3789;0.3975;0.3789;0.4286;0.4037;0.3851;0.4161;0.4534;0.4286;0.3540;0.3230;0.4099;0.3478;0.4161;0.3043;0.3230 |
| V1     | BRS71_original        | GradientBoosting  | 0;0;0;0;0;0;0;0;0;0;0;0;0;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1 | 0.0017;0.0025;0.0022;0.0017;0.0058;0.0031;0.0057;0.0035;0.0035;0.0058;0.0022;0.0034;0.0124;0.0094;0.0296;0.0250;0.0250;0.0250;0.0296;0.0133;0.0079;0.0250;0.0244;0.0244;0.0304;0.0250;0.0304;0.0296;0.0393;0.0244;0.0250;0.0296;0.0296;0.0250;0.0077;0.0053;0.0332;0.0304;0.0112;0.0013;0.0049 | 0.0017;0.0025;0.0021;0.0017;0.0058;0.0031;0.0057;0.0035;0.0035;0.0058;0.0021;0.0034;0.0124;0.0094;0.0296;0.0250;0.0250;0.0250;0.0296;0.0133;0.0079;0.0250;0.0244;0.0244;0.0304;0.0250;0.0304;0.0296;0.0392;0.0244;0.0250;0.0296;0.0296;0.0250;0.0077;0.0053;0.0332;0.0304;0.0112;0.0013;0.0049 |
| V1     | BRS71_original        | XGBoost           | 0;0;0;0;0;0;0;0;0;0;0;0;0;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1 | 0.6727;0.6727;0.6727;0.6727;0.6727;0.6727;0.6727;0.5179;0.5179;0.6727;0.6727;0.5179;0.8311;0.5179;0.7296;0.7342;0.7342;0.7342;0.7342;0.7342;0.5179;0.7342;0.7342;0.7342;0.7342;0.7342;0.7342;0.7296;0.7342;0.7342;0.7342;0.7342;0.7342;0.7342;0.5179;0.7296;0.7342;0.7342;0.7342;0.7342;0.7342 | 0.7315;0.7315;0.7315;0.7315;0.7315;0.7315;0.7315;0.5013;0.5013;0.7315;0.7315;0.5013;0.9669;0.5013;0.8160;0.8229;0.8229;0.8229;0.8229;0.8229;0.5013;0.8229;0.8229;0.8229;0.8229;0.8229;0.8229;0.8160;0.8229;0.8229;0.8229;0.8229;0.8229;0.8229;0.5013;0.8160;0.8229;0.8229;0.8229;0.8229;0.8229 |
| V1     | v3_panel_compact      | LogReg_l2         | 0;0;0;0;0;0;0;0;0;0;0;0;0;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1 | 0.9885;0.9854;0.9942;0.9896;0.9906;0.9866;0.9838;0.9936;0.9897;0.9906;0.9964;0.9956;0.9870;0.9980;0.9970;0.9989;0.9983;0.9982;0.9987;0.9995;0.9991;0.9992;0.9987;0.9993;0.9987;0.9966;0.9972;0.9978;0.9982;0.9985;0.9984;0.9989;0.9962;0.9990;0.9977;0.9973;0.9989;0.9965;0.9992;0.9999;0.9993 | 1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000 |
| V1     | v3_panel_compact      | LogReg_elasticnet | 0;0;0;0;0;0;0;0;0;0;0;0;0;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1 | 0.9937;0.9924;0.9977;0.9951;0.9953;0.9932;0.9920;0.9967;0.9956;0.9952;0.9981;0.9976;0.9945;0.9989;0.9987;0.9994;0.9991;0.9990;0.9993;0.9997;0.9995;0.9996;0.9993;0.9996;0.9994;0.9982;0.9986;0.9989;0.9991;0.9992;0.9992;0.9994;0.9982;0.9994;0.9989;0.9986;0.9994;0.9982;0.9995;0.9999;0.9995 | 1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000 |
| V1     | v3_panel_compact      | RandomForest      | 0;0;0;0;0;0;0;0;0;0;0;0;0;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1 | 0.5350;0.5450;0.5425;0.5350;0.5475;0.5500;0.5550;0.5575;0.5600;0.5300;0.5700;0.5850;0.6000;0.7000;0.6650;0.6575;0.6375;0.6850;0.6775;0.6900;0.6500;0.6700;0.6825;0.6825;0.7125;0.6625;0.6650;0.6650;0.6775;0.6700;0.6550;0.6875;0.6575;0.6800;0.6600;0.6425;0.6700;0.6625;0.6875;0.6075;0.6025 | 0.6875;0.7125;0.7062;0.6875;0.7188;0.7250;0.7375;0.7437;0.7500;0.6750;0.7750;0.8125;0.8500;1.0000;1.0000;0.9937;0.9437;1.0000;1.0000;1.0000;0.9750;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;1.0000;0.9875;1.0000;0.9937;1.0000;1.0000;0.9562;1.0000;1.0000;1.0000;0.8688;0.8563 |
| V1     | v3_panel_compact      | GradientBoosting  | 0;0;0;0;0;0;0;0;0;0;0;0;0;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1 | 0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000 | 0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000;0.0000 |
| V1     | v3_panel_compact      | XGBoost           | 0;0;0;0;0;0;0;0;0;0;0;0;0;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1 | 0.2023;0.2023;0.2023;0.2023;0.2023;0.2023;0.2023;0.2023;0.2023;0.2023;0.2023;0.2023;0.2023;0.2023;0.2023;0.2023;0.2023;0.2023;0.2023;0.2023;0.2023;0.2023;0.2023;0.2023;0.2023;0.2023;0.2023;0.2023;0.2023;0.2023;0.2023;0.2023;0.2023;0.2023;0.2023;0.2023;0.2023;0.2023;0.2023;0.2023;0.2023 | 0.1040;0.1040;0.1040;0.1040;0.1040;0.1040;0.1040;0.1040;0.1040;0.1040;0.1040;0.1040;0.1040;0.1040;0.1040;0.1040;0.1040;0.1040;0.1040;0.1040;0.1040;0.1040;0.1040;0.1040;0.1040;0.1040;0.1040;0.1040;0.1040;0.1040;0.1040;0.1040;0.1040;0.1040;0.1040;0.1040;0.1040;0.1040;0.1040;0.1040;0.1040 |

### ✅ `results/tables/v3_permutation_null.tsv`
_Rows: 1_

|   real_auc |   null_mean |   null_p5 |   null_p95 |   empirical_p |   n_perm |
|-----------:|------------:|----------:|-----------:|--------------:|---------:|
|      0.923 |       0.723 |     0.689 |      0.757 |      0.000999 |     1000 |


## 5. KNOWN LIMITATIONS + NEXT STEPS


### analysis_summary_v2.md · Remaining limitations

## Remaining limitations

- BRS71 is a proxy derived from TCGA-THCA, not the original published list.
- TierA67 hardcoded composite contains 67 category entries but 66 unique HGNC-normalized symbols because PAX8 appears in two categories.
- GSE27155 external labels are phenotype proxies inferred from histology/morphology, not direct molecular subtypes.
- GSE126698 external subset is very small after restricting to PTC vs FTC.
- TCGA mutation anchors depend on locally available open MAFs, not a controlled-access harmonized mutation callset.

### next_steps_v2.md (verbatim)

```markdown
# next_steps_v2

1. Replace the BRS71 proxy with the original source-verified supplementary list once the exact official table is pinned down.
2. Add a better external RNA-seq validation cohort with explicit follicular-patterned tumor labels rather than generic PTC-only labels.
3. Reprocess Affymetrix CEL files in R/Bioconductor for a cleaner microarray validation layer.
4. Expand TCGA anchors beyond BRAF V600E and RAS to fusion-aware subtype definitions when a consistent public callset is available.
5. Prepare a small Phase 1 pilot model package using only TDS16 and TierA67 on clearly labeled differentiated thyroid tumors.
```

## 6. BUILD STATE


### reports/html/_version.json

```json
{
  "build_time": "2026-04-24 11:29:10",
  "git_hash": "no-git",
  "n_samples": 1509,
  "n_datasets": 7,
  "n_genes": 51711,
  "v3_build_time": "2026-04-24 06:20:25 UTC",
  "v3_pages": [
    "19_honesty_audit.html",
    "20_panel_size.html",
    "21_three_class.html",
    "22_survival.html",
    "23_bethesda_sim.html",
    "24_multimodal.html"
  ],
  "v3_ext_pages": [
    25,
    26,
    27,
    28
  ],
  "v3_ext_build_time": "2026-04-24T06:36:07.810236Z",
  "v3_ext_primary_auc": 1.0,
  "v3_complete": true,
  "v3_ext_complete": true,
  "thyrai_rebrand": true,
  "last_update": "2026-04-24 16:22:51",
  "download_zip": "/reports/html/thyrai_repo_latest.zip",
  "total_pages": 32
}
```

- **pages/*.html count:** 46
- **figs_interactive/*.html count:** 105
- **Port 8012 /reports/html/index.html:** 200
