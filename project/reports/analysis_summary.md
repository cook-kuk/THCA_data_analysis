# analysis_summary

## Confirmed

- Usable datasets in this run: GSE213647, GSE126698, GSE27155, GSE76039, GSE97466, TCGA-THCA
- Direct histology annotations were available for: {"GSE213647": 625, "TCGA-THCA": 541, "GSE27155": 66, "GSE97466": 50, "GSE76039": 37, "GSE126698": 28}
- Partly inferred / weaker labels were present for: {"GSE97466": 91, "GSE27155": 33, "TCGA-THCA": 31, "GSE213647": 7}
- Bulk RNA-seq and microarray were processed separately. No pooled training across modalities was performed.
- MTC samples were preserved as separate labels and not merged with follicular-derived tumors.
- TCGA-THCA backbone model ran: yes.
- External validation on GSE213647 did not constitute a valid 2-class test because the follicular-neoplasm samples were annotated as normal rather than tumor in GEO metadata.

## Still Uncertain

- Exact source-verified BRS71 list was not recovered in this run; coverage and model use for BRS71 were skipped.
- Exact definition of the requested TierA67 panel was not supplied; TierA67 coverage and model use were skipped.
- GSE213647 contains a large patient-level bulk RNA-seq cohort but label alignment to TCGA is weaker than TCGA internal histology labels.
- The attempted TCGA-to-GSE213647 external validation was not reliable because only one effective tumor class remained after metadata filtering.
- Mutation calls were used only as auxiliary driver anchors. They were not used to finalize phenotype labels.

## Modality Integration

- Bulk RNA-seq and microarray should not be merged directly for supervised modeling at this stage.
- Platform effects and label non-equivalence remain large enough that microarray is treated as exploratory / directional validation only.
- Bulk RNA-seq external validation should be treated as provisional until a tumor-only follicular-patterned cohort with direct labels is secured.

## Coverage

```tsv
dataset	total_gene_count	TDS16_coverage_n	TDS16_coverage_frac	TDS16_covered_genes	BRS71_coverage_n	BRS71_coverage_frac	BRS71_covered_genes	TierA67_coverage_n	TierA67_coverage_frac	TierA67_covered_genes
TCGA-THCA	51711	16	1.0	DIO1,DIO2,DUOX1,DUOX2,FOXE1,GLIS3,NKX2-1,PAX8,SLC26A4,SLC5A5,SLC5A8,TG,THRA,THRB,TPO,TSHR	0			0		
GSE213647	51389	0	0.0		0			0		
GSE126698	57773	14	0.875	DIO1,DIO2,DUOX1,DUOX2,FOXE1,GLIS3,PAX8,SLC26A4,SLC5A5,TG,THRA,THRB,TPO,TSHR	0			0		
GSE27155	13237	14	0.875	DIO1,DIO2,DUOX1,DUOX2,FOXE1,NKX2-1,PAX8,SLC26A4,SLC5A5,TG,THRA,THRB,TPO,TSHR	0			0		
GSE76039	22880	16	1.0	DIO1,DIO2,DUOX1,DUOX2,FOXE1,GLIS3,NKX2-1,PAX8,SLC26A4,SLC5A5,SLC5A8,TG,THRA,THRB,TPO,TSHR	0			0
```

## sklearn baseline

```tsv
task	feature_set	model	n_train_samples	n_external_samples	cv_auc	cv_balanced_accuracy	cv_f1	external_auc	external_balanced_accuracy	external_f1	status
PTC_like_vs_follicular_patterned	TDS16	LogisticRegression	539.0	349.0	0.1824407936846596	0.2504053765735011	0.3201320132013201		0.0	0.0	ok
PTC_like_vs_follicular_patterned	TDS16	ElasticNetLogistic	539.0	349.0	0.1797098357158096	0.2469810113078728	0.3223684210526316		0.0	0.0	ok
PTC_like_vs_follicular_patterned	TDS16	LinearSVC	539.0	349.0	0.1829315126946874	0.2447194367399189	0.3273322422258592		0.0	0.0	ok
PTC_like_vs_follicular_patterned	TDS16	RandomForest	539.0	349.0	0.1602410923831875	0.3603797738425432	0.0745341614906832		0.0	0.0	ok
PTC_like_vs_follicular_patterned	TDS16	GradientBoosting	539.0	349.0	0.1844676765521655	0.3272989118839343	0.115079365079365		0.0	0.0	ok
PTC_like_vs_follicular_patterned	BRS71	LogisticRegression									skipped_insufficient_features
PTC_like_vs_follicular_patterned	BRS71	ElasticNetLogistic									skipped_insufficient_features
PTC_like_vs_follicular_patterned	BRS71	LinearSVC									skipped_insufficient_features
PTC_like_vs_follicular_patterned	BRS71	RandomForest									skipped_insufficient_features
PTC_like_vs_follicular_patterned	BRS71	GradientBoosting									skipped_insufficient_features
PTC_like_vs_follicular_patterned	TierA67	LogisticRegression									skipped_insufficient_features
PTC_like_vs_follicular_patterned	TierA67	ElasticNetLogistic									skipped_insufficient_features
PTC_like_vs_follicular_patterned	TierA67	LinearSVC									skipped_insufficient_features
PTC_like_vs_follicular_patterned	TierA67	RandomForest									skipped_insufficient_features
PTC_like_vs_follicular_patterned	TierA67	GradientBoosting									skipped_insufficient_features
PTC_like_vs_follicular_patterned	Top50Var	LogisticRegression	539.0	349.0	0.1615532323447834	0.2318647322381054	0.2857142857142857		1.0	1.0	ok
PTC_like_vs_follicular_patterned	Top50Var	ElasticNetLogistic	539.0	349.0	0.1573074461275869	0.220364838916151	0.2789915966386554		1.0	1.0	ok
PTC_like_vs_follicular_patterned	Top50Var	LinearSVC	539.0	349.0	0.1761681245999573	0.228440366972477	0.288107202680067		1.0	1.0	ok
PTC_like_vs_follicular_patterned	Top50Var	RandomForest	539.0	349.0	0.1229891188393428	0.2847877106891402	0.0946745562130177		0.0	0.0	ok
PTC_like_vs_follicular_patterned	Top50Var	GradientBoosting	539.0	349.0	0.1407403456368679	0.2712182632814167	0.1371428571428571		1.0	1.0	ok
```

## Most Realistic Next Step

- Source-verify the BRS71 and TierA67 panel definitions from the exact intended publications or internal spec.
- Add stricter TCGA subtype curation including tall-cell vs classical vs FVPTC splits and fusion-aware anchors where available.
- Replace the current GSE213647 external test with a cohort whose follicular-patterned samples are tumor-labeled and directly comparable to TCGA.
- If raw CEL reprocessing is required, install an R/Bioconductor stack and re-run microarray from CEL with RMA instead of relying on GEO processed matrices.
