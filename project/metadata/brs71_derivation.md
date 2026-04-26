# BRS71 proxy derivation

- Status: proxy
- Source cohort: TCGA-THCA
- Source expression: `project/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_log2.tsv`
- Mutation anchor source: open GDC masked somatic MAF files available locally
- Positive class: BRAF V600E only
- Negative class: NRAS/HRAS/KRAS mutant
- Exclusions: samples with both anchors, neither anchor, or unavailable MAF
- Differential expression implementation: Welch t-test on log2 normalized expression
- Threshold: abs(log2FC) > 1 and FDR < 0.05
- Selected genes: top 71 by significance among threshold-passing genes
- Important note: this is a BRS71 proxy derived internally from TCGA-THCA, not original Chakravarty 2011
- Mutation-anchor counts used: {'BRAF_like': 279, 'RAS_like': 54}
