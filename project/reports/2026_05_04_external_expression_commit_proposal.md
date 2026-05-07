# External expression validation — commit proposal (NOT executed)

Per spec, no commit was performed. Below are the proposed groups for the user to execute.

## Working tree footprint

- `M  .gitignore` — adds `project/results/p_external_expression_validation/raw/` to ignore list
- `??` `project/reports/2026_05_04_external_expression_access_check.md`
- `??` `project/reports/2026_05_04_external_expression_validation_report.md`
- `??` `project/reports/2026_05_04_external_expression_commit_proposal.md` (this file)
- `??` `project/results/p_external_expression_validation/` (whole subtree)
  - `raw/` (62 MB) — **gitignored**, will not appear in any commit
  - tracked artifacts: ~30 MB (4 gene-level matrices) + ~70 KB (tables) + ~1.1 MB (4 PNGs) + scripts

## Proposed commit groups

### Commit 1 — access check, report, scripts, gitignore
Lightweight, narrative + reproducibility.

**Add:**
```
.gitignore
project/reports/2026_05_04_external_expression_access_check.md
project/reports/2026_05_04_external_expression_validation_report.md
project/reports/2026_05_04_external_expression_commit_proposal.md
project/results/p_external_expression_validation/scripts/run_external_validation.py
project/results/p_external_expression_validation/scripts/make_figures.py
```

Suggested message:
> docs+infra: external expression validation sweep — access check, report, scripts (4 GPL570 cohorts; GSE126698 dropped at access)

### Commit 2 — processed matrices, scores, summary tables
Per-dataset gene-level matrices + score outputs. Total ~30 MB of `.tsv.gz`.

**Add:**
```
project/results/p_external_expression_validation/sample_metadata.tsv
project/results/p_external_expression_validation/external_gene_coverage.tsv
project/results/p_external_expression_validation/external_score_tests.tsv
project/results/p_external_expression_validation/external_direction_consistency.tsv
project/results/p_external_expression_validation/external_spearman.tsv
project/results/p_external_expression_validation/external_sample_scores.tsv.gz
project/results/p_external_expression_validation/GSE33630_expression_gene_log.tsv.gz
project/results/p_external_expression_validation/GSE29265_expression_gene_log.tsv.gz
project/results/p_external_expression_validation/GSE65144_expression_gene_log.tsv.gz
project/results/p_external_expression_validation/GSE53157_expression_gene_log.tsv.gz
```

Suggested message:
> data: external expression validation — gene-level expression matrices + score tests (4 GPL570 cohorts, n=205)

**Note:** ~30 MB of compressed TSV is tolerable for git but not ideal. Alternative: move the 4 `*_expression_gene_log.tsv.gz` to a gitignored cache (regenerable from `raw/` + `scripts/run_external_validation.py`) and keep only the score tables (~70 KB). User decision.

### Commit 3 — figures
Four PNGs, ~1.1 MB total.

**Add:**
```
project/results/p_external_expression_validation/external_rai_lineage_boxplots.png
project/results/p_external_expression_validation/external_dm1_nonoverlap_scatter_grid.png
project/results/p_external_expression_validation/external_direction_consistency_forest.png
project/results/p_external_expression_validation/external_dataset_qc_heatmap.png
```

Suggested message:
> figs: external expression validation — boxplots / scatter grid / direction-consistency forest / coverage QC

### NOT included in any commit
- `project/results/p_external_expression_validation/raw/` (62 MB) — Series Matrix `*.gz`, GPL570 annotation. Reproducible from public GEO URLs; gitignored.

## Forbidden-action audit (sweep total)

- ✅ No GEO accessions outside the user-supplied 8.
- ✅ No raw FASTQ alignment.
- ✅ No raw CEL bulk processing.
- ✅ No RunPod / GPU.
- ✅ No TCGA WSI touch.
- ✅ No H&E-DM1 retry.
- ✅ No TCGA methylation touch.
- ✅ No DepMap / CCLE / PRISM touch.
- ✅ No Paper 3 / Paper 4 work.
- ✅ No manuscript prose written.
- ✅ No voice-protected section drafted (Hook / Aim / Discussion / Limitations / Cover / Q9).
- ✅ No commit performed.
