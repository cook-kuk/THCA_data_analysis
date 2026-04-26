# UI Static Figure QA

Date: 2026-04-23

Scope: quick served-source validation against `http://127.0.0.1:8012/reports/html/index.html` and representative report pages. No dashboard HTML/CSS/JS edited.

## Result

- `index.html` served successfully (`200 OK`) and does not itself contain static figure-wall iframes.
- Static iframe figure walls do appear in served page source for representative content pages:
  - `pages/01_overview.html` (`overview_sankey`, `overview_timeline`, `overview_consort`)
  - `pages/03_sample_master.html` (`datasets_sample_bar`, `datasets_quality_grade`, `caveats_quality_grade`)
  - `pages/05_eda.html` (`cohort_compare_pca`, `cohort_batch_severity`, `sample_histology_dataset_heatmap`)
  - `pages/07_ml_baseline.html` (`panel_perf_curve`, `panel_cost_utility`, `panel_marginal_gain`, `shap_bar`, `shap_swarm`)
  - `pages/08_panel_comparison.html` (`scores_tds_violin`, `scores_brs_violin`, `ml_roc_overlay`)
  - `pages/10_gene_explorer.html` (`sample_histology_dataset_heatmap`, `scores_density`, `cohort_label_concordance`)
- Sample iframe targets also returned `200 OK` at the served `figs_interactive/*.html` URLs.

## Duplicate Check

- No obvious duplicate iframe embeds were found within the same page across the pages that currently contain iframe walls.
- One obvious cross-page duplicate remains: `../figs_interactive/sample_histology_dataset_heatmap.html` is embedded in both `pages/05_eda.html` and `pages/10_gene_explorer.html`.

## Verdict

Static iframe figure walls are present in served source for representative report pages. Obvious same-page duplicates were not found, but at least one cross-page duplicate figure remains.
