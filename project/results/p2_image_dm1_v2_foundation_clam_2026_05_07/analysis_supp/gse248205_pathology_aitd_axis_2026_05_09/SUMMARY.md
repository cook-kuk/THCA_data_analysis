# GSE248205 H&E-to-AITD spatial-axis scout

Fast handcrafted H&E feature screen across 16,985 spots and 8 thyroid Visium samples (Control n=2, GD n=3, HT n=3).

Verdict: `NO_GO_LOCAL_HE_CAVEAT_ONLY`. AP/TLS H&E prediction looks strong when pooled (rho 0.590), but it fails the sample-centered test (rho -0.090) and domain aggregation (sample-centered domain rho -0.164). HLA-II/AP, B/TLS, and thyrocyte axes are also not convincing after sample centering.

Use only as a negative control/caveat: not every thyroid autoimmune ST cohort supports a local H&E-to-spatial-RNA claim.

Primary artifacts:
- `GSE248205_PATHOLOGY_AITD_REPORT.md`
- `GSE248205_PATHOLOGY_AITD_SUMMARY.json`
- `gse248205_pathology_model_summary.tsv`
- `gse248205_pathology_domain_summary.tsv`
- `fig_gse248205_pathology_aitd_axis.png`
