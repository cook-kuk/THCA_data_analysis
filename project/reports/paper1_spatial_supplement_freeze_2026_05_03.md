# Paper 1 spatial supplement — FREEZE
**Date: 2026-05-03 · Status: frozen for manuscript marathon · Owner: Seungho Cook**

## Final verdict

**Verdict C — multi-cohort supplementary upgrade.** Three Visium ST cohorts (GSE250521 cancer progression PT→ATC; GSE230424 PTC + Hashimoto-overlap; GSE248205 autoimmune thyroid without cancer) provide supportive validation of the 8-gene RAI / DM1-like axis in Paper 1.

This is **not a main-figure-grade** result. The cancer progression direct evidence remains underpowered in GSE250521 alone (sample-mean Spearman ρ = −0.35, p = 0.18 in epithelial-enriched depth-corrected subset). External cohorts strengthen interpretability of the axis itself, not the stage-progression claim.

## Why not main-figure-grade
- GSE250521 (the only cancer-progression cohort) is a single 16-slide dataset with depth-corrected sample-mean Spearman p = 0.18 — direction-consistent but underpowered.
- GSE230424 has only one condition (PTC + Hashimoto-overlap) with no within-cohort progression contrast.
- GSE248205 has no cancer samples, so cannot speak to PT→ATC progression.
- Therefore primary cancer progression evidence remains a single underpowered cohort. The supplementary package supports the molecular validity of the score, not the stage trend.

## Risks addressed by the supplementary package

### (a) TDS overlap / score tautology concern
The TDS_overlap score shares 6 of 8 genes with RAI_8 (TG, TPO, TSHR, SLC5A5, DIO1, PAX8), so the strong DM1_like ↔ TDS_overlap anti-correlation in GSE250521 (ρ = −0.89) was partially tautological. We addressed this by **non-overlapping lineage cross-validation** using THYROID_NONOVERLAP (SLC26A4, IYD, DUOX1, DUOX2, TFF3, HHEX, GLIS3, DIO2; **zero gene overlap with RAI_8**). Across the 12 external slides, sample-mean depth-residualized epithelial top-25% Spearman ρ = −0.85 to −1.00 (p < 0.01 in both datasets, Pearson r = −0.99 to −0.996). This **supports** that the DM1-axis indexes a thyroid-lineage signal that survives outside the original 8-gene panel.

### (b) Inflammation artifact concern
GSE248205 (2 control, 3 Hashimoto's, 3 Graves' samples; no cancer) was used as a non-cancer autoimmune-inflammation control. Sample-mean DM1_like (depth-residualized, epithelial top 25%) was −0.029 in CONTROL, −0.370 in HT, −0.065 in GD; Mann-Whitney p > 0.10 for all pairwise comparisons (n = 2-3 per condition). This **argues against** a generic-inflammation interpretation of DM1_like — Hashimoto's and Graves' inflammation, on their own, do not raise DM1 above control-tissue baseline.

## Remaining limitation
Cancer progression direct evidence remains direction-consistent but underpowered: sample-mean Spearman against ordinal stage in GSE250521 epithelial-enriched depth-corrected subset is ρ = −0.35 (p = 0.18) for RAI_8 and ρ = +0.35 (p = 0.18) for DM1_like. The two external cohorts contribute non-overlapping lineage validation and an autoimmune negative control, but neither contains a PT→PTC→LPTC→ATC progression contrast. Independent confirmation of the stage trend will require additional spatial cancer cohorts in future work; this is deferred outside Paper 1 scope.

---

## Draft text (manuscript-ready, allowed-vocabulary)

### Results paragraph draft

> External Visium spatial transcriptomics cohorts provide supportive validation of the 8-gene RAI / DM1-like axis. In a 16-slide thyroid cancer cohort (GSE250521; 4 PT, 4 PTC, 4 LPTC, 4 ATC) the DM1-like score showed a direction-consistent but underpowered increase from PT to ATC in the epithelial top-50% depth-residualized subset (sample-mean Spearman ρ = +0.35, p = 0.18; n = 16). To address the partial gene overlap between the RAI_8 panel and conventional thyroid-differentiation scores, we performed non-overlapping lineage cross-validation using a zero-overlap thyroid-lineage gene set (SLC26A4, IYD, DUOX1, DUOX2, TFF3, HHEX, GLIS3, DIO2) in two additional Visium cohorts (GSE230424, n = 4 PTC + Hashimoto-overlap; GSE248205, n = 8 control / Hashimoto's / Graves'). Sample-mean DM1-like and the non-overlapping thyroid-lineage score anti-correlated strongly after depth correction, in both cohorts and at multiple epithelial-enrichment thresholds (epithelial top-50% sample-mean Spearman ρ = −1.00 in GSE230424, Pearson r = −0.99, p = 8 × 10⁻³, n = 4; epithelial top-25% sample-mean ρ = −0.88 in GSE248205, Pearson r = −0.996, p = 2 × 10⁻⁷, n = 8). In GSE248205 (no cancer), the DM1-like score was not raised above CONTROL baseline by Hashimoto's thyroiditis or Graves' disease alone (sample-mean DM1-like = −0.029 / −0.370 / −0.065 for CONTROL / HT / GD; Mann-Whitney p > 0.10 for all pairwise comparisons; n = 2 / 3 / 3).

### Methods sentence draft

> 8-gene RAI scoring and the depth-residualized variant were performed identically across all three external ST cohorts: per-spot log-normalized expression of each gene was residualized within sample against log(total_counts) + log(n_genes_by_counts), z-scored across spots within sample, and averaged. DM1_like = −RAI_8. Non-overlapping lineage cross-validation used SLC26A4, IYD, DUOX1, DUOX2, TFF3, HHEX, GLIS3, DIO2, scored identically. Stage and condition trends were tested at the sample-mean level by Spearman correlation (cancer progression in GSE250521) or Mann-Whitney U (autoimmune negative control in GSE248205); spot-level statistics were reported for transparency only.

### Discussion sentence draft

> Cross-validation against a non-overlapping thyroid-lineage gene set (SLC26A4, IYD, DUOX1, DUOX2, TFF3, HHEX, GLIS3, DIO2; zero gene overlap with RAI_8) supports the 8-gene RAI signal as a genuine thyroid-lineage / dedifferentiation read-out rather than a panel-internal artifact (sample-mean Spearman ρ < −0.85 in two independent Visium cohorts, depth-corrected, p < 0.01). The DM1-like score was not elevated by autoimmune thyroid inflammation alone (Hashimoto's, Graves') in a non-cancer cohort, which argues against a generic immune-microenvironment interpretation of the score; we note, however, that with two control / three Hashimoto's / three Graves' samples per condition, this comparison is direction-consistent but underpowered.

### Limitation sentence draft

> The cancer progression direct evidence (PT → ATC) remains underpowered to a single 16-slide cohort (GSE250521; sample-mean depth-corrected Spearman ρ = −0.35, p = 0.18 in epithelial-enriched subset). The two external cohorts (GSE230424, GSE248205) used in this study contribute non-overlapping lineage cross-validation and an autoimmune negative control, but neither contains a PT → PTC → LPTC → ATC progression contrast; independent confirmation of the stage trend in additional cancer ST cohorts is left to future work.

### Figure caption drafts

**Supp Fig X1.** GSE250521 spatial maps, PT → PTC → LPTC → ATC. 2 rows × 4 columns of representative slides (one per stage). Top row: RAI_8 score (within-sample z, log-normalized expression of TPO, DIO1, TSHR, PAX8, TG, FOXE1, NKX2-1, SLC5A5). Bottom row: DM1_like score (= −RAI_8). Overlays use scaled hires H&E and within-sample z-score color scaling. Single-cohort sample-mean Spearman against ordinal stage (PT = 0, ATC = 3) in epithelial top-50% subset: depth-corrected ρ = −0.35 for RAI_8 (p = 0.18, n = 16). Direction-consistent but underpowered at this n.

**Supp Fig X2.** Non-overlapping lineage cross-validation. Sample-mean DM1_like score versus sample-mean THYROID_NONOVERLAP score (SLC26A4, IYD, DUOX1, DUOX2, TFF3, HHEX, GLIS3, DIO2; zero gene overlap with RAI_8), both depth-residualized and computed on the epithelial top-25% subset within sample. n = 12 external slides (4 GSE230424 PTC + Hashimoto-overlap; 8 GSE248205 control / Hashimoto's / Graves'). Spearman ρ and Pearson r reported across all 12 slides. Per-cohort statistics: GSE230424 epithelial top-50% sample-mean Spearman ρ = −1.00 (Pearson r = −0.99, p = 8 × 10⁻³, n = 4); GSE248205 epithelial top-25% sample-mean Spearman ρ = −0.88 (Pearson r = −0.996, p = 2 × 10⁻⁷, n = 8). Markers: circles = GSE230424, squares = GSE248205. Color: condition (CONTROL / HT / GD / PTC + Hashimoto-overlap).

**Supp Fig X3.** Autoimmune thyroid (no cancer) control. GSE248205 sample-mean DM1_like score, depth-residualized, epithelial top-25% subset, by condition (CONTROL n = 2, HT n = 3, GD n = 3). Bar = mean ± SEM. Sample-mean DM1_like = −0.029 / −0.370 / −0.065 for CONTROL / HT / GD; Mann-Whitney p > 0.10 for all pairwise comparisons. Direction-consistent but underpowered; argues against a generic-inflammation interpretation of the DM1-axis without claiming proof.

### Supp table caption drafts

**Supp Table SX.** Sample-level scoring summary across 28 Visium slides from three cohorts. Rows: 16 GSE250521 (PT / PTC / LPTC / ATC), 4 GSE230424 (PTC + Hashimoto-overlap), 8 GSE248205 (CONTROL / HT / GD). Columns: dataset, condition, n_spots, sample-mean RAI_8 / DM1_like / TDS_overlap / THYROID_NONOVERLAP / Epithelial / Proliferation scores at multiple aggregation scopes (all spots, epithelial top 50%, epithelial top 25%) in raw within-sample-z and (where available) depth-residualized form. THYROID_NONOVERLAP scoring was performed only on the 12 external slides; the GSE250521 entry has the relevant cells empty by design and refers to the original cohort scoring run for raw scores plus the separate depth-residualization analysis.

**Supp Table SX+1.** Independent thyroid-lineage cross-validation and technical confounding. Combines (i) Spearman / Pearson correlations between sample-mean DM1_like and THYROID_NONOVERLAP across all-spots and epithelial-enriched subsets in GSE230424 + GSE248205; and (ii) per-cohort correlations of each score with log(total_counts) flagging |ρ| > 0.30. Raw-version scores carry substantial depth confounding in both external cohorts (|ρ| = 0.30 to 0.54); per-gene depth-residualization reduces this to |ρ| ≤ 0.05 across all scores including the THYROID_NONOVERLAP control.

---

## Exact file inventory (frozen)

```
project/supplementary/spatial_freeze_2026_05_03/
├── SuppFig_X1_GSE250521_spatial_PT_to_ATC.png
├── SuppFig_X1_GSE250521_spatial_PT_to_ATC.pdf
├── SuppFig_X2_DM1_thyroid_nonoverlap_crossvalidation.png
├── SuppFig_X2_DM1_thyroid_nonoverlap_crossvalidation.pdf
├── SuppFig_X3_GSE248205_autoimmune_negative_control.png
├── SuppFig_X3_GSE248205_autoimmune_negative_control.pdf
├── SuppTable_SX_28sample_score_summary.tsv          (28 rows × 18 cols)
├── SuppTable_SX1_independent_lineage_and_technical_confounding.tsv  (28 rows × 12 cols)
├── _make_supp_figs_x2_x3.py                         (regen-from-frozen-data script)
└── _make_supp_tables.py                             (regen-from-frozen-data script)
```

Source-of-truth pointers (read-only after this freeze; any change requires a new dated freeze):
- `project/results/01_spatial_score/all_spots_scored.tsv.gz` (GSE250521, 55,873 spots × 22 cols)
- `project/results/02_stage_trend/stage_trend_summary.csv`
- `project/results/02_stage_trend/depth_corrected_trend.csv`
- `project_external_st/results/scores/all_external_spots_scored.tsv.gz` (29,683 spots × 41 cols)
- `project_external_st/results/meta/sample_level_score_summary.tsv` (12 external samples × 80 cols)
- `project_external_st/results/meta/independent_lineage_validation.tsv`
- `project_external_st/results/meta/technical_confounding_summary.tsv`
- `project_external_st/results/meta/condition_comparison.tsv`
- `project_external_st/results/registry/dataset_sample_registry.tsv`

Background reports (frozen context, not in this supplement bundle):
- `project/reports/gse250521_dm1_spatial_validation_report.md`
- `project/reports/gse250521_paper1_inclusion_verdict_2026_05_03.md`
- `project_external_st/reports/external_st_triage_verdict_2026_05_03.md`

---

Spatial supplement frozen. Return to manuscript marathon.
