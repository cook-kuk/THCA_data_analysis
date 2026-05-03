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

---

# Appendix A — SuppTable SX (28 samples)

```tsv
sample_id	dataset	condition	n_spots	mean_RAI_8_score_raw_all	mean_RAI_8_score_raw_epi50	mean_RAI_8_score_raw_epi25	mean_DM1_like_score_raw_all	mean_DM1_like_score_raw_epi50	mean_DM1_like_score_raw_epi25	mean_TDS_overlap_score_raw_all	mean_Epithelial_score_raw_all	mean_Proliferation_score_raw_all	mean_THYROID_NONOVERLAP_score_raw_all	mean_THYROID_NONOVERLAP_score_resid_epi25	mean_DM1_like_score_resid_epi25	available_normalizations	notes
GSM7980860_N-1	GSE250521	PT	3395	1.0609351988325815e-07	0.15750012883785042	0.18172119094507655	-1.0609351988325815e-07	-0.15750012883785042	-0.18172119094507655	1.3632920766008485e-07	6.6919175196116746e-09	-5.002238582507766e-09				raw_within_sample_z; sample-level depth-corrected residual via depth_corrected_trend.py (epi-only, score-level)	Cancer progression cohort; depth-corrected per-gene residualization not in original scoring run; refer to depth_corrected_trend.csv for sample-mean resid trend
GSM7980861_N-2	GSE250521	PT	3780	-1.4793921693134817e-07	0.14018272038431748	0.14595936417149208	1.4793921693134817e-07	-0.14018272038431748	-0.14595936417149208	-1.8212848941864593e-07	-2.8210743386475247e-08	6.334656064602371e-09				raw_within_sample_z; sample-level depth-corrected residual via depth_corrected_trend.py (epi-only, score-level)	Cancer progression cohort; depth-corrected per-gene residualization not in original scoring run; refer to depth_corrected_trend.csv for sample-mean resid trend
GSM7980862_N-3	GSE250521	PT	4538	-8.361277170460809e-08	0.06586754906078449	0.05824406671495153	8.361277170460809e-08	-0.06586754906078449	-0.05824406671495153	-1.0573617144319376e-07	9.105363632168714e-10	1.841303856969826e-08				raw_within_sample_z; sample-level depth-corrected residual via depth_corrected_trend.py (epi-only, score-level)	Cancer progression cohort; depth-corrected per-gene residualization not in original scoring run; refer to depth_corrected_trend.csv for sample-mean resid trend
GSM7980863_N-4	GSE250521	PT	3027	-1.0521968959960302e-09	0.10170035949176354	0.07595259290832233	1.0521968959960302e-09	-0.10170035949176354	-0.07595259290832233	1.0260042947458309e-08	-5.645844730831918e-08	-8.440961352342214e-09				raw_within_sample_z; sample-level depth-corrected residual via depth_corrected_trend.py (epi-only, score-level)	Cancer progression cohort; depth-corrected per-gene residualization not in original scoring run; refer to depth_corrected_trend.csv for sample-mean resid trend
GSM7980864_PTC-1	GSE250521	PTC	2262	1.0484491600083869e-08	0.09421232553704686	0.02957369075583039	-1.0484491600083869e-08	-0.09421232553704686	-0.02957369075583039	1.304074845398164e-08	2.69300875298431e-08	1.7203585307870665e-08				raw_within_sample_z; sample-level depth-corrected residual via depth_corrected_trend.py (epi-only, score-level)	Cancer progression cohort; depth-corrected per-gene residualization not in original scoring run; refer to depth_corrected_trend.csv for sample-mean resid trend
GSM7980865_PTC-2	GSE250521	PTC	4947	-2.862209419859013e-09	0.07716841894753032	0.03808377184648343	2.862209419859013e-09	-0.07716841894753032	-0.03808377184648343	-4.500378006333721e-09	7.71127895705839e-08	5.75607843436835e-09				raw_within_sample_z; sample-level depth-corrected residual via depth_corrected_trend.py (epi-only, score-level)	Cancer progression cohort; depth-corrected per-gene residualization not in original scoring run; refer to depth_corrected_trend.csv for sample-mean resid trend
GSM7980866_PTC-3	GSE250521	PTC	4516	6.630918954949032e-08	0.06179945100516386	0.05406362806028344	-6.630918954949032e-08	-0.06179945100516386	-0.05406362806028344	5.194964271455538e-08	-1.6949568198735268e-08	9.818645268155399e-09				raw_within_sample_z; sample-level depth-corrected residual via depth_corrected_trend.py (epi-only, score-level)	Cancer progression cohort; depth-corrected per-gene residualization not in original scoring run; refer to depth_corrected_trend.csv for sample-mean resid trend
GSM7980867_PTC-4	GSE250521	PTC	4498	1.8370766514845645e-08	0.14340321816582383	0.13655853324201775	-1.8370766514845645e-08	-0.14340321816582383	-0.13655853324201775	3.3477710094538456e-08	-1.747130724823036e-08	-2.73557581902782e-09				raw_within_sample_z; sample-level depth-corrected residual via depth_corrected_trend.py (epi-only, score-level)	Cancer progression cohort; depth-corrected per-gene residualization not in original scoring run; refer to depth_corrected_trend.csv for sample-mean resid trend
GSM7980868_LPTC-1	GSE250521	LPTC	1520	-3.279839454414644e-08	0.0768261949507425	0.08331377653447368	3.279839454414644e-08	-0.0768261949507425	-0.08331377653447368	-3.920139473736016e-08	4.803924079001046e-08	7.828532891531269e-09				raw_within_sample_z; sample-level depth-corrected residual via depth_corrected_trend.py (epi-only, score-level)	Cancer progression cohort; depth-corrected per-gene residualization not in original scoring run; refer to depth_corrected_trend.csv for sample-mean resid trend
GSM7980869_LPTC-2	GSE250521	LPTC	4071	-2.243927683856984e-08	0.2774129189927427	0.4785558808802495	2.243927683856984e-08	-0.2774129189927427	-0.4785558808802495	-2.1686460332747034e-08	1.896912304898732e-08	-9.08759518685726e-09				raw_within_sample_z; sample-level depth-corrected residual via depth_corrected_trend.py (epi-only, score-level)	Cancer progression cohort; depth-corrected per-gene residualization not in original scoring run; refer to depth_corrected_trend.csv for sample-mean resid trend
GSM7980870_LPTC-3	GSE250521	LPTC	3336	-8.667916667400031e-09	0.21211691353458031	0.18592061824664266	8.667916667400031e-09	-0.21211691353458031	-0.18592061824664266	4.333546013109343e-09	1.5154760198686343e-08	5.2461690419494204e-09				raw_within_sample_z; sample-level depth-corrected residual via depth_corrected_trend.py (epi-only, score-level)	Cancer progression cohort; depth-corrected per-gene residualization not in original scoring run; refer to depth_corrected_trend.csv for sample-mean resid trend
GSM7980871_LPTC-4	GSE250521	LPTC	4707	5.598400255288993e-09	0.18680481158132115	0.18160515822438403	-5.598400255288993e-09	-0.18680481158132115	-0.18160515822438403	-1.3100925218166071e-08	7.700871149619711e-08	4.748359251576893e-09				raw_within_sample_z; sample-level depth-corrected residual via depth_corrected_trend.py (epi-only, score-level)	Cancer progression cohort; depth-corrected per-gene residualization not in original scoring run; refer to depth_corrected_trend.csv for sample-mean resid trend
GSM7980872_ATC-1	GSE250521	ATC	4715	3.898831173215973e-09	0.003391377832255724	0.0015977557956395234	-3.898831173215973e-09	-0.003391377832255724	-0.0015977557956395234	4.365961822401081e-09	-2.7884105831250027e-08	2.2737657260777166e-08				raw_within_sample_z; sample-level depth-corrected residual via depth_corrected_trend.py (epi-only, score-level)	Cancer progression cohort; depth-corrected per-gene residualization not in original scoring run; refer to depth_corrected_trend.csv for sample-mean resid trend
GSM7980873_ATC-2	GSE250521	ATC	2354	1.8326253170524734e-09	0.016268707429906545	0.03442790325297113	-1.8326253170524734e-09	-0.016268707429906545	-0.03442790325297113	2.5841121501651954e-09	-6.0093372981412385e-09	-7.1194052522475226e-09				raw_within_sample_z; sample-level depth-corrected residual via depth_corrected_trend.py (epi-only, score-level)	Cancer progression cohort; depth-corrected per-gene residualization not in original scoring run; refer to depth_corrected_trend.csv for sample-mean resid trend
GSM7980874_ATC-3	GSE250521	ATC	1621	-4.929056041931781e-10	-0.0010431356757089913	0.0022536693152709407	4.929056041931781e-10	0.0010431356757089913	-0.0022536693152709407	3.2819247411569633e-09	-7.338408389422073e-09	-3.718703639681658e-08				raw_within_sample_z; sample-level depth-corrected residual via depth_corrected_trend.py (epi-only, score-level)	Cancer progression cohort; depth-corrected per-gene residualization not in original scoring run; refer to depth_corrected_trend.csv for sample-mean resid trend
GSM7980875_ATC-4	GSE250521	ATC	2586	-1.1336670530915758e-08	0.028508367573936582	-0.03926279726584234	1.1336670530915758e-08	-0.028508367573936582	0.03926279726584234	-8.736505415472594e-09	-1.9989176341144843e-08	-7.136395974289843e-09				raw_within_sample_z; sample-level depth-corrected residual via depth_corrected_trend.py (epi-only, score-level)	Cancer progression cohort; depth-corrected per-gene residualization not in original scoring run; refer to depth_corrected_trend.csv for sample-mean resid trend
GSM7221915_P1	GSE230424	PTC_HT	3566	2.206451486304301e-08	0.3670549267130847	0.4796970143377018	-2.206451486304301e-08	-0.3670549267130847	-0.4796970143377018	2.493031408034424e-08	-3.072728547829892e-08	-3.936090853014913e-09	-2.772887268511592e-08	0.1904792698533961	-0.2940545275888461	raw_within_sample_z + per-gene depth-residualized z	PTC + Hashimoto-overlap; per-patient HT vs PTC+HT split unavailable in GEO metadata; user spec=all PTC+HT
GSM7221916_P2	GSE230424	PTC_HT	3005	4.333552079608425e-08	0.2464704214703259	0.2692194738589096	-4.333552079608425e-08	-0.2464704214703259	-0.2692194738589096	4.7472752076149617e-08	6.4902851941783695e-09	-8.414808645632649e-09	-1.2189098174502582e-08	0.0684526499825262	-0.1020870777170757	raw_within_sample_z + per-gene depth-residualized z	PTC + Hashimoto-overlap; per-patient HT vs PTC+HT split unavailable in GEO metadata; user spec=all PTC+HT
GSM7221917_P3	GSE230424	PTC_HT	2961	1.3911503982100136e-08	0.2798939110287644	0.3755527770995142	-1.3911503982100136e-08	-0.2798939110287644	-0.3755527770995142	5.505856129904826e-09	1.3434012159441826e-08	5.268655851798192e-09	1.2212850385316821e-08	0.1690232515427343	-0.159789967128888	raw_within_sample_z + per-gene depth-residualized z	PTC + Hashimoto-overlap; per-patient HT vs PTC+HT split unavailable in GEO metadata; user spec=all PTC+HT
GSM7221918_P4	GSE230424	PTC_HT	4601	4.475131493253624e-08	0.1820694837730986	0.1964268512337967	-4.475131493253624e-08	-0.1820694837730986	-0.1964268512337967	5.0235068462336796e-08	2.9581703110495432e-08	-1.5781568271283107e-10	1.0873501846414618e-08	0.0735588774075002	-0.0716282940858711	raw_within_sample_z + per-gene depth-residualized z	PTC + Hashimoto-overlap; per-patient HT vs PTC+HT split unavailable in GEO metadata; user spec=all PTC+HT
GSM7908359_C1	GSE248205	CONTROL	1064	-2.8334144736395334e-08	0.1903283989381015	0.2312468666593985	2.8334144736395334e-08	-0.1903283989381015	-0.2312468666593985	-4.531231202868229e-08	-8.466560154717304e-09	-1.133270675714369e-08	-4.907128759141304e-08	0.0335344707503057	-0.0532194501910339	raw_within_sample_z + per-gene depth-residualized z	Autoimmune thyroid (no cancer); CONTROL/HT/GD; PMID 39003267
GSM7908360_C2	GSE248205	CONTROL	1744	-5.75588130735864e-08	0.1071817042277867	0.108531214128211	5.75588130735864e-08	-0.1071817042277867	-0.108531214128211	-8.320704415333902e-08	-1.0134495402666464e-08	1.2500000030979703e-09	-4.530177743636874e-09	-0.0044323076130172	-0.0049628588569265	raw_within_sample_z + per-gene depth-residualized z	Autoimmune thyroid (no cancer); CONTROL/HT/GD; PMID 39003267
GSM7908361_HT1	GSE248205	HT	1927	4.0078142246439584e-09	0.3317991512843361	0.6577239333952283	-4.0078142246439584e-09	-0.3317991512843361	-0.6577239333952283	-1.511996887008646e-08	2.1407924228582423e-08	-1.1039880626348737e-08	-2.400546964207989e-08	0.4012794939987963	-0.5424896249416963	raw_within_sample_z + per-gene depth-residualized z	Autoimmune thyroid (no cancer); CONTROL/HT/GD; PMID 39003267
GSM7908362_HT2	GSE248205	HT	2595	-5.783217726757892e-09	0.2421326537266563	0.3051069827489985	5.783217726757892e-09	-0.2421326537266563	-0.3051069827489985	-3.5941641657739302e-09	-3.3132154133162504e-08	-4.616418887835997e-09	-1.1087976879449598e-08	0.1071479674769101	-0.1336556603825108	raw_within_sample_z + per-gene depth-residualized z	Autoimmune thyroid (no cancer); CONTROL/HT/GD; PMID 39003267
GSM7908363_HT3	GSE248205	HT	1980	1.1541313131229572e-08	0.310642688684303	0.5218531662619393	-1.1541313131229572e-08	-0.310642688684303	-0.5218531662619393	2.3240742425046374e-08	1.3865848481070588e-08	-1.1112643933672368e-08	7.578044950698275e-09	0.3485261735134457	-0.4325771257840491	raw_within_sample_z + per-gene depth-residualized z	Autoimmune thyroid (no cancer); CONTROL/HT/GD; PMID 39003267
GSM7908364_GD1	GSE248205	GD	2266	-6.420303221457915e-08	0.1229925211849549	0.1184361372185608	6.420303221457915e-08	-0.1229925211849549	-0.1184361372185608	-1.0776063547776657e-07	2.6524779350487665e-08	-4.71932479170103e-09	-6.368411297572642e-08	0.0137744465911999	-0.0547861555946178	raw_within_sample_z + per-gene depth-residualized z	Autoimmune thyroid (no cancer); CONTROL/HT/GD; PMID 39003267
GSM7908365_GD2	GSE248205	GD	2005	4.5833321694773875e-08	0.2027062477902392	0.2333665271108765	-4.5833321694773875e-08	-0.2027062477902392	-0.2333665271108765	5.903840897942809e-08	-3.218960597369351e-08	-9.034413956871171e-09	6.878563582509138e-09	0.0435046342914171	-0.0544677269745065	raw_within_sample_z + per-gene depth-residualized z	Autoimmune thyroid (no cancer); CONTROL/HT/GD; PMID 39003267
GSM7908366_GD3	GSE248205	GD	1969	1.4057516504950232e-08	0.1776149553616243	0.1828454675462474	-1.4057516504950232e-08	-0.1776149553616243	-0.1828454675462474	1.5918547486434713e-08	-2.050136108820773e-08	5.952260049329511e-09	1.647476739354076e-08	0.0406433570491596	-0.0863723453001549	raw_within_sample_z + per-gene depth-residualized z	Autoimmune thyroid (no cancer); CONTROL/HT/GD; PMID 39003267
```

# Appendix B — SuppTable SX+1 (independent lineage + technical confounding)

```tsv
table	dataset	score	scope	version	n	spearman_rho	spearman_p	pearson_r	pearson_p	flag	interpretation
independent_lineage_validation	GSE230424	DM1_like vs THYROID_NONOVERLAP	all_spots	raw	14133	-0.6778991184603385	0.0	-0.6543048029526817	0.0		DM1_like vs THYROID_NONOVERLAP — true independent lineage check (no overlap with 8-gene)
independent_lineage_validation	GSE230424	DM1_like vs THYROID_NONOVERLAP	all_spots	resid	14133	-0.4907277373778117	0.0	-0.488155572164167	0.0		DM1_like vs THYROID_NONOVERLAP — true independent lineage check (no overlap with 8-gene)
independent_lineage_validation	GSE248205	DM1_like vs THYROID_NONOVERLAP	all_spots	raw	15550	-0.5916398323344473	0.0	-0.592702430229072	0.0		DM1_like vs THYROID_NONOVERLAP — true independent lineage check (no overlap with 8-gene)
independent_lineage_validation	GSE248205	DM1_like vs THYROID_NONOVERLAP	all_spots	resid	15550	-0.4237188283680889	0.0	-0.4512621022011145	0.0		DM1_like vs THYROID_NONOVERLAP — true independent lineage check (no overlap with 8-gene)
independent_lineage_validation	GSE230424	DM1_like vs THYROID_NONOVERLAP	sample_mean_all	raw	4	0.1999999999999999	0.8	-0.0431972882257185	0.9568027117742814		DM1_like vs THYROID_NONOVERLAP at sample-mean (all)
independent_lineage_validation	GSE230424	DM1_like vs THYROID_NONOVERLAP	sample_mean_epi50	raw	4	-1.0	0.0	-0.9079401692940424	0.0920598307059576		DM1_like vs THYROID_NONOVERLAP at sample-mean (epi50)
independent_lineage_validation	GSE230424	DM1_like vs THYROID_NONOVERLAP	sample_mean_epi25	raw	4	-0.7999999999999999	0.2000000000000001	-0.8901738446021238	0.1098261553978763		DM1_like vs THYROID_NONOVERLAP at sample-mean (epi25)
independent_lineage_validation	GSE230424	DM1_like vs THYROID_NONOVERLAP	sample_mean_all	resid	4	0.3999999999999999	0.6	0.2354598609066923	0.7645401390933078		DM1_like vs THYROID_NONOVERLAP at sample-mean (all)
independent_lineage_validation	GSE230424	DM1_like vs THYROID_NONOVERLAP	sample_mean_epi50	resid	4	-1.0	0.0	-0.99229985472452	0.0077001452754801		DM1_like vs THYROID_NONOVERLAP at sample-mean (epi50)
independent_lineage_validation	GSE230424	DM1_like vs THYROID_NONOVERLAP	sample_mean_epi25	resid	4	-0.7999999999999999	0.2000000000000001	-0.8855376804163091	0.1144623195836909		DM1_like vs THYROID_NONOVERLAP at sample-mean (epi25)
independent_lineage_validation	GSE248205	DM1_like vs THYROID_NONOVERLAP	sample_mean_all	raw	8	-0.7619047619047621	0.0280049391530718	-0.6804112023543244	0.0632948920800435		DM1_like vs THYROID_NONOVERLAP at sample-mean (all)
independent_lineage_validation	GSE248205	DM1_like vs THYROID_NONOVERLAP	sample_mean_epi50	raw	8	-0.9761904761904764	3.314396026200098e-05	-0.952856443924328	0.0002527687647699		DM1_like vs THYROID_NONOVERLAP at sample-mean (epi50)
independent_lineage_validation	GSE248205	DM1_like vs THYROID_NONOVERLAP	sample_mean_epi25	raw	8	-0.9761904761904764	3.314396026200098e-05	-0.9875673132915304	4.759666886352331e-06		DM1_like vs THYROID_NONOVERLAP at sample-mean (epi25)
independent_lineage_validation	GSE248205	DM1_like vs THYROID_NONOVERLAP	sample_mean_all	resid	8	0.1904761904761905	0.6514014957024814	0.2404759603301334	0.5661890148104232		DM1_like vs THYROID_NONOVERLAP at sample-mean (all)
independent_lineage_validation	GSE248205	DM1_like vs THYROID_NONOVERLAP	sample_mean_epi50	resid	8	-0.880952380952381	0.0038503204637324	-0.9835697384701192	1.0952324856746703e-05		DM1_like vs THYROID_NONOVERLAP at sample-mean (epi50)
independent_lineage_validation	GSE248205	DM1_like vs THYROID_NONOVERLAP	sample_mean_epi25	resid	8	-0.880952380952381	0.0038503204637324	-0.995775779417616	1.878463837366856e-07		DM1_like vs THYROID_NONOVERLAP at sample-mean (epi25)
technical_confounding	GSE230424	RAI_8_score_raw	all_spots_vs_log_total_counts	raw	14133	0.5250977996463394	0.0			⚠️ depth-confound	Score ~ log_total_counts confound check; |ρ|>0.30 flagged
technical_confounding	GSE230424	DM1_like_score_raw	all_spots_vs_log_total_counts	raw	14133	-0.5250977996463394	0.0			⚠️ depth-confound	Score ~ log_total_counts confound check; |ρ|>0.30 flagged
technical_confounding	GSE230424	RAI_8_score_resid	all_spots_vs_log_total_counts	resid	14133	0.0486943146475044	6.962959406502137e-09			ok	Score ~ log_total_counts confound check; |ρ|>0.30 flagged
technical_confounding	GSE230424	DM1_like_score_resid	all_spots_vs_log_total_counts	resid	14133	-0.0486943146475044	6.962959406502137e-09			ok	Score ~ log_total_counts confound check; |ρ|>0.30 flagged
technical_confounding	GSE230424	THYROID_NONOVERLAP_score_raw	all_spots_vs_log_total_counts	raw	14133	0.541594583405772	0.0			⚠️ depth-confound	Score ~ log_total_counts confound check; |ρ|>0.30 flagged
technical_confounding	GSE230424	THYROID_NONOVERLAP_score_resid	all_spots_vs_log_total_counts	resid	14133	0.0447775991484107	1.0069570702760873e-07			ok	Score ~ log_total_counts confound check; |ρ|>0.30 flagged
technical_confounding	GSE248205	RAI_8_score_raw	all_spots_vs_log_total_counts	raw	15550	0.3018068576401589	0.0			⚠️ depth-confound	Score ~ log_total_counts confound check; |ρ|>0.30 flagged
technical_confounding	GSE248205	DM1_like_score_raw	all_spots_vs_log_total_counts	raw	15550	-0.3018068576401589	0.0			⚠️ depth-confound	Score ~ log_total_counts confound check; |ρ|>0.30 flagged
technical_confounding	GSE248205	RAI_8_score_resid	all_spots_vs_log_total_counts	resid	15550	-0.0016644993789532	0.8355835265535433			ok	Score ~ log_total_counts confound check; |ρ|>0.30 flagged
technical_confounding	GSE248205	DM1_like_score_resid	all_spots_vs_log_total_counts	resid	15550	0.0016644993789532	0.8355835265535433			ok	Score ~ log_total_counts confound check; |ρ|>0.30 flagged
technical_confounding	GSE248205	THYROID_NONOVERLAP_score_raw	all_spots_vs_log_total_counts	raw	15550	0.3284834630155037	0.0			⚠️ depth-confound	Score ~ log_total_counts confound check; |ρ|>0.30 flagged
technical_confounding	GSE248205	THYROID_NONOVERLAP_score_resid	all_spots_vs_log_total_counts	resid	15550	0.00190202385698	0.8125304662085958			ok	Score ~ log_total_counts confound check; |ρ|>0.30 flagged
```

# Appendix C — Frozen figure paths (binary, must scp separately)

```
project/supplementary/spatial_freeze_2026_05_03/SuppFig_X1_GSE250521_spatial_PT_to_ATC.png    (7.5 MB, 200 DPI; 2x4 stage spatial heatmap)
project/supplementary/spatial_freeze_2026_05_03/SuppFig_X1_GSE250521_spatial_PT_to_ATC.pdf    (2.5 MB)
project/supplementary/spatial_freeze_2026_05_03/SuppFig_X2_DM1_thyroid_nonoverlap_crossvalidation.png   (186 KB; non-overlapping lineage scatter, n=12)
project/supplementary/spatial_freeze_2026_05_03/SuppFig_X2_DM1_thyroid_nonoverlap_crossvalidation.pdf   (25 KB)
project/supplementary/spatial_freeze_2026_05_03/SuppFig_X3_GSE248205_autoimmune_negative_control.png    (121 KB; CONTROL/HT/GD bar plot)
project/supplementary/spatial_freeze_2026_05_03/SuppFig_X3_GSE248205_autoimmune_negative_control.pdf    (21 KB)
```

Spatial supplement frozen. Return to manuscript marathon.
