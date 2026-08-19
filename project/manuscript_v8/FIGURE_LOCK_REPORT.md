---
title: "Figure Lock Report — DM1 Paper 1"
date: 2026-07-30
auditor: Claude figure-lock agent
status: audit complete; no figures generated; no prose finalized
scope: Phase A inventory / Phase B architecture / Phase C missing panels / Phase D captions
---

# Figure Lock Report — DM1 Paper 1

## Phase A — Asset Inventory

### Locations searched

| Location | Exists | Finding |
|---|---|---|
| `project/manuscript_v8/figures/` | Yes | Only 2 files: `Fig7_dm1_mechanism.png/pdf` and `Fig8_epigenetic.png/pdf` — both belong to the **old** multi-figure v2 draft, not the 5-figure reviewer-first plan. |
| `project/manuscript_v8/assets/p1_onepage_audit/` | Yes | 30+ PNG files for web dossier only (fig01–fig18, deepdive_*, calibration, roc). None are journal-submission figures. |
| `project/manuscript_v8/figures_for_advisor/` | Does not exist. | |
| `project/dm1_story_web/public/figures/` | Yes | 40+ PNG; includes `Fig1_annotated.png`, `Fig3_annotated.png`, `Fig5_annotated.png`, `Fig6_annotated.png`, `EV_*` series, and web-only analysis charts. These are web-preview versions, not 600-dpi submission assets. |
| `project/dm1_story_web/dist/figures/` | Yes | Mirror of public/figures plus `manuscripts/figures/` subdirectory with `Fig1_discovery_axis.png` through `Fig6_reflex_translation.png` — same render as below. |
| `project/dm1_story_web/dist/manuscripts/figures/` | Yes | **Primary rendered main-figure PNGs (Fig1–Fig6) and PDFs** for the 6-figure NC architecture. This is the canonical render location. |
| `project/results/manuscript_v8_nc_main/` | Yes | **Most complete figure directory.** Contains Fig1–Fig6 (PDF+PNG), ED1–ED15 (PNG), build scripts (`fig1_discovery.py` through `fig6_reflex.py`, `build_ed.py`, `build_panels.py`, `annotate_figs.py`, `master_tcga.tsv`). |
| `project/results/` (fig*.png/pdf pattern) | Yes | Hundreds of analysis-intermediate figure files across dm1_robustness_v2026_05_08/rounds, dark_matter_phase2, v17_realfix, audit_2026_04_29. Useful as source panels but not final submission figures. |

### Canonical submission-ready figure set (as of 2026-07-30)

All 6 main figure composites exist as PDF + PNG in `project/results/manuscript_v8_nc_main/`:

- `Fig1_discovery_axis.pdf/.png`
- `Fig2_fusion_mechanism.pdf/.png`
- `Fig3_epigenetic.pdf/.png`
- `Fig4_sc_validation.pdf/.png`
- `Fig5_survival_portability.pdf/.png`
- `Fig6_reflex_translation.pdf/.png`

All 15 Extended Data figures exist as PNG in `project/results/manuscript_v8_nc_main/`:
ED1 through ED15.

**However, several individual panels within these composites are rendered with placeholder or synthetic data rather than real analysis outputs (see Phase C).**

### Orphan files requiring disposition

- `project/manuscript_v8/figures/Fig7_dm1_mechanism.png/.pdf` — belongs to the old v2 multi-figure architecture. No Fig. 7 exists in either the 5-figure reviewer-first plan. Status: REMOVE from submission package.
- `project/manuscript_v8/figures/Fig8_epigenetic.png/.pdf` — same issue. Status: EDITORIAL_DECISION_REQUIRED (content partially overlaps reviewer-first Fig. 3; author must decide whether any panels are promoted to the 5-figure composite or removed).

---

## Phase B — Architecture Conflict Resolution

### Comparison table

| Aspect | 5-figure reviewer-first plan (REVIEWER_FIRST_NC_REBUILD_2026_07_08.md) | 6-figure v2 NC draft (05_figure_captions_NC.md + NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md) |
|---|---|---|
| Main figure count | 5 | 6 |
| Panel count per figure | ~5 panels each | 7 panels each |
| Supplementary/ED figures | 6 Extended Data | 15 Extended Data |
| Story structure | State discovery → Biology → Methylation → Validation → Clinical | State discovery → Fusion/histology → Epigenetic → SC validation → Survival → Reflex/clinical |
| Fusion/kinase content | Fig1E (prevalence only) | Full Fig. 2 (7 panels, OR=7.41) |
| Single-cell content | Referenced but not a main figure | Full Fig. 4 (7 panels, SC UMAP + pseudotime) |
| Clinical emphasis | Harm avoidance; ATA uncertainty mosaic | Predictive biomarker (BRAF+/PFI interaction p=0.022) + selpercatinib waterfall |
| IHC 3-plex + power simulation | Deleted to supplement | Fig. 8 in v2 (4 panels; Fig. 8a–8d) |
| Interaction claim (Fig. 7) | Not present | Fig. 7a–7d: DM1 × BRAF interaction, head-to-head, multi-modal concordance |
| NC compliance fit | Conservative; strong fit for NC Article | Also fits NC Article; risk of reviewer attack on interaction claim without prospective validation |
| Reviewer risk | Lower — anticipates key attacks | Higher — "predictive biomarker" language without pretreatment RAI measurement is a direct reviewer vulnerability |
| Missing assets to generate | Fig3C, Fig5C (KM stack), Fig5F (time-ROC), Fig5G (calibration) + illustrative panels | Same missing assets + all of Fig.7 interaction panels + Fig.8 IHC/power + multiple SC panels in Fig.4 |
| Concordance with DM1_PROJECT_STATE.md | Directly cited as "figure architecture note" | Cited as "latest known draft"; conflicts called out in known narrative conflicts section |
| Current rendered output name | Does not exist as a compiled set | `project/results/manuscript_v8_nc_main/` Fig1–Fig6 follows this 6-figure plan |

### Architecture recommendation

**Adopt the 6-figure reviewer-first plan from `05_figure_captions_NC.md` with the five boundary constraints already locked in that document (no PRISM/DepMap, no TERT survival, no Mun proteome, no pan-cancer, no sub-B HLA) as the single authoritative architecture.** Retire the `REVIEWER_FIRST_NC_REBUILD_2026_07_08.md` 5-figure plan from active consideration.

Rationale: (1) The `manuscript_v8_nc_main/` build already renders 6 main figures and 15 ED figures matching the `05_figure_captions_NC.md` structure — implementing the 5-figure plan would require re-architecting all scripts, renaming panels, and rewriting captions from scratch. (2) The 6-figure plan is more defensible because Fig. 2 (fusion enrichment OR=7.41) and Fig. 4 (SC thyrocyte-intrinsic validation) are substantive independent claims, not just reviewer footnotes; excising them weakens the manuscript rather than tightening it.

The `REVIEWER_FIRST_NC_REBUILD_2026_07_08.md` document retains value as a reviewer-attack checklist but must no longer be treated as an alternative figure architecture. `DM1_PROJECT_STATE.md` must be updated to remove the phrase "5-figure reviewer-first plan" as the active architecture.

**Fig. 7 and Fig. 8** from the old `05_figure_captions.md` Cell Press companion build: Fig. 7 (interaction + replication) and Fig. 8 (IHC 3-plex + power) are referenced in the v2 results text but do NOT appear in `05_figure_captions_NC.md`. This is the primary unresolved numbering conflict. Resolution: treat all v2 text references to "Fig. 7" and "Fig. 8" as placeholder citations that map to content either already incorporated into Fig. 2 / Fig. 3 / Fig. 6 (for the interaction and IHC claims) or demoted to Extended Data. An explicit editorial log entry is required before the author's final copy-edit pass.

---

## Phase C — Missing Panels and Generation Plans

### Panel 1: Fig. 3C — Beta vs expression scatter for TPO, DIO1, TSHR, TG

**Status: REBUILD_REQUIRED (currently rendered with synthetic noise)**

Current state: `project/results/manuscript_v8_nc_main/fig3_epigenetic.py` lines 59–71 generate the 4-panel scatter using `expr = (1 - beta) + rng.normal(0, 0.08, len(beta))` — a synthetic expression proxy, not real RNA-seq values. This means the current `Fig3_epigenetic.png` panel C shows fabricated data.

Data needed:
- Beta (per-sample, per-gene): EXISTS at `project/results/audit_2026_04_30/round5/r5_2_sample_methylation_8gene.tsv` (n=503 samples, columns: `sample_short`, `DIO1`, `FOXE1`, `NKX2-1`, `PAX8`, `SLC5A5`, `TG`, `TPO`, `TSHR`, `mean_8g_beta`)
- RNA expression (per-sample, per-gene): EXISTS at `project/results/ncomm_push_2026_05_08/cbioportal_sweep/panel_expression_thpa_tcga_gdc.tsv` (n=513 samples, columns: `sampleId`, `DIO1`, `FOXE1`, `NKX2-1`, `PAX8`, `SLC5A5`, `TG`, `TPO`, `TSHR`, `RAI_8`, `DM1_like`, `DM_call`; values appear to be log2-transformed expression)
- Join key: `sample_short` (12-char TCGA barcode) from beta table vs `sampleId` (15-char, e.g. `TCGA-4C-A93U-01`) in expression table — truncate to 12 chars for join
- Paired n after join: expect approximately 490–503

Script action: Edit `fig3_epigenetic.py` lines 56–71 to load the expression table, join on truncated barcode, and plot real `(beta, log2_expr)` scatter per gene. Spearman ρ must be computed from real data, not re-used from the caption.

Estimated effort: **moderate** — data exists, join is straightforward, Spearman ρ needs recomputation and verification against the d-values in the caption.

### Panel 2: Fig. 3D (caption label F) — Per-driver-class mean 8-gene beta

**Status: REBUILD_REQUIRED (currently hardcoded placeholder values)**

Current state: `project/results/manuscript_v8_nc_main/fig3_epigenetic.py` lines 116–140 use hardcoded beta values `[0.37, 0.39, 0.38, 0.36, 0.27, 0.31]` for `BRAF V600E / RET fusion / NTRK fusion / ALK fusion / RAS-mutant / driver-neg`. These match the caption numbers but are not derived from a live computation joining the beta table to the driver-class column.

Data needed:
- Beta per sample: `project/results/audit_2026_04_30/round5/r5_2_sample_methylation_8gene.tsv`
- Driver class per sample: `project/results/manuscript_v8_nc_main/master_tcga.tsv` (column `driver_class`)
- Join: `sample_short` (consistent key in both)

The caption states "BRAF V600E β = 0.37 ≈ RET fusion β = 0.39 ≫ RAS-mutant β = 0.27." These specific numbers must be verified to emerge from the real joined table before locking the panel.

Estimated effort: **trivial** — data is joined, computation is a groupby mean. Risk is low but verification is required.

### Panel 3: Fig. 5C — Multi-cohort Kaplan-Meier stack (TCGA / MSK / pooled)

**Status: REBUILD_REQUIRED (script does not yet exist)**

Current state: The `fig5_survival.py` script exists but does not contain three-panel side-by-side KM code. The captions document (line 167) notes "C multi-cohort KM stack" as new build required and names the script `fig5c_multicohort_km_stack.py` — this script does not exist in `manuscript_v8_nc_main/`.

Data needed:
- TCGA OS: `project/results/manuscript_v8_nc_main/master_tcga.tsv` (columns `os_days`, `os_event`, `DM`)
- MSK-IMPACT OS: source file not confirmed in project; `project/results/audit_2026_04_30/round5/r5_1_msk_fusion.json` may contain summary statistics only, not per-patient survival data. **MSK per-patient OS table must be located before this panel can be generated.**

Estimated effort: **involved** if MSK per-patient data is missing; **moderate** if it exists as a TSV.

### Panel 4: Fig. 5F — Time-dependent ROC (1/3/5-year, IPCW)

**Status: REBUILD_REQUIRED (script does not yet exist)**

Script `fig5f_time_dependent_roc.py` noted in captions but absent from `manuscript_v8_nc_main/`. Data is available from `master_tcga.tsv`. Requires `lifelines` or `scikit-survival` IPCW implementation.

Estimated effort: **moderate**.

### Panel 5: Fig. 5G — Calibration plot (5-year, Hosmer-Lemeshow)

**Status: REBUILD_REQUIRED (script does not yet exist)**

Script `fig5g_calibration.py` noted in captions but absent. Same data source as 5F. TCGA has only 16 OS events — calibration with decile binning on 16 events is statistically questionable. **Editorial decision required on whether to keep this panel or replace with a different reliability visualization.**

Estimated effort: **moderate** if kept; author confirmation needed.

### Panel 6: Fig. 6A — Reflex-testing algorithm flowchart

**Status: REBUILD_REQUIRED (script does not yet exist)**

Script `fig6a_reflex_flowchart.py` noted in captions but absent. This is a diagram schematic, not a data-driven plot. Can be generated with matplotlib patches or a diagram library. The population yield estimate (48 selpercatinib-eligible per 1,000 PTC patients) is derived from TCGA proportions and requires explicit notation that it is a population projection, not a trial result.

Estimated effort: **trivial** (diagram only).

### Panel 7: Fig. 6D — Post-RAI refractory boxplot (GSE151179)

**Status: REBUILD_REQUIRED (panel in Fig6_reflex_translation.png is composite placeholder)**

An asset exists at `project/dm1_story_web/public/figures/EV_rai_lineage_box.png` that may serve this purpose, but it needs to be confirmed as the correct GSE151179 pre/post-RAI thyroid differentiation score comparison and re-rendered at submission resolution. Script `fig6d_post_rai_box.py` noted but absent from `manuscript_v8_nc_main/`.

Estimated effort: **trivial** if `EV_rai_lineage_box.png` is validated; **moderate** if GSE151179 data needs reprocessing.

### Panel 8: Fig. 6E — Decision-curve analysis

**Status: REBUILD_REQUIRED (script does not yet exist)**

Script `fig6e_decision_curve.py` noted but absent. DCA code (`fig_decision.png` and `fig_calibration.png`) exists in `project/manuscript_v8/assets/p1_onepage_audit/` — these may be salvageable with relabeling. Requires verification that the "DM1-stratified" strategy is properly defined (threshold range, endpoint, n).

Estimated effort: **moderate**.

### Missing source data (panels that cannot be generated without data acquisition)

| Panel | Missing data | Acquisition path |
|---|---|---|
| Fig. 2G | MSK-IMPACT per-patient fusion × DM table | Must confirm MSK-IMPACT data license and locate processed table in project |
| Fig. 4C | GSE241184 processed single-cell data | GSE241184 not found under `project/data/` or `project/results/`; requires GEO download |
| Fig. 4D | GSE232237 processed single-cell data | GSE232237 not found; requires GEO download |
| ED4 | Same as Fig. 4C/4D above | Same |
| ED11 | Mun 2025 proteogenomics (HRA004166) | `project/manuscript_biorxiv_2026_05_20/supp_data_bundle/04_Mun_2025_protein_unsup.tsv` exists (check); but the ED11 build in `build_ed.py` may not reference it |
| ED15 / Fig. 3G | GSE76039 (Landa 2016) processed expression matrix | `project/results/dm1_robustness_v2026_05_08/round4/landa_box.png` exists suggesting GSE76039 was processed; the raw matrix path used by that script must be traced to confirm availability |

**Note on Mun 2025:** `project/manuscript_biorxiv_2026_05_20/supp_data_bundle/04_Mun_2025_protein_unsup.tsv` was found during search — this may be the source for ED11. The `build_ed.py` script must be audited to confirm it references this file.

**Note on GSE76039 (Landa 2016):** The file `project/results/dm1_robustness_v2026_05_08/round4/landa_box.png` exists, suggesting the dataset was processed. The script `project/results/dm1_robustness_v2026_05_08/round4/` should be examined to recover the input path.

---

## Phase D — Caption Issues

The following gaps were identified in `05_figure_captions_NC.md`. Each caption is evaluated against four criteria: exact n, test+direction, error bars, source data citation.

### Figure 1

| Panel | n stated? | Test+direction? | Error bars? | Source data? | Gap |
|---|---|---|---|---|---|
| 1A | Yes (6 cohorts stated) | N/A (schematic) | N/A | N/A | No gap |
| 1B | No n for TIERA67 total | N/A (diagram) | N/A | Yoo 2016 cited | n for each TIERA67 category is missing |
| 1C | Yes (n=504) | KMeans k=2 stated | N/A (cluster colors) | Not specified | Source data file not named |
| 1D | Yes (n=504) | Tracks described | N/A | Not specified | Source data file not named; ATA 2015 annotation source not cited |
| 1E | Yes per driver | d and AUC stated; Mann-Whitney p for BRAF | N/A (bar plot) | Not specified | Not all driver p-values given (only BRAF p=0.57; TERT/KRAS/NRAS/HRAS AUC listed but p not given) |
| 1F | Yes (6 sizes) | ARI defined vs reference | N/A | Not specified | Hypergeometric p=3e-4 stated but one-sided vs two-sided not specified |
| 1G | Yes (n=143/361; 16 events) | Log-rank p stated; HR cross-ref to Fig5A | N/A (KM) | Not specified | Log-rank p value itself not given in caption (says "reported on plot") |
| Statistics block | Yes (n_TCGA=504) | Methods stated | N/A | N/A | Complete |

### Figure 2

| Panel | n stated? | Test+direction? | Error bars? | Source data? | Gap |
|---|---|---|---|---|---|
| 2A | Yes (180 dark-matter; 131 rescued) | % stated | N/A (Sankey) | Not specified | Missing Fisher test for Sankey flow |
| 2B | Yes (91 DM1; sub-A 72; sub-B 19) | Silhouette score 0.584 | N/A | Not specified | No p-value for sub-cluster structure |
| 2C | Yes (DM1: 63/82; DM2: 30.9%) | Fisher OR=7.41 CI stated; p=1.9e-13 | N/A (bar) | Source not named | Source data file not named |
| 2D | Yes (by partner class) | Counts given | N/A | cBioPortal cited | No test for partner distribution |
| 2E | Yes (27/33 = 81.8%) | Counts and rate stated | N/A | cBioPortal cited | No test for DM1 enrichment per fusion class |
| 2F | No exact n for FVPTC cells | Fisher OR=17.9 p=3e-31 | N/A | Not specified | n for FVPTC vs non-FVPTC subgroups missing |
| 2G | Yes (n=117 MSK) | Counts given | N/A | Not specified | No Fisher OR or p for MSK replication |

### Figure 3

| Panel | n stated? | Test+direction? | Error bars? | Source data? | Gap |
|---|---|---|---|---|---|
| 3A | Yes (n=503) | Cohen's d and Mann-Whitney p per gene | N/A (heatmap) | HM450 named | 2 control genes (DIO2, SLC26A4) referenced but no d/p reported for them |
| 3B | Yes (n=503) | DM1=0.385, DM2=0.253; direction stated | 95% bootstrap CI stated | HM450 named | "not_DM" group value (0.356) included in caption but n for not_DM not given |
| 3C | Yes (n=503 implied) | "Spearman ρ and within-DM-stratum slope reported" — but no values given | N/A (scatter) | Source not named | **Specific Spearman ρ per gene not in caption** (says "reported on plot"); source data file not named; **this is a BUILD-REQUIRED panel with synthetic data currently** |
| 3D | Yes (pooled n=1287; per-cohort n stated) | Spearman ρ per cohort; CI; pooled ρ=-0.327 | CI stated | Not named | I²=72.9% heterogeneity reason stated; Q p=0.005 given; complete |
| 3E | Yes (TCGA n=561; per-driver n) | Cohen's d Panel vs TDS; ΔAUC | N/A | Not specified | Source data file not named |
| 3F | Yes (BRAF, RET, RAS proportions implied) | Mean β values stated; direction stated | Not mentioned | Not specified | No CIs or SE for the per-driver-class mean β values; no n per driver class given |
| 3G | Yes (84 PDTC + 33 ATC) | 5/8 gene overlap stated | N/A (heatmap) | GSE76039 cited | No direction consistency test for the convergence; depends on visual comparison only |

### Figure 4

| Panel | n stated? | Test+direction? | Error bars? | Source data? | Gap |
|---|---|---|---|---|---|
| 4A | Yes (n=23 samples) | Direction stated (gradient present) | N/A | GSE193581 cited | No numerical score range or KS test for within-thyrocyte gradient |
| 4B | Yes (n=6 patients) | Spearman r range 0.798–0.886; p<10-10 Bonferroni | N/A | GSE184362 cited | n_cells per patient not given (caption says "n_cells annotation" on plot) |
| 4C | n not given for GSE241184 | Spearman r and BA residuals stated | N/A | GSE241184 cited | **n for GSE241184 is unlisted in caption** |
| 4D | n not given for GSE232237 | Direction stated | N/A | GSE232237 cited | **n for GSE232237 is unlisted** |
| 4E | Implied by 4A+4B | Regression slope and 95% CI | 95% CI stated | GSE184362+GSE193581 | Pooled n not explicit |
| 4F | Yes (TCGA n=513; Lee n=632) | Monotonic |ρ|≥0.95 stated | N/A | Lu 2023 reference cells n=67,678 stated | Complete |
| 4G | Cohort n stated (some unlisted) | Direction consistency ratio 6/7 | N/A | Not specified | GSE241184 and GSE232237 n still missing |

### Figure 5

| Panel | n stated? | Test+direction? | Error bars? | Source data? | Gap |
|---|---|---|---|---|---|
| 5A | Yes (n=504+117; 16+38 events) | HR with 95% CI; DL random-effects; I²=0% | Diamond = inv-variance | Not specified | Source data file not named |
| 5B | n implied from 5A | HR with CI per covariate | CI stated (forest) | Not specified | Model formula not fully stated (covariates listed but interaction with TERT unclear) |
| 5C | Yes (n=504 and 117) | Log-rank p; HR; CI | N/A (KM) | Not specified | **Script does not yet exist; MSK per-patient data not confirmed** |
| 5D | Yes (n=504+632+260+117) | KS p=0.44 for FFPE/FF | N/A (density) | Not specified | Source data file not named for each cohort |
| 5E | Yes (FFPE n=632, FF n=504) | KS p=0.44; Bland-Altman | N/A (density) | GSE213647 cited (FFPE) | Per-gene BA residuals mentioned but not quantified in caption |
| 5F | Not yet drafted | Not yet drafted | Not yet drafted | Not yet drafted | **Panel does not exist** |
| 5G | Not yet drafted | Not yet drafted | Not yet drafted | Not yet drafted | **Panel does not exist** |

### Figure 6

| Panel | n stated? | Test+direction? | Error bars? | Source data? | Gap |
|---|---|---|---|---|---|
| 6A | Yield estimate given (48/1,000) | N/A (flowchart) | N/A | TCGA proportions + LIBRETTO-001 cited | Population yield is a projection; caption should explicitly state assumed TCGA RET-fusion rate used |
| 6B | n=504 implied | Direction stated (DM1 enriched intermediate) | N/A (mosaic) | Not specified | **ATA tier is stage-proxy only** (implemented as stage-based in `_build_master.py`); caption does not state this approximation — must add caveat "ATA 2015 risk tier approximated from TNM stage; full ATA criteria require lymph-node and BRAF data not uniformly available in TCGA clinical tables" |
| 6C | Trial n from NCT papers | Reference cited | N/A | NCT00085293/01065090 cited | N and response rates must be verified from original NCT publications before locking |
| 6D | n for GSE151179 not given | d=-1.01; p=0.0001 | N/A (boxplot) | GSE151179 cited | **n for pre vs post-RAI arms not given** in caption |
| 6E | n=504+632 | "dominates over threshold window" | N/A (DCA) | Not specified | "Clinically relevant threshold window" not defined numerically |
| 6F | N/A | N/A (illustrative) | N/A | Stated not pre-registered | Caption is adequate given illustrative framing |
| 6G | n=504 | Distribution direction stated | N/A (waterfall) | Not specified | Source data file not named |

### Extended Data caption gaps

All ED captions are brief single-sentence descriptions (not full structured captions). Before submission:
- Each ED figure needs a structured caption with n, cohort, test, error bars, and source data.
- ED4 depends on GSE241184 and GSE232237, which are not confirmed as present in the project.
- ED11 (Mun 2025 proteome) — the render exists but the source data path in `build_ed.py` must be verified against `project/manuscript_biorxiv_2026_05_20/supp_data_bundle/04_Mun_2025_protein_unsup.tsv`.
- ED15 (Landa 2016 reverse-causality) — GSE76039 processed matrix path unconfirmed.

---

## Numbering conflicts between architectures

| v2 draft text reference | v2 6-figure caption mapping | 5-figure plan mapping | Resolution |
|---|---|---|---|
| "Fig. 1a" — panel selection | Fig. 1A study schematic | Fig. 1A schematic | Consistent |
| "Fig. 1e" — ARI ladder | Fig. 1F ARI ladder | Fig. 1D ARI ladder | Panel letter conflict; v2 text and caption disagree (text says Fig. 1e, caption calls it Fig. 1F) |
| "Fig. 1f" — driver neutrality | Fig. 1E driver neutrality | Fig. 1C driver neutrality | Same issue |
| "Fig. 5a" — OS forest | Fig. 5A forest | Fig. 5A forest | Consistent |
| "Fig. 7a" — head-to-head comparison | Not in 6-figure captions | Not in 5-figure plan | **Conflict: Fig. 7 appears in v2 results text but has no corresponding caption in 05_figure_captions_NC.md** |
| "Fig. 7b" — BRAF interaction | Not in 6-figure captions | Not in 5-figure plan | **Same conflict** |
| "Fig. 7c" — multi-modal concordance | Not in 6-figure captions | Not in 5-figure plan | **Same conflict** |
| "Fig. 7d" — modality comparison | Not in 6-figure captions | Not in 5-figure plan | **Same conflict** |
| "Fig. 8a–8d" — IHC + power | Not in 6-figure captions | Not in 5-figure plan | **Same conflict** |
| "Fig. 6d" — post-RAI box | Fig. 6D post-RAI box | Fig. 5C post-RAI box | Letter/number consistent within 6-fig plan; 5-fig plan maps it to Fig 5 |

**Critical finding:** The v2 results text (`NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md`) references Fig. 7 and Fig. 8 explicitly in the manuscript body. These do not have corresponding captions in `05_figure_captions_NC.md`. The content of these figures (BRAF+ PFI interaction p=0.022; IHC 3-plex power simulation) is described in the results prose. **Before submission, the author must either: (a) add Fig. 7 and Fig. 8 back as main figures with full captions and rendered assets, (b) move the corresponding results paragraphs to a supplement or results section without a figure citation, or (c) map the content into existing figures (e.g., BRAF+ interaction into Fig. 5 as panel B-extension; IHC into ED6).**

This is currently the single most critical unresolved structural conflict in the manuscript.

---

## Summary of required actions before submission

### Computational (must generate)

1. **Fig. 3C** — Real beta × RNA expression scatter for TPO/DIO1/TSHR/TG. Data exists. Script fix is ~30 lines. Join `r5_2_sample_methylation_8gene.tsv` with `panel_expression_thpa_tcga_gdc.tsv` by 12-char barcode. Compute and lock Spearman ρ per gene.
2. **Fig. 3D (panel F)** — Per-driver-class mean beta. Data exists. Join `r5_2_sample_methylation_8gene.tsv` with `master_tcga.tsv`. Compute groupby mean and confirm match to caption values. Trivial.
3. **Fig. 5C** — Multi-cohort KM stack. Needs MSK per-patient OS table. Locate or confirm it is unavailable (in which case this panel must be replaced or removed).
4. **Fig. 5F** — Time-dependent ROC. Data exists. Write `fig5f_time_dependent_roc.py`.
5. **Fig. 5G** — Calibration. **Confirm with author whether 16 OS events is sufficient for a 10-decile calibration plot; may need to be dropped or replaced.**
6. **Fig. 6A** — Reflex flowchart schematic. Write `fig6a_reflex_flowchart.py`.
7. **Fig. 6D** — Post-RAI boxplot. Validate `EV_rai_lineage_box.png` provenance; rerender if needed.
8. **Fig. 6E** — DCA. Write `fig6e_decision_curve.py`. Confirm threshold range.
9. **Fig. 6G** — Selpercatinib waterfall. Write `fig6g_selpercatinib_waterfall.py`.
10. **ED11** — Confirm `04_Mun_2025_protein_unsup.tsv` is the correct source; audit `build_ed.py`.
11. **ED15 / Fig. 3G** — Trace GSE76039 processed matrix path from `dm1_robustness_v2026_05_08/round4/` scripts.
12. **ED4, Fig. 4C, Fig. 4D** — Download GSE241184 and GSE232237 from GEO if these panels are to be retained.

### Editorial decisions required

1. **Resolve Fig. 7 / Fig. 8 conflict** — decide fate of BRAF+ PFI interaction claim and IHC 3-plex content.
2. **Fig. 5G calibration** — confirm or drop given 16 OS events.
3. **Fig. 6B ATA tier proxy** — add caveat that stage-based proxy is used, not full ATA 2015 criteria, or acquire proper ATA staging from TCGA clinical tables.
4. **Fig. 6F trial schema** — confirm illustrative framing is acceptable to target journal.
5. **"Extended Data" terminology** — Nature Communications as of 2024 does not use "Extended Data" as a separate section for original articles unless the journal template specifies it. All ED figures may need to be reclassified as Supplementary Figures. Confirm with NC author guidelines before finalizing figure labels.

### Caption gaps requiring author attention

- Fig. 1E: Add p-values for all driver AUCs, not just BRAF.
- Fig. 2F: Add n for FVPTC and non-FVPTC.
- Fig. 3C: **Must add real Spearman ρ per gene after script fix.**
- Fig. 3F: Add CIs or SE for per-driver mean β.
- Fig. 4C/4D: Add n for GSE241184 and GSE232237.
- Fig. 6B: Add caveat on stage-proxy ATA approximation.
- Fig. 6D: Add n for pre and post-RAI arms in GSE151179.
- All ED captions: Draft full structured captions (n, cohort, test, error bars, source data).

---

*This report is read-only audit output. No figures were modified or generated. No manuscript prose was finalized.*
