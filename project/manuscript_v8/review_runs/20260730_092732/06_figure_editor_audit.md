---
title: Figure editor audit
audit_date: 2026-07-30
auditor: manuscript-audit-agent (figure editor mode — no figure generation)
---

# 06 — Figure editor audit

## 1. Figure/panel references found in v2 draft text

All figure references found in `NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md`:

| Reference in text | Location (approx. line) |
|---|---|
| Fig. 1a,b | line 58 (panel gene curation) |
| Fig. 1c,d | line 60 (UMAP and clustering) |
| Fig. 1e | line 62 (ARI ladder) |
| Fig. 1f | line 62 (driver mRNA neutrality) |
| Fig. 2a | line 68 (dark matter Sankey) |
| Fig. 2b | line 70 (DM1 prevalence by cohort) |
| Fig. 3 | line 98 (differentiation signal reference) |
| Fig. 3a | line 246 caption (HM450 heatmap) |
| Fig. 3b | line 246 caption (mean beta) |
| Fig. 3c | line 246 caption (methylation-expression scatter) |
| Fig. 3d | line 246 caption (MAPK × Panel forest) |
| Fig. 4a,b | line 90 (Lu 2023 UMAP; Pu 2021 correlations) |
| Fig. 5a | line 72 (OS meta-analysis forest) |
| Fig. 5e | line 92 (FFPE vs FF concordance) |
| Fig. 6d | line 124 (post-RAI refractory alignment) |
| Fig. 7a | line 62 (head-to-head comparison) |
| Fig. 7b | lines 76, 78 (interaction result; subgroup HRs) |
| Fig. 7b, panel A | line 76 |
| Fig. 7b, panel B | line 78 |
| Fig. 7c | line 88 (external replication forest) |
| Fig. 7d | line 108 (multi-omics 5-way convergence) |
| Fig. 8a | line 98 (14 combination comparison) |
| Fig. 8b,c | line 100 (IHC 3-plex KM + comparison) |
| Fig. 8d | line 115 (Monte Carlo power simulation) |

**Total main figure references: Fig 1–8, with multiple sub-panels.**

---

## 2. Figure assets that actually exist

### In `project/results/manuscript_v8_nc_main/` (confirmed):

| Asset file | Corresponds to |
|---|---|
| `Fig1_discovery_axis.pdf/png` | Main Fig 1 |
| `Fig2_fusion_mechanism.pdf/png` | Main Fig 2 |
| `Fig3_epigenetic.pdf/png` | Main Fig 3 (partial — see missing panels below) |
| `Fig4_sc_validation.pdf/png` | Main Fig 4 |
| `Fig5_survival_portability.pdf/png` | Main Fig 5 |
| `Fig6_reflex_translation.pdf/png` | Main Fig 6 |
| `F3_driver_neutrality.pdf/png` (in `results/figures/`) | Fig 1 driver-neutrality panel |
| `F4_pangenome_robustness.pdf/png` (in `results/figures/`) | Fig 1 ARI ladder panel |

### In `project/dm1_story_web/public/figures/` and `dist/figures/` (confirmed):

| Asset file | Corresponds to |
|---|---|
| `fig_deep1_head_to_head.png` | v2 Fig 7a (head-to-head comparison) |
| `fig_deep2_prognostic_predictive.png` | v2 Fig 7b (interaction + subgroup HRs) |
| `fig_a2_replication.png` | v2 Fig 7c (external replication forest) |
| `fig_deep4_multiomics.png` | v2 Fig 7d (multi-omics convergence) |
| `fig_ihc3_killer.png` | v2 Fig 8 (IHC 3-plex + power simulation) |
| `fig_deep5_snubh_simulation.png` | v2 Fig 8d (Monte Carlo power simulation) |
| `fig_deep3_redifferentiation.png` | v2 Fig 6c (redifferentiation schematic?) |

### In `manuscript_v8/figures/` (only 2 files):

| Asset file | Corresponds to |
|---|---|
| `Fig7_dm1_mechanism.pdf/png` | Legacy Fig 7 (mechanism figure from v1 era) |
| `Fig8_epigenetic.pdf/png` | Legacy Fig 8 (epigenetic figure from v1 era) |

---

## 3. Orphaned references (cited in text but no confirmed asset)

| Citation | Asset status |
|---|---|
| Fig. 3c (methylation vs expression scatter, 4-gene grid: TPO/DIO1/TSHR/TG) | **MISSING** — confirmed not built. STATUS_INSILICO_FINALIZATION_2026_07_30.md lists this as remaining work. Script `fig3c_beta_expression_scatter.py` noted as needed. |
| Fig. 3d (per-driver-class beta bar: BRAF/RET/RAS) | **NEEDS REBUILD** — data exists (values in AUDIT_LOCKED_RESULTS.md M-12) but visualization asset not confirmed in nc_main directory. |
| Fig. 5d (ATA intermediate-risk mosaic) | **NEEDS BUILD** — STATUS file lists "Fig 5 ATA intermediate-risk mosaic rebuild" as remaining work. |
| Fig. 2b (DM1 prevalence by TCGA vs Korean cohorts) | Asset unclear — may be in Fig2_fusion_mechanism.png but panel letter b is specified as DM1 sub-A/sub-B silhouette in NC captions (05_figure_captions_NC.md line 39) not prevalence bar. Panel letter assignment conflict. |
| "Fig. 7c" (External replication forest) | `fig_a2_replication.png` exists but this is labeled as A2 replication; the full Lee+K2+GPL570+Landa forest referenced in v2 Fig 7c caption may be the EV_master_forest_annotated.png instead. Clarification needed. |

---

## 4. "Rebuild" and "missing" panels from REVIEWER_FIRST_NC_REBUILD_2026_07_08.md

The reviewer-first 5-figure plan explicitly lists:

| Rebuild/Build item | Status per current file inventory |
|---|---|
| Fig 1: "Rebuild simple schematic (1A)" | No schematic confirmed in nc_main |
| Fig 1: "Use F3_driver_neutrality.png (1C)" | PRESENT |
| Fig 1: "Use F4_pangenome_robustness.png, simplify to four bars (1D)" | PRESENT |
| Fig 2: "Rebuild stacked bar for DM1 prevalence by driver class (1E)" | Not confirmed in nc_main |
| Fig 3: "Rebuild 3C beta-expression scatter (missing mechanism bridge)" | **CONFIRMED MISSING** |
| Fig 3: "Rebuild 3D per-driver-class beta from v5A/per-driver class table" | Not confirmed in nc_main |
| Fig 4: "Use GSE76039 lineage plot/heatmap (2D)" | Not confirmed as Landa heatmap in nc_main |
| Fig 5: "Rebuild 5D GSE151179 post-RAI box" | `fig6d_post_rai_box.py` noted in captions as needed; asset not found |

---

## 5. Figure numbering conflict between v1 and v2 drafts

**Critical conflict:**

| Figure | v1 (reviewer-first 5-fig plan) | v2 draft (8-figure) | NC captions (6-figure) |
|---|---|---|---|
| Discovery (ARI, driver neutrality) | Fig 1 | Fig 1 | Fig 1 |
| Biological meaning / fusion | Fig 2 | Fig 2 | Fig 2 |
| Methylation mechanism | Fig 3 | Fig 3 | Fig 3 |
| External validation | Fig 4 | Fig 4 | Fig 4 |
| Clinical relevance / OS | Fig 5 | Fig 5 | Fig 5 |
| Reflex translation / GSE151179 | — (not in 5-fig) | Fig 6 | Fig 6 |
| Interaction + external replication | — (not in 5-fig) | **Fig 7 (NEW)** | — |
| IHC 3-plex + power simulation | — (not in 5-fig) | **Fig 8 (NEW)** | — |

**Specific conflict:** The v2 draft refers to "Fig. 7" (interaction result) and "Fig. 8" (IHC 3-plex) extensively. The NC captions file (05_figure_captions_NC.md) does NOT include Fig 7 or Fig 8 — it ends at Fig 6. The reviewer-first rebuild plan describes a 5-figure architecture with no Fig 7 or 8. The legacy figures in `manuscript_v8/figures/` are also named Fig7 and Fig8 but represent different content (mechanism and epigenetic schematic from v1 era).

**This conflict means that Fig 7 and Fig 8 in the v2 draft refer to figures that:**
- Are not in the NC captions file
- Are not in the reviewer-first rebuild plan
- Have the same names as legacy v1 figures with DIFFERENT content
- Exist only in `dm1_story_web/public/figures/` and `dist/figures/` rather than in the main `manuscript_v8_nc_main/` figure directory

**Resolution required:** An explicit editorial decision must determine whether the final submission will be a 5-figure, 6-figure, or 8-figure manuscript. Once decided, the figure numbering must be consistent across the draft text, all captions, and all figure files.

---

## 6. Extended Data terminology

The reviewer-first plan uses "Extended Data" terminology (e.g., "Extended Data Set" with 6 ED figures, numbered ED1-ED6). The NC captions file (05_figure_captions_NC.md) uses "Extended Data Fig. 1-15." The `manuscript_v8_nc_main/` directory uses "ED1-ED15" naming.

**NC journal policy:** Nature Communications does NOT use "Extended Data" — that term is specific to Nature, Nature Medicine, Nature Cancer, and other primary Nature journals. Nature Communications uses "Supplementary Figures" or "Supplementary Information." All "Extended Data" references in the manuscript must be changed to "Supplementary Figure" before submission.

**Action:** Replace all "Extended Data" and "ED" with "Supplementary Figure" throughout ALL manuscript files (v2 draft, captions, reviewer-first plan). Rename `ED1-ED15_*.png` files to `SFig1-SFig15_*.png`.

DM1_PROJECT_STATE.md (line 88) explicitly flags: "Extended Data terminology should be checked against current Nature Communications requirements and replaced with Supplementary Figures if needed." — This has NOT been acted on.

---

## 7. Summary of figure completeness

| Figure | Asset status | Issue |
|---|---|---|
| Fig 1 (a-g) | Partial — core panels present, schematic strip and pooled KM teaser may need build | MINOR |
| Fig 2 (a-g) | Largely present in nc_main | MINOR — panel letter assignment vs NC captions needs reconciliation |
| Fig 3 (a-g) | **Panel 3C MISSING; Panel 3D needs rebuild** | MAJOR |
| Fig 4 (a-g) | Present in nc_main | MINOR |
| Fig 5 (a-g) | Partial — time-dep ROC and calibration noted as needed; ATA mosaic needed | MAJOR |
| Fig 6 (a-g) | Partial — post-RAI box and reflex flowchart need build | MAJOR |
| Fig 7 (a-d) | Present in dm1_story_web but not in nc_main; not in NC captions | UNRESOLVED — architectural decision needed |
| Fig 8 (a-f) | Present in dm1_story_web but not in nc_main; not in NC captions | UNRESOLVED — architectural decision needed |
| ED / Supp Figs | 15 PNG files in nc_main labeled ED1-ED15; must be renamed to SFig | TERMINOLOGY ERROR |
