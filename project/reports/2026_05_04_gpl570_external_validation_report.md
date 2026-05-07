# GPL570 external validation pack — first-pass report
**Date:** 2026-05-04
**Scope:** Three GPL570 thyroid cancer microarray cohorts (GSE33630, GSE29265, GSE65144) used to extend the Landa 2016 / GSE76039 advanced-disease replication of the driver-orthogonal transcriptional differentiation axis (RAI lineage silencing readout).
**Discipline:** Processed Series Matrix only, within-cohort z-score per dataset, no TCGA classifier transfer, no batch-corrected pooling, no survival/mutation modelling, CPU-only.

---

## 1. Access summary

All three datasets passed access check (`2026_05_04_gpl570_validation_access_check.md`):

| Dataset | n | Histology | Platform | Series Matrix |
|---|---|---|---|---|
| GSE33630 | 105 | 49 PTC + 11 ATC + 45 Normal | GPL570 | ~30 MB |
| GSE29265 | 49 | 20 PTC + 9 ATC + 20 Normal | GPL570 | ~15 MB |
| GSE65144 | 25 | 12 ATC + 13 Normal | GPL570 | ~5.7 MB |
| **Total** | **179 new + 37 anchor (GSE76039) = 216 GPL570 thyroid samples** | | | |

GPL570 annotation (54,675 probes) reused from the Landa-2016 SOFT.

---

## 2. Sample / dataset summary

- **PTC** total across new cohorts: 69 (49 GSE33630 + 20 GSE29265). **+ 0 PTC** in the GSE76039 anchor.
- **ATC** total across new cohorts: 32 (11 + 9 + 12). **+ 20 ATC** in GSE76039 anchor.
- **Normal** total: 78 (45 + 20 + 13). **+ 0 Normal** in GSE76039 anchor.
- **PDTC**: 0 in the validation pack (GSE76039 remains the only PDTC source).

GSE29265 origin sub-stratification (sporadic vs Chernobyl PTC) is recorded in `sample_metadata.tsv` but **not used** as a contrast in this first pass.

---

## 3. Gene panel availability

22/22 panel genes mapped on GPL570 in all three datasets via best-mean-probe collapse:

```
RAI_8 (8):              SLC5A5, TPO, TG, TSHR, PAX8, NKX2-1 (alias TITF1), FOXE1, DIO1
THYROID_NONOVERLAP (8): SLC26A4, IYD, DUOX1, DUOX2, TFF3, HHEX, GLIS3, DIO2
TF_collapse (4):        FOXE1, NKX2-1, PAX8, HHEX
STAT3_AP1_DNMT (5):     STAT3, FOSL1, JUNB, DNMT1, DNMT3B
TROP2 (1):              TACSTD2
Panel union (deduped):  22 genes
```
Probe-level mapping persisted to `probe_to_gene_panel.tsv` (one row per gene per dataset).

---

## 4. Per-dataset results

Within-cohort z-score per gene → arithmetic mean per module → Mann–Whitney U (two-sided) and Cohen's d for between-group contrasts; Spearman across all samples for module-vs-module relationships.

### 4.1 GSE33630 (n=105 — full ladder)

**ATC vs Normal (n=11 vs 45)**

| axis | Cohen d | MW p |
|---|---|---|
| RAI_8 | −5.430 | 3.5e−7 |
| DM1_like | +5.430 | 3.5e−7 |
| THYROID_NONOVERLAP | **−6.532** | 3.5e−7 |
| TDS_like (16-gene) | −6.099 | 3.5e−7 |
| TF_collapse | −4.619 | 3.9e−7 |
| STAT3_AP1_DNMT | +1.707 | 1.9e−4 |
| TACSTD2 (TROP2) z | −0.039 | 0.73 (n.s.) |

**ATC vs PTC (n=11 vs 49)**

| axis | Cohen d | MW p |
|---|---|---|
| RAI_8 | −3.853 | 8.3e−7 |
| DM1_like | +3.853 | 8.3e−7 |
| THYROID_NONOVERLAP | −4.121 | 5.0e−7 |
| TDS_like | −4.122 | 6.2e−7 |
| TF_collapse | −4.082 | 3.7e−7 |
| STAT3_AP1_DNMT | +1.128 | 5.9e−3 |
| TACSTD2 z | −1.721 | 8.3e−4 |

**Spearman across all 105 samples**

- DM1_like vs THYROID_NONOVERLAP: ρ = **−0.941**, p = 2.8e−50
- TF_collapse vs DM1_like: ρ = −0.891, p = 4.2e−37
- STAT3_AP1_DNMT vs DM1_like: ρ = +0.428, p = 5.4e−6
- TACSTD2 vs DM1_like: ρ = +0.439, p = 2.9e−6

(PTC vs Normal contrasts also computed and persisted to `gpl570_score_tests.tsv`; not headlined here.)

### 4.2 GSE29265 (n=49 — full ladder, smaller)

**ATC vs Normal (n=9 vs 20)**

| axis | Cohen d | MW p |
|---|---|---|
| RAI_8 | −2.205 | 8.9e−4 |
| DM1_like | +2.205 | 8.9e−4 |
| THYROID_NONOVERLAP | −2.507 | 6.8e−5 |
| TDS_like | −2.381 | 8.3e−5 |
| TF_collapse | −2.480 | 1.0e−4 |
| STAT3_AP1_DNMT | +1.264 | 1.2e−2 |
| TACSTD2 z | −0.807 | 2.0e−2 |

**ATC vs PTC (n=9 vs 20)**

| axis | Cohen d | MW p |
|---|---|---|
| RAI_8 | −1.505 | 3.6e−2 |
| DM1_like | +1.505 | 3.6e−2 |
| THYROID_NONOVERLAP | −1.483 | 1.5e−2 |
| TDS_like | −1.526 | 2.5e−2 |
| TF_collapse | −2.106 | 3.7e−3 |
| STAT3_AP1_DNMT | +1.180 | 1.7e−2 |
| TACSTD2 z | −1.981 | 7.5e−4 |

**Spearman across all 49 samples**

- DM1_like vs THYROID_NONOVERLAP: ρ = **−0.846**, p = 2.1e−14
- TF_collapse vs DM1_like: ρ = −0.877, p = 1.3e−16
- STAT3_AP1_DNMT vs DM1_like: ρ = +0.306, p = 3.2e−2
- TACSTD2 vs DM1_like: ρ = +0.246, p = 8.9e−2 (n.s.)

### 4.3 GSE65144 (n=25 — ATC vs Normal only)

**ATC vs Normal (n=12 vs 13)**

| axis | Cohen d | MW p |
|---|---|---|
| RAI_8 | −2.946 | 3.2e−5 |
| DM1_like | +2.946 | 3.2e−5 |
| THYROID_NONOVERLAP | −3.227 | 3.2e−5 |
| TDS_like | −3.226 | 2.5e−5 |
| TF_collapse | −3.132 | 2.5e−5 |
| STAT3_AP1_DNMT | +1.103 | 1.1e−2 |
| TACSTD2 z | +0.253 | 0.61 (n.s.) |

**Spearman across all 25 samples**

- DM1_like vs THYROID_NONOVERLAP: ρ = **−0.904**, p = 6.0e−10
- TF_collapse vs DM1_like: ρ = −0.966, p = 5.0e−15
- STAT3_AP1_DNMT vs DM1_like: ρ = +0.589, p = 1.9e−3
- TACSTD2 vs DM1_like: ρ = −0.015, p = 0.94 (n.s.)

---

## 5. Direction-consistency summary

Headline ATC contrasts (ATC vs PTC and ATC vs Normal where available): 35 axis × dataset × contrast cells.

**Lineage axes (RAI_8, DM1_like, THYROID_NONOVERLAP, TDS_like, TF_collapse) AND inflammation/methylation axis (STAT3_AP1_DNMT):**

| axis | expected sign | observed direction-consistent cells | median d (signed) |
|---|---|---|---|
| RAI_8_score | − | **5/5** | −2.95 |
| DM1_like_score | + | **5/5** | +2.95 |
| THYROID_NONOVERLAP_score | − | **5/5** | −3.23 |
| TDS_like_score | − | **5/5** | −3.23 |
| TF_collapse_score | − | **5/5** | −3.13 |
| STAT3_AP1_DNMT_score | + | **5/5** | +1.18 |
| **Subtotal lineage + axis** | | **30/30** | |
| TACSTD2 (TROP2) z | + | **1/5** | −0.81 |

**Headline:** 30/30 (100%) direction-consistent for the differentiation/inflammation axis. **TACSTD2 (TROP2) at GPL570 microarray bulk level does NOT replicate the ATC-elevated direction seen in GSE76039 advanced-disease comparison** — only GSE65144 (ATC vs Normal) shows the expected positive direction (d = +0.25, n.s.). This is a meaningful negative on the TROP2 sub-claim and is not concealed.

**Module-vs-module orthogonality (Spearman ρ for DM1_like vs THYROID_NONOVERLAP):**

| Cohort | n | ρ | p |
|---|---|---|---|
| GSE76039 (anchor, ATC+PDTC) | 37 | −0.925 | 3.1e−16 |
| GSE33630 | 105 | −0.941 | 2.8e−50 |
| GSE29265 | 49 | −0.846 | 2.1e−14 |
| GSE65144 | 25 | −0.904 | 6.0e−10 |
| **All four GPL570 cohorts** | **216** | **strong negative across the board** | |

---

## 6. Comparison with GSE76039 anchor

- **Differentiation axis:** GSE76039 ATC vs PDTC d ≈ −3.47 (RAI_8) is now flanked by GSE33630 ATC vs PTC d ≈ −3.85 and ATC vs Normal d ≈ −5.43, GSE29265 ATC vs PTC d ≈ −1.51 and ATC vs Normal d ≈ −2.20, and GSE65144 ATC vs Normal d ≈ −2.95.
- **THYROID_NONOVERLAP module** (zero overlap with the RAI_8 readout) replicates with d = −1.48 to −6.53 across all five contrasts where both groups are present, confirming the anchor's most surprising finding (orthogonal-module replication).
- **Spearman ρ** for the within-cohort orthogonality is essentially identical across cohorts (−0.85 to −0.94, anchor −0.93).
- **STAT3 / AP-1 / DNMT axis** (mechanism-side) replicates direction in all 5 ATC contrasts (anchor d ≈ +1.72, validation pack +1.10 to +1.71).
- **TACSTD2 (TROP2)** does not replicate at GPL570 bulk level (anchor d ≈ +1.13 in ATC vs PDTC; validation pack 4/5 wrong direction). This is consistent with the known caveat that TROP2 elevation in advanced disease may be enriched at the spatial / clonal level and is not necessarily resolvable at bulk-microarray scale.

---

## 7. Verdict

| Claim | Verdict |
|---|---|
| Driver-orthogonal differentiation axis (RAI_8 / DM1_like) replicates in independent GPL570 thyroid cancer cohorts | **MAIN-FIGURE EXTENSION** — strong 5/5 direction-consistency, monotonic-magnitude in the ATC-vs-Normal contrast where Normal is available |
| THYROID_NONOVERLAP zero-overlap module replicates the same axis | **MAIN-FIGURE EXTENSION** — 5/5 direction-consistent; ρ ≈ −0.85 to −0.94 with DM1_like across cohorts |
| Inflammation / AP-1 / DNMT axis (STAT3_AP1_DNMT) elevates with ATC | **SUPPLEMENT** — 5/5 direction-consistent, smaller magnitude (median d ≈ +1.18); supportive but not primary |
| TROP2 (TACSTD2) bulk mRNA elevation in advanced disease | **INTERNAL ONLY** at GPL570 bulk level — does NOT replicate (only 1/5 expected direction). Do not headline as "TROP2 elevation in advanced disease" from this microarray pack. Spatial / clonal evidence remains a separate question. |

---

## 8. Risks

- **Microarray platform.** GPL570 is bulk Affymetrix; per-gene quantitation differs from RNA-seq, and our cross-platform calibration audit (`v17_korean_K2_calibration`) prohibits TCGA-trained absolute-form classifier transfer. We therefore use within-cohort z only — this is fit for direction and effect-size replication, **not** for shared-cutoff prediction.
- **Histology-label quality.** Labels are taken at face value from GEO `Sample_characteristics_ch1` / `Sample_title`. No central pathology re-review. PTC vs ATC is unambiguous in the labels but no inferred sub-classification is added beyond what the GEO record states.
- **Probe mapping.** Best-mean-probe collapse per gene; identical methodology as GSE76039. NKX2-1 mapped via TITF1 alias (Affy legacy). 22/22 panel genes found in all three datasets; this is the simpler probe-mapping case.
- **No survival / mutation / age / stage** in any of these GEO records. No correlative claims of that kind are made.
- **Sample sizes.** GSE65144 lacks PTC; GSE29265 ATC arm is small (n=9); GSE33630 ATC arm is moderate (n=11). Absolute-magnitude inference per dataset must respect the per-arm n.

---

## 9. Safe wording (use)

- "External GPL570 expression validation across 3 independent thyroid cancer cohorts."
- "Direction-consistent lineage-silencing readout across independent thyroid cancer datasets."
- "Within-cohort z-score per dataset; no cross-platform classifier transfer."
- "Replication of the orthogonality between the RAI_8 readout and a zero-overlap THYROID_NONOVERLAP module."

## 10. Forbidden wording (do NOT use)

- "Clinical validation"
- "Progression proven" / "PTC → ATC progression demonstrated"
- "Survival validated"
- "Fusion validated"
- "Monotonic 3-stage gradient" (Normal → PTC → ATC magnitudes are *direction-consistent*, not asserted as monotonic ordinal)
- "TROP2 elevation in advanced disease confirmed" (TACSTD2 negative at GPL570 bulk level — see §5/§7)

---

## 11. Files produced

- `project/results/p_gpl570_validation/sample_metadata.tsv` — 179 rows, dataset / sample_id / histology_raw / histology_clean / disease_group / platform / source_file / notes.
- `project/results/p_gpl570_validation/probe_to_gene_panel.tsv` — 66 rows (22 genes × 3 datasets), best-probe per gene.
- `project/results/p_gpl570_validation/{ACC}_expression_gene_log.tsv.gz` — panel-gene log expression matrix per dataset (3 files).
- `project/results/p_gpl570_validation/{ACC}_scores.tsv` — module scores per sample per dataset (3 files).
- `project/results/p_gpl570_validation/gpl570_score_tests.tsv` — all between-group MW + Cohen's d + Spearman tests, long format.
- `project/results/p_gpl570_validation/gpl570_meta_effect_summary.tsv` — coverage table, per-dataset axis summary, overall-per-axis summary, full headline table.
- `project/results/p_gpl570_validation/gpl570_rai_lineage_boxplots.png` — 7 axes × 3 datasets boxplot grid (Normal/PTC/ATC where available).
- `project/results/p_gpl570_validation/gpl570_dm1_nonoverlap_scatter_grid.png` — DM1_like vs THYROID_NONOVERLAP scatter, one panel per dataset, with Spearman ρ in title.
- `project/results/p_gpl570_validation/gpl570_direction_consistency_forest.png` — Cohen's d forest, color-coded for direction consistency vs expected sign.
- `project/notebooks_or_scripts/p_gpl570_validation_first_pass.py` — script that produces all of the above.

---

**End of report. No prose was written for any voice-protected manuscript section in this analysis.**
