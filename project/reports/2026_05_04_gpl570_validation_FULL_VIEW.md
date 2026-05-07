# GPL570 thyroid dedifferentiation validation pack — FULL VIEW
**Date:** 2026-05-04
**Anchor:** Landa 2016 / GSE76039 first-pass (PDTC + ATC, n=37)
**Validation pack:** GSE33630 (n=105) + GSE29265 (n=49) + GSE65144 (n=25) → +179 GPL570 thyroid samples
**Discipline:** CPU-only · within-cohort z-score per dataset · processed Series Matrix only (no CEL) · no TCGA-trained classifier transfer · no batch-corrected pooling · no survival/mutation/age/stage modelling
**Voice protection:** No manuscript Hook / Aim / Discussion / Limitations / Cover / Q9 prose drafted.

---

## 0. TL;DR

- **30/30** lineage-axis × dataset × ATC-contrast cells direction-consistent (RAI_8, DM1_like, THYROID_NONOVERLAP, TDS_like, TF_collapse, STAT3_AP1_DNMT).
- **Within-cohort orthogonality** (DM1_like ⊥ THYROID_NONOVERLAP) Spearman ρ ∈ [−0.846, −0.941] across the 3 validation cohorts (anchor −0.925).
- **TROP2 / TACSTD2 at GPL570 bulk level does NOT replicate** the ATC-elevated direction: 1/5 ATC contrasts in expected direction. Honest negative — supplement / spatial follow-up only.
- **22/22** panel genes mappable on GPL570 in all 3 datasets (best-mean-probe collapse).
- **Largest lineage effects:** GSE33630 ATC vs Normal — THYROID_NONOVERLAP d = −6.53; TDS_like d = −6.10; RAI_8 d = −5.43; DM1_like d = +5.43.
- **Verdict:** main-figure extension for the differentiation axis claim; supplement for STAT3_AP1_DNMT; internal-only for TROP2 bulk-microarray claim.

---

## 1. Access summary

All three datasets fully accessible via processed GEO Series Matrix; raw CEL not used.

| Dataset | Title (short) | Platform | n | Composition | Series Matrix | Histology source |
|---|---|---|---|---|---|---|
| **GSE33630** | Tomás 2011 — N / PTC / ATC | GPL570 | 105 | 49 PTC + 11 ATC + 45 Normal | ~30 MB | `pathological diagnostic` |
| **GSE29265** | Tomás 2012 — sporadic vs Chernobyl PTC + ATC | GPL570 | 49 | 20 PTC + 9 ATC + 20 Normal | ~15 MB | `Sample_title` |
| **GSE65144** | von Roemeling 2015 — ATC vs Normal | GPL570 | 25 | 12 ATC + 13 Normal | ~5.7 MB | `tissue type` |
| GSE76039 *(anchor)* | Landa 2016 — PDTC + ATC | GPL570 | 37 | 17 PDTC + 20 ATC | already on disk | `source_name` |

GPL570 SOFT annotation (54,675 probes) reused from `project/results/p_landa_2016/raw/GPL570_full.soft`. Raw Series Matrices stored under `project/results/p_gpl570_validation/raw/<ACC>/` and gitignored.

---

## 2. Cohort composition (full)

| Dataset | Normal | PTC | ATC | PDTC | Total |
|---|---|---|---|---|---|
| GSE33630 | 45 | 49 | 11 | 0 | 105 |
| GSE29265 | 20 | 20 | 9 | 0 | 49 |
| GSE65144 | 13 | 0 | 12 | 0 | 25 |
| GSE76039 | 0 | 0 | 20 | 17 | 37 |
| **Total GPL570 thyroid samples** | **78** | **69** | **52** | **17** | **216** |

PDTC absent from validation pack — GSE76039 remains the only PDTC anchor in the GPL570 family.

---

## 3. Gene panel (22 genes, all found on GPL570 in all 3 datasets)

```
RAI_8 (8):              SLC5A5, TPO, TG, TSHR, PAX8, NKX2-1 (alias TITF1), FOXE1, DIO1
THYROID_NONOVERLAP (8): SLC26A4, IYD, DUOX1, DUOX2, TFF3, HHEX, GLIS3, DIO2
TDS_like (16):          RAI_8 ∪ THYROID_NONOVERLAP
TF_collapse (4):        FOXE1, NKX2-1, PAX8, HHEX
STAT3_AP1_DNMT (5):     STAT3, FOSL1, JUNB, DNMT1, DNMT3B
TROP2 (1):              TACSTD2
Panel union (deduped): 22 genes
```
Best-mean-probe collapse per gene per dataset; identical methodology to GSE76039 anchor.

---

## 4. Headline ATC contrasts — Cohen's d

Two-sided Mann–Whitney U; within-cohort z-scored module/gene.

### 4.1 GSE33630 (n=105)

| Axis | ATC vs Normal (n=11 vs 45) | ATC vs PTC (n=11 vs 49) | PTC vs Normal (n=49 vs 45) |
|---|---|---|---|
| RAI_8_score | **−5.430** (p=3.5e−7) | −3.853 (p=8.3e−7) | −2.373 (p=1.3e−13) |
| DM1_like_score | **+5.430** (p=3.5e−7) | +3.853 (p=8.3e−7) | +2.373 (p=1.3e−13) |
| THYROID_NONOVERLAP | **−6.532** (p=3.5e−7) | −4.121 (p=5.0e−7) | −2.082 (p=6.5e−13) |
| TDS_like (16-gene) | −6.099 (p=3.5e−7) | −4.122 (p=6.2e−7) | −2.316 (p=1.3e−13) |
| TF_collapse | −4.619 (p=3.9e−7) | −4.082 (p=3.7e−7) | −1.537 (p=7.2e−10) |
| STAT3_AP1_DNMT | +1.707 (p=1.9e−4) | +1.128 (p=5.9e−3) | +0.635 (p=6.1e−4) |
| TACSTD2 (TROP2) z | −0.039 (p=0.73 n.s.) | **−1.721** (p=8.3e−4 — wrong dir) | +2.284 (p=2.5e−12, PTC>Normal) |

### 4.2 GSE29265 (n=49)

| Axis | ATC vs Normal (n=9 vs 20) | ATC vs PTC (n=9 vs 20) | PTC vs Normal (n=20 vs 20) |
|---|---|---|---|
| RAI_8_score | −2.205 (p=8.9e−4) | −1.505 (p=3.6e−2) | −1.295 (p=4.6e−4) |
| DM1_like_score | +2.205 (p=8.9e−4) | +1.505 (p=3.6e−2) | +1.295 (p=4.6e−4) |
| THYROID_NONOVERLAP | −2.507 (p=6.8e−5) | −1.483 (p=1.5e−2) | −1.288 (p=4.6e−4) |
| TDS_like | −2.381 (p=8.3e−5) | −1.526 (p=2.5e−2) | −1.388 (p=2.7e−4) |
| TF_collapse | −2.480 (p=1.0e−4) | −2.106 (p=3.7e−3) | −0.877 (p=2.0e−3) |
| STAT3_AP1_DNMT | +1.264 (p=1.2e−2) | +1.180 (p=1.7e−2) | +0.138 (p=0.84 n.s.) |
| TACSTD2 (TROP2) z | **−0.807** (p=2.0e−2 wrong dir) | **−1.981** (p=7.5e−4 wrong dir) | +1.342 (p=4.2e−4, PTC>Normal) |

### 4.3 GSE65144 (n=25, ATC vs Normal only)

| Axis | ATC vs Normal (n=12 vs 13) |
|---|---|
| RAI_8_score | −2.946 (p=3.2e−5) |
| DM1_like_score | +2.946 (p=3.2e−5) |
| THYROID_NONOVERLAP | −3.227 (p=3.2e−5) |
| TDS_like | −3.226 (p=2.5e−5) |
| TF_collapse | −3.132 (p=2.5e−5) |
| STAT3_AP1_DNMT | +1.103 (p=1.1e−2) |
| TACSTD2 (TROP2) z | +0.253 (p=0.61 n.s.) |

### 4.4 GSE76039 anchor (n=37, ATC vs PDTC) — for reference

| Axis | Cohen d | MW p |
|---|---|---|
| RAI_8_score | −3.469 | 4.6e−7 |
| DM1_like_score | +3.469 | 4.6e−7 |
| THYROID_NONOVERLAP | −3.317 | 6.3e−7 |
| TDS_like | −3.538 | 4.6e−7 |
| TF_collapse | −3.566 | 2.4e−7 |
| STAT3_AP1_DNMT | +1.720 | 3.6e−5 |
| TACSTD2 (TROP2) z | +1.131 | 2.7e−3 |

---

## 5. Direction-consistency summary (5 ATC contrasts × 7 axes = 35 cells)

| Axis | Expected sign | n cells consistent / total | Median signed d |
|---|---|---|---|
| RAI_8_score | − | **5/5** | −2.95 |
| DM1_like_score | + | **5/5** | +2.95 |
| THYROID_NONOVERLAP_score | − | **5/5** | −3.23 |
| TDS_like_score | − | **5/5** | −3.23 |
| TF_collapse_score | − | **5/5** | −3.13 |
| STAT3_AP1_DNMT_score | + | **5/5** | +1.18 |
| **Lineage + inflammation/methylation subtotal** | | **30 / 30** | |
| TACSTD2 (TROP2) z | + | **1/5** | −0.81 |

→ **Differentiation axis claim externally replicates in 30/30 cells. TROP2 bulk-microarray sub-claim fails — flagged honestly.**

---

## 6. Module-vs-module orthogonality — Spearman

### 6.1 DM1_like vs THYROID_NONOVERLAP (zero-overlap test)

| Cohort | n | ρ | p |
|---|---|---|---|
| **GSE33630** | 105 | **−0.941** | 2.8e−50 |
| **GSE29265** | 49 | **−0.846** | 2.1e−14 |
| **GSE65144** | 25 | **−0.904** | 6.0e−10 |
| GSE76039 *(anchor)* | 37 | −0.925 | 3.1e−16 |

Independent zero-overlap evidence the same axis is operative in every cohort.

### 6.2 TF_collapse vs DM1_like

| Cohort | n | ρ | p |
|---|---|---|---|
| GSE33630 | 105 | −0.891 | 4.2e−37 |
| GSE29265 | 49 | −0.877 | 1.3e−16 |
| GSE65144 | 25 | −0.966 | 5.0e−15 |
| GSE76039 *(anchor)* | 37 | −0.931 | 7.3e−17 |

### 6.3 STAT3_AP1_DNMT vs DM1_like

| Cohort | n | ρ | p |
|---|---|---|---|
| GSE33630 | 105 | +0.428 | 5.4e−6 |
| GSE29265 | 49 | +0.306 | 3.2e−2 |
| GSE65144 | 25 | +0.589 | 1.9e−3 |
| GSE76039 *(anchor)* | 37 | +0.681 | 3.5e−6 |

### 6.4 TACSTD2 vs DM1_like

| Cohort | n | ρ | p |
|---|---|---|---|
| GSE33630 | 105 | +0.439 | 2.9e−6 |
| GSE29265 | 49 | +0.246 | 8.9e−2 (n.s.) |
| GSE65144 | 25 | −0.015 | 0.94 (n.s.) |
| GSE76039 *(anchor)* | 37 | +0.443 | 6.1e−3 |

(The GSE76039 anchor showed a positive correlation; in the GPL570 validation pack only GSE33630 reproduces it. Consistent with bulk-microarray TROP2 weakness — the spatial / clonal question is separate.)

---

## 7. Comparison with the GSE76039 anchor

| Claim | Anchor (GSE76039) | Validation pack outcome |
|---|---|---|
| RAI_8 silencing in advanced disease | ATC vs PDTC d = −3.47 | flanked by ATC vs Normal d = −2.21 → −5.43 and ATC vs PTC d = −1.51 → −3.85 |
| THYROID_NONOVERLAP module replicates same axis | d = −3.32 (anchor) | d = −1.48 to −6.53 in 5/5 contrasts |
| Module orthogonality (DM1_like ⊥ NONOVERLAP) | ρ = −0.925 | ρ ∈ [−0.846, −0.941] in 3/3 cohorts |
| Inflammation / AP-1 / DNMT axis up in ATC | d = +1.72 | d = +1.10 to +1.71 in 5/5 contrasts |
| TROP2 (TACSTD2) bulk elevation in advanced disease | d = +1.13 (ATC vs PDTC) | 4/5 wrong direction at GPL570 bulk |

---

## 8. Verdict per claim

| Claim | Verdict | Action |
|---|---|---|
| Driver-orthogonal differentiation axis (RAI_8 / DM1_like) replicates externally | **MAIN-FIGURE EXTENSION** | use in main figure / Results |
| Zero-overlap THYROID_NONOVERLAP module replicates the same axis | **MAIN-FIGURE EXTENSION** | use in main figure / Results |
| Inflammation / AP-1 / DNMT axis up in ATC | **SUPPLEMENT** | supportive supp panel; smaller magnitude |
| Within-cohort module orthogonality (DM1_like ⊥ NONOVERLAP) | **REPRODUCED** | supports paragraph on independent-module validation |
| TROP2 (TACSTD2) bulk mRNA elevation in advanced disease | **INTERNAL ONLY at GPL570 bulk level — NEGATIVE** | do not headline; spatial / clonal question separate |

---

## 9. Risks / framing discipline

### 9.1 Platform & calibration
- GPL570 is bulk Affymetrix microarray. Within-cohort z only — no shared cutoff, no TCGA-trained absolute-form transfer (per `v17_korean_K2_calibration`).
- Best-mean-probe collapse per gene; identical methodology to anchor.
- 22/22 panel genes mapped in all three datasets — no missing-gene caveat.

### 9.2 Histology label quality
- Labels taken at face value from GEO `Sample_characteristics_ch1` / `Sample_title`. No central pathology re-review.
- No PDTC anywhere in the validation pack — GSE76039 remains the only PDTC source.
- GSE65144 lacks PTC entirely; GSE29265 ATC arm n=9; GSE33630 ATC arm n=11 — per-arm n must be respected for absolute-magnitude inference.

### 9.3 Metadata that does NOT exist in any of these GEO records
- Survival, mutation status, age (except GSE29265), stage, BRAF/RAS/fusion status. None claimed.

### 9.4 Probe mapping caveats
- NKX2-1 mapped via TITF1 alias (Affy legacy).
- Gene-symbol field may contain ' /// ' joins; we accept any probe whose gene-list contains the requested symbol or alias.

---

## 10. Safe wording (use)

- "External GPL570 expression validation across three independent thyroid cancer cohorts."
- "Direction-consistent lineage-silencing readout across independent cohorts."
- "Independent zero-overlap module replicates the same axis."
- "Within-cohort z-score per dataset; no cross-platform classifier transfer."

## 11. Forbidden wording (do NOT use)

- "Clinical validation"
- "Survival validated"
- "Fusion validated"
- "Progression proven" / "PTC → ATC progression demonstrated"
- "Monotonic 3-stage gradient"
- "TROP2 elevated in advanced disease confirmed" *(TACSTD2 negative at GPL570 bulk)*

---

## 12. Files produced

```
project/results/p_gpl570_validation/
├── sample_metadata.tsv                     179 rows (dataset / sample_id / histology_raw / histology / disease_group / source_file / notes)
├── probe_to_gene_panel.tsv                 66 rows (22 genes × 3 datasets, best-mean-probe collapse)
├── GSE33630_expression_gene_log.tsv.gz     panel-gene log expression, 22 × 105
├── GSE29265_expression_gene_log.tsv.gz     panel-gene log expression, 22 ×  49
├── GSE65144_expression_gene_log.tsv.gz     panel-gene log expression, 22 ×  25
├── GSE33630_scores.tsv                     per-sample module scores, n=105
├── GSE29265_scores.tsv                     per-sample module scores, n= 49
├── GSE65144_scores.tsv                     per-sample module scores, n= 25
├── gpl570_score_tests.tsv                  all MW + Cohen d + Spearman, long format
├── gpl570_meta_effect_summary.tsv          coverage + per-dataset axis + overall + headline tables
├── gpl570_rai_lineage_boxplots.png         Figure 1 — 7 axes × 3 datasets boxplot grid
├── gpl570_dm1_nonoverlap_scatter_grid.png  Figure 2 — DM1_like vs THYROID_NONOVERLAP scatter, with Spearman
├── gpl570_direction_consistency_forest.png Figure 3 — Cohen d forest, color-coded for direction consistency
└── raw/                                    GITIGNORED (~50 MB Series Matrices, reproducible from URLs)

project/notebooks_or_scripts/p_gpl570_validation_first_pass.py    analysis script
project/scripts/build_gpl570_validation_html.py                   HTML rendering script
project/papers_hub_2026_05_04/gpl570_validation_pack.html         self-contained web page (1.1 MB, figures inline base64)
project/reports/2026_05_04_gpl570_validation_access_check.md      access check report
project/reports/2026_05_04_gpl570_external_validation_report.md   first-pass report
project/reports/2026_05_04_gpl570_validation_FULL_VIEW.md         this consolidated view
```

---

## 13. Live web

- **Top-level:** http://40.82.129.113/gpl570_validation_pack.html
- **Hub subpath:** http://40.82.129.113/papers_hub_2026_05_04/gpl570_validation_pack.html

Page contains: 4 KPI tiles · cohort composition · 3 figures with dense captions · 35-row headline contrast table (OK/WRONG colour-coded) · Spearman summary including anchor reference rows · 4 verdict cards · risk panels · file manifest. Self-contained — figures embedded as base64; no external network needed once loaded.

---

## 14. Commit group proposal (from prior turn — still pending user approval)

| Group | Files |
|---|---|
| **G1 — script + analysis outputs** | `p_gpl570_validation_first_pass.py`, all `project/results/p_gpl570_validation/*.tsv` `*.tsv.gz` `*.png` |
| **G2 — reports** | `2026_05_04_gpl570_validation_access_check.md`, `2026_05_04_gpl570_external_validation_report.md`, this FULL_VIEW.md |
| **G3 — web deliverable** | `build_gpl570_validation_html.py`, `papers_hub_2026_05_04/gpl570_validation_pack.html` |
| **G4 — gitignore** | `.gitignore` (one new line: `project/results/p_gpl570_validation/raw/`) |

**Excluded from commit:** `project/results/p_gpl570_validation/raw/` (3 Series Matrix .txt.gz, ~50 MB; gitignored, reproducible).

Approval phrases:
- `commit G1-G4` — ship everything in 4 logical commits
- `commit G1+G2+G4 only, hold G3` — ship analysis + reports + gitignore, hold the HTML hub page
- `do not commit, audit X` — revisit before staging

---

**End of consolidated view. GPL570 external validation first-pass complete. No additional datasets acquired. No manuscript prose drafted.**
