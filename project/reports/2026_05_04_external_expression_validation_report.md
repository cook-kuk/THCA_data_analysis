# External expression validation — Paper 1 RAI-lineage / DM1-axis replication

**Date:** 2026-05-04
**Paper scope:** Paper 1 (DM1 molecular dark matter)
**Question:** Does the driver-orthogonal RAI-lineage / DM1-differentiation axis reproduce in independent thyroid expression cohorts?

Working directory: `project/results/p_external_expression_validation/`
Access check: `project/reports/2026_05_04_external_expression_access_check.md`
Wording mode: external expression validation • direction-consistent lineage-silencing • advanced-disease replication. **Not** clinical, survival, fusion, progression, or therapeutic-response validation.

---

## 1. Executive verdict

**REPLICATED.** Across 4 independent GPL570 cohorts (GSE33630 n=105, GSE29265 n=49, GSE65144 n=25, GSE53157 n=27), the RAI-lineage panel and the orthogonal `THYROID_NONOVERLAP` panel both drop in advanced disease (ATC, PDTC) relative to differentiated thyroid carcinoma (PTC/FVPTC/FTC) and to normal, with consistent direction and large effect sizes in three cohorts and the predicted (but underpowered) direction in the fourth. The DM1-like axis (`-RAI_8`) and the orthogonal lineage panel are tightly anti-correlated within each dataset (Spearman ρ between −0.84 and −0.94 in all four), supporting the claim that the axis is a real lineage-silencing signal rather than a panel artifact.

GSE126698 (Tier 1) is **dropped** at access — Series Matrix carries metadata only, raw is SRA-only, and raw FASTQ alignment is forbidden by spec.

---

## 2. Access summary

See `project/reports/2026_05_04_external_expression_access_check.md` for full table. Outcome:

| GEO | Tier | Status |
|-----|------|--------|
| GSE126698 | T1 (RNA-seq) | **DROP** — Series Matrix metadata-only (87 lines); supp files are DE summaries, no per-sample counts; raw=SRA-only; spec forbids alignment |
| GSE33630 | T2 GPL570 | **USED** — 105 samples (11 ATC / 49 PTC / 45 normal) |
| GSE29265 | T2 GPL570 | **USED** — 49 samples (9 ATC / 20 PTC / 20 paired-normal) |
| GSE65144 | T2 GPL570 | **USED** — 25 samples (12 ATC / 13 normal) |
| GSE53157 | T2 GPL570 | **USED** — 26 samples (5 PDTC / 8 FVPTC / 7 PTC / 4 FTC / 2 normal; 1 commercial pool dropped) |
| GSE53072 / GSE60542 / GSE120177 | T3 hold | not downloaded per user spec |

No additional GEO accessions queried, downloaded, or aligned.

---

## 3. Dataset / sample summary

`sample_metadata.tsv`:

| Dataset | n | ATC | PDTC | PTC | FVPTC | FTC | normal | dropped | Histology source |
|---------|---|-----|------|-----|-------|-----|--------|---------|------------------|
| GSE33630 | 105 | 11 | — | 49 | — | — | 45 | 0 | `Sample_source_name_ch1` ("non-tumor") + `Sample_characteristics_ch1` ("anaplastic"/"papillary") |
| GSE29265 | 49 | 9 | — | 20 | — | — | 20 | 0 | `Sample_title` ("Anaplastic …" / "Papillary …" / "Patient-matched non-tumor …") |
| GSE65144 | 25 | 12 | — | — | — | — | 13 | 0 | `Sample_source_name_ch1` ("Anaplastic Thyroid Carcinoma" / "Normal …") |
| GSE53157 | 27 | — | 5 | 7 | 8 | 4 | 2 | 1 commercial pool | `Sample_source_name_ch1` (`fresh-frozen <hist>`) |

Counts match the GEO landing-page declarations exactly.

---

## 4. Processed-data availability

All 4 used datasets have per-sample probe-level expression in the Series Matrix (~54,675 probes × N samples; GPL570 RMA-like log2 values). Auto-scale check in the script falls through (max < 50) — no re-log2 applied. Probe→gene collapse used the GPL570 GEO `annot.gz` (`Aug 09 2016`, 45,118 probes mapped). Multi-symbol probes (`A///B`) were expanded; per-sample, per-symbol expression was the **mean of all probes mapping to that symbol**. Result: 22,836 unique gene symbols across each dataset.

GSE126698 was **not** parsed for expression — only metadata.

---

## 5. Gene-panel coverage

`external_gene_coverage.tsv` — every panel gene was found in every used dataset:

| Panel | GSE33630 | GSE29265 | GSE65144 | GSE53157 |
|-------|----------|----------|----------|----------|
| RAI_8 (8 genes) | 8/8 | 8/8 | 8/8 | 8/8 |
| THYROID_NONOVERLAP (8) | 8/8 | 8/8 | 8/8 | 8/8 |
| TDS_TF (4) | 4/4 | 4/4 | 4/4 | 4/4 |
| MECHANISM (6) | 6/6 | 6/6 | 6/6 | 6/6 |
| OPTIONAL_TARGETS (12) | 12/12 | 12/12 | 12/12 | 12/12 |

Aliases used: `NKX2-1` ⇄ `TITF1`, `FOXE1` ⇄ `TTF2/FKHL15`, `SLC5A5` ⇄ `NIS` (none required for these four series — canonical symbols all present).

QC heatmap: `external_dataset_qc_heatmap.png`.

---

## 6. Per-dataset results

All scores are within-dataset z-scores (genes z'd within each dataset before averaging into panels). Key contrasts (full table: `external_score_tests.tsv`):

### 6.1 GSE33630 (n=105)
- **ATC vs PTC:** RAI_8 d=**−3.91**, p=6.8e-7; NONOVERLAP d=**−4.22**, p=4.6e-7; TDS d=−4.20; STAT3/AP1/DNMT d=+1.19, p=5.6e-3; TACSTD2 d=−1.58.
- **ATC vs normal:** RAI_8 d=**−5.48**, p=3.5e-7; NONOVERLAP d=**−6.72**, p=3.5e-7.
- **PTC vs normal:** RAI_8 d=−2.25, p=3.5e-13; NONOVERLAP d=−1.98 — partial RAI silencing already at PTC stage.
- DM1_like vs NONOVERLAP Spearman ρ=**−0.93**, p=9e-46 (n=105).

### 6.2 GSE29265 (n=49)
- **ATC vs PTC:** RAI_8 d=**−1.56**, p=2e-2; NONOVERLAP d=−1.51, p=1e-2; TDS d=−2.21; TACSTD2 d=−2.06; STAT3/AP1/DNMT d=+1.20.
- **ATC vs normal:** RAI_8 d=−2.29; NONOVERLAP d=−2.66 (p=4.5e-5).
- DM1_like vs NONOVERLAP Spearman ρ=−0.85, p=8e-15 (n=49).

### 6.3 GSE65144 (n=25)
- **ATC vs normal:** RAI_8 d=**−2.79**, p=4.0e-5; NONOVERLAP d=**−3.53**, p=3.2e-5; TDS d=−2.93; STAT3/AP1/DNMT d=+1.20, p=7e-3; TACSTD2 d≈0 (n.s.).
- DM1_like vs NONOVERLAP Spearman ρ=**−0.94**, p=7e-12 (n=25).

### 6.4 GSE53157 (n=26 after pool drop)
- **PDTC vs PTC:** RAI_8 d=−0.28 (n.s.); NONOVERLAP d=+0.15 (n.s.); TDS d=−0.20.
- **advanced (PDTC) vs DTC (PTC+FVPTC+FTC):** RAI_8 d=−1.05, p=0.09; NONOVERLAP d=−0.78, p=0.16; TDS d=−0.98 — direction consistent, **underpowered** (n=5 vs 19).
- **FVPTC vs PTC:** RAI_8 d=+1.50 (FVPTC retains higher RAI than classical PTC) — interesting but tangential to Paper 1; reported but not main.
- DM1_like vs NONOVERLAP Spearman ρ=−0.84, p=9e-8 (n=26).
- **Verdict:** weakest of the four; PDTC arm too small; assignment = supp/sensitivity, not main.

---

## 7. Direction consistency

`external_direction_consistency.tsv`, `external_direction_consistency_forest.png`.

Sign convention: positive d ⇒ higher in advanced (ATC/PDTC) than comparator.

| Score | advanced_vs_DTC sign across 4 datasets | advanced_vs_normal sign across 3 datasets (33630/29265/65144) |
|-------|----------------------------------------|------|
| RAI_8 | −, −, (no DTC), − | −, −, − |
| DM1_like | +, +, (no DTC), + | +, +, + |
| THYROID_NONOVERLAP | −, −, (no DTC), − | −, −, − |
| TDS / TF_collapse | −, −, (no DTC), − | −, −, − |
| STAT3_AP1_DNMT | +, +, (no DTC), − (53157, n.s.) | +, +, + |
| TACSTD2_z | −, −, (no DTC), − | + (33630, n.s.), −, ≈0 |

Direction is **consistent for the lineage panels** in every cohort where the contrast is testable. STAT3/AP1/DNMT trends positive (predicted direction) in 6/7 testable contrasts. TACSTD2 is heterogeneous (mostly negative — i.e. lower in advanced than DTC) and is **not claimed** as a positive marker for advanced disease in this report.

---

## 8. Comparison with GSE76039 (Landa 2016)

GSE76039 is the prior internal benchmark (PDTC/ATC vs PTC; cited in Discussion §3.1 — see `v17_landa2016_cite_save`). Conceptual link:

- Landa GSE76039 design: PDTC + ATC (advanced) RNA-seq.
- Present sweep adds **dichotomized expression-array replication** at the lineage-axis level, in 4 cohorts that are independent of GSE76039's cases.
- Effect sizes in the present sweep (RAI_8 d ≈ −1.6 to −5.5 in advanced contrasts) are concordant in sign with GSE76039 lineage signal and with the TCGA-derived 8-gene/RAI relationship described in the manuscript.

This is **expression replication only**. No survival/mutation claim is added or implied.

---

## 9. Which datasets are main / supp / internal / drop

| Dataset | Role | Reason |
|---------|------|--------|
| GSE33630 | **Main external replication panel** | Largest (n=105), full ATC + PTC + normal triad, strongest effect sizes |
| GSE65144 | **Main supporting** (ATC vs normal) | Cleanest 12 vs 13 ATC/normal contrast, RAI_8 d=−2.79 |
| GSE29265 | **Supp** (Sporadic vs Chernobyl) | Adds paired-normal architecture; RAI/NONOVERLAP both replicate but smaller effect sizes |
| GSE53157 | **Supp / sensitivity** | PDTC arm underpowered (n=5); direction-consistent only; not a main cohort |
| GSE126698 | **DROP** | Series Matrix metadata-only; raw FASTQ alignment forbidden |
| GSE53072 / GSE60542 / GSE120177 | **HOLD** (not used) | Per user spec — non-GPL570 / lymphoid confounding / Paper 9 scope |

---

## 10. Risks

- **Microarray vs RNA-seq:** all 4 used cohorts are GPL570 Affy U133+2 RMA values; the in-paper 8-gene index was derived on TCGA RNA-seq. Cross-platform agreement on direction is what is claimed; absolute coefficients are **not** transferable and were not transferred (per-dataset z-scoring was used instead of pooling).
- **Probe → gene collapse:** mean-of-probes was used. A few HGNC symbols (e.g. `NKX2-1`) historically map to legacy probe-set IDs under aliases (`TITF1`); the alias table was applied. No probe-set-vs-symbol curation was done — coverage is what GEO's Aug-2016 GPL570 annot.gz reports.
- **Histology label quality:** labels are taken verbatim from `Sample_source_name_ch1` / `Sample_title` / `Sample_characteristics_ch1`. No re-review of pathology was performed. GSE53157's 1 "commercial RNA pool" sample was dropped; otherwise no samples were re-classified. The 45 "patient-matched non-tumor control" samples in GSE33630 are treated as `normal`; some literature flags them as adjacent rather than histologically pristine — this caveat is carried forward.
- **No survival or mutation metadata** is being fetched or used. No survival, clinical-outcome, mutation, or therapeutic-response claim is made or implied from this sweep.
- **Cross-study batch effects:** datasets were never pooled before scoring. All effect sizes and tests are **within-dataset**. The "direction consistency" forest is a consensus *over independent studies*, not a meta-analytic pooled estimate.
- **GSE53157 PDTC underpower:** n=5 PDTC vs n=7 PTC; no contrast in this cohort survives multiple-testing correction. Reported as supportive only.
- **GSE126698:** spec-conformant SRA realignment would change the access list; it was **not** undertaken.

---

## 11. Safe wording to use in any downstream prose

- "external expression validation"
- "direction-consistent lineage-silencing"
- "advanced-disease replication"

## 12. Wording NOT to use (forbidden in this report and any derivative)

- "clinical validation"
- "survival validation"
- "fusion validation"
- "progression proven"
- "therapeutic response proven"

These claims are not supported by anything in this sweep.

---

## Output inventory

- `project/results/p_external_expression_validation/sample_metadata.tsv`
- `project/results/p_external_expression_validation/external_gene_coverage.tsv`
- `project/results/p_external_expression_validation/external_score_tests.tsv`
- `project/results/p_external_expression_validation/external_direction_consistency.tsv`
- `project/results/p_external_expression_validation/external_spearman.tsv`
- `project/results/p_external_expression_validation/external_sample_scores.tsv.gz`
- `project/results/p_external_expression_validation/{GSE33630,GSE29265,GSE65144,GSE53157}_expression_gene_log.tsv.gz`
- `project/results/p_external_expression_validation/external_rai_lineage_boxplots.png`
- `project/results/p_external_expression_validation/external_dm1_nonoverlap_scatter_grid.png`
- `project/results/p_external_expression_validation/external_direction_consistency_forest.png`
- `project/results/p_external_expression_validation/external_dataset_qc_heatmap.png`
- `project/results/p_external_expression_validation/scripts/run_external_validation.py`
- `project/results/p_external_expression_validation/scripts/make_figures.py`
- `project/results/p_external_expression_validation/raw/{GSE33630,GSE29265,GSE65144,GSE53157,GSE126698,GPL570}/…` (gitignored)
- `project/reports/2026_05_04_external_expression_access_check.md`
- `project/reports/2026_05_04_external_expression_validation_report.md` (this file)
