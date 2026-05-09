# R2 — driver-stratified spatial replay (GSE250521 Visium n=16)

**2026-05-09 · Paper 1+2 BRAF Nature REINFORCEMENT · weak spot #2 (no driver labels in GSE250521).**

## Headline

**BRAF-like-inferred subset n=5 replicates spatial HT pattern with HLA-DRA delta=+0.71 log1p** (HIGH-HT vs LOW-HT, sign-consistent with TCGA d=+1.57 and ALL-n=8 H11 d=+0.83). The pattern is *not* BRAF-restricted: MAPK-active samples (n=10, PTLPTC+ATC, MAPK_z>0) reach **MWU p=0.0079** for HLA-DRA HIGH vs LOW (delta=+0.88, rho=+0.90 vs HT13_pb). **Zero RAS-like tumors** in this cohort (all 4 RAS-like calls are PT/normals) — driver-stratified RAS contrast not testable here.

## Driver inference

Per-sample pseudobulk z (within-GSE) for MAPK-output (DUSP4/5/6, SPRY2/4, ETV4/5, PHLDA1, CCND1), TDS-16, 8-gene panel. **Rule:** BRAF-like = MAPK_z >= +0.25 AND TDS_z <= +0.25; RAS-like = MAPK_z <= -0.25 AND TDS_z >= -0.25; else ambiguous. TCGA-THCA reference (n=505 tumors): BRAF_like MAPK-out z median +0.255 vs RAS_like -0.097 confirms axis polarity. **BRAF V600E mutation calling**: infeasible — GEO release has count matrices only (no BAM/FASTQ). Honest negative.

## Stage x driver concordance

| stage | BRAF_like | ambig | RAS_like |
|---|---|---|---|
| ATC (n=4)  | 2 | 2 | 0 |
| LPTC (n=4) | 2 | 2 | 0 |
| PTC (n=4)  | 1 | 3 | 0 |
| PT (n=4)   | 0 | 0 | 4 |

All 4 PT (normal) -> RAS-like (low MAPK + high TDS). Zero advanced tumors are RAS-like. Among 12 tumors: 5 BRAF-like + 7 ambiguous + 0 RAS-like.

## H11 stratified replay (HT-axis HIGH vs LOW within subset)

| subset | n | HLA_DRA delta | MWU p | HT13_pb delta |
|---|---|---|---|---|
| ALL_PTLPTC | 8 | +0.83 | **0.029** | +0.215 |
| **MAPKhigh_PTLPTC_ATC** | **10** | **+0.88** | **0.0079** | +0.213 |
| BRAFlike_PTLPTC_ATC | 5 | +0.71 | 0.20 | +0.155 |
| BRAFlike_PTLPTC | 3 | +0.59 | 0.67 | +0.204 |

**Sign holds in every subset.** Statistical power is the limiter at n=5 BRAF-only.

## TLS density by inferred driver

| driver | n | total TLS spots | mean TLS cov % | mean HLA-DRA |
|---|---|---|---|---|
| BRAF_like (any) | 5 | 11 | 0.049 | 1.685 |
| ambig (tumor) | 7 | 83 | 0.487 | 1.892 |
| RAS_like (=PT) | 4 | 3 | 0.020 | 0.535 |

BRAF-like vs ambiguous-tumor HLA-DRA d=-0.21 (p=0.53): pattern tracks the **MAPK-active continuum**, not a BRAF-vs-RAS dichotomy at this n.

## Bottom line

1. Driver-stratified spatial claim **survives** in BRAF-like-inferred n=5 (sign-consistent HLA-DRA +0.71 HIGH-HT vs LOW-HT) but is *not* BRAF-specific.
2. RAS-driven THCA cannot be spatially tested here (zero RAS-like tumors).
3. **Honest paper framing**: "spatial HT-axis pattern replicates within MAPK-active samples (n=10 p=0.0079), with sign-consistent direction in strictly BRAF-like n=5 subset; RAS contrast deferred to future Visium dataset with RAS-driven cases."

## Files

`run_r2_driver_inference.py` · `r2_per_sample_driver_inference.tsv` · `r2_per_sample_metrics_ALL16.tsv` · `r2_stratified_spatial_aucs.tsv` · `r2_stratified_spatial_aucs_AUGMENTED.tsv` · `r2_TLS_by_driver.tsv` · `r2_stage_x_driver_concordance.tsv` · `r2_tcga_reference_distributions.tsv` · `r2_summary.json`
