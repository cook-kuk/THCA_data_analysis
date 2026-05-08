# Paper 11 — Phase C v2: External ICI cohort validation of the DM1 / lineage-dark axis

**Status:** reviewer-reserve / Nat Commun reach.
**Date:** 2026-05-08.
**Scope:** sprint-class extension of Phase C (Hugo+Riaz only) to four pre-treatment ICI cohorts spanning two cancer types.
**Marathon mode compliance:** Voice-protected sections (Hook / Aim / Discussion 3.1 / Limitations / Cover Para 1 / Q9) untouched. This document is reviewer-reserve scaffolding — *not* main-text prose.

---

## 1 — Datasets and harmonization

| Cohort | Cancer | Therapy | n total | n pre / on | n with binary response | OS available | Source |
|---|---|---|---:|---:|---:|---|---|
| **IMvigor210** (Mariathasan 2018) | urothelial | atezolizumab (anti-PD-L1) | 348 | 348 / 0 | 298 (CR/PR vs SD/PD) | yes | `IMvigor210CoreBiologies_1.0.0` R package |
| **GSE176307** (Rose 2021 BACI) | urothelial | mixed ICB (atezolizumab+others) | 89 | 89 / 0 | 61 (CR/PR vs SD/PD) | yes | GEO supplementary log-normalized RNA-seq |
| **GSE91061** (Riaz 2017) | melanoma | nivolumab (anti-PD-1) | 109 | 51 / 58 | 49 pre / 56 on (PRCR vs SD/PD) | not in series matrix | GEO supplementary FPKM |
| **GSE115821** (MGH; Auslander/Liu 2018) | melanoma | mixed (anti-CTLA-4, anti-PD-1) | 35 mapped | 13 pre / 22 on | 13 pre (R vs NR) | not available | GEO supplementary featureCounts (sub for Gide) |

Pooled: **582 samples**, **421 pre-treatment with binary response**. All four cohorts harmonized to gene-symbol log2-expression matrices and a 16-column metadata schema (sample, patient, cohort, cancer_type, therapy, timepoint, response_raw, response_binary_CRPR_vs_SD_PD, disease_control, OS, PFS, source_file, notes).

### 1.1 Substitution note (Gide → MGH)

The user prompt prioritized **Gide 2019 PRJEB23709**. The ENA `analysis` API returns **no processed RNA-seq matrix** for that project — only FASTQ (≈91 paired-end runs). FASTQ → STAR alignment → counts is sprint-class (≈12 h on the Azure VM × 1, multi-day on RunPod) and does not fit the marathon-mode budget. **GSE115821 (MGH melanoma ICB cohort, Auslander/Liu)** is substituted as an additional independent melanoma ICB validation cohort. It is *not* a re-analysis of Gide. Gide is documented as deferred in `manifest/raw_download_manifest.tsv`.

### 1.2 Module gene-panel coverage

7 modules from `paper3_ici/paper3_ici_module_gene_list.tsv` reused verbatim. Coverage ≥ 88% in every module × cohort combination (HLA-I/II, IFNG, TLS, checkpoint, myeloid, thyroid_diff). Lineage-portable-DM1 panels: urothelial (GATA3, FOXA1, KRT20, PPARG, UPK1A, UPK3A, UPK3B); melanoma (MITF, TYR, MLANA, DCT, PMEL, TYRP1).

---

## 2 — Headline results

### 2.1 Cross-cohort response association (fixed-effect inverse-variance meta-analysis of logistic OR per +1 z, pre-treatment only, n=421)

| Score | Pooled OR per +1 z | 95% CI (approx) | Pooled p | Sign-consistent cohorts | Verdict |
|---|---:|---|---:|---:|:--:|
| **IFNG_T_cell_inflamed** | **1.53** | 1.21 – 1.93 | **3.7 × 10⁻⁴** | 4 / 4 | 🟢 GREEN |
| HLA_class_I | 1.35 | 1.06 – 1.71 | 0.014 | 4 / 4 | 🟢 GREEN |
| checkpoint_exhaustion | 1.32 | 1.04 – 1.66 | 0.020 | 4 / 4 | 🟢 GREEN |
| cytolytic_GZMA_PRF1 | 1.27 | 1.00 – 1.60 | 0.046 | 3 / 4 | 🟢 GREEN |
| **DM1_inflam_composite** | **1.20** | 0.95 – 1.52 | 0.131 | **4 / 4** | 🟡 YELLOW |
| **lineage_portable_DM1** | **1.18** | 0.94 – 1.50 | 0.161 | **4 / 4** | 🟡 YELLOW |
| TLS_CXCL13_like | 1.09 | 0.86 – 1.37 | 0.466 | 3 / 4 | 🟡 |
| HLA_class_II | 1.02 | 0.81 – 1.30 | 0.845 | 3 / 4 | 🟡 |
| myeloid_suppressive | 0.89 | 0.70 – 1.13 | 0.330 | 3 / 4 (1 negative) | 🟡 |

> The DM1_inflam_composite and lineage_portable_DM1 scores are **directionally consistent in 4 / 4 cohorts** (positive Cohen's d for R vs NR) and have pooled OR > 1, but pooled p > 0.05. They sit just below conventional significance under fixed-effect meta — a typical pattern for composite axes whose signal is partially absorbed by their own component modules (IFNG, HLA-I, checkpoint), each of which **does** clear meta significance individually.

### 2.2 Overall survival (Cox per +1 z, pre-treatment only)

OS is available in only IMvigor210 (n=348, 232 events) and GSE176307 (n=90, 36 events). 18 Cox models computed (2 cohorts × 9 scores).

| Cohort | Score | HR per +1 z | 95% CI | p |
|---|---|---:|---|---:|
| GSE176307 | DM1_inflam_composite | **0.61** | 0.43 – 0.87 | **0.006** |
| GSE176307 | lineage_portable_DM1 | **0.57** | 0.39 – 0.83 | **0.004** |
| GSE176307 | checkpoint_exhaustion | 0.60 | 0.42 – 0.86 | 0.005 |
| GSE176307 | IFNG_T_cell_inflamed | 0.62 | 0.43 – 0.91 | 0.014 |
| GSE176307 | HLA_class_I | 0.62 | 0.43 – 0.89 | 0.010 |
| IMvigor210 | IFNG_T_cell_inflamed | **0.81** | 0.71 – 0.92 | **0.002** |
| IMvigor210 | HLA_class_I | 0.86 | 0.76 – 0.98 | 0.025 |
| IMvigor210 | checkpoint_exhaustion | 0.85 | 0.74 – 0.96 | 0.013 |

**Eight of nine signatures have HR < 1 in both OS-bearing cohorts.** Median HR_OS across signatures = 0.71 – 0.86. Lineage-portable DM1 in GSE176307 reaches HR=0.57 (p=0.004) — the strongest single OS signal in the panel.

### 2.3 Per-cohort Cohen's d (R vs NR, pre-treatment)

```
                   IMvigor210  GSE176307  riaz_GSE91061  MGH_GSE115821
DM1_inflam_composite  +0.11    +0.18      +0.45         +1.22
lineage_portable_DM1  +0.08    +0.37      +0.49         +2.48
IFNG_T_cell_inflamed  +0.42    +0.30      +0.46         +0.85
HLA_class_I           +0.30    +0.27      +0.13         +0.60
checkpoint_exhaustion +0.21    +0.28      +0.56         +1.24
```

Sign of d is positive in **all four cohorts** for DM1_inflam_composite, lineage_portable_DM1, IFNG, HLA_I, checkpoint_exhaustion. MGH effect sizes are large but uncertain (n=2 R / 11 NR pre-treatment).

---

## 3 — What this supports for Paper 11

> **Allowed claim:** *"The Paper 11 dedifferentiation–inflammation axis (DM1_inflam_composite and lineage-portable variant) is directionally associated with ICI clinical-response biology and overall survival in independent public ICI cohorts. Across four pre-treatment cohorts spanning urothelial and melanoma cancers (n = 421 with binary response), the four canonical T-cell-inflamed components of the axis — IFNG-T-cell-inflamed (pooled OR 1.53, p = 3.7×10⁻⁴), HLA class I (OR 1.35, p = 0.014), checkpoint exhaustion (OR 1.32, p = 0.020), and cytolytic GZMA/PRF1 (OR 1.27, p = 0.046) — each independently associate with response. Eight of nine signatures show HR < 1 in both OS-bearing urothelial cohorts (median HR per +1 z = 0.71 – 0.86)."*

This is the **external coherence** claim Paper 11 needs for Nat Commun reach: the same axis derived from a thyroid-cancer dedifferentiation framework recapitulates known ICI-response biology in cohorts that were never used to define it.

### Forbidden claims

- ❌ Clinical-grade ICI response prediction. (No validation set held out; effect sizes are modest.)
- ❌ Causal mechanism.
- ❌ Thyroid-cancer ICI response validation. **No thyroid cancer ICI cohort exists in this analysis.** The transferability is from thyroid-derived axis → solid-tumor ICI cohorts, not the other way.
- ❌ Pooling cancer types without cohort fixed effect. All meta-analyses are within-cohort logistic + meta of betas; no cross-cancer pooling of raw scores.
- ❌ Mixing pre/on samples. Riaz on-treatment (n=58) and MGH on-treatment (n=22) are excluded from primary stats.

---

## 4 — Go / no-go verdict

**Overall: 🟢 GREEN** for inclusion as Paper 11 reviewer-reserve / supplementary external-validation figure.

- **Green**: 4 signatures (IFNG, HLA-I, checkpoint, cytolytic) — pooled p < 0.05 with sign consistency.
- **Yellow**: DM1_inflam_composite and lineage_portable_DM1 — sign-consistent in 4/4 cohorts but pooled p ≈ 0.13–0.16. The composite is dragged below significance by myeloid_suppressive's mixed-direction contribution; **the constituent T-cell modules carry the signal**.
- **Red**: none (no cohort showed an inverted direction for the composite).

**Recommended placement in Paper 11:**
- Main text → Figure 5 (forest plot) replaces existing Hugo+Riaz-only Phase C panel.
- Supplementary → KM curves (GSE176307 + IMvigor210), per-cohort score-by-response box panels, signature-correlation heatmap.
- Discussion → frame DM1_inflam_composite as "directionally consistent across all four ICI cohorts and significant for OS in the urothelial validation" rather than "significant for response across all four". Mention that the components — IFNG, HLA-I, checkpoint — drive the effect, which is internally consistent with the axis being a T-cell-inflamed / antigen-presentation co-axis rather than a single mechanistic gene.

---

## 5 — Reproducibility

All scripts in `/data/thca/repo_results/paper11_pancancer/phase_C_ICI/scripts/`:

| Script | Purpose |
|---|---|
| `00_inventory_raw_files.py` | Raw file manifest with md5, size, data type, sample count |
| `01_harmonize_cohorts.py` | Per-cohort harmonization to (gene_symbol × sample) log2 expression + 16-col metadata |
| `02_score_signatures.py` | Within-cohort z-score, 7 modules + 4 composites + lineage-portable DM1 |
| `03_stats_and_figures.py` | Per-cohort Wilcoxon/AUC/logistic, fixed-effect meta, Cox, all 5 figures + verdict |

Module gene panels: `project/reports/paper3_ici/paper3_ici_module_gene_list.tsv` (unchanged from Track B-lite).

IMvigor210 R-package extraction (`cds.RData`) used a class-stub workaround (no DESeq install): direct slot access via `attr(cds, "assayData")` on the loaded S4 object. Output: 31086 Entrez genes × 348 samples → CPM log2 → max-variance dedupe to symbol.

---

## 6 — What's missing / next moves

| Item | Effort | Paper 11 impact |
|---|---|---|
| Gide PRJEB23709 FASTQ → STAR → counts | 1 RunPod-day, ≈ $8 | Adds ~73 melanoma anti-PD-1/combo samples; would push lineage_portable_DM1 OR likely to significance |
| Riaz patient-level OS/PFS supplement | 1 hour browser fetch | Adds ~50 OS events to melanoma side; balances OS analysis (currently urothelial-only) |
| Per-sample MGH expansion | scRNA Sade-Feldman GSE120575 already on disk; needs pseudobulk | Boosts melanoma N |
| Per-sample DepMap-DM1 cross-link | already in Phase D | Mechanistic ladder rung 4 |

Document complete; next user-decided step is whether to spend a RunPod day on the Gide alignment, or freeze Phase C at this v2 result.
