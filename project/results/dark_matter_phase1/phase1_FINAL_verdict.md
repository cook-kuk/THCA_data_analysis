# Dark Matter — Phase 1 FINAL Verdict
**Date:** 2026-04-29 (PM single-day sweep)
**Cohorts analyzed:** TCGA-THCA (482 with mutation calls), GSE241184 (3 samples / 1 patient sc)
**Cohorts queued:** Yoo 262 / K2 (mutation calls absent), 분당 (outreach), GSE33630 (not in warehouse)

---

## TL;DR

**GO** with the **reframed paper** (Frame B — Molecular Taxonomy), **NOGO** for the original prognostic Frame A.

- **Headline:** within driver-negative ("Dark Matter") thyroid cancer, an 8-gene transcriptional axis recapitulates the cPTC↔FVPTC histological dichotomy at single-cell resolution (**r = 0.905**), with DICER1/EIF1AX/PPM1D as the genomic anchor for the FVPTC-like subgroup (DM2; bulk Fisher **p = 0.0175**).
- **Drop:** prognostic claim. TCGA-THCA event scarcity (PFI 9 / OS 6 in DM cohort) makes survival HR untestable; multivariate cluster HR ≈ 1.03 with PFI.
- **Target venue:** Cell Reports Medicine / JCI Insight (primary). Reach: Nat Comm if Korean cohort + adult sc validation succeed.
- **Timeline:** 6-month bioRxiv → 9-month submission still feasible under Frame B.

---

## Cohort & data summary

| Cohort | Expr | BRAF/RAS calls | TERT | DICER1/EIF1AX | Fusion | 8-gene cluster | Survival |
|---|---|---|---|---|---|---|---|
| **TCGA-THCA (482)** | ✅ | ✅ | ✅ 36 | ⚠ 7 (mutation_genes) | ⚠ proxy | ✅ DM1 n=89 / DM2 n=55 | ✅ via Liu 2018 TCGA-CDR (PFI 50 / OS 16 events whole cohort) |
| K2 / PRJEB11591 (260) | ⚠ panel only (kallisto h5 in `/data/thca/PRJEB11591_quant_se/`) | ❌ | ❌ | ❌ | ❌ | ✅ predicted | ❌ |
| 분당 SNUH | ❌ outreach | ❌ | ❌ | ❌ | ❌ | ❌ | ❌ |
| GSE33630 | ❌ | — | — | — | — | — | — |
| GSE241184 (sc, 1 patient) | ✅ 30,493 cells post-QC | n/a | n/a | n/a | n/a | per-cell 8-gene score | n/a |

Master table: `tcga_dm_master_with_pfi.tsv` (482 rows, joined: mutation_groups + sample_master_v17 + clinical_extended + Liu 2018 CDR).

---

## Final per-step results

### Step 1 — Dark Matter prevalence ✅ PASS

| Cohort | n | BRAF V600E+ | RAS hotspot+ | **DM (BRAF−/RAS−)** | DM % |
|---|---|---|---|---|---|
| TCGA-THCA | 482 | 285 | 60 | **137** | **28.4%** |

Threshold ≥5% kill rule: **PASS**. Multi-cohort DM% pending K2 mutation calling.

### Step 2 — Within-DM survival HR ❌ FAIL (endpoint-limited)

| Endpoint | n | events | HR(DM2/DM1) | 95% CI | p |
|---|---|---|---|---|---|
| OS | 136 | 6 | 0.83 | 0.15-4.51 | 0.82 |
| PFI (Liu 2018) | 136 | 9 | 1.33 | 0.35-4.95 | 0.68 |
| DFI | 103 | 3 | 3.16 | 0.29-34.98 | 0.35 |
| DSS | 134 | 3 | 0.73 | 0.07-8.06 | 0.80 |
| **PFI multivariate (+age + stage)** | 133 | 8 | **1.03** | 0.25-4.16 | 0.97 |
| Within Xing low-risk (DM2 vs DM1) | 137 | 10 | 1.10 | 0.31-3.91 | 0.88 |

Cluster does NOT prognostically stratify within DM at any endpoint. TCGA-THCA structural limitation (whole-cohort 50 PFI events / 482 patients). Korean cohort survival required to revive prognostic claim.

### Step 3 — Alt-driver enrichment ✅ STRONG PASS

| Driver | DM1 (n=81) | DM2 (n=55) | OR | Fisher p |
|---|---|---|---|---|
| **DICER1 / EIF1AX / PPM1D** | **1 (1.2%)** | **6 (10.9%)** | **0.10** | **0.0175** |
| Any non-BRAF/RAS anchor | 4 (4.9%) | 7 (12.7%) | 0.36 | 0.12 (trend) |
| TERT promoter | 4 (4.9%) | 1 (1.8%) | 2.81 | 0.65 |

DM2 is **9× enriched** for DICER1/EIF1AX vs DM1.

### Step 3b — Continuous phenotype validation ✅ MASSIVE EFFECT

| Variable | DM1 (n=89) | DM2 (n=55) | Welch p | Cohen's d |
|---|---|---|---|---|
| RAI uptake gene score | 7.73 | 9.19 | **4.6e-16** | **1.54** |
| TDS16 (differentiation) | 6.93 | 8.11 | **4.0e-16** | similar |
| DICER1+ vs others (RAI) | — | — | **0.0023** | — |

DM2 is more differentiated and has stronger RAI-uptake signature than DM1.

### Step 3c — Histology ✅ TAXONOMY CONFIRMED

| Histology | DM1 | DM2 |
|---|---|---|
| cPTC | **64 (72%)** | 19 |
| FVPTC | 14 | **32 (58%)** |
| unknown | 11 | 4 |

DM1 = cPTC-enriched, DM2 = FVPTC-enriched within driver-negative tumors.

### Step 4 — sc Intra-tumor heterogeneity ✅ PASS

| Metric | Value |
|---|---|
| Total variance (8-gene score, thyrocytes n=7,259) | 0.133 |
| Mean within-sample variance | 0.130 |
| Between-sample variance | 0.038 |
| **Within / Between ratio** | **3.48×** |
| Tumor bimodality coef (n=2,427) | 0.32 (< 0.555) → **unimodal continuous gradient** |

ITH dominates inter-sample variation. The 8-gene panel is a continuous differentiation gradient at sc level, not a binary clusterer.

### Step 4b — sc headline finding (NEW, Figure 5D) ✅ STRONGEST RESULT

| Within tumor thyrocytes (n=2,427) | r |
|---|---|
| **8-gene score vs FVPTC signature** | **+0.905** (p < 1e-300) |
| 8-gene score vs cPTC signature | -0.363 |
| 8-gene score vs DICER1 pathway proxy | +0.001 |
| 8-gene score vs EIF1AX pathway proxy | -0.046 |

Replicated across samples: r = 0.900 (Normal) / 0.905 (Tumor) / 0.884 (LN_Met). The 8-gene panel captures a **universal thyrocyte differentiation axis** that aligns to FVPTC/cPTC histology. DICER1/EIF1AX bulk enrichment is co-occurrence at the genomic level — not a transcriptional pathway readout.

### Step 5 — Trajectory direction ✅ descriptive

| Sample | 8-gene | FVPTC | cPTC |
|---|---|---|---|
| Normal | 1.043 | 0.928 | 0.237 |
| Tumor | 0.769 | 0.619 | 0.587 |
| LN_Met | 0.669 | 0.539 | 0.602 |

`Normal → Tumor → LN_Met`: 8-gene ↓, FVPTC ↓, cPTC ↑ — classical dedifferentiation. Formal pseudotime requires multi-patient sc cohort.

### Step 6 — Xing axis rescue ⚠ size only, not prognostic

- Xing low-risk (BRAF−/TERT− ∪ BRAF+/TERT−) = 446 patients
- Of those, DM2 = 78 patients = **17.5%** "rescue rate" (count)
- **BUT** PFI events: DM1 6/82 (7.3%) vs DM2 4/55 (7.3%) → identical rate, logrank p=0.88

Counting metric passes; prognostic interpretation does not. Do not use Step 6 as prognostic evidence in the paper.

---

## Honest framing decision

| Frame | Claim | Evidence | Verdict |
|---|---|---|---|
| **A — Prognostic** (Nat Cancer ambition) | DM2 = high-risk substratifier in driver-negative thyroid cancer | None (PFI/OS/DFI/DSS all NS, multivariate HR≈1) | **DROP** |
| **B — Molecular Taxonomy** (Cell Rep Med / JCI Insight) | 8-gene axis recapitulates cPTC↔FVPTC at sc level (r=0.91); DICER1/EIF1AX = genomic anchor for FVPTC-like subgroup | r=0.91 sc, p=0.0175 DICER1 enrichment, p=4.6e-16 RAI/TDS, 28% prevalence, 17.5% Xing rescue size, ITH 3.48× | **GO** |
| C — RAI prediction | DM2 = high-RAI-uptake = treatable; DICER1/EIF1AX = treatable subtype | Same evidence as B | merge into B as clinical hook |

DM1 = "true Dark Matter" (cPTC-architectured driver-negative, unexplained, n=89). DM2 = "FVPTC look-alike with alternative drivers" (DICER1/EIF1AX-anchored, n=55). The mystery population is DM1, not DM2 — paper should foreground that asymmetry.

---

## Meeting talking points (3)

1. **Frame B is more defensible AND more interesting.** "8-gene reads FVPTC↔cPTC histology at r=0.91 in driver-negative thyroid cancer; DICER1/EIF1AX is the genomic anchor of the FVPTC-like cluster." Stronger than the prognostic story would have been even if it had worked.
2. **Prognostic drop is the right move.** PFI (9 events) + OS (6) + multivariate (HR=1.03) + Xing-rescue (logrank p=0.88) all consistent. TCGA-THCA cannot support the survival claim regardless of cluster definition. Korean cohort follow-up data is the only way back to a prognostic angle.
3. **Clinical hook for translation.** FNA Bethesda III/IV with ambiguous histology → 8-gene RNA score predicts cPTC vs FVPTC architecture without sequencing. Targets the 28% driver-negative subgroup where current molecular tests are non-informative.

---

## Figure 5 panel candidates (paper-ready, sc)

| Panel | File | Description |
|---|---|---|
| 5A | `fig5_sc/umap_overview.png` | All 30,493 cells, colored by sample × celltype × 8-gene score |
| 5B | `fig5_sc/hist_8gene_thyrocyte.png` | Density of 8-gene score in thyrocytes by sample (Normal>Tumor>LN_Met) |
| 5C | `fig5_sc/umap_thyrocytes.png` | Thyrocyte sub-cluster UMAP |
| **5D** | `fig5_sc/fig5D_headline_8gene_vs_fvptc.png/.pdf` | **HEADLINE — 8-gene vs FVPTC scatter colored by cPTC, r=0.905** |
| 5D-supp | `fig5_sc/fig5D_supp_by_sample.png` | Same scatter split by Normal/Tumor/LN_Met (universal r≈0.9) |
| 5E | derive from variance numbers | Within vs between bar (0.130 vs 0.038) |
| 5F | text annotation | Trajectory schematic Normal→Tumor→LN_Met |

---

## Files index (`project/results/dark_matter_phase1/`)

**Reproducible scripts:**
- `run_steps_1_to_6.py` — Step 1, 3, 6 (TCGA mutation × cluster)
- `step2_redo_pfi.py` — Step 2 with Liu 2018 CDR PFI/DFI/DSS/OS
- `step2_advanced.py` — multivariate Cox + Xing-rescue + alt-driver PFI
- `step2_rai_pivot.py` — Step 3b/3c (RAI score, TDS, histology by cluster)
- `sc_analysis.py` — Step 4-5 sc pipeline (load → QC → norm → annotate → score → ITH → corr)
- `make_fig5d.py` — Figure 5D headline scatter generator

**Tables:**
- `tcga_dark_matter_master.tsv` — 482 TCGA tumors with mutation + cluster
- `tcga_dm_master_with_pfi.tsv` — same + Liu 2018 CDR endpoints
- `step3_dm_cluster_alt_driver_enrichment.tsv` — Fisher table per alt-driver
- `step6_xing_rescue.tsv` — Xing-group × cluster rescue assignments
- `step2_dm_histology.tsv` — DM cluster × histology subtype
- `sc_cell_metadata.tsv` — 30,493 cells × per-cell scores + celltype + cluster

**Summaries:**
- `step1_to_6_summary.json` — D0 machine-readable
- `step2_pfi_results.json` — PFI/DFI/DSS Cox results
- `step2_advanced.json` — multivariate + Xing + alt-driver
- `step2_rai_pivot.json` — RAI/TDS pivot
- `sc_step4_5_summary.json` — sc machine-readable

**External data:**
- `tcga_cdr.xlsx` — Liu 2018 *Cell* PanCanAtlas Clinical Data Resource (downloaded today from GDC)

**Source documents:**
- `viability_dark_matter.md` — D0 deliverable
- `viability_dark_matter_step4_5.md` — D1-D2 sc deliverable
- `data_inventory.md` — cohort × variable matrix + blockers
- `phase1_FINAL_verdict.md` — **this document, single source of truth**

---

## Phase 2 priorities (D3+ ordered by impact)

1. **K2 RNA-seq mutation calling** (kallisto h5 in `/data/thca/PRJEB11591_quant_se/`) → BRAF/RAS Korean validation + cohort 2 DM% confirmation. **2-4 days** if matched normals exist.
2. **Adult multi-patient sc dataset** (GSE193581 Lu et al, or GSE164289 Wang) → external sc validation of r=0.91 FVPTC↔8-gene claim. **1-2 weeks.**
3. **Yoo 2016 supplementary follow-up** — search Yoo SK *Mol Ther* 2016 supplementary for Korean recurrence data. Author contact if needed.
4. **분당 outreach status check** — emails drafted but unsent per memory. Decide send/escalate this week.
5. **Full driver map (성공 4 viability)** — TCGA fusion classes × clinical phenotype. Can reuse `tcga_dark_matter_master.tsv`. **1-2 days.**
6. **Korean NRG1 (성공 3)** — separate trajectory; KoGES summary statistics if no germline access. Defer to Phase 3.

---

## Single-line conclusion

> **The 8-gene panel does not predict survival in driver-negative thyroid cancer (TCGA-THCA cannot test it), but it captures the FVPTC↔cPTC histological axis at single-cell resolution with r = 0.91 and identifies DICER1/EIF1AX as the genomic anchor for the FVPTC-like subgroup. Reframe the paper as molecular taxonomy with clinical translation through histology prediction.**
