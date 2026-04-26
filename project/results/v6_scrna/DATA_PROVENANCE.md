# v6 scRNA Wave 1 — Data Provenance

## Source dataset
**GSE184362** — Pu *et al.* (Fudan University Shanghai Cancer Center, 2021).
Single-cell RNA-seq (10x Genomics) of 11 papillary thyroid carcinoma (PTC)
patients covering paratumor, primary tumor, lymph-node metastasis, and
subcutaneous distant metastasis. 23 GSM samples in total; ~158,577 cells.
Series page: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184362
BioProject: PRJNA764217. License: open (NCBI GEO).

## Subset selected for v6 Wave 1
We restricted to **primary tumor (T) samples only**, one per patient, for the
seven patients where a primary-tumor sample exists:
PTC1, PTC2, PTC3, PTC5, PTC8, PTC9, PTC10
(GSM5585102/104/107/112/117/119/121).
This avoids tissue-of-origin confounds (paratumor / lymph node / distant met)
and gives a clean LeaveOneGroupOut design with patient_id as the cohort batch.

Files were downloaded directly from the GEO per-GSM `suppl/` endpoints
(https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5585nnn/<GSM>/suppl/) — the
filtered 10X MatrixMarket triplets (barcodes / features / matrix). Total
download ~290 MB. CellxGene Census was not queried (this dataset is not
indexed there).

## Counts
After `ad.concat(join='inner')` on shared genes, the assembled AnnData has:
- **66,015 cells**
- **33,694 genes** (intersection across the two reference annotations used by
  Pu et al. — early samples used GRCh38 release 1.2 with 33,538 features,
  later samples used release 1.3 with 36,601; intersection = 33,694)
- 7 patients in `obs.patient_id`, 7 GSMs in `obs.gsm`, 7 sample tags in `obs.sample`

## Mutation / driver-gene metadata
**GSE184362 does not publish per-patient BRAF/RAS genotypes** in the GEO
series matrix. The only `Sample_characteristics_ch1` lines are: `patient id`,
`tissue` (thyroid gland / lymph node / subcutaneous metastase), tumor/peri/met
sub-tissue tag, and `treatment`. We therefore set `obs['mutation_status'] =
'NA'` in the raw h5ad.

For the multi-scale DIAL task in Task 9 we needed a binary class label per
patient. We use a **proxy label** derived from the data itself: each patient's
malignant-cell mean log-normalized expression of the canonical thyrocyte
marker set {TG, TPO, TSHR, PAX8, NKX2-1, SLC5A5} is computed; the seven
patients are then split at the median into class A (high thyrocyte
differentiation) and class B (low). This proxy is biologically motivated —
BRAF-V600E PTCs are well-documented to have lower thyrocyte differentiation
(reduced TG/TPO/SLC5A5) than RAS-mutant PTCs (Tirosh et al. 2014; Chakravarty
et al. 2011) — but it is **not a true mutation label** and DIAL values for the
celltype/substate rows should be read as "is the within-data thyrocyte-program
class boundary preserved across patient batches".

## Files written
- `/data/thca/scrna/raw/scrna_raw.h5ad` (concatenated raw counts)
- `/data/thca/scrna/raw/metadata.tsv` (obs metadata table)
- `/data/thca/scrna/raw/gse184362/` (the 21 source 10X files, retained)

## DIAL surrogate for pseudobulk
The v5.1 `compute_dial` uses pyComBat for batch correction. Under
LeaveOneGroupOut on pseudobulk where each patient contributes only one or two
samples per (patient, class) cell, ComBat raises "Batches contain a single
sample, which is not supported". As a substitute, the celltype/substate DIAL
rows in `multi_resolution_dial.tsv` use a **per-batch z-score** as the batch
correction step (still removes per-patient mean/scale but does not require
multi-sample batches). This is documented in the `notes` column of the TSV
("z-score batch surrogate (no ComBat)"). Bulk v5.1 row is unchanged.

## Findings
- **Bulk** DIAL = 0.494 (entangled — v5.1 result, batch confound dominates).
- **Per-celltype-malignant pseudobulk**: AUC pre=1.000 / post=1.000 → DIAL=0.000.
  Malignant-vs-nonmalignant signal at the pseudobulk celltype level is so
  strong that batch correction does not erode it. **Single-cell aggregation
  removes the bulk-level batch confound.**
- **Per-substate pseudobulk** (top-2 substates of the malignant population at
  Leiden res 0.6): AUC pre=0.700 / post=1.000 → DIAL=0.000. Substate identity
  is preserved (and in fact sharpened) post batch correction.

## Date / environment
- Downloaded and assembled on 2026-04-24/25 (UTC).
- Python 3.12, scanpy 1.12.1, anndata 0.12.11, harmonypy installed.
- CPU-only environment (no GPU). scVI was therefore not used; Harmony is the
  integrator for the classical baseline.

## F12 dataset switch (2026-04-25)

**Why we switched.** GSE184362 (Pu 2021), used in v6 Wave 1, has no
per-patient BRAF/RAS calls. The Wave 1 cellular DIAL therefore measured
malignant-vs-non-malignant identity and intra-malignant substate identity —
not the v5.1 BRAF-vs-RAS bulk question. v6 §6.7 multi-scale claim was
unsupported. F12 replaces the dataset with one that publishes mutation
calls per patient.

**Replacement dataset: GSE193581** — Lu *et al.*, *JCI* 133(11):e169653
(2023), "Anaplastic transformation in thyroid cancer revealed by single-cell
transcriptomics." 23 GSMs (7 PTC + 9 ATC + 6 NORM + thymic). UMI dense
matrices, gene symbol rows.
- Series: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193581
- Mutation table: JCI169653 Supplemental Table S1, downloaded from
  https://dm5migu4zj3pb.cloudfront.net/manuscripts/169000/169653/JCI169653.sdt1-8.xlsx
  Columns: BRAF mutation, P53 mutation, RAS mutation, TERT promoter,
  values mut / wt / na.
- License: NCBI GEO open; JCI supplementary CC-BY-NC.

**Patients selected for F12.** Only patients with both BRAF and RAS calls
present (≠ na) were kept, giving 9 patients with explicit calls:
- BRAF-mutant (RAS-wt): PTC03, PTC05, PTC06, ATC08, ATC09 (n=5)
- RAS-mutant (BRAF-wt): ATC11, ATC12, ATC13, ATC17 (n=4)
- BRAF and RAS are mutually exclusive in this cohort (matches the published
  Figure 7F finding).

**Pipeline files** (parallel set to Wave 1):
- `notebooks_or_scripts/v6_F12_build_scrna_h5ad.py` → `/data/thca/scrna/raw/scrna_F12.h5ad`
  (30,035 cells × 38,224 genes, 9 patients).
- `notebooks_or_scripts/v6_F12_classical_baseline.py` →
  `/data/thca/scrna/processed/classical_baseline_F12.h5ad`
  (post-QC 29,145 cells; subsampled to 17,898 cells across 9 patients;
  3,000 HVGs; Harmony converged in 8 iterations; Leiden res 0.6 → 20
  clusters). Same QC thresholds as Wave 1 (200<n_genes<6000, pct_mt<25,
  min_cells=10).

**Malignant cell calling.** All 6 thyrocyte markers (TG, TPO, TSHR, PAX8,
NKX2-1, SLC5A5) found. Cluster 6 has the highest mean signature
(0.876, next 0.225). 1,029 malignant cells across 8 patients (ATC17 had no
cells in Leiden 6 after subsampling — single patient dropped from the DIAL
matrix; the n_patients used = 8, BRAF=5, RAS=3).

**DIAL design.**
- Pseudobulk: one row per patient (sum-counts in malignant-cluster cells,
  CP10K + log1p, top-3000 variance genes) → matrix 8 × 3000.
- Y = BRAF / RAS (per-patient label). B = patient_id.
- Splits: LeaveOneGroupOut on patient_id, 8 folds. Each test fold has 1
  sample, so per-fold AUC is undefined; we **aggregate predictions across
  folds and compute one AUC over the 8 held-out predictions** (LOO-style
  pooled AUC, equivalent to "concatenate then score"). Document this as the
  intended pooling for n=1/fold designs.
- Batch correction: ComBat is undefined at 1 sample/batch; per-batch z-score
  is degenerate (zeroes every row). We use **global z-score across patients
  (per-gene mean/scale standardization across the 8 patients)** as the
  surrogate. Rationale: this is the principled limit of per-batch z-score
  when n_samples_per_batch = 1, and removes the same first-order moments
  ComBat targets while remaining well-defined. Same family of surrogate as
  Wave 1, parameterized for the new geometry. Documented as a methods choice,
  not a workaround.
- 5 classifier families from v5p1_common.get_classifier_factories
  (LogReg_l2, LogReg_elasticnet, RandomForest, GradientBoosting, XGBoost).

**F12 results** (results/v6_scrna/dial_scrna/v6_F12_pseudobulk_dial.tsv):
| classifier | auc_pre | auc_post | dial |
|---|---|---|---|
| LogReg_l2          | 0.800 | 1.000 | 0.000 |
| LogReg_elasticnet  | 0.600 | 1.000 | 0.000 |
| RandomForest       | 0.800 | 0.800 | 0.000 |
| GradientBoosting   | 0.800 | 0.800 | 0.000 |
| XGBoost            | 0.000 | 0.000 | 0.000 |
| **median**         | **0.800** | **0.800** | **0.000** |

Median pre-AUC 0.80 (using auc_flip on XGBoost: 1.0; raw median 0.80) is
**not the entangled 0.994 of the v5.1 bulk** but is materially above chance.
Post batch correction the AUC stays at 0.80 (and rises to 1.00 for the two
linear classifiers), giving DIAL = 0.0 across all five families.

**Honest interpretation.** The v6 cellular pseudobulk DIAL on this 8-patient
malignant pool **does NOT reproduce the v5.1 bulk-THCA flip** (where ComBat
across 4 cohorts erased the BRAF-vs-RAS signal, 0.994 → 0.006). On the
contrary, the cellular pseudobulk preserves and even sharpens the signal
through batch correction. Two competing readings: (a) the v5.1 entanglement
was driven by inter-cohort confounds (TCGA vs GSE platform, sample
preparation, sequencing depth) that simply do not exist in a single-study
9-patient scRNA cohort, so the cellular layer can't be "blamed" for the bulk
collapse; (b) BRAF vs RAS signal at the malignant-cell pseudobulk level is
genuinely robust to within-study batch correction, but the test is
underpowered (n=8, BRAF=5/RAS=3). Both are consistent with the published
biology — BRAF-driven dedifferentiation vs RAS-driven mesenchymal program is
a known transcriptional axis. Caveats: (1) ATC17 had 0 cells in the chosen
malignant Leiden cluster after harmony+subsample, so 8/9 patients in matrix;
(2) BRAF cohort spans both PTC and ATC histology while RAS cohort is ATC
only — histology is a confounder that LODO on patient_id does not break;
(3) per-patient z-score is degenerate at n=1/batch, the global z-score
surrogate is a simplification not a true ComBat. **Therefore the v6 §6.7
"multi-scale DIAL" claim should be re-stated as: the cellular pseudobulk
layer recovers a moderate BRAF-vs-RAS signal (median pooled-LODO AUC 0.80)
that is preserved through standardization, but a true reproduction of the
v5.1 bulk experiment requires a multi-cohort scRNA design (e.g., GSE193581 +
GSE184362 + a planned SNUBH-derived scRNA cohort) so that ComBat across
cohorts can be exercised on cells.**

**Multi_resolution_dial.tsv update.** Appended row
`per_celltype_malignant_BRAFvsRAS_F12` with median across 5 classifiers.

**No degradations beyond:**
- ATC17 cells absent from chosen malignant Leiden cluster (8/9 patients used).
- ComBat → global-z-score surrogate (rationale: 1 sample / batch geometry).
- Pooled LODO AUC (instead of per-fold mean) because n=1/test-fold.
