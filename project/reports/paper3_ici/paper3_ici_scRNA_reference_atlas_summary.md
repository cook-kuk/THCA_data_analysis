# Paper 3 ICI — scRNA Reference Atlas Plan

**Working title:** HLA loss, neoantigen architecture, and immune ecotypes define ICI vulnerability in molecularly dark thyroid cancer
**Track:** A — design only. No integration runs. No Seurat / scanpy / scVI / scANVI / Harmony invocation.
**Author:** Seungho Cook
**Status date:** 2026-05-04

**Claim guard:** Atlas defines cell states for *signature derivation* and bulk projection. State frequencies in patients are descriptive, not predictive of ICI response in thyroid until DIAL audit + (where possible) thyroid ICI raw RNA-seq supports such a claim.

---

## 1. Cohorts & expected scale

Per dataset registry (`paper3_ici_dataset_registry.md` §2):

| Accession | Subtype skew | Expected cells (rough order) | Atlas role |
|---|---|---|---|
| GSE184362 | PTC + ATC | 30–80K | Anchor — broadest cell-type coverage |
| GSE193581 | PTC | 20–60K | Immune sub-states |
| GSE232237 | mixed | 20–50K | Integration validation |
| GSE191288 | PTC | 10–40K | Replication |
| GSE148673 | mixed | 10–30K | Replication |
| Han 2024 JCI Insight | dedifferentiation | 10–40K | Dedifferentiation anchor |

**Total estimate:** 100–300K cells after QC. Of these, PDTC + ATC fraction is expected to be small (<5K) — a critical caveat.

**Atlas-side risk register:**
- **R1.** PDTC/ATC scRNA cells <5K total → de novo dedifferentiation trajectory inference is underpowered. Mitigation: project pseudobulk-derived dedifferentiation signature from `paper3_ici_signature_registry` Module 7 onto epithelial cells; do not run Slingshot/Monocle3 on ATC alone.
- **R2.** Sample-level batch effects dominate disease-level effects in thyroid scRNA. Mitigation: scVI with `batch_key=sample_id`, treat patient-level as donor-level; never use disease as integration covariate.
- **R3.** Thyrocyte epithelial cells often contaminated by adjacent normal in surgical specimens. Mitigation: copy-number-based malignant-cell calling (inferCNV / CopyKAT) before downstream epithelial-state analysis.

---

## 2. Reference atlas dependencies (citation only — no runs)

Cell-type annotation will lean on published reference atlases via **transfer learning, not de novo discovery**. Each is referenced read-only — no re-derivation.

| Reference | Use | Tool |
|---|---|---|
| Zheng et al. 2021 (pan-cancer T-cell atlas) | T-cell sub-state labels | scANVI/celltypist |
| Andreatta et al. 2022 ProjecTILs | T-cell sub-state QC | ProjecTILs |
| Mulder et al. 2021 (myeloid pan-tissue) | Myeloid sub-state labels | celltypist |
| Lambrechts / Salcher pan-cancer NSCLC | Stromal CAF sub-states | celltypist |
| HCA thyroid (if available) | Thyrocyte normal reference | scANVI |

If HCA thyroid is unavailable, normal-thyroid cells from each cohort serve as the within-study normal reference.

---

## 3. Integration strategy

**Default:** scVI → scANVI (label-aware fine-tune for cell types) → Harmony as cross-check.

Sequence:
1. Per-cohort QC (scanpy default + `scrublet` doublet removal).
2. Concatenate raw counts; subset to top-3000 highly variable genes per cohort, take union (~5–8K HVGs).
3. scVI training — `n_layers=2`, `n_latent=30`, `batch_key=sample_id`. Single seed first; replicate with 3 seeds for stability.
4. scANVI fine-tune with reference cell-type labels from celltypist seeded annotation.
5. Harmony on the same HVG matrix as alternative integration; UMAP-level concordance reported.
6. Cluster (Leiden, multiple resolutions 0.3 / 0.6 / 1.0); pick resolution by silhouette × biological coherence.

**Acceptance bars:**
- kBET batch metric — improvement vs raw PCA must exceed 0.3.
- Marker-gene coherence — top 20 cluster markers must contain ≥3 known canonical markers for the assigned coarse cell type.
- LISI / iLISI — sample-mixing within cell-type bounds.

---

## 4. Cell-type taxonomy (locked at design)

**Coarse (Level 1) — cluster-then-annotate:**
- Epithelial (thyrocyte normal / malignant)
- T cell (CD4 / CD8 / Treg)
- B cell + plasma
- Myeloid (mono / macro / DC)
- Mast cell
- Endothelial
- Fibroblast / CAF
- Pericyte / smooth muscle

**Fine (Level 2) — within each Level 1, define operational sub-states:**

T-cell:
- Naive / central memory CD8 (TCF7+ LEF1+)
- Effector memory CD8 (GZMK+)
- Effector / tissue-resident CD8 (GZMB+ PRF1+ ITGAE+)
- Exhausted CD8 (PDCD1+ HAVCR2+ TOX+ — possibly CXCL13+ subset = `scCXCL13_CD8` — tumor-reactive proxy)
- CD4 Th1 (CXCR3+ TBX21+)
- CD4 Tfh-like (CXCL13+ ICOS+ BCL6+ — `scCXCL13_TFH`)
- Treg (FOXP3+)
- Proliferating T (MKI67+)

B-cell + plasma:
- Naive B (IGHD+)
- Memory B (CD27+)
- Germinal-center B (BCL6+ AICDA+ — `scTLS_BCELL_GC`) — *referenced for atlas; deeper BCR analysis is Paper 2 territory*
- Plasma cell (MZB1+ JCHAIN+)

Myeloid:
- Classical mono (S100A8+ S100A9+)
- Non-classical mono (FCGR3A+)
- Inflammatory TAM (CXCL9+ CXCL10+ — `scTAM_INFLAM`)
- M2-like / TREM2-high TAM (SPP1+ MRC1+ TREM2+ — `scTAM_M2_LIKE`)
- cDC1 (CLEC9A+ XCR1+)
- cDC2 (CD1C+)
- pDC (LILRA4+)
- LAMP3+ mature DC (mregDC)
- Mast (TPSAB1+)

CAF:
- Inflammatory CAF (CXCL12+ IL6+ — `scCAF_INFLAM`)
- Myo-CAF (ACTA2+)
- Antigen-presenting CAF (HLA-DR+) — flag as candidate Class II source distinct from myeloid

Epithelial (post-malignancy call):
- Normal thyrocyte
- Malignant — well-differentiated (TG+ TPO+)
- Malignant — dedifferentiated (TG-low, TDS-low — `scDEDIFF_EPI`)
- Malignant — IFN-responsive (HLA-DR-high, ISG-high)

---

## 5. Derived single-cell signatures → bulk projection

For each Level-2 sub-state, derive a signature by:
1. Wilcoxon DE of that sub-state vs all other cells of the same Level-1 type.
2. Take top 50 markers with adj. p < 0.01 and logFC > 0.5.
3. Filter against thyroid normal contamination markers (TG, TPO, TSHR top genes excluded).
4. Filter against Y-chromosome / mitochondrial / ribosomal genes.
5. Lock the gene list into `paper3_ici_signature_registry.md` §8 (currently placeholders).

Bulk projection method:
- ssGSEA on bulk cohorts using locked gene lists.
- Validation — predicted state fraction (CIBERSORTx with sc-derived signature matrix) correlated with ssGSEA score across samples; r > 0.5 required.

---

## 6. Headline scRNA outputs

- **Atlas UMAP** — Level 1 + Level 2 cluster overlay, dataset-source overlay, malignancy overlay.
- **State-frequency table** — per patient, fraction of each Level 2 state, separated by PTC vs PDTC vs ATC and by BRAF/RAS/fusion/dark-matter group.
- **Tumor-reactive CD8 (CXCL13+ exhausted) frequency** — per patient, by dedifferentiation tertile (TDS).
- **TLS-like B/Tfh co-occurrence** — per patient, evidence for organized B/Tfh niches; cross-validated against spatial cohorts.
- **TAM polarization landscape** — `scTAM_INFLAM` vs `scTAM_M2_LIKE` per patient.
- **Locked sc-derived signatures** populated into bulk registry.

---

## 7. What this atlas does *not* do

- **Does not** infer ICI response trajectories. State frequencies are descriptive.
- **Does not** re-analyze Paper 2's GSE286332 PTC+HT BCR / TLS / AICDA story. That dataset may appear as one cohort in the atlas, but state-level B/plasma analysis defers to Paper 2.
- **Does not** define DM1 / DM2 molecular clusters — that is Paper 1.
- **Does not** establish neoantigen-T-cell mapping at single-cell level (would require TCR-seq paired with WES; not consistently available across cohorts).

---

## 8. Outputs (filenames to be produced in Track B)

- `project/results/paper3_ici/scrna_atlas/atlas_obs.parquet` — cell-level metadata.
- `project/results/paper3_ici/scrna_atlas/atlas_markers.tsv` — Level 2 markers.
- `project/results/paper3_ici/scrna_atlas/sc_signatures_locked.json` — gene lists for §5.
- `project/results/paper3_ici/scrna_atlas/state_frequency_per_patient.tsv` — for ecotype linkage.
- `project/results/paper3_ici/scrna_atlas/figures/` — UMAPs, heatmaps.

---

Track A completed. No marathon violation.
