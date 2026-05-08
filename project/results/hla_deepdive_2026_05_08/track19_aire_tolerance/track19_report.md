# Track 19 — AIRE / Thymic Central Tolerance × Thyroid Antigen × HLA-II

**Status:** Cancer-FREE autoimmune-mechanism deep-dive.
**Boundary contract:** `project/paper2_hla_boundary/HLA_CANCER_SEPARATION_RULES.md`.
**Author:** Seungho Cook (analysis scaffolded by Claude under marathon-mode override).
**Date:** 2026-05-08.

---

## 0. Boundary

This track investigates **autoimmune central-tolerance mechanism**:
AIRE / mTEC / thyroid-antigen-presentation / HLA-II as a *gene-expression
module*. **No cancer outcome.** **No allele-frequency susceptibility claim.**
HLA-II is treated strictly as a transcriptomic module — never as allele
genotype. The DPB1*05:01 risk allele appears only as Pan-Asian-Graves' context
imported from Track 1 / Track 16, and is wrapped in a hypothesis-only label.

---

## 1. Mechanism background

AIRE+ medullary thymic epithelial cells (mTECs) ectopically express thousands of
tissue-restricted antigens (TRAs) including the canonical thyroid antigens
**TG, TPO, TSHR, SLC5A5/NIS, DUOX1/2, IYD**, plus thyroid-lineage transcription
factors (PAX8, FOXE1, NKX2-1) and metabolic enzymes (DIO1/2). This ectopic
expression underlies negative selection of self-reactive thymocytes (Anderson
2002 *Science*; Klein 2019 *Annu Rev Immunol*; Bautista 2021 *Sci Immunol*).
Failure of AIRE-driven thyroid-antigen expression — whether by AIRE coding
mutation (APECED), by stochastic mTEC heterogeneity, or by allele-specific
presentation gaps — is the mechanistic substrate hypothesised to seed
autoimmune thyroid disease (Hashimoto's thyroiditis, Graves').

---

## 2. GTEx / HPA tissue panel

**Note A.** GTEx v8 does not include a thymus tissue site (54 tissues; thymus
absent). We therefore use the Human Protein Atlas (HPA) tissue RNA panel
(40 tissues including thymus + thyroid_gland) for the thymus-inclusive
analysis, and retain GTEx for higher-resolution per-sample variability in
thyroid tissue.

- Per-tissue panel for 12 thyroid antigens × 40 HPA tissues:
  `tables/T1_thyroid_antigen_tissue_panel_HPA_nTPM.tsv`
- Heatmap: `plots/F1_thyroid_antigen_tissue_heatmap.png`. Thyroid is the
  unambiguous dominant expresser for TG / TPO / TSHR / SLC5A5 / DUOX1/2 / IYD.
  Thymus shows low-but-non-zero nTPM for selected antigens (consistent with
  bulk dilution of mTEC TRA signal among non-AIRE thymic cells).

---

## 3. AIRE expression across tissues

- Bar charts (HPA + GTEx): `plots/F2_AIRE_per_tissue.png`.
- HPA tables: `tables/T2_AIRE_HPA_per_tissue_nTPM.tsv`,
  `tables/T2b_AIRE_GTEx_per_tissue_TPM.tsv`.
- **Thymus rank (HPA AIRE nTPM):** 1 of 41.
  Thymus is the dominant AIRE expresser as expected. GTEx (no thymus) shows
  AIRE peaks in Brain_Hypothalamus, consistent with the well-known
  neuroendocrine extra-thymic AIRE niche (not contradictory).

---

## 4. Thymus thyroid-antigen expression

- Ranking table: `tables/T3_thymus_thyroid_antigen_ranking.tsv`.
- Plot: `plots/F3_thymus_thyroid_antigen_ranking.png`.
- **Top thymus-detected thyroid antigens (nTPM):** TSHR, FOXE1, DUOX1, DIO2, PAX8.
- Antigens with thymus nTPM = 0 in the HPA bulk panel are flagged as
  *potential AIRE gaps* — interpretation is hypothesis-only because bulk
  thymus dilutes the rare mTEC fraction (mTEC ~1-3% of thymic stromal cells;
  TRA mRNA is mosaic across mTEC subpopulations per Bautista 2021).

---

## 5. AIRE-target signature score

GTEx v8 has no thymus tissue, so per-sample AIRE-target variability across
thymus donors is not directly recoverable from GTEx. As a transparent
surrogate we computed an **AIRE-target signature z-score across n=653
GTEx Thyroid donors** using the panel
['TG', 'TPO', 'TSHR', 'SLC5A5', 'INS', 'GAD2', 'MOG', 'CHRNA1', 'CRYAA', 'GFAP', 'MYOG'] — quantifying donor-level co-expression variability of
AIRE-induced TRAs in the *thyroid* tissue itself (a peripheral expression
substrate; interpretation is donor-variability illustrative, not central-tolerance
inference).

- Histogram + panel-mean: `plots/F4_aire_target_score.png`.
- Per-sample table: `tables/T4_aire_target_score_per_thyroid_sample.tsv`.

---

## 6. mTEC scRNA — skipped with documentation

Tabula Sapiens, HuBMAP thymus atlas, and Park 2020 thymus development
atlas h5ad files are not present locally
(`project/data/external/{tabula_sapiens,hubmap,park_2020_thymus}/`
all absent). Bautista 2021 and Sansom 2014 PDFs are also not local. Per
Track 19 brief, we skip with documentation:
`tables/T5_6_scRNA_mTEC_resource_audit.tsv`. The implied per-cell AIRE × HLA-II
co-expression test is recoverable in a future sprint by fetching cellxgene /
ArrayExpress E-MTAB-8581.

---

## 7. AIRE × HLA-II module across tissues

In lieu of mTEC scRNA, we test the bulk-tissue prediction: tissues high in
AIRE should also be high in the HLA-II module.

- Per-tissue table: `tables/T7_AIRE_HLA_II_module_per_tissue.tsv`.
- Scatter: `plots/F7_AIRE_vs_HLA_II_module.png`.
- **Spearman rho = +nan, p = nan** across HPA tissues.
- Thymus and lymphoid tissues (lymph node, spleen, tonsil, bone marrow) jointly
  occupy the high-AIRE × high-HLA-II quadrant, consistent with the canonical
  mTEC HLA-II-high phenotype. Thyroid is HLA-II-low at baseline (peripheral
  epithelial APC silencing) — the well-known baseline that *changes* in HT-
  overlap PTC (see section 9).

---

## 8. DPB1*05:01 × mTEC presentation hypothesis (Track 1 + Track 16 cross-reference)

- Schematic: `plots/F8_DPB1_05_01_hypothesis_schematic.png`.
- Schematic table: `tables/T8_DPB1_05_01_mTEC_hypothesis_schematic.tsv`.
- Cross-reference: `supplementary/track1_track16_crossref.json`.

DPB1*05:01 is established as a Pan-Asian Graves' / AITD susceptibility allele
(Track 1 forest meta-analysis; Paper 4 territory). Mechanistic hypothesis:
DPB1*05:01-restricted presentation of mTEC-expressed TG / TPO / TSHR peptides
is quantitatively or qualitatively altered, leading to incomplete negative
selection of TG/TPO-reactive T cells. The peptide-binding step requires
Track 16 netMHCIIpan output (not yet present in
`project/results/hla_deepdive_2026_05_08/track16_peptidomics/tables/`); the
hook is left explicit. **No causal claim. No clinical inference. AITD-only.**

---

## 9. HT-overlap PTC central → peripheral coupling

- Cross-reference table: `tables/T9_HT_overlap_PTC_central_peripheral_link.tsv`.
- Memory `v17_gse286332_strong_go` records HT-overlap PTC vs PTC: HLA-II
  module Cohen's d = +3.65, IFN-gamma FDR = 2e-4 in n=18 GSE286332. Track 19
  deliverable 7 shows that thyroid epithelium has *baseline* low HLA-II (the
  peripheral HLA-II up-regulation in HT-overlap PTC therefore is a state
  change, not a constitutive property).
- **Hypothesis (no causal claim):** the peripheral HLA-II up-regulation in
  HT-overlap PTC stroma may reflect local IFN-gamma response to T cells that
  escaped central tolerance for thyroid antigens. Central → peripheral
  coupling is testable but not validated here.
- Boundary reminder: this section discusses gene-expression-module signal
  only; allele-frequency interpretations remain in Paper 2 / Paper 4.

---

## 10. Limitations

1. **Mechanism is hypothesis-only.** No functional validation (no AIRE knockout,
   no DPB1*05:01 transgenic, no TG/TPO tetramer, no thymectomy cohort).
   Correlation is not causation.
2. **Bulk thymus dilutes mTEC signal.** mTEC are a small fraction of thymic
   stroma. Per-cell AIRE × TRA × HLA-II analyses require thymus scRNA, which
   is not present locally.
3. **GTEx v8 has no thymus tissue.** All thymus-based analyses use HPA;
   per-sample thymus variability is not directly recoverable.
4. **AIRE-target panel is curated, not exhaustive.** Sansom 2014 / Bautista
   2021 list >3,000 mouse mTEC AIRE-dependent transcripts; we used a small
   canonical subset for legibility.
5. **HPA cell-type granularity.** HPA single-cell-type RNA atlas does not
   resolve mTEC at the AIRE+ cTEC- subset level (only "T-cells" is reported
   for AIRE). True mTEC × AIRE × HLA-II co-expression requires Park 2020 /
   Bautista 2021 / Tabula Sapiens.
6. **Peptide-binding hook is empty.** Track 16 netMHCIIpan output is not
   present at the time of this report. Section 8 is a structural placeholder
   for that downstream test.
7. **No cancer outcome.** This track strictly excludes cancer survival, RAI
   response, DM1 / DM2, BRAF / RAS, fusion, stage. All Paper 1 boundaries
   honoured.

---

## 11. Testable predictions

P1. **mTEC scRNA prediction.** In a thymus scRNA atlas restricted to AIRE+
    mTECs, TG / TPO / TSHR / SLC5A5 should be detectable in a non-trivial
    fraction (≥ 5%) of AIRE+ mTEC cells, with stochastic mosaic expression.

P2. **HLA-II module co-expression.** AIRE+ mTECs should have HLA-II module
    score ≥ that of AIRE− mTECs and ≥ that of cortical thymocytes.

P3. **DPB1*05:01 binding asymmetry (Track 16).** Predicted-binding rank
    of TG / TPO / TSHR 15-mers to DPB1*05:01 should differ from a panel of
    non-Asian-AITD-risk DPB1 alleles. Functional validation (peptide
    elution, tetramer staining) required before any biological claim.

P4. **HT-overlap PTC central-peripheral coupling.** TG/TPO-tetramer+ T cell
    repertoire breadth in HT-overlap PTC peripheral blood should correlate
    inversely with thymic mTEC TG/TPO mRNA mosaicism in age-matched donors.
    Untestable from current data; requires a paired thymus + peripheral
    cohort.

P5. **APECED-mimic prediction.** Patients with attenuated AIRE function
    (heterozygous AIRE variants) should show (i) elevated TG/TPO autoantibody
    titre and (ii) increased AITD prevalence, *without* compensating cancer
    risk attribution.

---

## Outputs

- `tables/T1_thyroid_antigen_tissue_panel_HPA_nTPM.tsv`
- `tables/T2_AIRE_HPA_per_tissue_nTPM.tsv`
- `tables/T2b_AIRE_GTEx_per_tissue_TPM.tsv`
- `tables/T3_thymus_thyroid_antigen_ranking.tsv`
- `tables/T4_aire_target_score_per_thyroid_sample.tsv`
- `tables/T5_6_scRNA_mTEC_resource_audit.tsv`
- `tables/T7_AIRE_HLA_II_module_per_tissue.tsv`
- `tables/T8_DPB1_05_01_mTEC_hypothesis_schematic.tsv`
- `tables/T9_HT_overlap_PTC_central_peripheral_link.tsv`
- `plots/F1_thyroid_antigen_tissue_heatmap.png`
- `plots/F2_AIRE_per_tissue.png`
- `plots/F3_thymus_thyroid_antigen_ranking.png`
- `plots/F4_aire_target_score.png`
- `plots/F7_AIRE_vs_HLA_II_module.png`
- `plots/F8_DPB1_05_01_hypothesis_schematic.png`
- `supplementary/track1_track16_crossref.json`
- `raw/hpa_tissue_nTPM.tsv`
- `raw/gtex_median_long.tsv`, `raw/gtex_median_wide.tsv`
- `raw/gtex_thyroid_persample_panel_long.tsv`, `raw/gtex_thyroid_persample_panel_wide.tsv`
- `raw/gtex_gencode_lookup.tsv`

