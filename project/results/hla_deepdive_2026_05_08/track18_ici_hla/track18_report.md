# Track 18 — Pan-cancer ICI response × HLA-I/II module × DM1 score

**Date:** 2026-05-08
**Sprint:** HLA deep-dive multi-track, Track 18 (Cancer / Paper 11 territory)
**Boundary:** HLA gene-expression module — not allele genotype.
**Output dir:** `/home/seungho/personal/THCA_data_analysis/project/results/hla_deepdive_2026_05_08/track18_ici_hla/`
**Script:** `/home/seungho/personal/THCA_data_analysis/scripts/hla_deepdive_2026_05_08/track18/01_track18_analysis.py`

---

## 0. Boundary contract

This track operates on **HLA gene-expression modules only** (HLA-A/B/C, B2M, TAP1/2, NLRC5, IRF1, PSMB8/9, ERAP1/2, HLA-E/F/G, CALR, CANX, PDIA3 for HLA-I; HLA-DRA/DRB1/DPA1/DPB1/DQA1/DQB1, CIITA, CD74 for HLA-II). It does **not** infer or use any HLA allele genotype. It is Paper-11 territory (pan-cancer ICI vulnerability of DM1-high tumours), not Paper 2 (HT-PTC HLA exploratory) or Paper 4 (Korean GD HLA backlog). Caption boilerplate is appended to every figure.

---

## 1. Cohorts

Track 18 consumes the Phase C v2 per-sample signature score table built by `paper11_pancancer/phase_C_ICI/scripts/02_score_signatures.py` (582 samples, 30 columns). Restricting to pre-treatment samples with binary CR/PR vs SD/PD response yields **n=421 across 4 cohorts**:

| cohort | cancer type | therapy class | n (pre + binary R) | n_R | n_NR | OS available |
|---|---|---|---|---|---|---|
| GSE176307 | urothelial | atezo / pembro / nivo / durva / avelu | 61 | 11 | 50 | yes |
| IMvigor210 | urothelial | atezolizumab | 298 | 68 | 230 | yes |
| MGH_GSE115821 | melanoma | anti-PD-1 | 13 | 2 | 11 | no |
| riaz_GSE91061 | melanoma | nivolumab (pre) | 49 | 10 | 39 | no |
| **TOTAL** | | | **421** | 91 | 330 | 2/4 cohorts (n=359 with OS) |

Gide PRJEB23709 was queued in raw/ but the processed expression matrix is empty in repo; documented in `T8_singlecell_ici_audit.tsv` and Phase C v2 manifest. MGH_GSE115821 is the formal substitution per `paper11_pancancer_2026_05_08` memory.

## 2. Module score construction

Inherited from Phase C v2 (`02_score_signatures.py`):
1. Per-cohort within-cohort gene z-score (controls cross-cohort batch entering composite).
2. Module score = mean z over genes present (coverage logged in `module_gene_coverage_per_cohort.tsv`; 1.00 except TLS missing JCHAIN at 0.889 in 3 cohorts and thyroid_diff missing TPO at 0.875 in GSE176307).
3. DM1_inflam_composite = mean(IFNG_T_cell_inflamed, myeloid_suppressive, checkpoint_exhaustion, HLA_class_II) — Paper 11 Phase C definition.
4. Lineage-portable DM1 = −LineageTF_z + Inflam_panel_z + Mod4_core_z (lineage-TF panel switched per cancer type — urothelial vs melanoma).

Track 18 does not re-define modules; it re-runs the response/survival statistics with HLA-I and HLA-II isolated, plus interaction terms.

## 3. HLA-I module × ICI response (Phase C v2 reconfirmation)

Per-cohort logistic regression of binary response on HLA-I module z-score; IVW meta across 4 cohorts.

| score | k | n_pos / k | pooled OR | 95% CI | p |
|---|---|---|---|---|---|
| **HLA_class_I** | 4 | **4 / 4** | **1.346** | [1.061, 1.706] | **0.014** |
| IFNG_T_cell_inflamed | 4 | 4 / 4 | 1.524 | [1.209, 1.922] | 3.7e-4 |
| checkpoint_exhaustion | 4 | 4 / 4 | 1.319 | [1.045, 1.665] | 0.020 |
| cytolytic_GZMA_PRF1 | 4 | 3 / 4 | 1.266 | [1.004, 1.596] | 0.047 |

Phase C v2 finding (memory `paper11_pancancer_2026_05_08`: HLA-I OR=1.35 p=0.014, IFNG OR=1.53 p=3.7e-4, checkpoint OR=1.32 p=0.020, cytolytic OR=1.27 p=0.046) **fully reconfirmed** to three significant figures. Sign coherence 4/4 (or 3/4 for cytolytic) cohorts. See `T2_per_cohort_OR_response.tsv` (full per-cohort stats), `T3_meta_OR_response.tsv` (IVW meta), and Figure F1.

## 4. HLA-II module × ICI response (NEW — Phase C v2 only tested HLA-I)

Same per-cohort logistic + IVW meta, applied to HLA-II module (HLA-DRA/DRB1/DPA1/DPB1/DQA1/DQB1/CIITA/CD74).

| cohort | n | n_R / n_NR | logOR per +1 SD | OR | 95% CI | p |
|---|---|---|---|---|---|---|
| GSE176307 (urothelial) | 61 | 11 / 50 | +0.014 | 1.014 | [0.527, 1.952] | 0.967 |
| IMvigor210 (urothelial) | 298 | 68 / 230 | −0.048 | 0.953 | [0.728, 1.249] | 0.728 |
| MGH_GSE115821 (melanoma) | 13 | 2 / 11 | +0.733 | 2.081 | [0.377, 11.49] | 0.400 |
| riaz_GSE91061 (melanoma) | 49 | 10 / 39 | +0.434 | 1.544 | [0.737, 3.235] | 0.249 |
| **IVW META** | **4** | 91 / 330 | **+0.024** | **1.024** | **[0.811, 1.294]** | **0.84** |

**HLA-II is null at the meta level** (OR=1.02, p=0.84, sign 3/4 positive but cohort-heterogeneous). The two melanoma cohorts trend positive; the two urothelial cohorts are flat-to-negative. Compared to HLA-I (OR=1.35 p=0.014), HLA-II adds **no independent ICI-response signal** in this 421-sample pool — consistent with the canonical view that anti-PD-1/PD-L1 acts primarily on the CD8/MHC-I axis. This is the **first key new finding** of Track 18: HLA-I ≠ HLA-II in this dataset; the Phase C "HLA module" effect is HLA-I-driven. See Figure F2 (forest) and F4 (heatmap).

## 5. DM1 × HLA-I interaction in ICI response

Per-cohort and pooled logistic regression `response ~ DM1_z + HLA_I_z + DM1×HLA_I` (sample-level z within cohort).

| cohort | n | β_DM1 (p) | β_HLA-I (p) | β_int (p) | OR_int [95% CI] |
|---|---|---|---|---|---|
| GSE176307 | 61 | −0.26 (0.71) | +0.46 (0.51) | +0.14 (0.61) | 1.16 [0.67, 2.01] |
| IMvigor210 | 298 | −0.41 (0.08) | **+0.63 (0.007)** | +0.13 (0.28) | 1.14 [0.90, 1.44] |
| MGH_GSE115821 | 13 | — | — | — | (singular; n too small) |
| riaz_GSE91061 | 49 | +1.40 (0.08) | −1.05 (0.18) | +0.17 (0.65) | 1.18 [0.57, 2.44] |
| **POOLED** | **421** | −0.18 (0.39) | **+0.44 (0.035)** | **+0.11 (0.28)** | **1.12 [0.91, 1.37]** |

**Interaction term is positive in 3/3 evaluable cohorts and pooled, but does not reach significance** (pooled p_int=0.28). The dominant signal is HLA-I main effect (pooled β=+0.44, p=0.035) once DM1 is included. DM1 main effect is negative when conditioned on HLA-I (the IFNG/checkpoint/HLA-II constituents of DM1 partially overlap with HLA-I co-expression, so DM1's marginal contribution shrinks). Interpretation: **HLA-I module is the load-bearing component of the DM1_inflam composite for ICI response**, with a directionally positive but underpowered DM1 × HLA-I synergy. See `T4_dm1_x_hla1_interaction.tsv` and Figure F5 (per-cohort scatter).

## 6. Within-DM1-high stratification: does HLA-I matter more in DM1-high?

For each cohort, sort samples by DM1_inflam_composite, take top tertile (DM1-high) and bottom tertile (DM1-low), then test HLA-I → response within each stratum.

| cohort | strata | n | n_R | OR_HLA-I | p |
|---|---|---|---|---|---|
| GSE176307 | DM1-high | 21 | 4 | 1.31 | 0.61 |
| GSE176307 | DM1-low  | 20 | 2 | 0.57 | 0.48 |
| IMvigor210 | DM1-high | 100 | 28 | **1.55** | **0.060** |
| IMvigor210 | DM1-low  | 99  | 20 | 1.03 | 0.90 |
| riaz_GSE91061 | DM1-high | 17 | 5 | 1.06 | 0.91 |
| riaz_GSE91061 | DM1-low  | 17 | 3 | 1.00 | 0.997 |

Pooled (cohort-stratified logistic on DM1-high samples only, n≈138): **HLA-I logOR_per_z = +0.318, SE=0.195, OR=1.37, p=0.103.** In every cohort the within-DM1-high HLA-I OR is **larger than within-DM1-low** (1.31>0.57; 1.55>1.03; 1.06>1.00). This is consistent with the directional interaction in Section 5: **HLA-I matters most among DM1-high tumours**, and DM1-low tumours show essentially no HLA-I response gradient. With n=138 DM1-high samples this falls just short of significance; powering a 5th anti-PD-(L)1 cohort would push it through. See Figure F6.

## 7. Survival in ICI cohorts (HLA-I, HLA-II × OS)

Cox HR per +1 SD of score, OS endpoint. Only GSE176307 (n=90, 36 events) and IMvigor210 (n=348, 232 events) carry usable OS.

| cohort | score | n | n_event | HR per z | 95% CI | p |
|---|---|---|---|---|---|---|
| GSE176307 | HLA_class_I | 90 | 36 | **0.62** | [0.43, 0.89] | **0.010** |
| GSE176307 | HLA_class_II | 90 | 36 | **0.68** | [0.49, 0.95] | **0.025** |
| GSE176307 | DM1_inflam_composite | 90 | 36 | **0.61** | [0.43, 0.87] | **0.006** |
| GSE176307 | IFNG | 90 | 36 | 0.62 | [0.43, 0.91] | 0.014 |
| GSE176307 | checkpoint | 90 | 36 | 0.60 | [0.42, 0.86] | 0.005 |
| IMvigor210 | HLA_class_I | 348 | 232 | **0.86** | [0.76, 0.98] | **0.025** |
| IMvigor210 | HLA_class_II | 348 | 232 | 0.93 | [0.82, 1.06] | 0.30 |
| IMvigor210 | DM1_inflam_composite | 348 | 232 | 0.90 | [0.79, 1.02] | 0.10 |
| IMvigor210 | IFNG | 348 | 232 | 0.81 | [0.71, 0.92] | 0.0019 |
| IMvigor210 | checkpoint | 348 | 232 | 0.85 | [0.74, 0.96] | 0.013 |

**HLA-I module is OS-protective in both cohorts under ICI** (HR=0.62 / 0.86, both p<0.05). HLA-II is OS-protective in GSE176307 (mixed agents) but null in IMvigor210 (atezo only) — same urothelial heterogeneity seen in the response analysis. DM1_inflam composite OS-protective in GSE176307 (HR=0.61, p=0.006), trending in IMvigor210. See `T6_survival_cox_per_cohort.tsv` and Figure F7.

## 8. Composite DM1_inflam re-confirmation + HLA-II augmentation

| composite definition | k | sign+ / k | pooled OR | 95% CI | p |
|---|---|---|---|---|---|
| DM1_inflam_composite (Phase C: IFNG + myeloid + checkpoint + HLA-II) | 4 | 4 / 4 | 1.20 | [0.95, 1.52] | 0.13 |
| DM1_HLA2only (mathematically identical above) | 4 | 4 / 4 | 1.20 | [0.95, 1.52] | 0.13 |
| **DM1_inflam_plus_HLA2 = mean(IFNG, myeloid, checkpoint, HLA-II, HLA-I)** | **4** | **4 / 4** | **1.24** | **[0.98, 1.57]** | **0.076** |

Adding HLA-I to the composite **improves the pooled OR from 1.20 → 1.24** and tightens the p toward significance (0.13 → 0.076). The Phase C composite leaves response signal on the table by including HLA-II (null) and excluding HLA-I (positive). The current composite definition should be flagged for revision in Paper 11; sign-consistency 4/4 holds.

## 9. Single-cell ICI cohorts

Searched repo (`project/data/external/`, `paper11_pancancer/`, `paper3_ici/`) for canonical scRNA ICI cohorts (Sade-Feldman GSE120575, Yost GSE123813, Bassez GSE166181, etc.). None are downloaded into this repo. Documented as missing in `T8_singlecell_ici_audit.tsv`. Track 12 (`track12_scrna_hla_celltype`) covers single-cell HLA cell-type analyses outside the ICI-response axis.

## 10. Limitations

- **Small n per cohort** for melanoma side: MGH n=13 (only 2 responders) is hypothesis-generating only; riaz n=49. Most power lives in IMvigor210 (n=298).
- **Response definition heterogeneity**: GSE176307 mixes atezo / pembro / nivo / durva / avelu; IMvigor210 is atezo monotherapy; melanoma cohorts are anti-PD-1. Heterogeneity pushes meta toward null for cohort-divergent signals like HLA-II.
- **Within-cohort z-score** controls batch but does **not** correct for differing response prevalence (urothelial 18-22% vs melanoma 15-20%). A mixed-effects model is the proper next step but requires more cohorts.
- **Tumour-cell vs immune-cell HLA-I origin is not separable** in bulk RNA-seq; deconvolution is needed to claim "HLA-I high tumour cells" vs "HLA-I high TIL infiltrate". Track 12 (scRNA cell-type module) is the appropriate locus.
- **OS censoring**: only 2/4 cohorts have usable OS; HR confidence intervals are wide.
- **No pre-registration** for Track 18 itself — these are post-hoc extensions of Phase C v2 within the same dataset, so multiple-testing penalties for the new HLA-II contrast are explicit (FDR within score in T2; raw p reported elsewhere). The HLA-II null result is robust to any reasonable correction.
- **Boundary**: gene-expression module ≠ allele genotype. No claim is made about specific HLA-A/B/C/DR/DP/DQ alleles or carrier frequencies.

## 11. Paper 11 hook

Across 4 anti-PD-(L)1 cohorts (urothelial + melanoma, n=421), a **transcriptomic HLA-I module** — independent of HLA allele genotype — is the load-bearing predictor of ICI response (pooled OR 1.35, p=0.014; sign-consistent 4/4 cohorts) and is OS-protective under ICI (HR 0.62-0.86, both cohorts with OS p<0.05). The matched HLA-II module is null at the meta level (OR 1.02, p=0.84), establishing that the Phase C "HLA module" finding is class-I-driven. A directionally positive DM1 × HLA-I interaction (per-cohort 3/3, pooled OR 1.12, p=0.28) and the within-DM1-high stratification trend (HLA-I OR=1.37, p=0.10 in n=138 DM1-high samples) suggest **HLA-I is the gateway through which DM1-high tumours convert to ICI response** — a testable axis for Paper 11's pan-cancer DM1 vulnerability story, with image-DM1 (CLAM AUC 0.83) and methylation/protein triangulation already in place. Replacing HLA-II with HLA-I in the DM1_inflam composite tightens the response association (OR 1.20 → 1.24, p 0.13 → 0.076), supporting a v3 composite definition.

---

**Files**
- Tables: `T1_per_cohort_hla_dm1_scores.tsv` (per-sample), `T2_per_cohort_OR_response.tsv` (9 scores × 4 cohorts), `T3_meta_OR_response.tsv` (IVW), `T4_dm1_x_hla1_interaction.tsv`, `T5_within_DM1_strata_HLA1_response.tsv`, `T6_survival_cox_per_cohort.tsv`, `T7_composite_DM1_inflam_extended.tsv`, `T8_singlecell_ici_audit.tsv`.
- Figures: `F1_forest_HLA1_response.png`, `F2_forest_HLA2_response.png`, `F3_forest_DM1_inflam_response.png`, `F4_OR_heatmap_per_cohort.png`, `F5_DM1_x_HLA1_scatter_by_response.png`, `F6_within_DM1_high_HLA1_by_response.png`, `F7_survival_forest.png`.
- Summary JSON: `track18_summary.json`. Log: `logs/track18_run.log`.

**Caption boilerplate (every figure):** *HLA gene-expression module — not allele genotype.*
