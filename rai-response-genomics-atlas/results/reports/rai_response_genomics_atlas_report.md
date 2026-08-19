# RAI Response Genomics Atlas — initial build report

**Compiled**: 2026-05-21
**Workspace**: `/home/seungho/personal/THCA_data_analysis/rai-response-genomics-atlas/`
**Data**: `/data2/rai_atlas/` (Premium SSD, symlinked as `data/`)

This first build delivers the scaffolding, eight-gene panel, dataset inventory, scripts, Figure 1 (multi-gate RAI model), and a working first validation against GSE151179.

## 1. What was created

```
rai-response-genomics-atlas/
├── README.md                              # how-to + workspace overview
├── requirements.txt                       # GEOparse, scikit-learn, etc.
├── config/
│   ├── datasets.yaml                      # 10 priority datasets, Tier-mapped
│   └── eight_gene_panel.yaml              # canonical Paper 1 panel (DIO1 included)
├── docs/
│   ├── literature_review.md               # 10-section scientific background
│   ├── dataset_inventory.md               # human-readable inventory
│   ├── rai_response_label_taxonomy.md     # 5-tier label hierarchy (Tier 1 → Tier 5)
│   ├── manuscript_strategy.md             # 3 candidate titles + figure plan + journal ladder
│   ├── reviewer_risk_register.md          # 10 anticipated reviewer concerns with defenses
│   ├── data_request_boucai.md             # Tier 1 data-request draft
│   └── data_request_mu_hra004166.md       # Tier 2 controlled-access request draft
├── scripts/
│   ├── 00_setup_environment.py            # env audit → logs/env.txt
│   ├── 01_search_geo_metadata.py          # print priority datasets
│   ├── 02_download_geo_dataset.py         # GEO download with GEOparse + FTP fallback
│   ├── 03_parse_geo_metadata.py           # parses RAI / iodine tokens + key:value pairs
│   ├── 04_score_gene_panel.py             # within-cohort z-score of 8 panel genes
│   ├── 05_validate_rai_labels.py          # MWU, Cohen d, ROC AUC, Kruskal
│   ├── 06_make_figures.py                 # boxplot + ROC + heatmap
│   ├── 07_build_dataset_inventory.py      # emit dataset_inventory.csv
│   ├── _extract_gpl23159_probe2gene.py    # Clariom S annotation helper
│   └── figure1_multigate_rai_model.py     # Nature Medicine-style Fig 1
├── notebooks/
│   ├── 01_gse151179_validation.ipynb      # DONE (results documented)
│   ├── 02_gse299988_validation.ipynb      # placeholder + plan
│   ├── 03_tcga_thca_discovery_plan.ipynb  # placeholder (re-use Paper 1 outputs)
│   └── 04_redifferentiation_dataset_scan.ipynb  # placeholder
├── data → /data2/rai_atlas/               # symlink
│   ├── raw/GSE151179/                     # SOFT family + series matrix + GPL annot
│   ├── interim/GSE151179_metadata.tsv     # 52 samples × 61 metadata cols, RAI fields parsed
│   ├── interim/GSE151179_expression.tsv   # 27,189 probes × 52 samples
│   └── processed/
│       ├── GSE151179_panel_score.tsv      # panel z + n_genes_available
│       └── GSE151179_label_joined.tsv     # panel + metadata + parsed labels
└── results/
    ├── tables/
    │   ├── dataset_inventory.csv          # 10 datasets, all columns per spec
    │   └── GSE151179_test_results.tsv     # MWU + AUC + Kruskal
    ├── figures/
    │   ├── figure1_multigate_rai_model.png  # Nature Medicine-style Fig 1
    │   ├── figure1_multigate_rai_model.pdf
    │   ├── GSE151179_boxplot_panel_by_label.png
    │   └── GSE151179_roc_refractory_vs_avid.png
    └── reports/
        ├── GSE151179_validation_brief.md
        └── rai_response_genomics_atlas_report.md  ← this file
```

## 2. Eight-gene panel (canonical, Paper 1 confirmed)

`SLC5A5 · TPO · TG · TSHR · PAX8 · NKX2-1 · FOXE1 · DIO1`

- Verified against `project/results/p3_gse286332/8gene_panel_per_sample.tsv` and `project/submission/npj/tables/AMP4_8gene_model_coefficients.tsv`.
- **DIO1 is in, DUOX2 is not** (corrected from the early specification).
- Polarity rule (from Paper 1 R17 memory): panel_z high = preserved RAI / DM2; panel_z low = silenced / DM1.

## 3. Dataset coverage at this build

10 priority datasets indexed; one fully run end-to-end:

| Dataset | Tier | Status |
|---|---|---|
| GSE151179 | 2 (avidity) | **Pipeline complete**. Tumor vs non-neoplastic AUC=0.96, p=9.5e-7. Refractory vs avid underpowered (avid n=4). |
| GSE299988 | 2 (avidity) | Notebook placeholder + plan. Run pipeline next. |
| Boucai 2023 CCR | 1 (RECIST) | Data-request email drafted. Supplement scan recommended first. |
| Mu 2024 HRA004166 | 2 (4-class avidity, gray-zone) | Controlled-access plan + email draft. Supplement scan recommended first. |
| TCGA-THCA | 4 (proxy) | Re-use Paper 1 outputs (panel z + HM450 β + R17 zone). No new compute needed. |
| GSE112202 | 5 (redifferentiation) | In plan. |
| GSE184362 (Pu 2021 sc) | 2 / cell-state | Plan-only (scRNA, heavy). |
| Redifferentiation trials | 5 | Supplement-scan plan. |
| GSE151180 (miRNA) | 2 | Deferred. |
| GSE76039 / GSE286332 / Mun 2025 / Lu 2023 | 4 | Re-use from Paper 1 R17 layers. |

## 4. GSE151179 validation results (this run)

| Comparison | n_test | n_ref | Cohen d | MWU p | AUC | Tier |
|---|---|---|---|---|---|---|
| tumor vs non-neoplastic | 39 | 13 | **−1.77** | **9.5×10⁻⁷** | **0.96** | sanity |
| refractory vs avid (tumor only) | 35 | 4 | −0.24 | 0.49 | 0.61 | 2 |
| no uptake vs uptake (tumor only) | 20 | 19 | −0.27 | 0.68 | 0.54 | 2 |
| persistence vs remission (tumor only) | 35 | 4 | −0.24 | 0.49 | 0.61 | 3 |
| Kruskal by lesion class (4 driver groups) | 4 groups (BRAFV600E 15, Fusion 9, pTERT 3, WT 11) | — | — | 0.44 | — | 2 |

**Honest summary:** The panel cleanly separates tumor from non-neoplastic thyroid (sanity passed, large effect). Within tumors, label imbalance (avid n=4) makes the refractory-vs-avid binary test underpowered, but the direction matches expectation (refractory has lower panel z). The strength of the manuscript will come from **multi-cohort integration** (Boucai + Mu + GSE299988) rather than from GSE151179 alone.

## 5. Figure 1 (Nature Medicine style)

`results/figures/figure1_multigate_rai_model.png` and `.pdf` — four-panel composite:

- **a. Clinical unmet need**: DTC → I-131 → RAI-avid remission vs RAI-refractory persistent/metastatic; 10-yr survival drop > 90% → 10–14% in distant RAIR-DTC.
- **b. Multi-gate biology**: blood → SLC5A5/NIS (gate 1) → PAX8/NKX2-1/FOXE1/TSHR lineage program (gate 2) → TG/TPO/DUOX1/2 organification (gate 3) → colloid retention + radiation killing (gate 4). Footer: "NIS alone is insufficient — coordinated thyroid differentiation across all four gates is required."
- **c. Three molecular states**: RAI-avid differentiated (z = +1.6, all gates open, blue) → molecular gray zone (z = 0, partial gates, gray-purple) → RAI-refractory dedifferentiated (z = −1.6, collapsed gates, red). BRAF / RAS / TERT / TP53 / fusion shown as modifier badges, not the main axis.
- **d. Evidence ladder**: Discovery (TCGA + integrated, Tier 4 proxy) → Panel definition (8 genes) → Label-anchored validation (GSE151179 / GSE299988 / Boucai 2023 / Mu 2024, Tier 1/2 anchor) → Clinical stratification (3 groups).

Style: muted biomedical palette (blue / red / gray-purple / beige), thin lines, sans-serif labels, no 3D effects.

## 6. Manuscript framing locked in

Working title: **"A Thyroid Differentiation-Silencing Axis Identifies Radioiodine-Refractory Thyroid Cancer Across Independent Cohorts"** (see `docs/manuscript_strategy.md` for alternates).

**Tier-respecting wording rule** (project-wide):
- Never write "predicts RAI response" against Tier 4 data.
- Use "associated with" / "stratifies" / "captures a differentiation-linked RAI failure axis" / "consistent with RAI refractoriness".
- Every figure caption tags its dataset tier.

## 7. Reviewer-risk defense already drafted

10 entries in `docs/reviewer_risk_register.md`, including the predictable ones:
- R1 "TCGA has no RAI label" → strict Tier-4 framing.
- R2 "Eight genes are cherry-picked" → LOGO + random-panel permutation + TDS-16 comparison built into `05_validate_rai_labels.py`.
- R3 "RAI refractoriness is not binary" → Mu 2024 4-class gray zone framing.
- R4 "BRAF already explains this" → driver-stratified analyses standard.
- R8 "Sample sizes are small" → cross-cohort meta-analysis + effect-size emphasis.

## 8. Round 2 additions (2026-05-21 follow-up)

The four "다 고고" follow-ups completed:

### 8.1 GSE299988 pipeline run (Tier 2 supportive) — DONE

Verified citation: GSE299988 = "Gene Expression Analysis of Papillary Thyroid Carcinoma with Lymph Node Metastasis and Radioactive Iodine Refractive". Platform: Agilent SurePrint G3 V3 (GPL21185).

Cohort = **5 LN-negative RAI-avid PTC + 5 LN-positive non-RAI-avid PTC + 4 adjacent normal** = 14 total. All 8 panel genes mapped (GPL21185 has clean GENE_SYMBOL column).

Results:

| Comparison | n_test | n_ref | Cohen d | MWU p | AUC | Direction |
|---|---|---|---|---|---|---|
| tumor vs non-neoplastic | 10 | 4 | **−1.61** | **0.002** | **1.00** | ✅ as expected |
| refractive vs avid (tumor only) | 5 | 5 | **+1.13** | 0.15 | 0.20 | ⚠️ reversed |

**Honest caveat**: GSE299988's 5-vs-5 binary went *opposite* to the GSE151179 direction. The dataset selects LN-positive vs LN-negative as the "non-avid vs avid" axis — likely a selection confound rather than a refutation of the panel. Reported as caveat in `notebooks/02_gse299988_validation.ipynb` and Figure 3b. Cannot be used as a clean Tier-2 supportive replication.

### 8.2 Boucai 2023 supplement scan — DONE

PMC fetch confirmed:
- **PMID 36780190 · PMC10106408** · *Clin Cancer Res* 2023.
- Data: "available on request" — no GEO / SRA / dbGaP accession.
- Cohort: **n=8 ER vs n=16 NR** (RECIST v1.1, 1:2 matched on histology + stage).
- Signatures: 16-gene canonical TDS + **64-gene enhanced TDS (eTDS)** + 71-gene BRAF-RAS Score + ERK transcriptional output.
- Supp Table S4 (64-gene eTDS components) is directly downloadable from the PMC supplement → high-yield first move (allows external cross-check of our 8-gene panel against the broader iodide axis without raw data).
- Email draft locked in `docs/data_request_boucai.md`.

### 8.3 Mu 2024 supplement scan — DONE

PMC fetch confirmed:
- **PMID 38060243 · PMC11031230** · *J Clin Endocrinol Metab* 2024 Apr 19;109(5):1231-1240.
- NGDC accession: **HRA004166** (controlled access).
- 220 patients enrolled, **214 analyzed**; 4 RAI uptake patterns with verified counts:
  - **I-RAIR n = 80** (initially RAI refractory)
  - **C-RAIA n = 48** (continually RAI avid)
  - **G-RAIR n = 19** (gradually RAI refractory)
  - **P-RAIR n = 10** (partly RAI refractory)
- Driver frequencies (per PMC main text):
  - BRAF V600E 61.1 % of mutated I-RAIR / 50.7 % TERT promoter in I-RAIR.
  - RAS more frequent in I-RAIA.
  - TP53 associated with I-RAIR.
  - Late-hit (TERT/TP53/PIK3CA) **50 % I-RAIR vs 27 % I-RAIA**.
- **Gray zone = P-RAIR + G-RAIR = n = 29 (14 % of 214)** — usable for Figure 4 framing without per-patient access.
- Email + DAC application draft locked in `docs/data_request_mu_hra004166.md`.

### 8.4 Figure 2 + Figure 3 built — DONE

- **Figure 2 (discovery context)**: 6 sub-panels from Paper 1 R17 outputs (TCGA bulk RNA, Lu 2023 sc, Mun 2025 proteome, TCGA HM450 β, per-zone 8-gene profile, TERT × zone). `results/figures/figure2_discovery_context.png` (718 KB) + `.pdf`.
- **Figure 3 (label-anchored validation)**: a — GSE151179 boxplot (normal/avid/refractory); b — GSE299988 boxplot with honest reversed-direction caveat; c — Mu 2024 4-class pictorial summary with gray-zone bracket. `results/figures/figure3_label_anchored_validation.png` (image rendered) + `.pdf`.

## 9. Next priority moves

1. Fetch Boucai 2023 **Supp Table S4** (64-gene eTDS) → compute overlap with our 8-gene panel; report concordance.
2. Email Boucai + submit Mu HRA004166 DAC application (drafts ready, mailbox + ORCID need to be inserted).
3. **GSE112202 redifferentiation pipeline run** (Tier 5 direction-of-effect).
4. Author-keyboard sections (Hook / Disc 3.1 / Limit / Q9 / Cover ¶1) — voice-protected per marathon mode.

## 9. Honest limitations of this initial build

- No Tier-1 (RECIST) anchor yet — Boucai 2023 access pending.
- Mu 2024 4-class gray-zone data not yet pulled.
- GSE151179 refractory-vs-avid binary is underpowered (avid n=4).
- The 8-gene panel inherits Paper 1's selection — for this manuscript we should report a re-derivation null (random-panel permutation) and a LOGO sensitivity within the RAI-labeled datasets, not only within TCGA.

## 10. Workspace integrity

- `data → /data2/rai_atlas` symlink verified; no large outputs land on root.
- `requirements.txt` lists all libs; `00_setup_environment.py` confirms env OK.
- All scripts are idempotent and accept `<ACCESSION>` arg.
- Voice-protected manuscript sections (per marathon mode rule) are **not** generated in this workspace.
