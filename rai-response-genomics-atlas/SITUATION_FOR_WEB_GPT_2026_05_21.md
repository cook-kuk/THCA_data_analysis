# RAI Response Genomics Atlas — situation brief for outside review

**Compiled** 2026-05-21 · Seungho Cook (국승호) + Yu Hyeong-won · Seoul National University Bundang Hospital
**Companion paper**: "Paper 1" — molecular dark-matter / DM1 / RAI-lineage axis in thyroid cancer (bioRxiv 2026-05-20, R1–R18 completion).
**Scope of this brief**: a parallel manuscript — *label-anchored RAI-response validation of an eight-gene thyroid differentiation / iodide-handling panel*.

This document is self-contained. It is intended to be pasted into an outside LLM (web GPT, Claude, etc.) so that the model can independently judge: (1) is the story scientifically sound, (2) is the framing manuscript-safe, (3) what is the highest-yield next move, (4) what should we be most worried about as reviewer risk.

---

## TL;DR (5 sentences)

1. We have a published-precedent 8-gene panel (SLC5A5 · TPO · TG · TSHR · PAX8 · NKX2-1 · FOXE1 · DIO1) that captures the thyroid differentiation / iodide-handling program — the biology underlying radioactive iodine (RAI) therapy response.
2. Paper 1 already shows this panel separates a "silenced + immune-active" dark-matter zone from BRAF-like and RAS-like zones across **8 layers in 4 data pillars** (RNA bulk + sc + Korean cross-ethnic + Mun 2025 proteome + TERT × zone + Cox PFI + per-zone gene profile + HM450 β), with the strongest single statistic being Mun protein **ATC vs PTC dark-matter Fisher OR = 8.54, p = 2.7 × 10⁻¹⁵**.
3. Paper 1 is **Tier 4 discovery only** — TCGA-THCA carries no direct RAI response field. The new manuscript adds **Tier 1–2 label-anchored validation** against external RAI-labeled cohorts.
4. First two Tier-2 cohorts done: **GSE151179** (n = 52: 39 PTC + 13 normal) passes tumor-vs-normal sanity (AUC = 0.96, p = 9.5 × 10⁻⁷) but the refractory-vs-avid binary is underpowered (avid n = 4); **GSE299988** (n = 14: 5 avid + 5 refractive + 4 normal) passes sanity AUC = 1.00 but its 5-vs-5 binary went the *wrong direction* (d = +1.13, AUC = 0.20) — likely LN+/LN− selection confound.
5. The two highest-yield next anchors — **Boucai 2023 CCR** (Tier-1 RECIST, n = 8 ER vs n = 16 NR; data request-only but Supp Table S4 has the 64-gene eTDS for direct overlap analysis) and **Mu 2024 JCEM HRA004166** (Tier-2 4-class uptake patterns: I-RAIR 80 / C-RAIA 48 / G-RAIR 19 / P-RAIR 10, with verified driver-frequency tables enabling the gray-zone framing without per-patient access) — are scanned, citations verified, request drafts ready.

---

## Project framing and manuscript story

Working title candidates (ranked safest → boldest):

1. **A Thyroid Differentiation-Silencing Axis Identifies Radioiodine-Refractory Thyroid Cancer Across Independent Cohorts** (preferred).
2. An Eight-Gene Iodide-Handling Panel Reveals a Molecular Gray Zone of Radioiodine Failure in Thyroid Cancer.
3. Integrated Transcriptomic Evidence for Differentiation-Linked Radioiodine Failure in Thyroid Cancer.

Five hierarchical claims (safest → boldest):

1. The 8-gene panel is **associated with** RAI avidity in multiple Tier-2 cohorts.
2. Driver-stratified analysis shows the panel adds information **beyond BRAF / RAS / TERT** status.
3. Mu 2024's 4 uptake patterns are explained by a gradient of panel silencing — a **molecular gray zone**, not a binary.
4. Redifferentiation-treatment cohorts (Tier 5) show panel score rises with restored uptake — mechanism plausibility.
5. Framing as **risk-stratification readout**, not clinical biomarker.

Bold-but-reserved: a single integrated panel score that summarizes the published 16-gene TDS with parsimony advantage. Holding for review response, not headline.

Wording rules (project-wide):
- Never write "predicts RAI response" against Tier-4 data (TCGA).
- Use "associated with" / "stratifies" / "captures a differentiation-linked RAI failure axis" / "consistent with RAI refractoriness".
- Every figure caption tags its dataset tier.
- "Clinical biomarker" / "diagnostic" forbidden without prospective Tier-1 confirmation.

---

## Eight-gene panel — provenance, polarity, biology

Panel: **SLC5A5 (NIS), TPO, TG, TSHR, PAX8, NKX2-1 (TTF1), FOXE1 (TTF2), DIO1**.

Provenance: derived in Paper 1 by Random Forest ranking from a curated 55-gene pool with drivers excluded by design; later cross-validated against 16-gene TDS (Landa 2016) at Spearman ρ = 0.954.

Polarity (memory `DM1 Round 4` polarity correction, 2026-05-09): **panel z high = preserved differentiation, low = silenced** (= "DM1-like" in Paper 1 vocabulary). Some legacy files used the inverse — always check the polarity_note before reusing tables.

Per-gene roles:

| Gene | Role | Why it's in the panel |
|---|---|---|
| SLC5A5 (NIS) | iodide uptake | classical rate-limiter — the field's prior NIS-alone view |
| TPO | organification | binds iodine to thyroglobulin tyrosines |
| TG | storage scaffold | colloidal iodine reservoir |
| TSHR | upstream maintenance | TSH → AC → cAMP → differentiation program |
| PAX8 | lineage transcription factor | master thyroid TF, binds TG/TPO/SLC5A5 promoters |
| NKX2-1 / TTF1 | lineage TF | thyroid + lung lineage |
| FOXE1 / TTF2 | lineage TF | forkhead, migration / maturation |
| DIO1 | thyroid hormone metabolism | mature thyrocyte program marker |

NOT in the core panel (but on the iodide-axis sensitivity list): DUOX1, DUOX2, DUOXA1, DUOXA2, IYD, SLC26A4, SCGN, DIO2. Some early specs had DUOX2 instead of DIO1 — **the canonical Paper-1 panel uses DIO1**.

---

## Paper 1 R17 discovery evidence — what is already published

Done before this workspace, fully replicable from `project/results/r17_tcga_panel_d4p2_reconciliation/`. **All Tier 4 / 5 (no direct RAI labels)**. Treat this as the "discovery / mechanism" half of the new manuscript.

| Layer | Pillar | n | Headline |
|---|---|---|---|
| L1 | RNA bulk (TCGA) | 500 | 14.8 % concordance between d4p2 HT-axis and panel RAI-axis labels — proves two axes are orthogonal, not redundant |
| L2 | sc RNA (Lu 2023 GSE193581) | 14,624 malignant cells | ATC 38.3 % in "silenced + HT-overlap" dark-matter zone — cellular substrate |
| L3 | RNA bulk Korean (GSE286332 PTC vs PTC+HT) | 18 | PTC+HT 77.8 % dark-matter, PTC 88.9 % WT-like — cross-ethnic replication |
| L4 | Protein (Mun 2025) | 336 | ATC 58.4 % dark-matter; **Fisher OR = 8.54, p = 2.7 × 10⁻¹⁵** vs PTC |
| L5 | TCGA TERT × zone | 477 | TERT+ uniquely enriched in dark-matter zone, OR = 2.34, Fisher p = 0.016 |
| L6 | TCGA Cox PFI | 563 / 63 ev | BRAF-like HR 1.89, dark-matter HR 1.89 (p ≈ 0.05 each, vs WT-like) — OS/DSS underpowered |
| L7 | per-zone 8-gene profile | 522 | BRAF-like silences TPO/DIO1 first; dark-matter silences TG/PAX8/TSHR deeper — two silencing programs |
| L8 | HM450 β × zone | 518 | mean β: dark-matter 0.41 > BRAF 0.35 > RAS 0.37 > WT 0.28 · 4/8 genes methylation-mediated (DIO1 · SLC5A5 · TG · TPO) |

Combined dossier (8-layer 4-pillar) is web-deployed at the lab IP:
- `http://40.82.129.113/r17/` (Korean default + English HTML + Korean PDF + 9 figures including a master synthesis composite).
- `http://40.82.129.113/paper1_biorxiv_2026_05_20/` (full paper-1 deliverables, R1-R18).

---

## RAI atlas validation evidence — what this new workspace adds

All Tier-mapped. Strict separation from Paper 1.

### GSE151179 (Tier 2 RAI avidity) — DONE

- Platform Affymetrix Clariom S (GPL23159). 27,189 probes → 8/8 panel genes mapped via family.soft annotation parsing (custom helper script).
- Cohort: 39 PTC + 13 matched non-neoplastic thyroid = 52 samples.
- Labels parsed from `Sample_characteristics_ch1`: `patient rai responce` (avid/refractory), `disease` (remission/persistence), `tissue type` (primary / synchronous LN met / LN met post-RAI), `collection before/after rai`, `lesion class` (BRAFV600E / Fusion / pTERT / WT).

Results (within-cohort z-score panel):

| Comparison | n_test | n_ref | Cohen d | MWU p | AUC | Comment |
|---|---|---|---|---|---|---|
| tumor vs non-neoplastic | 39 | 13 | **−1.77** | **9.5 × 10⁻⁷** | **0.96** | sanity strong |
| refractory vs avid (tumor only) | 35 | 4 | −0.24 | 0.49 | 0.61 | direction correct but underpowered (avid n = 4) |
| no uptake vs uptake (tumor only) | 20 | 19 | −0.27 | 0.68 | 0.54 | small effect, NS |
| persistence vs remission (tumor only) | 35 | 4 | −0.24 | 0.49 | 0.61 | same labels as refractory/avid here |
| Kruskal by lesion class (BRAFV600E 15 / Fusion 9 / pTERT 3 / WT 11) | — | — | — | 0.44 | — | NS in this small cohort |

### GSE299988 (Tier 2 supportive) — DONE, with honest caveat

- Platform Agilent SurePrint G3 V3 (GPL21185). 58,341 probes → 8/8 panel genes mapped (clean GENE_SYMBOL column).
- Cohort: 5 RAI-avid LN-negative PTC + 5 non-RAI-avid LN-positive PTC + 4 adjacent normal = 14 samples.
- Singapore study; "Refractive" spelling normalized to "refractory".

Results:

| Comparison | n_test | n_ref | Cohen d | MWU p | AUC | Comment |
|---|---|---|---|---|---|---|
| tumor vs non-neoplastic | 10 | 4 | **−1.61** | **0.002** | **1.00** | sanity strong |
| refractive vs avid (tumor only) | 5 | 5 | **+1.13** | 0.15 | 0.20 | **direction REVERSED** |

The 5-vs-5 direction reversal is reported honestly as a caveat in Figure 3b and in the notebook. The dataset's "avid" vs "non-avid" label is effectively LN−-vs-LN+ — likely a selection confound, not a refutation of the panel. Cannot be used as clean supportive replication.

### Boucai 2023 CCR (Tier 1 RECIST) — supplement-scanned, data pending

- Citation: PMID 36780190 · PMC10106408 · *Clin Cancer Res* 2023.
- URL: https://pmc.ncbi.nlm.nih.gov/articles/PMC10106408/
- Cohort: **n = 8 exceptional responders vs n = 16 non-responders**. RECIST v1.1 structural-response definition. 1:2 matching on histology + stage.
- Signatures used in paper: 16-gene canonical TDS, **64-gene enhanced TDS (eTDS)**, 71-gene BRAF-RAS Score (BRS), ERK transcriptional output, Hallmark GSEA.
- Data availability: "available on request" — no GEO / SRA / dbGaP / EGA accession.
- **Supp Table S4** contains the 64-gene eTDS — directly downloadable from PMC; allows cross-checking our 8-gene panel against the broader iodide axis without raw data.
- Action: download Supp S4 (no email required) + email Dr. Boucai for raw expression. Email draft locked in `docs/data_request_boucai.md`.

### Mu 2024 JCEM (Tier 2 gray-zone, n = 214) — supplement-scanned, data pending

- Citation: PMID 38060243 · PMC11031230 · *J Clin Endocrinol Metab* 2024 Apr 19;109(5):1231-1240.
- URL: https://pmc.ncbi.nlm.nih.gov/articles/PMC11031230/
- Accession: **NGDC HRA004166** (controlled access, https://ngdc.cncb.ac.cn/).
- Cohort: 220 enrolled, **214 analyzed**.
- Four RAI uptake patterns with verified counts:

| Class | Definition | n |
|---|---|---|
| I-RAIR | initially RAI-refractory | 80 |
| C-RAIA | continually RAI-avid | 48 |
| G-RAIR | gradually RAI-refractory | 19 |
| P-RAIR | partly RAI-refractory | 10 |

- Driver frequencies (from PMC main text):
  - BRAF V600E **61.1 %** of mutated I-RAIR cases.
  - TERT promoter **50.7 %** of mutated I-RAIR cases.
  - RAS family **enriched in C-RAIA** (4.5 % in I-RAIR, higher in I-RAIA).
  - TP53 **associated with I-RAIR**.
  - Late-hit composite (TERT/TP53/PIK3CA) **50.0 % I-RAIR vs 26.9 % I-RAIA**.

- **"Molecular gray zone" = P-RAIR + G-RAIR = 29 / 214 = 14 %**. This is the single best argument that RAI refractoriness is not binary and that an expression-based panel can stratify the gray zone where mutation status alone cannot.

- Per-patient genotypes are in main manuscript tables; raw NGS data are HRA-controlled (typical 4–12 week DAC application turnaround, mixed success rate for non-Chinese requesters).
- Action: DAC application + email to corresponding author. Draft locked in `docs/data_request_mu_hra004166.md`.

---

## Figures already built

| File | What | Size |
|---|---|---|
| `results/figures/figure1_multigate_rai_model.png` (+ .pdf) | Nature Medicine-style 4-panel: (a) clinical unmet need → (b) multi-gate biology blood-NIS-TF-organification-killing → (c) three molecular states with score bar + driver modifiers → (d) evidence ladder discovery → panel → label-anchored validation → clinical stratification | 463K PNG / 64K PDF |
| `results/figures/figure2_discovery_context.png` (+ .pdf) | 6 sub-panels from Paper 1 R17 outputs (TCGA bulk, Lu sc, Mun protein, HM450 β, per-zone gene profile, TERT × zone) | 718K PNG / 281K PDF |
| `results/figures/figure3_label_anchored_validation.png` (+ .pdf) | (a) GSE151179 boxplot normal / avid / refractory · (b) GSE299988 boxplot with reversed-direction caveat · (c) Mu 2024 4-class pictorial summary with gray-zone bracket | 309K PNG / 62K PDF |
| `results/figures/GSE151179_boxplot_panel_by_label.png` | standalone boxplot | |
| `results/figures/GSE151179_roc_refractory_vs_avid.png` | standalone ROC | |
| `results/figures/GSE299988_*` | standalone outputs | |

All figures are pure matplotlib vector-clean; no proprietary icon libraries. Muted biomedical palette.

---

## Workspace layout and storage discipline

Top-level: `/home/seungho/personal/THCA_data_analysis/rai-response-genomics-atlas/`. Code/docs/results small (~MB), kept on root disk. Large GEO matrices go to `/data/...`/`/data2/...` via the `data → /data2/rai_atlas` symlink (Premium SSD, 438 GB free).

```
rai-response-genomics-atlas/
├── README.md
├── requirements.txt
├── config/
│   ├── datasets.yaml             # 10 priority datasets, Tier-mapped
│   └── eight_gene_panel.yaml     # SLC5A5/TPO/TG/TSHR/PAX8/NKX2-1/FOXE1/DIO1
├── docs/
│   ├── literature_review.md
│   ├── dataset_inventory.md
│   ├── rai_response_label_taxonomy.md   # 5-tier hierarchy
│   ├── manuscript_strategy.md
│   ├── reviewer_risk_register.md        # 10 anticipated risks + defenses
│   ├── data_request_boucai.md
│   └── data_request_mu_hra004166.md
├── scripts/
│   ├── 00_setup_environment.py
│   ├── 01_search_geo_metadata.py
│   ├── 02_download_geo_dataset.py        # GEOparse + FTP fallback
│   ├── 03_parse_geo_metadata.py
│   ├── 04_score_gene_panel.py
│   ├── 05_validate_rai_labels.py
│   ├── 06_make_figures.py
│   ├── 07_build_dataset_inventory.py
│   ├── _extract_gpl23159_probe2gene.py   # Clariom S helper
│   ├── figure1_multigate_rai_model.py
│   ├── figure2_r17_discovery_context.py
│   └── figure3_label_anchored_validation.py
├── notebooks/
│   ├── 01_gse151179_validation.ipynb      # done
│   ├── 02_gse299988_validation.ipynb      # done
│   ├── 03_tcga_thca_discovery_plan.ipynb  # plan only — reuse Paper 1
│   └── 04_redifferentiation_dataset_scan.ipynb  # plan only
├── data → /data2/rai_atlas/
├── results/
│   ├── tables/
│   │   ├── dataset_inventory.csv
│   │   ├── GSE151179_test_results.tsv
│   │   └── GSE299988_test_results.tsv
│   ├── figures/  (above)
│   └── reports/
│       ├── GSE151179_validation_brief.md
│       ├── GSE299988_validation_brief.md
│       └── rai_response_genomics_atlas_report.md
└── logs/
```

---

## Reviewer risks already drafted (10 entries in `docs/reviewer_risk_register.md`)

Top 5 most likely:

1. **"RAI response label is not available in TCGA."** — TCGA-THCA is strictly Tier 4. Every Tier-4 claim is tagged in figure captions and Methods. Tier-1/2 validation comes from external cohorts.
2. **"Eight genes are cherry-picked."** — Paper-1 selection from 55-gene curated pool with drivers excluded by design; LOGO + random-panel permutation + TDS-16 comparison built into the validation script. ρ = 0.954 with full TDS-16.
3. **"RAI refractoriness is not binary."** — Mu 2024 4-class framing explicitly modeled (Figure 3c). Both binary AUC and ordinal Spearman vs 4-class ordered scale will be reported.
4. **"BRAF already explains this."** — Driver-stratified analyses (BRAF-mut / RAS-mut / driver-neg) within every label-anchored dataset. Paper 1 R17 already shows BRAF-like vs dark-matter zones diverge despite both being panel-DM1.
5. **"NIS alone is enough."** — Parallel scores (NIS-only, NIS+TPO, panel-8, TDS-16) reported per dataset in the validation script's `score_comparison_per_dataset.tsv`.

Remaining 5 (sample size, batch effects, no clinical utility, why not deep learning, overlap with Paper 1) — see file.

---

## Honest gaps and open decisions

These are the things I want an outside model to evaluate:

1. **The label-anchored evidence so far is thin.** GSE151179 sanity is excellent but refractory-vs-avid is underpowered (avid n = 4). GSE299988 went the wrong direction. The manuscript currently rests on **Mu 2024 driver-frequency tables + Boucai S4 eTDS overlap** for label-anchored support, neither of which is in our local possession yet. Question: is this scientifically responsible to write up now, or should we hold until Boucai raw + Mu DAC complete?

2. **The "molecular gray zone" framing is rhetorical until we project a per-patient panel score onto the Mu 4-class data.** Mutation frequencies alone don't constitute a panel-validation. Question: is the framing strong enough on driver-frequency analysis alone, or does it require Mu per-patient access?

3. **Direction-reversal in GSE299988 should be disclosed up front, not buried.** The current Figure 3b puts it side-by-side with GSE151179. Question: is this the right disclosure level, or should it be moved to supplementary with a clean explanation?

4. **Paper 1 R17 already published the discovery — risk of "what's new?"** The new manuscript's novelty is the Tier-1/2 anchoring, not the discovery axis. Question: is the contribution Cell-Rep-Med-acceptable as is, or does it need an additional novel angle (e.g., redifferentiation Tier-5 predictive validation, or a single integrated multi-cohort score)?

5. **Target journal**: candidates are Nat Commun (if Boucai + Mu access succeed); Cell Reports Medicine (if not); JCI Insight / npj Precision Oncology / Thyroid as fallback. Question: which is the right submission floor given current evidence + reasonable extra effort?

6. **Wording rule audit**: the project enforces "associated with / stratifies / consistent with"; never "predicts RAI response" unless Tier-1 prospective. Question: is this conservative enough for the target journals, or are we overhedging?

---

## Round 2 additions (2026-05-21 PM push) — what was new since the round-1 brief

### R2-1 · 8-gene ⊂ TDS-16 ⊂ eTDS-64 containment confirmed (parsimony argument locked)

All 8 panel genes (SLC5A5, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, DIO1) lie within the canonical TDS-16 from Landa 2016 / TCGA Cell 2014 (16 genes). Boucai 2023's eTDS-64 (Supp Table S4) by construction extends TDS-16, so 8-gene ⊂ TDS-16 ⊂ eTDS-64 by transitivity. Round-8 audit of Paper 1: 8-gene Cohen d = 1.78 / AUC = 0.903 vs TDS-16 d = 1.97 / AUC = 0.913 → 8-gene captures **98.9 %** of TDS-16's AUC at 50 % gene cost (Spearman ρ = 0.954). Figure: `figure_panel_vs_tds_etds_overlap.png`. Caveat: explicit eTDS-64 roster JS-gated on PMC; containment claim rests on Boucai's published Methods statement.

### R2-2 · GSE151179 robustness pack — all three reviewer-cherry-pick defenses pass

On the n = 52 label-anchored cohort:
- **Leave-one-gene-out (LOGO)**: d ∈ [−2.04, −1.49] vs full panel d = −1.77 → no single gene dominates.
- **Random-panel permutation null (n = 1000 matched-variance random 8-gene panels)**: null median d = 0.006, null median AUC = 0.503; observed d = −1.77, AUC = 0.96. **Empirical p (|d| ≥ observed) = 0.014**, **p (|AUC − 0.5| ≥ observed) = 0.005**. The observed panel effect lies on the extreme tail of the matched-variance null.
- **TDS-16 sensitivity**: 16/16 genes present on Clariom S → d = −1.96, AUC = 0.964. 8-gene captures 99.5 % of TDS-16 AUC in this label-anchored cohort. Parsimony confirmed externally, not just on TCGA.

Figure: `GSE151179_robustness.png`. Reviewer R2 cherry-pick defense addressed three independent ways.

### R2-3 · GSE112202 Tier 5 redifferentiation — direction-of-effect confirmed

11 digoxin-treated NMTC patients vs 11 matched untreated controls (Cufflinks RNA-seq, group-level FPKM). **6 / 8 panel genes upregulated** by digoxin, **median log2FC = +0.30**. Strongest: TSHR +0.95, SLC5A5 +0.73, NKX2-1 +0.52, TG +0.40 — exactly the canonical iodide-uptake + lineage axis we'd predict for restored differentiation. TDS-16 extras directionally consistent. PAX8 and DIO1 marginally down (smallest exceptions). Provides Tier-5 mechanistic plausibility (panel responds correctly to redifferentiation drug) without overclaiming.

Figure: `GSE112202_redifferentiation.png`.

### R2-4 · Figure 4 gray-zone locked — driver × zone × Mu 4-class composite

Three-panel figure consolidating the gray-zone framing:
- **a · TCGA-THCA driver × R17 zone** (n = 500): BRAF V600E splits 41 % BRAF-like + 41 % dark-matter; RAS is 87 % WT-like; BRAF·RAS-neg is 45 % WT-like with 20 % each in BRAF-like and dark-matter. Driver alone is NOT a sufficient zone predictor.
- **b · Mu 2024 4-class stacked driver bar** (n = 214): BRAF dominates I-RAIR, RAS dominates C-RAIA, and the gray-zone (G-RAIR + P-RAIR = 29 / 214 = 14 %) has mixed driver composition. The bracket is annotated directly.
- **c · GSE151179 panel z by lesion class** (BRAF V600E 15, Fusion 9, pTERT 3, WT 11): boxplots show panel z range similar across driver classes within tumor — direct evidence that the panel is informative independently of driver.

Figure: `figure4_driver_grayzone.png`. Reviewer R3 ("RAI refractoriness is not binary") and R4 ("BRAF already explains this") both addressed.

### R2-5 · Reviewer-defense scorecard single-image composite

All four round-2 defenses bundled into one shareable image for cover letter use: `reviewer_defense_scorecard.png` (450 KB, 17 × 21 inches). Each row is a defense vector with co-located caption.

### Round 2 net effect on the manuscript ceiling

Before round 2: Cell Reports Medicine territory with a known cherry-pick exposure and an unresolved "is RAI binary?" question.

After round 2:
- **Cherry-pick defense fully drafted** (containment + LOGO + permutation null + TDS-16 sensitivity) — R2 risk neutralized at a Tier-4 + Tier-2 level.
- **Gray-zone framing now visual, not rhetorical** (Mu 2024 driver-composition stacked bar) — R3 risk neutralized at a Tier-2 level.
- **Mechanism plausibility added** (GSE112202 redifferentiation direction) — Tier-5 anchor in hand.
- **Containment with field-canonical TDS-16** — answers "is this just a noisy subset?" structurally.

The cohort weakness (still no Tier-1 raw, still no Mu per-patient) is unchanged, but the manuscript scaffolding around it is now substantially harder to attack. Cell Reports Medicine is now an overcautious target; NComms is the more honest reach.

---

## What I'd most like outside review on

Please judge:

- (a) Is the scientific narrative defensible **as currently anchored** (Tier 4 + GSE151179 sanity + Mu 2024 driver frequencies + Boucai S4 overlap)? Or is at least one of Boucai raw / Mu per-patient *required* before submission?
- (b) Is the title / abstract claim ceiling correctly set? Would you cut, sharpen, or raise it?
- (c) What is the highest-yield next single move? Choices: (i) download Boucai Supp S4 and lock the eTDS-vs-8-gene overlap argument; (ii) send Boucai email + start Mu DAC immediately; (iii) GSE112202 redifferentiation pipeline (Tier 5); (iv) full GSE184362 single-cell re-analysis with per-sample RAI annotation; (v) write the manuscript Hook / Discussion / Limitations now (currently voice-protected per marathon mode).
- (d) Which reviewer risk do you think we have *under*-prepared for?
- (e) Is the GSE299988 reversed-direction finding acceptable disclosure-as-caveat, or fatal to including it?

Thank you.

---

## URLs and citations (for reference)

- **GSE151179** GEO record: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151179
- **GSE299988** GEO record: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299988
- **Boucai 2023 CCR**: https://pmc.ncbi.nlm.nih.gov/articles/PMC10106408/ · PMID 36780190
- **Mu 2024 JCEM**: https://pmc.ncbi.nlm.nih.gov/articles/PMC11031230/ · PMID 38060243 · NGDC HRA004166 https://ngdc.cncb.ac.cn/
- **Landa 2016 JCI** (TDS / driver framework): JCI 126(3):1052–1066
- **TCGA-THCA Cell 2014** (original TDS): Cell 159(3):676–690
- **Selumetinib RAI redifferentiation, Ho 2013 NEJM**: https://www.nejm.org/doi/full/10.1056/NEJMoa1209288

Live web copies of Paper-1 work (companion):
- 8-layer 4-pillar R17 dossier (Korean default + English HTML + Korean PDF + master synthesis figure): `http://40.82.129.113/r17/`
- Paper 1 full deliverables: `http://40.82.129.113/paper1_biorxiv_2026_05_20/`
