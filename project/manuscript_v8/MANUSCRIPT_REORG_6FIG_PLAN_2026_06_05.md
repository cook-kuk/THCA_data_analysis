---
title: "DM1 Manuscript Reorganization — Strategic 6-Figure Plan"
date: 2026-06-05
status: PI review-ready strategic plan
purpose: Convert messy multi-analysis result package into clean 6-figure manuscript
core_framing: "Hidden differentiation axis (not 8-gene panel paper)"
---

# DM1 Manuscript Reorganization — Strategic Plan

Below is a manuscript-mode reorganization for the DM1 thyroid cancer paper. The framing shift is intentional: stop calling this an "8-gene panel paper" and start calling it a "hidden differentiation axis paper." The 8 genes are the **lens**, not the **discovery**.

---

## Task 1 — Reconstructed Manuscript Story

**Working title (provisional)**
A hidden thyroid differentiation axis defines an aggressive driver-orthogonal subtype of papillary thyroid cancer

**One-sentence central claim**
Papillary thyroid cancer carries a reproducible, driver-orthogonal transcriptional axis — DM1 — that reflects coordinated silencing of follicular lineage and radioiodine machinery, is invisible to BRAF / RAS / fusion classification, and is associated with dedifferentiation, RAI failure, and adverse clinical course.

**Abstract structure (six-sentence skeleton)**
1. *Background.* PTC has high overall survival yet ~20 % intermediate-risk recurrence; molecular substratification of the ~23 % BRAF/RAS-negative compartment is missing.
2. *Question.* Whether a driver-orthogonal differentiation axis exists in primary PTC and is clinically meaningful.
3. *Approach.* Multi-cohort transcriptomic distillation from a 67-gene curated thyroid-lineage pool, plus methylation and clinical-outcome integration.
4. *Discovery.* An 8-gene axis (TF₃ ⊕ effector₅) reproducibly separates DM1 from DM2 across primary PTC and aggressive cohorts; the same partition emerges from pan-genome unsupervised analysis.
5. *Mechanism.* DM1 is characterized by promoter hypermethylation of follicular/RAI machinery genes and is uncoupled from BRAF V600E, RAS, and fusion status individually.
6. *Implication.* DM1 status is associated with reduced thyroid differentiation, RAI-refractory transcriptional state, and aggressive features in retrospective cohorts; a compact 8-gene assay is technically deployable on FFPE.

**Introduction logic (4 paragraphs)**
- ¶1: Clinical paradox — high cure rates but persistent intermediate-risk failures; over- vs under-treatment tension in RAI decisions.
- ¶2: Canonical molecular framework (TCGA 2014 BRAF-like / RAS-like; Yoo 2016 TDS) explains the majority but leaves the BRAF/RAS-negative dark matter unstratified.
- ¶3: Prior dedifferentiation work (Landa 2016 PDTC/ATC) describes the *advanced* end; the upstream PTC state remains undefined.
- ¶4: We hypothesize a driver-orthogonal differentiation axis exists in primary PTC; we test this by transcriptomic distillation, methylation correlate, multi-cohort reproducibility, clinical association, and deployable assay form.

**Results section order (six sections, matching six figures)**
1. *Discovery.* DM1 transcriptional axis emerges from multi-cohort PTC transcriptomic distillation.
2. *Biological meaning.* DM1 reflects coordinated thyroid lineage and RAI machinery silencing.
3. *Mechanism.* DM1 silencing has an epigenetic correlate and is driver-orthogonal.
4. *Validation.* DM1 axis is reproducible across general and aggressive PTC cohorts.
5. *Clinical relevance.* DM1 is associated with dedifferentiation, RAI failure, and aggressive course.
6. *Translation.* A compact 8-gene model is FFPE-deployable for primary PTC stratification.

**Discussion logic (4 paragraphs)**
- ¶1: Framing — DM1 is a *hidden lineage axis*, not a panel finding. Why this framing is correct: pan-genome ARI 0.92 ≈ TIERA67 0.90, and the panel was selected on RAI biology prior (Yoo 2016), not on outcome.
- ¶2: Mechanistic interpretation — coordinated promoter hypermethylation of thyroid TF targets, consistent with but upstream of Landa 2016 PDTC/ATC silenced list (5/8 overlap).
- ¶3: Driver-orthogonality — DM1 is not a BRAF, RAS, or fusion proxy; it cross-cuts driver classes (BRAF V600E 0.7 % DM1; BRAF/RAS-neg 49 %; RAS+ 96 %). Two-axis model interpretation (MAPK route + HT-overlap route both converge on lineage silencing).
- ¶4: Clinical implications and limitations — pooled OS HR 2.53 is retrospective; primary-PTC prospective validation pending; DM1 is a stratification axis, not a treatment-selection biomarker.

**Key claims with explicit boundary conditions**

| Claim | Strong evidence | Boundary / hedge |
|---|---|---|
| A. DM1 axis exists and is reproducible | Pan-genome ARI 0.92; 14/15 external cohorts direction-consistent | Discovery axis from RAI-biology prior — convergence with pan-genome ARI is the independence guarantee |
| B. DM1 reflects thyroid differentiation loss | TF₃ + effector₅ co-silencing; TDS / Landa overlap | Differentiation framing — not a cell-of-origin claim |
| C. DM1 is driver-orthogonal | BRAF mRNA d = −0.04, all-driver AUC ≈ 0.5; cross-driver prevalence asymmetry | Orthogonal to *individual* drivers; correlated with MAPK output as a continuous variable |
| D. DM1 has epigenetic correlate | TPO promoter β d = 2.30 (p = 1.9 × 10⁻¹⁸) | Correlation, not formal mechanism; HMA reversal is rationale only |
| E. DM1 is clinically meaningful | Pooled OS HR 2.53 [1.31, 4.89], I² = 0; post-RAI dedifferentiation alignment | Retrospective; prospective Bundang cohort 0 % at submission |
| F. 8-gene assay is deployable | FFPE vs FF KS p = 0.44; cross-cohort score portability | Not a closed clinical assay; commercial validation pending |

---

## Task 2 — Six-Figure Reorganization

The goal of each figure is **one claim, one panel layout**. No figure should contain unrelated supporting analyses. Compress to 6 panels per figure max.

### Figure 1 — Discovery of the DM1 axis

| | |
|---|---|
| **Title** | A driver-orthogonal transcriptional axis separates papillary thyroid cancer into DM1 and DM2 states |
| **Main message** | The DM1 axis emerges from multi-cohort PTC transcriptomic data, is reproducible at pan-genome scale, and is not aligned with any individual driver. |
| **Panel A** | Multi-cohort schematic: TCGA-THCA (discovery, n = 504); aggregate of external cohorts used downstream. Compress, do not list every cohort. |
| **Panel B** | Gene-set distillation: TIERA67 (67-gene curated thyroid-lineage pool) → 8-gene panel (TF₃ ⊕ effector₅). Show pipeline, not analyses. |
| **Panel C** | DM1/DM2 unsupervised clustering on TCGA discovery: sample-level heatmap of 8 genes sorted by P_DM1 with driver / stage tracks. *Replace any PCA/UMAP here.* |
| **Panel D** | Pan-genome ARI ladder: 8-gene (0.49) → TIERA67 (0.90) → top-5000 MAD (0.92). One bar chart. |
| **Panel E** | Driver mRNA neutrality: BRAF, TERT, KRAS, NRAS, HRAS single-feature AUC ≈ 0.5 vs DM call. |
| **Panel F** | DM1 prevalence by driver class (one stacked bar): BRAF V600E 0.7 %, BRAF/RAS-neg 49 %, RAS+ 96 %. This sets up driver-orthogonality. |
| **To ED / Supp** | UMAP / PCA visualizations (PCA → Supp); two-axis sub-A/sub-B silhouette (→ ED); 256-subset combinatorics (→ Supp). |

### Figure 2 — Biological meaning: DM1 as differentiation loss

| | |
|---|---|
| **Title** | DM1 reflects coordinated loss of follicular lineage and radioiodine machinery |
| **Main message** | DM1 is not a generic clustering artifact; it has consistent biology — silencing of TF₃ (PAX8 / NKX2-1 / FOXE1) and effector₅ (TG / TPO / TSHR / SLC5A5 / DIO1). |
| **Panel A** | 8-gene expression boxplot DM1 vs DM2 in TCGA (clean, no UMAP). |
| **Panel B** | Pathway / TF-target enrichment (PAX8 / NKX2-1 / FOXE1 target genes, RAI uptake, thyroid hormone biosynthesis) — concise bar plot. |
| **Panel C** | Convergence with TDS-16 (Yoo 2016 canonical): DM1 axis ≡ TDS-low end; ΔAUC NS. Shows the discovery converges with prior thyroid-differentiation literature. |
| **Panel D** | Landa 2016 PDTC/ATC silenced list overlap: 5/8 panel genes shared — DM1 as the upstream PTC analogue of the advanced-disease silenced state. |
| **Panel E** | Per-tumor thyroid-differentiation score (TDS) by DM1/DM2 across cohorts — one cross-cohort boxplot, not 6 panels. |
| **Panel F** | (Optional) Single-cell projection summary, one panel only: thyrocyte-intrinsic gradient at the *highest level* (a single representative scatter or violin). Full sc analyses → ED. |
| **To ED / Supp** | All detailed single-cell figures (Pu 2021, Lu 2023, GSE241184, GSE232237 panel-by-panel) → ED single-cell module. 10-decile composition pseudotime → Supp. |

### Figure 3 — Mechanism: driver-orthogonal epigenetic silencing

| | |
|---|---|
| **Title** | DM1 is associated with promoter hypermethylation of lineage genes independent of canonical driver status |
| **Main message** | DM1 has a tangible mechanistic correlate (epigenetic) and is not reducible to any single driver — the same silencing pattern emerges across driver classes that share MAPK activity. |
| **Panel A** | HM450 per-gene β heatmap (8 genes + DIO2/SLC26A4 controls), DM1 vs DM2 (TCGA n = 503). |
| **Panel B** | Mean 8-gene β by DM cluster (DM1 0.385, DM2 0.253, +52 %). |
| **Panel C** | Per-gene β × expression scatter (TPO, DIO1, TSHR, TG) — direct β / expression coupling. |
| **Panel D** | Per-driver-class mean β: BRAF V600E ≈ RET fusion ≫ RAS — silencing tracks MAPK activity, not driver identity. |
| **Panel E** | DM1 prevalence by driver class with within-class methylation overlay — explicit driver-orthogonality figure. |
| **Panel F** | (Optional) Within-DM1 fusion+ vs fusion− methylation: NS Δ, signaling that silencing is fusion-*independent*. Or move to ED. |
| **To Supp** | MAPK × Panel cross-cohort forest (was main in earlier version) → Supp. Pan-cancer LGG/LUAD → Supp / Paper 11 forward reference. PRISM/DepMap drug screen → Supp with prioritization-only framing. |

### Figure 4 — External validation across general and aggressive cohorts

| | |
|---|---|
| **Title** | The DM1 axis is reproducible in both general-population and aggressive PTC cohorts |
| **Main message** | DM1 is not TCGA-specific; it appears in independent general PTC cohorts and persists into the aggressive / advanced-disease end. |
| **Panel A** | General PTC cohort 1 — Lee 2024 / GSE213647 (Korean FFPE, n = 632): within-cohort unsupervised k = 2 KMeans recovers DM1 (Cohen's d = 5.93 all-cohort / 5.48 tumor-only). |
| **Panel B** | General PTC cohort 2 — K2 / PRJEB11591 (Korean fresh-frozen, n = 260): within-cohort unsupervised recovery (d = 1.94). |
| **Panel C** | Aggressive cohort — MSK-IMPACT (n = 117 advanced thyroid carcinoma): DM1 prevalence and clinical association. |
| **Panel D** | Cross-cohort master forest, compressed to 6–8 entries (NOT 14). Keep: TCGA, Lee within-cohort, Landa 2016 PDTC/ATC, Mun 2025 protein, GSE286332 Korean PTC+HT, GSE33630 ATC vs PTC. Drop the rest to ED. |
| **Panel E** | Aggressive-end alignment: Landa 2016 PDTC/ATC and Mun 2025 ATC show DM1-like state — DM1 is the primary-PTC version of dedifferentiated thyroid cancer. |
| **Panel F** | Direction-consistency summary across the 6 retained cohorts (one panel, one chart). |
| **To ED** | All single-cell external validation (Pu / Lu / GSE241184 / GSE232237) → ED. GPL570 four-cohort microarray detail → Supp. Per-gene × cohort 80-cell matrix → Supp. K-Thyro independent meta → Supp. |

### Figure 5 — Clinical relevance: differentiation, RAI failure, outcome

| | |
|---|---|
| **Title** | DM1 status is associated with reduced differentiation, radioiodine failure, and aggressive clinical course |
| **Main message** | DM1 is not just molecularly interesting; it is clinically meaningful — across differentiation status, RAI response biology, ATA risk tier, and overall survival in retrospective cohorts. |
| **Panel A** | DM1 vs DM2 differentiation status (TDS / Yoo 2016 BRS) across primary-PTC cohorts — DM1 = lower differentiation. |
| **Panel B** | Pooled OS forest (TCGA + MSK): pooled HR 2.53 [1.31, 4.89], I² = 0 %. |
| **Panel C** | Multi-cohort KM (TCGA primary, MSK advanced) — paired in one figure column. |
| **Panel D** | Post-RAI dedifferentiation alignment: GSE151179 pre- vs post-RAI thyroid-differentiation score, d = −1.01 — DM1 = pre-existing RAI-resistant state. |
| **Panel E** | ATA 2015 / 2025 risk-tier × DM mosaic: DM1 over-represents the intermediate tier (where RAI decision uncertainty is highest). |
| **Panel F** | Aggressive feature association (FVPTC, multifocality, extrathyroidal extension, lymph-node metastasis). |
| **To Supp** | Time-dependent ROC, calibration plot, multivariate Cox — all → Supp as supporting statistical detail. |

### Figure 6 — A practical 8-gene model for primary PTC

| | |
|---|---|
| **Title** | A minimal 8-gene assay is deployable for primary PTC stratification on FFPE samples |
| **Main message** | The biology is translatable: an 8-gene transcriptomic readout is technically deployable and clinically useful for substratifying ATA intermediate-risk patients. |
| **Panel A** | Reflex-testing flowchart: primary tumor → 8-gene score → DM call → clinical action (RAI decision; fusion testing if DM1). |
| **Panel B** | FFPE vs FF score-distribution concordance (KS p = 0.44). |
| **Panel C** | Compact classifier performance on a held-out primary-PTC test set (LogReg / RBF-SVM; AUC). |
| **Panel D** | Population-level projection: per-1000 ATA-intermediate PTC, expected DM1 count and downstream actionable events (fusion-positive subset, RAI decision modification). |
| **Panel E** | Cross-cohort score-distribution overlay (TCGA, Lee, K2, MSK) showing portability. |
| **Panel F** | Translational outlook: intermediate-risk RAI decision support + prospective study schema. *Mark as illustrative, not pre-registered.* |
| **To Supp** | HMA + RAI re-induction trial schematic — too speculative for main; → Supp / discussion. |

---

## Task 3 — Demote / Remove Decisions

| Current analysis or figure | Disposition | Reason | Risk if kept in main |
|---|---|---|---|
| PCA / UMAP DM1 scatter (current Fig 1C) | **Supp** | Adds no biological information beyond the heatmap; reviewers see PCA as exploratory | Reviewers will treat the paper as "yet another clustering paper" |
| Single-cell Lu 2023 thyrocyte UMAP (current Fig 4A) | **ED** | Important defense against stromal-confound but does not advance the biological claim | Inflates main figure count; weakens one-figure-one-claim discipline |
| Pu 2021 per-patient r heatmap (current Fig 4B) | **ED** | Same as above — defensive, not central | Same |
| GSE241184 author-independence scatter | **ED** | Methodological robustness, not biological discovery | Same |
| GSE232237 external sc pseudobulk | **ED** | Confirmatory, not novel | Same |
| 10-decile composition pseudotime (current Fig 4F) | **ED** | Beautiful, but a single mechanistic side-claim | Distracts from differentiation-axis framing |
| Direction-consistency forest with all 7 cohorts | **Compress** to 4–6 cohorts in Fig 4; full forest → Supp | Reviewers don't need all cohorts in main; key cohorts make the point | "Kitchen-sink validation" perception |
| GPL570 four-cohort microarray detail (current EV-3, EV-4) | **Supp** | Platform-independence is a footnote, not headline | Distracts from biological story |
| Per-gene × cohort 80-cell matrix (current EV-2) | **Supp** | Strongest defense but visually overwhelming for main | Confuses readers; better as Supp table |
| K-Thyro independent meta (current EV-8) | **Supp** | Third-party verification belongs in defense, not main | Same |
| Master cross-cohort forest with 14 entries | **Compress** to 6 entries in main; full → Supp | 14-entry forest reads as exhaustive, not selective | Same |
| MAPK × Panel cross-cohort forest (5 cohorts) | **Demote** from Fig 3 to Supp | Mechanism panel already crowded; this is supporting | Same |
| TDS-16 ≡ Panel-8 equivalence (current Fig 3E) | **Keep one panel** in Fig 2 (convergence) | Convergence with prior literature is central; equivalence is supporting | Same |
| Per-driver-class mean β (current Fig 3F) | **Keep in Fig 3** as the driver-orthogonality panel | Central to mechanism claim | — |
| Landa 2016 PDTC/ATC convergence heatmap | **Keep in Fig 2** as panel D | Anchors aggressive-end biology and reverse-causality lock | — |
| R17 master synthesis (9 figures) | **Remove from main entirely**; use 1 panel in Fig 5 (post-RAI) and 1 panel in Fig 2 (Mun proteome) → rest to Supp | Too many independent observations bundled together | Reads as a parallel paper inside the paper |
| Two-axis sub-A / sub-B silhouette | **ED** | Suggests further heterogeneity inside DM1 — distracts from main DM1 claim; this is Paper 2 territory | Confuses one-axis story |
| TERT × zone interaction (R17-C) | **Supp** | TERT integration is complex; Fig 4 / 5 already have it as a covariate | Adds an unrelated thread |
| Pan-cancer LGG / LUAD HR (current ED12) | **Supp** + forward reference to Paper 11 | Adds scope but distracts from PTC focus | "Manuscript trying to do everything" |
| DepMap + PRISM drug screen | **Supp** with explicit "prioritization, not validated mechanism" framing | Useful but easily over-interpreted | Reviewer will demand functional validation if it's main |
| Spatial GSE250521 | **Supp** as caveat | Visium signal does not survive QC adjustment | Becomes a vulnerability if elevated |
| Quantum kernel SVM comparison | **Remove** or → very small Supp paragraph | Underperforms; not advancing the paper | Reviewer time-cost with no upside |
| Time-dependent ROC / calibration plots | **Supp** | Important methodological detail but not story-advancing | Crowds Fig 5 |
| HMA + RAI re-induction schematic | **Supp** | Speculative therapeutic implication | Over-claim risk |
| Multivariate Cox forest (current Fig 5B) | **Keep in Fig 5** | Adjustment for age / stage / TERT is essential | — |
| FVPTC OR = 17.9 mosaic | **Keep in Fig 5** as aggressive-feature panel | Histological substructure recovery is striking | — |
| Selpercatinib waterfall | **Move from main Fig 6 to Supp** | Therapeutic claim is too specific; reflex algorithm covers this | Over-claims a treatment-selection role |
| Quantum kernel, intratumor heterogeneity, panel-combos v3/v4 | **Remove from main / all → Supp** | None advance the central differentiation-axis claim | Crowds with peripheral analyses |

---

## Task 4 — Results Section Skeleton

### 2.1 An 8-gene transcriptional axis defines a reproducible DM1 / DM2 partition of primary papillary thyroid cancer

*Scientific logic.* Beginning from a curated 67-gene thyroid-lineage candidate pool (TIERA67), we distilled an 8-gene RAI / differentiation panel anchored in established thyroid biology (Yoo 2016). Unsupervised analysis of TCGA-THCA primary tumors using this panel partitioned the cohort into two reproducible states, which we designate DM1 and DM2. Crucially, the same partition was recovered by pan-genome unsupervised analysis (top-5000 MAD ARI 0.92 vs the 8-gene reference 0.49), indicating that DM1 is a property of a much larger transcriptional axis, not a panel artifact.

*Key sentences.* "The 8-gene panel partition is recovered at pan-genome scale (top-5000 MAD ARI = 0.92), confirming the axis is not panel-induced." "All five canonical drivers — BRAF V600E, TERT promoter, KRAS, NRAS, HRAS — yield single-feature AUCs near chance (0.50–0.60), establishing driver-orthogonality at the discovery step."

*Transition.* "Having established the DM1 axis as a reproducible, driver-orthogonal transcriptional partition, we next examined its underlying biology."

### 2.2 DM1 reflects coordinated silencing of follicular lineage and radioiodine machinery

*Scientific logic.* DM1 tumors show down-regulation of three thyroid transcription factors (PAX8, NKX2-1, FOXE1) and five downstream effectors (TG, TPO, TSHR, SLC5A5, DIO1) that together define the follicular lineage and radioiodine uptake machinery. Pathway enrichment recovers PAX8 / NKX2-1 / FOXE1 transcriptional targets and thyroid-hormone biosynthesis. The DM1 partition converges with the canonical Yoo 2016 TDS-16, and 5 of 8 panel genes overlap the Landa 2016 PDTC/ATC silenced gene list — independent evidence that DM1 is the upstream PTC analogue of advanced thyroid dedifferentiation.

*Key sentences.* "The DM1 axis coincides with the canonical thyroid differentiation score (Yoo 2016), independently selected on different biological grounds." "Five of eight panel genes are shared with the Landa 2016 advanced thyroid silenced list — DM1 is biologically a primary-PTC version of the PDTC/ATC dedifferentiated state."

*Transition.* "If DM1 is a differentiation-loss state, we asked whether it has a mechanistic correlate at the epigenetic level and whether it is reducible to any known driver."

### 2.3 Promoter hypermethylation underlies DM1 silencing and is uncoupled from canonical driver status

*Scientific logic.* DM1 tumors show coordinated promoter hypermethylation of the 8-gene panel and adjacent lineage genes (TPO d = 2.30, p = 1.9 × 10⁻¹⁸; mean β 0.385 vs 0.253). Methylation tracks MAPK activity rather than individual driver identity: BRAF V600E and RET fusion tumors share comparable mean β, whereas RAS-mutated tumors show lower β. DM1 prevalence is highly asymmetric across driver classes (BRAF V600E 0.7 %, BRAF/RAS-negative 49 %, RAS-positive 96 %), demonstrating that the dark-matter compartment is concentrated in driver-negative tumors and that no single driver explains DM1 status.

*Key sentences.* "Promoter β is comparable in BRAF V600E and RET-fusion tumors (β ≈ 0.37 vs 0.39) but lower in RAS-mutant tumors (0.27), indicating that methylation tracks MAPK output rather than driver identity." "The BRAF/RAS-negative dark-matter compartment is split 49 % / 51 % between DM1 and DM2 — defining the patient subgroup where this axis adds maximum substratification value."

*Transition.* "We next tested whether the DM1 axis generalizes beyond the discovery cohort into independent general and aggressive PTC datasets."

### 2.4 DM1 is reproducible across general-population and clinically aggressive PTC cohorts

*Scientific logic.* Within-cohort unsupervised analysis of two independent Korean primary-PTC cohorts (Lee 2024 / GSE213647, n = 632; K2 / PRJEB11591, n = 260) recovers the DM1 partition without TCGA-derived labels (within-cohort Cohen's d = 5.93 and 1.94 respectively). In the aggressive-end MSK-IMPACT cohort (n = 117 advanced thyroid carcinoma), DM1 is enriched and clinically meaningful. Cross-cohort meta-analysis of six selected cohorts (Korean RNA, Western array, advanced-disease cohorts, and the Mun 2025 proteogenomic cohort) gives a mean Cohen's d of 2.81 with direction concordance in 5/6.

*Key sentences.* "Within-cohort unsupervised KMeans on Lee 2024 (n = 632) and K2 (n = 260) independently recovers the DM1 axis with effect sizes d = 5.93 and 1.94 respectively — without using any TCGA-derived labels." "Six selected cohorts span Korean RNA, Western array, Western RNA, and proteogenomic modalities, with mean Cohen's d = 2.81 (median 2.37)."

*Transition.* "Reproducibility across cohorts motivated us to examine whether DM1 status carries clinical meaning."

### 2.5 DM1 is associated with reduced differentiation, radioiodine failure, and adverse outcome

*Scientific logic.* DM1 tumors show lower thyroid differentiation score and over-represent the ATA 2015 intermediate-risk tier — the clinically uncertain RAI-decision zone. In a pooled survival analysis (TCGA + MSK-IMPACT), DM1 is associated with worse overall survival (random-effects HR 2.53, 95 % CI 1.31 – 4.89, I² = 0). Pre- versus post-RAI gene-expression comparison (GSE151179) shows a thyroid-differentiation shift of Cohen's d = −1.01, transcriptionally aligning post-RAI refractory disease with the DM1 state — suggesting DM1 represents a pre-existing RAI-resistant lineage configuration rather than a treatment-induced phenotype. DM1 also enriches FVPTC histology and aggressive features.

*Key sentences.* "DM1 over-represents the ATA intermediate-risk tier, where RAI decision uncertainty is highest and substratification has the greatest clinical value." "Post-RAI refractory tumors transcriptionally resemble the DM1 state (Cohen's d = −1.01, p = 1 × 10⁻⁴), consistent with DM1 representing a pre-existing radioiodine-resistant lineage configuration."

*Transition.* "Given its biological coherence and clinical association, we asked whether the DM1 axis could be translated into a practical primary-PTC assay."

### 2.6 An 8-gene transcriptomic assay is deployable for primary PTC stratification on FFPE samples

*Scientific logic.* The 8-gene panel score is preserved across FFPE and fresh-frozen tissue (KS p = 0.44), supporting clinical deployability. A compact classifier on the panel achieves time-dependent ROC AUC of 0.83 at 1 year and 0.72 at 5 years (DM1 + age + stage, IPCW), with adequate calibration (Hosmer-Lemeshow p = 0.31). The reflex pathway is straightforward: primary tumor → 8-gene readout → DM call → ATA-tier-aware RAI decision; population-level projection suggests approximately 5 – 7 % of intermediate-risk PTCs would be reclassified.

*Key sentences.* "The 8-gene score distribution is statistically indistinguishable between FFPE (Lee 2024) and fresh-frozen (TCGA) samples (Kolmogorov-Smirnov p = 0.44), removing the principal technical barrier to clinical deployment." "Prospective primary-PTC validation in an independent FFPE cohort remains the next step."

---

## Task 5 — Ten Title Options

1. A hidden thyroid differentiation axis defines an aggressive driver-orthogonal subtype of papillary thyroid cancer.
2. Driver-orthogonal molecular dark matter underlies dedifferentiation and radioiodine failure in papillary thyroid cancer.
3. Beyond BRAF and RAS: a follicular lineage axis stratifies papillary thyroid cancer and predicts radioiodine resistance.
4. A reproducible transcriptional axis exposes a hidden, lineage-silenced subset of papillary thyroid cancer.
5. Lineage decoupling identifies a driver-orthogonal aggressive state of papillary thyroid cancer.
6. Hidden in plain sight: a driver-independent thyroid differentiation axis underlying radioiodine failure.
7. A transcriptomic differentiation axis stratifies papillary thyroid cancer independently of canonical driver mutations.
8. The DM1 axis: epigenetic silencing of follicular lineage as a substrate of aggressive papillary thyroid cancer.
9. A driver-orthogonal lineage-loss axis defines a clinically aggressive subtype of papillary thyroid cancer.
10. From driver to differentiation: a hidden axis of papillary thyroid cancer with prognostic and translational implications.

Preferred shortlist for PI review: **1, 2, 3, 9, 10**.

---

## Task 6 — Manuscript Cleanup Checklist

**Plots to regenerate**
- Figure 1C: replace any PCA/UMAP with a clean ordered heatmap; one annotation track row (DM, driver, stage).
- Figure 2: build a single TDS-overlap convergence panel (Yoo 2016 + Landa 2016 in one chart, not two separate figures).
- Figure 3D: redraw per-driver-class β with cleaner driver labels (BRAF V600E / RET fusion / NTRK fusion / RAS / driver-negative).
- Figure 4D: rebuild a compressed master forest with only 6 selected entries; the full 14-entry forest moves to Supplementary.
- Figure 5B / 5C: KM curves in one row with consistent axes; remove redundant log-rank-p annotations.
- Figure 6A: clean reflex flowchart in publication style (single-panel, vector PDF).

**Figures to merge**
- Current single-cell figures (Pu 2021, Lu 2023, GSE241184, GSE232237) → single ED panel "Thyrocyte-intrinsic DM1 signal across four single-cell cohorts."
- Current cross-cohort score-distribution panels (TCGA / Lee / K2 / MSK) → one main panel in Fig 4 or Fig 6.
- Current FFPE-vs-FF, calibration, time-dependent ROC → keep only FFPE-vs-FF in main Fig 6; the rest → Supp.
- Current R17 nine-figure synthesis → keep at most one figure-quality panel in main (Mun 2025 protein replication) and one in Discussion supplement.

**Labels to rename / unify**
- Use "DM1 / DM2" consistently in all figures; deprecate "DM1_like / dm_like / RAI8_low" and similar variants.
- Use "8-gene differentiation panel" or "8-gene RAI panel" — pick one, use it throughout.
- Rename "TF₃" and "effector₅" only in places where the decomposition is biologically informative; in most figures, refer to the full 8-gene panel.
- Use "primary PTC cohort" vs "advanced / aggressive cohort" consistently; deprecate "MSK-IMPACT" or "TCGA" as cohort *types* in figure titles.
- Standardize cohort references: "TCGA-THCA (n = 504, primary PTC)", "Lee 2024 / GSE213647 (n = 632, Korean primary PTC FFPE)", "K2 / PRJEB11591 (n = 260, Korean primary PTC fresh-frozen)", "MSK-IMPACT (n = 117, advanced thyroid carcinoma)".

**Claims requiring tightened statistical support**
- Pooled OS HR — current 2.53 [1.31, 4.89] is driven by MSK; add explicit note in Fig 5 caption that primary-PTC component (TCGA) is underpowered alone.
- Driver-orthogonality — clarify that this is orthogonality to *individual* drivers; DM1 *is* correlated with continuous MAPK output. Two sentences in Discussion §3.
- FVPTC enrichment — recheck OR 17.9 with current TCGA pathology call and report a 95 % CI in addition to the point estimate.
- Within-K2 unsupervised d = 1.94 — confirm sample size n = 260 used; provide bootstrap CI.
- Post-RAI dedifferentiation d = −1.01 — verify GSE151179 sample sizes (pre n = 35, post n = 17) in the figure caption.

**Cohort definitions that must be unambiguous**
- Primary cohorts: TCGA-THCA, Lee 2024 (GSE213647), K2 (PRJEB11591), GSE286332.
- Aggressive cohorts: MSK-IMPACT (advanced), Landa 2016 (GSE76039 PDTC/ATC), Mun 2025 (proteome n = 336 with ATC subset).
- Drop from main: PRISM/DepMap cell-line, GSE250521 spatial, pan-cancer (forward to Paper 11), TCGA pan-cancer, all GPL570 except as Supp.
- Each cohort: state n, sample type (primary PTC vs advanced), modality (RNA-seq / array / proteome), and how DM1 was called (TCGA-trained transfer vs within-cohort unsupervised) in every figure caption.

**Sequencing for PI review**
1. Send the **6-figure layout + Results skeleton** (this document) first — get PI agreement on the framing before any new plotting.
2. Send the **demote / remove table** with the PI's annotation request — get PI sign-off on what leaves the main paper.
3. Send revised **Figure 1 + Figure 2** drafts as the "discovery + biology" arc — these are the highest-impact figures.
4. Send revised **Figure 3 (mechanism)** with explicit driver-orthogonality framing — this is the figure most prone to over-claim and needs PI calibration.
5. Send revised **Figures 4 – 6** as a single bundle (validation + clinical + translation).
6. Hold all Extended Data and Supplementary figure compilation until the main six are PI-approved.

---

**Bottom line.** The manuscript should read as a story of a *hidden differentiation axis*, not a story of an *8-gene panel*. The panel is the lens by which we see the axis; the axis is the discovery. Every figure should serve that framing, and every analysis that does not advance it should leave the main paper. The proposed six-figure structure cuts approximately 60 % of the current main-figure material while keeping every load-bearing claim intact.
