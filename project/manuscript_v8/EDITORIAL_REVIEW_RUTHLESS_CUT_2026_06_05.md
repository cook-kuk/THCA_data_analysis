---
title: "Editorial review — DM1 manuscript rescue plan (ruthless cut)"
date: 2026-06-05
role: Senior Nature Communications editor perspective
purpose: Remove analyses, not add them. Convert 38-figure archive → 5-figure NC paper.
core_principle: When uncertain between keep vs remove, REMOVE.
companion_doc: MANUSCRIPT_REORG_6FIG_PLAN_2026_06_05.md (strategic 6-figure plan)
---

# Editorial review — DM1 manuscript rescue plan

I am writing as a senior reviewer. My job is to tell you what to remove. Every recommendation below is biased toward cutting, because that is what this manuscript needs.

---

## STEP 1 — Inventory and ruthless classification

I classify each item against the only three questions that matter:
- Q1 Does DM1 exist?
- Q2 What biological state does it represent?
- Q3 Why does it matter clinically?

If a figure does not directly answer one of these, it leaves the main paper.

### Discovery layer

| Item | Verdict | Editor's reasoning |
|---|---|---|
| 8-gene heatmap sorted by P_DM1 (annotation tracks) | **Essential** | The single visual that proves DM1/DM2 exist as a structure. |
| Driver-neutrality bar (BRAF / TERT / RAS AUC ≈ 0.5) | **Essential** | One panel, sets up driver-orthogonality before any cohort claim. |
| Pan-genome ARI ladder | **Essential** | The one defense against "your panel created the cluster." |
| PCA scatter (any version) | **Remove** | Reviewers do not need PCA to believe a heatmap. PCA reads as exploratory. |
| UMAP scatter (any version) | **Remove from main; minimal supp** | The heatmap shows the same partition with more information. |
| TIERA67 → 8-gene Sankey | **Useful, supplementary** | Process diagram, not a result. |
| Two-axis sub-A vs sub-B silhouette | **Remove** | This is Paper 2 territory and weakens the one-axis story. |
| 256-subset combinatorics | **Remove** | Defensive, exhausting, no biological information. |
| Quantum-kernel SVM comparison | **Remove** | Underperforms; serves no editorial purpose. |
| Cluster-stability bootstrap (97.4 %) | **One Methods sentence**, not a figure | Belongs in text. |

### Biological meaning layer

| Item | Verdict | Editor's reasoning |
|---|---|---|
| 8-gene expression boxplot DM1 vs DM2 | **Essential** | Direct visualization of the biology. |
| TF / pathway enrichment (PAX8 / NKX2-1 / FOXE1 / RAI uptake) | **Essential** | Anchors the differentiation-loss interpretation. |
| TDS-16 (Yoo 2016) convergence | **Essential** | The single most important "we are not making up biology" panel. |
| Landa 2016 PDTC/ATC overlap heatmap | **Essential** | Connects DM1 to the canonical advanced-disease literature. |
| MAPK output × Panel-8 cross-cohort forest (5 cohorts) | **Remove from main** | Mechanism support, not headline; move to supp. |
| TDS-16 ≡ Panel-8 per-driver-class equivalence | **Supplement** | Defensive only. |
| Cell-composition deconvolution (NNLS / nu-SVR / Ridge / LR-clip) | **Supplement** | Methodological; not biological. |

### Mechanism layer

| Item | Verdict | Editor's reasoning |
|---|---|---|
| HM450 per-gene β heatmap | **Essential** | The mechanism evidence. |
| Mean β bar DM1 vs DM2 | **Essential** | One-number summary the reader will remember. |
| β × expression scatter (4 genes) | **Essential** | Shows the methylation–silencing link. |
| Per-driver-class mean β (BRAF ≈ RET ≫ RAS) | **Essential** | The driver-orthogonality panel. |
| DM1 prevalence by driver class (BRAF 0.7 % / RAS+ 96 %) | **Essential** | Defines the patient subgroup. |
| Within-DM1 fusion+/− methylation | **Supplement** | Fusion-independence is defensive support. |
| Pan-cancer LGG / LUAD HR | **Remove** | Belongs in a different paper. |
| DepMap + PRISM drug screen | **Remove** | Over-interpreted; reviewer will demand functional validation. |
| Spatial GSE250521 | **Remove** | Does not survive QC adjustment; becomes a vulnerability if shown. |

### External validation layer

| Item | Verdict | Editor's reasoning |
|---|---|---|
| TCGA-THCA (discovery) | **Essential** | Origin cohort. |
| Lee 2024 / GSE213647 (n = 632 Korean PTC FFPE) | **Essential** | Largest independent PTC cohort, within-cohort recovery d = 5.93. |
| Landa 2016 / GSE76039 (PDTC + ATC) | **Essential** | Anchors aggressive end. |
| Mun 2025 proteogenomic (n = 336) | **Essential** | The only cross-modality (protein) validation. |
| K2 / PRJEB11591 (n = 260) | **Supplement** | Calibration mismatch becomes a distraction at the main level. |
| MSK-IMPACT advanced | **Keep ONLY in Fig 5** for survival; do not duplicate in validation | One use, not two. |
| GSE286332 (Korean PTC + HT, n = 18) | **Supplement** | n = 18; mention as supporting cohort only. |
| GSE33630 / GSE29265 / GSE65144 / GSE53157 (GPL570 array, n = 205) | **Supplement** | Platform-agnosticism is a footnote, not a headline. |
| K-Thyro independent meta | **Supplement** | Defensive only. |
| Master cross-cohort forest with 14 entries | **Compress to 4** in main; full → supp | The 14-entry forest reads as exhaustive, not selective. |
| Per-gene × cohort × contrast 80-cell matrix | **Supplement** | Strong defense, weak figure aesthetics. |
| External direction-consistency forest | **Supplement** | Replicated by the compressed master forest. |
| Nonoverlap scatter grid (4 cohorts) | **Supplement** | Defensive reverse-causality argument. |
| Korean applicability bar | **Discussion sentence**, not a figure | One number. |

### Single-cell validation layer

| Item | Verdict | Editor's reasoning |
|---|---|---|
| Lu 2023 thyrocyte UMAP | **Extended Data (one panel)** | Sufficient to defend against stromal-confound. |
| Pu 2021 per-patient r heatmap | **Remove or merge** | Adds nothing beyond Lu. |
| GSE241184 author-independence | **Supplement** | Methodological. |
| GSE232237 external pseudobulk | **Supplement** | Confirmatory. |
| Multisite pooled scatter | **Remove** | Redundant. |
| 10-decile composition pseudotime | **Remove from main; supp** | Beautiful but a side-claim. |
| Direction-consistency forest across 7 sc cohorts | **Supplement** | Reviewers cannot evaluate seven sc cohorts in detail. |

### Clinical layer

| Item | Verdict | Editor's reasoning |
|---|---|---|
| Pooled OS forest (TCGA + MSK) HR = 2.53 | **Essential** | The headline survival statistic. |
| TCGA-THCA KM | **Supplement** (or merge into one panel showing pooled curves) | Underpowered alone, dominates ≠ leads. |
| MSK advanced KM | **Supplement** | Same. |
| Multivariate Cox forest (DM1 + age + stage + TERT) | **Essential** | Independence claim. |
| Post-RAI dedifferentiation GSE151179 box | **Essential** | Cleanest "DM1 has clinical relevance" panel. |
| ATA risk-tier × DM mosaic | **Essential** | Defines where the assay would matter. |
| Time-dependent ROC (1y/3y/5y) | **Supplement** | Methodological. |
| 5-year calibration plot | **Supplement** | Methodological. |
| Decision-curve analysis | **Supplement** | Net-benefit framing is reviewer-controversial. |
| Selpercatinib waterfall | **Remove** | Over-claims a treatment-selection role. |
| HMA + RAI re-induction schematic | **Supplement** | Speculative. |
| Prospective trial schema | **Remove from main; Discussion paragraph** | Imaginary trial does not belong in a results figure. |
| FFPE vs FF concordance (KS p = 0.44) | **Essential** | The deployability statement. |

### R17 master synthesis (9 figures)

| Item | Verdict | Editor's reasoning |
|---|---|---|
| R17-A master synthesis | **Remove** | Reads as a parallel paper inside the paper. |
| R17-B two-axis reconciliation | **Discussion sentence** | Concept, not figure. |
| R17-C TERT × differentiation zone | **Supplement** | Side thread. |
| R17-D / E / F zone forest / gene profile / HM450 β | **Remove** | Three separate restatements of the same idea. |
| R17-G Mun 2025 proteome | **Use as Fig 4 panel** | Only R17 figure that earns its keep. |
| R17-H GSE286332 two-axis | **Supplement** | n = 18, supporting evidence only. |
| R17-I sc Lu two-axis | **Remove** | Subsumed by Lu sc panel in ED. |

### Extended Data set (current 15)

Reduce to ~6 ED items:
- ED1 P_DM1 heatmap detail with driver tracks
- ED2 Pan-genome ARI fine ladder
- ED3 Single-cell thyrocyte-intrinsic gradient (one panel, Lu only)
- ED4 Within-DM1 fusion+/− methylation
- ED5 Mun 2025 proteogenomic per-gene effect
- ED6 Time-dep ROC + calibration

Everything else → supplementary tables / text.

---

## STEP 2 — The smallest viable Nature Communications manuscript

Five main figures is the right number. Six only if Translation cannot be accommodated in Discussion. I recommend you target five.

### Figure 1 — Discovery

| | |
|---|---|
| **Title** | A driver-orthogonal transcriptional axis separates papillary thyroid cancer into two reproducible states |
| **One-sentence claim** | An 8-gene readout, anchored in thyroid lineage biology, partitions PTC into DM1 and DM2 — a partition independent of canonical drivers and recovered at pan-genome scale. |
| **Panels** | A: cohort schematic (3 cohorts only). B: 8-gene heatmap sorted by P_DM1, annotation tracks (DM / driver / stage). C: driver-neutrality bar. D: pan-genome ARI ladder. |
| **Delete** | PCA. UMAP. Sankey. Sub-A/sub-B. ARI 6-bar ladder → 4 bars only. |

### Figure 2 — Biological meaning

| | |
|---|---|
| **Title** | DM1 represents coordinated loss of thyroid lineage and radioiodine machinery |
| **One-sentence claim** | DM1 is a differentiation-loss state — recapitulating the canonical thyroid differentiation score and overlapping with the advanced-disease silenced lineage program. |
| **Panels** | A: 8-gene expression boxplot DM1 vs DM2. B: TF target / RAI pathway enrichment. C: TDS-16 (Yoo 2016) convergence. D: Landa 2016 silenced-list overlap. |
| **Delete** | Single-cell here; cell-composition deconvolution; per-cohort TDS panels (one boxplot, not six). |

### Figure 3 — Mechanism

| | |
|---|---|
| **Title** | DM1 silencing has an epigenetic correlate and is uncoupled from individual driver mutations |
| **One-sentence claim** | DM1 tumors show coordinated promoter hypermethylation of lineage genes that tracks MAPK activity rather than driver identity, with no canonical driver alone explaining DM1 status. |
| **Panels** | A: HM450 per-gene β heatmap (8 genes + 2 controls). B: mean β bar DM1 / DM2. C: β × expression scatter (4 representative genes). D: per-driver-class mean β. E: DM1 prevalence by driver class (BRAF V600E 0.7 % / BRAF-RAS-neg 49 % / RAS+ 96 %). |
| **Delete** | MAPK cross-cohort forest, TDS-equivalence, fusion-vs-non-fusion methylation comparison, pan-cancer cameo, drug screens. |

### Figure 4 — External validation

| | |
|---|---|
| **Title** | The DM1 axis is reproducible in independent PTC, advanced-disease, and proteogenomic cohorts |
| **One-sentence claim** | DM1 recovers from an independent Korean primary PTC cohort, persists into advanced thyroid carcinoma, and is detectable at the protein level. |
| **Panels** | A: Lee 2024 within-cohort unsupervised recovery (d = 5.93). B: Landa 2016 PDTC/ATC alignment. C: Mun 2025 proteogenomic d. D: compressed cross-cohort forest (4 entries only — TCGA / Lee / Landa / Mun). |
| **Delete** | K2 detail, GPL570 four-cohort, MSK as separate panel here (it goes to Fig 5), single-cell external cohorts, K-Thyro meta, 80-cell matrix, nonoverlap scatter grid. |

### Figure 5 — Clinical relevance and deployability

| | |
|---|---|
| **Title** | DM1 is associated with adverse outcome and radioiodine failure, and is detectable on FFPE samples |
| **One-sentence claim** | DM1 status is associated with reduced differentiation, post-RAI refractory transcriptional state, and adverse survival in retrospective cohorts — and is technically deployable on FFPE primary-PTC samples. |
| **Panels** | A: pooled OS forest (TCGA + MSK), HR 2.53. B: multivariate Cox forest (DM1 + age + stage). C: post-RAI dedifferentiation (GSE151179) d = −1.01. D: ATA risk-tier × DM mosaic. E: FFPE vs FF concordance (KS p = 0.44). |
| **Delete** | Selpercatinib waterfall, trial schema, decision-curve, HMA schematic, time-dep ROC, calibration, ATA + reflex flowchart, ATA + population projection. Move these to Discussion paragraphs and supplement. |

### (Optional) Figure 6 — Translation

I recommend **not** including Figure 6. Replace it with a single Discussion paragraph titled "Translational outlook" and one supplementary figure (compact-assay schematic + FFPE concordance restated). A separate Figure 6 will not strengthen the manuscript and will invite reviewer demands for prospective validation that you do not yet have.

---

## STEP 3 — Duplicated evidence

| Duplication | What you have now | What to keep | What to remove |
|---|---|---|---|
| **Cohort proliferation** | 19 external cohorts | TCGA + Lee 2024 + Landa 2016 + Mun 2025 (4 cohorts) | K2 (→ supp), GSE286332 (→ supp), GPL570 (×4, → supp), K-Thyro meta (→ supp), Pu 2021 / Lu 2023 / GSE241184 / GSE232237 (single-cell → 1 ED panel), GSE151179 used once in Fig 5 |
| **Dimensionality plots** | PCA, UMAP scatters, 8-gene heatmap | 8-gene heatmap (Fig 1B) | All PCA / UMAP |
| **Cluster proofs** | 6 different "DM1/DM2 separates" panels | One in Fig 1 + ARI ladder | UMAP, sub-A/B, 256-subset, quantum kernel, two-axis silhouette |
| **Survival plots** | TCGA KM, MSK KM, pooled forest, multivariate Cox, time-dep ROC, calibration, decision curve | Pooled forest (Fig 5A) + multivariate Cox forest (Fig 5B) | KMs individually (→ supp), time-dep ROC (→ supp), calibration (→ supp), decision-curve (→ supp) |
| **Single-cell** | Pu 2021, Lu 2023, GSE241184, GSE232237, pseudotime, direction forest | Lu 2023 only, in ED | All other sc analyses |
| **Mechanism panels** | HM450 heatmap, mean β, β-expr scatter, MAPK forest, TDS equivalence, per-driver β, Landa convergence, fusion-vs-non methylation | First 3 + per-driver β + Landa in Fig 2/3 | MAPK forest (→ supp), TDS equivalence (→ supp), fusion comparison (→ supp) |
| **Forest plots** | Master cross-cohort 14-entry, MAPK 5-cohort, direction-consistency 7-cohort | One compressed 4-entry forest in Fig 4 | All other forests (→ supp) |
| **Heatmaps** | Expression heatmap, β heatmap, β × expression heatmaps, per-driver-class β heatmap | One in Fig 1, one in Fig 3 | Per-gene × cohort 80-cell heatmap (→ supp) |
| **R17 synthesis** | 9 figures | None | All 9 (one panel → ED for Mun 2025 only) |

Editor's rule: **if two figures make the same point, the second figure is automatically supplementary.**

---

## STEP 4 — Reviewer risks and minimum defenses

### Figure 1 (Discovery)
- **Risk:** "Your panel created the cluster."
  **Defense:** Pan-genome ARI ladder is in Fig 1D; explicitly state pan-genome top-5000 MAD ARI = 0.92.
- **Risk:** "Why 8 genes specifically?"
  **Defense:** One supplementary panel showing AUC stability with 4, 8, 12, 16-gene variants + LOO sensitivity.
- **Risk:** "PCA would show separation more clearly."
  **Defense:** One supplementary UMAP panel — single image, no text.

### Figure 2 (Biology)
- **Risk:** "TDS overlap is circular because you selected on TDS."
  **Defense:** One Methods sentence + supplementary panel showing the 8-gene panel was distilled from TIERA67 (curated thyroid-lineage pool) before any TDS comparison.
- **Risk:** "Pathway enrichment is hand-picked."
  **Defense:** One supplementary table with full enrichment output across MSigDB Hallmark + KEGG + Reactome.

### Figure 3 (Mechanism)
- **Risk:** "Methylation correlation is not mechanism."
  **Defense:** Acknowledge explicitly in caption; provide supplementary panel showing per-CpG β stratified by clinical variables.
- **Risk:** "Driver-orthogonality is sampling-driven."
  **Defense:** Multivariate Cox in Fig 5B adjusts for drivers; supplementary panel with DM1 prevalence × driver class confidence intervals.

### Figure 4 (Validation)
- **Risk:** "Only one Korean primary cohort?"
  **Defense:** One supplementary panel with K2 (PRJEB11591) within-cohort unsupervised recovery (d = 1.94) acknowledging calibration limitations.
- **Risk:** "Single-cell validation should be in main."
  **Defense:** ED panel with Lu 2023 thyrocyte gradient.
- **Risk:** "Compressed forest is selective."
  **Defense:** Supplementary table with all 14 cross-cohort entries.

### Figure 5 (Clinical)
- **Risk:** "Pooled HR is dominated by advanced disease (MSK)."
  **Defense:** State this explicitly in Fig 5A caption and supplement; primary-PTC component is underpowered alone.
- **Risk:** "ATA tier mosaic is descriptive."
  **Defense:** Multivariate Cox in Fig 5B handles this.
- **Risk:** "Post-RAI sample sizes are small."
  **Defense:** Supplementary table with GSE151179 n breakdown.
- **Risk:** "FFPE concordance only one cohort."
  **Defense:** Acknowledge; commit to prospective validation in Discussion.

Total supplementary figures needed for defense: **6 supplementary figures + ~5 supplementary tables**. No more.

---

## STEP 5 — Final paper structure

```
Figure 1 · Discovery
   The DM1 axis exists and is driver-orthogonal.

Figure 2 · Biological meaning
   DM1 is differentiation loss.

Figure 3 · Mechanism
   DM1 is epigenetically silenced and not reducible to any single driver.

Figure 4 · External validation
   DM1 reproduces in independent PTC, advanced-disease, and protein cohorts.

Figure 5 · Clinical relevance and deployability
   DM1 marks adverse course and post-RAI failure, and is FFPE-deployable.

(No Figure 6.)

Extended Data: 6 panels.
Supplementary: 5 tables + ~6 figures + Methods.
```

That is the entire paper. Anything not on this list is either supplementary or removed.

---

## STEP 6 — Realistic 2-week submission timeline

The goal is clarity, not perfection. Use only existing data.

### Week 1 — Decisions and figure assembly (no new analyses)

**Day 1 (Mon)** — Editorial alignment with PI
- Send PI: title shortlist (preferred: "A driver-orthogonal differentiation axis defines an aggressive subtype of papillary thyroid cancer").
- Send PI: 5-figure layout (this document).
- Lock the framing: "DM1 axis paper," not "8-gene panel paper."
- Approve the **remove list** before any plotting work begins.

**Day 2 (Tue)** — Figure 1 and 2
- Pull existing assets: 8-gene heatmap, driver-neutrality bar, pan-genome ARI bars, TDS convergence, Landa overlap, pathway enrichment.
- Compose Fig 1 and Fig 2 from existing PNGs (do not re-run analyses; use what is in `papers_hub_2026_05_04/assets/paper1_nc/` and `manuscript_v8_nc_main/`).
- Single layout pass each.

**Day 3 (Wed)** — Figure 3 and 4
- Compose Fig 3 from existing HM450 heatmap, mean β bar, β × expression scatter (4-gene grid), per-driver-class β, driver-class prevalence.
- Compose Fig 4: Lee 2024 within-cohort recovery, Landa 2016 alignment, Mun 2025 protein effect, compressed 4-entry forest.

**Day 4 (Thu)** — Figure 5 and Extended Data
- Compose Fig 5: pooled OS forest, multivariate Cox, GSE151179 post-RAI box, ATA mosaic, FFPE vs FF density.
- Assemble Extended Data (6 panels max): P_DM1 heatmap detail, fine ARI ladder, single-cell thyrocyte gradient (Lu only), fusion-vs-non methylation, Mun proteogenomic per-gene, time-dep ROC + calibration.

**Day 5 (Fri)** — Write Results (Section 2)
- Six paragraphs total, one per figure plus one transition pair (Fig 5 absorbs the discovery → biology → mechanism → validation → clinical handoffs).
- Target: ~2,500 words for Results.

**Day 6–7 (Sat–Sun)** — Buffer / co-author review of figures
- Send Figs 1–5 plus ED to all co-authors as a single PDF.
- Collect only structural feedback at this stage, not wording changes.

### Week 2 — Text and finalization

**Day 8 (Mon)** — Introduction
- Four paragraphs: clinical paradox, canonical molecular framework, gap (driver-orthogonal axis missing), preview of findings.
- Target: ~700 words.

**Day 9 (Tue)** — Discussion
- Four paragraphs: differentiation-axis framing, mechanism interpretation, driver-orthogonality and clinical implication, limitations and translational outlook.
- Target: ~1,000 words. No Fig 6; clinical translation is one Discussion paragraph plus one supplementary figure.

**Day 10 (Wed)** — Abstract, cover letter, key resources table
- Six-sentence structured abstract.
- Cover letter (200 words) framing DM1 as a hidden differentiation axis with translational implications, not as a panel paper.

**Day 11 (Thu)** — Methods + Supplement
- STAR Methods compressed: cohort summary table (one), score computation (one paragraph), statistical methods (one paragraph), data/code availability.
- Supplement: 5 tables + 6 figures only.

**Day 12 (Fri)** — Internal scientific review
- One scientific reviewer (co-author or external colleague) reads the full draft cold.
- One editorial reviewer (you) reads only for clarity, redundancy, and over-claim.

**Day 13 (Sat)** — Final edits
- Tighten claims. Add boundary language. Remove any reviewer-attack surfaces still in the text.
- Verify every statistic cited in text matches every figure.

**Day 14 (Sun)** — Submission packet
- Final PDF assembly. ORCID, data deposition, code release tag. Submit.

---

## Editorial summary

You currently have approximately 38 figures across main, extended data, and supplementary positioning. The manuscript needs **5 main figures, 6 extended data figures, and 5 supplementary tables**. The reduction from 38 to ~16 figures (∼ 60 % cut) is what converts this from a result archive into a publishable Nature Communications paper.

The single most important framing change is: stop introducing the manuscript with the panel. Introduce it with the axis. The 8 genes are how you measure DM1. DM1 is what you discovered. Everything else — fusion, methylation, survival, FFPE — is downstream of that one finding.

When you are uncertain about a figure, the answer is **remove it**.
