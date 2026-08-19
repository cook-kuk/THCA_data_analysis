# Reviewer risk register — RAI Response Genomics Atlas

Every entry maps an anticipated reviewer concern to a specific defensive response, the figure/table/data that supports the defense, and the action item.

---

### R1. "RAI response label is not available in TCGA."

**Strength of concern**: high (universal first-line reviewer trap).

**Response**: TCGA-THCA is positioned strictly as the **discovery / molecular-context** cohort (Tier 4 per `rai_response_label_taxonomy.md`). All label-anchored validation is performed against Tier 1 (Boucai 2023 CCR exceptional responders) and Tier 2 (GSE151179 RAI-avidity, GSE299988, Mu 2024 4-class uptake). Every TCGA-derived claim is explicitly tagged "Tier 4 proxy" in figure captions and text.

**Supporting evidence**: explicit Tier mapping in `docs/dataset_inventory.md`; manuscript Methods § "Label tier conventions".

**Action**: enforce wording rule via lint check in scripts (Method `04_score_gene_panel.py` emits warnings if "predict" appears in any plot title).

---

### R2. "Eight genes are cherry-picked."

**Strength**: high (always asked of small panels).

**Response**: The panel is a **published Paper 1 panel** (Cook & Yu 2026) selected from a 55-gene curated pool by Random-Forest ranking with drivers explicitly excluded by design (see Paper 1 memory `v17 8-gene audit 2026-04-29`). We perform:

- Leave-one-gene-out (LOGO) — show no single gene dominates the score (per-gene mean |d| range 0.72–1.33 in Paper 1 DM1 robustness work).
- Random-panel permutation (n=1000 random 8-gene panels from matched expression range) — show observed panel z-score effect at the > 99th percentile of nulls.
- TDS-16 comparison — Spearman ρ=0.954 with full Landa 16-gene TDS (Paper 1 R8 audit), confirming our 8-gene captures 90% of TDS-16 information at 50% gene cost.

**Supporting evidence**: `results/figures/fig_logo.png`, `results/figures/fig_random_panel_null.png`, `results/tables/panel_vs_tds16.tsv`.

**Action**: script `05_validate_rai_labels.py` runs LOGO + random null automatically.

---

### R3. "RAI refractoriness is not binary."

**Strength**: high (a Mu 2024-aware reviewer will press this).

**Response**: We model RAI refractoriness as a **tiered, gray-zone phenotype** by leveraging:

- Mu 2024 JCEM's 4-class uptake patterns (initially refractory / continually avid / gradually refractory / partly refractory).
- GSE151179 sample-type stratification (primary tumor / LN met / distant met with avidity annotation).
- Panel-score continuum (z-score), not binarized.

We report **both** binary AUC and ordinal Spearman correlation with the 4-class scale.

**Supporting evidence**: Figure 4d gray-zone diagram; Methods § "Ordinal modeling of RAI uptake patterns".

**Action**: `05_validate_rai_labels.py` emits both binary (Mann-Whitney + ROC) and ordinal (Kruskal-Wallis + Spearman vs ordered class) tests.

---

### R4. "BRAF already explains this."

**Strength**: high.

**Response**: We perform driver-stratified analyses **within** BRAF-mutant, RAS-mutant, and driver-negative groups separately. Paper 1 R17 shows the panel separates "BRAF-like silenced + no HT" zone from "dark-matter silenced + HT" zone — two BRAF-positive populations with different prognostic and therapeutic implications. Within driver-negative BRAF·RAS-neg subset (the most BRAF-confounded subset removed), the panel retains effect.

**Supporting evidence**: R17 L1 cross-tab + L7 per-zone gene profile; Figure 4a driver × panel-zone heatmap.

**Action**: every panel-vs-label plot is also drawn within each driver stratum (`06_make_figures.py` produces per-driver panel sets).

---

### R5. "NIS alone is enough."

**Strength**: medium (older reviewers; less common in 2025+).

**Response**: We compare:

- panel score using all 8 genes
- NIS (SLC5A5) alone
- NIS + TPO (functional pair)
- TDS-16 (Landa)

across each label-anchored cohort. Hypothesis: the 8-gene panel matches or modestly outperforms NIS-alone with much smaller variance, and is well within TDS-16 performance with 50% gene cost.

**Supporting evidence**: `results/tables/score_comparison_per_dataset.tsv`; Figure 3 supplement.

**Action**: `04_score_gene_panel.py` produces parallel scores for NIS-only, NIS+TPO, panel-8, TDS-16 when source data permits.

---

### R6. "No clinical utility."

**Strength**: high (especially at translational journals).

**Response**: We frame the panel as a **risk-stratification readout** that complements driver mutation status, not as a clinical biomarker or diagnostic. The proposed downstream uses are:

- Pre-RAI triage for redifferentiation therapy candidates.
- Trial enrichment for RAI-refractory studies.
- Functional readout in MAPK-inhibitor redifferentiation trials.

No claim of immediate clinical deployment. A prospective study would be the next step.

**Supporting evidence**: `docs/manuscript_strategy.md` § "Decision diagram"; Figure 5c.

**Action**: language audit before submission — no "biomarker", "diagnostic", or "predicts RAI response" without Tier-1 prospective validation.

---

### R7. "The 4-pillar (RNA + sc + protein + methylation) is already in Paper 1 — what's new here?"

**Strength**: medium (only relevant if Paper 1 and this paper overlap reviewers).

**Response**: Paper 1's 4-pillar establishes **molecular biology of the differentiation-silencing axis**. This paper newly contributes **Tier 1/2 label-anchored RAI-response validation** (GSE151179 + Boucai 2023 + Mu 2024) plus **redifferentiation cohort direction-of-effect** (Tier 5). The two papers are complementary, not redundant.

**Supporting evidence**: explicit cross-reference table in Discussion § "Relationship to companion paper".

**Action**: cite Paper 1 R17 results in Methods, keep all RAI-label-specific analyses in this manuscript.

---

### R8. "Sample sizes are small."

**Strength**: high.

**Response**: Acknowledge directly. Three responses:

- We pool across cohorts where designs permit (random-effects meta-analysis of standardized mean differences).
- We emphasize effect sizes and confidence intervals over p-values for small-n datasets.
- We use Mu 2024's larger 220-patient cohort (mutation-frequency-only) as the **driver-overlap anchor**.

**Supporting evidence**: forest plot of effect sizes with 95% CI in Figure 3 supplement; sample-size table in every figure caption.

---

### R9. "Cross-cohort batch effects."

**Strength**: medium.

**Response**: Panel scores are **within-cohort z-scored** before any cross-cohort comparison. We never combine raw expression. For each dataset we report normalization status and platform. ComBat or rank-based harmonization is used only for the meta-analysis forest, not for primary tests.

**Supporting evidence**: Methods § "Cross-cohort harmonization"; per-dataset normalization audit in `data/interim/<dataset>_normalization_audit.txt`.

---

### R10. "Why not deep learning?"

**Strength**: low (but mentioned).

**Response**: This is an **interpretable, reproducible 8-gene readout** built for clinical translation context. Deep models on small RAI-labeled cohorts (n ≤ 50) would overfit. A parsimonious panel with documented per-gene biology is the right unit for this disease and sample-size regime.

**Action**: response only — no figure needed.
