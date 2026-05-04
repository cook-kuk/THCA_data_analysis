# 8-gene curated vs de novo / driver-map strategy — final decision memo

**Date:** 2026-05-04
**Owner:** Seungho Cook
**Mode:** Strategy / report-only. **No new analysis. No download. No RunPod. No Paper 3/4 touch. No voice-protected prose.**
**Scope:** Paper 1 (DM1 molecular dark matter, per `paper_numbering_2026_05_04`) — central-axis decision for the 5/4–6/13 marathon and bioRxiv 6/13.
**Authoritative anchors used (read-only):**

- `project/reports/2026_05_04_paper1_dm1_full_molecular_only_lock.md` — Paper 1 manuscript-safe envelope (16 claims, forbidden language)
- `project/reports/2026_05_04_paper1_dm1_full_conflict_audit.md` — claim-by-claim source verification
- `project/reports/2026_05_04_image_dm1_final_nogo_decision.md` — H&E-DM1 closure (NO-GO)
- `project/reports/PAPER1_MARATHON_AUDIT_BUNDLE_2026_05_04.md` — numerical / xref / bib audit
- `project/reports/2026_05_03_results_R1_R5_prose.md` — fact-only Results scaffolding
- `project/reports/2026_05_03_manuscript_v8_OUTLINE.md` — pre-2026-05-04 5-pillar outline (now sliced into Paper 1 / Paper 2 / Paper 4)
- `project/reports/paper1_spatial_supplement_freeze_2026_05_03.md` — ST supplement freeze (verdict C)
- `project/reports/2026_04_30_taskB_8gene_audit.md` — 8-gene forensic audit (TIERA67 / TDS_core subset)
- Memory: `v17_8gene_audit_2026_04_29` · `v17_8gene_pangenome_robustness` · `v17_dark_matter_pivot_2026_04_29` · `v17_marathon_mode_post_pillar1` · `paper_numbering_2026_05_04` · `v18_paper2_HT_isolated`

---

## 1. Executive verdict

**Demote 8-gene out of the title and out of the discovery framing. Keep it as a compact RAI-lineage / dedifferentiation readout that operationalises an axis that is genuinely robust at the broader transcriptional level (TIERA67 ARI ≈ 0.90 ≈ pan-genome top-5000 ARI ≈ 0.92).** The 8-gene panel is not novel as a panel (8/8 overlap with Yoo 2016 TDS-16) and ARI=0.49 alone is too weak to support headline framing. The novel claim is the **mechanistic axis** that the score reads out: lineage-TF collapse (FOXE1↓ NKX2-1↓ PAX8↓) → DNMT1↑ / STAT3↑ / FOSL1↑ → TROP2 re-expression, with prognostic association independent of BRAF/RAS/age/stage and supportive spatial validation across 28 ST slides. Therefore: **lineage/RAI-silencing axis = the lead; 8-gene = the score; driver landscape = the orthogonality argument.** The "BRAF/RAS-negative dark matter" sub-population claim is **honest-negative supp only** (N1 n=156, HR=1.20, p=0.80, NS) — not a headline.

Recommended one-line position: **"a transcriptional differentiation-axis subtype, anchored in canonical RAI biology and orthogonal to canonical driver mutations, identifies a TROP2-targetable thyroid-cancer population."**

---

## 2. Why BRAF/RAS alone is not enough

The advisor feedback ("BRAF/RAS dichotomy is exhausted; novelty lives beyond it") is consistent with the published literature:

1. **TCGA Cancer Genome Atlas Research Network 2014 (Cell)** already operationalised BRAF-like vs RAS-like with the BRAF-RAS Score (BRS) and a thyroid differentiation score (TDS). Any paper that retreads this dichotomy without a mechanism on top will read as derivative.
2. **Yoo SK 2016 (PLOS Genet, SNU-GMI Korean PTC)** built the 16-gene transcriptional differentiation panel that this work's 8-gene is a subset of. The 16-gene panel is the canonical RAI/lineage axis. Our 8-gene therefore inherits the framework — it does not establish it.
3. **Landa 2016 (JCI, PDTC/ATC)** established that dedifferentiation along the lineage axis is the dominant axis of clinical aggressiveness in advanced disease, with BRAF/RAS-+TERT co-mutations marking the high-risk end. This sets the precedent: **mutation status × differentiation axis** is the productive frame, not mutation status alone.
4. **Krishnamoorthy 2025 (Nat Commun, ATC)** reinforces dedifferentiation/TF-axis collapse as the mechanism in advanced disease (note: per memory `v17_landa2016_cite_save`, our Discussion §3.1 earlier mis-attributed several specific facts to Krishnamoorthy and was corrected to Landa 2016).
5. **Pan/Ge 2025 (Nat Commun, ATC proteogenomic)** identified RRP9, C5AR1, etc. as ATC-specific subtypes — again moving beyond BRAF/RAS into mechanism.

**Implication for our framing.** A "Paper 1 = BRAF-like vs RAS-like" framing has already been done. A "Paper 1 = a new 8-gene panel" framing reads as a smaller TDS-16. The defensible novelty has to come from one of:

- (i) a **mechanistic axis** that is orthogonal to BRAF/RAS and downstream-actionable (TROP2 vulnerability),
- (ii) a **molecular subtype within BRAF/RAS-negative tumors** (genuine "dark matter" — but our N1 n=156 outcome is NS, so this can only be supp),
- (iii) **fusion / TERT / EIF1AX / DICER1 driver expansion** beyond the BRAF/RAS axis (not in our current data envelope),
- (iv) **lineage-TF collapse → DNMT → STAT3 → TROP2** four-step framework as a transcriptional re-expression hypothesis (we have transcriptional evidence for this; methylation is absent).

We have (i) and (iv) cleanly. We do not have (iii). We have (ii) only as an honest negative.

---

## 3. Evidence supporting the 8-gene panel (and its limits)

### 3.1 External literature support for the score (not for the panel)

- **Yoo 2016 TDS-16 / Riesco-Eizaguirre & Santisteban 2014 (RAI biology)** — 8/8 of our 8-gene panel are TDS-16 members. The score reads out canonical iodide-uptake / thyroid-lineage biology.
- **TCGA 2014** — the differentiation axis was already an axis of the dataset; we are operationalising, not discovering.
- **Landa 2016 / Krishnamoorthy 2025** — dedifferentiation along lineage TFs is the recognised mechanism in advanced disease.

### 3.2 Internal data support (from molecular_only_lock §2 envelope, source-verified)

| # | claim | role of 8-gene |
|---|---|---|
| 1 | TCGA bulk DM1_like vs THYROID_NONOVERLAP r=−0.885, n=561, p=5×10⁻¹⁸⁸ | **score-derived** ; non-overlap lineage cross-validation (8-gene-free) |
| 2 | TCGA DM1 cluster vs DM1_like score (DM1=0.00 vs DM2=−0.66) | **score = anchor** for DM1 cluster |
| 3 | PFI HR=2.04 [1.15–3.61] p=0.015; DFI multivariate HR=1.41 p=0.025 (BRAF/RAS/age/stage adj., n=461) | **score-dichotomised** |
| 4 | TF activity collapse — FOXE1 −1.21 / NKX2-1 −0.65 / STAT3 +1.55 / FOSL1 +0.98 / DNMT1 +0.74 (n=261) | mechanism downstream of the score |
| 5 | TROP2/TACSTD2 Δmean +4.33, FDR=1.2×10⁻³² (n=261) | DE in DM1-high vs DM1-low (score) |
| 6 | tumor vs normal — DM1 p=1.8×10⁻²⁰; TROP2 1.4× tumor>normal (n=517) | ddx separation |
| 7 | Hallmark GSEA — IL6_JAK_STAT3 +8.0; OXPHOS −5.5; EMT, IFN-γ, TNF-α/NF-κB up (n=261) | mechanism |
| 8 | TF network — FOXE1-PAX8 r=+0.55 etc. (n=261) | within-axis coordination |
| 9 | Bootstrap 95% CI 12-slide [−0.997, −0.916] | spatial sample-mean robustness |
| 10 | 8-gene drop-one-out r ∈ [−0.987, −0.981] | **panel robustness — Supp referee Q** |
| 11 | Moran's I 0.36 mean, 26/28 perm p<0.05 | spatial structure |
| 12 | Bivariate Moran (DM1 × NONOVERLAP) all 12 negative | spatial anti-correlation, 8-gene-free |
| 13 | Margin gradient — 23/28 negative ρ | spatial validation |
| 14 | Pan-cancer THCA r=−0.892 outlier (33 cancers) | thyroid-specific framing |
| 15 | Dark-matter subset (BRAF-/RAS-, n=156) HR=1.20 p=0.80 NS | **honest negative — Supp** |
| 16 | Harmony 28-slide UMAP | ST QC backbone |

Plus from `paper1_spatial_supplement_freeze_2026_05_03.md`:

- **Non-overlap lineage cross-validation** (SLC26A4, IYD, DUOX1, DUOX2, TFF3, HHEX, GLIS3, DIO2; **zero gene overlap** with 8-gene): sample-mean Spearman ρ = −0.85 to −1.00 (p<0.01), Pearson r = −0.99 to −0.996, across 12 external slides — **the score reads out a thyroid-lineage signal that survives outside the 8-gene panel**. This is the strongest defence of the score.
- **Autoimmune-without-cancer negative control** (GSE248205, n=2 ctrl / 3 HT / 3 GD): DM1_like not raised above ctrl by HT/GD inflammation alone (MW p>0.10) — **score is not generic-inflammation-driven**.

### 3.3 Pan-genome cluster robustness (from `v17_8gene_pangenome_robustness` memory + R-PRO §R4)

- **8-gene panel alone:** ARI = 0.49 (modest)
- **TIERA67 67-gene pool:** ARI = 0.903
- **Pan-genome top-5000 MAD (unrestricted):** ARI = 0.918
- **Driver_anchor 12 genes alone:** ARI = −0.007 (random)
- **TIERA67 enrichment in pan-genome top-100:** hypergeometric p = 3 × 10⁻⁴

This is the **single most important piece of evidence in this memo**: the **DM1/DM2 cluster is a property of the transcriptome at large**, not of the 8-gene panel. The 8-gene panel is a *clinically interpretable readout of that axis*, not the axis itself. This permits the rhetorical move: "the cluster is robust at the broad transcriptome level; the 8-gene panel is the compact readout we recommend for translation."

### 3.4 Forensic 8-gene audit (`2026_04_30_taskB_8gene_audit.md`)

- 8-gene = **TDS_core (16) sub-category subset**, hardcoded in `v17_ULTIMATE_common.py`
- TIERA67 candidate pool **does include** BRAF/NRAS/HRAS/KRAS/RET/TERT in `[Driver_anchor]` — drivers were **not explicitly excluded**, they were *out-ranked* in univariate Cohen's d (Driver_anchor ranks #52, #56, #63, #66, #67 vs 8-gene ranks #4, #5, #18, #19, #25, #38, #40, #50)
- ΔAUC 8 vs 16 = 0.013 (CI overlap, NS) — **8-gene retains ~99% of TDS-16 discriminative power at lower panel cost**
- The earlier framing "drivers excluded by design" is **inaccurate** and should be replaced with "drivers were in the candidate pool but ranked at the bottom of univariate Cohen's d" — already corrected in R-PRO §R3 and the audit Q&A.

### 3.5 Limits of the panel (must appear in Limitations / Discussion)

- **8-gene panel alone ARI=0.49** — modest. Defended by clinical-interpretability trade-off, not by maximum statistical power.
- **N1 BRAF-/RAS- dark-matter HR=1.20, p=0.80 (NS, n=156)** — the score *is* prognostic in the full cohort but **not** within the BRAF-/RAS- subset. So "dark matter has worse outcome" is not a Paper 1 claim. The N1 result is honest-negative supp.
- **Q3 TROP2 spot-level co-localisation ρ=−0.016** — **TROP2 is a *tumor-level* vulnerability, not a *DM1-spot-level* vulnerability**. Forbidden language: "DM1-high spots are TROP2-high."
- **H&E → DM1 NO-GO** — closure battery 2026-05-04 verdict D. **Cannot revisit during marathon.** Forbidden language: "H&E-inferable molecular subtype," "morphology-derived DM1."
- **Methylation absent.** TF→DNMT→TROP2 is "consistent with a 4-step framework," not "demonstrates."
- **Functional validation absent** (no organoid, no IHC, no IC50). TROP2 ADC actionability is a hypothesis, not a demonstrated dependency.

---

## 4. Candidate-axis tournament table

Scoring rubric (0–5 per axis): 1 = unsupported / unsafe; 5 = strongly supported / safe. Total /50. **Higher = better fit for Paper 1 main axis.** Read column headers as: N=Novelty beyond BRAF/RAS · D=Data availability today · CL=Clinical association · MV=Multivariate independence · M=Mechanistic coherence · TX=Therapeutic actionability · EX=External / ST / scRNA support · S=Reviewer simplicity · BD=Boundary safety w/ Paper 2/3/4 · OC=Overclaim risk inverted (5 = low risk).

| ID | Strategy | N | D | CL | MV | M | TX | EX | S | BD | OC | **Total** | Rank |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **A** | 8-gene-curated *as the headline panel* | 1 | 5 | 4 | 4 | 3 | 3 | 4 | 5 | 4 | 2 | **35** | 5 |
| **B** | Broader TDS-16 / RAI / lineage-axis lead, 8-gene as compact readout | 3 | 5 | 4 | 4 | 4 | 3 | 5 | 4 | 5 | 4 | **41** | **2 (TIE)** |
| **C** | Full driver-map first (BRAF/RAS/TERT + RET/NTRK/ALK fusions) | 3 | 1 | 3 | 3 | 2 | 3 | 1 | 5 | 4 | 3 | **28** | 8 |
| **D** | BRAF/RAS-neg dark-matter fusion / epigenetic | 4 | 1 | 1 | 1 | 2 | 2 | 1 | 3 | 4 | 2 | **21** | 9 |
| **E** | Epigenetic / DNMT / TF-collapse (mechanism layer) | 4 | 4 | 4 | 4 | 4 | 3 | 3 | 4 | 4 | 3 | **37** | 4 |
| **F** | Therapeutic vulnerability — TROP2 / NAMPT / KCNN4 / LYN | 4 | 4 | 3 | 3 | 3 | 5 | 3 | 4 | 4 | 3 | **36** | 6 (TIE) |
| **G** | De novo discovery (NMF / consensus / WGCNA / sparse PCA) | 4 | 2 | 3 | 3 | 3 | 3 | 2 | 3 | 4 | 4 | **31** | 7 |
| **H** | Metabolic reprogramming-only (OXPHOS / glycolysis / hypoxia) | 2 | 3 | 2 | 2 | 3 | 2 | 2 | 3 | 4 | 3 | **26** | — |
| **I** | Immune / ICI-readiness | 3 | 3 | 3 | 3 | 3 | 4 | 3 | 4 | **1** | 3 | **30** | (Paper 3) |
| **J** | Korean germline / NRG1 / somatic | 3 | 2 | 2 | 2 | 2 | 2 | 2 | 3 | **1** | 3 | **22** | (Paper 4) |
| **★** | **Lineage/RAI-silencing axis + driver orthogonality + TF→DNMT→STAT3→TROP2 mechanism + multivariate Cox + ST validation; 8-gene = compact readout** (B + E + F integrated) | 4 | 5 | 4 | 4 | 4 | 4 | 5 | 4 | 5 | 4 | **43** | **1** |

Notes on scoring:

- **Strategy A (8-gene-headline)** loses points on novelty (8/8 = TDS-16 subset; not novel as a panel) and overclaim risk (positioning a non-novel panel as the discovery puts the paper at reviewer risk).
- **Strategy C (full driver map first)** loses on data: we have BRAF/RAS/TERT mutation status and transcript expression but **no fusion calls (RET, NTRK, ALK)**. The "driver map" framing requires data we do not have today.
- **Strategy D (BRAF/RAS-neg dark-matter as headline)** is the **memory-default** framing (`v17_dark_matter_pivot_2026_04_29`) but is contradicted by N1 NS (n=156, HR=1.20, p=0.80). Headline-grade only if we had fusion / methylation / external survival validation in this subset, which we don't.
- **Strategy I (Immune/ICI)** belongs to Paper 3 (`v19_paper3_ici_track_a` FROZEN). Boundary-violation if pulled into Paper 1.
- **Strategy J (Korean germline)** belongs to Paper 4 backlog. Boundary-violation if pulled into Paper 1.
- **Strategy ★ (integrated)** is the best total because it stacks B (axis robustness) + E (mechanism coherence) + F (translation hook) into a single argument while keeping Paper 2/3/4 boundaries clean.

---

## 5. Recommended final architecture

**Driver-orthogonality first → broader axis robustness → compact readout → mechanism → outcome → spatial validation → therapeutic vulnerability.**

| Layer | Claim | Source / envelope row |
|---|---|---|
| 1. Driver landscape orthogonality | BRAF/RAS/TERT mutation × transcript: |d|<0.4; AUC<0.7; Driver_anchor cluster ARI=−0.007 | R3, lock-claim — derived from `t_driver_features.tsv`, `m_master_regulator_volcano.png` |
| 2. DM1/DM2 cluster axis | Pan-genome top-5000 MAD ARI=0.918 ≈ TIERA67 ARI=0.903 ≫ Driver_anchor ARI=−0.007 | R4, memory `v17_8gene_pangenome_robustness` |
| 3. Compact 8-gene readout | 8/8 ⊂ TDS-16; ΔAUC vs 16 = 0.013 (NS); drop-one-out r ∈ [−0.987, −0.981] | Audit §A4–A5; lock-claim 10 |
| 4. Mechanism (TF → DNMT → STAT3/AP-1 → TROP2) | TF activity collapse FOXE1 −1.21 / NKX2-1 −0.65 / STAT3 +1.55 / FOSL1 +0.98 / DNMT1 +0.74; Hallmark IL6_JAK_STAT3 +8.0 / OXPHOS −5.5 | lock-claims 4, 7 |
| 5. Outcome | PFI HR=2.04 [1.15–3.61] p=0.015; DFI multivar HR=1.41 p=0.025 (BRAF/RAS/age/stage adj.) | lock-claim 3 |
| 6. Spatial validation | Bulk r=−0.885 (n=561); 8-gene-free non-overlap lineage ρ ≤ −0.85 (12 slides); Moran 26/28 sig; bivariate Moran all 12 negative; autoimmune neg ctrl HT/GD ≯ ctrl | lock-claims 1, 11, 12, 13 + spatial freeze §a, §b |
| 7. Tumor-level therapeutic vulnerability | TROP2 Δmean +4.33 FDR=1.2×10⁻³²; tumor vs normal 1.4×; framed strictly as **tumor-population vulnerability** (Q3 spot-level NEG → no co-localisation claim) | lock-claims 5, 6 |
| 8. Honest negatives (Supp + Limitations) | N1 BRAF-/RAS- subset NS (n=156, HR=1.20); H&E→DM1 NO-GO closure battery; methylation absent; functional validation pending | lock-claim 15; closure battery |

**De novo discovery does not enter the body of Paper 1.** It can appear as a one-paragraph methods-supp showing pan-genome top-5000 MAD ARI = 0.918 (already computed; that *is* the de novo arm) and the hypergeometric enrichment of TIERA67 in pan-genome top-100 (p = 3×10⁻⁴). No new analysis required.

---

## 6. Curated 8-gene strategy — strengths, weaknesses, sentences

### Strengths to use
- **Compact** — 8 transcripts, RT-qPCR-feasible, panel cost ~half of TDS-16.
- **Robust** — drop-one-out r ∈ [−0.987, −0.981]; ΔAUC vs TDS-16 = 0.013 NS.
- **Anchored** — every gene is canonical RAI/iodide-uptake or thyroid-lineage TF (Yoo 2016, TCGA 2014).
- **Mechanistically coherent** — score correlates with TF-collapse, DNMT, STAT3, TROP2 axes, and supports an interpretable 4-step framework.
- **Spatially validated by zero-overlap genes** — the axis the score reads out survives in non-overlapping thyroid-lineage transcripts (SLC26A4, IYD, DUOX1, DUOX2, TFF3, HHEX, GLIS3, DIO2).

### Weaknesses to acknowledge
- **Not novel as a panel** — 8/8 ⊂ Yoo 2016 TDS-16. Cannot be presented as a "new biomarker discovery."
- **ARI=0.49 alone is modest** — must be defended via TIERA67 / pan-genome ARI ≈ 0.92.
- **Drivers were not "explicitly excluded"** — they were in the candidate pool but ranked at the bottom of univariate Cohen's d. Earlier "driver-excluded by design" framing is wrong (Task B audit) and is being corrected in Methods.
- **N1 BRAF-/RAS- subset HR=1.20 NS** — the score does not improve risk-stratification within the BRAF-/RAS- subset at n=156.

### Required sentences (template, fact-only — voice-protected wording belongs to author)
- "An eight-gene panel (SLC5A5/NIS, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, DIO1) provides a compact readout of canonical thyroid radioiodide-uptake biology (Yoo et al., 2016)."
- "The DM1/DM2 cluster definition is robust to candidate-pool restriction: pan-genome top-5000 MAD ARI = 0.918 ≈ TIERA67 ARI = 0.903 (vs Driver_anchor alone ARI = −0.007)."
- "The eight-gene readout retains ≈99% of the discriminative power of the canonical 16-gene differentiation panel (ΔAUC 0.013, 95% CI overlap)."
- "DM1_like supports an independent prognostic association with PFI (HR 2.04 [1.15–3.61], p=0.015) and DFI (multivariate HR 1.41, p=0.025) after adjustment for BRAF, RAS, age, and stage."

### Forbidden sentences (per molecular_only_lock §5)
- "An 8-gene panel proves..." → use "supports."
- "All risks resolved." → false.
- "DM1-high spots are TROP2-high." → Q3 NEG.
- "H&E-inferable molecular subtype." → closure battery NO-GO.
- "Cancer Cell-ready" / venue-probability percentages.

---

## 7. De novo discovery strategy — design only, no execution

**Purpose.** Pre-empt reviewer Q: "Did you cherry-pick the panel from a curated pool?"

**Feature universe (3-arm).**
1. Biology-curated thyroid-lineage / RAI genes (TIERA67) — already done.
2. Transcriptome-wide top-N MAD (top-100, top-200, top-1000, top-5000) — already done at top-5000 (ARI 0.918) and top-1000 (ARI 0.902) and top-200 (ARI 0.864).
3. Driver / epigenetic / therapeutic genes (Driver_anchor) — already done (ARI −0.007).

**Methods.**
- KMeans k=2 on z-scored expression (current implementation).
- ARI vs original DM1/DM2 labels.
- Hypergeometric enrichment of TIERA67 in pan-genome top-100 by Cohen's d (p = 3×10⁻⁴, already computed).

**Validation.**
- Already covered: drop-one-out (lock-claim 10), bootstrap 95% CI (lock-claim 9), non-overlap lineage cross-val (spatial freeze §a).
- Not yet done and **not required for Paper 1**: NMF / WGCNA module recovery in TCGA bulk; consensus clustering with stability index. Can be a 2-paragraph methods-supp scaffold without rerunning anything if user wants.

**Possible outcomes (all already known; none requires a new run).**
- A. de novo top-N converges on TDS_core ∪ MAPK_output ∪ Aggressive_marker → **already observed; reinforces 8-gene as the compact subset**. Use this.
- B. de novo top-N is dominated by immune / stromal genes → rules out a tumor-intrinsic axis. **Not observed** (TIERA67 is enriched p=3×10⁻⁴; immune/stromal genes ranked outside top-100).
- C. de novo top-N is dominated by drivers → would invalidate orthogonality. **Not observed** (Driver_anchor ARI=−0.007).
- D. de novo top-N is dominated by metabolic genes → metabolism-only main axis. **Not observed** (OXPHOS down at NES=−5.5 is mechanism, not the cluster definer).

**Marathon discipline.** Write the de novo arm as a 1-page Supp Methods + 1 supp figure (already exists as F4 panel A). **Do not run NMF / WGCNA / consensus clustering during the marathon.** Defer to pre-revision response.

---

## 8. Fusion strategy — what we can and cannot say

### Driver landscape we have
- BRAF V600E mutation status (n=273 carriers / n=182 WT)
- HRAS / NRAS / KRAS hotspot mutation status
- TERT promoter mutation status (recovered v2: n=36 carriers; HR=4.33 lifelines per `v17_tert_recovery_v2`)
- Driver transcript expression for all of the above (Cohen's d, AUC)
- Driver_anchor 12-gene cluster (ARI=−0.007)

### Driver landscape we do NOT have
- RET / NTRK / ALK / kinase fusion calls (not in Paper 1 envelope; not in molecular_only_lock §2)
- DICER1 / EIF1AX point mutation calls beyond what TCGA published
- BRAF+TERT and RAS+TERT joint analysis at this scope (TERT v2 was a separate npj submission package; whether it enters Paper 1 v8 is open — see §12)

### Recommended fusion-strategy framing (safe)
- **Drivers act as the orthogonality argument** ("DM1/DM2 axis is independent of canonical mutation status") — Layer 1 of the architecture.
- **Drivers do not act as a driver-map figure** in Paper 1 because:
  - We have point-mutation status only, not fusion calls.
  - Driver_anchor cluster ARI=−0.007 means "drivers cannot define the axis," not "drivers define a different axis we are now showing."
  - Adding a half-driver-map invites reviewer Q "where are the fusions?" with no good answer.
- **TERT** is the one driver where we have survival data with effect size (HR=4.33). **Decision pending (§12)**: include in Paper 1 as multivariate co-variate (already in Cox), or hold for a TERT-focused future paper. Recommended: keep TERT as a co-variate in the multivariate Cox (we already do), do not give it its own figure.

### Sentences to use (safe)
- "BRAF V600E, RAS-hotspot, and TERT promoter status do not propagate to discriminative driver-transcript signal (|d|<0.4, AUC<0.7); driver-only clustering yields ARI=−0.007."
- "The DM1/DM2 axis is orthogonal to canonical driver mutation status."

### Sentences to avoid
- "We provide a comprehensive driver map" — we don't.
- "RET/NTRK/ALK fusion-positive tumors..." — we did not call fusions.
- "Driver_anchor + 8-gene combined classifier" — adds nothing (ARI 0.918 is already saturated by transcriptome).

---

## 9. Figure plan

### Main figures (envelope-supported, source files verified in conflict_audit §2)

| panel | content | claim |
|---|---|---|
| **F1A** | Driver mRNA × mutation status (BRAF/HRAS/NRAS/KRAS) | R3, claim outside envelope but supported by R-PRO |
| **F1B** | Single-feature driver AUC for DM cluster | R3 |
| **F1C** | TIERA67 univariate Cohen's d ranking, with Driver_anchor highlighted at #52–67 | R3 |
| **F1D** | Driver_anchor-only cluster ARI = −0.007 | R3 |
| **F2A** | Pan-genome top-N MAD ARI ladder (top-100 / -200 / -1000 / -5000) vs TIERA67 vs 8-gene vs Driver_anchor | R4 |
| **F2B** | TIERA67 enrichment in pan-genome top-100 (hypergeometric p=3×10⁻⁴) | R4 |
| **F3A** | TCGA bulk DM1_like vs THYROID_NONOVERLAP scatter (r=−0.885, n=561) — **non-overlap validation, 8-gene-free** | claim 1 |
| **F3B** | TCGA DM1 cluster vs DM1_like score boxplot | claim 2 |
| **F4A** | TF activity volcano (FOXE1 / NKX2-1 / PAX8 down; STAT3 / FOSL1 / DNMT1 up) | claim 4 |
| **F4B** | Hallmark GSEA bars (IL6_JAK_STAT3 / OXPHOS / EMT / IFN-γ / TNF-α) | claim 7 |
| **F4C** | TROP2 / TACSTD2 DE volcano + drug-target overlay | claim 5 |
| **F4D** | Tumor vs normal — DM1 score + TROP2 panel | claim 6 |
| **F5A** | KM curve PFI by DM1_like dichotomy | claim 3 |
| **F5B** | Multivariate Cox forest (BRAF / RAS / age / stage / DM1_like) | claim 3 |
| **F5C** | DM1_like dichotomy across N1 BRAF-/RAS- subset (HR=1.20 NS — **honest negative inset** or move to Supp) | claim 15 |
| **F6** | Spatial — 4-panel: Moran's I per slide, bivariate Moran (DM1 × NONOVERLAP), 8-gene drop-one-out per slide, autoimmune neg ctrl (GSE248205) | claims 9, 11, 12; spatial freeze |

Five main figures (F1–F5) is conservative; F6 spatial may be split into Supp depending on space.

### Supplementary candidates

- 28-slide Harmony UMAP (claim 16)
- 8-gene drop-one-out per slide (claim 10) — referee-Q answer
- Pan-cancer THCA outlier (claim 14)
- TF network FOXE1-PAX8 / PAX8-DIO1 (claim 8)
- Bootstrap 95% CI 12-slide (claim 9)
- Margin gradient (claim 13)
- N1 dark-matter subset KM + forest (claim 15) — **if not in main F5C**
- GSE250521 cancer-progression PT→ATC (spatial freeze §a, direction-consistent ρ=−0.35 p=0.18 underpowered)
- GSE230424 PTC + Hashimoto-overlap non-overlap lineage cross-val (spatial freeze §a)
- GSE248205 autoimmune-no-cancer negative control (spatial freeze §b)
- Closure battery negative-feasibility 1-page Methods supp (Phase A NO-GO, no manuscript-facing claim)

### Internal-only — must NOT enter submission package

- §9 / §11 of `PAPER1_DM1_FULL_2026_05_04.md` (RunPod G1/G2/G3 plan — superseded by closure)
- venue-probability framing
- "Cancer Cell-ready" / "Nat Cancer reach" wording
- TROP2 spot-level co-localisation (Q3 NEG)
- H&E-DM1 / WSI-projection in any active-claim form
- Fusion-map figure (we do not have fusion calls)

---

## 10. Title variants — ranked

| # | title | pros | cons | rank |
|---|---|---|---|---|
| **T1** | "**A transcriptional differentiation axis stratifies thyroid cancer orthogonally to canonical driver mutations and identifies a tumor-population TROP2 vulnerability**" | lineage/RAI + BRAF/RAS-orthogonality + translation hook; novelty without overclaim; no panel name in title | longer | **1** |
| T2 | "**Beyond BRAF and RAS: a thyroid-lineage silencing axis defines a TROP2-targetable papillary thyroid cancer subtype**" | "Beyond BRAF/RAS" hook; concise | "silencing" without methylation may be challenged | 2 |
| T3 | "**Lineage transcription-factor collapse defines a dedifferentiation subtype of papillary thyroid cancer with TROP2 vulnerability**" | mechanism-led | "subtype" needs molecular definition spelled out elsewhere | 3 |
| T4 | "**Driver-mutation-independent dedifferentiation in thyroid cancer: an RNA-defined axis with prognostic and therapeutic implications**" | safe and accurate; matches data | bland; weaker hook | 4 |
| T5 | "**An eight-gene transcriptional axis stratifies papillary thyroid cancer beyond canonical driver mutations**" | 8-gene-forward | 8/8 ⊂ Yoo 2016 ⇒ panel novelty thin | **avoid as headline** |
| T6 | "**A BRAF/RAS-negative dark-matter subtype of papillary thyroid cancer**" | "dark matter" hook | N1 NS (HR=1.20 p=0.80, n=156) ⇒ contradicts data | **avoid** |
| T7 | "**Spatially-resolved differentiation axis in papillary thyroid cancer**" | spatial-led | spatial is supportive, not lead | avoid as headline |

**Recommended:** T1 or T2.

---

## 11. Claims to use vs avoid

### Manuscript-safe (cite-allowed)
- "Supports an independent prognostic association" (claim 3, multivariate-adj.)
- "Consistent with a 4-step framework: TF collapse → DNMT silencing → STAT3 / AP-1 activation → TROP2 re-expression" (claim 4 + 5 + 7)
- "Tumor-level / tumor-population TROP2 vulnerability" (claim 5)
- "Supportive spatial validation in 28 ST slides" (claims 9–13)
- "Pan-genome cluster robustness ARI ≈ 0.92 ≈ TIERA67 ARI ≈ 0.90" (R4)
- "Eight-gene panel retains ≈99% of TDS-16 discriminative power (ΔAUC 0.013, NS)" (audit §A5)
- "Functional validation remains required (TROP2 IHC, sacituzumab IC50, organoid)" — must appear in Limitations.

### Risky (only with hedging)
- "Beyond BRAF/RAS" — fine in Discussion, fine in title; **not** "BRAF-/RAS-negative subset has worse outcome" (N1 NS).
- "Dark matter" — fine as a *narrative* term; not as a *prognostic* claim within BRAF-/RAS-.
- "TROP2-targetable" — fine as hypothesis with "remains to be validated"; not "actionable in DM1 patients."

### Forbidden (per molecular_only_lock §5)
- `proves` / `establishes` / `all risks resolved`
- `Cancer Cell-ready` / venue-probability percentages
- `H&E-inferable` / `H&E predicts DM1` / `pathology-AI triage of DM1` / `morphology-derived DM1`
- `TROP2 spatially colocalizes with DM1-high spots` / `DM1-high spots are TROP2-high`
- `tile-level DM1 inference` / `WSI-validated DM1`
- `RunPod G1/G2/G3` / `image-DM1 Phase B/C`

---

## 12. Final recommendation

### 12.1 Final title (recommended)
**"A transcriptional differentiation axis stratifies thyroid cancer orthogonally to canonical driver mutations and identifies a tumor-population TROP2 vulnerability"**

(Backup: T2 "Beyond BRAF and RAS: a thyroid-lineage silencing axis defines a TROP2-targetable papillary thyroid cancer subtype")

### 12.2 Final main claim
The DM1/DM2 transcriptional axis stratifies papillary thyroid cancer along a thyroid-lineage / RAI-uptake silencing dimension that is (i) **robust at the unrestricted transcriptome level** (pan-genome ARI ≈ 0.92), (ii) **orthogonal to canonical BRAF / RAS / TERT mutation status** (Driver_anchor ARI = −0.007), (iii) **mechanistically coherent** with lineage-TF collapse + DNMT1↑ + STAT3 / FOSL1 activation, (iv) **independently prognostic** in multivariate Cox (DFI HR=1.41, p=0.025), (v) **operationalisable as a compact eight-gene readout** (≈99% of TDS-16 power; drop-one-out r ∈ [−0.987, −0.981]), and (vi) **associated with a tumor-population TROP2 / TACSTD2 expression window** that motivates evaluation of TROP2-directed antibody-drug-conjugate strategies in this subtype.

### 12.3 Final first-paragraph logic (skeleton, voice-protected — author keyboard)
- Anatomic risk-stratification gap in DTC.
- Lineage TF / RAI / TDS framework already operational since Yoo 2016 / TCGA 2014.
- Recurrence prediction beyond BRAF/RAS gap remains.
- We define a transcriptionally-anchored axis, validate it spatially with non-overlapping lineage genes, link it to TF-collapse and TROP2, and demonstrate independence from canonical drivers.
- (Author keyboard for the actual prose.)

### 12.4 What 8-gene becomes in the manuscript
- **Title:** removed.
- **Abstract:** appears in one sentence as "compact readout (eight RAI-uptake genes)" — not as the discovery.
- **Introduction:** mentioned in §2 (existing molecular landscape) as the TDS-16 / Yoo 2016 lineage panel that we restrict to a clinically interpretable 8-gene subset.
- **Methods:** Gene-panel selection rationale paragraph (Task B audit §3.1 draft).
- **Results:** appears in R3 (TIERA67 ranking) and R4 (ARI), and as the score behind R5/R6/R7 mechanism / outcome / TROP2.
- **Figures:** drop-one-out as Supp; not a main figure.
- **Discussion:** "compact readout, not a discovery panel; rationale = clinical translation" (audit §3.3 draft).
- **Limitations:** ARI=0.49 alone modest; TDS-16 inheritance acknowledged; N1 BRAF-/RAS- subset NS.

### 12.5 Cell Reports Medicine-safe version (the default ladder rung)
- T1 title.
- 5 main figures (F1 driver orthogonality / F2 ARI robustness / F3 axis + non-overlap validation / F4 mechanism + TROP2 / F5 outcome).
- F6 spatial in Supp.
- TROP2 framed as tumor-population vulnerability, **no functional claim**.
- Limitations explicit on (a) ARI 0.49 alone, (b) N1 NS, (c) functional validation pending, (d) H&E-DM1 NO-GO disclosed, (e) methylation absent.
- bioRxiv 6/13.

### 12.6 Nature Cancer reach version
- Same body. **Adds**: independent external survival validation of DM1_like in a second molecular cohort (e.g., Lee 2024 GSE213647 dichotomised by DM1_like; *not currently in Paper 1 envelope* — would need a small post-marathon analysis or Paper-2-overlap cohort). + functional TROP2 IHC on a small TCGA-FFPE panel (wet-lab; not available now). Without these, **stay at CRM**.

### 12.7 Cancer Cell aim — what is required (not promised)
- (a) Functional perturbation (FOXE1 / NKX2-1 KO restoring or eliminating RAI uptake in thyroid organoid).
- (b) Methylation evidence for DNMT1-mediated TROP2 re-expression.
- (c) Sacituzumab-govitecan IC50 in DM1-high vs DM1-low organoid / cell-line panel.
- (d) Prospective cohort with DM1_like-dichotomised follow-up.
- **None of these are obtainable in 5/4–6/13.** Cancer Cell is not a marathon target. This must NOT appear in any cover letter or planning doc as a stated venue.

### 12.8 Next concrete action (this week)
Per `2026_05_04_marathon_state_after_image_dm1_nogo.md` §5 priority order:

1. **Verify bib-m3m4 completion** (`2026_05_04_bib_m3m4_straggler_cleanup.md`). If incomplete, finish.
2. **Otherwise → voice-protected drafting (author keyboard only):** Hook ¶1 / Aim ¶4 / Discussion §3.1 / §3.4 Limitations / Cover letter ¶1 / Reviewer Q9.
3. **Otherwise / parallel safe-to-defer (Claude OK):** Methods scaffolding paragraph for the 8-gene-curated framing (audit §3.1 draft already exists; transcribe into Methods M-PRO or M-SCAF as a single-paragraph addition, with explicit "drivers were in the candidate pool but ranked at the bottom of univariate Cohen's d" correction). 1-page Supp Methods entry for the closure-battery negative-feasibility result. Update the v8 Outline title to T1.
4. **Bib citation pass for de novo / Yoo 2016 / Landa 2016 / TCGA 2014 / Krishnamoorthy 2025 / Pan 2025** (verify they are in `2026_05_03_references.bib`). The cite-save memo `v17_landa2016_cite_save` is the source of truth for the Krishnamoorthy → Landa correction.

### 12.9 What this memo does NOT do
- Does not run any new analysis.
- Does not modify any voice-protected section.
- Does not generate Hook / Aim / Discussion §3.1 / Limitations / Cover ¶1 / Q9 prose.
- Does not touch Paper 2 / 3 / 4.
- Does not re-attempt H&E-DM1.
- Does not commit any file.
- Does not update memory.

### 12.10 Single-line final recommendation
**Keep 8-gene as a compact RAI-lineage silencing readout, not as the sole discovery premise. Lead with the lineage-axis + driver-orthogonality + TF/DNMT/STAT3/TROP2 mechanism. Demote 8-gene out of the title. Keep "BRAF/RAS-negative dark matter" out of headline framing — it is honest-negative supplement only.**

---

## Appendix — quick rebuttal table for the 7 most likely reviewer challenges

| # | Reviewer challenge | Defence |
|---|---|---|
| 1 | "8-gene is just TDS-16 minus 8 genes." | Correct. 8-gene retains ≈99% of TDS-16 discriminative power (ΔAUC 0.013 NS) at lower panel cost. We frame it as compact readout, not novel discovery. |
| 2 | "ARI=0.49 alone is too weak." | The cluster definition is robust at the broader transcriptome level: pan-genome top-5000 MAD ARI=0.918, TIERA67 ARI=0.903 (vs Driver_anchor ARI=−0.007). 8-gene is a clinical-interpretability trade-off. |
| 3 | "Why not include RET/NTRK/ALK fusions in the driver map?" | Out of scope for this study; mutation/transcript status are presented to establish orthogonality, not as a comprehensive driver map. Future cohorts with fusion calls. |
| 4 | "TROP2 colocalises with DM1 spots?" | We tested — Q3 spot-level ρ=−0.016. We claim only tumor-population TROP2 vulnerability. Spot-level colocalization is explicitly disclaimed. |
| 5 | "Where is methylation?" | Not in this work. We frame the TF→DNMT→TROP2 model as "consistent with" rather than "demonstrates." |
| 6 | "Where is functional validation?" | Limitations explicit: TROP2 IHC, sacituzumab IC50, organoid perturbation are stated as future work. The paper makes no functional claim. |
| 7 | "Did you try H&E?" | Yes — closure battery 2026-05-04: NO-GO at r=0.022, AUROC=0.511. We pre-tested and explicitly disclaim image-based DM1 inference in this study. |

---

*Memo authored 2026-05-04 by Claude (Opus 4.7) under marathon-mode discipline. Read-only audit + design-only strategy. No new analysis, no download, no RunPod, no voice-protected prose generated, no Paper 3/4 touch, no manuscript file modified, no commit. All claims sourced from molecular_only_lock §2 envelope or pre-existing reports/memories listed in the header.*
