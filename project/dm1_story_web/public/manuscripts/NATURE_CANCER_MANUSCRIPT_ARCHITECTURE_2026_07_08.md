---
title: "Paper 1 Nature Cancer-style manuscript architecture"
date: 2026-07-08
author: Codex scaffold for Seungho Cook
status: scaffold only; voice-protected prose slots preserved
source_context:
  - project/dm1_story_web/
  - project/manuscript_v8/DM1_MASTER_BRIEF_FOR_GPT_AND_MEETING_2026_06_25.md
  - project/manuscript_v8/04_results.md
  - project/manuscript_v8/05_figure_captions.md
  - project/papers_hub_2026_0509/paper1_nature_cancer_board.html
voice_boundary:
  - Do not use this file to overwrite Hook, Aim, Discussion 3.1, Limitations 3.4, cover-letter paragraph 1, or Q9.
---

# Nature Cancer-Style Manuscript Architecture

## Editorial Position

This should not be framed as an 8-gene biomarker discovery paper. The high-impact version is:

> A compact thyroid-lineage state reveals a clinically actionable, driver-orthogonal route to iodine-handling collapse in thyroid cancer, linking dedifferentiation biology to avoidance of ineffective high-dose radioiodine.

The manuscript has to read as a lineage-state biology paper with a clinical decision consequence, not as a model paper. The score is the lens; the state is the discovery.

## Best Title Candidates

1. **A thyroid-lineage state marks iodine-handling failure in papillary thyroid cancer**
2. **A compact iodine-handling axis reveals radioiodine-poor thyroid cancer states**
3. **Driver-orthogonal thyroid differentiation collapse defines radioiodine-poor thyroid cancer**
4. **An 8-gene thyroid-lineage axis resolves radioiodine-poor dark matter in papillary thyroid cancer**
5. **A clinically translatable thyroid differentiation axis identifies iodine-handling-low thyroid cancer**

Recommended: title 1 for Nature Cancer style. It avoids sounding like a diagnostic-kit paper while keeping the biology central.

## One-Sentence Claim

A literature-anchored 8-gene thyroid differentiation and iodine-handling axis resolves a driver-orthogonal, radioiodine-poor molecular state that is reproducible across public cohorts, enriched for kinase-fusion biology and promoter methylation of thyroid machinery, and clinically positioned to reduce ineffective high-dose RAI exposure.

## Claim Boundary

Use this exact hierarchy throughout the paper:

| Claim tier | Safe claim | Avoid |
|---|---|---|
| Discovery | The axis identifies iodine-handling-low thyroid-lineage state | The panel predicts RAI response definitively |
| Biology | DM1 reflects thyroid differentiation collapse across RNA, methylation, proteome, and single-cell evidence | DM1 is proven causal mechanism |
| Mechanism | Fusion/MAPK and methylation evidence motivate a mechanism model | MAPK/methylation causes all DM1 |
| Clinical | The state may support RAI harm-avoidance and reflex testing | The state selects treatment today |
| Translation | SNUBH/FFPE/IHC validation is the next clinical lock | Public cohorts alone establish clinical utility |

## Abstract Skeleton

Do not paste this as final prose without author review. It is a structured content scaffold.

**Background:** differentiated thyroid cancer has excellent survival, but a clinically important subset receives radioiodine despite poor iodine-handling biology. Current molecular risk frameworks emphasize canonical drivers and do not directly measure thyroid-lineage functional state.

**Approach:** define an 8-gene iodine-handling / thyroid-differentiation axis using TG, TPO, TSHR, SLC5A5/NIS, DIO1, PAX8, NKX2-1/TTF-1, and FOXE1/TTF-2; evaluate the axis across TCGA-THCA, Korean cohorts, single-cell datasets, methylation, proteomic, and external expression validation.

**Findings:** the axis resolves iodine-handling-low DM1 and iodine-handling-high DM2 states; the partition is reproduced by pan-genome clustering (ARI 0.92), is not explained by driver-gene expression, shows strong cross-cohort direction consistency, and is supported by thyroid machinery promoter methylation and multi-modality differentiation loss. DM1 is enriched for kinase-fusion biology and has retrospective outcome association, but clinical response prediction remains prospective-validation territory.

**Interpretation:** the 8-gene axis is best positioned as a compact readout of thyroid-lineage collapse and a candidate harm-avoidance tool for identifying patients unlikely to benefit from repeated high-dose RAI.

## Results Spine

### 1. A compact thyroid-lineage axis resolves iodine-handling-low and iodine-handling-high states

Purpose: establish the biological lens and defend against "arbitrary eight genes."

Must show:

| Evidence | Main number |
|---|---|
| TCGA-THCA discovery | n = 504 |
| 8-gene panel | TG, TPO, TSHR, SLC5A5, DIO1, PAX8, NKX2-1, FOXE1 |
| Pan-genome recovery | ARI = 0.92 |
| TIERA67 recovery | ARI = 0.90 |
| Driver-only failure | ARI = -0.007 |
| BRAF mRNA neutrality | d = -0.044, p = 0.57 |

Professor-level paragraph direction:

> We began from the premise that radioiodine responsiveness is ultimately constrained by thyroid-cell function rather than by driver identity alone. We therefore treated iodine uptake, organification, hormone synthesis and lineage transcription factors as a single biological circuit, and asked whether a compact readout of this circuit could resolve clinically meaningful states in PTC.

### 2. The axis recovers clinically relevant dark matter not explained by canonical drivers

Purpose: make DM1/DM2 matter clinically before going deep into mechanism.

Must show:

| Evidence | Main number |
|---|---|
| Xing dark-matter rescue | 131 / 180, 73% |
| TCGA DM1 prevalence | 28.4% |
| Korean dark-matter enrichment | 37.8% |
| TCGA + MSK pooled OS | HR 2.53 [1.31, 4.89], I2 = 0% |

Boundary: retrospective outcome association, not prospective risk model.

### 3. External validation shows this is a lineage state, not a TCGA artifact

Purpose: reviewer defense against cherry-pick / platform / ancestry concerns.

Must show:

| Evidence | Main number |
|---|---|
| Master cross-cohort forest | mean d = 2.81, median d = 2.37 |
| Direction consistency | 80 / 80 cells direction-consistent |
| Lee 2024 Korean PTC | n = 632, d = 5.93 |
| GPL570 4-cohort validation | 4 / 4 rho <= -0.84 |
| Mun 2025 proteomics | n = 336, 7 / 7 protein direction-consistent |
| Pu 2021 single-cell | per-patient r = 0.798-0.886 |
| FFPE compatibility | KS p = 0.44 |

### 4. DM1 is linked to fusion/MAPK biology but not reducible to driver identity

Purpose: elevate from biomarker to mechanism-facing biology.

Must show:

| Evidence | Main number |
|---|---|
| DM1 fusion positivity | 63 / 82, 76.8% |
| DM1 vs DM2 fusion enrichment | OR 7.41, p = 1.9e-13 |
| RET-fusion capture | 27 / 33, 81.8% |
| SV missingness | MAR chi-square p = 0.56 |
| MAPK x Panel-8 pooled | rho = -0.327 [ -0.376, -0.278 ], n = 1,287 |

Boundary: MAPK/fusion evidence supports a mechanistic model but does not prove causality.

### 5. DM1 carries epigenetic silencing of thyroid differentiation machinery

Purpose: mechanism layer and Nature Cancer biological depth.

Must show:

| Evidence | Main number |
|---|---|
| TCGA HM450 | n = 503 |
| Mean 8-gene beta | DM1 0.385 vs DM2 0.253 |
| TPO promoter methylation | d = 2.30, p = 1.9e-18 |
| DIO1 | d = 1.24 |
| TSHR | d = 1.20 |
| PAX8 | d = 0.97 |
| Landa overlap | 5 / 8 genes |

Boundary: use "consistent with epigenetic silencing" and "motivates" language. Do not state causal methylation proof.

### 6. Clinical translation: harm avoidance before treatment selection

Purpose: make the paper clinically mature.

Main message:

> The most compelling clinical use is not "send patients to systemic therapy faster" as the first claim. It is "avoid repeated ineffective high-dose RAI in tumors whose lineage program already indicates iodine-handling collapse."

Must show:

| Evidence | Main number |
|---|---|
| GSE151179 post-RAI refractory similarity | d approx -1.0, MW p approx 1e-4 |
| ATA uncertainty zone | intermediate-risk RAI decision gap |
| FFPE compatibility | KS p = 0.44 |
| SNUBH roadmap | NGS / RNA / IHC validation pending |
| TSO500-only caveat | 3-gene AUC 0.76, not enough alone |
| Compact add-on | TSO500 + TPO + DIO1 AUC 0.940 |

## Main Figure Architecture

For Nature Cancer-style readability, reduce the main paper to 5-6 main figures and push audit layers to Extended Data.

| Figure | Title | Job |
|---|---|---|
| Fig. 1 | A compact iodine-handling axis resolves thyroid-lineage states | Discovery, why 8 genes, pan-genome defense |
| Fig. 2 | The iodine-handling-low state is clinically and ancestrally relevant | dark-matter rescue, Korean enrichment, survival association |
| Fig. 3 | The state generalizes across platforms and cellular resolution | external forest, per-gene consistency, single-cell/protein |
| Fig. 4 | DM1 is enriched for fusion/MAPK biology | fusion stack, RET capture, MAPK correlation |
| Fig. 5 | DM1 carries epigenetic silencing of thyroid machinery | HM450, Landa overlap, TDS-16 canonicality |
| Fig. 6 | Translation roadmap for RAI harm avoidance | FFPE, TSO500 caveat, SNUBH/IHC validation pathway |

Extended Data should carry stress tests: calibration mismatch, spatial negative/caveat panels, PRISM boundary, residualization, deconvolution, subset sensitivity, random modules.

## Nature Cancer Probability Audit

| Axis | Current strength | Editorial read |
|---|---|---|
| Novel state biology | Strong | A compact axis reading lineage collapse is compelling |
| Cross-cohort validation | Strong | 19 cohorts / multi-modality is a real asset |
| Clinical endpoint | Medium | Retrospective, event-limited; needs humble framing |
| Mechanism causality | Medium-low | Methylation/MAPK are coherent but correlative |
| Own cohort / wet-lab | Missing | Main reason Nature Cancer would be difficult |
| Translational specificity | Medium | RAI harm avoidance is better than treatment-selection overclaim |

Practical assessment: Nature Cancer is a stretch unless SNUBH/IHC or perturbation data arrives. The same package is much stronger for Nature Communications / JCI Insight / npj Precision Oncology if the claim boundary is kept clean.

## Professor-Level Manuscript Voice Rules

Use:

- "lineage state"
- "iodine-handling collapse"
- "compact readout"
- "driver-orthogonal"
- "clinically positioned"
- "retrospective association"
- "prospective validation"
- "harm avoidance"

Avoid:

- "predicts RAI response" unless a true response-labeled cohort is being discussed
- "causes"
- "proves"
- "treatment selection biomarker"
- "Nature Cancer ready" without the missing validation layer
- "DM1 is just fusion-positive" framing

## Author-Keyboard Slots

These must remain for Seungho Cook to write or approve:

| Slot | File / section | Instruction |
|---|---|---|
| Hook | 03_introduction.md / 04_intro_1_1_hook.md | Leave as author voice |
| Aim | Introduction Aim | Leave as author voice |
| Discussion 3.1 | mechanism interpretation | Use this scaffold only as factual input |
| Limitations 3.4 | limitation voice | Keep honest boundary, author-written |
| Cover letter paragraph 1 | 08_cover_letter.md | Author-written |
| Q9 | 09_reviewer_qa.md | Preserve "motivates rather than confirms" |

## Immediate Writing Plan

1. Convert `04_results.md` into the six-figure spine above.
2. Move older Fig. 7/8 mechanism density into Extended Data and reviewer reserve.
3. Rename DM1/DM2 in the paper-facing layer as iodine-handling-low / iodine-handling-high, with DM1/DM2 retained as internal labels.
4. Build a one-page "claim boundary table" for Methods or Supplement.
5. Prepare SNUBH/IHC validation language as future/ongoing validation, not as completed evidence.

