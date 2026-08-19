# RAI response label taxonomy

A single project-wide taxonomy that every dataset, every plot, and every claim in this workspace must be mapped to. The manuscript must be **explicit** about which tier each dataset belongs to. This avoids the canonical reviewer trap of "you used TCGA as RAI response data".

---

## Tier 1 — Direct clinical RAI response

The strongest possible label. Requires post-RAI structural or functional re-imaging.

- **RECIST v1.1 response** after RAI therapy (CR / PR / SD / PD).
- **Exceptional responder** vs **non-responder** (Boucai 2023 CCR convention).
- **Structural tumor shrinkage** after RAI on cross-sectional imaging.
- **Progression despite RAI** (radiologic progression within 6–12 mo of RAI).

Datasets that *currently* qualify: Boucai 2023 CCR (controlled / supplement-level), some redifferentiation trials with post-RAI response endpoint.

---

## Tier 2 — RAI avidity / uptake

The next strongest tier. Measures whether the tumor took up RAI at all, regardless of whether tumor killing followed.

- **Metastatic RAI uptake yes/no** (I-131 or I-124 whole-body scan or PET).
- **RAI-avid vs non-avid** at the lesion level.
- **Uptake-pattern subgroups** (e.g., Mu 2024 JCEM's 4-class system: initially refractory / continually avid / gradually refractory / partly refractory).
- **Pre-/post-rhTSH-stimulated uptake** if available.

Datasets that qualify: GSE151179, GSE151181 (SuperSeries), GSE299988, Mu 2024 JCEM HRA004166. The 4-class Mu schema specifically captures the **gray zone** that this project's framing emphasizes.

---

## Tier 3 — Clinical disease status after RAI

Composite endpoints that mix RAI response with subsequent therapy and follow-up.

- **Remission** (excellent response, ATA 2015).
- **Persistent disease** (biochemical incomplete, structural incomplete).
- **Recurrent disease** after initial remission.
- **Indeterminate response** (ATA 2015).

Useful but noisier than Tier 2 because remission can be achieved by surgery + suppression even without RAI.

---

## Tier 4 — Molecular proxy

No clinical RAI label; use molecular surrogates that have been independently linked to RAI response.

- **Thyroid differentiation score (TDS)** — Landa 2016, TCGA-THCA Cell 2014.
- **Iodide-handling gene expression** — this project's 8-gene panel.
- **MAPK output / BRS** (BRAF-RAS score).
- **MAPK activation signature** (e.g., MAPK16-gene from TCGA-THCA Cell 2014).

Datasets: TCGA-THCA, GSE213647 (Lee 2024), GSE286332, Mun 2025 proteome, Lu 2023 sc.

These are **discovery / mechanism** anchors. Never claim they "predict RAI response" on their own.

---

## Tier 5 — Treatment-induced redifferentiation

Pre-/post-treatment paired expression with a re-RAI step.

- **Restored uptake after MEK / BRAF inhibitor** (selumetinib, dabrafenib, trametinib, vemurafenib).
- **Increased iodide-handling gene expression** after redifferentiation treatment.
- **MERAIODE / SEL-I-METRY trial supplements**.

Datasets: GSE112202 (digoxin), redifferentiation trial supplements. Tier 5 is mechanism + therapeutic relevance, not validation of the panel as a direct biomarker.

---

## Mapping rule for this workspace

Every dataset entry in `config/datasets.yaml` must declare `label_tier:` using one of the five strings above. Every figure caption and report must include the tier in parentheses, e.g., "validated against RAI avidity labels (Tier 2)". Every claim must respect the tier:

| Tier | Allowed claim language |
|------|------------------------|
| 1 | "associated with RAI response", "stratifies RAI responders from non-responders" |
| 2 | "associated with RAI avidity / uptake", "stratifies RAI-avid from refractory lesions" |
| 3 | "associated with persistence after RAI" |
| 4 | "captures a differentiation-linked RAI failure axis", "consistent with RAI refractoriness" |
| 5 | "consistent with redifferentiation biology", "panel score increases with restored iodide-handling" |

**Forbidden claim language** (per-project rule):

- "predicts RAI response" — only allowed when Tier 1 data + held-out validation + ROC AUC reported.
- "clinical biomarker" — only allowed after prospective or independent multi-center Tier 1 confirmation.
- "diagnostic" — never; this is risk stratification, not diagnosis.
