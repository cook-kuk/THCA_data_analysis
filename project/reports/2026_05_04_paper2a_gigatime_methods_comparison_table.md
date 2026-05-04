# Paper 2A — GigaTIME vs our H&E negative-feasibility test (Methods supplement comparison table)

**Date:** 2026-05-04 · **Owner:** Seungho Cook
**Target file:** `project/manuscript_v8/07_star_methods.md` (Author copies into the H&E negative-feasibility supp sub-section, immediately following the multimodal-AI-context paragraph already drafted in `2026_05_04_paper1_supp_methods_he_negative_feasibility.md`).

**Authority:** `2026_05_04_paper2a_gigatime_methods_comparison_table.md` (this file) supplements `2026_05_04_paper1_supp_methods_he_negative_feasibility.md` (commit `cf14656` + `2ffe929`) with a structured side-by-side comparison. All numerical cells are sourced from `2026_05_04_gigatime_thca_relevance_scan.md` §1/§4 (confirmed metadata only) or our own committed result files.

**Voice rule:** Methods (factual, neutral). NOT voice-protected. Author copies and refines wording as needed.

---

## Drop-in markdown (copy below into `07_star_methods.md`)

```markdown
### Comparison to recent multimodal-AI baselines for H&E and tumor biology

To position our pre-registered negative result against the current
state of multimodal AI in tumor pathology, the table below summarizes
the principal scope and scale differences between our Visium-paired
H&E test and the most recent population-scale multimodal AI work
(GigaTIME, Valanarasu et al., 2026 Cell 189[2]:386–400.e19).

| Dimension | Our H&E negative-feasibility test | GigaTIME (Valanarasu et al., 2026) |
|---|---|---|
| Task class | Continuous regression of a depth-residualized 8-gene transcriptional axis (DM1_like_score_resid) from H&E tile embeddings | Cross-modal image-to-image translation: H&E → virtual multiplex-IF (mIF) protein-channel images |
| Output modality | Scalar transcriptional axis score per spot (RNA-derived ground truth) | Multi-channel protein expression image per cell (mIF-derived ground truth) |
| Model class | Frozen ResNet50 ImageNet1K_V2 (2,048-dim) + per-fold StandardScaler + Ridge / ElasticNet | Custom multimodal foundation model (architecture details in original paper) |
| Training paired data | None for the predictor (frozen ImageNet weights); ground-truth labels from 16-slide Visium ST | ~40 × 10⁶ cells with paired H&E ↔ mIF (~21 proteins) |
| Test cohort | GSE250521, 16 thyroid Visium slides, 3,200 spot-centered tiles (200/sample, balanced PT/PTC/LPTC/ATC) | 14,256 patients across 51 hospitals + ~1,000 clinics, Providence Health 7 US states |
| Cancer scope | Thyroid only (PT → ATC dedifferentiation trajectory) | 24 cancer types, 306 subtypes (thyroid inclusion not confirmed in public metadata) |
| Image modality | Spot-aligned hires Visium PNG (≈110 µm at 224 px tile, default tissue_hires_scalef) | Whole-slide histology + paired mIF (resolution and scanner specifications per original) |
| Validation strategy | Leave-one-slide-out CV (16 folds), fold-isolated scaler; pre-registered tile-size improvement gate (Δr ≥ +0.05 vs 224 px) | Independent 10,200 TCGA patient validation; 1,234 statistically significant protein-biomarker-staging-survival associations |
| Negative controls | Housekeeping panel (8 genes) + 30 random 8-gene panels (residualized) + 50 random panels (raw) | Not directly comparable (different task class) |
| Primary metric | LOSO pooled Spearman r and DM1_high (q75) AUROC on continuous score | Image-translation fidelity + association recovery |
| Result | Pooled Spearman r = 0.022, AUROC = 0.511; real signal indistinguishable from random gene panels (max 0.031) | 1,234 protein-biomarker-staging-survival associations recovered at population scale |
| Pre-registered gate met | No (Spearman < 0.20 BORDERLINE; AUROC < 0.65 BORDERLINE) | Not applicable to a different task |
| Scope of inference | Bounds H&E-alone recoverability of the depth-residualized 8-gene transcriptional axis at hires Visium 224 px tile resolution under ImageNet feature transfer | Establishes feasibility of population-scale H&E → mIF translation given paired training at scale |

The two studies are therefore complementary rather than competing:
GigaTIME demonstrates that H&E carries substantial predictive signal
for the tumor immune microenvironment when paired with molecular
ground-truth at the order of 10⁷ cells; our pre-registered test
establishes the depth-residualized 8-gene transcriptional axis as a
target that is not recoverable from H&E alone at our scale and
resolution under frozen ImageNet feature transfer. Bridging the two —
RNA-paired H&E training on a thyroid cohort at scale — is identified
as the principal future-work direction for any image-based triage of
the molecular subtype (see 5 re-entry conditions, Discussion).
```

---

## Notes for author review

- **Word count:** ~360 words (table + 2 framing sentences). Trimmable to ~220 by collapsing rows that author judges redundant.
- **All numerical cells are source-verified:**
  - GigaTIME side: from `2026_05_04_gigatime_thca_relevance_scan.md` §1 (confirmed metadata only); never from §5 (unsupported claims)
  - Our side: from `closure_battery_metrics.tsv`, `loso_metrics_resnet50.tsv`, `negative_controls_summary.tsv`, `negative_controls_raw_summary.tsv` (all committed)
- **Table placement:** immediately after the GigaTIME multimodal-AI-context paragraph that was integrated into the supp draft via commit `2ffe929`. Together they form a 2-paragraph + 1-table block.
- **Do NOT add:** any specific GigaTIME architecture details, AUC numbers, panel composition, or THCA-specific results until §5 of the scan memo is upgraded with verified data.
- **Citation:** `\citep{Valanarasu2026}` (BibTeX entry already in `03_intro_references.bib` via commit `2ffe929`).

## Cross-references

- `2026_05_04_gigatime_thca_relevance_scan.md` (S3 scan, this session)
- `2026_05_04_paper1_supp_methods_he_negative_feasibility.md` (S1 base supp draft, commits `cf14656` + `2ffe929`)
- `project/manuscript_v8/03_intro_references.bib` (`@Valanarasu2026`)
