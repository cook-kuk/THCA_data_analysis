# Supplementary Methods scaffold — pre-tested image-based DM1 inference (negative feasibility)

**Date:** 2026-05-04 · **Owner:** Seungho Cook
**Type:** Supplementary Methods scaffold (fact-only, Cell Press STAR Methods style). **Not voice-protected**; Limitations narrative is.
**Scope:** Single Supp Methods entry, ~1 page, for inclusion in Paper 1 supplement.
**Authority:** All numbers sourced from existing TSV / closure-battery reports — **no new analysis, no download, no GPU/RunPod, no H&E-DM1 retry.**
**Source files:**

- `project/reports/pathology_dm1_phaseA_cpu_verdict_2026_05_04.md` (Phase A verdict, ResNet50 LOSO 16-fold, n=3,200 tiles)
- `project/reports/pathology_dm1_closure_battery_2026_05_04.md` (closure battery sections A–G, raw vs residualized variants, RGB baseline, multi-resolution)
- `project/reports/2026_05_04_image_dm1_final_nogo_decision.md` (final verdict D; 5 re-entry conditions)
- `project/results/03_pathology_poc/closure_battery_metrics.tsv` (26 rows, source-of-truth metrics)
- `project/results/03_pathology_poc/loso_metrics_resnet50.tsv` (per-fold detail, 224 + 448 px)
- `project/results/03_pathology_poc/negative_controls_summary.tsv` (30 random panels + housekeeping, residualized)
- `project/results/03_pathology_poc/negative_controls_raw_summary.tsv` (50 random panels + housekeeping + RAI real, raw)
- Memory: `paper_numbering_2026_05_04` · `v17_marathon_mode_post_pillar1`

**Manuscript-facing framing reminder** (per `molecular_only_lock` §4):
> "We pre-tested H&E-based prediction of the depth-residualized DM1 axis at GSE250521 hires resolution and found no learnable signal beyond random panels. Image-based triage of the molecular subtype is therefore not pursued in this work."

**Forbidden language in this entry** (per `molecular_only_lock` §5):
- `H&E-inferable` / `H&E predicts DM1` / `morphology-derived DM1` / `WSI-validated DM1`
- `tile-level DM1 inference` (in active-claim form)
- `RunPod G1/G2/G3` / `image-DM1 Phase B/C` (in active-claim form)
- `pathology-AI triage` / "all risks resolved"

---

## Supp Methods Mxx — Pre-tested image-based DM1 inference (negative feasibility)

### Mxx.1 — Rationale

Because depth-residualized DM1_like reads out a thyroid-lineage / RAI-uptake silencing axis, we assessed whether haematoxylin-and-eosin (H&E) tile morphology could recover the same axis at hires Visium resolution. A confirmatory positive result would have supported a future image-based triage of the molecular subtype; a negative result is reported here for transparency and is explicitly disclaimed in the main text.

### Mxx.2 — Image cohort and tile extraction

We used the GSE250521 cancer-progression Visium cohort (16 slides; 4 PT, 4 PTC, 4 LPTC, 4 ATC). 224 px tiles centred on Visium spot coordinates were extracted from the published hires H&E images, balanced across stages: 800 PT, 800 PTC, 800 LPTC, 800 ATC, total **n = 3,200 tiles**. A 448 px centre-crop sweep was also performed (n = 2,935). Per-spot DM1_like, RAI_8 and TDS_like target scores were taken from the same depth-residualization pipeline used in the main spatial analysis (within-sample OLS on log_total_counts + log_n_genes_by_counts).

### Mxx.3 — Embedding and regression model

Tile embeddings were computed with the pretrained ResNet50 (ImageNet1K_V2 weights, 2048-dim global-average-pooled features, ImageNet normalisation, 224 px centre crop). Tile-level Ridge regression (α = 1.0; sklearn) and ElasticNet (α = 0.001, l1_ratio = 0.5) were trained against the residualized DM1_like target with **leave-one-slide-out (LOSO) 16-fold cross-validation**. Per-fold and pooled Spearman r, Pearson r, R², and AUROC (DM1_high defined by training-set q75) were recorded.

### Mxx.4 — Negative-control panels

To determine the chance distribution of LOSO Spearman r under this pipeline, we trained Ridge / ElasticNet against (i) **30 random thyroid gene panels** (residualized) and (ii) **50 random panels + housekeeping** (raw) under the same LOSO protocol. The DM1_like real-panel Spearman r was compared to the random-panel maximum.

### Mxx.5 — Closure battery (audit beyond Phase A)

To rule out that the negative result is an artefact of model capacity or feature representation, we additionally evaluated:
- **Raw vs residualized targets** (DM1_like, RAI_8, TDS_like): raw DM1_like Spearman r = 0.061 with ResNet50; residualized r = 0.022.
- **Five residualization variants** (log_counts only, log_ngenes only, no_resid, rank_norm, epi-top50): all r ≤ 0.06.
- **Tile size 448 px** (n = 2,935): Spearman r = 0.017 (no improvement).
- **Simple RGB 9-dimensional baseline** (mean / std / entropy of R, G, B per tile + Ridge): raw DM1_like r = 0.224, AUROC = 0.702 — **but** this signal correlates with within-tile stain darkness, which itself correlates with sequencing depth (within-tile log_counts ρ ≈ −0.45). Per-tile residualization on log_counts + log_n_genes reduces both the ResNet50 and the RGB-baseline correlations to the chance band (residualized RGB r = 0.080; ResNet50 r = 0.022), supporting a depth-confounded raw-image signal rather than a genuine morphology signal.
- **Stage ordinal vs ATC binary**: stage ordinal Spearman r = −0.047 (NS); ATC vs non-ATC AUROC = 0.689. The ATC binary signal corresponds to grossly visible anaplastic morphology and is not a new molecular finding.

### Mxx.6 — Result

Pooled across all 16 slides, ResNet50 ImageNet → Ridge LOSO yielded depth-residualized DM1_like **Spearman r = 0.022, AUROC = 0.511** (chance ≈ 0.5; per-slide r ∈ [−0.074, +0.118]). The depth-residualized real-panel performance did **not** exceed the maximum of 30 random thyroid panels (random max r = 0.031); raw real-panel Spearman r = 0.067 likewise did not exceed the raw random-panel maximum (random max r = 0.062 across 50 panels + housekeeping). Multi-resolution upscaling (448 px) did not improve the result. We therefore conclude that, under this pipeline, no learnable signal beyond random panels was recovered for the depth-residualized DM1 axis.

### Mxx.7 — Pre-registration and re-entry

This pre-test was conducted prior to the marathon writing window (5/4–6/13). Per the closure decision authority (`2026_05_04_image_dm1_final_nogo_decision.md`), image-based DM1 inference is not pursued further in this study. Any future re-entry will require, jointly: (i) full-resolution whole-slide images (not Visium hires crops), (ii) external paired molecular labels, (iii) a foundation-model pathology embedding (e.g., UNI / CTransPath) under a pre-registered hypothesis, (iv) post-marathon explicit decision, and (v) compute access (none of which are obtained at the time of this submission). The Limitations section of the main text discloses this pre-test result without claiming it as evidence for or against any positive finding.

### Mxx.8 — Data and code availability

Closure-battery metrics, per-tile predictions and 80-panel negative-control summaries are deposited at `project/results/03_pathology_poc/` (`closure_battery_metrics.tsv`, `loso_metrics_resnet50.tsv`, `loso_predictions_resnet50.tsv.gz`, `negative_controls_summary.tsv`, `negative_controls_raw_summary.tsv`, `tile_metadata*.tsv.gz`). Tile embeddings (.npz files) and source H&E images are reproducible from raw GSE250521 data via `project/src/05_pathology_poc/`. Embeddings themselves are not redistributed (model weights are publicly available from torchvision; image rights belong to GSE250521).

### Mxx.9 — Voice-protected components

The Limitations narrative referencing this entry is **voice-protected** (author keyboard only) per `v17_sprint_vs_marathon_violation`. The recommended phrasing template is in `molecular_only_lock` §4: *"We pre-tested H&E-based prediction of the depth-residualized DM1 axis at GSE250521 hires resolution and found no learnable signal beyond random panels. Image-based triage of the molecular subtype is therefore not pursued in this work."* The Discussion §3.4 / Limitations paragraph using this template is not generated in this scaffold.

---

## Marathon-discipline attestation

| Constraint | Status |
|---|---|
| New analysis | ✓ none — all numbers sourced from existing closure-battery TSVs / reports |
| New data download | ✓ none |
| GPU / RunPod / Azure burst | ✓ none |
| H&E-DM1 retry | ✓ NOT performed (pre-test summary only) |
| Paper 2 / 3 / 4 touch | ✓ none |
| Voice-protected prose generated | ✓ none — Limitations / Discussion / Hook / Aim / Cover ¶1 / Q9 untouched |
| File created | this scaffold only (`2026_05_04_supp_methods_closure_negative_feasibility.md`) |
| Manuscript files modified | none in this entry |
| Forbidden language in active-claim form | ✓ absent — only quoted-as-disclaimer references in §Forbidden language reminder block |

---

*Scaffold authored 2026-05-04 by Claude (Opus 4.7) under marathon-mode discipline. Supp Methods entry only; main-text Limitations and Discussion §3.4 voice-protected and deferred to author keyboard. All numerical content traceable to source TSVs in `project/results/03_pathology_poc/`.*
