# Paper 1 — STAR Methods supplementary draft: H&E negative feasibility

**Date:** 2026-05-04 · **Owner:** Seungho Cook
**Status:** Claude-drafted scaffolding for author copy-into-manuscript review.
**Target file:** `project/manuscript_v8/07_star_methods.md` (new sub-section under "Quantification and Statistical Analysis" or new "Pre-registered failed hypotheses" sub-section).
**Authority:** `2026_05_04_paper2_post_image_dm1_nogo_status.md` §5 explicitly suggested this Methods supplementary slot. `2026_05_04_paper1_dm1_full_molecular_only_lock.md` §3 / §6 require negative-result framing be Methods-supplementary-only (NOT main text, NOT Limitations).

**Voice rule:** This is Methods (technical, factual, no rhetorical phrasing). NOT voice-protected. Claude scaffolding allowed. Author refines wording.

---

## Drop-in markdown (copy below into 07_star_methods.md)

```markdown
### H&E predictability of the depth-residualized DM1 axis (pre-registered negative feasibility)

To explicitly bound what histology can recover at hires Visium resolution, we
pre-tested whether spot-aligned H&E tiles can predict the depth-residualized
8-gene DM1_like score (DM1_like_score_resid; per-sample OLS residual on
log(total_counts) and log(n_genes_by_counts)) in GSE250521 (16 thyroid
Visium slides, 3,200 spot-centered tiles balanced 200/sample). Tiles
(224 px in hires-image space, ≈110 µm at the default tissue_hires_scalef)
were embedded with a frozen torchvision ResNet50 ImageNet1K_V2 backbone
(2,048-dim features, 224 px center crop, ImageNet normalization). Ridge
regression (α = 1.0) and ElasticNet (α = 0.001, l1_ratio = 0.5) with
fold-isolated StandardScaler were trained under leave-one-slide-out
cross-validation (16 folds). DM1_high binary discrimination used the
train-fold q75 threshold; pooled metrics also report the global q75
sensitivity. Negative controls: a housekeeping panel (ACTB, GAPDH, B2M,
HPRT1, PPIA, RPL13A, RPLP0, TBP) and 30 random 8-gene panels drawn from
genes detected in ≥30% of spots in every sample (excluding RAI_8 and
housekeeping) were processed through the identical residualization +
LOSO-Ridge pipeline. A 50-panel raw-target arm (no residualization)
quantified the noise floor in ResNet50 feature space.

ResNet50-on-DM1_like_score_resid yielded pooled Spearman r = 0.022 and
DM1_high AUROC = 0.511 (n = 3,200 held-out predictions); fold-level
Spearman ranged [−0.07, +0.12]. Random 8-gene panels under the same
residualized pipeline gave mean Spearman = +0.002 (max +0.031); the
housekeeping panel returned −0.008. The real DM1_resid signal was thus
indistinguishable from the random null. A 9-dimensional simple-RGB
baseline (per-tile mean, standard deviation, and Shannon entropy of
each channel) outperformed the 2,048-dim ResNet50 on the raw target
(r = 0.224 vs 0.061) but not on the residualized target (r = 0.080
vs 0.022), confirming that the apparent raw signal is captured almost
entirely by tile-level color statistics — a stain-darkness ↔
sequencing-depth artifact that within-sample residualization is
designed to remove. Larger context (448 px tiles, n = 2,935) gave
DM1_resid Spearman = 0.017 (slightly worse than 224 px); 672 px
extraction was prepared but the embed step was cancelled per the
pre-registered improvement gate (Δr ≥ +0.05 vs 224 px not met).
ResNet50 ImageNet did discriminate gross anaplastic thyroid carcinoma
(ATC vs non-ATC LogReg AUROC = 0.689) but failed on continuous
within-tumor molecular gradients (stage ordinal Ridge Spearman =
−0.047; epithelial and proliferation raw scores Spearman ≤ 0.03).

We therefore conclude that the depth-residualized DM1 axis does not
have a recoverable tile-level morphological correlate at GSE250521
hires resolution under the tested model class. The Hashimoto-resolution
limit is consistent with a published expected gain of +0.05–0.15
Spearman for histology foundation models (UNI, CONCH, Virchow2) over
ResNet50 ImageNet on tile-level molecular tasks, which falls short of
the +0.28 gap to our pre-specified GO threshold (Spearman ≥ 0.30). The
state-of-the-art for population-scale H&E + tumor biology multimodal AI
(GigaTIME; Valanarasu et al., 2026) achieves protein-level inference
via cross-modal H&E → virtual multiplex-IF translation trained on
~40 × 10⁶ paired cells across 14,256 patients and 24 cancer types,
demonstrating that H&E carries substantial information when paired
with molecular ground-truth at scale; our pre-registered test, by
contrast, asks a frozen ImageNet ResNet50 to recover a depth-
residualized 8-gene transcriptional axis from spot-aligned hires
Visium tiles — a different task class (continuous molecular
regression rather than image-to-image protein translation) at
radically smaller scale (3,200 tiles). The negative result therefore
bounds what is recoverable from H&E **alone** for this transcriptional
axis at our resolution, and motivates RNA-paired H&E training at
scale (rather than ImageNet-frozen tile embedding) as the only
realistic path to image-DM1 inference. H&E-based triage of the
molecular subtype is therefore not pursued in this work. All numerical
results, per-fold tables, negative-control panels, and pre-registered
gates are archived under `project/results/03_pathology_poc/`
(closure_battery_metrics.tsv, loso_metrics_resnet50.tsv,
negative_controls_summary.tsv, negative_controls_raw_summary.tsv,
closure_battery_summary.png, pred_vs_obs_resnet50.png).
```

---

## Notes for author review

- **Voice / register:** Methods-paper neutral. No "we propose", no "we demonstrate". Past-tense factual reporting of a negative pre-registered test.
- **Word count:** ~440 words. Trimmable to ~250 if STAR Methods page budget is tight (drop the simple-RGB baseline detail; keep the 224 px main result + the 448 px context test + the conclusion).
- **Numerical claims (all from committed files):**
  - 16 slides / 3,200 tiles 224 px, 2,935 tiles 448 px, 2,386 tiles 672 px (extracted only)
  - DM1_resid Spearman 0.022, AUROC 0.511 — `closure_battery_metrics.tsv` row 1 + `loso_metrics_resnet50.tsv`
  - Random 8-gene mean +0.002, max +0.031 — `negative_controls_summary.tsv`
  - Housekeeping −0.008 — same
  - 9-dim RGB raw 0.224 vs ResNet50 raw 0.061 — `closure_battery_metrics.tsv` F section + A section
  - 9-dim RGB resid 0.080 vs ResNet50 resid 0.022 — same
  - 448 px 0.017 vs 224 px 0.022 — `loso_metrics_resnet50.tsv`
  - ATC AUROC 0.689 — `closure_battery_metrics.tsv` D2 row
  - stage ordinal Spearman −0.047 — `closure_battery_metrics.tsv` D1 row
  - epithelial r 0.026, proliferation r 0.028 — `closure_battery_metrics.tsv` D3 rows
- **Citations to add (author):**
  - UNI: Chen et al. 2024 Nat Med (foundation model)
  - CONCH: Lu et al. 2024 Nat Med
  - Virchow2: Vorontsov et al. 2024 (paige-ai)
  - **GigaTIME / multimodal AI upper bound**: Valanarasu et al. 2026 Cell 189(2), DOI 10.1016/j.cell.2025.11.016 — already added to `project/manuscript_v8/03_intro_references.bib` as `@Valanarasu2026`
  - Visium platform: 10x Genomics SOP
  - GSE250521: Lu et al. spatial thyroid progression dataset

## Cross-reference

- `project/reports/2026_05_04_paper1_dm1_full_molecular_only_lock.md` §3, §5, §6
- `project/reports/2026_05_04_image_dm1_final_nogo_decision.md` (5 re-entry conditions)
- `project/reports/2026_05_04_paper2_post_image_dm1_nogo_status.md` §5 (suggested home for this Methods entry)
- `project/reports/pathology_dm1_closure_battery_2026_05_04.md` (long-form)
- `CLOSURE_BATTERY_2026_05_04.md` (top-level)

---

*Methods-supplementary draft only. Not main text. Not Limitations. Author copies into 07_star_methods.md after review.*
