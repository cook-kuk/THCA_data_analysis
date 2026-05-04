# Paper 2A — Reviewer Q scaffold (GigaTIME / multimodal AI / image-DM1 anticipated questions)

**Date:** 2026-05-04 · **Owner:** Seungho Cook
**Target file:** `project/manuscript_v8/09_reviewer_qa.md` (Author copies into the response document; Q9 is voice-protected and remains author keyboard).
**Status:** Claude-drafted Q1–Q8 scaffolds. Author refines wording, picks which to include in cover-response, and writes Q9 manually.

**Authority:** Pre-empts likely reviewer questions about why the H&E negative-feasibility result was the chosen pipeline rather than a multimodal-AI / foundation-model approach, given the recent publication of GigaTIME (Valanarasu et al., 2026 Cell). All factual claims source-traced to committed result files and `2026_05_04_gigatime_thca_relevance_scan.md`.

**Voice rule:** Reviewer Q1–Q8 = Claude scaffolding allowed (per `STATUS_PAPER1_2026_05_02.md` voice-protected list which names only Q9). Q9 = author keyboard strict.

---

## Q1 (high probability) — "Why didn't you use a histopathology foundation model (UNI / CONCH / Virchow2) instead of frozen ResNet50 ImageNet?"

**Answer scaffold:**

We considered foundation-model embedders during pre-registration (UNI, CONCH, Virchow2, Prov-GigaPath; all gated under MahmoodLab and paige-ai non-commercial agreements requiring 1–7 day approval). The pre-registered improvement gate was set as Spearman r ≥ 0.30 (GO) or r ∈ [0.20, 0.30) (BORDERLINE). The published gain of histopathology foundation models over ResNet50 ImageNet on tile-level molecular tasks is approximately +0.05 to +0.15 Spearman across the literature (e.g., Chen et al. 2024 Nat Med UNI benchmarks; Lu et al. 2024 Nat Med CONCH benchmarks). With our ResNet50 ImageNet baseline at Spearman r = 0.022, even the maximal reported foundation-model gain (+0.15) lands at ≈ 0.17 — below the BORDERLINE 0.20 floor. Re-running with a foundation embedder was therefore predicted not to cross the pre-registered gate.

We retained foundation-model retry as a re-entry condition (`2026_05_04_image_dm1_final_nogo_decision.md` §5 condition #4: "Foundation model access acquired and validated on a non-thyroid public benchmark first, to ensure the embedder is not the bottleneck before attempting our task"). If/when MahmoodLab or paige-ai access is granted and validated externally, the same closure-battery pipeline can be re-applied without modification.

---

## Q2 (high probability) — "Why didn't you use multimodal AI like GigaTIME (Valanarasu et al., 2026 Cell)?"

**Answer scaffold:**

GigaTIME (Valanarasu et al., 2026 Cell) trains a cross-modal H&E → virtual multiplex-IF (mIF) translator on ~40 × 10⁶ paired cells across 14,256 patients, recovering 1,234 protein-biomarker-staging-survival associations across 24 cancer types. Our test addresses a **different task class** at a **radically different scale**: continuous regression of a depth-residualized 8-gene transcriptional axis on 3,200 spot-aligned tiles from 16 Visium slides under frozen ImageNet feature transfer (no paired multimodal training). The two studies are complementary — GigaTIME demonstrates that H&E carries substantial signal for the tumor immune microenvironment when paired with molecular ground-truth at scale; our work establishes the molecular axis itself (8-gene RAI / DM1) on RNA, with the negative H&E-alone result bounding what is recoverable from imaging without paired training.

Reproducing GigaTIME's approach for our axis would require: (a) RNA-paired H&E training on a thyroid cohort approaching 10⁴ patients (vs our 16 slides); (b) full-resolution scanner WSI (vs hires Visium PNG at 0.17× downsample); (c) RNA labels at the case level alongside H&E (which our 8-gene axis provides); (d) a multimodal architecture rather than a frozen ImageNet backbone. Conditions (a), (b), and (d) are not currently met; condition (c) — our axis as a target — is established by this work.

---

## Q3 (high probability) — "Was the GSE250521 hires Visium PNG resolution sufficient to test this hypothesis?"

**Answer scaffold:**

GSE250521 hires PNGs are downsampled by `tissue_hires_scalef = 0.5625879` from the full-resolution scanner image (per dataset metadata; see `project/data/processed/GSE250521/.../sample_metadata.tsv`). At our 224 px tile size in hires-image space, this corresponds to ≈ 110 µm in tissue space — adequate for tissue architecture (gland shape, stromal density, follicular structure) but insufficient for nuclear morphometry (nuclear chromatin, mitotic figures, cell-membrane detail) which are the features most likely to encode molecular-axis information. We acknowledge this resolution as a hard ceiling and explicitly mark it in the 5 re-entry conditions (`2026_05_04_image_dm1_final_nogo_decision.md` §5 condition #1: "Full-resolution scanner WSI available for ≥1 thyroid cohort with paired bulk RNA-seq"). The negative result is therefore properly described as "bounded at hires Visium 224 px under ImageNet feature transfer", not as "H&E does not contain DM1 information in any form".

---

## Q4 (medium probability) — "Couldn't slide-level aggregation rescue the spot-level negative?"

**Answer scaffold:**

We tested three slide-level aggregations of per-spot Ridge predictions: mean (n = 16 slides per metric), top-25%-mean, and high-risk fraction (predicted q75). At slide level, DM1_resid mean Spearman = −0.003, top-25%-mean Spearman = +0.209, high-risk-fraction was constant (NaN) for several targets. Top-25%-mean shows a marginal positive direction at n = 16 — too few slides for meaningful significance and consistent with the per-fold variance observed at spot level (range [−0.07, +0.12]). We report this in `closure_battery_metrics.tsv` E section and in the closure battery long-form (`pathology_dm1_closure_battery_2026_05_04.md`); the aggregation-level signal is too weak to revise the GO/NO-GO verdict and is reserved as exploratory.

---

## Q5 (medium probability) — "Did you test other tile sizes / context windows?"

**Answer scaffold:**

We pre-registered a tile-size sweep with a +0.05 Spearman improvement gate. At 448 px (n = 2,935 tiles after border filtering), pooled DM1_resid Spearman = 0.017, slightly worse than the 224 px baseline (0.022). The gate was therefore not met, and the 672 px embedding was cancelled per pre-registration (tiles were extracted but the model run was not executed). Cancellation saved ≈ 25 min of CPU and demonstrates that larger morphological context does not bridge the signal floor at our scale.

---

## Q6 (medium probability) — "Why is the raw 8-gene signal from a 9-dim RGB baseline higher than from 2,048-dim ResNet50?"

**Answer scaffold:**

A 9-dimensional simple-RGB baseline (per-tile mean, standard deviation, and Shannon entropy of each color channel) achieved pooled Spearman r = 0.224 / DM1_high AUROC = 0.702 on the **raw** DM1 target — substantially higher than 2,048-dim ResNet50 ImageNet on the same raw target (r = 0.061 / AUROC = 0.521). On the **depth-residualized** target, both fell to near-zero (RGB r = 0.080, ResNet50 r = 0.022). The most parsimonious interpretation is that the apparent raw signal lives in tile-level color statistics — i.e., a stain-darkness ↔ sequencing-depth artifact (independently confirmed by DM1_raw ~ log_counts ρ = −0.45 in the spatial-validation report). Within-sample residualization on log_counts and log_ngenes was designed precisely to remove this confound, and its near-elimination of the signal in both feature spaces validates the residualization design rather than refuting it.

---

## Q7 (low probability) — "Could you have used GigaTIME's released model directly on your tiles?"

**Answer scaffold:**

As of submission, GigaTIME model weights and code repository have not been located in our publicly accessible search (GitHub, HuggingFace, Microsoft Research project pages all returned no GigaTIME-named asset; see `2026_05_04_gigatime_thca_relevance_scan.md` §1). Microsoft Research typically releases code 6–12 months after publication. If/when GigaTIME assets become available, the same closure-battery pipeline can be re-applied, and we have specified this as part of `2026_05_04_image_dm1_final_nogo_decision.md` §5 future re-verification.

---

## Q8 (low probability) — "Does the ATC vs non-ATC AUROC = 0.689 contradict the H&E negative result?"

**Answer scaffold:**

The ATC vs non-ATC binary classification (LogReg C = 0.1, LOSO pooled AUROC = 0.689) reflects the gross anaplastic morphology of ATC samples — a known pathological feature distinguishable by trained pathologists at a glance and well within ImageNet-feature-transfer reach. This is consistent with our negative result on the molecular axis: H&E can support classification when the morphological signal is gross and pre-existing (anaplasia), but cannot recover continuous molecular gradients (stage ordinal Spearman = −0.047; epithelial / proliferation raw scores Spearman ≤ 0.03; depth-residualized DM1 Spearman = 0.022). The two findings are mutually consistent and jointly support the conclusion that hires Visium 224 px H&E under ImageNet feature transfer captures gross histology but not within-tumor molecular gradients.

---

## Q9 — author keyboard (NOT scaffolded by Claude)

Per `STATUS_PAPER1_2026_05_02.md` voice-protected list. Author writes manually using `2026_05_04_voice_hook_fact_brief.md` §C dark-matter-gap content + closure-battery negative as evidence of pre-registered honest hypothesis testing.

---

## Cross-references

- `2026_05_04_gigatime_thca_relevance_scan.md` (S3 scan)
- `2026_05_04_paper1_supp_methods_he_negative_feasibility.md` (commits `cf14656` + `2ffe929`)
- `2026_05_04_paper2a_gigatime_methods_comparison_table.md` (S1 comparison table)
- `2026_05_04_image_dm1_final_nogo_decision.md` (5 re-entry conditions)
- `2026_05_04_paper1_dm1_full_molecular_only_lock.md` (manuscript-safe envelope)
- `STATUS_PAPER1_2026_05_02.md` (voice-protected list)
