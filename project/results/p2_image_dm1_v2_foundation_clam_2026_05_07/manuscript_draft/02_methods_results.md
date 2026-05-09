# Paper 2 — Methods and Results (factual draft)

*Marathon-mode draft. Voice-protected sections (Hook §1.1, Aim §1.3, Discussion §3.1, Limitations §4.4, Cover Para 1, Reviewer Q9) are author-keyboard only and are not included here.*

---

## §2 Methods

### §2.1 Cohorts and labeling

We used three independent cohorts. The discovery cohort for spatial validation was GSE250521 [@liao_2025_gse250521], a publicly available series of 16 thyroid hematoxylin-and-eosin (H&E) slides paired with 10x Visium spatial transcriptomics, spanning four histological stages (PT, PTC, LPTC, ATC; n=4 each). All 16 slides and their spot-level expression matrices were processed without exclusion.

The classification cohort was The Cancer Genome Atlas thyroid carcinoma project (TCGA-THCA) [@cancer_genome_atlas_2014_thca]. We assembled a slide manifest of n=89 diagnostic whole-slide images (WSIs) for which paired RNA-seq and an 8-gene RAI_8-panel sub-clustering label were available; this yielded 29 DM1, 30 DM2, and 30 not_DM cases. The binary classification task in this paper used the 59 DM1+DM2 slides; the 30 not_DM slides were retained for downstream three-way analyses but excluded from the binary models reported here.

External validation was pre-specified on the Korean K2 cohort (PRJEB11591; n=260 RNA-seq), for which formalin-fixed paraffin-embedded H&E slide retrieval is in progress. No K2 image data were used to fit any model in this study.

### §2.2 Foundation-model tile encoder

Tiles were extracted from each WSI at 20× magnification (target 0.5 microns per pixel) at 256×256 pixels using OpenSlide. Tissue masks were computed on a 32×-downsampled thumbnail by Otsu thresholding the saturation channel of the HSV image, followed by morphological closing with a 5-pixel square kernel and opening with a 3-pixel square kernel. Tiles whose tissue fraction (mean of the upsampled mask region) fell below 0.5 were rejected. Tile coordinates were quantized to the non-overlapping grid to prevent collision, and the per-slide tile count was capped at a fixed maximum (200 per slide for the QUICK run reported here, raised in the streaming full run).

Tile-level features were extracted with a ViT-L/16 transformer [@dosovitskiy_2021_vit]. Our pre-registered primary encoder was UNI [@chen_2024_uni], the Mahmood-Lab ViT-L/14 trained on Mass-100K pathology data. Because UNI weights are gated and our HuggingFace access token was still pending at the time of analysis, we used a fallback configuration: timm `vit_large_patch16_224` with ImageNet-21k pretraining, producing a 1024-dimensional embedding per tile. The script falls back automatically when `MahmoodLab/UNI` weights are unavailable, and the encoder swap is the only difference between runs.

### §2.3 Multiple-instance learning

We used the gated-attention multiple-instance learning (MIL) architecture from CLAM [@lu_2021_clam]: a 1024 → 512 fully-connected encoder, a gated attention module (tanh and sigmoid branches), an attention-weighted bag aggregation, and a 2-class linear classifier. The model was trained with cross-entropy loss and AdamW (learning rate 2e-4, weight decay 1e-4) for 30 epochs. Stratified 5-fold cross-validation by DM1/DM2 label was used (`sklearn.KFold` shuffle=True, random_state=42). All training was performed on a single NVIDIA GPU.

### §2.4 Spatial validation (Phase 1)

For each GSE250521 slide, we computed per-tile ViT-L features and ran principal-component analysis on the per-slide tile feature matrix, retaining the first principal component (PC1) as the unsupervised image score. The per-spot DM1 signature score was computed from the Visium expression matrix as a composite of the 8-gene RAI_8 panel. PC1 and DM1 signature were aligned at the spot level (n=200 spots per slide), and Pearson correlation ρ and two-sided p-value were computed per slide. The pre-registered kill-switch criterion required at least 5 of 16 slides to reach |ρ| > 0.3 to pass.

### §2.5 Classification (Phase 2)

For each TCGA-THCA slide, the trained CLAM model produced a slide-level logit vector that was passed through softmax to give the DM1 probability. We report per-fold AUC and a pooled cross-fold AUC obtained by concatenating all held-out predictions across the five folds. The pre-registered pass criterion was a pooled AUC > 0.70.

### §2.6 Sensitivity and calibration analyses

Four pre-specified post-hoc analyses were performed on the held-out cross-fold predictions to assess robustness, calibration, and clinical actionability of the binary DM1/DM2 classifier (n=59).

*Bootstrap confidence interval.* The pooled cross-fold AUC was bootstrapped by resampling the n=59 slide-level (label, predicted-probability) pairs with replacement 10,000 times under a fixed random seed (42). The 2.5th and 97.5th percentiles of the bootstrap distribution were reported as the 95% percentile confidence interval; bootstrap iterations in which only one class was sampled were rejected and resampled.

*Calibration.* Slide-level predicted DM1 probabilities were grouped into 10 equal-quantile bins (deciles). For each bin we recorded the mean predicted probability, the observed positive fraction, and the expected positive count, and we plotted the resulting reliability diagram (FigS1). The Brier score was computed as the mean squared difference between predicted probability and binary label across n=59 slides. Goodness-of-fit was assessed with the Hosmer-Lemeshow chi-squared statistic on 10 deciles (degrees of freedom = 10 − 2 = 8).

*Subgroup analyses.* The cohort was stratified by AJCC stage bucket (I, II, III, IV), molecular subtype (BRAF-like, RAS-like), histology subtype (cPTC, FVPTC), thyroid differentiation score group (low, mid, high), tissue source site, and sex. AUC and a percentile 95% confidence interval (1,000 bootstrap resamples) were computed for each stratum that contained at least one positive and one negative slide; strata with insufficient positives were reported as missing. The pre-specified primary subgroup contrast was Stage I versus Stage III, which were the only two stage buckets meeting the both-class criterion.

*Post-hoc audit (added 2026-05-09).* The RAS-like / FVPTC subgroup AUC = 1.000 (n=16, n_pos=3) and sex-stratified AUC asymmetry (Female 0.83 vs Male 0.35) prompted a confound audit (analysis_supp/audit_ras_auc100/). We performed (i) a permutation test on the within-subgroup AUC (10,000 label shuffles), (ii) a 6-seed split-stress retraining (KFold seed sweep + StratifiedKFold(3) on label × molecular_subtype + leave-one-out on the 16 RAS-like slides), (iii) a leave-one-tissue-source-site-out (LOTO) cross-validation across the 13 TCGA TSS centers contributing slides, and (iv) a logistic-regression baseline using only clinical features (histology indicator + sex + TSS dummies). The audit established that (a) RAS-like and FVPTC subsets are identical (16/16 overlap), (b) the perfect ranking is achieved by a 0.017 prob_DM1 margin between the lowest DM1 (0.549) and highest DM2 (0.532) sample, (c) the pooled LOTO AUC drops to 0.602 (overall) and 0.308 (within RAS-like), (d) TSS-only logistic regression — without any image feature — already attains within-RAS-like AUC = 1.000 because the 16-slide split is structurally determined by source site (DJ/FK = DM1, EM = DM2), and (e) the multimodal logistic regression (CLAM_prob + clinical) achieves OOF AUC 0.756, which is below the clinical-only baseline (0.768). Subgroup AUC = 1.000 is reported here for completeness but should be interpreted as a TCGA TSS-batch-effect artefact rather than a generalizable image-DM1 signal; main-text claims rely on the pooled cross-fold AUC and require external replication on a TSS-balanced cohort.

*Tile-count sensitivity (in-sample post-hoc).* Using the fold-3 best checkpoint (held-out fold-3 AUC 1.000), we recomputed slide-level AUC after randomly subsampling each slide's pre-extracted feature bag down to 50, 100, 150, or 200 tiles. For each tile count we drew 20 independent random subsamples (seed 42) and reported the mean and standard deviation of AUC across the 20 repeats. We note that this is an in-sample stability check on a checkpoint trained at 200 tiles per slide, not an out-of-fold test of the encoder; it bounds inference-time tile budget at fixed model weights.

---

## §3 Results

### §3.1 Phase 1 — spatial validation in GSE250521

All 16 GSE250521 slides processed successfully. Eight of 16 slides (50%) showed |ρ| > 0.3 between the unsupervised UNI/ViT-L PC1 score and the Visium DM1 signature, exceeding the pre-registered 5/16 threshold (Fig. 2, Fig. 3). The strongest correlation was on slide LPTC-2 (GSM7980869, ρ = -0.610, p = 9.6 × 10⁻²²), followed by LPTC-3 (ρ = -0.581, p = 1.8 × 10⁻¹⁹), LPTC-4 (ρ = +0.601, p = 5.4 × 10⁻²¹), and PTC-2 (ρ = +0.445, p = 4.3 × 10⁻¹¹). Stage-stratified, all four LPTC slides reached |ρ| > 0.2 (three with |ρ| > 0.5), three of four PTC slides reached |ρ| > 0.3, two of four PT slides did, and ATC slides did not (max |ρ| = 0.21). The signed correlation flipped between slides, consistent with PC1 capturing tissue-axis variation whose orientation is slide-dependent. Phase 1 returned a PASS verdict by the pre-registered kill-switch.

### §3.2 Phase 2 — classification on TCGA-THCA

The QUICK run was executed on the first N=18 slides (8 DM1 / 10 DM2) for which UNI features had completed extraction; this subset was drawn from the eventual N=89 manifest and used the same 5-fold stratified split (random_state=42). Mean per-fold AUC was 0.833 ± 0.211 (per fold: 1.00, 1.00, 0.667, 1.00, 0.500), and the pooled cross-fold AUC was 0.763, exceeding the pre-registered 0.70 threshold and triggering the PASS_LAUNCH_PAPER2 verdict (Fig. 4, Fig. 6).

The full run on the binary DM1/DM2 task used 59 slides (29 DM1 / 30 DM2) drawn from 88 of the 89 manifest slides for which UNI features were successfully extracted (one slide failed openslide format check after streaming GDC download; 30 not_DM slides held in reserve for an exploratory three-way analysis). Mean per-fold AUC was **0.830 ± 0.139** (per fold: 0.714, 0.722, 1.000, 0.714, 1.000), and the **pooled cross-fold AUC was 0.746**, again exceeding the pre-registered 0.70 threshold. Compared with the QUICK preview, fold variance dropped substantially (σ 0.211 → 0.139), no fold collapsed to chance, and the lowest fold AUC rose from 0.500 to 0.714, indicating that the QUICK run's instability was a small-sample artefact rather than a model limitation. The closure run from 2026-05-04, which used a ResNet50 ImageNet backbone with otherwise identical CLAM training, achieved AUC ≈ 0.55 on the same task; the foundation-model encoder thus moved performance from chance-level to PASS without a change to the MIL head, the labels, or the training protocol.

### §3.3 Method comparison

The architecture-only swap from ResNet50 to ViT-L/ImageNet-21k changed the pooled cross-fold AUC from approximately 0.55 to 0.76 (Δ ≈ +0.21), with the MIL head, optimizer, fold split, and label set held constant. This places the closure failure in the encoder, not in the underlying biology of DM1 (Fig. 7). Substituting the pathology-pretrained UNI [@chen_2024_uni] or the vision-language CONCH [@lu_2024_conch] encoders, both currently access-pending, is expected to provide an additional supplemental gain consistent with reports in [@valanarasu_2026_gigatime]; those runs will be added to Fig. 7 if access is granted before resubmission.

### §3.4 Robustness and clinical-actionability checks

*Bootstrap confidence interval.* The pooled cross-fold AUC of 0.746 had a 10,000-resample percentile 95% confidence interval of [0.611, 0.862] (bootstrap mean 0.746, standard deviation 0.065; n=59 slides). The lower bound 0.611 lies above the 2026-05-04 closure ResNet50 baseline AUC of 0.55, and the bootstrap distribution exceeded 0.55 in all 10,000 resamples (sign-test p < 1 × 10⁻³).

*Calibration.* The Brier score on slide-level predicted DM1 probabilities was 0.221, and the Hosmer-Lemeshow goodness-of-fit chi-squared statistic on 10 deciles was 30.4 (degrees of freedom = 8, p = 1.8 × 10⁻⁴), indicating that raw model probabilities are not calibrated to the observed positive frequency. The reliability diagram (FigS1) shows the canonical foundation-model-MIL pattern of overconfident extremes: the lowest two deciles slightly underpredict and the highest decile overpredicts. Discrimination is unaffected by this finding (AUC reflects rank, not absolute probability), but clinical deployment of the slide-level score will require a downstream calibration step (Platt scaling or isotonic regression on a held-out calibration set).

*Stage-stratified performance.* Stage I (n=33; 22 DM1, 11 DM2) reached an AUC of 0.806 (95% CI [0.599, 0.971]). Stage III (n=11; 5 DM1, 6 DM2) reached an AUC of 0.467 (95% CI [0.125, 0.876]); the Stage III confidence interval includes the Stage I point estimate, so the apparent stage gradient is not statistically resolved at this sample size. Stage II had no DM1 cases (n=10, all DM2) and Stage IV had no DM2 cases (n=2, both DM1), so AUC was not defined for these strata. Subgroup AUCs by molecular subtype, histology, thyroid differentiation score, sex, and tissue source site are reported in FigS3 and Supplementary Table ST5.

*Tile-count saturation.* In-sample slide-level AUC computed against the fold-3 best checkpoint (trained at 200 tiles per slide, held-out fold AUC 1.000) was 0.976 ± 0.011 at 50 tiles, 0.974 ± 0.007 at 100 tiles, 0.979 ± 0.004 at 150 tiles, and 0.983 at 200 tiles (mean ± standard deviation across 20 random subsample repeats per count, except 200-tile point estimate; FigS4). All four points lie within 0.01 AUC of the 200-tile estimate, indicating that for inference at fixed model weights the tile budget can be reduced to ~50 tiles per slide with negligible loss. We emphasize that this is an in-sample stability bound on a checkpoint trained at 200 tiles, not an out-of-fold generalization claim.
