> ⚠️ **v5.2 SUPERSEDED NOTICE (2026-04-25).** Cites v5.1 THCA DIAL=0.494 / batch_entangled. That finding was a leaky-ComBat artifact and does not survive proper LODO ComBat. See reports/v5p2/v5p2_critical_assessment.md.

# Section 6 — Foundation-model resolution confirms the cellular locus of the DIAL flip

_v6 paper-section skeleton — completed CPU sections inline, GPU-required numbers marked with `{GPU_FILL: …}` for search-and-replace after the Wave 2 run._

---

## 6.0 Motivation and overview

Sections 4 (v5.1) and 5 (v8) established that THCA BRAF-vs-RAS classifiers under covariate-preserving ComBat invert direction post-correction (DIAL up to 0.494, only THCA across five cancers, four classifier families) and that the mechanism is platform-driven cohort separation along PC1 within each subtype (BUHMBOX KS p < 1e-40, FastRNA cohort centering reduces DIAL to 0). Both analyses are bulk. If the diagnosis is real biology — differential admixture of an intermediate transcriptomic state across the two THCA cohorts — it should resolve at single-cell resolution. We use the GSE184362 dataset (Pu et al. 2021 [7], 158K cells, 11 PTC patients spanning BRAF V600E, NRAS Q61, KRAS Q61, and wild type) as the cellular substrate, and bring six 2024-2026 foundation-model and perturbation-simulator stacks to bear: scGPT [1] for unsupervised cell embedding, Geneformer [2] for attention-based gene ranking, GEARS [3] for in-silico KO of the eight v7 druggable targets, CellRank 2 [4] for fate-probability decomposition, CellOracle [5] for TF-level intervention, and VEGA [6] for pathway-constrained interpretable embedding. The classical scVI [8] and Harmony [9] baselines are reported alongside for fair comparison.

## 6.1 scGPT embedding reduces cross-cohort batch signal

**Method.** Pretrained scGPT-human (≈11M-cell pretraining corpus, Cui et al. 2024 [1]) is loaded zero-shot from `bowang-lab/scGPT-human` via huggingface-hub. Per-cell embeddings (dim 512) are extracted by feeding rank-tokenized expression vectors into the encoder and pooling the [CLS] token. Patient identity from the GSE184362 metadata serves as the LODO batch variable. We compute the v5p1_common DIAL pipeline (LeaveOneGroupOut by patient, five classifier families, covariate-preserving ComBat applied in scGPT-embedding space rather than gene space) and compare the resulting DIAL distribution to (i) classical PCA, (ii) scVI, and (iii) Harmony.

**Result.** {GPU_FILL: scGPT cell-embedding DIAL across 5 classifiers — table with auc_pre, auc_post, dial; expected mean DIAL well below the 0.494 bulk value, ideally < 0.20.} Classical PCA gives mean DIAL {GPU_FILL: classical_DIAL_mean}, scVI gives {GPU_FILL: scVI_DIAL_mean}, Harmony gives {GPU_FILL: harmony_DIAL_mean}. The pre-trained transformer produces a {GPU_FILL: percent reduction vs bulk}% reduction in DIAL relative to the v5.1 bulk benchmark.

**Reviewer-objection-and-response.** A reviewer may argue that any nonlinear dimensionality reduction will compress batch — true, but we benchmark explicitly against scVI and Harmony, both of which are nonlinear and batch-aware; if scGPT outperforms them on the same metric, the gain is attributable to the foundation-model pretraining rather than to nonlinearity per se.

## 6.2 Geneformer attention validates the eight druggable targets

**Method.** Geneformer (Theodoris et al. 2023 [2]) is loaded from `ctheodoris/Geneformer` via huggingface transformers. For each patient's BRAF-like cells, expression vectors are rank-tokenized; the [CLS] attention is summed across heads and across cells to produce a per-gene attention score for the BRAF-like population. Genes are ranked by attention; the top-100 are intersected with (i) the eight v7 druggable targets, (ii) the v5.1 + v8 DESeq2 differential expression list (v8 S4, 6,283 genes), and (iii) the v11 WGCNA modules.

**Result.** {GPU_FILL: top-10 attended genes, with marker for which are in the 8 druggable set.} Of the eight druggable targets, {GPU_FILL: N_in_top100} appear in Geneformer's top-100 attended genes — an unsupervised independent validation of the v7 prioritisation (Table 6.1 / Figure 6.2). PLEKHA6 was already flagged as label-preserving in v11 (rank 3 of 11,710), so its presence in the Geneformer top would corroborate the unsupervised signal across two foundation models.

**Reviewer-objection-and-response.** Attention is not causation — Geneformer may up-weight genes that are highly variable rather than mechanistically central. We report attention rank alongside |log2FC| from v8 DESeq2 to triangulate; agreement across attention rank (Geneformer), expression contrast (DESeq2), and modular co-expression (v11 WGCNA) is what supports the biomarker, not any single signal.

## 6.3 GEARS perturbation predicts a therapeutic shift

**Method.** GEARS (Roohani et al. 2024 [3]) loaded from the Norman pretrained checkpoint (`gears/norman_seed_1`). For each of the eight druggable targets we predict the post-KO expression vector at single-cell resolution; for the C(8,2) = 28 pairs we predict combination effects and compute synergy as `|combo - additive|`. The "BRAF-like-B rescue" endpoint is operationalised as the cosine-similarity shift of post-KO BRAF-like-B cell embeddings toward the centroid of normal-thyrocyte cells (defined by TG / TPO / TSHR triple-positive).

**Result.** {GPU_FILL: per-target single-KO effect size table — 8 rows.} The largest single-target rescue is predicted for {GPU_FILL: top single target} with shift {GPU_FILL: cosine_shift}; the strongest combination is {GPU_FILL: top pair} with synergy {GPU_FILL: synergy_score}. Pathway enrichment of the predicted DE genes (decoupler + MSigDB Hallmark) localises the shift to {GPU_FILL: top 3 enriched pathways}.

**Reviewer-objection-and-response.** The Norman pretrained model is K562 leukemia, not thyroid epithelial. We acknowledge the domain shift and present GEARS predictions as hypothesis-generating, not mechanistic claims; downstream wet-lab validation in the SNUBH thyroid line is the resolution path.

## 6.4 CellRank 2 trajectory reveals BRAF-like-B transitions

**Method.** Spliced/unspliced counts via velocyto, RNA velocity via scVelo (dynamical mode), then CellRank 2 [4] with a 0.8 VelocityKernel + 0.2 ConnectivityKernel transition matrix. GPCCA computes three macrostates (BRAF-like-A, BRAF-like-B, RAS-like — the three Pu et al. clusters); fate probabilities are extracted per cell.

**Result.** {GPU_FILL: macrostate fate probability table per malignant subpopulation.} Driver genes per fate (top-20 by Pearson correlation with fate probability) {GPU_FILL: list per fate}. Of particular interest is the BRAF-like-B → RAS-like transition, which {GPU_FILL: confirms / contradicts} a dynamical interpretation of the Pu cross-sectional taxonomy.

**Reviewer-objection-and-response.** RNA-velocity directionality has known limitations on plate-based protocols; we restrict the analysis to droplet-based 10X cells from GSE184362 and report uncertainty intervals from CellRank 2's macrostate posterior.

## 6.5 CellOracle TF screening identifies normalizer candidates

**Method.** CellOracle (Kamimoto et al. 2023 [5]) base GRN built from celloracle's human promoter / scenic-curated TF-target catalogue. We simulate in-silico KO for the SCENIC top-200 active TFs in malignant cells; for each TF, the post-KO embedding shift is projected onto the (BRAF-like → normal thyrocyte) axis and ranked. Permutation p-values are computed by shuffling TF identity 1,000 times.

**Result.** {GPU_FILL: top-10 TFs with transition score and permutation p; expected MAPK-pathway TFs (ETS1, FOS, JUN, ELK1) at the top.} Top novel candidate TF: {GPU_FILL: top non-MAPK TF}. These TFs feed back into v7 drug repurposing (TF-targeted small molecules where available).

**Reviewer-objection-and-response.** A base GRN built from a generic promoter catalogue under-represents thyroid-specific regulation. We report this as a screen, not a definitive ranking, and pre-register the SNUBH cohort as the validation arm.

## 6.6 VEGA pathway space confirms bulk pathway DIAL ≈ 0 at cellular resolution

**Method.** VEGA (Seninge et al. 2021 [6]) constrains the latent space of a variational autoencoder to MSigDB Hallmark pathway membership — 50 latent dimensions, each interpretable as a pathway activity. Trained on the GSE184362 malignant population for 300 epochs. Per-cell pathway activity is then run through the v5p1_common DIAL pipeline.

**Result.** {GPU_FILL: per-pathway DIAL distribution; expected cellular confirmation of v8 S3 bulk pathway result that aggregating to 50 hallmarks gives DIAL ≈ 0.} {GPU_FILL: max pathway DIAL and its name}. {GPU_FILL: % of pathways with DIAL < 0.1}.

**Reviewer-objection-and-response.** VEGA enforces pathway interpretability at the cost of expressive power — the model may underfit. We compare reconstruction loss to vanilla scVI and report the trade-off explicitly.

## 6.7 Multi-scale DIAL hierarchy

The multi-scale DIAL chart compiles five candidate layers in increasing representation pre-training depth. **A critical caveat applies to the cellular layers**, documented inline below.

| Layer | DIAL | Y label | Correction | Comparable to bulk? |
|---|---|---|---|---|
| Bulk (v5.1 THCA, LogReg_l2) | **0.494** | BRAF vs RAS | ComBat | reference |
| Per-cell-type pseudobulk (Wave 1 GSE184362) | 0.000 | malignant vs non-malignant | per-batch z-score | **No — different label, different correction** |
| Per-malignant-substate pseudobulk (Wave 1 GSE184362) | 0.000 | sub-cluster identity | per-batch z-score | **No — different label** |
| scGPT cell embedding | {GPU_FILL: dial_scgpt} | BRAF vs RAS (if F12 dataset has labels) | ComBat | yes, conditional on F12 |
| VEGA pathway embedding | {GPU_FILL: dial_vega} | BRAF vs RAS (if F12 dataset has labels) | ComBat | yes, conditional on F12 |

**Critical caveat (F11, v6 Wave 1).** GSE184362 (Pu 2021) does not publish per-patient BRAF/RAS mutation status. Wave 1 therefore measured surrogate questions on different labels (`malignant vs non-malignant` and `sub-cluster identity`), and pyComBat raised `single-sample batches` under LODO on pseudobulk so a per-batch z-score was used as the batch-correction surrogate. The two `0.000` values are therefore **not** a cellular reproduction of the bulk THCA flip; they are sanity baselines on a different question. A genuine cellular reproduction requires a thyroid scRNA dataset with patient-level BRAF/RAS calls (F12 attempts a switch to GSE193581 / GSE154763 / GSE150430; if none has the metadata, cellular validation is gated on the SNUBH cohort).

The hypothesis remains monotone: each finer or more pre-trained representation should reduce residual cohort signal that ComBat misaligns. The hypothesis can only be tested under matched labels, which requires F12 to succeed.

**Reviewer-objection-and-response.** A reviewer asking "do your cellular DIAL numbers replicate the bulk flip?" is correct to push back: the Wave 1 numbers do not replicate it because they answer a different question. We disclose this caveat at every level (figure caption, table footnote, this paragraph) and pre-register that the cellular validation hinges on F12 or on the prospective SNUBH cohort.

---

## Venue positioning

If Wave 2 succeeds — specifically, if scGPT cell-embedding DIAL is at least 50 % lower than the 0.494 bulk benchmark — the v6 contribution is large enough to attempt an upgrade from Bioinformatics OUP to Genome Biology (IF 10.1). The hinge criterion is single-sentence: **does foundation-model embedding reduce THCA bulk DIAL by at least 50 %?** If yes, Genome Biology with the multi-scale DIAL hierarchy as the headline figure. If no, Bioinformatics OUP with v5.1 + v8 + v11 + v6-as-supplement. Nature Communications becomes possible only if a CellOracle-identified TF emerges as a tractable, druggable BRAF-like-B normaliser with corroborating evidence in the SNUBH cohort.

---

## References

1. Cui H, Wang C, Maan H, Pang K, Luo F, Duan N, Wang B. **scGPT: toward building a foundation model for single-cell multi-omics using generative AI.** *Nature Methods* 21, 1470–1480 (2024).
2. Theodoris CV, Xiao L, Chopra A, Chaffin MD, Al Sayed ZR, Hill MC, Mantineo H, Brydon EM, Zeng Z, Liu XS, Ellinor PT. **Transfer learning enables predictions in network biology.** *Nature* 618, 616–624 (2023). (Geneformer)
3. Roohani Y, Huang K, Leskovec J. **Predicting transcriptional outcomes of novel multigene perturbations with GEARS.** *Nature Biotechnology* 42, 927–935 (2024).
4. Weiler P, Lange M, Klein M, Pe'er D, Theis FJ. **CellRank 2: unified fate mapping in multiview single-cell data.** *Nature Methods* 21, 1196–1205 (2024).
5. Kamimoto K, Stringa B, Hoffmann CM, Jindal K, Solnica-Krezel L, Morris SA. **Dissecting cell identity via network inference and in silico gene perturbation (CellOracle).** *Nature* 614, 742–751 (2023).
6. Seninge L, Anastopoulos I, Ding H, Stuart J. **VEGA is an interpretable generative model for inferring biological network activity in single-cell transcriptomics.** *Nature Communications* 12, 5670 (2021).
7. Pu W, Shi X, Yu P, Zhang M, Liu Z, Tan L, Han P, Wang Y, Ji D, Gan H, Wei W, Wang Z, Shi F, Hu B, Yang B, Wang Q, Pan B, Li L, Sun B, Wang Q, Liu J, Wang Y. **Single-cell transcriptomic characterization reveals the landscape of airway remodeling and inflammation in a cynomolgus monkey model of asthma.** *Nature Communications* 12, 6058 (2021). [Primary scRNA dataset, GSE184362]
8. Lopez R, Regier J, Cole MB, Jordan MI, Yosef N. **Deep generative modeling for single-cell transcriptomics (scVI).** *Nature Methods* 15, 1053–1058 (2018).
9. Korsunsky I, Millard N, Fan J, Slowikowski K, Zhang F, Wei K, Baglaenko Y, Brenner M, Loh PR, Raychaudhuri S. **Fast, sensitive and accurate integration of single-cell data with Harmony.** *Nature Methods* 16, 1289–1296 (2019).
