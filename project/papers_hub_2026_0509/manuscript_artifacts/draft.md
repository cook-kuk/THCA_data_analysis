# A leakage-aware post-TESLA neoantigen vaccine benchmark reveals source bias and enables HLA-aware candidate prioritization

**Authors**: Seungho Cook *et al.*
**Corresponding**: kukshomr@gmail.com
**Date**: 2026-05-07
**Status**: DRAFT — under internal review.

---

## Abstract

**Background.** Neoantigen vaccine design increasingly depends on computational immunogenicity prediction, but post-TESLA datasets differ in source composition, HLA distribution, curation practice, and leakage risk relative to the corpora used to train commonly cited predictors. The widely-reported benchmark performance of biophysics- and protein-language-model-based predictors has not been systematically tested under strict cross-source generalization.

**Methods.** We integrated 13 post-TESLA neoantigen, T-cell, TCR, and vaccine-related resources into a leakage-aware benchmark with explicit `test_set_safety` annotation. We systematically evaluated peptide and HLA representations, including biophysical features, ESM2-35M protein-language-model embeddings, and shortcut baselines, under five split designs: random stratified k-fold, peptide GroupKFold, peptide-HLA GroupKFold, Leave-One-Source-Out (LOSO), and Leave-One-HLA-Out. We performed shortcut tests (source-only / HLA-only / label-shuffle within-source-or-HLA) and HLA-stratified feature confounding analysis. Candidate prioritization was re-run under the strict-generalization model with explicit autoimmunity-self-similarity safety flags.

**Results.** In the current verified analysis the cleaned benchmark contains **84,549 labeled records** spanning **46,731 unique peptides** and **235 HLA alleles** across **11 distinct curation sources** (after dedup; 119,237 raw master records). Random 5-fold cross-validation gave Logistic-regression AUROC ≈ 0.96 for biophysics+source+HLA features but the same model under Leave-One-Source-Out collapsed to AUROC ≈ 0.48–0.58, near or below chance. Source-only and HLA-only shortcut classifiers reached random-CV AUROC of 0.926 and 0.745 respectively but dropped to 0.500 and 0.517 under LOSO — direct evidence of source/HLA memorization. ESM2-35M mean-pooled embeddings were the only representation that maintained meaningful cross-source signal under LOSO (AUROC 0.63–0.74 across LR/RF), partially rescuing the worst-affected source (NEPdb leave-out 0.40 → 0.48). Per-HLA stratified analysis preserved the aromatic-fraction effect direction but attenuated or sign-flipped several other pooled biophysical signals, consistent with HLA-composition confounding. Candidate re-prioritization under the strict-generalization Logistic-regression model preserved KRAS G12D peptide families (top: VVVGADGVG, score 0.69) but with reduced absolute scores compared to leakage-prone training.

**Conclusions.** Our results indicate that conventional random cross-validation substantially overestimates neoantigen immunogenicity prediction performance in the presence of source and HLA-composition bias. Representation choice and evaluation design dominate over architecture for strict generalization. We provide a leakage-aware benchmark resource, an HLA-aware candidate prioritization layer with explicit autoimmunity flags, and a preclinical validation roadmap. We do not claim clinical efficacy, validated immunogenicity, or production-grade pMHC structural prediction.

---

## 1 · Introduction

The TESLA consortium (Wells *et al.*, *Cell* 2020) catalyzed open benchmarking of neoantigen immunogenicity predictors. Subsequent resources — IEDB, NEPdb, IMPROVE, McPAS-TCR, Neodb, dbPepNeo, TSNAdb, CEDAR, BigMHC, NeoRanking, VDJdb — have grown rapidly but vary in scope, curation, HLA representation, and the assays used to assign immunogenicity labels. Many published predictors (NetMHCpan, MHCflurry, PRIME, BigMHC) train on IEDB, so benchmarking on IEDB-derived peptides causes test-set contamination.

Two underlying problems remain insufficiently addressed:

1. **Evaluation design.** Random k-fold cross-validation on a multi-source corpus implicitly assumes independence between training and test peptides at the source/HLA level. When source curators differ in their positive/negative ratios, in their HLA representation, or in their assay choice, random k-fold can be solved by source memorization.

2. **HLA-composition confounding.** Biochemical features that appear to discriminate immunogenic from non-immunogenic peptides at the pooled level can be Simpson's-paradox artifacts of HLA-distribution imbalance. Without per-HLA stratification, attribution methods (SHAP, permutation importance) can highlight features that are biologically uninterpretable.

In this work we ask: under the strictest evaluation protocols, what actually determines neoantigen immunogenicity prediction performance? Inspired by recent systematic-benchmark work in proteomics-based drug-response prediction (which showed that pretrained representation quality, not architecture or fusion strategy, explained statistical improvements once strict splits were enforced), we translate that philosophy to neoantigen prediction: split by source/peptide/HLA, audit shortcut signals, stratify feature interpretation by HLA, and ask whether protein-language-model representations rescue cross-source generalization.

---

## 2 · Results

### 2.1 · A leakage-aware 13-source benchmark with explicit test-set safety

We integrated 13 post-TESLA neoantigen / T-cell / TCR / vaccine resources (Methods). After QC, normalization, and dedup, the cleaned master contained **119,237** records, of which **84,549** carried a clear immunogenic-positive or immunogenic-negative label (Figure 1, Table 1). The benchmark covered **46,731 unique peptides**, **235 normalized HLA alleles**, and **11 distinct sources** (after merging dataset variants).

We assigned each record a `test_set_safety` flag based on per-source rules (Methods, Table S2). **TRAINING_OVERLAP** records (n=61,656) — those known or strongly suspected to be in the training corpora of commonly used predictors (IEDB v3 T-cell, BigMHC el_train) — were excluded from benchmark evaluation. **HELD_OUT** records (n=15,765, primarily NEPdb), **EXTERNAL_TEST** records (n=1,197, primarily TESLA), and **PARTIAL_OVERLAP** records (n=4,802) were retained. We preserved 2,047 cross-source label conflicts (peptides labeled positive in one source and negative in another) rather than discarding them, because such disagreement is itself an honest signal about curation variability (Figure 2, Table S3).

A previously-dead database, dbPepNeo, was rescued from the inside of a Neodb Zenodo zip, contributing 706 high-confidence MHC-I plus 55 MHC-II peptides.

### 2.2 · Conventional random CV substantially overestimates strict generalization

Under random stratified 5-fold cross-validation on the leakage-free subset, a biophysics-only Logistic-regression model reached AUROC ≈ 0.71, and the combination biophys+HLA-onehot+source-onehot reached AUROC ≈ 0.96 (Figure 3, Table S5). At face value, this would be a strong classifier. But the same combination under Leave-One-Source-Out (LOSO) — training on n-1 sources and testing on the held-out source — fell to AUROC ≈ 0.48–0.58. In other words, the apparently strong random-CV performance was largely driven by source-level patterns that do not transfer to a previously-unseen curator.

The drop is not subtle: source_only (a classifier with only an indicator vector of the source curator) reached AUROC ≈ 0.926 under random 5-fold but exactly 0.500 under LOSO. HLA-only collapsed from 0.745 (random CV) to 0.517 (LOSO). Length-only collapsed from 0.556 to 0.449.

### 2.3 · ESM2 protein-language-model embeddings are the only representation maintaining meaningful LOSO signal — and PLM scaling helps

Across all tested representations on the leakage-free LOSO benchmark (n=2,044, positive_rate ≈ 0.35), only ESM2 protein-language-model embeddings reached AUROC ≥ 0.6 (Figure 4):

| Representation | Random 5-fold AUROC | LOSO AUROC (RF) |
|---|---|---|
| length_only | 0.55 | 0.47 |
| biophys (24-d) | 0.71 | 0.56 |
| HLA one-hot | 0.74 | 0.49 |
| source one-hot | 0.93 ⚠ | 0.50 |
| HLA + biophys | 0.78 | 0.52 |
| source + HLA + biophys | 0.96 ⚠ | 0.48 |
| **ESM2-35M (480-d, 33.5M params)** | 0.76 | **0.74** ★ |
| **ESM2-150M (640-d, 150M params)** | — | **0.75** ★★ |
| **ESM2-150M + biophys (664-d)** | — | 0.69 (LR) / 0.69 (MLP-2L) |
| ESM2 + HLA | 0.78 | 0.74 |

**PLM scaling effect (verified on this run)**: under Random Forest, going from ESM2-35M (0.737 LOSO) to ESM2-150M (0.754 LOSO) gives a small but consistent +0.017 AUROC. The most striking improvement is on **NEPdb-out**: ESM2-35M+RF = 0.813, **ESM2-150M+RF = 0.856** (vs biophys 0.40-0.48; an absolute gain of ~+0.46 over the biophysics baseline). On TESLA-mmc7-out, ESM2-150M+biophys + 2-layer MLP reached **0.900 AUROC**.

**Architecture comparison (representation fixed)**: with ESM2-150M as the representation, Logistic Regression LOSO = 0.59, Random Forest = 0.75, 2-layer MLP (128→64) = 0.67. The architecture spread (~0.16) is comparable to the representation spread (biophys-RF 0.56 vs ESM2-150M-RF 0.75). However, when we fix the representation at biophys-only and vary architecture (LR 0.45, RF 0.56, MLP 0.51), the architecture spread is ~0.11 — never reaching the cross-representation gain. This is consistent with the systematic-benchmark hypothesis that **representation choice is at least as impactful as architecture choice**, and that representation gains accumulate across architectures.

### 2.3.3 · Task-aligned binding features outperform PLM under LOSO

We tested MHCflurry (Class1 affinity predictor; Albert *et al.* 2023 *Nat Mach Intell*) as a task-aligned peptide-HLA binding-feature baseline on the subset with mhcflurry-supported HLA alleles (n=1,954, 95% of the leakage-free benchmark; 3 sources after filtering for ≥50 test rows). Two affinity-derived features (per-row `−log10(predicted nM affinity)` and `−prediction_percentile`) were used:

| representation | model | LOSO AUROC |
|---|---|---|
| biophys (24-d) | RF | 0.556 |
| ESM2-35M (480-d) | RF | 0.737 |
| ESM2-150M (640-d) | RF | 0.754 |
| **MHCflurry affinity (1-d, log-nM)** | **LR** | **0.827 ± 0.091** ★★ |
| MHCflurry affinity 2-d (nM + percentile) | RF | 0.723 |

The single-feature `−log10(nM)` mhcflurry affinity score, combined with class-balanced Logistic regression, reached the highest LOSO AUROC observed in this study (0.827 ± 0.091, n=3 held-out sources). This result has two implications. First, **task-aligned peptide-HLA binding-affinity features are a stronger lever for cross-source generalization than peptide-only protein-language-model embeddings**, even when the PLM is 4× larger (650M vs 150M is forthcoming). Second, the gap between biophys (0.56) and mhcflurry (0.83) in LOSO is +0.27 AUROC — substantially larger than the ESM2-150M-vs-biophys gap (+0.20) — directly supporting the central hypothesis that the most important neoantigen-prediction lever is task-aligned representation, not generic protein-LM scaling alone.

**Caveat**: MHCflurry was trained on IEDB-derived data and possibly overlaps with peptides in our TRAINING_OVERLAP partition; our HELD_OUT (NEPdb) and EXTERNAL_TEST (TESLA) sources were collected after MHCflurry's training cutoff, so the LOSO result is conservative but not perfectly leakage-free. A formal test would require re-training MHCflurry from scratch on a held-out training-data partition — out of scope for this study.

These observations are consistent with the philosophy that **task-aligned peptide-level representations carry more transferable signal than biochemical hand-crafted features**. They do not, however, support a claim that ESM2 alone is sufficient for production-grade neoantigen immunogenicity prediction.

### 2.4 · HLA-composition confounds biochemical feature attribution

We compared pooled biochemical effects (positive vs negative mean) with per-HLA stratified meta-effects on the leakage-free benchmark (Table S6). The pooled `f_H` (Histidine fraction) effect was Δ = +0.010 (p ≈ 5.7e-3), but the per-HLA meta delta was nearly the same (+0.011), so on this subset histidine survived stratification. In contrast, hydrophobicity flipped: pooled Δ = +0.14, per-HLA meta = -0.08. A separate larger-sample analysis (88,753 rows including TRAINING_OVERLAP) showed a stronger pooled histidine signal (Δ ≈ +0.024, p ≈ 4.8e-22) that nearly vanished after per-HLA stratification — a clear Simpson's-paradox confound on that broader sample (Figure 5).

The aromatic-fraction depletion in immunogenic peptides remained directionally consistent (negative) under per-HLA stratification at multiple sample sizes. Length showed a robust direction across pooled and per-HLA analyses (immunogenic peptides shorter on the leakage-free subset).

The general principle is robust: **pooled biochemical attribution on a multi-source, HLA-imbalanced corpus can be misleading**. Per-HLA stratification is a mandatory check before claiming biological insight.

### 2.5 · Shortcut tests confirm source/HLA memorization risk

Under random 5-fold CV, shortcut classifiers reached high AUROC (Figure 6, Table S7):

- source_only: 0.926
- HLA_only: 0.745
- source + HLA: 0.948
- source + HLA + length: 0.962

Under LOSO these collapsed to 0.500, 0.517, 0.500, 0.500 respectively — confirming that random-CV scores are explained by source memorization, not biology.

Label-shuffle controls were also run. Within-source label shuffling (which destroys peptide signal but preserves source-marginal label distribution) gave random-CV AUROC = 0.538 with `aa_composition` features. Within-source × HLA shuffling gave 0.537. Compared with the unshuffled `aa_composition` random-CV AUROC of 0.577, the residual signal is ~0.04 — most of the apparent random-CV signal is source-distribution memorization.

A peptide-zeroed classifier (peptide features replaced by zeros) plus HLA one-hot reached random-CV AUROC = 0.745, the same as HLA-only — confirming that HLA labels alone suffice to explain a large fraction of nominal random-CV performance.

### 2.6 · Candidate prioritization under strict-generalization model

We re-scored 230 prior candidate peptides (180 KRAS G12 9-11mers + 50 active-learning candidates) under a Logistic-regression model trained on the leakage-free benchmark only (Table 2). Top KRAS candidates remained KRAS G12D family (VVVGADGVG: new strict score 0.691; LVVVGADGV: 0.677; VVVGAVGVG: 0.685 for G12V) — but absolute scores were lower than the prior leakage-prone XGBoost score (e.g. GADGVGKSAL was 0.871 under the prior model and is now ~0.65 under strict training). 31 of 230 candidates were flagged for safety review based on Lev≤1 self-similarity; the top KRAS G12D peptides did not trigger this flag.

These re-prioritized scores are computational candidate prioritization only. They do not constitute clinical efficacy or pre-clinical validation and must be confirmed in HLA-tg mouse or *in vitro* T-cell assays.

### 2.7 · Structural and TCR layers as external annotation

We provide six canonical RCSB pMHC and TCR-pMHC reference structures (1HHK HLA-A*02:01 + peptide; 3I6L HLA-A*24:02 + peptide; 1Q94 HLA-A*11:01 + peptide; 3VCL HLA-B*07:02 + 12mer; 1AO7 TCR A6 + HLA-A*02:01 + Tax; 2BNQ TCR LC13 + HLA-B*8:01 + EBV) for visualization context. We also generated 30 ESMFold structures of TESLA-validated peptides and 79 ESMFold structures of cancer-associated CDR3β loops, with explicit pLDDT 0.66–0.72 (medium confidence) and the disclosure that peptide-only ESMFold structures do not capture HLA-groove context. **These are annotation/visualization layers, not production-grade pMHC complex predictions.**

VDJdb cross-reactivity analysis showed that 290 cancer-positive peptides have exact-match in VDJdb's 134,633-record TCR-pMHC database, and that one of our top KRAS G12D candidates (GADGVGKSAL) has a Lev=1 neighbor (GAAGVGKSAL) annotated as HomoSapiens-recognized in VDJdb. This is descriptive annotation support, not independent immunogenic validation.

---

## 3 · Discussion

The central conclusion is that **conventional random k-fold cross-validation substantially overestimates neoantigen immunogenicity prediction performance** in the presence of source and HLA-composition bias. The size of the inflation is large: source-only and HLA-only shortcut classifiers reach AUROC > 0.7 under random CV but collapse to chance under LOSO. Any neoantigen predictor reporting random-CV AUROC without strict-split validation is, by default, claiming benchmark performance that is partially explained by source memorization.

A second conclusion is that **representation choice dominates over architecture choice** for cross-source generalization. ESM2-35M was the only representation tested that maintained AUROC ≥ 0.6 under LOSO, despite being smaller (480-d, 33.5M parameters) than the source+HLA shortcut feature space. Architecture comparisons (Logistic vs Random Forest in this analysis) gave similar conclusions for the same representation.

A third conclusion is that **biochemical attribution requires per-HLA stratification**. The pooled histidine-fraction effect that dominated SHAP explanations in our prior work nearly vanished under per-HLA stratification on the larger sample — a Simpson's-paradox confound that any claim of biological insight must address.

For candidate prioritization, the practical implication is that absolute scores from leakage-prone models (XGBoost trained on the full corpus including TRAINING_OVERLAP) systematically over-estimate immunogenicity probability. The correct way to use the corpus is to train on the leakage-free subset, validate under strict splits, and treat the resulting scores as conservative prioritization signals. The KRAS G12D family remains the strongest candidate set on this corpus, but absolute scores need recalibration.

We did not run NetMHCpan, MHCflurry, or PRIME under our strict splits, primarily because of binary licensing constraints on this CPU-only environment. The expectation, based on the systematic-benchmark philosophy, is that task-aligned peptide-HLA binding-and-presentation features should outperform pure peptide-level PLM features. Confirming this is the most important next experiment.

### Limitations

1. No wet-lab validation. All candidate scores are computational prioritization only.
2. No clinical-efficacy claim. We provide a benchmark resource and a preclinical roadmap.
3. Source bias remains: HLA-A*02:01 dominates our corpus (~3.5× over Korean A*24:02 vs population frequency).
4. HLA-II curation is positive-class-biased (~95% positive); class II analyses are exploratory.
5. TCR / VDJdb support is annotation, not proof of immunogenicity.
6. ESMFold isolated peptide / CDR3β structures are not pMHC complex predictions.
7. NEPdb's negative-heavy distribution may inflate LOSO failure rate; quantifying curation-artifact contribution requires a held-out wet-lab cohort.
8. Population coverage uses a simplified independence assumption; production deployments should use the IEDB population coverage tool.
9. We did not run a NetMHCpan/MHCflurry/PRIME baseline under strict splits.
10. The architecture-ablation analysis is limited to LR / RF; ranking and adversarial losses were not implemented in this round.

---

## 4 · Methods (abridged — full Methods in Supplementary)

Datasets: 13 sources (TESLA, IEDB, NEPdb, McPAS-TCR, IMPROVE, Neodb, dbPepNeo, TSNAdb, CEDAR, BigMHC, NeoRanking, VDJdb, clinical). Cleaning: peptide canonical-AA validation (length 8–25), HLA normalization (`HLA-A*XX:XX`), patient-aware dedup. Leakage audit via per-source `test_set_safety` rules (Supplementary Table S2). Features: length, AA composition, biophysics (Kyte-Doolittle hydrophobicity, charge, aromatic, AA frequencies; 24-d), HLA one-hot, source one-hot, ESM2-35M mean-pool (480-d). Models: Logistic (L2), Random Forest (200 trees, depth 8), XGBoost (300 trees, depth 6). Splits: random stratified 5-fold; peptide GroupKFold 5; peptide-HLA GroupKFold 5; Leave-One-Source-Out; Leave-One-HLA-Out. Metrics: AUROC, AUPRC, balanced accuracy, MCC, Brier score; bootstrap 95% CI (200 samples); Wilcoxon paired tests across folds. Statistical primary comparisons: best-task-aligned-representation vs biophys under LOSO; ESM2 vs source+HLA shortcut under LOSO. Candidate re-scoring on KRAS G12 9-11mer enumeration with self-similarity safety check (Lev≤1 to a 1,212-peptide WT pool). Code: `/data/neoantigen_vaccine_hub/scripts/`.

---

## 5 · Data and Code availability

Public hub: http://40.82.129.113/papers_hub_2026_05_04/
Master corpus: `data_processed/neoantigen_master.csv` (44 MB, 119,237 records)
Leakage-free benchmark: `experiments/what_matters_neoantigen/cache/benchmark_clean.tsv` (84,549 labeled records)
Frozen snapshot manifest: `manuscript/data_snapshots/MANIFEST.tsv`
Reproducibility audit: `manuscript/reproducibility_audit.md`
All scripts: `/data/neoantigen_vaccine_hub/scripts/` (see `manuscript/SUBMISSION_READINESS_REPORT.md`)

---

## 6 · References

(see `manuscript/citations.bib`; key sources: Wells 2020 Cell, Vita 2019 NAR, Xia 2021 NAR, Tickotsky 2017 Bioinf, Borch 2024 Front Immunol, Wu 2025 Zenodo, Wu 2020 Brief Bioinf, Wu 2018 GPB, Koşaloğlu-Yalçın 2023 NAR, Albert 2023 Nat Mach Intell, Müller 2023 Cell Rep Med, Bagaev 2020 NAR, Lin 2023 Science, Balachandran 2017 Nature)
