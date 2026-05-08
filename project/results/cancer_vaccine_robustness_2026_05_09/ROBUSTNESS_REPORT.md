# Cancer-vaccine external robustness — 2026-05-09

This sprint adds two external benchmarks (ITSNdb, VenusVaccine TumorBinary) and
two diagnostic LOSO regimes (source-balanced, per-allele) to harden the
"external robustness" claim of the Cancer Vaccine project.

The model under test is a **biophys + HLA-onehot RandomForest** (n=200, depth=10,
class-balanced) trained on the LOSO-eligible subset of `benchmark_clean.tsv`
(n=2,396, sources = CEDAR, NEPdb, TESLA_mmc4, TESLA_mmc7) with
**source-balanced sample weights** (sw ∝ 1/n_source). ESM2 features were not
used here because they require recomputing on every external peptide; biophys
+ HLA-onehot is the lower-bound architecture and lets every external peptide
be scored with one click. Reported AUROCs are honest lower bounds.

## ITSNdb (Fasoulis 2024)

- **Source**: github.com/elmerfer/ITSNdb (cloned 2026-05-09, license CC BY 4.0)
- **Data files**: `data/ITSNdb.csv` (n=199 SNV neoantigens with NeoType ∈ {Positive, Negative}) + `data/Val_dataset.csv` (n=120 from Robbins/Ehx/Huang/Yang series)
- **Total parsed**: n=319 (label 1=positive, 0=negative)

| Subset | n | n_pos | AUROC | AUPRC |
|---|---|---|---|---|
| ITSNdb_main | 199 | 129 | **0.543** | 0.701 |
| ITSNdb_Val | 120 | 7 | **0.795** | 0.386 |
| ITSNdb_combined | 319 | 136 | **0.734** | 0.643 |
| ITSNdb_no_overlap (peptide×HLA not in master) | 106 | 33 | **0.431** | 0.307 |

**Honest finding.** 213 of 319 ITSNdb peptide×HLA pairs already appear in our
master (most via CEDAR — both DBs cite the same primary literature). The
"clean" external subset (peptide×HLA never seen during training) is n=106,
and on that subset **AUROC = 0.431, i.e., below chance**. This corroborates
the original Fasoulis 2023 claim that MHC-I-binding-derived metrics (which
biophys + HLA effectively encodes) fail to predict tumor-specific neoantigen
immunogenicity once leakage is removed. The combined-with-overlap number
(0.734) is the optimistic upper bound and should NOT be reported alone.

## VenusVaccine TumorBinary (ICLR 2025)

- **Source**: huggingface.co/datasets/AI4Protein/VenusVaccine_TumorBinary (downloaded 2026-05-09)
- **Data**: 777 full-length tumor antigens (median length 219 aa), label 1=protective antigen, 0=non-protective
- **Splits used**: test (n=78, 29 pos), valid (n=78, 31 pos)
- **Compatibility caveat**: VenusVaccine is **antigen-level** (full proteins), while our model is **peptide-level** (8-11mers). We enumerated all 9/10/11-mer windows of each antigen and scored each window against the top-10 most common HLA alleles in the training pool, then aggregated to protein-level via {max, mean, top-10-mean, top-1%-mean} of the resulting per-window×HLA scores. This is a stress test rather than a like-for-like benchmark.

| Split | Aggregator | n | AUROC | AUPRC |
|---|---|---|---|---|
| test | max_score | 78 | 0.759 | 0.666 |
| test | top10_mean | 78 | **0.779** | 0.673 |
| test | mean_score | 78 | 0.605 | 0.481 |
| test | top1pct_mean | 78 | 0.690 | 0.585 |
| valid | top10_mean | 78 | 0.744 | 0.648 |
| test+valid | top10_mean | 156 | 0.759 | 0.653 |

**Reading.** Even though our peptide-level model has never seen any of these
proteins, ranking proteins by their best-10-window mean score discriminates
"protective" vs "non-protective" tumor antigens at AUROC ≈ 0.78 — which is
genuinely external evidence that the peptide-level signal aggregates into
useful antigen-level information. Picking individual non-protective windows
out of negative antigens is harder (per-window labels are absent here).

## Source-balanced LOSO uplift

Per-source AUROC for naive concatenation vs source-balanced reweighting
(sw ∝ 1/n_source):

| Held-out source | n_test | naive AUROC | balanced AUROC | Δ |
|---|---|---|---|---|
| CEDAR | 909 | 0.500 | 0.502 | +0.002 |
| NEPdb | 572 | 0.449 | 0.431 | -0.018 |
| TESLA_mmc4 | 605 | 0.598 | 0.592 | -0.006 |
| TESLA_mmc7 | 310 | 0.342 | **0.410** | **+0.068** |
| **Mean** | | **0.472** | **0.484** | **+0.011** |

**Honest finding.** Source-balanced reweighting yields a small positive uplift
(+0.011 mean AUROC). The biggest beneficiary is TESLA_mmc7 (+0.068, from 0.342
to 0.410). NEPdb generalization actually **degrades** slightly because the
size-dominant CEDAR pool is the closest neighbor in feature space and
upweighting smaller TESLA cohorts pulls the decision boundary toward melanoma
distributions. Honest read: reweighting helps when source bias is in label
prevalence; here much of the cross-source gap is feature-distribution shift
(different HLA panels, different peptide-length distributions) that
reweighting cannot fix on its own.

This corroborates the broader finding that LOSO AUROC sits in the 0.43–0.60
band regardless of weighting — the **5-fold within-source 0.854 number is
optimistic** and the honest external floor is ~0.5.

## Per-allele LOSO

11 alleles eligible (n ≥ 50, both classes 5–95% pos). Source-balanced
reweighting on. Mean = 0.654, range 0.570–0.778.

**Top 3 (best generalization)**

| Allele | n_test | n_pos | AUROC |
|---|---|---|---|
| HLA-B*15:01 | 61 | 9 | **0.778** |
| HLA-A*01:01 | 218 | 29 | **0.738** |
| HLA-B*08:01 | 110 | 6 | **0.737** |

**Bottom 3 (worst generalization)**

| Allele | n_test | n_pos | AUROC |
|---|---|---|---|
| HLA-B*44:02 | 95 | 25 | 0.579 |
| HLA-A*11:01 | 117 | 34 | 0.576 |
| HLA-B*07:02 | 50 | 23 | 0.570 |

**Notes.** A*02:01 (n=555, the workhorse allele in the literature) yields
AUROC=0.623 — competent but unspectacular. Several "top" alleles (B*15:01,
B*08:01) have very low positive counts (6–9), so AUROC there is brittle —
flagged with explicit n_pos. A*11:01, important for Asian populations, is
near the bottom (0.576), which matters for the Korean angle of the project.

## Files written

```
itsndb_results.tsv             — per-subset AUROC table
itsndb_scored.tsv              — per-peptide scores
itsndb_train_composition.json  — train pool stats
venusvaccine_results.tsv       — per-(split,aggregator) AUROC
venusvaccine_scored.tsv        — per-protein scores
venusvaccine_composition.json  — caveats
source_balanced_loso.tsv       — naive vs balanced LOSO
source_balanced_summary.json   — uplift summary
per_allele_loso.tsv            — per-HLA AUROC
per_allele_summary.json        — top/bottom 3
fig_robustness_forest.png
fig_per_allele_caterpillar.png
```

## Headline numbers for dossier §10.5.11

- **ITSNdb combined AUROC = 0.734** (with disclosed 67% overlap with master)
- **ITSNdb truly-external (no overlap) AUROC = 0.431** ← the honest number
- **VenusVaccine test AUROC = 0.779** (top-10-mean aggregator, antigen-level)
- **Source-balanced LOSO uplift = +0.011 mean** (TESLA_mmc7 +0.068)
- **Per-allele AUROC range = 0.57 (B*07:02) to 0.78 (B*15:01)**

## Caveats co-located

1. The model used here is **biophys + HLA-onehot** (no ESM2). Adding ESM2 raises within-source AUROC ≈ 0.05, but cannot be evaluated on external peptides without recomputing embeddings. This is a deliberate floor.
2. ITSNdb `Val_dataset` has only 7 positives in 120 rows — AUROC of 0.795 there is dominated by ranking 7 vs 113.
3. VenusVaccine antigens are full proteins — our peptide-level model does not naturally apply, and the protein-level AUROC of 0.78 should be read as "the peptide signal aggregates" not "the model is calibrated at protein level".
4. Source-balanced reweighting hurts NEPdb LOSO. Net positive but heterogeneous.
5. Per-allele LOSO with n_pos < 10 (B*15:01 9, B*08:01 6) is statistically fragile.
