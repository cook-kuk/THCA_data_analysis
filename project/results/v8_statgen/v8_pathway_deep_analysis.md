# Deep analysis: why does the 50-Hallmark pooled DIAL collapse to 0.000?

> ⚠️ **FALLS under v5.2 retraction.** This analysis was framed as
> *"does pathway aggregation eliminate the v5.1 gene-level flip?"*
> Under v5.2 (proper LODO ComBat fit on train only), the gene-level
> THCA flip itself is 0 (auc_post 0.99, DIAL 0.000) — there is no
> flip to eliminate. The pathway-level DIAL = 0 finding still
> reproduces under v5.2 (because both gene and pathway DIAL are 0
> under leak-safe ComBat), but the *contrast* this analysis was built
> on no longer exists. The biological observation (29/50 Hallmark
> pathways are coherent BRAF-vs-RAS biology) survives as a stand-alone
> result; the "gene-level artefact, not biology reversal" framing is
> retracted because the gene-level "flip" was a leak artefact, not a
> high-dimensional residual-batch artefact as originally claimed.

_v8 supplement to Task 3. Input: `v8_pathway_activity_THCA.tsv` (392 samples ×
50 MSigDB Hallmark pathways, ssGSEA). Output:
`v8_pathway_perpath_dial.tsv`. Generated 2026-04-24._

## The question

Task 3 ran the full v5.1 LODO DIAL pipeline on the 50-pathway matrix and got
DIAL = **0.000 for every one of the five classifiers**, while the 3000-gene
pipeline had DIAL up to 0.494. Two non-exclusive explanations:

1. **Aggregation-dilutes.** Some individual pathways still flip, but the
   pooled classifier averages flipped and non-flipped pathways and the net
   DIAL cancels.
2. **Coherent biology.** No individual pathway flips; the 3000-gene flip was
   a high-dimensional residual-batch artefact that pathway aggregation removes
   by stabilising each sample's cohort-mean in a 180-gene-neighbourhood
   average.

Discriminator: run LODO DIAL **once per pathway** (each column as a single
feature, same ComBat-preserve step applied to the joint 50-D matrix) and
inspect the per-pathway DIAL distribution.

## Method

- Feature matrix: 50-pathway ssGSEA activity from Task 3 (392 × 50).
- ComBat: `_combat_preserve(X, Y, B)` on the full 50-D matrix, identical to
  v8 Task 3. (ComBat on a 1-D column is numerically unstable — `inmoose`
  returns NaN — so we decompose the **same** 50-D ComBat output column-wise.)
- Per pathway j: LogReg_l2 (C=1.0) on the single-column feature under LODO on
  B (TCGA-THCA ↔ GSE27155). AUC_pre from raw column, AUC_post from ComBat
  column, DIAL_j = max(0, AUC_pre − max(AUC_post, 1 − AUC_post)).
- Rules match v5.1: DIAL ≥ 0.3 → `batch_entangled`; AUC_pre ≥ 0.7 and
  AUC_post ≥ 0.7 → `true_biology`; AUC_pre < 0.6 → `no_signal`; else
  `ambiguous`.

## Result

| Metric                         | Value  |
|--------------------------------|-------:|
| n pathways tested              |     50 |
| DIAL median / mean / max       | 0.000 / 0.005 / 0.085 |
| `batch_entangled` (DIAL ≥ 0.3) |  **0** |
| `true_biology`                 |     29 |
| `no_signal`                    |      8 |
| `ambiguous`                    |     13 |
| flipped_direction (post < 0.5) |  4 / 50 |

**No pathway individually crosses the batch_entangled threshold.** The
4 that technically invert post-ComBat all have AUC_pre < 0.5 (null
pre-correction) — these are noise-floor inversions, not biology
reversals: Cholesterol Homeostasis, Hedgehog Signaling, Unfolded
Protein Response, mTORC1 Signaling (each with DIAL = 0.000 because
AUC_pre_flipped ≈ AUC_post_flipped).

### Top-10 pathways by DIAL (all well below the 0.3 threshold)

| Pathway                         | AUC_pre | AUC_post | DIAL   | Interpretation |
|---------------------------------|--------:|---------:|-------:|----------------|
| Oxidative Phosphorylation       | 0.897   | 0.813    | 0.085  | true_biology   |
| DNA Repair                      | 0.639   | 0.616    | 0.023  | ambiguous      |
| TGF-beta Signaling              | 0.807   | 0.788    | 0.020  | true_biology   |
| Adipogenesis                    | 0.787   | 0.770    | 0.017  | true_biology   |
| p53 Pathway                     | 0.845   | 0.830    | 0.015  | true_biology   |
| Angiogenesis                    | 0.854   | 0.840    | 0.015  | true_biology   |
| Androgen Response               | 0.720   | 0.705    | 0.014  | true_biology   |
| Coagulation                     | 0.954   | 0.944    | 0.010  | true_biology   |
| Fatty Acid Metabolism           | 0.832   | 0.823    | 0.008  | true_biology   |
| Hypoxia                         | 0.787   | 0.780    | 0.007  | true_biology   |

The highest-DIAL pathway is Oxidative Phosphorylation (0.085) — a well-known
RAS-high / BRAF-low signature in THCA (Landa *et al.* 2016 *Cell*). AUC_pre
0.897 → AUC_post 0.813 is a mild ComBat over-correction on a real biology
axis, not a flip. DIAL 0.085 is well below the 0.3 batch_entangled bar and
below even the 0.1 "tolerance" line shown in v8 S1 figures.

## Interpretation

**Explanation 2 (coherent biology) wins.** The pooled DIAL=0.000 is
*not* an artefact of averaging flipped and non-flipped pathways. There
are no flipped pathways to average: the maximum per-pathway DIAL is
0.085, 29/50 pathways are clean `true_biology`, and the remaining
17/50 are null or ambiguous — none is `batch_entangled`.

The 3000-gene flip is therefore a **high-dimensional residual-batch
phenomenon**, not a reversal of BRAF-vs-RAS pathway biology. Each
Hallmark score is a rank-weighted mean over ≈ 180 member genes; this
average stabilises the GSE27155 cohort distribution relative to
TCGA-THCA *before* ComBat runs, so the covariate-preserving shift
leaves BRAF vs RAS intact.

### Implication for v5.1 interpretation

The pathway result is **constructive**, not "just dilution":

- The BRAF-vs-RAS signal is real and biologically coherent (29/50
  Hallmark pathways achieve AUC ≥ 0.7 pre- and post-ComBat).
- The v5.1 gene-level DIAL 0.494 is a noise-level artefact of ComBat
  in the 3000-gene space, not a biology reversal.
- **Recommended v5.1 downstream protocol:** report pathway-level AUCs
  alongside gene-level DIAL. For the published paper this means the
  "batch_entangled" THCA classifier row in Table 1 should carry a
  pathway-level addendum — e.g. "gene-level DIAL = 0.494, pathway-level
  DIAL = 0.000, pathway-level AUC = 0.96 (LogReg_l2)".

## Caveats

- Hallmark 2020 is coarse (50 sets, median 180 genes). Running this
  analysis with C2.CP.REACTOME (≈ 1 600 sets, median 23 genes) would
  partially bridge back toward the 3000-gene regime; sparse small
  gene-sets may reintroduce per-set flipping. Planned sensitivity
  check (not yet run).
- LODO has only 2 folds because THCA has only 2 cohorts in v5.1; a
  GSE33630 / GSE29265 refresh would give 3- or 4-fold.
- LogReg_l2 is the only single-pathway-feature classifier tested.
  Tree ensembles on 1-D features collapse to a univariate threshold
  and would reproduce the LogReg result up to tie-handling; the
  choice is deliberate.

## References

- Haenzelmann S *et al.* GSVA. *BMC Bioinformatics* 14:7 (2013).
- Barbie DA *et al.* ssGSEA. *Nature* 462:108 (2009).
- Liberzon A *et al.* Hallmark gene set collection.
  *Cell Systems* 1:417 (2015).
- Landa I *et al.* Integrative genomic characterization of
  poorly-differentiated and anaplastic thyroid cancer.
  *Cell* 169:803 (2016). (OxPhos/RAS-high signature)
- v8 Task 3 report: `reports/v8/v8_section_S3_pathway.md`.
