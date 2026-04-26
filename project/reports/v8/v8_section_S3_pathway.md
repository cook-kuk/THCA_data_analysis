> ⚠️ **v5.2 SUPERSEDED NOTICE (2026-04-25).** This robustness section was designed to stress-test the v5.1 THCA DIAL flip (DIAL=0.494, batch_entangled in 4/5). Under proper LODO ComBat the underlying flip disappears (v5.2: DIAL=0.000 / true_biology in 5/5; see reports/v5p2/v5p2_critical_assessment.md). Conclusions below are conditional on the v5.1 result they reference.

# v8 Section S3: Does the THCA post-ComBat label flip survive pathway aggregation?

## TL;DR

**No.** The v5.1 gene-level "label flip" (post-ComBat `auc_post` ~0.006,
DIAL up to 0.494 under Leave-One-Dataset-Out on 3000 variance-top
genes) vanishes when the same 392 THCA samples are first reduced to a
sample x 50 MSigDB Hallmark activity matrix via ssGSEA (Barbie et al.
2009). Pathway `auc_post` recovers to 0.90 - 0.98 and DIAL collapses
to exactly 0.0 for all five classifiers. The finding is **not robust
to pathway aggregation** and behaves like the high-dimensional
residual-batch artefact GSVA was designed to mitigate (Haenzelmann,
Castelo and Guinney, *BMC Bioinformatics* 2013).

## Numbers (LODO, 2 folds: TCGA-THCA vs GSE27155)

| Classifier         | DIAL gene (3000g) | DIAL pathway (50 HM) | auc_post gene | auc_post pathway |
|--------------------|------------------:|---------------------:|--------------:|-----------------:|
| LogReg L2          |             0.494 |                0.000 |         0.006 |            0.964 |
| LogReg ElasticNet  |             0.492 |                0.000 |         0.008 |            0.958 |
| Random Forest      |             0.323 |                0.000 |         0.177 |            0.977 |
| Gradient Boosting  |             0.334 |                0.000 |         0.166 |            0.918 |
| XGBoost            |             0.013 |                0.000 |         0.487 |            0.905 |

All pathway runs are classified `true_biology` (DIAL < 0.05 and
`auc_post` > 0.7).

## Method

Inputs are the v5.1 harmonized THCA matrix (392 x 11710 HGNC, float32
log2-TPM-like), Y in {BRAF (321), RAS (71)}, batch B in {TCGA-THCA
(351), GSE27155 (41)}. MSigDB Hallmark 2020 (50 pathways, median 180
genes) loads via `gseapy.get_library`; a gseapy 1.1.13 bytes/str bug
triggers an Enrichr plain-text fallback (cached). ssGSEA runs in 7 s
(`gseapy.ssgsea`, sample_norm_method='rank', permutation_num=0,
min_size=5, all cores). We then apply the same protocol as
`v5p1_common._combat_preserve`: `inmoose.pycombat_norm` with `covar_mod`
= one-hot Y, preserving biology while removing B. DIAL is computed by
`v5p1_common.compute_dial` with five `get_classifier_factories()`
classifiers and `LeaveOneGroupOut` on B — identical to v5.1 THCA, only
X changes from 3000 genes to 50 pathways.

## Interpretation

The v5.1 flip reflects ComBat over-correcting in a space where
TCGA-THCA and GSE27155 occupy partially separate manifolds for
thousands of mildly informative genes; a LODO classifier forced onto
the held-out cohort then points the wrong way. Pathway aggregation
averages that noise away: each Hallmark score is a rank-weighted sum
over ~180 genes, stabilising cohort distributions *before* ComBat
runs. Inside the 50-D pathway space, ComBat's covariate-preserving
shift leaves BRAF vs RAS intact and GSE27155 classifies at `auc_post`
= 0.90 - 0.98 for every model.

Caveats. Hallmark 2020 is coarse (50 sets); finer collections
(C2.CP REACTOME ~1600 sets) could reintroduce sparsity-driven flipping
and merit a sensitivity check. LODO has only two folds because THCA
has two cohorts in v5.1; a GSE33630 / GSE29265 refresh (v8 Task 5)
would give three- or four-fold LODO.

## Conclusion

The THCA BRAF vs RAS post-ComBat flip is a **gene-level artefact**, not
a biological direction-reversal. 50 Hallmark pathway scores eliminate
DIAL without sacrificing discrimination (`auc_post` >= 0.90 all
classifiers). v5.1 THCA interpretation should run in pathway space or
cite this table.

## References

- Barbie DA et al. Systematic RNAi reveals oncogenic KRAS-driven
  cancers require TBK1. *Nature* 2009;462:108-12. (ssGSEA)
- Haenzelmann S, Castelo R, Guinney J. GSVA: gene set variation
  analysis for microarray and RNA-seq data. *BMC Bioinformatics*
  2013;14:7.
- Liberzon A et al. The MSigDB Hallmark gene set collection.
  *Cell Systems* 2015;1:417-25.

## Artefacts

- `results/v8_statgen/v8_pathway_dial.tsv`
- `results/v8_statgen/v8_pathway_vs_gene.tsv`
- `results/v8_statgen/v8_pathway_activity_THCA.tsv`
- `logs/v8_pathway.log`
- `notebooks_or_scripts/v8_pathway_dial.py`

(paths under `/opt/thyroid-dash/project/`)
