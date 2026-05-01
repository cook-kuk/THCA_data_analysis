# Lu et al. 2023 (GSE193581) — v17p35 8-gene panel transfer (SLIM)

**Generated:** 2026-04-28
**Cohort:** GSE193581 — Lu et al. *J Clin Invest* 2023;133:e169653
**Design:** 10x Genomics 3' scRNA-seq, 23 samples (NORM + PTC + ATC)
**Strategy:** memory-slim — extract only 8 panel rows post-normalization, no full HVG/scale/PCA/UMAP

---

## Summary

| metric                              | value                                    |
|-------------------------------------|------------------------------------------|
| Samples used                        | 23 ({'ATC': 9, 'NORM': 6, 'PTC': 6}) |
| Total cells post-QC                 | 67,678                        |
| Malignant cells used                | 15,330                    |
| Malignant by histology              | NORM=675, PTC=8,621, ATC=6,034 |
| Panel genes found                   | 8/8                |
| Panel genes missing                 | none |

## Cross-cohort statistics

* Kruskal-Wallis 3-group (['NORM', 'PTC', 'ATC']) DM_score: H=**11175.3**, p=**0.00e+00**
* Mann-Whitney PTC vs ATC: U=51270308, p=**0.00e+00**
* chi-square histology × DM_class: chi²=9724.8, p=**0.00e+00**
* ARI(histology, DM_class) = **0.550**

NORM median DM_score = +1.355; PTC median = +0.291; ATC median = -0.594.

## Per-sample
```
       histology  n_cells  median_DM  mean_DM  pct_DM1
sample                                                
ATC08        ATC      213     -0.594   -0.584    0.005
ATC10        ATC       67     -0.594   -0.590    0.000
ATC11        ATC      174     -0.594   -0.592    0.000
ATC12        ATC     1296     -0.594   -0.593    0.000
ATC13        ATC     3308     -0.594   -0.565    0.008
ATC14        ATC       15     -0.594   -0.357    0.133
ATC15        ATC       56     -0.594   -0.200    0.196
ATC17        ATC        8     -0.594   -0.594    0.000
ATC09        ATC      897     -0.524   -0.514    0.000
NORM19      NORM      261      0.979    1.031    0.954
NORM03      NORM      136      1.314    1.188    0.824
NORM07      NORM      123      1.479    1.426    0.919
NORM21      NORM        8      1.747    1.637    1.000
NORM20      NORM       24      2.069    2.406    1.000
NORM18      NORM      123      2.136    2.090    1.000
PTC04        PTC     1197      0.113    0.181    0.629
PTC03        PTC       86      0.194    0.138    0.605
PTC05        PTC      629      0.197    0.169    0.715
PTC01        PTC     3565      0.279    0.269    0.825
PTC07        PTC     2727      0.362    0.372    0.897
PTC06        PTC      417      0.378    0.350    0.849
```

## Files

* `/opt/thyroid-dash/project/results/v17_lu2023/GSE193581_*.tsv` + `.json`
* `/opt/thyroid-dash/project/reports/html/figs_interactive/v17/v17_lu2023_*.html` (4 figures)
* `/opt/thyroid-dash/project/logs/v17_lu2023_slim.log`

*author: Seungho Cook · seed=42*
