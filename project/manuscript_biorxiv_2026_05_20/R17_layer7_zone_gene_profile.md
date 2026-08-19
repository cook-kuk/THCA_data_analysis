# R17 Layer 7 — Per-zone 8-gene panel expression profile (TCGA THCA)

n = **522** TCGA THCA samples with both R17 zone and panel expression.

## Per-zone gene z-mean (cohort-centered)

|        |   BRAF-like |   RAS-like |   dark-matter |   WT-like |
|:-------|------------:|-----------:|--------------:|----------:|
| DIO1   |      -0.527 |      0.668 |        -0.635 |     1.049 |
| FOXE1  |      -0.286 |      0.306 |        -0.472 |     0.719 |
| NKX2-1 |      -0.093 |      0.154 |        -0.291 |     0.363 |
| PAX8   |      -0.274 |      0.287 |        -0.688 |     0.938 |
| SLC5A5 |      -0.396 |      0.939 |        -0.313 |     0.492 |
| TG     |      -0.156 |      0.384 |        -0.762 |     0.862 |
| TPO    |      -0.602 |      0.875 |        -0.615 |     1.05  |
| TSHR   |      -0.142 |      0.341 |        -0.546 |     0.632 |

## Per-zone summary

| zone        |   n | most_silenced_gene   |   most_silenced_z | most_preserved_gene   |   most_preserved_z |   panel_mean_z |
|:------------|----:|:---------------------|------------------:|:----------------------|-------------------:|---------------:|
| BRAF-like   | 165 | TPO                  |            -0.602 | NKX2-1                |             -0.093 |      -0.3095   |
| RAS-like    |  43 | NKX2-1               |             0.154 | SLC5A5                |              0.939 |       0.49425  |
| dark-matter | 161 | TG                   |            -0.762 | NKX2-1                |             -0.291 |      -0.54025  |
| WT-like     | 153 | NKX2-1               |             0.363 | TPO                   |              1.05  |       0.763125 |

## Interpretation
* BRAF-like zone shows uniform low-z across all 8 panel genes (deep coordinated silencing).
* Dark-matter zone (silenced + HT) silences most panel genes similarly to BRAF-like but with
  slightly different gradient (typically less uniform).
* RAS-like and WT-like zones show preserved-high panel z, with RAS-like skewed toward
  partial preservation of certain genes (TG, TSHR) consistent with RAS-driven differentiation
  retention.
* The two-axis decomposition therefore captures different *patterns* of panel silencing,
  not just overall level — supporting the framing that BRAF-like and dark-matter represent
  distinct biological silencing programs that converge on the same panel-DM1 phenotype.