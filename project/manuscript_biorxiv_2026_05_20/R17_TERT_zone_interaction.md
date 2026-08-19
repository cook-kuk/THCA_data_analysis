# R17 × TERT — TCGA THCA two-axis zone × TERT promoter status interaction

Samples with known TERT status: **n = 477** of 522 R17-classified TCGA samples.
Overall TERT+ rate: **7.5%**.

## TERT+ fraction by R17 zone

| zone        |   n |   n_TERT_pos |   pct_TERT_pos |
|:------------|----:|-------------:|---------------:|
| BRAF-like   | 155 |            9 |            5.8 |
| dark-matter | 150 |           18 |           12   |
| WT-like     | 134 |            9 |            6.7 |
| RAS-like    |  38 |            0 |            0   |

## Fisher exact (zone vs rest)

| zone        |   TERT+_in_zone |   TERT-_in_zone |   TERT+_rest |   TERT-_rest |    OR |      p |
|:------------|----------------:|----------------:|-------------:|-------------:|------:|-------:|
| BRAF-like   |               9 |             146 |           27 |          295 | 0.674 | 0.36   |
| dark-matter |              18 |             132 |           18 |          309 | 2.341 | 0.0156 |
| WT-like     |               9 |             125 |           27 |          316 | 0.843 | 0.847  |
| RAS-like    |               0 |              38 |           36 |          403 | 0     | 0.1    |

## PFI event rate by zone × TERT status

| zone        | tert_bool   |   n |   events |   event_rate_pct |
|:------------|:------------|----:|---------:|-----------------:|
| BRAF-like   | False       | 146 |       12 |              8.2 |
| BRAF-like   | True        |   9 |        2 |             22.2 |
| RAS-like    | False       |  38 |        2 |              5.3 |
| WT-like     | False       | 125 |        7 |              5.6 |
| WT-like     | True        |   9 |        3 |             33.3 |
| dark-matter | False       | 132 |       18 |             13.6 |
| dark-matter | True        |  18 |        5 |             27.8 |

## Interpretation
* TERT+ is an established late event in dedifferentiation. The R17 zones predict the histologic
  / molecular state along two orthogonal axes (RAI silencing × HT overlap), so the question is
  whether TERT+ enriches in any single zone.
* The dark-matter zone (silenced + HT) is the predicted high-aggression cell-state intersection.
  TERT+ enrichment here would tighten Paper 1's HR=7.57 TERT survival claim into a specific
  cellular context.
* PFI event rate by zone × TERT row shows the joint effect; n is small in some cells so treat
  as exploratory.

## Caveats
* TCGA THCA is overwhelmingly low-risk PTC, so TERT+ rate is low overall (~10-13%); large-effect
  zone enrichment may be hard to detect with this prior.
* PFI event count is small per zone × TERT cell; survival inference within zones is underpowered.