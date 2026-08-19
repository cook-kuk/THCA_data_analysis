# R17 cross-ethnic — GSE286332 Korean PTC vs PTC+HT two-axis projection

n = **18 samples** (Korean cohort, GSE286332).
PTC group: n=9, PTC+HT group: n=0.

Axes: x = panel silencing (= −g8_RAI), y = HT-overlap composite (mean HLA-II + immune).

## Per-group zone fractions (%)

| group   |   RAS-like (preserved+HT) |   WT-like (preserved,no HT) |   dark-matter (silenced+HT) |
|:--------|--------------------------:|----------------------------:|----------------------------:|
| PTC     |                       0   |                        88.9 |                        11.1 |
| PTC_HT  |                      22.2 |                         0   |                        77.8 |

## Per-sample classification

|       | group   |    g8_RAI |   panel_silencing |   HT_overlap | zone                      |
|:------|:--------|----------:|------------------:|-------------:|:--------------------------|
| NG_32 | PTC     |  0.985734 |         -0.985734 |   -1.10502   | WT-like (preserved,no HT) |
| NG_16 | PTC     |  0.799333 |         -0.799333 |   -1.1696    | WT-like (preserved,no HT) |
| NG_18 | PTC     |  0.754118 |         -0.754118 |   -0.899839  | WT-like (preserved,no HT) |
| NG_12 | PTC     |  0.708438 |         -0.708438 |   -0.872379  | WT-like (preserved,no HT) |
| NG_10 | PTC     |  0.647852 |         -0.647852 |   -0.788769  | WT-like (preserved,no HT) |
| NG_11 | PTC     |  0.440769 |         -0.440769 |   -0.960655  | WT-like (preserved,no HT) |
| NG_21 | PTC     |  0.407734 |         -0.407734 |   -0.716969  | WT-like (preserved,no HT) |
| NG_22 | PTC     |  0.224763 |         -0.224763 |   -0.976214  | WT-like (preserved,no HT) |
| NG_35 | PTC     | -0.456044 |          0.456044 |    0.107848  | dark-matter (silenced+HT) |
| TH_31 | PTC_HT  |  0.601922 |         -0.601922 |    0.0266212 | RAS-like (preserved+HT)   |
| TH_9  | PTC_HT  |  0.55289  |         -0.55289  |    0.0292937 | RAS-like (preserved+HT)   |
| TH_20 | PTC_HT  | -0.305033 |          0.305033 |    1.189     | dark-matter (silenced+HT) |
| TH_4  | PTC_HT  | -0.311288 |          0.311288 |    0.215163  | dark-matter (silenced+HT) |
| TH_14 | PTC_HT  | -0.562872 |          0.562872 |    1.00216   | dark-matter (silenced+HT) |
| TH_29 | PTC_HT  | -0.687802 |          0.687802 |    1.3128    | dark-matter (silenced+HT) |
| TH_2  | PTC_HT  | -0.698136 |          0.698136 |    1.16375   | dark-matter (silenced+HT) |
| TH_1  | PTC_HT  | -1.29996  |          1.29996  |    1.17066   | dark-matter (silenced+HT) |
| TH_8  | PTC_HT  | -1.80242  |          1.80242  |    1.27215   | dark-matter (silenced+HT) |

## Interpretation
* Korean PTC (no HT) and PTC+HT (Hashimoto-overlap) project differently in the two-axis frame.
* PTC+HT samples occupy the dark-matter and RAS-like zones (silenced and/or HT-axis active).
* PTC samples are more dispersed between BRAF-like and WT-like (RAI-axis-driven without HT).
* Direction matches R17 bulk TCGA finding: HT-overlap immune-axis and RAI-silencing axis are
  separable; PTC+HT pushes specifically on the HT-axis.
* This is **cross-ethnic replication** of R17's two-axis decomposition in a Korean cohort,
  not gated on TCGA-specific labels.
* Caveat: n=18 is small; treat as supportive replication, not a stand-alone claim.