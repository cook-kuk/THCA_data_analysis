# 추가 public data addendum — 2026-05-08

사용자 요청에 따라 11-topic pilot 중 추가 public 데이터를 더 가져와서 바로 붙일 수 있는 검증 레이어를 확인했다.

## 새로 가져온 데이터

| Source | What | Local file | Use |
|---|---|---|---|
| cBioPortal TCGA-THCA Cell 2014 | sample clinical, ARMDRIVER_CN, ARM_SCNA_CLUSTER, selected-gene GISTIC | `data_search/cbio_thca_tcga_pub_*` | T00 CNV residual class independent validation |
| cBioPortal TCGA-GDC 2025 | selected-gene GISTIC, hg38 GDC mirror | `data_search/thpa_tcga_gdc_*` | T00 consistency check |
| cBioPortal GATCI 2024 ATC | n=190 ATC GISTIC + clinical | `data_search/thyroid_gatci_2024_*` | T11 advanced CNA progression context |
| cBioPortal MSK 2016 PDTC/ATC | n=117 targeted-panel GISTIC + clinical | `data_search/thyroid_mskcc_2016_*` | negative/sparse targeted-panel context |
| Web snapshots | JCI Insight 2026 atlas + GSE301163 OmicsDI page | `data_search/external_source_snapshots.json` | prior-art guardrail |

## T00 CNV residual class — cBioPortal validation

| Metric | DM2 rate | DM1 rate | OR | p |
|---|---:|---:|---:|---:|
| cbio_ARMDRIVER_signature_any | 0.354 | 0.0132 | 41.1 | 1.6e-07 |
| cbio_ARMDRIVER_12q | 0.25 | 0 | inf | 4.4e-06 |
| cbio_ARMDRIVER_7q | 0.292 | 0.0132 | 30.9 | 4.7e-06 |
| cbio_ARMDRIVER_16p | 0.229 | 0 | inf | 1.3e-05 |
| ARM_SCNA_CLUSTER=Many SCNA | 0.271 | 0.0132 | 27.9 | 1.4e-05 |
| cbio_ARMDRIVER_7p | 0.25 | 0.0132 | 25 | 4.0e-05 |
| cbio_ARMDRIVER_2p | 0.208 | 0 | inf | 4.0e-05 |
| cbio_ARMDRIVER_16q | 0.188 | 0 | inf | 0.000118 |
| cbio_ARMDRIVER_2q | 0.167 | 0 | inf | 0.000343 |
| ARM_SCNA_CLUSTER=Quiet | 0.5 | 0.724 | 0.382 | 0.0136 |

해석: local arm-CNV table뿐 아니라 cBioPortal의 TCGA-THCA sample-level `ARMDRIVER_CN`에서도 같은 7-arm signature가 DM2에 강하게 몰린다. T00은 `B+`에서 `A-`로 상향 가능하다.

## T00 selected-gene GISTIC support within Class6

| Gene | Event | DM2 rate | DM1 rate | OR | p |
|---|---|---:|---:|---:|---:|
| DCN | gain_or_amp | 0.255 | 0 | inf | 8.5e-06 |
| MET | gain_or_amp | 0.298 | 0.0147 | 28.4 | 1.0e-05 |
| CDK6 | gain_or_amp | 0.277 | 0.0147 | 25.6 | 2.8e-05 |
| BRAF | gain_or_amp | 0.298 | 0.0294 | 14 | 5.5e-05 |
| HAVCR2 | gain_or_amp | 0.213 | 0 | inf | 6.9e-05 |
| TERT | gain_or_amp | 0.213 | 0 | inf | 6.9e-05 |
| EGFR | gain_or_amp | 0.255 | 0.0147 | 23 | 7.6e-05 |
| TPO | loss_or_del | 0.191 | 0 | inf | 0.000194 |
| COL1A1 | gain_or_amp | 0.213 | 0.0147 | 18.1 | 0.000519 |
| CCL5 | gain_or_amp | 0.213 | 0.0147 | 18.1 | 0.000519 |
| PAX8 | loss_or_del | 0.17 | 0 | inf | 0.000532 |
| LGALS9 | gain_or_amp | 0.191 | 0.0147 | 15.9 | 0.00131 |
| DICER1 | gain_or_amp | 0.149 | 0 | inf | 0.00144 |
| NKX2-1 | gain_or_amp | 0.149 | 0 | inf | 0.00144 |
| DIO2 | gain_or_amp | 0.149 | 0 | inf | 0.00144 |
| CXCR3 | gain_or_amp | 0.128 | 0 | inf | 0.00382 |
| SLC5A5 | gain_or_amp | 0.128 | 0 | inf | 0.00382 |
| ALK | loss_or_del | 0.191 | 0.0294 | 7.82 | 0.00704 |

해석: 7q/7p representative genes (MET, CDK6, BRAF, EGFR)와 12q/16p 관련 selected genes가 Class6-DM2에서 gene-level GISTIC gain으로도 반복된다. 단 gene-level은 arm-level의 proxy라서 main claim은 arm-CNV로 유지.

## Advanced thyroid cancer external context

| Study | Top selected CNA gains | Interpretation |
|---|---|---|
| thyroid_gatci_2024 | EGFR 38.7%, MET 34.3%, TERT 32.1%, TG 29.2%, CDK6 27.0%, BRAF 25.5%, NKX2-1 23.4%, COL1A1 20.4% | ATC에서 EGFR/MET/TERT/CDK6/BRAF gains가 흔해 advanced-CNA progression context를 제공. |
| thyroid_mskcc_2016 | NKX2-1 4.3%, TERT 0.9%, ALK 0.9%, CDK6 0.0%, BRAF 0.0%, DICER1 0.0%, EGFR 0.0%, MET 0.0% | Targeted-panel GISTIC가 매우 sparse; negative control/coverage caveat로만 사용. |
| thpa_tcga_gdc | TERT 5.5%, MET 5.1%, BRAF 5.1%, CDK6 5.0%, CCL5 5.0%, COL1A1 5.0%, HAVCR2 5.0%, LGALS9 4.8% | TCGA-GDC mirror에서 TCGA selected-gene gain frequency가 비슷한 방향으로 확인. |

## 반영된 판정

| Topic | Old | New | Reason |
|---|---:|---:|---|
| T00 | B+ | A- | cBioPortal ARMDRIVER_CN independently supports 7-arm Class6-DM2 enrichment: 35.4% vs 1.3%, Fisher p=1.6e-7. |
| T01 | B+ | B+ | GSE301163 already included; extra cBio data adds context but not independent DICER1/DGCR8 class validation. |
| T11 | B- | B | GATCI 2024 ATC GISTIC shows frequent EGFR/MET/TERT/CDK6/BRAF gains, supporting an advanced-CNA progression layer, but not the same PTC residual-class test. |

## Files

- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/external_cna_selected_gene_summary.tsv`

- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/topic_grade_updates_after_extra_data.tsv`

- `/home/seungho/personal/THCA_data_analysis/project/reports/2026_05_08_high_impact_extra_data_addendum.md`
