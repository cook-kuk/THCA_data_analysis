# 11개 고IF 후보 주제 pilot screen — 2026-05-08

범위: 현재 CNV residual class 포함 + 추가 10개 후보. 기존 원고/figure는 수정하지 않았고, 로컬 public-data-derived 결과만 사용.

## 최종 순위

| Rank | Topic | Grade | Score | 한 줄 판정 | Next gate |
|---:|---|---:|---:|---|---|
| 1 | T02 · T02 tumor-specific TROP2 spatial niche | A- | 90 | Bulk-negative TROP2 is rescued as a tumor-specific spatial niche rather than a whole-tumor expression marker. | Keep as spatial-niche biology, not ADC vulnerability; add GeoMx/IF/IHC if pursuing high-tier. |
| 2 | T00 · T00 CNV-defined true-driver-negative residual class | B+ | 84 | True driver-negative thyroid cancers contain a recurrent arm-CNV residual class concentrated in DM2/FVPTC-like tumors. | Validate from independent GDC/cBioPortal GISTIC source; adjust for FVPTC/Hurthle histology and purity. |
| 3 | T01 · T01 miRNA-biogenesis-defective progression route (DICER1/DGCR8 GeoMx) | B+ | 84 | DICER1/DGCR8-associated lesions trace a compartment-resolved route from microPTC to PDTC with lineage loss and immune/stromal divergence. | Patient-stratified model; separate DICER1 from DGCR8 and histology; add CNV/DM overlay if possible. |
| 4 | T04 · T04 CAF/ECM ligand remodeling as a suppressor context for thyroid-lineage/RAI programs | B+ | 82 | A conserved CAF/ECM ligand program tracks spatial thyroid-lineage loss and bulk RAI suppression in driver-negative tumors. | Disentangle stromal abundance from tumor-cell intrinsic loss; validate with GeoMx VIM vs PanCK compartments. |
| 5 | T07 · T07 immune/TLS protective paradox in thyroid cancer | B+ | 82 | Immune-rich/TLS-like thyroid cancers may represent a protective indolence axis rather than a generic aggressive inflammation axis. | Validate in independent survival cohort and separate Hashimoto-like immune context from tumor immune escape. |
| 6 | T03 · T03 Hashimoto-overlap PTC TLS / B-cell / IGHV spatial niche | B | 78 | Hashimoto-overlap PTC contains spatially organized TLS/B-cell/IGHV niches visible directly in tissue. | Independent PTC+HT cohort or local IF/IHC; maintain no causal HT-to-cancer language. |
| 7 | T06 · T06 NIS/SLC5A5 non-methylation silencing exception | B | 75 | RAI-lineage silencing splits into methylation-driven hormone-biosynthesis loss and non-methylation NIS/SLC5A5 loss. | Add histone/miRNA/TF motif evidence; this is mechanism-follow-up, not standalone clinical biomarker yet. |
| 8 | T08 · T08 DICER1/EIF1AX rare-driver FVPTC-like residual class | B- | 68 | DICER1/EIF1AX-mutant PTC forms a rare FVPTC-like differentiated residual class that may bridge genomics and miRNA-biogenesis biology. | Pool TCGA + K2 + GeoMx DICER1 lesion data; do not pitch as standalone until n improves. |
| 9 | T05 · T05 spatial immune-checkpoint ligand-receptor routes in aggressive thyroid cancer | B- | 68 | ATC/LPTC progression may route through spatially localized LGALS9-HAVCR2, CXCL10-CXCR3, and CCL5-CCR5 immune-checkpoint circuits. | Cross-check against the 2026 JCI Insight atlas and avoid ICI-response claims without treated cohorts. |
| 10 | T11 · T11 ATC spatial coherence collapse with single-sample super-organization | B- | 68 | ATC progression collapses thyroid differentiation spatial coherence while creating stress/immune/hypoxia super-organized regions. | Keep supplementary unless an independent ATC spatial cohort or GeoMx PDTC/ATC validation supports it. |
| 11 | T10 · T10 morphology-invisible RAI/DM molecular axis | negative-use-only | 60 | Routine morphology does not recover the RAI/DM transcriptomic axis, arguing for molecular readouts rather than H&E surrogates. | Use as negative-control/methods paper or supplement; not a wet-lab biology headline. |
| 12 | T09 · T09 TERT-only triple-negative high-risk microclass | C | 55 | A tiny TERT-only residual class may carry disproportionate recurrence risk. | External TERT-rich cohort or MSK/GENIE-level expansion required; hold as section only. |

## Topic별 pilot 요약

### T02 — T02 tumor-specific TROP2 spatial niche

**Headline:** Bulk-negative TROP2 is rescued as a tumor-specific spatial niche rather than a whole-tumor expression marker.

**Pilot result:** PTC/LPTC Visium slides show high TROP2 spatial autocorrelation versus normal/autoimmune controls.

**Next gate:** Keep as spatial-niche biology, not ADC vulnerability; add GeoMx/IF/IHC if pursuing high-tier.


| Metric | n/rate/effect | p | Note |
|---|---|---:|---|
| TROP2 Moran PTC/LPTC vs negative controls | n_a=8; n_b=12; mean_a=0.433; mean_b=0.0398; diff_a_minus_b=0.393; cohen_d_a_minus_b=4.99 | 1.6e-05 | 28-slide Visium corpus |
| TROP2 niche threshold rate >0.25 | n_a=8; n_b=12; rate_a=1; rate_b=0 | 7.9e-06 | Moran's I threshold |

### T00 — T00 CNV-defined true-driver-negative residual class

**Headline:** True driver-negative thyroid cancers contain a recurrent arm-CNV residual class concentrated in DM2/FVPTC-like tumors.

**Pilot result:** CNV signature marks Class6-DM2 far more than Class6-DM1; top arms include 7p/7q gain, 2p/2q loss, 16p/16q gain, 12q gain.

**Next gate:** Validate from independent GDC/cBioPortal GISTIC source; adjust for FVPTC/Hurthle histology and purity.


| Metric | n/rate/effect | p | Note |
|---|---|---:|---|
| Class6-DM2 vs Class6-DM1 CNV-signature rate | n_a=45; n_b=65; rate_a=0.311; rate_b=0.0154 | 9.5e-06 | 7-arm quick signature |
| 16p_Gain DM2 vs DM1 | rate_a=0.222; rate_b=0 | 7.9e-05 | top Class6 arm CNV |
| 16q_Gain DM2 vs DM1 | rate_a=0.2; rate_b=0 | 0.0002 | top Class6 arm CNV |
| 2q_Loss DM2 vs DM1 | rate_a=0.2; rate_b=0 | 0.0002 | top Class6 arm CNV |
| 2p_Loss DM2 vs DM1 | rate_a=0.2; rate_b=0 | 0.0002 | top Class6 arm CNV |
| 7q_Gain DM2 vs DM1 | rate_a=0.227; rate_b=0.0154 | 0.0005 | top Class6 arm CNV |
| 7p_Gain DM2 vs DM1 | rate_a=0.222; rate_b=0.0154 | 0.0006 | top Class6 arm CNV |
| 12q_Gain DM2 vs DM1 | rate_a=0.178; rate_b=0 | 0.0006 | top Class6 arm CNV |

### T01 — T01 miRNA-biogenesis-defective progression route (DICER1/DGCR8 GeoMx)

**Headline:** DICER1/DGCR8-associated lesions trace a compartment-resolved route from microPTC to PDTC with lineage loss and immune/stromal divergence.

**Pilot result:** GeoMx n=78 ROI supports a PDTC-vs-microPTC lineage drop and HLA/DM-axis shifts, but genotype/histology are partially confounded.

**Next gate:** Patient-stratified model; separate DICER1 from DGCR8 and histology; add CNV/DM overlay if possible.


| Metric | n/rate/effect | p | Note |
|---|---|---:|---|
| PDTC PanCK+ vs microPTC PanCK+ DM1_axis | n_a=9; n_b=5; mean_a=3.98; mean_b=4.39; diff_a_minus_b=-0.411; cohen_d_a_minus_b=-3.3 | 0.0010 | GSE301163 GeoMx |
| PDTC PanCK+ vs microPTC PanCK+ RAI_8 | n_a=9; n_b=5; mean_a=4.46; mean_b=4.83; diff_a_minus_b=-0.37; cohen_d_a_minus_b=-4.54 | 0.0010 | GSE301163 GeoMx |
| PDTC PanCK+ vs microPTC PanCK+ Eight_gene_DM1 | n_a=9; n_b=5; mean_a=2.16; mean_b=2.23; diff_a_minus_b=-0.0658; cohen_d_a_minus_b=-0.547 | 0.5185 | GSE301163 GeoMx |
| PDTC PanCK+ vs microPTC PanCK+ HLA_II | n_a=9; n_b=5; mean_a=1.98; mean_b=2.84; diff_a_minus_b=-0.865; cohen_d_a_minus_b=-6.05 | 0.0010 | GSE301163 GeoMx |
| PDTC PanCK+ vs microPTC PanCK+ Cell_cycle | n_a=9; n_b=5; mean_a=2.14; mean_b=2.2; diff_a_minus_b=-0.0567; cohen_d_a_minus_b=-0.637 | 0.2977 | GSE301163 GeoMx |
| PDTC PanCK+ vs microPTC PanCK+ Hypoxia | n_a=9; n_b=5; mean_a=3.16; mean_b=3.37; diff_a_minus_b=-0.202; cohen_d_a_minus_b=-3.13 | 0.0010 | GSE301163 GeoMx |
| PDTC PanCK+ vs microPTC PanCK+ EMT | n_a=9; n_b=5; mean_a=2.7; mean_b=2.97; diff_a_minus_b=-0.271; cohen_d_a_minus_b=-4.18 | 0.0010 | GSE301163 GeoMx |
| Eight_gene_DM1 microPTC PanCK+ vs Normal PanCK+ | n_a=5; n_b=6; mean_a=2.23; mean_b=2.55 | 0.0173 | precomputed GeoMx stat |
| DM1_axis PDTC PanCK+ vs microPTC PanCK+ | n_a=9; n_b=5; mean_a=3.98; mean_b=4.39 | 0.0010 | precomputed GeoMx stat |
| RAI_8 PDTC PanCK+ vs microPTC PanCK+ | n_a=9; n_b=5; mean_a=4.46; mean_b=4.83 | 0.0010 | precomputed GeoMx stat |
| Thyroid_TF PDTC PanCK+ vs microPTC PanCK+ | n_a=9; n_b=5; mean_a=4.73; mean_b=4.28 | 0.0070 | precomputed GeoMx stat |

### T04 — T04 CAF/ECM ligand remodeling as a suppressor context for thyroid-lineage/RAI programs

**Headline:** A conserved CAF/ECM ligand program tracks spatial thyroid-lineage loss and bulk RAI suppression in driver-negative tumors.

**Pilot result:** DCN/COL1A2/SFRP2/COL1A1 are top spatial ligands; TCGA driver-negative CAF/M1_M2 modules correlate negatively with RAI.

**Next gate:** Disentangle stromal abundance from tumor-cell intrinsic loss; validate with GeoMx VIM vs PanCK compartments.


| Metric | n/rate/effect | p | Note |
|---|---|---:|---|
| spatial ligand DCN | n=16; median_r=0.286 | NA | NicheNet-lite spatial ligand screen |
| spatial ligand COL1A2 | n=16; median_r=0.279 | NA | NicheNet-lite spatial ligand screen |
| spatial ligand SFRP2 | n=16; median_r=0.229 | NA | NicheNet-lite spatial ligand screen |
| spatial ligand COL1A1 | n=16; median_r=0.182 | NA | NicheNet-lite spatial ligand screen |
| spatial ligand COL3A1 | n=16; median_r=0.178 | NA | NicheNet-lite spatial ligand screen |
| spatial ligand LUM | n=16; median_r=0.175 | NA | NicheNet-lite spatial ligand screen |
| spatial ligand BGN | n=16; median_r=0.167 | NA | NicheNet-lite spatial ligand screen |
| spatial ligand SFRP4 | n=15; median_r=0.144 | NA | NicheNet-lite spatial ligand screen |
| TCGA NEG CAF vs RAI | n=132; spearman_r=-0.285 | NA | BRAF/RAS-negative bulk replication |
| TCGA NEG M1_M2 vs RAI | n=132; spearman_r=-0.479 | NA | BRAF/RAS-negative bulk replication |
| TCGA NEG TLS vs RAI | n=132; spearman_r=-0.253 | NA | BRAF/RAS-negative bulk replication |
| TCGA NEG CD8 vs RAI | n=132; spearman_r=-0.267 | NA | BRAF/RAS-negative bulk replication |

### T07 — T07 immune/TLS protective paradox in thyroid cancer

**Headline:** Immune-rich/TLS-like thyroid cancers may represent a protective indolence axis rather than a generic aggressive inflammation axis.

**Pilot result:** DSS HRs for M1_M2/CD8 are protective in TCGA; TLS trends protective but needs morphology-level validation.

**Next gate:** Validate in independent survival cohort and separate Hashimoto-like immune context from tumor immune escape.


| Metric | n/rate/effect | p | Note |
|---|---|---:|---|
| DSS Cox M1_M2 | n=566; HR=0.347 | 0.0050 | univariate TCGA module Cox |
| DSS Cox CD8 | n=566; HR=0.277 | 0.0111 | univariate TCGA module Cox |
| DSS Cox TLS | n=566; HR=0.425 | 0.0737 | univariate TCGA module Cox |
| DSS Cox CAF | n=566; HR=1.15 | 0.7331 | univariate TCGA module Cox |
| DSS Cox RAI | n=566; HR=1.7 | 0.2578 | univariate TCGA module Cox |
| DSS Cox age+stage M1_M2 | n=564; HR=0.53 | 0.0362 | adjusted TCGA module Cox |
| DSS Cox age+stage CD8 | n=564; HR=0.618 | 0.1321 | adjusted TCGA module Cox |
| DSS Cox age+stage TLS | n=564; HR=0.658 | 0.1676 | adjusted TCGA module Cox |

### T03 — T03 Hashimoto-overlap PTC TLS / B-cell / IGHV spatial niche

**Headline:** Hashimoto-overlap PTC contains spatially organized TLS/B-cell/IGHV niches visible directly in tissue.

**Pilot result:** 4/4 GSE230424 slides show organized HLA-II/B-cell/TLS/IGHV Moran patterns; n=4 is the main limitation.

**Next gate:** Independent PTC+HT cohort or local IF/IHC; maintain no causal HT-to-cancer language.


| Metric | n/rate/effect | p | Note |
|---|---|---:|---|
| GSE230424 HLA_II Moran summary | n=4 | NA | 4 PTC+HT slides |
| GSE230424 B_cell Moran summary | n=4 | NA | 4 PTC+HT slides |
| GSE230424 TLS Moran summary | n=4 | NA | 4 PTC+HT slides |
| GSE230424 IGHV_AICDA Moran summary | n=4 | NA | 4 PTC+HT slides |
| PTC_HT immune-Moran mean vs others | n_a=16; n_b=96; mean_a=0.567; mean_b=0.35; diff_a_minus_b=0.217 | NA | axis-level summary |

### T06 — T06 NIS/SLC5A5 non-methylation silencing exception

**Headline:** RAI-lineage silencing splits into methylation-driven hormone-biosynthesis loss and non-methylation NIS/SLC5A5 loss.

**Pilot result:** TPO/TG/TSHR/PAX8 show strong methylation effects; SLC5A5/NIS is the outlier with nonsignificant promoter methylation.

**Next gate:** Add histone/miRNA/TF motif evidence; this is mechanism-follow-up, not standalone clinical biomarker yet.


| Metric | n/rate/effect | p | Note |
|---|---|---:|---|
| SLC5A5 promoter methylation DM1-DM2 | n_a=90; n_b=54; diff_a_minus_b=0.026; cohen_d_a_minus_b=0.225 | 0.4209 | HM450 per-gene methylation |
| TPO promoter methylation DM1-DM2 | n_a=90; n_b=54; diff_a_minus_b=0.415; cohen_d_a_minus_b=2.3 | 1.9e-18 | HM450 per-gene methylation |
| TG promoter methylation DM1-DM2 | n_a=90; n_b=54; diff_a_minus_b=0.129; cohen_d_a_minus_b=0.856 | 2.2e-06 | HM450 per-gene methylation |
| TSHR promoter methylation DM1-DM2 | n_a=90; n_b=54; diff_a_minus_b=0.117; cohen_d_a_minus_b=1.2 | 9.8e-12 | HM450 per-gene methylation |
| PAX8 promoter methylation DM1-DM2 | n_a=90; n_b=54; diff_a_minus_b=0.045; cohen_d_a_minus_b=0.97 | 4.5e-08 | HM450 per-gene methylation |
| NKX2-1 promoter methylation DM1-DM2 | n_a=90; n_b=54; diff_a_minus_b=0.053; cohen_d_a_minus_b=0.631 | 8.9e-07 | HM450 per-gene methylation |
| FOXE1 promoter methylation DM1-DM2 | n_a=90; n_b=54; diff_a_minus_b=0.042; cohen_d_a_minus_b=0.835 | 1.0e-05 | HM450 per-gene methylation |
| DIO1 promoter methylation DM1-DM2 | n_a=90; n_b=54; diff_a_minus_b=0.227; cohen_d_a_minus_b=1.24 | 6.5e-11 | HM450 per-gene methylation |
| TPO methylation effect minus SLC5A5 exception | - | NA | mechanism split |

### T08 — T08 DICER1/EIF1AX rare-driver FVPTC-like residual class

**Headline:** DICER1/EIF1AX-mutant PTC forms a rare FVPTC-like differentiated residual class that may bridge genomics and miRNA-biogenesis biology.

**Pilot result:** TCGA n=7 is small but heavily DM2-like with high RAI/TDS and low methylation; K2 rare-driver evidence exists but needs harmonization.

**Next gate:** Pool TCGA + K2 + GeoMx DICER1 lesion data; do not pitch as standalone until n improves.


| Metric | n/rate/effect | p | Note |
|---|---|---:|---|
| Class4 vs Class6 DM2 rate | n_a=7; n_b=125; rate_a=0.857; rate_b=0.384 | 0.0186 | TCGA driver class |
| Class4 vs Class6 rai_score_v17 | n_a=7; n_b=125; mean_a=9.03; mean_b=8.33; diff_a_minus_b=0.704; cohen_d_a_minus_b=0.589 | 0.1437 | TCGA driver class |
| Class4 vs Class6 tds16_score_v17 | n_a=7; n_b=125; mean_a=7.96; mean_b=7.41; diff_a_minus_b=0.551; cohen_d_a_minus_b=0.574 | 0.1203 | TCGA driver class |
| Class4 vs Class6 mean_8g_beta | n_a=7; n_b=123; mean_a=0.237; mean_b=0.336; diff_a_minus_b=-0.0994; cohen_d_a_minus_b=-1.13 | 0.0034 | TCGA driver class |
| Class4 vs Class6 age | n_a=7; n_b=122; mean_a=52.6; mean_b=46.4; diff_a_minus_b=6.23; cohen_d_a_minus_b=0.361 | 0.5923 | TCGA driver class |

### T05 — T05 spatial immune-checkpoint ligand-receptor routes in aggressive thyroid cancer

**Headline:** ATC/LPTC progression may route through spatially localized LGALS9-HAVCR2, CXCL10-CXCR3, and CCL5-CCR5 immune-checkpoint circuits.

**Pilot result:** ATC dominates the top ligand-receptor activity rows, but this is spatial-signature only and sample-limited.

**Next gate:** Cross-check against the 2026 JCI Insight atlas and avoid ICI-response claims without treated cohorts.


| Metric | n/rate/effect | p | Note |
|---|---|---:|---|
| ATC LGALS9→HAVCR2 | n=4 | NA | spatial ligand-receptor activity |
| LPTC LGALS9→HAVCR2 | n=4 | NA | spatial ligand-receptor activity |
| ATC CCL5→CCR5 | n=4 | NA | spatial ligand-receptor activity |
| LPTC CCL5→CCR5 | n=4 | NA | spatial ligand-receptor activity |
| ATC CXCL10→CXCR3 | n=4 | NA | spatial ligand-receptor activity |
| PTC LGALS9→HAVCR2 | n=4 | NA | spatial ligand-receptor activity |
| PTC CXCL13→CXCR5 | n=4 | NA | spatial ligand-receptor activity |
| ATC CD86→CTLA4 | n=4 | NA | spatial ligand-receptor activity |
| PTC CXCL9→CXCR3 | n=4 | NA | spatial ligand-receptor activity |
| PT CCL5→CCR5 | n=4 | NA | spatial ligand-receptor activity |
| PTC CD86→CTLA4 | n=4 | NA | spatial ligand-receptor activity |
| ATC IL10→IL10RA | n=4 | NA | spatial ligand-receptor activity |

### T11 — T11 ATC spatial coherence collapse with single-sample super-organization

**Headline:** ATC progression collapses thyroid differentiation spatial coherence while creating stress/immune/hypoxia super-organized regions.

**Pilot result:** Spatial v12 shows differentiation collapse in ATC and an outlier-like highly organized ATC sample, but n=4 ATC and n=1 outlier caveats dominate.

**Next gate:** Keep supplementary unless an independent ATC spatial cohort or GeoMx PDTC/ATC validation supports it.


| Metric | n/rate/effect | p | Note |
|---|---|---:|---|
| ATC vs LPTC stage-gradient Thyroid_differen | n_a=1; n_b=1; mean_a=0.124; mean_b=0.612; diff_a_minus_b=-0.488 | 1.0000 | 16-slide GSE250521 spatial v12 |
| ATC vs LPTC stage-gradient Cell_cycle | n_a=1; n_b=1; mean_a=0.437; mean_b=0.161; diff_a_minus_b=0.276 | 1.0000 | 16-slide GSE250521 spatial v12 |
| ATC vs LPTC stage-gradient T_cell | n_a=1; n_b=1; mean_a=0.419; mean_b=0.182; diff_a_minus_b=0.237 | 1.0000 | 16-slide GSE250521 spatial v12 |
| ATC vs LPTC stage-gradient Macrophage_TAM | n_a=1; n_b=1; mean_a=0.346; mean_b=0.214; diff_a_minus_b=0.132 | 1.0000 | 16-slide GSE250521 spatial v12 |
| ATC vs LPTC stage-gradient Hypoxia | n_a=1; n_b=1; mean_a=0.496; mean_b=0.526; diff_a_minus_b=-0.03 | 1.0000 | 16-slide GSE250521 spatial v12 |
| ATC vs LPTC stage-gradient EMT_program | n_a=1; n_b=1; mean_a=0.429; mean_b=0.507; diff_a_minus_b=-0.078 | 1.0000 | 16-slide GSE250521 spatial v12 |
| ATC vs LPTC stage-gradient Cellular_stress | n_a=1; n_b=1; mean_a=0.569; mean_b=0.598; diff_a_minus_b=-0.0288 | 1.0000 | 16-slide GSE250521 spatial v12 |
| top ATC super-organization sample | - | NA | case-like outlier screen |

### T10 — T10 morphology-invisible RAI/DM molecular axis

**Headline:** Routine morphology does not recover the RAI/DM transcriptomic axis, arguing for molecular readouts rather than H&E surrogates.

**Pilot result:** Closure battery stays near random (max |Spearman| < 0.1; AUROC near 0.55).

**Next gate:** Use as negative-control/methods paper or supplement; not a wet-lab biology headline.


| Metric | n/rate/effect | p | Note |
|---|---|---:|---|
| Closure battery best absolute Spearman | max_abs_spearman=0.0929; max_auroc=0.547 | NA | H&E/Visium closure battery |
| A raw_DM1_like_score ridge | spearman_r=0.0606 | NA | negative result |
| A raw_DM1_like_score enet | spearman_r=-0.0929 | NA | negative result |
| A raw_RAI_8_score ridge | spearman_r=0.0606 | NA | negative result |
| A raw_RAI_8_score enet | spearman_r=-0.0929 | NA | negative result |
| A raw_TDS_like_score ridge | spearman_r=0.079 | NA | negative result |
| A raw_TDS_like_score enet | spearman_r=-0.0792 | NA | negative result |
| B B1_resid_log_counts_only ridge | spearman_r=0.0118 | NA | negative result |
| B B2_resid_log_ngenes_only ridge | spearman_r=0.0113 | NA | negative result |

### T09 — T09 TERT-only triple-negative high-risk microclass

**Headline:** A tiny TERT-only residual class may carry disproportionate recurrence risk.

**Pilot result:** PFI event rate appears high but TCGA n=5 is too small for a paper.

**Next gate:** External TERT-rich cohort or MSK/GENIE-level expansion required; hold as section only.


| Metric | n/rate/effect | p | Note |
|---|---|---:|---|
| Class5 TERT-only vs Class6 PFI event rate | n_a=5; n_b=125; rate_a=0.4; rate_b=0.048 | 0.0304 | TCGA sparse pilot |
| Class5 vs Class6 rai_score_v17 | n_a=5; n_b=125; mean_a=6.79; mean_b=8.33; diff_a_minus_b=-1.54; cohen_d_a_minus_b=-1.26 | 0.0131 | TCGA sparse pilot |
| Class5 vs Class6 mean_8g_beta | n_a=5; n_b=123; mean_a=0.356; mean_b=0.336; diff_a_minus_b=0.0197; cohen_d_a_minus_b=0.223 | 0.6402 | TCGA sparse pilot |
| Class5 vs Class6 age | n_a=5; n_b=122; mean_a=54.4; mean_b=46.4; diff_a_minus_b=8.02; cohen_d_a_minus_b=0.477 | 0.2892 | TCGA sparse pilot |

## 산출 파일

- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/topic_scorecard.tsv`

- `/home/seungho/personal/THCA_data_analysis/project/results/high_impact_topic_pilots_2026_05_08/pilot_metrics_long.tsv`

- `/home/seungho/personal/THCA_data_analysis/project/reports/2026_05_08_high_impact_11_topic_pilots.md`
