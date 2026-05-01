# v17 HLA analysis — does HLA reveal NEW biology beyond DM1/DM2 + BRAF/RAS/TERT?

**Author**: Seungho Cook  **Date**: 2026-04-28  **random_state**: 42

## Hypotheses

- **H1** HLA Class I (HLA-A/B/C + B2M + TAP1/2 + NLRC5) suppressed in BRAF V600E (literature: MAPK -> HLA-I down, immune evasion).
- **H2** HLA Class I/II higher in DM1 (immune-hot) than DM2 (immune-cold).
- **H3** HLA Class II tracks lymphocytic infiltration (Hashimoto/autoimmune-PTC axis).
- **H4** HLA scores decline Normal -> PTC -> PDFP -> ATC in Korean GSE213647.
- **H5** HLA-A 4-digit genotype distribution differs by molecular subtype (Thorsson 2018 OptiType calls).

## Data sources

- TCGA-THCA log2 expression: `/data/thca/data_processed/bulk_rnaseq_v3/TCGA-THCA_v3_log2.tsv` (572 samples, all 14 HLA genes present)
- DM1/DM2 labels: `/opt/thyroid-dash/project/results/v17p3/tables/A2_dm_score_full_cohort.tsv` (DM1=403, DM2=110)
- BRAF/TERT calls: `/opt/thyroid-dash/project/results/v17_tert_recovery/v3/v3_external_combined_BRAF_TERT.tsv` (V600E=282, RAS-only=54, TripleNeg=161)
- Korean GSE213647 counts: `/data/thca/v17_korean/GSE213647/expression_counts.tsv.gz` (632 samples, all 14 HLA ENSG present)

## Results

### H1 — HLA Class I in BRAF V600E (TCGA-THCA, n=569)

| group | n | median HLA-I z |
|---|---|---|
| BRAF V600E+ | 319 | **+0.231** |
| BRAF V600E- (RAS / TripleNeg) | 250 | -0.552 |

Mann-Whitney U=54,755, **p=2.1e-14, Cohen's d=+0.634**. Class II same direction (d=+0.602, p=8.1e-12).

**Verdict**: H1 is REJECTED in its literature-classical form. In TCGA-THCA, BRAF V600E PTCs have HIGHER (not lower) HLA expression. This is the well-known "BRAF-PTC = immune-infiltrated, lymphocyte-rich" phenotype — in differentiated PTC, BRAF V600E recruits immune cells and the bulk HLA signal reflects the infiltrate, not tumor-cell-intrinsic HLA. To rescue the literature signal would need single-cell or HLA-loss-of-heterozygosity inference.

### H2 — HLA in DM1 vs DM2 (TCGA-THCA, n=517)

| score | DM1 median | DM2 median | p | Cohen's d |
|---|---|---|---|---|
| HLA-I  | +0.237 | -0.952 | **1.6e-34** | **+1.53** |
| HLA-II | +0.306 | -1.031 | **8.1e-37** | **+1.75** |

Spearman rho(prob_DM1, HLA-I) = **0.588** (p=5e-49); rho(prob_DM1, HLA-II) = **0.619** (p=1e-55).

**Verdict**: H2 STRONGLY SUPPORTED — HLA expression is one of the largest known single-axis separators of DM1 vs DM2 (d~1.5-1.75). But this is **redundant** with the existing DM1=immune-hot definition: HLA-I/II loads on the same axis as the 8-gene panel anti-correlation in DM1 vs DM2.

### H3 — HLA Class II vs lymphocyte signature (TCGA-THCA, n=572)

| pair | Spearman rho | p |
|---|---|---|
| HLA-II vs lymph (CD8A,B,CD3DEG,GZMA/B,PRF1,IFNG,CXCL9/10) | **0.795** | 1e-125 |
| HLA-I  vs lymph | 0.784 | 6e-120 |

**Verdict**: HLA expression is ~80% explained by infiltrating lymphocyte abundance. HLA-II does NOT cleanly separate "Hashimoto-autoimmune-PTC" from "tumor-rejection cytotoxic infiltrate" using bulk RNA — both look identical. To carve out a true autoimmune-PTC subgroup orthogonal to DM, we would need (a) thyroglobulin/TPO autoantibody titre or (b) lymphocytic thyroiditis pathology annotation, neither in TCGA freeze.

### H4 — Korean GSE213647 by histology (n=632)

| histology | n | median HLA-I z | median HLA-II z |
|---|---|---|---|
| Normal | 262 | -0.458 | -0.588 |
| PTC | 353 | +0.243 | +0.397 |
| PDFP | 9 | -0.400 | +0.007 |
| UTC/ATC | 8 | -0.014 | -0.042 |

Kruskal-Wallis: HLA-I p=1.7e-21, HLA-II p=4.6e-28. PTC > Normal (p~1e-23 / 1e-29). Normal vs PDFP/ATC NOT significant (small n=9, n=8).

**Verdict**: H4 monotonic dedifferentiation pattern NOT seen — bulk HLA goes UP in PTC vs Normal (lymphocyte infiltration), then plateau in PDFP/ATC. The 8-gene thyroid-differentiation panel (TG, SLC5A5, TPO, ...) gives the clean monotonic loss; HLA does not.

### H5 — HLA-A/B/C 4-digit genotype

DEFERRED. Thorsson 2018 PanCanAtlas OptiType calls (`TCGA.HLA.4digit.byPatient.txt`) live on Synapse syn4602499 and require Synapse credentials. Public cBioPortal `thca_tcga_pan_can_atlas_2018` clinical export does not include them (verified — only IMMUNE_SUBTYPE-style attributes absent for THCA). Recommended next step: register Synapse account, re-run Fisher exact for HLA-A*02:01 / DRB1*03:01 (Hashimoto risk allele) against DM1/DM2 and BRAF V600E.

## Figures

- `reports/html/figs_interactive/v17/v17_hla_braf_split.html` — violin HLA-I by BRAF V600E
- `reports/html/figs_interactive/v17/v17_hla_dm1dm2.html` — split-violin HLA-I/II by DM1/DM2
- `reports/html/figs_interactive/v17/v17_hla_korean_histology.html` — box HLA-I/II by histology (Korean)

## Conclusion — is HLA worth a paper section?

**HONEST verdict: REDUNDANT with DM1/DM2 + lymphocyte axis. Not a standalone new finding.**

1. The HLA-I and HLA-II composite scores correlate r~0.6 with the DM1 probability and r~0.8 with a generic CD8/T-cell signature. They live on the same "immune-hot vs cold" axis already covered by DM1/DM2 and FIX5 single-cell breakdown.
2. The classical BRAF -> HLA-I-loss immune-evasion result REVERSES in bulk TCGA-THCA because BRAF V600E PTCs are lymphocyte-rich. Stating this would contradict v17's "BRAF-MAPK immune-cold" narrative unless we explicitly separate tumor-intrinsic vs infiltrate HLA. Doing that requires single-cell HLA deconvolution (out of scope for v17 timeline).
3. HLA Class II does NOT carve out an autoimmune-PTC subgroup orthogonal to DM — it is collinear with the cytotoxic infiltrate signal. A real autoimmune-PTC carve-out requires anti-Tg/anti-TPO serology or pathology lymphocytic thyroiditis annotation, neither in the freeze.
4. Genotype-side analysis (HLA-A*02:01, DRB1*03/04/05 Hashimoto alleles) is genuinely orthogonal to expression-side and could yield a new finding, but Synapse credentialed access is the gating step.

**Recommendation for the v17 npj submission**: keep HLA out of the main text. Use 1 supplementary paragraph + 1 supplementary figure (the DM1/DM2 split-violin) to demonstrate that "DM1/DM2 captures the HLA-immune axis" — i.e. cite HLA as VALIDATION of the DM1/DM2 story, not as a new axis. Re-open the genotype angle once Synapse access lands, target a follow-on Korean-cohort paper.

## Suggested 1-paragraph paper text (Korean)

> HLA Class I (HLA-A/B/C, B2M, TAP1/2, NLRC5) 및 Class II (HLA-DRA/DRB1/DPA1/DPB1/DQA1/DQB1, CIITA) 발현 점수는 TCGA-THCA에서 DM1 vs DM2 dark-matter cluster를 강하게 분리하였다 (Class I Cohen's d=1.53, p=1.6e-34; Class II d=1.75, p=8.1e-37). 그러나 두 점수는 모두 종양 침윤 림프구 지표(CD8A/B, CD3, GZMA/B, IFNG, CXCL9/10)와 매우 강한 상관(Spearman rho=0.78-0.80)을 보여, HLA 신호는 본 8-유전자 패널과 DM1/DM2 축이 이미 포착하고 있는 면역-hot/cold 축의 일부로 해석된다. BRAF V600E PTC에서 HLA Class I이 오히려 상승하는 결과는 잘 알려진 BRAF-PTC의 림프구 침윤 표현형을 반영하는 것으로, 종양세포 내재적 HLA 손실은 본 bulk RNA 데이터로 분리되지 않는다. HLA-A/B/C 4-digit 유전형 기반의 분석(Hashimoto 연관 DRB1*03/04/05 등)은 PanCanAtlas Synapse 접근이 확보되는 대로 후속 보고할 예정이다.
