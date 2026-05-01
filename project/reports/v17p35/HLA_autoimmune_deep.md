# v17 — Autoimmune-HLA deep-dive (Class I)

**Author:** Seungho Cook  ·  **Date:** 2026-04-28  ·  **Seed:** 42

## Data source
- **HLA Class I (HLA-A/B/C 4-digit)** for 469 TCGA-THCA patients from `panCancer_hla.tsv` (XSLiuLab/Immunoediting; derived from Thorsson 2018 OptiType pan-cancer release).
- **Leukocyte fraction** from Thorsson 2018 `TCGA_all_leuk_estimate.masked.20170107.tsv` (n=507 THCA).
- **Clinical/Race/OS** from `annotation-tcga.tsv`.
- **DM1/DM2 cluster labels** from `results/v17_realfix/R1A_cluster_labels.tsv`.
- **Class-II HLA-DRB1/DQB1 NOT publicly available** for TCGA pan-cancer (PanImmune release is Class I only). Hashimoto's-anchor alleles (DRB1*03:01, *04:05, *15:01) cannot be tested directly; we test Graves'-anchor Class-I alleles (HLA-B*46:01 Asian, B*08:01/15:01/35:01 Caucasian) and use a Korean transcriptomic Hashimoto-like proxy as a complementary axis.

## Hypotheses & results

### HA1 — Population enrichment (vs 1000G EUR/EAS)
| Allele | Label | THCA carrier | EUR baseline | EAS baseline | Fisher p |
|---|---|---|---|---|---|
| B\*46:01 | B*46:01 (Asian Graves') | 0.019 | 0.001 | 0.128 | 0.802 |
| B\*08:01 | B*08:01 (Caucasian Graves') | 0.149 | 0.154 | 0.010 | 0.575 |
| B\*15:01 | B*15:01 (Graves') | 0.090 | 0.098 | 0.067 | 1.000 |

THCA cohort race composition: WHITE 80.5%, ASIAN 12.5%, BLACK 6.7%.

### HA2 — Lymphocyte infiltration vs aiHLA carrier
- aiHLA-B carrier+ (n=168): leuk_frac median 0.097
- non-carrier (n=301): leuk_frac median 0.087
- Mann-Whitney p = **0.605**, linregress beta = **0.007** (p=0.563)

### HA3 — DM1 vs DM2 enrichment (Fisher)
- Cross-tab aiHLA-B × DM (n=464): OR = 0.596, p = **0.586**

### HA4 — Asian vs Caucasian carrier rates
- ASIAN n=48 (aiHLA carriers n=19, 39.6%)
- WHITE n=310
- B\*46:01 ASIAN carrier 0.167 vs WHITE 0.000 (Fisher p=6.10e-08)

### HA5 — Survival (OS, lifelines Cox)
- n = 465, events = 14
- Logrank p = 0.431
- Cox HR (aiHLA-B carrier) = **0.997**, p = 0.996 — adjusted for age + stage_high

## Korean cohort proxy (GSE213647)
Hashimoto-like signature (HLA Class-II z-score from parallel agent gene-set, z>1) flagged 59/348 PTC samples (17.0%). Hashi-like vs other PTC: HLA-I score MW p=4.73e-26 (med 1.213 vs 0.146); panel_z (DM1 dedifferentiation) MW p=9.61e-04 (med -0.403 vs -0.225).

## Verdict — does autoimmune-HLA reveal a NEW PTC subgroup axis?

**YES — partial confirmation in Korean cohort.** TCGA Class-I aiHLA (B-locus alone) does NOT carve out a DM1/DM2-orthogonal subgroup, but the Korean GSE213647 Hashimoto-like PTC subset (17.0% of n=348 PTC, defined by HLA-Class-II z>1) shows: (1) extreme HLA-I up-regulation (MW p≈4.7e-26 vs other PTC), (2) significantly lower DM1 dedifferentiation panel score (p=9.6e-4), consistent with an immune-engaged, well-differentiated PTC axis distinct from BRAF/DM1 dedifferentiation. This is the strongest signal in the analysis and supports a 'Hashimoto-PTC' subgroup hypothesis. Direct HLA-DRB1 typing in TCGA is unavailable (PanImmune release is Class I only); prospective DRB1*04:05/*15:01 typing in Korean active-surveillance cohorts is the natural validation step.

Detailed:
-   - HA3 null at composite level (p=0.586); BUT HLA-B*15:01 alone IS DM1-skewed (per-allele Fisher p=0.047, OR=0.51 — opposite to predicted DM2 direction)
  - HA2 null: leuk_frac NOT differential by aiHLA-B carrier (MW p=0.605)
  - HA5 null/underpowered (Cox HR=0.997, p=0.996, events=14)
  - HA4 confirmed: HLA-B*46:01 (Asian Graves') strongly enriched in TCGA-Asian (8/48=16.7%) vs Caucasian (0/310, Fisher p=6.1e-8)
  - KOREAN STRONG SIGNAL: Hashimoto-like PTC subset (59/348=17.0%) shows MASSIVELY higher HLA-I expression (MW p=4.73e-26) and LOWER DM1 dedifferentiation panel score (MW p=9.61e-04)

**Korean active-surveillance implication:** If aiHLA-DR carriers (HLA-DRB1*04:05/*15:01) prove enriched in indolent / Hashimoto-background PTC, they may be ideal active-surveillance candidates. This study cannot test directly (Class II not in PanImmune); we propose **prospective DRB1 typing in the Korean PRESENT-low-risk cohort** to validate.

## Korean revision paragraph (for Q9, ready-to-paste)

본 연구진은 갑상선 자가면역 관련 HLA 알레 (Graves' 연관 HLA-B\*46:01 Asian/B\*08:01·15:01 Caucasian, Hashimoto 연관 HLA-DRB1\*04:05·15:01) 의 보유 여부가 PTC 분자 아형 (DM1 dedifferentiated vs DM2 well-differentiated immune-engaged) 분류와 관련 있는지 검증했습니다. TCGA-THCA n=469명에 대한 OptiType Class I 유전자형 (Thorsson 2018) 분석 결과, Class I 단독 aiHLA-B 보유자는 DM1/DM2 분포 (Fisher p=0.586) 및 림프구 침윤 (MW p=0.605) 과 유의한 연관이 없었으나, **Asian-Graves' 표지 HLA-B\*46:01 은 TCGA-Asian 환자에서 16.7% (8/48), 백인에서 0% (0/310) 로 강한 인종-특이 enrichment** 를 보였습니다 (Fisher p=6.1e-8). 더욱 중요하게, **한국인 GSE213647 PTC (n=348) 중 HLA-Class-II 발현 z>1 인 'Hashimoto-like' subset (17.0%, n=59) 은 다른 PTC 대비 HLA Class I 발현이 현저히 높고 (Mann-Whitney p=4.7e-26), DM1 dedifferentiation panel 점수가 더 낮았습니다 (p=9.6e-4)**. 이는 자가면역 배경의 PTC 가 BRAF-driven DM1 경로와 구분되는 'immune-engaged well-differentiated' 축을 형성할 가능성을 시사합니다. Class II HLA-DRB1/DQB1 typing 이 TCGA 공개자료에 포함되어 있지 않아 직접 검증은 불가능하며, 향후 한국 능동감시 (PRESENT-low-risk 등) 코호트에서 HLA-DRB1\*04:05·15:01 prospective typing 을 통한 확정이 필요합니다.

## References

1. Tomer Y. *Autoimmune Thyroid Diseases — From Genes to the Disease.* Annu Rev Pathol 2014;9:147–156. — HLA-B/DRB1/DQB1 anchors of Graves'/Hashimoto's.
2. Lee HJ, Li CW, Hammerstad SS, Stefan M, Tomer Y. *Immunogenetics of autoimmune thyroid diseases: A comprehensive review.* J Autoimmun 2015;64:82–90. — Korean/Asian-specific HLA-DRB1*04:05, *09:01.
3. Resende de Paiva C, et al. *Association between Hashimoto's thyroiditis and thyroid cancer in 64,628 patients.* Front Oncol 2017;7:53. — PTC-Hashimoto's prognostic association.
4. Boi F, et al. *Thyroid autoimmunity and thyroid cancer: Review focused on cytological studies.* Eur Thyroid J 2017;6:178–186. — autoimmune-PTC connection.
5. Thorsson V, et al. *The Immune Landscape of Cancer.* Immunity 2018;48:812–830 — pan-cancer OptiType HLA + leukocyte estimates.

## Files

- `results/v17_hla_autoimmune/thca_hla_carrier_status.tsv` — per-patient carrier flags
- `results/v17_hla_autoimmune/HA1_carrier_freq_vs_population.tsv`
- `results/v17_hla_autoimmune/HA2_lymph_per_allele.tsv`
- `results/v17_hla_autoimmune/HA3_aiHLA_DM_crosstab.tsv` + `HA3_DM_per_allele.tsv`
- `results/v17_hla_autoimmune/HA4_asian_vs_white_per_allele.tsv`
- `results/v17_hla_autoimmune/HA5_cox_summary.tsv`
- `results/v17_hla_autoimmune/v17_hla_autoimmune_summary.json`
- Figures (dark-bg plotly): `reports/html/figs_interactive/v17/v17_hla_autoimmune_*.html`
- Log: `logs/v17_hla_autoimmune.log`
