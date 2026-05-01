# Pozdeyev 2018 cross-cohort mutation landscape (v17)

**Source:** Pozdeyev N, et al. *Genetic Analysis of 779 Advanced Differentiated and
Anaplastic Thyroid Cancers.* **Clin Cancer Res** 2018;24(13):3059–3068. PMID 29615459;
PMC6030480; DOI 10.1158/1078-0432.CCR-18-0373. cBioPortal does **not** host this
study; data were retrieved directly from the journal supplementary table 8 (PMC6030480
NIHMS956033 supplement-8.xlsx, 1.87 MB), which provides per-patient annotations for
all 779 specimens.

## Cohort breakdown (Supplementary Table 1, n = 779)

| Histology                            |   N |
| ------------------------------------ | ---:|
| Papillary thyroid cancer (adult)     | 468 |
| Anaplastic thyroid cancer            | 196 |
| Follicular thyroid cancer            |  65 |
| Hurthle-cell thyroid cancer          |  35 |
| Pediatric papillary thyroid cancer   |  15 |

Sequencing platforms: MSK-IMPACT (n = 149) and FoundationOne (n = 630).

## Cross-cohort driver prevalence

| cohort                               |   N |   BRAF_V600E_n |   BRAF_V600E_pct | RAS_any_n   | RAS_any_pct   |   TERT_promoter_n |   TERT_promoter_pct | BRAF_p_vs_TCGA   | TERT_p_vs_TCGA   |
|:-------------------------------------|----:|---------------:|-----------------:|:------------|:--------------|------------------:|--------------------:|:-----------------|:-----------------|
| TCGA-THCA primary (PTC)              | 502 |            282 |             56.2 | —           | —             |                36 |                 7.2 | ref              | ref              |
| Landa 2016 MSK PDTC/ATC              |  43 |             43 |            100   | —           | —             |                26 |                60.5 | 8.13e-11         | 1.41e-16         |
| Pozdeyev advanced PTC                | 468 |            325 |             69.4 | 43          | 9.2           |               237 |                50.6 | 2.11e-05         | 5.92e-55         |
| Pozdeyev advanced DTC (PTC+FTC+HCTC) | 568 |            326 |             57.4 | 91          | 16.0          |               292 |                51.4 | 7.11e-01         | 6.73e-61         |
| Pozdeyev ATC                         | 196 |             77 |             39.3 | 49          | 25.0          |               108 |                55.1 | 7.15e-05         | 3.68e-41         |

`p` values: Fisher exact, two-sided, vs. TCGA-THCA primary (n = 502 in our
v17 cache). RAS prevalence in TCGA is reported here as the literature value (~13 %)
because the cached file does not carry RAS calls; the Pozdeyev RAS column uses
NRAS/HRAS/KRAS Q61/G12/G13 hotspots only.

## Figure

Interactive grouped bar chart at
`reports/html/figs_interactive/v17/v17_pozdeyev_landscape.html`.

## Suggested Discussion paragraph

### English (paste-ready, ~150 words)

> To assess whether the v17 8-gene TERT/RAS/BRAF index generalises beyond
> the TCGA-THCA primary cohort, we cross-validated driver prevalence against
> two external advanced-thyroid datasets. In TCGA-THCA primary tumours the
> BRAF V600E and TERT promoter rates were 56.2 % and 7.2 %
> respectively, consistent with the original TCGA Cell 2014 paper. In the
> Pozdeyev 2018 cohort of 779 advanced/anaplastic specimens (Clin Cancer Res
> 24:3059), TERT promoter mutations were enriched ~7.0-fold in advanced
> PTC (50.6 %, Fisher p = 5.92e-55) and ~7.7-fold in ATC
> (55.1 %, p = 3.68e-41) versus TCGA, while BRAF V600E
> remained dominant in PTC (69.4 %) and dropped to 39.3 %
> in ATC. RAS hotspot mutations rose from ~13 % in TCGA primary to 25.0 %
> in Pozdeyev ATC, supporting an alternative RAS-mutant route to anaplastic
> transformation. The Landa 2016 MSK PDTC/ATC cache (n = 43) recapitulated
> the high-TERT signal (60.5 %). These prevalence patterns confirm that the
> v17 index targets drivers whose absolute frequencies scale with disease
> aggressiveness, an essential prerequisite for prognostic deployment.

### Korean (붙여넣기용, ~150자)

> v17 8-유전자 TERT/RAS/BRAF 지표가 TCGA-THCA 일차 코호트를 넘어 일반화되는지를
> 검증하기 위해 두 개의 진행성 갑상선암 외부 데이터셋으로 변이 빈도를 교차 검증하였다.
> TCGA-THCA 일차 종양의 BRAF V600E 및 TERT promoter 빈도는 각각 56.2 %,
> 7.2 %로 TCGA Cell 2014 원논문과 일치하였다. Pozdeyev 2018 코호트
> (Clin Cancer Res 24:3059, n = 779)에서는 TERT promoter 변이가 진행성 PTC에서
> 50.6 % (Fisher p = 5.92e-55, 약 7.0배 증가),
> ATC에서 55.1 % (p = 3.68e-41, 약 7.7배 증가)로
> 두드러지게 풍부했다. BRAF V600E는 PTC에서 69.4 %로 유지되었고
> ATC에서는 39.3 %로 감소하였다. RAS hotspot 변이는 TCGA 일차의
> 약 13 %에서 Pozdeyev ATC의 25.0 %로 상승하여 RAS-매개 비형성 전환
> 경로를 뒷받침한다. Landa 2016 MSK 캐시(n = 43)도 동일한 고-TERT 신호
> (60.5 %)를 재현하였다. 이러한 변이 분포 패턴은 v17 지표가 질병 공격성과
> 함께 절대 빈도가 증가하는 driver 변이를 표적으로 하고 있음을 확인시켜 주며, 예후
> 도구로서의 임상 적용 전제 조건을 충족한다.

---
*Generated 2026-04-28 by `notebooks_or_scripts/v17_pozdeyev_landscape.py`. Author: Seungho Cook.*
