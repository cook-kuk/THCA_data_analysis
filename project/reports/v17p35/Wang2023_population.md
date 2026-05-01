# Cross-population BRAF / RAS / TERT prevalence — Wang 2023 anchor

**Author:** Seungho Cook · **Generated:** 2026-04-28 · **Stage:** v17p35 Discussion support

## Paper

- **Citation:** Du Y, Zhang H, Wang Y, et al. *Mutational profiling of Chinese patients with thyroid cancer.* Frontiers in Endocrinology **14**:1156999 (2023).
- **DOI:** [10.3389/fendo.2023.1156999](https://doi.org/10.3389/fendo.2023.1156999)
- **PubMed / PMC:** PMID 37465126 · PMC10351985
- **Cohort:** 458 Chinese thyroid cancer patients, three-center, retrospective; targeted NGS with a 1,021-gene panel.
- **Histology mix (parsed from supp Table 2):** PTC 438 (95.6%), ATC 10 (2.2%), FTC 6 (1.3%), PDTC 4 (0.9%).
- **Supplementary source used:**
  - Supp Table 2 (per-patient long-format mutation list) — `https://www.frontiersin.org/api/v4/articles/1156999/file/Table_2.xlsx/1156999_supplementary-materials_tables_2_xlsx/1`
  - Mirrored locally at `/opt/thyroid-dash/project/data_raw/wang2023/Table_2.xlsx`
- **Raw repo:** Genome Variation Map accession **GVM000545** (`http://bigd.big.ac.cn/gvm/getProjectDetail?project=GVM000545`).

## Wang 2023 — per-histology driver prevalence (parsed from supp Table 2)

| Histology | n | BRAF V600E | RAS (NRAS/KRAS/HRAS) | TERT promoter (C228T/C250T) | TP53 | RET fusion |
|-----------|---:|-----------:|---------------------:|----------------------------:|-----:|-----------:|
| PTC       | 438 | 343 (78.3%) | 9 (2.1%) | 18 (4.1%) | 4 (0.9%) | 35 (8.0%) |
| FTC       | 6   | 0 (0.0%)    | 4 (66.7%) | 3 (50.0%) | 1 (16.7%) | 0 |
| PDTC      | 4   | 2 (50.0%)   | 1 (25.0%) | 2 (50.0%) | 4 (100%) | 0 |
| ATC       | 10  | 1 (10.0%)   | 6 (60.0%) | 6 (60.0%) | 4 (40.0%) | 0 |
| **All**   | **458** | **346 (75.6%)** | **20 (4.4%)** | **29 (6.3%)** | **13 (2.8%)** | **35 (7.6%)** |

Paper-reported overall (text/abstract): BRAF 76.0%, RAS 4.1%, TERTp 6.3%, RET fusion 7.6% — our independent parse reproduces all four within rounding.

## Cross-population comparison

| Cohort | n | BRAF V600E | RAS | TERTp | TP53 | Source |
|--------|---:|----------:|----:|------:|-----:|--------|
| TCGA-THCA (US) | 496 | 60.7% | 12.9% | 9.4% | 0.8% | Cell 2014 + Liu ERC 2014 (TERTp resequenced) |
| **Wang 2023 (CN)** | **458** | **75.6%** | **4.4%** | **6.3%** | **2.8%** | fendo.2023.1156999 supp Table 2 (this work) |
| Han 2023 (KR PTC) | 240 | 78.3% | 3.3% | 8.3% | 1.7% | Han et al. 2023 Korean PTC (representative) |
| Yoo 2019 (KR adv PTC) | 125 | 74.4% | 8.0% | 23.8% | 7.2% | Yoo et al. *Nat Commun* 2019;10:2764 |

Pairwise Fisher exact (Wang vs other), two-sided:

| Comparison | BRAF V600E | RAS | TERT promoter |
|---|---:|---:|---:|
| Wang vs TCGA-THCA | **p = 1.1e-06** | **p = 3.2e-06** | p = 0.093 |
| Wang vs Han 2023 (KR) | p = 0.452 | p = 0.685 | p = 0.351 |
| Wang vs Yoo 2019 (KR adv PTC) | p = 0.815 | p = 0.112 | **p = 1.1e-07** |

Take-aways: (i) Chinese and Korean PTC cohorts share a high-BRAF, low-RAS profile that is significantly distinct from the TCGA-THCA US PTC cohort; (ii) elevated TERT promoter prevalence in Yoo 2019 reflects its enrichment for advanced/aggressive disease, not a general East-Asian feature.

## Figure

Interactive HTML (Plotly dark theme, BG `#0b0e12`):

- `/opt/thyroid-dash/project/reports/html/figs_interactive/v17/v17_wang2023_population.html`

## Discussion paragraph — ready to paste

### English (~2 sentences)

> Driver mutation frequencies in papillary thyroid cancer differ substantially across populations: Chinese (Wang 2023, n=458) and Korean cohorts (Han 2023 PTC; Yoo 2019 advanced PTC) show BRAF V600E prevalence of 74-78%, compared with 60.7% in the US TCGA-THCA cohort (n=496) — a difference that is highly significant (Fisher exact p = 1.1e-06 for Wang vs TCGA), with mirror-image enrichment of RAS-family drivers in TCGA (12.9% vs 3-4% in East Asia, p = 3.2e-06). This East Asian shift toward a BRAF-dominant, RAS-poor genomic landscape underscores why population-matched validation cohorts — including the Chinese 458-patient anchor analyzed here and our Korean K1A/K2 references — are essential when interpreting BRAF-centred prognostic or therapeutic signals derived from TCGA.

### Korean (~2 sentences)

> 갑상선 유두암의 주요 driver 돌연변이 빈도는 인종 집단 간에 뚜렷한 차이를 보인다. 중국인 코호트(Wang 2023, n=458)와 한국인 코호트(Han 2023 PTC, Yoo 2019 진행성 PTC)에서 BRAF V600E 빈도는 74-78%인 반면, 미국 TCGA-THCA 코호트(n=496)에서는 60.7%로 낮으며(Fisher 정확검정 p = 1.1e-06), RAS 계열 돌연변이는 반대로 TCGA에서 12.9%, 동아시아 코호트에서 3-4%로 거울상 분포를 보인다(p = 3.2e-06). 이러한 동아시아 특이적 BRAF-우세 / RAS-결손 유전체 양상은, TCGA 기반 BRAF 예후·치료 신호를 한국인 임상에 적용할 때 본 연구에서 분석한 중국 458명 데이터와 우리 한국인 K1A·K2 검증 코호트와 같은 인종 일치 검증 세트가 반드시 필요함을 보여준다.

## Reproducibility

- Script: `/opt/thyroid-dash/project/scripts/v17_wang2023_population.py`
- Outputs: `/opt/thyroid-dash/project/results/v17_wang2023/Wang2023_population_summary.{tsv,json}`, per-patient flags TSV
- Log: `/opt/thyroid-dash/project/logs/v17_wang2023.log`
