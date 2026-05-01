---
title: "references.bib audit (manuscript_v8)"
date: 2026-05-01
purpose: "manuscript_v8/03_intro_references.bib 의 incomplete entries 보고. Voice-protected 영역이므로 content 변경 X — 본인 키보드 처리 의제."
---

# references.bib audit — `manuscript_v8/03_intro_references.bib`

⚠️ **본 audit 은 report only**. `manuscript_v8/` 는 voice-protected 영역 (Hook/Aim/Disc 3.1/Limitations/Cover Para 1/Q9 + bib entry content). Entry content 자체는 본인 키보드 처리.

## 1. 전체 통계

| 항목 | 값 |
|---|---|
| Total entries | 17 |
| Incomplete (missing doi / journal=TBD / title placeholder) | **6** |
| Chu 2018 (Pillar 1 backbone) status | ✅ found, but incomplete |

## 2. Incomplete entries (6 / 17)

| Key | Year | Type | Issues |
|---|---|---|---|
| **Chu2018** ⭐ | 2018 | article | doi missing, journal=TBD, title placeholder |
| Liu2017 | 2017 | article | doi missing, journal=TBD, title placeholder |
| Wang2024 | 2024 | article | doi missing, title placeholder |
| SEER_PTC | 2024 | article | doi missing, journal=TBD |
| Lu2023 | 2023 | article | doi missing, title placeholder |
| Bradley2010 | 2010 | article | doi missing, journal=TBD |

## 3. Chu 2018 entry — Pillar 1 backbone 정확 정보

본 brief 에서 verified 정보 (PMC 6161647 web fetch 2026-05-03):

```bibtex
@article{Chu2018,
  author    = {Chu, X. and Pan, C-M. and Zhao, S-X. and others},
  title     = {Han Chinese Graves' disease HLA fine-map (4-digit SNP2HLA imputation against Pan-Asian reference panel)},
  journal   = {Journal of Medical Genetics},
  volume    = {55},
  number    = {10},
  pages     = {685--692},
  year      = {2018},
  doi       = {10.1136/jmedgenet-2017-105146},
  pmcid     = {PMC6161647},
  note      = {Paper 2 Pillar 1 backbone. n=1,468 GD vs 1,490 ctrl Han Chinese. DPB1*05:01 OR=1.90 [1.69, 2.14] published.}
}
```

→ Manuscript_v8 voice-protected 이므로 본인이 위 entry 로 갱신.

## 4. 다른 incomplete entries — 본인 verify 자료

본 audit 은 manuscript_v8 prose 영역이라 content 채우지 않음. 본인 verify 시 다음 사항 도움될 수 있음:

### 4.1 Wang 2024 (Shanghai n=2,844)
- 분야: PTC mutation landscape, Asian population
- 메모리에서 referenced: "Wang 2024 Shanghai n=2,844 (mutation only)"
- 본인 verify 필요

### 4.2 Liu 2017 (Asian PTC n=583)
- 메모리에서 referenced: "Liu 2017 Asian n=583"
- 본인 verify 필요

### 4.3 SEER_PTC 2024
- SEER database PTC incidence/mortality 통계
- 본인 verify 필요

### 4.4 Lu 2023 (single-cell)
- 메모리: "Lu 2023 GSE193581 (n=23 samples)"
- GEO accession + 정확 cite 본인 verify

### 4.5 Bradley 2010
- 본인 verify 필요 (context: PTC 진단 history?)

## 5. 정상 entries (11 / 17)

DOI / journal / title 모두 채워진 entries 11개는 audit 통과. 이들은 voice-protected 이지만 이미 정확 — 본인이 추가 작업 불필요.

## 6. References for re-verification

| Entry | Suggested verification source |
|---|---|
| Chu2018 | PMC 6161647 (verified 2026-05-03 web fetch) |
| Wang2024 | PubMed search "Wang 2024 thyroid cancer Shanghai" |
| Liu2017 | PubMed search "Liu 2017 Asian thyroid 583" |
| SEER_PTC | https://seer.cancer.gov/ |
| Lu2023 | GEO GSE193581 + linked publication |
| Bradley2010 | 본인 origin 자료 확인 |

---

*Generated 2026-05-01. Audit only — content 갱신은 본인 키보드. references.bib 는 manuscript_v8 voice-protected 영역.*
