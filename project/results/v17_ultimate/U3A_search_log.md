# U3A — Search log

_Date: 2026-04-27_

| # | Source | URL / query | Outcome |
|---|--------|-------------|---------|
| 1 | PubMed Mu 2024 | `WebSearch query 'Mu Z 2024 thyroid radioiodine avidity JCEM 109 1231'` | FOUND PMID 38060243 |
| 2 | Oxford Academic JCEM | `https://academic.oup.com/jcem/article/109/5/1231/7460630` | Found article landing page (paywalled body) |
| 3 | PMC full text | `https://pmc.ncbi.nlm.nih.gov/articles/PMC11031230/` | Open-access PMC version available; Data Availability extracted |
| 4 | PubMed | `https://pubmed.ncbi.nlm.nih.gov/38060243/` | Confirmed PMID + abstract |
| 5 | Ovid mirror | `https://www.ovid.com/journals/jceme/fulltext/10.1210/clinem/dgad697...` | Paywalled |
| 6 | NGDC GSA-Human HRA004166 | `https://ngdc.cncb.ac.cn/gsa-human/browse/HRA004166` | CONTROLLED ACCESS; DAC contact: linys@pumch.cn |
| 7 | GEO GSE151181 | `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151181` | OPEN; n=99 SuperSeries |
| 8 | GEO GSE151179 | `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151179` | OPEN; n=52 mRNA |
| 9 | GEO GSE151180 | `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151180` | OPEN; n=47 miRNA |
| 10 | GEO GSE190966 | `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190966` | OPEN but bovine — exclude |
| 11 | GEO GDS landscape | `https://www.ncbi.nlm.nih.gov/gds/?term=thyroid+AND+radioiodine+AND+expression` | 15 thyroid+RAI series enumerated |
| 12 | Sci-Hub / preprint | `(skipped per ethics)` | NOT ATTEMPTED |

## Notes

- All lookups executed via the agent's WebFetch / WebSearch tools on 2026-04-27.
- Sci-Hub / preprint paywall-bypass routes were intentionally skipped.
- HRA004166 (Mu 2024 NGS) is **controlled access** — DAC application required; no scraping.
- GEO landscape query returned 15 thyroid+RAI series; the GSE151179/180/181 trio is the strongest open RAI-avid vs RAI-refractory ground truth (n=99, PMID 33198784).
