# 07 lit_enrich summary

_Generated: 2026-05-02 01:08_

_Output dir: `/home/seungho/personal/THCA_data_analysis/project/manuscript_p2_brief/lit_enrich_2026_05_02`_


Multi-source literature enrichment across 12 free academic-data APIs.

## Tasks

- **verify_references**: {'entries': 17, 'incomplete': 8}
- **competitive_landscape**: {'claims': 8}
- **hla_frequencies**: {'alleles': 5, 'populations': 4}
- **misattribution_check**: {'checked': 11, 'issues': 2}
- **clinical_landscape**: {'queries': 5, 'results_per_query': {'papillary thyroid cancer immunotherapy': 3, 'RET fusion thyroid': 12, 'BRAF V600E thyroid': 16, 'anaplastic thyroid carcinoma': 30, 'thyroid cancer Hashimoto': 4}}
- **full_text_discovery**: {'entries': 17, 'oa_count': 10}

## Sources used

| # | Source | Purpose | Auth |
|---|---|---|---|
| 1 | OpenAlex | Broad academic + abstracts + citations graph | polite-pool email |
| 2 | CrossRef | DOI authority + metadata + references | polite-pool email |
| 3 | PubMed E-utilities | Biomedical authority | polite-pool email |
| 4 | Europe PMC | Full-text search for OA papers | none |
| 5 | Semantic Scholar | Influential citations ranking | none (rate limited) |
| 6 | arXiv | Preprint search | none |
| 7 | bioRxiv (via E-PMC) | 2024-2026 biomedical preprints | none |
| 8 | ICite (NIH) | RCR citation impact | none |
| 9 | AFND | HLA allele frequency (HTML scrape) | none |
| 10 | ClinicalTrials.gov v2 | Trial landscape | none |
| 11 | Unpaywall | OA discovery | polite-pool email |
| 12 | CORE | OA full-text (200M) | free API key (skipped if absent) |

## How to re-run

```bash
project/.venv/bin/python project/notebooks_or_scripts/v17_lit_enrich_run.py
```

Outputs land back in this directory. Edit `claims`, `alleles`, `trial_queries` in the runner to broaden scope.
