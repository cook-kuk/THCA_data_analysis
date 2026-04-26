# Patent landscape scan — status (NOT completed)

_Date: 2026-04-25 · Author: Seungho Cook._

## Why deferred

All four free public patent search endpoints failed automated access from this run:

| Endpoint | Result |
|---|---|
| Google Patents (`patents.google.com/?q=...`) | SPA / JS-rendered; WebFetch returned empty body |
| Lens.org search list | HTTP 403 (anti-bot) |
| WIPO PatentScope | HTTP 403 (anti-bot) |
| EPO OPS REST (`ops.epo.org/3.2/...`) | HTTP 403, fair-use policy violation (requires registered key) |
| Europe PMC `SRC:PAT` filter | 0 hits across all relevant queries — coverage too thin to use as substitute |

This is an environmental access limitation, not a substantive null result.

## What needs to happen before paper submission

A **manual ~10-minute** patent search by Seungho on each of these (recommended order, copy-paste queries):

### Google Patents (https://patents.google.com)

1. `BRAF TROP2 thyroid antibody-drug conjugate` (filter: 2018–2026, English)
2. `TACSTD2 thyroid sacituzumab stratification`
3. `BRAF V600E papillary thyroid TROP2 biomarker`
4. `sacituzumab govitecan thyroid carcinoma`
5. `TROP2 antibody drug conjugate thyroid cancer`

### Espacenet (https://worldwide.espacenet.com)

1. `BRAF AND TROP2 AND thyroid` (CPC: A61K, A61P)
2. `TACSTD2 AND thyroid AND ("antibody drug conjugate" OR ADC)`

### USPTO Patents Public Search (https://ppubs.uspto.gov)

Same queries, US-only filter.

## What to look for

- **Composition-of-matter** patents on TROP2 ADC variants for thyroid indication (especially Daiichi Sankyo, Gilead/Immunomedics, MSD, AstraZeneca portfolios).
- **Method-of-use** patents claiming BRAF-stratified TROP2 ADC dosing or BRAF-mutant patient selection.
- **Diagnostic** patents claiming TROP2 IHC + BRAF genotype combination as a companion diagnostic.

## Decision rule

- **No relevant hits** → submission proceeds; cite this null result in the IP-disclosure section if required by venue.
- **Any composition-of-matter hit** for thyroid TROP2 ADC by major sponsor → not a scoop on our paper (we make no IP claim), but worth a sentence in the discussion noting active commercial development.
- **Any method-of-use hit** claiming BRAF-stratified TROP2 ADC dosing → reframe v14 as "independent academic corroboration of an industry-staked stratification hypothesis," cite the patent in §6.X.

## Time budget

10 minutes manual. Can be done the morning of submission.

## Why this scan is low priority for the manuscript

We are submitting a **computational corroboration paper**, not an IP filing. Patent landscape is reference material for the Discussion, not a novelty determinant. The journal-literature + preprint scans (Tasks 1–2 + this preprint pass) are the load-bearing scoop checks; the patent pass is hygiene.
