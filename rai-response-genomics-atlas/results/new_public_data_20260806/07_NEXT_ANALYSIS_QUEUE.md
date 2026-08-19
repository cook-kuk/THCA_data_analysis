# 07 — Next analysis queue

Ordered by value per unit of effort, given that no new patient-level molecular data entered
the project this round.

## Queue

| # | Task | Blocked by | Effort | Value |
|---|---|---|---|---|
| 1 | **Draft the Capdevila request and hold it** for author approval | nothing — done, `docs/decision_minimum_data_request.md` | done | the only route to a 57-patient untreated comparator |
| 2 | **Position the redifferentiation pooled result in the manuscript** — one Discussion sentence in Paper 1, one figure in the companion | author decision on Paper 1 vs companion split | low | draws the exact boundary of the driver-axis null |
| 3 | **Fold the endpoint hierarchy into Methods and Limitations** | nothing — `results/rai_endpoint_hierarchy.tsv` exists | low | pre-empts the "why is your endpoint so indirect" review |
| 4 | **Add the six-portal single-cell negative to Limitations** with names and counts | nothing | low | converts an absence into a stated, checkable search |
| 5 | **Resolve Figure 6D** — remove or move to supplement | author decision | low | P0 blocking item; a reviewer reproduces the confound in ten minutes |
| 6 | Ask Netea-Maier (Radboud) for the GSE112202-successor cohort | author approval to send | low | cheapest outreach with a deposit precedent |
| 7 | Confirm the Barcelona/Madrid methylation cohort actually carries RAI-avidity annotation before requesting it | author email or a full-text read of Rodríguez-Lloveras 2025 | low | 127 samples with 17 paired metastases; the richest untapped European asset |
| 8 | iProX file listing for Zhang 2026 | **manual browser action 1** | low once unblocked | would let us test the panel in the one cohort with full RAI annotation |
| 9 | Chirra 2026 Methods sentence | **manual browser action 2** | trivial once unblocked | decides whether an n = 1,348 cohort has a contrast at all |
| 10 | DECISION raw-read streaming quantification (download → pseudoalign → delete, per sample) | **must not start** until per-patient PFS is obtained | 25–40 CPU-hours, ~8 GB transient disk | without endpoints it answers nothing |
| 11 | Zheng Baidu matrix | **manual browser action 3** | low | metabolomic layer only; cannot test the panel |

## What is explicitly not in the queue, and why

**Do not start the DECISION raw download.** It is technically feasible — 125 samples,
~7.8 GB each, streaming so only ~8 GB of transient disk is ever needed, roughly 25–40
CPU-hours of pseudoalignment. But every patient in DECISION was already refractory at entry,
so it cannot address who becomes refractory. And the per-patient PFS that would make the
placebo arm useful was **destroyed by spreadsheet formatting** in the published supplement
(nine-digit subject IDs rendered as `1E+08`). Downloading a terabyte to compare score
distributions across arms would answer nothing. This stays blocked on item 1.

**Do not pool Zheng and Wang.** Same institution, same laboratory, same corresponding author,
same refractory definition, overlapping period. Treat as one group's programme
(`04_WANG_ACAC_SCREENING.md`).

**Do not attempt an eight-gene analysis in ERRITI.** There is no molecular matrix — only a
panel gene list and eligibility criteria. Its role is the endpoint-hierarchy reference
(`03_ERRITI_ENDPOINT_AUDIT.md`).

**Do not re-run the public-repository sweep.** Six waves have now covered GEO, SRA/ENA,
ArrayExpress/BioStudies, EGA, dbGaP, GSA, iProX, Zenodo, Figshare, GitHub, nine single-cell
and spatial portals, ClinicalTrials.gov structured results, five national biobanks, the
NIS-theranostics literature outside thyroid, and ~300 of the 600 filtered TCGA-THCA citers.
The remaining ~300 citers are the only unexhausted seam, and page 3 onward was already almost
entirely reviews and protocols. Further sweeping has a low expected yield; the marginal
return has moved to author contact.

## Standing constraints honoured this round

No manuscript file was modified. No commit or push was made. No email was sent. No raw FASTQ
was downloaded. No login, CAPTCHA or paywall was circumvented. Total download 7.20 MB against
a 5 GB soft budget. `/data/rai_atlas` was not used for new work and its existing contents were
neither moved nor deleted.
