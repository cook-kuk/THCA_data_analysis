# 06 — Manual actions requiring a browser or an institutional login

Nothing here was automated or bypassed. Each item is blocked by a login, a paywall or a
CAPTCHA, and the standing instruction is to record rather than circumvent.

## Action 1 — iProX IPX0011848000 (Zhang 2026, n = 113 proteomics)

**Why it cannot be scripted.** The entire iProX API is behind a login session.
`download.iprox.cn/IPX0011848000/` returns HTTP 403 — the path exists, only the listing is
gated.

**What is needed.** A free account, then the project page. Capture the **file listing only**.
Priority order: processed abundance matrix, then sample annotation, then technical/QC
annotation. Raw mass-spec files are not needed.

**Why it is the highest-value manual item.** Zhang 2026 (*Cell Rep Med*, PMID 41794039) is
the richest cohort found in six waves: n = 113 advanced DTC with WES/panel sequencing, mass
spectrometry, and single-cell/spatial data, annotated with **RAI-avid versus refractory
status, biochemical and structural response, and PFS**. Its DNA accessions are GSA HRA011340
and HRA001107; the proteomics is IPX0011848000. We have already used the published
proteomic-subtype counts (CC1 21% → CC2 57% → CC3 87% refractory, P = 1.4 × 10⁻⁷) but not
the underlying matrix.

## Action 2 — Chirra 2026, *Clin Cancer Res* (Caris, n = 1,348)

**Why it cannot be scripted.** Fully paywalled. Unpaywall, Semantic Scholar and institutional
repositories have no copy; the AACR figshare mirror route that worked for Boucai and ERRITI
does not expose this article.

**What is needed.** SNU library AACR subscription. **One sentence from the Methods** answers
the whole question: *is RAI-refractory a per-patient coded variable, or an inclusion
criterion applied to the whole cohort?* If it is an inclusion criterion, the dataset has no
contrast and drops out.

## Action 3 — Zheng 2024 metabolomics, Baidu repository (new this round)

**Why it cannot be scripted.** The raw matrix is stated to sit in a Baidu repository with the
password given in the paper text. That requires a browser session. No bypass was attempted.

**What is needed.** Confirm only whether a **sample-level metabolite matrix** actually
exists there. Do not invest further before that is confirmed.

**Important caveat before spending any time on it.** This is **serum metabolomics, not
tumour transcriptomics** — the eight-gene panel cannot be tested with it under any
circumstance. Its only role would be a metabolic layer for the atlas and possible pathway
convergence with the Liu proteomics. It is the lowest-priority of the three.

**And if you do contact the author**, ask in the same message for the Wang tissue
metabolomics matrix (*The Oncologist* 2024, PMID 38760956). Same institution, same
laboratory, same corresponding author (Zhu Xin) — see `04_WANG_ACAC_SCREENING.md`. It is one
request, not two.

## Not a manual action — closed

- **Liu 2024 medRxiv** — resolved. The author downloaded the bundle personally; it is in
  `external_data/liu_rrptc_2024/` with a SHA-256 manifest.
- **Rodríguez-Lloveras 2025** — opened via CORE, then confirmed a dead lead and locked as an
  exclusion in `audit/rai_integration_20260806/16_SCREENED_EXCLUSIONS.tsv`.
- **DECISION supplement** — obtained through the AACR figshare mirror; no manual step needed.
  What it does *not* contain is documented in `01_DECISION_READINESS.md`.

## Author requests drafted but deliberately not sent

Per the standing instruction (*"이메일 발송은 하지 마라"*), the following exist as drafts only:

| Recipient | Ask | File |
|---|---|---|
| Jaume Capdevila, Vall d'Hebron | DECISION normalised transcript matrix, sample-ID-to-Supplementary-Table-1 mapping, per-patient PFS time and event | `docs/decision_minimum_data_request.md` |
| Four earlier recipients (Boucai/Fagin, Jordà/Robledo, Weber, Muzza) | cohort-specific asks, figures updated to current numbers | `docs/data_requests_drafts_2026_08_06.md` |

Romana Netea-Maier (Radboud UMC) is the newest and cheapest outreach target identified — a
verified address, an explicit request-based data-availability statement, and a real precedent
of depositing an RAI-sensitive-versus-refractory RNA-seq series (GSE112202). No draft has
been written for her yet.
