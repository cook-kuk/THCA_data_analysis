# Task 1 — BRAF + TROP2 thyroid prior-art assessment

_Date: 2026-04-25 · Author: Seungho Cook · Source: PubMed E-utilities, queries A–H, 65 unique records, all top hits manually reviewed._

## What is genuinely novel in our v13/v14 work

1. **Cell-line-level (CCLE) corroboration of the BRAF-axis TROP2 up-regulation hypothesis.**
   The 2018 IHC reports (PMID 31949805, PMID 29228520) established BRAF V600E ↔ TROP2 high in primary PTC at the protein level. Our v14 result (top-3 TACSTD2 lines = all V600E in 13-line CCLE panel) is the first **cell-line transcriptomic** restatement, and matters because it lets the link travel into in-vitro pharmacology (PRISM, drug screens, xenograft). It is not a new discovery — it is an **orthogonal substrate** for an existing claim.

2. **PRISM-level SN-38-class selectivity for BRAF-like thyroid lines.**
   No prior PubMed record (queries G/H/T2c/T2d, 33 records reviewed) reports topotecan/irinotecan ΔLFC stratified by BRAF status in thyroid cell lines. The closest hit, PMID 34560251 (2021, *Pharmacological Research*), describes TOP1 modulation during melanoma BRAF-i resistance — this is a **resistance** narrative in melanoma, not a **primary-sensitivity** narrative in thyroid. The thyroid-specific PRISM finding is — to our reading — a first specific report.

3. **Two-gate biomarker proposal (BRAF-like × TACSTD2-high) for sacituzumab govitecan stratification.**
   Both currently active trials (NCT06235216, NCT07521670) enrol unselected for BRAF status (and NCT07521670 explicitly waives TROP2 IHC). The *proposal to stratify* is therefore actionable; the underlying biology is older.

## What is incremental but supported by literature

- **TROP2 as an ADC target in thyroid cancer.** PMID 35158847 (2022, *Cancers*) and PMID 33289434 (2021, *Int J Surg Pathol*) both name TROP-2 as a viable ADC target in **anaplastic** thyroid carcinoma; PMID 35477165 (2022, *Pathobiology* TMA on 18,563 tumours) maps TROP2 expression pan-tumour and includes thyroid. Our PTC focus is incremental relative to the established ATC framing.
- **CCLE-thyroid → BRAF subtype mapping.** Standard reference; v8 BRS52 cell-line failure (66.7 %) is acknowledged in v14_impact_summary.md.

## What was already known

| Finding | Established by | Year |
|---|---|---|
| BRAF V600E ↔ TROP2 over-expression in PTC (IHC) | PMID 31949805 (Liu et al., *Int J Clin Exp Pathol*) | 2018 |
| TROP2 IHC predicts BRAF mutation in PTC | PMID 29228520 (Bychkov et al., *J Pathol Transl Med*) | 2018 |
| TROP2 IHC in follicular-patterned thyroid neoplasms | PMID 35121239 (*Ann Diagn Pathol*) | 2022 |
| TROP2 / Nectin-4 / GPNMB / B7-H3 as ADC targets in ATC | PMID 35158847 (*Cancers*) | 2022 |
| TROP-2 / 5hmC / IDH1 in ATC | PMID 33289434 (*Int J Surg Pathol*) | 2021 |
| Pan-tumour TROP2 IHC landscape (18,563 tumours) including thyroid | PMID 35477165 (*Pathobiology*) | 2022 |
| Sacituzumab govitecan clinical activity in TNBC, urothelial | PMID 39067901 (*Lancet*, TROPiCS-02) | 2024 |
| BRAF-mutant mCRC standard chemo backbones include irinotecan (FOLFOXIRI, FOLFIRINOX, encorafenib + cetuximab vs FOLFIRI) | PMID 41852303, 41505697 (multiple) | 2024–2026 |

## Recommendation: framing for paper

The accurate framing is **orthogonal corroboration plus stratification proposal**, not de novo discovery. Concretely:

- Open §5.Y (CCLE/PRISM) by **citing the 2018 IHC papers** explicitly: "Liu 2018 and Bychkov 2018 reported BRAF V600E ↔ TROP2 up-regulation in primary PTC by IHC; here we ask whether the relationship persists at the cell-line transcriptomic level (CCLE) and whether it travels into pharmacological vulnerability (PRISM)."
- Frame the SN-38 / topotecan / irinotecan finding as **first thyroid-specific PRISM read-out**, while flagging the melanoma TOP1-resistance literature (PMID 34560251) and the established BRAF-mCRC irinotecan exposure as related but distinct.
- For the stratification claim: anchor it to the trial-design facts (NCT06235216/NCT07521670 are BRAF-unselected) rather than to the biomarker novelty.
- Do **not** claim "first report that BRAF drives TROP2 expression" — that claim was made in 2018.

## Caveats on this query

- PubMed coverage only; preprints (bioRxiv, medRxiv), conference abstracts (ASCO, ESMO, AACR), and patent filings not searched.
- E-utilities queries used Boolean term match; broader literature on "thyroid TROP2 IHC + BRAF correlation" may exist under non-MeSH wording.
- One auto-flag: the 2026 *Cancer Medicine* paper (PMID 41852303, encorafenib + cetuximab vs FOLFIRI in BRAF V600E mCRC) again exposes BRAF-mutant patients to topo-I poison — keep this on the radar in §5.Y discussion.
