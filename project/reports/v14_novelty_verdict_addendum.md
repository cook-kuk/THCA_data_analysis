# v14 Novelty Verdict — Addendum

_2026-04-25. Two additional prior-art references surfaced by the 8-query pipeline run (`notebooks_or_scripts/v14_03_priorart.py`) that are not already cited in `reports/v14_novelty_verdict.md`. Both extend, but do not invalidate, the verdict's "Option A — conservative corroboration" framing._

## A1. Kalfert et al. 2024 — Pathology Research and Practice (PMID 38696857)

**"BRAF mutation, selected miRNAs and genes expression in primary papillary thyroid carcinomas and local lymph node metastases."** Kalfert D, Ludvikova M, Pesta M et al.

60 primary PTC tumours + 40 paired lymph node metastases. Measures BRAF mutation status alongside mRNA expression of four genes (LGALS3, NKX2-1, **TACSTD2**, TPO) and four miRNAs (miR-21, miR-34a, miR-146b, miR-9). Authors' own words: *"to the best of our knowledge, this is the first integrated study of BRAF mutational status, the expression levels of mRNAs of selected genes and miRNAs in primary PTC, and paired LNM."*

**Impact on our novelty claim.** This extends the 2018 Liu (PMID 31949805) + Bychkov (PMID 29228520) IHC line of evidence into the **transcriptomic axis** at tissue level. It is *strictly closer prior art* than either 2018 IHC paper for a tissue-level BRAF × TACSTD2 mRNA claim. It does NOT propose ADC repurposing, does NOT stratify by RAI-refractoriness, and does NOT extend to cell lines.

**Required change to the main verdict and drafts.**
- Move this to **Section B (already-established) row-2 of the verdict**, between the 2018 IHC primaries and the 2022 Pathobiology TMA.
- In `reports/v13/v13_trop2_standalone_draft.md §3 (Literature review)`, add explicit sentence: *"Kalfert et al. (2024) extend this observation to the mRNA axis in 60 primary PTC + 40 paired LNM, establishing BRAF status × TACSTD2 transcript relationships in tissue — the transcriptomic foundation on which our cell-line and pharmacogenomic evidence rest."*
- **Do not** let our cell-line finding (v14 §5.Y.3: top-3 TACSTD2 CCLE lines all V600E) be framed as "first transcriptomic link". Frame it instead as "first extension to cell-line models" — which is accurate.

## A2. Nieto-Jiménez et al. 2023 — Clinical and Translational Medicine (PMID 37740463)

**"Uncovering therapeutic opportunities in the clinical development of antibody-drug conjugates."** Nieto-Jiménez C, Sanvicente A, Díaz-Tejeiro C et al.

Systematic review matching ADC clinical landscape to public genomic TAT datasets. Key sentence directly relevant to us: *"Sacituzumab govitecan (anti-TROP2) in **pancreatic, gastric, thyroid** or endometrial cancer, among others"* — explicitly lists thyroid cancer as a "non-explored niche indication" for sacituzumab govitecan at review level.

**Impact on our novelty claim.** Any language like "we propose sacituzumab govitecan for thyroid cancer" or "we identify thyroid as a TROP2-ADC opportunity" is scooped by this 2023 review. Our actual contribution is narrower and must be stated as such:
- Not "first to propose SG for thyroid" — scooped (Nieto-Jiménez 2023)
- Not "first to identify TROP2 as thyroid ADC target" — scooped (Cancers 2022 PMID 35158847; IJSP 2021 PMID 33289434 in ATC; Nieto-Jiménez 2023)
- **DOES remain** "first to propose **BRAF V600E sub-stratification** of the already-proposed thyroid TROP2 ADC opportunity" — novel, and supported by the trial-design gap (NCT06235216, NCT07521670 both BRAF-unselected).

**Required change to the main verdict and drafts.**
- Add to Section B of verdict as a new row.
- In `reports/v13/v13_trop2_standalone_draft.md` Abstract and §5 (Proposed correlative design), replace any "we propose" / "we identify" TROP2-ADC-for-thyroid language with "we sub-stratify the previously-proposed thyroid TROP2-ADC opportunity by BRAF V600E status."
- Cite Nieto-Jiménez 2023 in the Introduction alongside the 2018 IHC priors.

## Combined implication

Three tiers of priority reduction now stacked:
1. **2018 IHC (Liu, Bychkov)** — BRAF V600E ↔ TROP2 at IHC level in PTC. Priors.
2. **2024 Kalfert** — BRAF × TACSTD2 at mRNA level in PTC + LNM. Closer priors.
3. **2023 Nieto-Jiménez** — SG for thyroid proposed at review level. Closer prior for the therapy-repurposing framing.

Genuinely preserved novelty after all three:

| Our claim | Still novel after A1 + A2? | Why |
|---|---|---|
| Cell-line (CCLE 13-line) transcriptomic corroboration | **Yes** | Kalfert's 60 primary PTC tissues do not cover cell lines. |
| PRISM SN-38-class BRAF-selectivity in thyroid | **Yes** | Neither Kalfert nor Nieto-Jiménez ran pharmacogenomic screens in thyroid. |
| BRAF V600E sub-stratification for ongoing thyroid SG trials (NCT06235216, NCT07521670) | **Yes** | Neither paper proposes this. Both trials remain BRAF-unselected in 2026-04-25 eligibility. |
| Two-gate biomarker (BRAF-like classifier + TACSTD2 z-score) | **Yes** | No prior paper combines transcriptomic classifier with TROP2 expression as a single enrichment signal. |
| DIAL-guided end-to-end pipeline | **Yes** | Methodology-level; no prior paper integrates DIAL → 8-target panel → multi-modality matrix → trial biomarker in thyroid. |

## Recommended action

- Treat the main verdict's **Option A** framing as still correct — nothing in this addendum changes the recommendation.
- Insert Kalfert 2024 + Nieto-Jiménez 2023 into the Introduction paragraph of any paper draft before submission. Both belong in the "previously established — on which we build" paragraph, not in the "we contribute" paragraph.
- Re-check the rest of the prior-art run (`results/v14_ccle/v14_priorart/all_queries_results.tsv`) once more before submission to see if any other MEDIUM-relevance 2023+ paper deserves promotion.

## Trial-design confirmation (unchanged from main verdict, but explicit)

- **NCT06235216** — Recruiting since 2024-09, n=42, sponsor Grupo Espanol de Tumores Neuroendocrinos. Zero BRAF or TROP2 mention in eligibility. Primary completion 2027-12.
- **NCT07521670 (STRAP)** — Not yet recruiting, planned start 2026-05, n=68, sponsor National Cancer Centre Singapore. TROP2 testing **explicitly waived** ("not required for enrollment"), no BRAF criterion. Primary completion 2028-10.

Both trials collect archival tumour tissue for retrospective translational studies, leaving room for a biomarker-only correlative sub-study overlay as proposed in `v13_trop2_standalone_draft.md §5`.

---

_Provenance: 8 PubMed queries (62 records) + CT.gov v2 API (2 trials) executed 2026-04-25 via `v14_03_priorart.py`. Raw output at `results/v14_ccle/v14_priorart/`. Addendum integrates with, does not overwrite, `reports/v14_novelty_verdict.md`._
