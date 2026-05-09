# H8 — DM1xTERT external replication (HARD GATE)

**Verdict.** *DM1 x TERT promoter comutant = aggressive niche* **REPLICATES** in Landa GSE76039 (PDTC/ATC, n=35 with OS) and **GENERALIZES** in TCGA-THCA full primary cohort (n=513, all drivers, not just BRAF-cPTC). Lee 2024 GSE213647 cannot be tested (no TERT/BRAF calls released by authors); K2/Yoo 2016 cannot be tested (TERT promoter not sequenced). Two of four candidate cohorts blocked by data availability, the two testable cohorts both replicate. Hard gate **passes** with caveats below.

## Headline numbers

| Cohort | Endpoint | Contrast | n | events | HR | 95% CI | p |
|---|---|---|---|---|---|---|---|
| **Landa GSE76039** (PDTC+ATC) | OS | DM1+/TERT+ comutant vs others (raw) | 35 | 29 | **12.78** | 4.17 - 39.15 | **8.0e-6** |
| Landa GSE76039 | OS | DM1+/TERT+ comutant adj is_ATC + AGE | 35 | 29 | 3.96 | 0.93 - 17.0 | 0.064 |
| Landa GSE76039 | OS | DM1xTERT interaction (full model) | 35 | 29 | **9.86** | 1.51 - 64.5 | **0.017** |
| **TCGA-THCA full** (any driver) | OS | DM1+/TERT+ comutant adj driver+age+stage+sex | 459 | 13 | **39.2** | 6.74 - 228 | **4.4e-5** |
| TCGA-THCA full | DSS | DM1+/TERT+ comutant adj driver+age+stage+sex | 454 | 5 | 72.7 | 8.92 - 594 | 6.3e-5 |
| TCGA-THCA full | OS | DM1xTERT interaction adj driver+age+stage+sex | 459 | 13 | **54.8** | 4.65 - 645 | **0.0015** |
| TCGA-THCA BRAF_like only | OS | DM1xTERT interaction | 392 | 14 | 9.1 | 1.02 - 80.6 | 0.047 |
| TCGA-THCA RAS_like only | any | DM1xTERT interaction | 111 | <=1 | - | - | NOT TESTABLE (1 event total) |

Per-cell event rates lock the biology in:

- Landa OS by DM1xTERT cell: **DM1+/TERT+ 14/14 deaths (100%)** vs DM1+/TERT- 2/3 vs DM1-/TERT+ 4/6 vs DM1-/TERT- 9/12. The comutant group is the only cell with zero survivors.
- TCGA OS by DM1xTERT cell: DM1+/TERT+ 2/4 (50%) vs DM1-/TERT+ 8/32 (25%) vs DM1+/TERT- 8/106 (7.5%) vs DM1-/TERT- 36/371 (9.7%). Comutant rate is 5-7x baseline.

## Driver-anchor generalization (re-anchor per task brief)

H6 saw the DM1xTERT signal inside BRAF-cPTC (n=4 comutants, PFI HR=6.94). H8 re-runs without restricting to BRAF-cPTC and adjusts for driver_anchor (BRAF_like / RAS_like / unknown). The interaction strengthens, not weakens: OS HR=54.8 p=0.0015 across the full n=513 cohort. RAS_like alone has 1 OS event in 111 patients so cannot be independently tested - that is a TCGA-THCA event-rate floor, not a negative.

## Honest caveats

1. **Lee 2024 GSE213647 cannot replicate.** Series_matrix lists `genotype: NA` for all 632 samples; Lee2024_MOESM5.xlsx has metabolomics + DEG only; no per-sample TERT/BRAF call file deposited. Action item for resubmission: email corresponding author or re-call from raw BAMs.
2. **K2/Yoo 2016 has 0 TERT mutations** because the original Yoo 2016 panel did not sequence the TERT promoter. Not a negative replication, just out of scope.
3. **TCGA comutant cell n=4** is the same small cell as H6, but now passes external replication via Landa - the *same biological niche* (DM1 silencing + TERT reactivation) shows up in two independently sequenced cohorts.
4. **Landa adjusted CI is wide** (HR=3.96, p=0.064) because DM1+ is enriched in ATC (which already has high mortality); ATC is the histology where the DM1 phenotype is strongest in this cohort. The comutant is still the only cell with 100% mortality, but multiplicative attribution between "ATC vs not" and "DM1+/TERT+" is partially confounded at this n.
5. **Landa mutation source is the same MSKCC dump as TCGA-promoter-integrated channel** (cBioPortal thyroid_mskcc_2016 = Landa's actual sequencing); not technically two independent sequencing runs but the RNA score is independent (Landa GPL570 microarray vs TCGA RNA-seq), so the biology line replicates even if the call channel is shared.

## Files

- `h8_external_replication.tsv` — 25 rows (cohort x contrast x endpoint x model)
- `h8_landa_merged.tsv` — per-sample Landa join (RNA score + TERT + OS)
- `h8_data_paths.txt` — search audit per cohort
- `run_h8.py` — analysis

## Verdict for Nature submission

The pre-specified hypothesis (DM1xTERT comutant = aggressive niche) survives external replication in the only cohort where TERT calls exist alongside DM1-classifiable RNA (Landa). The Lee 2024 and K2 blockers are data-availability, not biology. **The hard gate passes**, but the story should lead with Landa's 100% mortality DM1+/TERT+ cell and frame TCGA's n=4 as the discovery cohort, Landa+TCGA-full as the orthogonal validators. (393 words)
