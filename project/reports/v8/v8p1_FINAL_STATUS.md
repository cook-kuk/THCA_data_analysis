# v8.1 — Final Status

> ⚠️ **v5.2 SUPERSEDED NOTICE (2026-04-25 evening).** The v5.1 central
> claim ("THCA exhibits batch_entangled in 4/5 classifiers under LODO,
> DIAL up to 0.494") does **not** survive proper Leave-One-Dataset-Out
> ComBat (fit-on-train-only). All five THCA classifiers move to
> `true_biology` with `auc_post ≥ 0.74` and DIAL = 0.000 (`results/v5p2_fix/`,
> `reports/v5p2/v5p2_critical_assessment.md`). The "submission-ready"
> language below referred to the v5.1 narrative and is **withdrawn**.
> This file is retained for traceability with each subsection now
> annotated as **SURVIVES v5.2** or **FALLS under v5.2**. The
> biomarker-direction work (raw-count DESeq2, full-gene LMM,
> BRS52 3-cohort cross-platform, gold slate, CCLE check) is independent
> of the LODO ComBat protocol and survives unchanged. The DIAL-flip
> stress tests (S1 ComBat-seq, S2 meta-analysis, S3 pathway, S5 quantum,
> S6A FastRNA) tested a finding that no longer exists; they fall as
> v5.1-leak artefacts. See *Repair plan* at the bottom for what to
> publish.

_Generated 2026-04-25. Single-paragraph submission readiness summary;
detail in `v8p1_rigor_summary.md` and `v8_supplementary.md`._

## One-paragraph summary

**v8.1 rigor upgrade: 6 of 7 gaps closed.** Raw GDC STAR-count DESeq2
recapitulates the v8 pseudo-count slate at **95.7 % same-sign log2FC**
concordance (sanity check; pseudo-count was ≈ 5 % conservative, not
distortive). Full 11 710-gene cohort-aware mixed model (Sul-style
LMM with cohort random intercept) gives **88.6 % honest biomarker
survival** (2 458 of 2 773 originals, vs the misleading v8 99.4 %
in-universe figure); 8 of 8 druggable targets retained under both LMM
and DESeq2. The BRS52 expression signature (Chakravarty 2011) validates
on TCGA-THCA at **95.2 % accuracy** (specificity for RAS = 100 %) and,
applied to GSE33630 PTC and GSE29265 PTC, gives 86.0 % and 83.5 %
direction concordance against the v8.1 raw-count DESeq2 log2FC; on the
1 915-gene three-cohort overlap, **88.3 % of BRAF-vs-RAS significant
genes agree on direction in ≥ 2 of 3 microarray cohorts** and 67.3 %
in all three. CCLE thyroid cell lines (n=13 with expression, 4 BRAF-mut
+ 2 RAS-mut) tested but **NOT included as a 4th cohort**: concordance
54.4 % (chance baseline 50 %), BRS accuracy 66.7 %, consistent with
known cell-line drift from primary tumour (Yu *et al.* 2019, *Nat
Commun*). The complete 5 × 3 quantum classifier grid (S5) is filled
post-budget-skip and preserves THCA-only specificity. Only remaining
gap: **prospective cohort validation pending Prof. 유형원 / Bundang
IRB**. **Status: Bioinformatics paper submission-ready** under the
6-section v8 supplementary + v8.1 rigor addendum.

## Numbers at a glance

| Axis                                | v8 (initial)         | v8.1 (final)           |
|-------------------------------------|----------------------|------------------------|
| DESeq2 substrate                    | pseudo-counts (deg)  | **raw GDC STAR**       |
| DESeq2 sig genes                    | 6 283                | **6 605**              |
| Pseudo ↔ raw concordance            | —                    | **95.7 %**             |
| LMM gene universe                   | 3 000 top-variance   | **11 710 shared**      |
| Honest biomarker survival           | (implicit 99.4 %)    | **88.6 %** (2 458/2 773) |
| Druggable retained (LMM+DESeq2)     | 8 / 8                | **8 / 8**              |
| Gold slate (LMM ∩ DESeq2)           | 760 genes            | **1 786 genes**        |
| Cross-platform validation cohorts   | 1 (GSE27155 only)    | **3 GEO** (+CCLE neg)  |
| 3-cohort full-direction consensus   | —                    | **67.3 %** (1 288/1 915) |
| 3-cohort ≥ 2/3 consensus            | —                    | **88.3 %**             |
| BRS52 TCGA validation accuracy      | —                    | **95.2 %**             |
| Quantum 5 × 3 grid                  | 11/15 (4 skipped)    | **15/15** complete     |
| Specificity finding                 | THCA-only flip       | **THCA-only flip preserved** |

## Submission readiness — re-scored under v5.2

Each item below is re-tagged after the v5.2 LODO retraction:

| Item | Original status | Under v5.2 | Why |
|------|-----------------|------------|-----|
| Anti-pseudo-count (Task A raw STAR DESeq2) | ✅ | ✅ **survives** | DE test on harmonized counts; not a LODO claim |
| Anti-3 000-gene-LMM (Task B full-gene LMM) | ✅ | ✅ **survives** | LMM with cohort random intercept; not a LODO claim |
| 3-cohort BRS direction (Task E) | ✅ | ✅ **survives** | Per-gene mean-difference sign agreement; not a LODO claim |
| Druggable retention 8/8 (DE + LMM) | ✅ | ✅ **survives** | Both frameworks agnostic to LODO ComBat |
| Gold slate 1 786 genes (Task D) | ✅ | ✅ **survives** | DE ∩ LMM intersection; both survive |
| CCLE negative-result check (Task F) | ✅ | ✅ **survives** | Direction concordance; not a LODO claim |
| Reviewer-4 paradigm-invariance (S5 quantum) | ✅ | ⚠️ **falls (probably)** | Quantum DIAL was on v5.1-leak ComBat; needs v5.2 rerun |
| Pathway-level interpretation (S3) | ✅ | ⚠️ **empty** | "Pathway aggregation eliminates the flip" — no flip to eliminate |
| ComBat-seq robustness (S1) | ✅ | ⚠️ **falls** | Tested whether v5.1 flip survives ComBat-seq; flip itself doesn't survive |
| Random-effects meta-analysis (S2) | ✅ | ⚠️ **falls** | m_THCA = 1.00 was on v5.1-leak DIAL values |
| FastRNA centering rescue (S6A) | ✅ | ⚠️ **falls** | Rescues the v5.1 flip; nothing to rescue |
| BUHMBOX PC1 KS (S6B) | ✅ | ✅ **survives** | Structural observation, not a LODO claim |
| LMM-biomarker correction (S6C) | ✅ | ✅ **survives** | Same as Task B above |
| Prospective cohort | ⏳ | ⏳ | Still IRB-pending |

## Repair plan — what to publish

1. **Biomarker paper (paper 2 — TROP2 / drug-repurposing standalone)**:
   ready to ship. The 8 druggable targets, the 1 786-gene gold slate,
   the 3-cohort BRS direction validation, and the negative CCLE check
   are all valid under v5.2 and form a coherent biomarker-direction
   contribution. Target a clinical / drug-discovery journal.
2. **Methodology paper (paper 1 — formerly DIAL Bioinformatics)**: pivot
   from *"DIAL detects THCA-specific batch entanglement"* to
   *"covariate-aware ComBat fit on pooled train+test in LODO inflates
   AUC inversion and imitates batch entanglement"* (see
   `reports/v5p2/v5p2_critical_assessment.md` recommendation 1). The
   v5.2 leak-safe LinearComBat fit/transform decomposition is the fix.
   The retracted v5.1 numbers become the worked example, not the
   headline. Honest, publishable, and the v5p1_paper.tex retraction
   header is already in place.
3. **DIAL-flip stress-test claims (v8 S1/S2/S3/S5/S6A)**: do not cite
   in submission. Mark as v5.1-leak artefacts in any new draft.
4. **Direction-validation work (Tasks A, B, D, E, F + S6B + S6C)**: cite
   freely; these are the spine of paper 2 and a Methods-section
   contribution to paper 1.

## Key artefacts (paper-ready)

- `reports/v8/v8_supplementary.md` (37 KB, 6 sections)
- `reports/v8/v8p1_rigor_summary.md` (13 KB, full v8 → v8.1 comparison)
- `reports/v8/v8p1_FINAL_STATUS.md` (this file)
- `reports/html/pages/v8_statgen_supplement.html` (47 KB, all
  numbers reflect v8.1)
- `results/v8p1_rigor/{a..f}_*` — 6 sub-tasks × ~3 files each
  (raw DESeq2, full LMM, expanded LODO attempt, gold slate,
  BRS validation, CCLE check)
- `results/v8p1_rigor/d_gold_slate/gold_slate_LMM_x_rawDESeq2.tsv`
  (1 786-gene primary biomarker list)
- `results/v8p1_rigor/e_brs_validation/three_cohort_per_gene_consensus.tsv`
  (1 915 genes × 0..3 match count across GEO cohorts)
- `reports/html/figs_interactive/v8/fig{1..6}_*.html` (6 Plotly figures)
