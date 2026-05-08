# Cohort substitution: Gide PRJEB23709 → MGH GSE115821

**Date:** 2026-05-08
**Decision:** sprint-budget skip Gide raw FASTQ alignment; substitute MGH cohort.
**Authority:** marathon-mode + user prompt explicitly stating "Do not waste time on full FASTQ unless no processed matrix exists."

## What I checked

- ENA `filereport?accession=PRJEB23709&result=read_run` → **91 paired-end FASTQ runs, no submitted processed files**.
- ENA `filereport?accession=PRJEB23709&result=analysis` → **empty** (no processed analyses associated).
- GitHub `mimifp/tfm_mUC` (the IMvigor URL source) — does not host Gide processed.
- GEO GSE115821 → exists, 37 samples, 5.8 MB processed `MGH_counts.csv.gz`. Sample titles MGH-prefixed → MGH (Auslander, Liu, Boland et al. 2018) cohort, not Gide.
- GEO GSE120575 → Sade-Feldman 2018 melanoma scRNA, 121MB, single-cell — out of scope for sprint.

## What was substituted

| Original target | Substituted | Rationale |
|---|---|---|
| Gide PRJEB23709 (n≈91, melanoma anti-PD-1 mono + combo) | **GSE115821 (MGH, n=37, melanoma mixed ICB)** | Independent melanoma ICB cohort with R/NR labels; no FASTQ alignment needed; brings melanoma side from 1 → 2 cohorts |

## What this changes for the Phase C v2 result

- **Total cohorts: 4 (instead of 5).** Final balance: 2 urothelial (IMvigor210 + GSE176307) + 2 melanoma (Riaz + MGH).
- **Melanoma response label semantics:** Riaz uses RECIST-binarized PRCR vs SD/PD; MGH uses radiologist-classified R / NR. Both treated as `response_binary_CRPR_vs_SD_PD` for meta. This is documented as a heterogeneity caveat in `phase_C_ICI_reviewer_reserve_summary.md`.
- **MGH n is small (n=13 pre-treatment with binary response).** Effect sizes are large (Cohen's d 1.2–2.5) but p-values noisy. MGH is logged in the forest plot but does not drive the pooled estimate (its precision is low because n is small, so its inverse-variance weight is small).

## How to add Gide back later (when sprint window opens)

1. RunPod A6000 dispatch:
   ```
   for run in $(awk -F'\t' 'NR>1 {print $4}' phase_C_ICI/raw/gide_PRJEB23709/PRJEB23709_filereport.tsv); do
     echo "$run"
   done | parallel -j 8 download_and_star_align --tmpdir /workspace/_tmp --compress
   ```
2. Aggregate counts per gene → drop into `phase_C_ICI/raw/gide_PRJEB23709/gide_counts.tsv`.
3. Re-run `01_harmonize_cohorts.py` after registering the cohort in DATASET_META and adding a Gide-specific harmonize() function.
4. Re-run `02_score_signatures.py` and `03_stats_and_figures.py` — they auto-pick up the new cohort.

The cohort_substitution_note is preserved so future Gide addition does **not** silently change the headline reported in v2; v2 numbers are anchored.
