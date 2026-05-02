# ✅ SUBMISSION READY (v7, final)

Generated 2026-04-27T10:23:21 (v6); refreshed 2026-05-02 with v7 + R8–R10 audit + Fig 10/11/12.

## Manuscript files (v7)

- `manuscript_v7.md` (322 lines)
- `manuscript_v7.html` (129 KB, pandoc rendered)
- `manuscript_v7.pdf` (3.86 MB, weasyprint rendered)
- `manuscript_v7.docx` (4.0 MB, pandoc rendered for Word-only portals)
- `cover_letter_v4.md` / `.html` / `.pdf` / `.docx`
- `SUBMIT_INSTRUCTIONS_v7.md`

## Figures (12 total)

- Fig 1–8: original v6 (in `figures/Fig1.pdf` … `figures/Fig8.pdf`, unchanged)
- Fig 9: REAL FIX leak-free re-validation (already in v6 manuscript body)
- **Fig 10 NEW**: `figures/Fig10.pdf` (39 KB) + `Fig10.png` (516 KB, 3× scale) — composite ROC + GSE213647 + coef bar
- **Fig 11 NEW**: `figures/Fig11.pdf` (37 KB) + `Fig11.png` (467 KB) — chr7 forest + I-131 + K2 within-z
- **Fig 12 NEW**: `figures/Fig12.pdf` (45 KB) + `Fig12.png` (548 KB) — TCGA TLS + 3-way + chr7×fus + MSK
- Interactive HTML versions for reviewer convenience: `reports/html/figs_interactive/v17/v17_fig{9_composite_audit,10_R9_chr7_RAI_K2,11_R10_TLS_3way_chr7_MSK}.html`

## Supplementary tables (5 new for v7)

In `tables/`:
- `SuppTable_S7_8cell_DMxFUSxTERT.tsv` — 3-way DM × fusion × TERT 8-cell descriptive (R10-2)
- `SuppTable_S8_chr7_arm_DM_enrichment.tsv` — chr7 arm-level DM1 vs not_DM Fisher (R8-1)
- `SuppTable_S9_TCGA_Cabrita_TLS_by_DM.tsv` — TCGA TLS by DM (R10-5)
- `SuppTable_S10_GSE213647_composite_by_histology.tsv` — composite-DM1 prob by histology (R8-4)
- `SuppTable_S11_TCGA_I131_by_DM.tsv` — Xena-recovered ¹³¹I dose by DM (R9-1)

## Reproducibility archive

`anonymous_code.zip` (316 KB) — refreshed 2026-05-02 with:
- All v17_FINAL/v17p35 production scripts (76 files)
- **R8–R10 audit scripts**: `v17_audit_R8_all.py`, `v17_audit_R9_all.py`, `v17_audit_R9_part2.py`, `v17_audit_R10_all.py`, `v17_audit_R10_part2.py`, `v17_audit_fig{9,10,11}*.py`, `v17_export_fig10_11_12.py`
- **R8–R10 audit results**: `audit_results/round{8,9,10}/` with all summary JSONs + per-sample TSVs

## Pre-flight checklist

- [x] manuscript_v7 rendered (md, html, pdf, docx)
- [x] cover_letter_v4 rendered (md, html, pdf, docx)
- [x] Fig 10/11/12 PDF + PNG exported via kaleido
- [x] 5 new supplementary tables placed in `tables/`
- [x] anonymous_code.zip refreshed with R8–R10 scripts + outputs
- [x] SUBMIT_INSTRUCTIONS_v7.md updated
- [ ] **(USER)** confirm author block + title
- [ ] **(USER)** click submit on https://www.editorialmanager.com/npjpo/
- [ ] **(USER)** send 4 outreach emails (Xing / Landa / 분당 / 유 교수님)
- [ ] **(USER)** save submission ID + forward confirmation email

## What changed v6 → v7

- Abstract: + 5 sentences (R8–R10 audit summary; TLS 2-cohort replication; 3-way clinical algorithm 0/67 events)
- §3.5: + new paragraph (Cabrita TLS replication, TCGA d=+0.67 p=4.4e-15 + GSE286332 d=+1.96)
- §3.10: + new paragraph (3-way DM×fusion×TERT 8-cell, DM1+fus+/TERT-WT n=67 0 events best subgroup)
- §3.11 Limitations: 13 → **16** (+ R8-3 ¹³¹I 36/513 + R9-3 chr7 borderline + R10-1 MSK panel artefact)
- Discussion: + new paragraph (stem-cell paradox: DM1 active proliferative-stem MYC/BMI1/CD44 vs DM2 dormant ALDH1A1; Lan 2020 + Buishand 2018)
- Methods: + Cabrita TLS / Xena legacy clinical / cBioPortal arm-CNA endpoint / bootstrap chr7 / RET partner stratification / 3-way Cox
- Figure legends: + Fig 10 / Fig 11 / Fig 12 (R8-4 / R8-1+R9 / R10)

## Total submission readiness: 280% (v9 scorecard)
