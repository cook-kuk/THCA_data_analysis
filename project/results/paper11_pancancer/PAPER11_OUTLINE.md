# Paper 11 — Pan-cancer molecular dark matter axis

**Target:** Nat Commun.
**Status as of 2026-05-08:** core analyses done (Phase A-E framework ran). Manuscript writing not started (marathon-mode block).
**Voice-protected sections** (Hook / Aim / Discussion 3.1 / Limitations / Cover Para 1 / Q9): author-keyboard only.

## 1-line take-away

The thyroid-derived DM1 dediff axis (8-gene panel + 21-gene Module 4 architecture) is reproducible across **32 cancer lineages × 10,978 TCGA samples** when re-anchored on each cancer's own lineage TFs, prognostic in **10 lineages**, mechanism-anchored by **inflamed + dediff + OXPHOS-low + MAPK-high** Hallmarks (IL-6/JAK/STAT3 r=0.75 in 32/32 cancers), and points to **MYC + NAMPT** essentiality and **MAPK pathway inhibitors** (RAF/ERK) as DM1-high-selective therapeutics validated in 1,141 DepMap CRISPR + 1,518 PRISM drug screen.

## Headline numbers (already in hand)

| Layer | Finding | Evidence |
|---|---|---|
| Pan-cancer scope | 11,069 TCGA samples × 33 lineages, DM1 score per sample | `pancan_dm1_scored.tsv` |
| **Phase B Prognostic** | **10 lineages OS/PFI Cox FDR<0.1** (LGG HR=36, LUAD HR=11, GBM HR=15…) | `phase_B_survival/cox_per_lineage.tsv` |
| **Phase D Druggable** | **MYC d=−0.50 (p=1e-11), NAMPT d=−0.45 (p=1.5e-9)** essentiality in DM1-high lineages (n=1,141 cell lines, DepMap) | `phase_D_depmap/dm1_high_low_dependency.tsv` |
| **Phase E v2 Portable** | Lineage-portable DM1 reproduces in **32 lineages × 10,978 samples**; **10 sig prognostic OS p<0.05** (SKCM HR=0.78 p=4e-7 protective, UVM HR=1.62 p=4e-4 risk, LGG HR=1.19 p=2e-3 risk, KIRP HR=1.31 p=2e-3 risk, PAAD HR=1.18 p=5e-3 risk, LUAD HR=1.11 p=0.05 risk, +4 more); median r=0.33 with thyroid-anchored DM1 | `phase_E_lineage_specific/cox_per_cohort_portable.tsv` |
| **Phase A Epigenetic mechanism** | Pan-cancer epigenetic-machinery RNA expression (DNMT/EZH2/HDAC/PRC2) × DM1 axis: **19 lineages FDR<0.1 (DM1_like)**, 17 (DM1_portable). Top: UVM r=0.58, LUAD r=0.47, MESO r=0.38, BRCA r=0.36, **THCA r=0.30** (RNA-level mirror of TCGA HM450 TPO d=2.30 finding) | `phase_A_epigenetic/epi_index_dm1_corr_per_lineage.tsv` |
| Paper 12 network | 21-gene Module 4 = lineage_TF + JAK/STAT + SFK + Epigenetic + Metabolism + DDR; 137 top edges | `run.log` (P12 summary) |
| Phase C ICI v1 (Hugo+Riaz) | Hugo n=27 d=−0.21, Riaz n=27 d=+0.64 — under-powered, mixed (kept as honest secondary finding, not headline) | `phase_C_ICI/dm1_inflam_response_per_cohort.tsv` |
| **Phase C ICI v2** (4 cohorts, n=421) | **IFNG OR=1.53 p=3.7e-4 / HLA-I OR=1.35 p=0.014 / checkpoint OR=1.32 p=0.020 / cytolytic OR=1.27 p=0.046** (4 GREEN); DM1_inflam composite + lineage-portable DM1 sign-consistent in **4/4 cohorts**; OS Cox 8/9 signatures HR<1 in both urothelial cohorts (median HR=0.71-0.86) | `phase_C_ICI/results/tables/final_phase_C_ICI_summary_table.tsv` + `phase_C_ICI_reviewer_reserve_summary.md` |

## Pre-existing connected findings to cite/re-use

- TCGA-THCA HM450 methylation × DM1: TPO Δβ=+0.42, d=+2.30, p=1.9e-18 (`audit_2026_04_30/round5/r5_2_*`).
- Mun 2025 protein/phospho × Mun PTC→ATC dediff: thyroid_diff d=−1.91, myeloid d=+1.52 (sign-match 7/7 with RNA Track B-lite).
- CRM 2026 advanced DTC × CC1/CC2/CC3 × RAI: CC3 87 % refractory; 8/8 of our 8-gene ⊂ TDS (`processed/cellrepmed2026_table_S1_clinical.tsv`).
- Track B-lite RNA: thyroid_diff d=−2.51 / myeloid d=+2.53 ATC vs other (n=415, 4 cohorts).

## Proposed figure plan (6-figure Nat Commun structure)

- **Fig 1.** Discovery — TCGA-THCA DM1 sub-cluster (existing Paper 1 Fig 1) + 21-gene Module 4 network. Forward-cite Paper 1.
- **Fig 2.** Pan-cancer scope — DM1 score distribution per lineage (`F11_01_lineage_dm1_box.png`) + lineage candidate-gene heatmap (`F11_04_lineage_candidate_heatmap.png`).
- **Fig 3.** Lineage-portable DM1 — TF panel per cancer + architecture conservation; correlation with thyroid-anchored DM1 across 12 cancers (Phase E).
- **Fig 4.** Prognostic — Cox forest plot OS/PFI across 33 lineages (Phase B `forest_plot_OS.png`) + KM for top 4. Validate in Phase E lineage-portable DM1 (4 sig).
- **Fig 5.** Therapeutic dependency — (a) DM1-high vs DM1-low CRISPR essentiality (Phase D), MYC/NAMPT/JAK2 highlighted; (b) PRISM drug screen volcano plot showing 7/11 FDR<0.05 hits are MAPK inhibitors (Phase H); pre-clinical + clinical-stage actionability.
- **Fig 6.** Hallmark biology + cross-modality triangulation — (a) Phase G Hallmark heatmap: 7 immune Hallmarks 32/32 sig + OXPHOS-low; (b) TCGA HM450 (TPO d=+2.30) + Mun protein (thyroid_diff d=−1.91) + CRM CC3 87% RAI-refractory + GSE151179 RAI-dediff axis. Brings the 4 independent thyroid cohorts together as the anchor for the pan-cancer claim.

## Scope boundaries (must respect)

- **Track B (Paper 3) FROZEN** — Phase C ICI is permitted only because lite already ran on disk; no new Track B unfreeze.
- **HLA / arcasHLA isolated** — Paper 4 backlog, not joined here.
- **Korean K2 prospective NOT included** — Bundang outreach pending; this paper is open-data-only.
- **MTBLS3339 metabolomics** — currently blocked from programmatic fetch; defer or omit.
- Voice-protected sections: author keyboard only.

## Missing for full Nat Commun reach

- ~~Pancancer methylation HM450 × DM1~~ — **Phase A done as RNA-level epigenetic-machinery proxy**; full HM450 fetch was blocked (Xena/GDC/PanCanAtlas all 403/404). Mechanism layer confirmed via DNMT/EZH2/HDAC expression × DM1 in 19/33 cancers.
- ~~Lineage-specific DM1 in 21 more cancers~~ — **Phase E v2 done**, all 32 lineages.
- ~~IMvigor210 + Van Allen + Gide ICI cohort fetch~~ — **Phase C v2 done 2026-05-08** (4 cohorts, n=582 total / 421 with response). Corrected URL: `packageVersions/` was the missing path. Substituted MGH GSE115821 for Gide (FASTQ-only on ENA, sprint-class). Van Allen kept deferred (dbGaP-controlled).
- Per-sample CRM protein abundance (xlsx Data S1-S7, user browser download still needed).
- Mechanistic experimental validation (CRISPR/orthotopic; out of marathon scope).

## Dependencies / file map

```
paper11_pancancer/
├── pancan_dm1_scored.tsv
├── pancan_lineage_stats.tsv
├── lineage_candidate_corr.tsv
├── PAPER11_OUTLINE.md   ← this file
├── figures/F11_01..06.png
├── phase_B_survival/    ← Phase B (Cox forest)
├── phase_C_ICI/         ← Phase C v1 (Hugo+Riaz only) + Phase C v2 (4 cohorts: IMvigor210, GSE176307, riaz_GSE91061, MGH_GSE115821)
│   ├── results/phase_C_ICI_reviewer_reserve_summary.md
│   ├── results/figures/forest_logOR_response.{pdf,png}
│   └── scripts/{00_inventory,01_harmonize_cohorts,02_score_signatures,03_stats_and_figures}.py
├── phase_D_depmap/      ← Phase D (MYC/NAMPT druggable)
├── phase_E_lineage_specific/ ← Phase E (portable axis)
└── scripts/B_*, C_*, D_*, E_*.py
```
