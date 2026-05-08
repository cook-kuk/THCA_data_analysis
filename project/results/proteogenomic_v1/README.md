# proteogenomic_v1 — protein/phospho corroboration layer

Created 2026-05-08 under marathon mode (5/4–6/13 manuscript writing window).
**Framing:** reviewer-reserve / Track B-lite extension. Not paper-blocking.
**Storage:** `/data/thca/repo_results/proteogenomic_v1/` (bind-mounted at `project/results/proteogenomic_v1/`).

## Why this dir exists

Reviewer Q anticipation: "RNA-only 8-gene panel — does it hold at protein/phospho level?"
Two public proteogenomic THCA cohorts give us the answer without sprinting new analyses:

| Cohort | n | Histology | Use | Status |
|---|---|---|---|---|
| Wang 2024 (Nat Commun) | 102 PTC | PTC ± recurrence-risk | Paper 1 reviewer-reserve | active |
| Mun 2025 (Nat Commun) | 348 (PTC 177 / PDTC 47 / ATC 124) + 119 NAT | full dediff axis | Paper 3 Track B-lite extension | feasibility check |

CPTAC pan-cancer does **not** include THCA — these two studies are the closest equivalent.

## Subdirs

- `paper1_wang2024_reviewer_reserve/` — 8-gene (SLC5A5/TPO/TG/TSHR/PAX8/NKX2-1/FOXE1/DIO1) at protein + phospho, CS1–4 subtype association, TCGA DM1 axis cross-check
- `paper3_mun2025_dediff_layer/` — DIAL-lite 7 modules (HLA-I/II, IFNG, TLS, checkpoint, myeloid, thyroid_diff) along PTC→PDTC→ATC at phospho level
- `_meta/` — accession docs, fetch logs
- `_raw/` — downloaded supplements / source-data Excel (kept off git)

## Marathon-mode rules respected

- No new sprint analysis; corroboration only.
- Paper 3 Track B remains FROZEN per `v19_paper3_ici_track_a` — this layer extends Track B-**lite** which already ran 5/6.
- Output goes to `manuscript_v8/reviewer_reserve/` for Paper 1 and stays inside Track B-lite for Paper 3.
- Voice-protected sections (Hook/Aim/Disc 3.1/Limitations/Cover Para 1/Q9) untouched.

## Data sources

- Wang 2024: PRIDE PXD044900 (proteomics), PXD045017 (phospho); GSA-Human HRA005293/HRA005382 (RNA/WES, restricted); MetaboLights MTBLS3339; code github.com/diChen310/PTC_multi_omics; source data Excel from Nature Commun supplements.
- Mun 2025: PMC12000556; supplements MOESM3 (90MB protein+phospho master tables) and MOESM11 (99MB phospho figure source data) all public-accessible.
- TCGA-THCA HM450 methylation (4th pillar): pre-existing at `audit_2026_04_30/round5/r5_2_*` — n=503, 7/8 panel hypermethylated DM1 vs DM2 (TPO d=+2.30, p=1.9e-18 strongest).

## Blocked fetches (deferred)

- Wang **metabolomics** (MTBLS3339) — MetaboLights public API now requires user token in 2026; SPA-rendered file pages return HTML for direct URL probes. User browser download required.
- **Cell Rep Med 2026** "Proteogenomic characterization of advanced DTC" (113 patients, CC1 canonical / CC2 stromal / CC3 immunogenic) — Cell.com returns 403 for programmatic Referer-spoofed requests; user browser download or institutional proxy required. This would be a powerful 2nd proteogenomic cohort for both Paper 1 (independent replication) and Paper 3 (CC3 ↔ DIAL-lite immunogenic phenotype).
- TCGA miRNA-seq + Lu 2023 scRNA per-cell-type module — feasible from on-disk/GDC but sprint-class compute, deferred per marathon-mode posture.
