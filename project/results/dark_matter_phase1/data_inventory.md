# Dark Matter Phase 1 — Data Inventory
**Date:** 2026-04-29

## Cohort × variable matrix

| Cohort | Expression | Clinical | BRAF V600E | RAS hotspot | TERT prom | DICER1/EIF1AX | Fusion | 8-gene cluster |
|---|---|---|---|---|---|---|---|---|
| **TCGA-THCA** | ✅ 572 × 51,711 (log2 TPM) | ✅ 506 OS, 16 events | ✅ `tcga_thca_mutation_groups.tsv` | ✅ same | ✅ 36 from `FINAL_tert_status_integrated.tsv` | ⚠ in `mutation_genes` col, 7 calls | ⚠ `dark_matter_fusion_overlay.tsv`, sparse | ✅ `v17_dark_cluster` DM1/DM2 |
| **GSE33630** | ❌ NOT IN WAREHOUSE | ❌ | ❌ | ❌ | ❌ | ❌ | ❌ | ❌ |
| **K2 / PRJEB11591 / Yoo 262** | ⚠ 8-gene panel only | ⚠ SRA run metadata only | ❌ | ❌ | ❌ | ❌ | ❌ | ✅ 260 predictions |
| **분당 SNUH** | ❌ OUTREACH stage | ❌ | ❌ | ❌ | ❌ | ❌ | ❌ | ❌ |
| **GSE241184 (sc)** | ✅ 30,660 cells × 38,349 genes | ⚠ 3 patients only (1T/1N/1LN) | ❌ | ❌ | ❌ | ❌ | ❌ | ⚠ per-cell DM score, no sample call |
| GSE126698 | ✅ 28 × 57,773 | ✅ partial | ❌ | ❌ | ❌ | ❌ | ❌ | ⚠ score only |
| GSE213647 | ✅ 632 × 51,389 | ✅ stage/age/sex/outcome | ❌ | ❌ | ❌ | ❌ | ❌ | ✅ predictions |
| GSE76039 | ✅ 37 × 22,880 (microarray) | ✅ partial | ❌ | ❌ | ❌ | ❌ | ❌ | ✅ predictions |

## Critical files used in Phase 1 D0
- `project/results/tables/tcga_thca_mutation_groups.tsv` — BRAF/RAS calls (468 samples, has_braf_v600e/has_ras_mut as bool)
- `project/results/v17/tables/sample_master_v17_full.tsv` — cluster + driver_anchor_v17 + mutation_genes (1510 rows multi-cohort)
- `project/results/v17_tert_recovery/v2/FINAL_tert_status_integrated.tsv` — 36 TERT promoter mutations
- `project/results/tables/tcga_thca_clinical_extended.tsv` — OS data (506 samples)

## Blockers for full multi-cohort Dark Matter test
1. **K2 BRAF/RAS calls** — NOT EXTRACTED. Path: re-process K2 RNA-seq with variant calling pipeline (`bcftools` from STAR-aligned BAMs), or pull from Yoo 2016 supplementary tables. **2-4 days work** if matched normals exist.
2. **GSE33630** — NOT IN WAREHOUSE. Either acquire from GEO or drop from analysis. **1 day** to acquire.
3. **분당 cohort** — OUTREACH (emails drafted, not sent per memory). **No timeline yet.**
4. **TCGA-CDR PFI** — NOT INTEGRATED. Liu 2018 *Cell* SuppTable provides ~30+ PFI events for THCA. **<1 day** to pull and integrate. **Highest priority follow-up.**
5. **DICER1/EIF1AX systematic calling** — only 7 calls aggregated in TCGA `mutation_genes`. Underlying MAFs are at `data_raw/gdc/TCGA-THCA/` — could be re-parsed for completeness. **1 day** if MAFs are present.

## Working assumption for Phase 1
TCGA-only analysis is sufficient to make a **GO / NOGO call on the Dark Matter framing**. Multi-cohort validation is a Phase 2 milestone, not a Phase 1 viability blocker.
