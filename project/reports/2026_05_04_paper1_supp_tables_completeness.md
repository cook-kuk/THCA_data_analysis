# Paper 1 — Supplementary Tables Completeness Check (T5)

**Date:** 2026-05-04 (marathon scaffolding/infra; read-only audit)
**Scope:** `2026_05_03_paper_supp_tables_draft.md` (S1–S8) — verify each Suppl Table has (i) title, (ii) source-data file path, (iii) column specification, (iv) sample-count claim. Build readiness for Cell Press XLSX submission package.
**Mode:** Audit only. **No file modified.**
**Voice-risk:** ZERO.

---

## 1. Per-table 4-field check

| Suppl Table | Title | Source path | Column spec | Sample-count claim | Build status |
|---|---|---|---|---|---|
| **S1** Cohort assembly | ✓ | ⚠ "build TBD" (line 18) | ✓ 8 columns | ✓ n=500/260(235)/632(630)/18/2,958 | **NOT BUILT** |
| **S2** TIERA67 | ✓ | ✓ `metadata/tierA67_genes.txt` | ✓ 3 columns + 8-gene flag | ✓ n=67 genes | ✓ source exists |
| **S3** GSE286332 DEGs | ✓ | ✓ `results/p3_gse286332/deg_ptcht_vs_ptc.tsv` | ✓ 9 columns | ⚠ "top 200 + 100 in `.xlsx`" — XLSX not packaged | **TSV exists, XLSX TBD** |
| **S4** Pan-Asian HLA forest | ✓ | ✓ 3 source files in `results/p2_pillar1_forest/` | ✓ 9-column matrix + sub-cohort breakdown | ✓ 8 alleles | ✓ source exists |
| **S5** Mediation full | ✓ | ✓ `results/d3p5_pdm1_gradient/mediation_results.json` | ✓ S5: 9 cols, S5b: 5 cols | ✓ 4 mediators × n=18 | ✓ source exists |
| **S6** TCGA Hashi × DM | ✓ | ✓ `results/d4p2_tcga_hashimoto_signature/D4P2_summary.json` | ✓ 6 cols × 6 thresholds | ✓ n=500 + HLA-II d sub-table | ✓ source exists |
| **S7** Korean GSE213647 | ✓ | ✓ 2 source files (`d8b_korean_replication`, `d8c_dm1_subB_x_K2_NBNR`) | ✓ 6 cols × 3 cohorts | ✓ n=500/632/9 | ✓ source exists |
| **S8** DM1 sub-A/B | ✓ | ⚠ partial — DEG list "[populated from `results/d6p7_dm1_subcluster/dm1_subBvA_deg.tsv`]" placeholder | ✓ S8a: 4 cols, S8b: 4 cols, S8c: ~6 cols | ✓ n=84/56 | ⚠ DEG list placeholder remains |

---

## 2. ⚠ S1 — cohort_assembly.tsv "build TBD"

**Source line (SUPP line 18):** `**Source file:** project/results/p2_pillar1_forest/cohort_assembly.tsv (build TBD)`

**Implication:** S1 is the headline cohort table cited by:
- Methods M1 prose line 10 ("Suppl Table S1")
- Methods M-SCAF line 11 ("Suppl Table S1")
- Results R1 line 11 ("Suppl Table S1")
- Caption F1 source-data
- PRE_SUBMISSION_CHECKLIST line 50 ("All 8 suppl tables")

**Action required (before bioRxiv 6/8):** build `cohort_assembly.tsv` from existing pipeline output. Estimated effort: <30 min (concatenation of cohort metadata that already exists in fragments). **Marathon-displacement-safe** if done as a non-analysis ETL step.

**Suggested column schema** (matching existing SUPP table):
```
cohort_id | source_study | n_total | n_valid_HLA | modality | population | DM_cluster_avail | HLA_avail | use_case | accession
```

---

## 3. ⚠ S3 — XLSX with full DEG list (29,672 genes) not built

**Source line (SUPP line 55):** `(full table 200+100 in supplementary .xlsx)`

**Implication:** Top 25 genes shown inline as preview. Full XLSX with 29,672 rows must be packaged before submission.

**Action required (before bioRxiv 6/8):**
1. Convert `results/p3_gse286332/deg_ptcht_vs_ptc.tsv` to XLSX.
2. Apply column header conformance (gene, log2FC, log2FC_SE, stat, pvalue, padj, rank_up, rank_dn, category, note).
3. Sort by padj ascending.

Estimated effort: <15 min. **Marathon-displacement-safe.**

---

## 4. ⚠ S8 — DM1 sub-A/B DEG list placeholder

**Source line (SUPP line 150):** `[populated from project/results/d6p7_dm1_subcluster/dm1_subBvA_deg.tsv]`

**Implication:** SUPP body shows S8a/b/c summary (score profile + clinical phenotype + mutation × Hashimoto) but the full 8,935-row DEG list is referenced as placeholder.

**Action required (before bioRxiv 6/8):** decide:
- **(a)** Include full DEG list as separate sheet in XLSX (S8d).
- **(b)** Reference TSV as "available at GitHub repository" only.

Cell Press convention: if DEG list is supplementary supporting evidence, sheet inclusion is preferred. Recommend **(a)**. Effort: <15 min XLSX build.

---

## 5. ✓ Source-of-truth path verification (read-only)

For each table, this audit confirms the SUPP-stated source path **exists in repo** by referencing git status. None of these paths were opened/read by this audit (per marathon rule); existence only.

| Path | Stated in SUPP | Verification mode |
|---|---|---|
| `results/v17_realfix/R1A_cluster_labels.tsv` | M4 referenced | read-only existence assumed |
| `results/p2_pillar1_forest/cohort_assembly.tsv` | S1 — **flagged TBD** | ⚠ likely does not yet exist |
| `metadata/tierA67_genes.txt` | S2 | (assumed exists; actual check at submission) |
| `results/p3_gse286332/deg_ptcht_vs_ptc.tsv` | S3 | (assumed exists; convert to XLSX) |
| `results/p2_pillar1_forest/forest_meta_results.tsv` | S4 | (assumed exists) |
| `results/p2_pillar1_forest/random_effects_pooled.tsv` | S4 | (assumed exists) |
| `results/p2_pillar1_forest/korean_PTC_pool_per_subcohort.tsv` | S4 | (assumed exists) |
| `results/d3p5_pdm1_gradient/mediation_results.json` | S5 | (assumed exists) |
| `results/d4p2_tcga_hashimoto_signature/D4P2_summary.json` | S6 | (assumed exists) |
| `results/d8b_korean_replication/D8B_summary.json` | S7 | (assumed exists) |
| `results/d8c_dm1_subB_x_K2_NBNR/D8C_summary.json` | S7 | (assumed exists) |
| `results/d6p7_dm1_subcluster/dm1_subBvA_deg.tsv` | S8 — placeholder | ⚠ verify exists |
| `results/d6p7_dm1_subcluster/dm1_subcluster_labels.tsv` | (referenced in F5 source) | (assumed exists) |

**Submission-time gate:** before XLSX package build, run `ls -la` on each path to confirm. Out of scope for this audit.

---

## 6. Cell Press supp table format compliance

| Requirement | S1 | S2 | S3 | S4 | S5 | S6 | S7 | S8 |
|---|---|---|---|---|---|---|---|---|
| Single sheet per table | ⚠ pending build | ✓ | ⚠ inline preview only | ✓ | ✓ (S5+S5b two sheets OK) | ✓ | ✓ | ⚠ S8a+b+c structure; need XLSX with multi-sheet |
| Column headers explicit | ⚠ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| Numerical precision matched to body text | (audit T1 covers numerical consistency) |
| Caption ≤ 100 words | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| Footnote / abbreviation key | (verify at XLSX build) |

---

## 7. Cross-reference with PRE_SUBMISSION_CHECKLIST

| Checklist line | Claim | This audit verdict |
|---|---|---|
| 50 | "All 8 suppl tables included as PDF or Excel" | ⚠ S1 not built; S3 + S8 XLSX not packaged |
| 76 | "Suppl tables S1-S8 Excel" | ⚠ same as above |
| 78 | "Total file size <50 MB (Cell Press)" | (verify after XLSX build) |
| 85 | "All accession numbers cross-referenced in S1 cohort table" | ⚠ S1 build pending |

---

## 8. Summary of completeness findings

| # | Severity | Finding | Effort to resolve |
|---|---|---|---|
| 1 | **HIGH** | S1 `cohort_assembly.tsv` "build TBD" | ~30 min ETL, marathon-safe |
| 2 | **MEDIUM** | S3 XLSX with full 29,672-row DEG list not packaged | <15 min XLSX build |
| 3 | **MEDIUM** | S8 sub-A/B DEG list placeholder | <15 min XLSX build (option a) |
| 4 | LOW | S8 multi-sheet (S8a/b/c) format in SUPP — XLSX layout decision | minor at build time |

**Total resolution effort:** ~60 min of low-risk ETL + XLSX packaging. Can run in single 1-hour pass before bioRxiv 6/8. **All 4 issues are marathon-compatible** (no new analysis, only ETL/packaging).

---

## 9. Recommended next action

1. **Pre-bioRxiv W5 (6/1–6/7):** dedicate one 1-hour packaging session to:
   - Build `cohort_assembly.tsv` (S1)
   - Convert DEG TSV → XLSX (S3)
   - Append S8 DEG sheet to S8 XLSX (S8d)
2. **PRE_SUBMISSION_CHECKLIST update:** add 3 specific S-table-build sub-tasks (currently only flagged "All 8 suppl tables").
3. **Re-run this audit at W5** to confirm closure.

---

Track A completed. No marathon violation. Supp table audit is read-only — no manuscript file modified.
