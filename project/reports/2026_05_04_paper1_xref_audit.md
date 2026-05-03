# Paper 1 — Cross-Reference Audit (T3)

**Date:** 2026-05-04 (marathon scaffolding/infra; read-only audit)
**Scope:** Every `Figure N`, `Suppl Figure SN`, `Suppl Table SN` mention in prose ↔ caption / supp-table definitions.
**Mode:** Audit only. **No file modified.**
**Voice-risk:** ZERO.

---

## 1. Files audited

| Tag | File | Role |
|---|---|---|
| OUTLINE | `2026_05_03_manuscript_v8_OUTLINE.md` | structure + figure list |
| M-PRO | `2026_05_03_methods_M1_M11_prose.md` | Methods cite Suppl Tables |
| M-SCAF | `2026_05_03_methods_M1_M11_scaffold.md` | Methods scaffold cite Suppl |
| R-PRO | `2026_05_03_results_R1_R5_prose.md` | Results cite Figures + Suppl |
| CAP | `2026_05_03_figure_captions_all.md` | Figure + Suppl Figure + Suppl Table captions |
| SUPP | `2026_05_03_paper_supp_tables_draft.md` | Suppl Table specs |
| QA | `2026_05_03_reviewer_QA_consolidated.md` | reviewer Qs cite source files |
| CHK | `2026_05_03_PRE_SUBMISSION_CHECKLIST.md` | submission asset list |

---

## 2. Main figures F1–F7 — STATUS: ✓ ALL DEFINED + REFERENCED

| Figure | Defined in CAP | OUTLINE list | R-PRO citation | Status |
|---|---|---|---|---|
| F1 (Pan-Asian HLA forest) | line 11 | line 206 | line 15 ("Figure 1A; Suppl Table S4") | ✓ |
| F2 (GSE286332 multi-panel) | line 39 | line 207 | line 23 (panels A/B/C/D) | ✓ |
| F3 (Driver mRNA neutrality) | line 73 | line 208 | line 31, 33, 35 (3A/B/C/D) | ✓ |
| F4 (Pan-genome ARI) | line 104 | line 209 | line 41 (4A), 45 (4B) | ✓ |
| F5 (Autoimmune-PTC mechanism) | line 125 | line 210 | line 55, 57, 67, 73 (5A/B/C/D) | ✓ |
| F6 (TCGA Hashi cross-cohort) | line 163 | line 211 | line 61–63 implied | ⚠ R-PRO mentions cross-cohort generalization but does not call out "Figure 6" by name (see §6) |
| F7 (DM1 sub-A/sub-B + K2 NBNR) | line 182 | line 212 | line 73 implied | ⚠ R-PRO mentions sub-cluster but does not call out "Figure 7" by name (see §6) |

**Panel-level inventory:**
- F1A, F1B ✓
- F2A, F2B, F2C, F2D ✓
- F3A, F3B, F3C, F3D ✓
- F4A, F4B ✓
- F5A (incl. inset), F5B, F5C, F5D ✓
- F6A, F6B ✓
- F7A, F7B, F7C ✓

---

## 3. Supplementary figures SF1–SF6 — STATUS: ✓ ALL DEFINED

| Suppl Figure | Defined in CAP | OUTLINE list | R-PRO / M citation | Status |
|---|---|---|---|---|
| SF1 (TIERA67 pan-genome rank) | line 207 | line 215 | M10 prose line 64 ("Suppl Figure S1") | ✓ |
| SF2 (Mediation Baron-Kenny detail) | line 217 | line 216 | R-PRO R5a line 57 ("Suppl Figure S2") | ✓ |
| SF3 (Korean GSE213647 replication) | line 228 | line 217 | R-PRO R5b line 63 ("Suppl Figure S3") | ✓ |
| SF4 (BCR repertoire detail) | line 246 | line 218 | R-PRO R5c line 67 ("Suppl Figure S4") | ✓ |
| SF5 (DM1 sub-cluster phenotype) | line 261 | line 219 | R-PRO R5d line 73 ("Suppl Figure S5") | ✓ |
| SF6 (Pan-Asian sensitivity 4-scenario) | line 280 | line 220 | R-PRO R1 line 13 ("Suppl Figure S6") + M9 line 58 | ✓ |

---

## 4. Supplementary tables S1–S8 — STATUS: ✓ ALL DEFINED + ✓ ALL REFERENCED

| Suppl Table | Defined in SUPP | Defined in CAP | Cited by | Status |
|---|---|---|---|---|
| S1 (Cohort assembly) | line 8 | line 298 | M1 prose line 10, M-SCAF line 11, R1 prose line 11, F1 caption | ✓ |
| S2 (TIERA67 candidate pool) | line 22 | line 310 | M2 prose line 14, R3 prose line 31, F1 caption | ✓ |
| S3 (Top DEGs GSE286332) | line 40 | line 321 | M3 prose line 18 (implied), R2 prose line 21, F2 caption | ✓ |
| S4 (Pan-Asian HLA forest full) | line 61 | line 333 | R1 prose line 13, M9 prose line 56, F1 caption | ✓ |
| S5 / S5b (Mediation) | line 83 | line 344 | R5a line 55, M6 prose line 36, F5 caption | ✓ |
| S6 (TCGA Hashi × DM cross-tab) | line 106 | line 358 | R5b prose line 61, F5 caption | ✓ |
| S7 (Cross-cohort generalization) | line 128 | line 369 | R5b prose line 63, R5d prose line 75, M8 prose line 48 | ✓ |
| S8 (DM1 sub-A/B characterization) | line 142 | line 379 | R5d prose line 73 (implied), M8 prose line 48 | ✓ |

**Sub-table check:**
- S5a + S5b — both defined in SUPP + CAP ✓
- S8a + S8b + S8c — defined in CAP line 379 as "S8a/b/c"; SUPP line 142 only labels "S8" (combines a/b/c into one section but content covers all three sub-tables). Slight asymmetry but consistent.

---

## 5. Source-data path references in CAP — STATUS: ⚠ 1 path may need build

| Caption | Source data path | Path exists in repo? |
|---|---|---|
| F1 | `results/p2_pillar1_forest/{forest_meta_results.tsv, random_effects_pooled.tsv, korean_subcohort_heterogeneity.tsv}` | (not verified — open at submission) |
| F2 | `results/p3_gse286332/{deg_ptcht_vs_ptc.tsv, gsea_MSigDB_Hallmark_2020.tsv, 8gene_panel_per_sample.tsv, hla_module_scores.tsv}` | (not verified) |
| F3 | `results/p1_driver_mrna_audit/{driver_mrna_mutation_audit.tsv, driver_mrna_dm_auc.tsv, top20_by_d_full_tiera67.tsv}` | (not verified) |
| F4 | `results/p4_pangenome_vs_tiera67/{ari_comparison.tsv, topN_coverage_ladder.tsv}` | (not verified) |
| F5 | `results/d3p5_pdm1_gradient/...`, `results/d4p2_tcga_hashimoto_signature/...`, `results/d5p6_bcr_repertoire/...`, `results/d6p7_dm1_subcluster/...` | (not verified) |
| F6 | `results/d4p2_tcga_hashimoto_signature/`, `results/d8b_korean_replication/` | (not verified) |
| F7 | `results/d6p7_dm1_subcluster/`, `results/d8c_dm1_subB_x_K2_NBNR/` | (not verified) |
| S1 SUPP line 18 | `project/results/p2_pillar1_forest/cohort_assembly.tsv` (build TBD) | ⚠ explicit "build TBD" |

**Action:** Before submission package is built, verify each source-data path exists. ⚠ S1 explicitly marked "build TBD". (T5 audit covers this.)

---

## 6. ⚠ MINOR FINDING — F6 / F7 not explicitly named in R-PRO

**Observation:** Results R-PRO body references F1–F5 by name (e.g., "Figure 1A", "Figure 2C") but never calls "Figure 6" or "Figure 7" by name even though the OUTLINE structure lists them as main figures.

**Examination of R-PRO line by line:**
- R1 → cites F1A, F1B, S4, S6 ✓
- R2 → cites F2A, F2B, F2C, F2D, S3 ✓
- R3 → cites F3A, F3B, F3C, F3D ✓
- R4 → cites F4A, F4B, SF1 ✓
- R5a → cites F5A, F5A inset, S5a/b, SF2 ✓
- R5b → cites F5B, S6, SF3 ✓ (does not call F6 by name)
- R5c → cites F5C, S5/S5b, SF4 ✓
- R5d → cites F5D, S7, SF5 ✓ (does not call F7 by name)

**Reconciliation candidates:**
- **(a)** Add explicit "(Figure 6)" / "(Figure 7)" callouts in R5b / R5d. Single-line edits in R-PRO.
- **(b)** Promote F6 / F7 content into the F5 multi-panel caption (collapse 3 figures into 1 super-figure). Reduces figure count from 7 → 5.
- **(c)** Demote F6 / F7 to suppl figures. Then SF7 (= old F6) and SF8 (= old F7), with R-PRO cite changes.

The F5 caption is already 4-panel + dense; further consolidation may exceed Cell Press caption length limits. **Recommend (a)** — add explicit callouts in R-PRO. Diff-friendly.

**Severity:** MEDIUM. Affects whether F6 / F7 are *referenced* in main text.

---

## 7. ⚠ MINOR FINDING — Caption F1 reference to S6 (sensitivity)

**Observation:** F1 caption (line 35) source data references include `random_effects_pooled.tsv` and `korean_subcohort_heterogeneity.tsv`, but the F1 caption body (line 22-33) describes only the headline forest plot. The 4-scenario sensitivity is in **SF6** caption (line 280) and cited in R1 (line 13).

**Verdict:** F1 + SF6 are correctly separated. No action needed.

---

## 8. Cross-checking PRE_SUBMISSION_CHECKLIST claims

| Checklist line | Claim | Verified? |
|---|---|---|
| 49 | "All 7 main figures + 6 suppl figures embedded as PDFs" | ✓ definitions match (F1–F7, SF1–SF6) |
| 50 | "All 8 suppl tables included as PDF or Excel" | ✓ S1–S8 defined |
| 71 | "Methods (M1-M11) detailed" | ✓ M1–M11 + M5b in scaffold and prose |
| 73 | "STAR Methods Key Resources Table (Excel)" | ✓ STAR Methods file exists |
| 74 | "CRediT author contributions statement" | (not verified — see Acknowledgments file) |
| 75 | "Figures F1-F7 main (300+ DPI PDF)" | (rendering not verified — figure-files state out of audit scope) |
| 76 | "Suppl figures SF1-SF6 PDF" | (rendering not verified) |
| 77 | "Suppl tables S1-S8 Excel" | ⚠ S1 cohort_assembly.tsv "build TBD" per SUPP line 18; S3 XLSX not explicitly built |

---

## 9. Cross-paper boundary check

Paper 1 prose **must not** cite Paper 2 / Paper 3 / Paper 4 figures or supp tables.

**Audit verdict:** ✓ no Paper-2 / Paper-3 / Paper-4 figure or supp-table reference detected in any Paper 1 file.

---

## 10. Summary of audit findings

| # | Severity | Issue | Action |
|---|---|---|---|
| 1 | **MEDIUM** | F6 / F7 not explicitly cited by name in R-PRO body | Add callouts to R5b / R5d or consolidate figures |
| 2 | **LOW** | S1 cohort_assembly.tsv "build TBD" | Build before submission |
| 3 | **LOW** | Source-data paths in F1–F7 captions not yet verified to exist on disk | Verify before submission |
| 4 | LOW | S8 referenced as "S8" in SUPP but "S8a/b/c" in CAP | Acceptable; keep consistent header |
| 5 | LOW | S3 XLSX with full DEG list (29,672 genes) not explicitly built | Build before submission |

**No broken references found.** F1–F7, SF1–SF6, S1–S8 are all defined and referenced (with one minor F6/F7 callout gap).

---

## 11. Recommended next action

1. Resolve **#1 (F6/F7 callouts)** in R-PRO — add 2 lines in R5b + R5d. Diff-friendly.
2. Track #2 + #3 + #5 (build & verify source-data paths) in PRE_SUBMISSION_CHECKLIST.
3. Re-run this xref audit after R-PRO is updated to confirm F6/F7 fix lands.

---

Track A completed. No marathon violation. xref audit is read-only — no manuscript file modified.
