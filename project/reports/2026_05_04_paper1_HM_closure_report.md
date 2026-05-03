# Paper 1 — HIGH/MEDIUM Issue Closure Report (2026-05-04)

**Author:** Seungho Cook
**Date:** 2026-05-04 (marathon scaffolding/infra session, Day 1 of 6-week sprint to bioRxiv 6/13)
**Mode:** Diff-friendly batch edit per audit `PAPER1_MARATHON_AUDIT_BUNDLE_2026_05_04.md`.
**Scope:** HIGH (H1, H2) + MEDIUM (M1, M2, M5, M6) closures. LOW deferred to cosmetic pass.

---

## 1. Decisions applied (user-approved)

| ID | Issue | Decision | Status |
|---|---|---|---|
| **H1** | DM1 sub-B "96%" mutation-negative inconsistency | candidate (a): use 53/56 = 94.6%; remove all "96%"; SUPP T8 counts maintained | ✓ closed |
| **M1** | sub-A 69% vs 61% denominator ambiguity | denominator disambiguation: prose long form + caption short form | ✓ closed |
| **H2** | S1 cohort_assembly.tsv "build TBD" | build TSV from existing metadata fragments (no analysis) | ✓ closed |
| **M2** | F6 / F7 not explicitly named in R-PRO body | add (Figure 6) to R5b, (Figure 7) to R5d | ✓ closed |
| **M5** | S3 full DEG XLSX not packaged | TSV → XLSX (29,672 rows, sorted by padj) | ✓ closed |
| **M6** | S8d sub-A/B DEG XLSX sheet not packaged | TSV → multi-sheet XLSX (S8a/b/c/d) | ✓ closed |

Total HIGH: 2/2. Total MEDIUM: 4/6 closed; **M3 + M4 deferred** (require GEO/PubMed lookup, marathon-displacement-safe to defer).

---

## 2. Files modified (changed) — 5 files

| File | Lines changed | Edits |
|---|---:|---|
| `project/reports/2026_05_03_manuscript_v8_OUTLINE.md` | 2 lines | H1 Abstract, H1+M1 R5d structure |
| `project/reports/2026_05_03_results_R1_R5_prose.md` | 2 paragraphs | H1+M1 R5d (sub-A 69%/61% disambiguation, sub-B 53/56 (94.6%)), M2 R5b (Figure 6 callout), M2 R5d (Figure 7 callout) |
| `project/reports/2026_05_03_figure_captions_all.md` | 2 panels | H1+M1 F5 panel D, H1+M1 F7 panel A |
| `project/reports/2026_05_03_reviewer_QA_consolidated.md` | 2 places | H1+M1 Q12 source-evidence bullets, H1 Q12 answer body |
| `project/reports/2026_05_03_paper_supp_tables_draft.md` | 3 source-path lines | H2 S1 source ("build TBD" → built), M5 S3 XLSX path, M6 S8 multi-sheet XLSX path |

**Voice-protected sections untouched:** Hook ¶1, Aim ¶4, Discussion §3.1, Discussion §3.4 Limitations, Cover letter ¶1, Reviewer Q9 — all verified untouched in this diff. ✓

**`git diff --stat HEAD`:**
```
project/reports/2026_05_03_figure_captions_all.md      | 11 ++++++-----
project/reports/2026_05_03_manuscript_v8_OUTLINE.md    |  4 ++--
project/reports/2026_05_03_paper_supp_tables_draft.md  | 17 ++++++++++++-----
project/reports/2026_05_03_results_R1_R5_prose.md      |  4 ++--
project/reports/2026_05_03_reviewer_QA_consolidated.md |  6 +++---
5 files changed, 25 insertions(+), 17 deletions(-)
```

Net change: +8 lines across 5 files. Diff-friendly. ✓

---

## 3. Files created (new) — 3 artifacts

| Path | Type | Size | Source |
|---|---|---:|---|
| `project/results/p2_pillar1_forest/cohort_assembly.tsv` | TSV | 1,830 bytes | H2 build from existing metadata (5 cohorts × 10 columns) |
| `project/results/p3_gse286332/SuppTable_S3_GSE286332_DEGs.xlsx` | XLSX | 2,236,160 bytes | M5 packaging from `deg_ptcht_vs_ptc.tsv` (29,672 rows, sorted by padj asc) |
| `project/results/d6p7_dm1_subcluster/SuppTable_S8_DM1_subcluster.xlsx` | XLSX | 3,088,130 bytes | M6 packaging — 4 sheets (S8a 7×5, S8b 2×5, S8c 2×11, S8d 51,711×5) |

---

## 4. Per-issue change detail

### H1 — DM1 sub-B "96%" → "53/56 (94.6%)" propagation

**Replacement strings used:**

| Location | Old | New |
|---|---|---|
| OUTLINE Abstract | `DM1 sub-B = 96% mutation-negative NBNR cluster.` | `DM1 sub-B = 53/56 (94.6%) mutation-negative NBNR cluster.` |
| OUTLINE R5d (line 170) | `sub-B (n=56, 96% mut-neg)` | `sub-B (n=56; 53/56 mut-neg = 94.6%)` |
| R-PRO R5d | `Sub-B contained only 1 BRAF+ and 2 RAS+ of 49 samples (96% mutation-negative)` | `Sub-B contained only 1 BRAF+ and 2 RAS+ of 56 samples (53/56, 94.6% mutation-negative)` — denominator corrected from 49 → 56 per SUPP T8 |
| CAP F5 panel D | `53 mutation-negative (96%)` | `53 mutation-negative (53/56, 94.6%)` |
| CAP F7 panel A | `sub-B (n=56) 96% mutation-negative` | `sub-B (n=56) 53 mutation-negative (53/56, 94.6%)` |
| QA Q12 source | `sub-B: 2/49 RAS+, 1/49 BRAF+ (96% mutation-negative)` | `sub-B: 2/56 RAS+, 1/56 BRAF+, 53/56 mutation-negative (94.6%)` — denominator corrected |
| QA Q12 answer | `Sub-B (n=56) ... — 96% mutation-negative,` | `Sub-B (n=56) ... — 53/56 (94.6%) mutation-negative,` |

**SUPP T8 cell** `53 (94%)` left unchanged (rounding consistent with itself; counts maintained). User decision: "Count는 SUPP T8 유지: BRAF+ 1, RAS+ 2, mut-neg 53, total 56." ✓

### M1 — sub-A 69% vs 61% denominator disambiguation

**Replacement strings used (per user's exact phrasing):**

| Location | Old | New |
|---|---|---|
| OUTLINE R5d (line 170) | `sub-A (n=84, 69% RAS+)` | `sub-A (n=84; 51/74 mut-tested RAS+ = 69%, 51/84 = 61% of total)` |
| R-PRO R5d | `sub-A contained 51 of 74 RAS-positive samples (69%) and 1 BRAF V600E-positive` | `51/74 mutation-tested sub-A samples were RAS-positive (69%; 51/84, 61% of total sub-A) and 1 was BRAF V600E-positive` |
| CAP F5 panel D | `sub-A: 51 RAS+ (69%), 1 BRAF+` | `sub-A: 51 RAS+ (69% of mutation-tested; 61% of total n=84), 1 BRAF+` |
| CAP F7 panel A | `sub-A (n=84) 69% RAS+ classical FVPTC` | `sub-A (n=84) 51 RAS+ (69% of mutation-tested; 61% of total n=84) classical FVPTC` |
| QA Q12 source | `sub-A: 51/74 RAS+, 1/74 BRAF+ (69% RAS+ FVPTC core)` | `sub-A: 51/74 mutation-tested RAS+, 1/74 BRAF+ (69% of mutation-tested; 51/84 = 61% of total sub-A; FVPTC core)` |

### H2 — `cohort_assembly.tsv` build (no new analysis)

**Built from existing metadata fragments only** (SUPP T1 + STAR Methods Key Resources + manuscript Methods M1).

```
cohort_id   source_study                                        n_total  n_valid_HLA  modality          population            DM_cluster_avail        HLA_avail                       use_case                     accession
TCGA-THCA   Cancer Genome Atlas Research Network 2014 Cell      500      NA           RNA-seq Illumina  EUR-dominant + mixed  DM1=140; DM2=360        HLA-LA imputed                  Discovery + classifier ...   dbGaP phs000178
K2          Yoo SK et al. 2016 Mol Cell Biol; SNU-GMI            260      235          RNA-seq Illumina  Korean                DM1=14; DM2=246         arcasHLA 4-digit (v3.44.0)      Korean PTC validation        ENA PRJEB11591
Lee2024     Lee Y et al. 2024; Macrogen                         632      630          RNA-seq Illumina  Korean                DM proxy via 8-gene     arcasHLA 4-digit (v3.44.0)      Korean PTC replication       GEO GSE213647
GSE286332   Lim DW et al. 2025 Dongguk University               18       9            RNA-seq NovaSeq X Korean                DM2 100% (all 18)       arcasHLA 4-digit (PTC arm)      PTC vs PTC+HT discovery      GEO GSE286332
Chu2018     Chu X et al. 2018 J Med Genet 55(10):685-692        2958     NA           SNP2HLA           Han Chinese           NA                      summary statistics published    Pan-Asian forest replication doi:10.1136/jmedgenet-2017-105146
```

5 cohorts × 10 columns. **No statistical re-computation.** Pure metadata aggregation. ✓

### M2 — Figure 6 / Figure 7 explicit callouts in R-PRO

| Location | Change |
|---|---|
| R-PRO R5b — Korean replication paragraph | added `(Figure 6;` to Suppl Figure S3 / Suppl Table S7 callout |
| R-PRO R5d — Korean GSE213647 sub-B-like replication | added `(Figure 7; ` to Suppl Table S7 callout, and `Figure 7` after `Figure 5D` in mutation-status sentence |

### M5 — S3 XLSX packaging

- Source: `project/results/p3_gse286332/deg_ptcht_vs_ptc.tsv` (29,672 rows, 6 cols: gene, baseMean, log2FoldChange, stat, pvalue, padj)
- Output: `SuppTable_S3_GSE286332_DEGs.xlsx` (2.2 MB; single sheet `S3_DEG_PTC+HT_vs_PTC`; sorted by padj ascending)
- **No DEG re-computation.** Pure TSV → XLSX. ✓

### M6 — S8 multi-sheet XLSX packaging

- Output: `SuppTable_S8_DM1_subcluster.xlsx` (3.1 MB; 4 sheets):
  - **S8a_score_profile** (7 rows × 5 cols) ← `subcluster_score_profile.tsv`
  - **S8b_clinical** (2 rows × 5 cols) ← `subcluster_clinical.tsv`
  - **S8c_mut_x_hashi** (2 rows × 11 cols) ← assembled from SUPP T8 + reconciliation (counts unchanged from SUPP T8; added denominator-disambiguation columns RAS_pos_pct_of_mut_tested vs RAS_pos_pct_of_total)
  - **S8d_subBvA_DEGs** (51,711 rows × 5 cols) ← `dm1_subBvA_deg.tsv` sorted by padj ascending
- **No new DEG analysis.** TSV/JSON aggregation only. ✓

---

## 5. Marathon discipline check

| Constraint | Verified? |
|---|---|
| No new data download | ✓ — only existing TSVs read |
| No new statistical analysis | ✓ — only TSV → XLSX packaging + metadata aggregation |
| No Paper 3/4 file touch | ✓ — `project/reports/paper3_ici/` and Paper 4/GD memories untouched |
| No voice-protected prose modification | ✓ — Hook ¶1, Aim ¶4, Discussion §3.1, Discussion §3.4 Limitations, Cover ¶1, Reviewer Q9 all unchanged |
| Diff-friendly | ✓ — 5 file diff with +25/−17 lines, individual single-line changes |
| H&E-DM1 closure NOT re-attempted | ✓ — Paper 2 territory, untouched |

**Marathon-safe: no new analysis.** ✓

---

## 6. Remaining open issues (this audit cycle)

### Deferred MEDIUM (require external lookup; marathon-displaceable)

| ID | Issue | Why deferred |
|---|---|---|
| M3 | `Lee2024GSE213647` bib author placeholder ("Y and others") | Requires GEO / PubMed lookup, ~15 min, schedule pre-bioRxiv W6 |
| M4 | `Lim2025GSE286332` bib title ("Vitamin D" framing) vs PTC-vs-PTC+HT use | Verify against PMID 41113708, ~5 min, same window as M3 |

### Deferred LOW (cosmetic; W6 batch pass)

| ID | Issue |
|---|---|
| L1 | TIERA67 / pan-genome ARI rounding (0.90 vs 0.903) |
| L2 | TCGA Hashi top-30% Fisher p (6e-10 vs 6.4e-10) |
| L3 | BRAF mRNA Cohen d (-0.04 vs -0.044) |
| L4 | references.bib header "~22" vs actual 29 |
| L5 | ★ Unicode in Chu2018JMG note |
| L6 | Pan2025NatComm pages format |
| L7 | 4 unused bib entries (Tuttle2019, Yi2016, Chen2024, Kim2014) |
| L8 | "★ exceptional" verb in F2D caption — voice polish |
| L9 | Possible Krishnamoorthy 2025 bib gap (depends on §3.1 voice draft) |

### Active cross-paper boundary status

- Paper 3 ICI: Track A FROZEN. **No touch in this session.** ✓
- Paper 4 GD HLA: backlog. **No touch.** ✓
- Paper 2 H&E-DM1: own marathon track. **No touch from Paper 1 closure side.** ✓

---

## 7. Recommended next action

1. **Re-run numerical consistency audit** (`2026_05_04_paper1_numerical_consistency_audit.md`) to confirm H1/M1 fully closed across all 9 audited files. Expected: 0 occurrences of "96%" remaining in Paper 1 prose.
2. **Schedule M3 + M4** (bib author/title verification) as one ~20 min session. Marathon-safe, single-file edit.
3. **W5 (6/1–6/7) cosmetic pass** for L1–L9 in single batch, single 30-min session.
4. **Pre-bioRxiv W6 (6/8–6/13)** voice-protected drafting (user keyboard) — Hook, Aim, Discussion §3.1, Limitations, Cover ¶1, Q9.

---

## 8. Verification commands (post-edit smoke test)

```bash
# Confirm "96%" removed from Paper 1 prose (excluding intentional source notes):
grep -n "96%" project/reports/2026_05_03_manuscript_v8_OUTLINE.md \
                project/reports/2026_05_03_results_R1_R5_prose.md \
                project/reports/2026_05_03_figure_captions_all.md \
                project/reports/2026_05_03_reviewer_QA_consolidated.md
# Expected: zero matches (or only audit/closure-report self-references).

# Confirm Figure 6 / Figure 7 callouts present in R-PRO:
grep -n "Figure 6\|Figure 7" project/reports/2026_05_03_results_R1_R5_prose.md
# Expected: at least 1 line each.

# Confirm cohort_assembly.tsv exists and parses:
python3 -c "import pandas as pd; df=pd.read_csv('project/results/p2_pillar1_forest/cohort_assembly.tsv', sep='\t'); print(df.shape); print(df['cohort_id'].tolist())"
# Expected: (5, 10) ['TCGA-THCA', 'K2', 'Lee2024', 'GSE286332', 'Chu2018']

# Confirm XLSX files readable:
python3 -c "import pandas as pd; print(pd.ExcelFile('project/results/p3_gse286332/SuppTable_S3_GSE286332_DEGs.xlsx').sheet_names)"
python3 -c "import pandas as pd; print(pd.ExcelFile('project/results/d6p7_dm1_subcluster/SuppTable_S8_DM1_subcluster.xlsx').sheet_names)"
# Expected: ['S3_DEG_PTC+HT_vs_PTC'] and ['S8a_score_profile', 'S8b_clinical', 'S8c_mut_x_hashi', 'S8d_subBvA_DEGs']
```

---

Marathon-safe: no new analysis. No marathon violation. Paper 3/4 untouched. Voice-protected sections untouched.
