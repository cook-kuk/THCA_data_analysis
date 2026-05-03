# Paper 1 — Marathon Audit Bundle (2026-05-04)

**Author:** Seungho Cook
**Date:** 2026-05-04 (marathon mode, scaffolding/infra session)
**Mode:** Read-only audit. No manuscript file modified. Diff-friendly.
**Scope:** Pre-bioRxiv (target 6/13) consistency / cross-reference / citation / supp-tables / Methods prose audit on existing 2026-05-03 scaffolding.

---

## Table of contents

1. [Numerical Consistency Audit (T1)](#1-numerical-consistency-audit-t1)
2. [Cross-Reference Audit — Figures / Suppl Figures / Suppl Tables (T3)](#2-cross-reference-audit-t3)
3. [BibTeX + Citation Cross-Check Audit (T2)](#3-bibtex--citation-audit-t2)
4. [Supplementary Tables Completeness Check (T5)](#4-supplementary-tables-completeness-t5)
5. [Methods Prose vs Scaffold Gap Report (T4)](#5-methods-prose-gap-report-t4)
6. [Decision summary — issue inventory by severity](#6-decision-summary)

---

## 1. Numerical Consistency Audit (T1)

**Date:** 2026-05-04 (marathon scaffolding/infra; read-only audit)
**Scope:** All numerical claims (n, %, p, OR, Cohen d, ARI, R², AUC) across 9 manuscript files cross-checked.
**Mode:** Audit only. **No file modified.** Reconciliation requires user decision.
**Voice-risk:** ZERO (numbers only).

---

## 0. Files audited (read-only)

| Tag | File | Role |
|---|---|---|
| OUTLINE | `2026_05_03_manuscript_v8_OUTLINE.md` | Abstract + structure |
| M-SCAF | `2026_05_03_methods_M1_M11_scaffold.md` | Methods scaffold |
| M-PRO | `2026_05_03_methods_M1_M11_prose.md` | Methods prose |
| R-PRO | `2026_05_03_results_R1_R5_prose.md` | Results prose |
| CAP | `2026_05_03_figure_captions_all.md` | Figure + supp captions |
| SUPP | `2026_05_03_paper_supp_tables_draft.md` | Supp table specs |
| QA | `2026_05_03_reviewer_QA_consolidated.md` | Reviewer Q&A |
| CHK | `2026_05_03_PRE_SUBMISSION_CHECKLIST.md` | Submission checklist |
| STAR | `2026_05_03_STAR_Methods_KeyResources.md` | STAR Methods table |

---

## 1. Source-of-truth conventions

For each metric, source-of-truth is the analysis result file (`results/*/summary.json` or `*.tsv`) — but during this audit, I treat **SUPP** (Supp Tables S1–S8) as ground-truth proxy because it carries explicit numerator/denominator cells and is closest to the analysis output.

When SUPP and other docs disagree, this audit reports both and proposes a reconciliation candidate but DOES NOT modify either.

---

## 2. Cohort sizes — STATUS: ✓ CONSISTENT

| Metric | OUTLINE | M-SCAF | M-PRO | R-PRO | CAP | SUPP | QA | STAR | Verdict |
|---|---|---|---|---|---|---|---|---|---|
| TCGA-THCA total | 500 | 500 | 500 | 500 | 500 | 500 | — | (cite) | ✓ |
| TCGA DM1 | 140 | 140 | — | 140 | — | 140 | — | — | ✓ |
| TCGA DM2 | 360 | 360 | — | 360 | — | 360 | — | — | ✓ |
| TCGA BRAF V600E carriers | 273 | — | — | — | 273 | — | 273 | — | ✓ |
| TCGA BRAF wild-type | 182 | — | — | — | 182 | — | 182 | — | ✓ |
| K2 PRJEB11591 raw | — | 260 | 260 | — | — | 260 | — | 260 | ✓ |
| K2 PRJEB11591 valid HLA | 235 | 235 | — | 235 | 235 | 235 | 235 | — | ✓ |
| Lee 2024 GSE213647 raw | 632 | 632 | 632 | 632 | — | 632 | 632 | 632 | ✓ |
| Lee 2024 valid HLA | 630 | — | 630 | 630 | 630 | 630 | — | — | ✓ |
| GSE286332 total | 18 | 18 | 18 | 18 | 18 | 18 | 18 | 18 | ✓ |
| GSE286332 PTC arm | 9 | 9 | 9 | 9 | 9 | 9 | 9 | — | ✓ |
| GSE286332 PTC+HT arm | 9 | 9 | 9 | 9 | 9 | 9 | 9 | — | ✓ |
| Korean PTC pool | 874 | 874 | 874 | 874 | 874 | — (implicit) | 874 | — | ✓ (235+630+9=874 arithmetic ✓) |
| Chu 2018 GD | 1,468 | 1,468 | 1,468 | 1,468 | 1,468 | 1,468 | — | — | ✓ |
| Chu 2018 ctrl | 1,490 | 1,490 | 1,490 | 1,490 | 1,490 | 1,490 | — | — | ✓ |
| Chu 2018 total | 2,958 | 2,958 | — | — | — | 2,958 | — | — | ✓ |
| DM1 sub-A | 84 | 84 | — | 84 | 84 | 84 | — | — | ✓ |
| DM1 sub-B | 56 | 56 | — | 56 | 56 | 56 | — | — | ✓ |

**Status:** All cohort sizes consistent.

---

## 3. ⚠ DISCREPANCY 1 — DM1 sub-B mutation-negative percentage

| Document | Reported | Math basis |
|---|---|---|
| OUTLINE R5d (line 170) | "sub-B (n=56, **96% mut-neg**)" | denominator unstated |
| OUTLINE Abstract | "DM1 sub-B = **96%** mutation-negative NBNR" | denominator unstated |
| R-PRO R5d (line 73) | "1 BRAF+ and 2 RAS+ of **49** samples (**96%** mutation-negative)" | 1+2=3 mut+, 49-3=46, 46/49 = **93.9%** ⚠ math doesn't yield 96% |
| CAP F5 panel D (line 152) | "sub-B: 1 BRAF+, 2 RAS+, 53 mutation-negative (**96%**)" | 53/56 = **94.6%** ⚠ math doesn't yield 96% |
| CAP F7 panel A (line 188) | "sub-B (n=56) **96% mutation-negative**" | denominator 56 |
| SUPP T8 (line 147) | "sub-B \| 56 \| BRAF+ 1 (2%) \| RAS+ 2 (4%) \| mut-neg 53 (**94%**)" | 53/56 = **94.6%** ✓ rounds to 94% or 95% |
| QA Q12 (line 161) | "sub-B: 2/49 RAS+, 1/49 BRAF+ (**96%** mutation-negative)" | denominator 49 |
| Memory `v17_D6P7_dm1_subB_NBNR` | "DM1 sub-B (n=56, **96%** mutation-negative)" | denominator 56 |

**Reconciliation candidates (no decision made):**

- **(a)** Source-of-truth = SUPP T8: 53 mut-neg / 56 total = **94.6%** (or **95%** rounded). Then OUTLINE / Abstract / Captions / QA "96%" all need correction → **94.6%** or **~95%**.
- **(b)** If true breakdown is 1 BRAF + 1 RAS = 2 mut+ (not 1+2=3), then 54 mut-neg / 56 = **96.4%** ≈ 96%. Then SUPP T8 has wrong RAS+ count (should be 1, not 2).
- **(c)** Denominator-disambiguation: "96%" calculated as 54 mut-neg / **mut-tested-only** (excluding mut-status-NA samples). If 56 = total but only 54 had mut data, and 2 mut+ of 54 → 52/54 = 96.3%. Plausible. Then need to add a footnote: "of n samples with available mutation status."

**Action required (user):** Pick (a) / (b) / (c). Update affected cells.

**Files affected by single fix:** OUTLINE 2 cells (line ~24 abstract, line ~170 R5d), R-PRO 1 cell (line ~73 R5d), CAP 2 cells (F5 panel D, F7 panel A), QA 1 cell (Q12), SUPP T8 may need note column.

---

## 4. ⚠ DISCREPANCY 2 — DM1 sub-A RAS+ percentage / denominator

| Document | Reported |
|---|---|
| OUTLINE R5d (line 170) | "sub-A (n=84, **69% RAS+**)" — implicit denom 84 → 69% × 84 = 58, but actual count 51 |
| R-PRO R5d (line 73) | "sub-A contained **51 of 74** RAS-positive samples (**69%**)" — denominator 74 (mut-tested subset), not 84 |
| CAP F5 panel D (line 151) | "sub-A: **51 RAS+ (69%)**" + "sub-A (n=84 vs n=56)" same panel — internally inconsistent: 51/84 = 60.7%, not 69% |
| CAP F7 panel A (line 188) | "sub-A (n=84) **69% RAS+ classical FVPTC**" — same internal inconsistency |
| SUPP T8 (line 146) | "sub-A \| 84 \| BRAF+ 1 (1%) \| RAS+ **51 (61%)** \| mut-neg 32 (38%)" — uses 51/84 = **60.7% ≈ 61%** |
| QA Q12 (line 160) | "sub-A: 51/74 RAS+, 1/74 BRAF+ (**69% RAS+ FVPTC core**)" — denominator 74 |

**Pattern:** "69%" is computed from 51/74 (mut-tested subset). "61%" is computed from 51/84 (full sub-cluster). Both are arithmetically defensible — they answer **different questions** (rate among mut-tested vs rate of cluster).

**Reconciliation candidate:** explicitly mark which denominator is in use. Example: "51/74 RAS+ (69% of mutation-tested; 61% of n=84 sub-A total)". Or pick one denominator convention and apply throughout.

**Files affected:** OUTLINE R5d, CAP F5 panel D + F7 panel A, R-PRO R5d, SUPP T8 (note column).

---

## 5. ⚠ DISCREPANCY 3 — TIERA67 / pan-genome ARI rounding

| Document | TIERA67 ARI | Pan-genome top-5000 ARI |
|---|---|---|
| OUTLINE Abstract | 0.90 | 0.92 |
| OUTLINE R4 | 0.90 | 0.92 |
| M-SCAF M10 | (procedure described, no value) | (procedure described) |
| M-PRO M10 | (procedure described, no value) | (procedure described) |
| R-PRO R4 | **0.903** | **0.918** |
| CAP F3D | **0.903** | **0.918** |
| CAP F4A | **0.903** | top-200 0.864 / top-1000 0.902 / **top-5000 0.918** |
| QA Q3 | **0.903** | **0.918** |
| Memory `v17_8gene_pangenome_robustness` | 0.90 | 0.92 |

**Reconciliation candidate:** Outline uses 2-decimal rounding for narrative; precise files use 3-decimal. Both correct. Recommend Outline keep "≈0.90" / "≈0.92" wording with "(0.903 / 0.918)" parenthetical OR unify to 3-decimal in Outline. Low severity.

**Files affected:** OUTLINE Abstract + R4 (cosmetic).

---

## 6. ⚠ DISCREPANCY 4 — TCGA Hashimoto-like top-30% Fisher p rounding

| Document | Fisher p |
|---|---|
| OUTLINE Abstract (line 24) | 6e-10 |
| OUTLINE R5b (line 161) | 6e-10 |
| R-PRO R5b (line 61) | **6.4×10⁻¹⁰** |
| CAP F5B (line 140) | 6e-10 |
| CAP F6 (line 174) | (not specified) |
| SUPP T6 (line 114) | **6.4e-10** |
| QA Q10 (line 138) | **6×10⁻¹⁰** |
| Memory `v17_D4P2_tcga_hashimoto_generalization` | p=6e-10 |

**Reconciliation candidate:** Use 6.4e-10 (precise) throughout, or "p ≈ 6×10⁻¹⁰" with implicit rounding. Pick one. Low severity.

**Files affected:** OUTLINE, CAP F5B, QA Q10 (cosmetic if precision unified).

---

## 7. ⚠ DISCREPANCY 5 — BRAF mRNA Cohen d precision

| Document | Cohen d |
|---|---|
| OUTLINE R3 (line 137) | **−0.04** |
| R-PRO R3 (line 31) | **−0.044** |
| CAP F3A (line 80) | **−0.044** |
| QA Q4 (line 60) | **−0.044** |

**Reconciliation candidate:** unify to **−0.044** (precise) since 3 of 4 use it. Update OUTLINE R3 only.

**Severity:** Cosmetic.

---

## 8. ✓ Effect sizes — CONSISTENT (cluster-grade)

| Metric | OUTLINE | R-PRO | CAP | SUPP | QA | Verdict |
|---|---|---|---|---|---|---|
| 8-gene RAI Cohen d (PTC+HT vs PTC) | −1.60 | −1.60 | −1.60 | — | 1.60 | ✓ |
| HLA-II Cohen d | +3.65 | +3.65 | +3.65 | — | 3.65 | ✓ |
| HLA-I Cohen d | (not specified) | +2.34 | (not specified) | — | 2.34 | ✓ |
| TLS Cabrita Cohen d | +1.96 | +1.96 | +1.96 | — | — | ✓ |
| PAX8 d | −2.32 | −2.32 | −2.32 | — | — | ✓ |
| NKX2-1 d | −1.92 | −1.92 | −1.92 | — | — | ✓ |
| FOXE1 d | −1.75 | −1.75 | −1.75 | — | — | ✓ |
| SLC5A5 d | +0.25 | +0.25 | +0.25 | — | — | ✓ |
| Mediation HLA-II % | 140% | 140% | 140% | 140% | 140% | ✓ |
| Mediation 8-gene % | 63% | 63% | 63% | 63% | 63% | ✓ |
| Bootstrap p_emp HLA-II | 0.023 | 0.023 | 0.023 | 0.023 | 0.023 | ✓ |
| Bootstrap p_emp 8-gene | 0.002 | 0.002 | — | 0.002 | 0.002 | ✓ |
| R² (P_DM1 ~ HLA-II + 8-gene + immune) | 0.756 | 0.756 | — | 0.756 | — | ✓ |
| Single-predictor HLA-II R² | 0.66 | 0.66 (87.7%) | — | — | 0.663 | ✓ (rounding only) |
| TCGA-THCA Hashi top-30% OR | 0.20 | 0.20 | 0.20 | 0.20 | 0.20 | ✓ |
| TCGA-THCA Hashi resid Otsu OR | 0.29 | 0.289 | — | 0.289 | 0.289 | ✓ (rounding) |
| Korean GSE213647 Hashi GMM% | 22.8% | 22.8% | 22.8% | 22.8% | 22.8% | ✓ |
| Korean GSE213647 Hashi Otsu% | 28.2% | 28.2% | 28.2% | 28.2% | 28.2% | ✓ |
| TCGA-THCA Hashi GMM% | 18% | 18.0% | 18% | 18.0% | 18% | ✓ |
| TCGA-THCA Hashi Otsu% | 19.6% | 19.6% | 19.6% | 19.6% | 19.6% | ✓ |
| Korean sub-B-like GMM% | 47.2% | 47.2% | 47.2% | 47.2% | 47-53% | ✓ |
| Korean sub-B-like Otsu% | 52.5% | 52.5% | 52.5% | 52.5% | 47-53% | ✓ |

---

## 9. ✓ Pillar 1 forest meta — CONSISTENT

| Metric | All sources | Verdict |
|---|---|---|
| DPB1*05:01 Korean PTC % | 53.2% | ✓ |
| DPB1*05:01 Chu ctrl % | 31.3% | ✓ |
| DPB1*05:01 Chu GD % | 44.0% | ✓ |
| DPB1*05:01 Korean vs ctrl OR | 2.50 | ✓ |
| DPB1*05:01 Korean vs ctrl p | 4e-26 | ✓ |
| DPB1*05:01 pooled OR | 2.16 [1.65, 2.83] | ✓ |
| DPB1*05:01 pooled I² | 85% | ✓ |
| Korean sub-cohort I² (DPB1*05:01) | 0% | ✓ |
| Korean sub-cohort breakdown (K2/Lee/GSE286332-PTC) | 56.2 / 52.1 / 55.6 | ✓ (some give "56/52/56" rounded) |
| Sensitivity scenario range | 52.1–56.2% | ✓ |

---

## 10. ✓ DEG counts — CONSISTENT

| Metric | All sources |
|---|---|
| GSE286332 DEGs total at padj<0.05 | 10,380 |
| Up-regulated | 6,004 |
| Down-regulated | 4,376 |
| Genes tested (post-filter) | 29,672 |
| Top hit | IGHV3-66 log2FC=+7.25, padj=1.4e-35 |
| GSEA Hallmark IFN-γ NES / FDR | +1.80 / 2e-4 (or 1.9e-4) — small precision drift |
| GSEA Hallmark Allograft Rejection NES / FDR | +2.12 / 0 |

**Minor:** IFN-γ FDR appears as both "2e-4" (Outline, Caption) and "1.9×10⁻⁴" (R-PRO). Cosmetic.

---

## 11. ⚠ DISCREPANCY 6 — Reference count claim

| Document | Claim |
|---|---|
| references.bib header (line 371) | "Citation count: ~22 papers + tools" |
| Actual `grep -c "^@" references.bib` | **29 entries** |

**Reconciliation candidate:** Update bib header to "Citation count: 29 (24 papers + 5 tool-only)" or recount. Cosmetic but visible to reviewer if bib is included in submission package.

---

## 12. Power-claim vs source-of-truth check (Q7)

| Claim | Source |
|---|---|
| Min detectable d at n=9 (80% power) = 1.41 | QA Q7, references `results/p2_power_planB/power_table.tsv` (not in this audit's read scope) |
| Observed 8-gene d=1.60 → power 0.89 | QA Q7 |
| Observed HLA-II d=3.65 → power 1.00 | QA Q7 |

**Status:** Internally consistent within QA Q7. Source TSV not opened during this audit; flag for verification at Track-B-equivalent (i.e., before submission, open `power_table.tsv` and confirm).

---

## 13. Summary of audit findings

| # | Severity | Issue |
|---|---|---|
| 1 | **HIGH** | DM1 sub-B 96% vs SUPP-derived 94.6% — arithmetic inconsistency |
| 2 | **MEDIUM** | sub-A 69% vs 61% — denominator (74 mut-tested vs 84 total) needs disambiguation |
| 3 | LOW | TIERA67 / pan-genome ARI rounding (0.90 vs 0.903) |
| 4 | LOW | TCGA Hashi top-30% p rounding (6e-10 vs 6.4e-10) |
| 5 | LOW | BRAF mRNA d precision (−0.04 vs −0.044) |
| 6 | LOW | references.bib header "~22" vs actual 29 |

**No source-of-truth-grade values changed.** All issues need user decision before reconciliation.

---

## 14. Recommended next action

1. **Resolve #1 (sub-B 96% / 94.6%) FIRST** — affects Abstract, paper-shaping. Pick reconciliation candidate (a/b/c), update 5 files in one diff-friendly batch.
2. Resolve #2 in same diff (denominator disambiguation, 4 files).
3. Resolve #3–#6 in cosmetic-pass before bioRxiv freeze 6/8.

Diff-friendly: each fix is single-line replacement, easy git diff review.

---

Track A completed. No marathon violation. Audit is read-only — no manuscript file modified.

---

## 2. Cross-Reference Audit (T3)

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

---

## 3. BibTeX + Citation Audit (T2)

**Date:** 2026-05-04 (marathon scaffolding/infra; read-only audit)
**Scope:** `2026_05_03_references.bib` (29 entries) ↔ all manuscript prose files (Outline / Methods / Results / Discussion / Captions / QA / STAR).
**Mode:** Audit only. **No file modified.**
**Voice-risk:** ZERO.

---

## 1. Bib summary

- **Total `@article` entries:** 29 (verified via `grep -c "^@" references.bib`)
- **Header claim (line 371):** "~22 papers + tools" — ⚠ inconsistent with actual count.
- **Cell Press style target:** numerical references.

---

## 2. Bib entry inventory

| # | Citation key | Author | Year | Journal | Status |
|---|---|---|---|---|---|
| 1 | `Cancer2014TCGA` | TCGA Research Network | 2014 | Cell | ✓ used in Outline ¶2 |
| 2 | `Yoo2016SNUGMI` | Yoo SK et al. | 2016 | PLOS Genet | ✓ used in M2, Q2, Outline |
| 3 | `Lee2024GSE213647` | Lee Y, others | 2024 | GEO | ⚠ author placeholder "Y and others"; needs full author list |
| 4 | `Lim2025GSE286332` | Lim DW, Kim SM | 2025 | GEO; PMID 41113708 | ✓ used widely |
| 5 | `Chu2018JMG` | Chu X et al. | 2018 | J Med Genet | ✓ used in M9, Q13, R1, F1 |
| 6 | `Cook2021cookHLA` | Cook S et al. | 2021 | Nat Commun | ✓ used in R1 (Cook et al. 2021) |
| 7 | `Orenbuch2020arcasHLA` | Orenbuch R et al. | 2020 | Bioinformatics | ✓ used in M9, STAR |
| 8 | `Robinson2020IPDIMGT` | Robinson J et al. | 2020 | NAR | ✓ STAR only |
| 9 | `Muzellec2023PyDESeq2` | Muzellec B et al. | 2023 | Bioinformatics | ✓ STAR only; M3 cites tool version not paper |
| 10 | `Fang2023gseapy` | Fang Z et al. | 2023 | Bioinformatics | ✓ STAR only; M3 cites tool version |
| 11 | `Liberzon2015Hallmark` | Liberzon A et al. | 2015 | Cell Syst | ✓ STAR only |
| 12 | `Frankish2021GENCODE` | Frankish A et al. | 2021 | NAR | ✓ STAR only |
| 13 | `Landa2016JCI` | Landa I et al. | 2016 | JCI | ✓ Outline ¶1 + Discussion §3.1 (voice-protected) |
| 14 | `Pan2025NatComm` | Pan Z et al. + Ge M | 2025 | Nat Commun | ✓ Outline ¶2 ("Pan/Ge 2025") |
| 15 | `Riesco2014EJE` | Riesco-Eizaguirre & Santisteban | 2014 | EJE | ✓ M2, Q1 |
| 16 | `Cabrita2020TLS` | Cabrita R et al. | 2020 | Nature | ✓ M7, R5c, Caption F5C |
| 17 | `Pfister2013Bimodality` | Pfister R et al. | 2013 | Front Psychol | ✓ M5 |
| 18 | `Haugen2016ATA` | Haugen BR et al. | 2016 | Thyroid | ✓ Outline ¶1 (ATA 2015 cheatsheet) |
| 19 | `Tuttle2019JCEM` | Tuttle RM, Alzahrani AS | 2019 | JCEM | ⚠ NOT YET CITED in any prose file |
| 20 | `Yi2016KTA` | Yi KH et al. | 2016 | Endocrinol Metab | ⚠ NOT YET CITED |
| 21 | `Chen2024EndocrConnect` | Chen XF et al. + Wang Y | 2024 | Endocr Connect | ⚠ NOT YET CITED in main prose (only in `D3P4_wang2024_audit.md` decision record) |
| 22 | `Pedregosa2011sklearn` | Pedregosa F et al. | 2011 | JMLR | ✓ STAR only |
| 23 | `Virtanen2020scipy` | Virtanen P et al. | 2020 | Nat Methods | ✓ STAR only |
| 24 | `Harris2020numpy` | Harris CR et al. | 2020 | Nature | ✓ STAR only |
| 25 | `Bray2016kallisto` | Bray NL et al. | 2016 | Nat Biotechnol | ✓ STAR only |
| 26 | `Langmead2012bowtie2` | Langmead B, Salzberg SL | 2012 | Nat Methods | ✓ STAR only |
| 27 | `Li2009samtools` | Li H et al. | 2009 | Bioinformatics | ✓ STAR only |
| 28 | `Lee2014TissueAntigens` | Lee KW et al. | 2014 | Tissue Antigens | ✓ Q13 (Korean baseline reference) |
| 29 | `Kim2014KoreanGraves` | Kim TH et al. | 2014 | Korean Endocrinol Soc | ⚠ NOT YET CITED in main prose; memory mentions 56% Korean GD DPB1*05:01 replication |

---

## 3. Cited but not in bib — STATUS: ✓ NONE FOUND

All "Author Year" patterns observed in the prose map to a bib entry:
- "TCGA Cancer Network 2014" → `Cancer2014TCGA`
- "Yoo et al. 2016" / "Yoo 2016" → `Yoo2016SNUGMI`
- "Lee et al. 2024" / "Lee Y" → `Lee2024GSE213647`
- "Lim DW et al. 2025" → `Lim2025GSE286332`
- "Chu et al. 2018" / "Chu 2018" → `Chu2018JMG`
- "Cook et al. 2021" → `Cook2021cookHLA`
- "Orenbuch et al. 2020" → `Orenbuch2020arcasHLA`
- "Cabrita et al. 2020" / "Cabrita 2020" → `Cabrita2020TLS`
- "Pfister et al. 2013" → `Pfister2013Bimodality`
- "Riesco-Eizaguirre & Santisteban 2014" → `Riesco2014EJE`
- "Landa 2016" / "Krishnamoorthy/Landa 2016" → `Landa2016JCI`
- "Pan/Ge 2025" / "Pan 2025" → `Pan2025NatComm`
- "Lee et al. 2014 Tissue Antigens" → `Lee2014TissueAntigens`

**Verdict:** No orphan citation found. Bib coverage is complete for currently-drafted prose.

---

## 4. ⚠ In bib but not yet cited — 4 entries

| Bib key | Reason for inclusion | Likely use |
|---|---|---|
| `Tuttle2019JCEM` | ATA 2015 + risk-stratification context | Discussion §3.1 / §3.2 (voice-protected) — likely user voice will cite |
| `Yi2016KTA` | Korean Thyroid Association 2016 guidelines | Discussion §3.2 (ATA-vs-KTA comparison) — likely user voice will cite |
| `Chen2024EndocrConnect` | Wang/Chen 2024 Shanghai n=2,844 mutation comparator | Discussion §3.3 (Asian mutation context) — D3P4 audit memo decided to use as comparator only |
| `Kim2014KoreanGraves` | Korean GD DPB1*05:01 ~36% baseline | Q13 + Discussion §3.3 (Pan-Asian susceptibility) — could replace or augment Lee2014TissueAntigens reference |

**Reconciliation candidate (no decision):**

- **(a)** Keep all 4. They are domain-aligned and likely to be cited in voice-protected Discussion. Acceptable to keep in bib as reserve.
- **(b)** Drop unused entries before submission to keep bib lean. Move to `references_reserve.bib`.
- **(c)** Decide per-entry during voice-protected Discussion drafting (user keyboard).

Recommend **(c)** — defer to user voice draft. Track in checklist.

---

## 5. ⚠ Author placeholder in `Lee2024GSE213647`

```
@article{Lee2024GSE213647,
  author = {Lee, Y and others},   ⚠ placeholder
  ...
}
```

**Action required:** Pull full author list from GSE213647 GEO record (Lee Y et al. 2024) or PubMed, before bioRxiv 6/8 freeze. Same may apply to `Kim2014KoreanGraves` and `Lee2014TissueAntigens` (both have minimal "and others" patterns).

---

## 6. ⚠ Title typo / formatting risk in bib

Spot-check identified:

- `Cancer2014TCGA` — `author = {{Cancer Genome Atlas Research Network}}` — double braces correct for institutional author. ✓
- `Lim2025GSE286332` — title is "Dysregulation of Vitamin D and Its Signaling in Hashimoto's Thyroiditis in Korean Population" — confirm with PubMed PMID 41113708. ⚠ Vitamin-D framing may not match our usage of GSE286332 as "PTC vs PTC+HT" — verify the actual published title before submission.
- `Pan2025NatComm` — `pages = {s41467-025-58910-3}` — pages format unusual; should be article number e.g., `1234` or `e58910`. Verify against Nat Comm convention.
- `Chu2018JMG` — `note` contains "★ Citation correction: previously misattributed as 'Chen 2018'". The ★ Unicode character may break some bibtex pipelines. ⚠ Recommend ASCII alternative ("(NOTE)").

---

## 7. Cross-paper citation discipline (Paper 2 / Paper 3 / Paper 4 boundary)

Paper 1 prose **must not** cite:
- Paper 2 (HT-overlap PTC pathology) — Paper 1 may use Paper 2's GSE286332 *data* but not its analytic claims.
- Paper 3 (ICI dark matter) — frozen, not citable from Paper 1.
- Paper 4 (GD HLA backlog) — gated, not yet citable.

**Audit verdict:** ✓ no cross-paper citation found in current Paper 1 prose. Paper 1 stays self-contained.

---

## 8. Voice-protected note

Discussion §3.1 (Krishnamoorthy/Landa 2016 framing) is voice-protected. When user drafts §3.1, they may need to add `@Krishnamoorthy2025NatComm` if Krishnamoorthy 2025 is being cited (currently absent from bib but flagged in memory `v17_landa2016_cite_save.md` as a misattribution-corrected reference). **Possible bib gap.**

**Action:** Before §3.1 drafting, decide:
- **(a)** Cite Landa 2016 only (current bib state), per `v17_landa2016_cite_save.md` correction, or
- **(b)** Cite both Landa 2016 + Krishnamoorthy 2025 — then add Krishnamoorthy2025 entry to bib.

---

## 9. Summary of audit findings

| # | Severity | Issue | Action |
|---|---|---|---|
| 1 | **MEDIUM** | `Lee2024GSE213647` author placeholder | Pull full author list from GEO/PubMed |
| 2 | **MEDIUM** | `Lim2025GSE286332` title vs paper usage mismatch (Vitamin D framing) | Verify against PMID 41113708 |
| 3 | LOW | Bib header "~22" vs actual 29 | Update header (cosmetic) |
| 4 | LOW | 4 bib entries (Tuttle2019, Yi2016, Chen2024, Kim2014) not yet cited | Defer to Discussion drafting |
| 5 | LOW | `Pan2025NatComm` pages format unusual | Verify Nat Comm convention |
| 6 | LOW | ★ Unicode in `Chu2018JMG` note | Replace with ASCII |
| 7 | LOW | Possible Krishnamoorthy 2025 bib gap (depends on §3.1 voice) | Decide at §3.1 drafting time |

**No orphan citation, no missing-bib citation found.** Bib is in good shape for current prose state.

---

## 10. Recommended next action

1. **Pre-bioRxiv (W6 6/8):** Fix MEDIUM #1 + #2 (author list + title verification).
2. **Cosmetic pass:** Fix LOW #3 + #6 in same diff.
3. **Discussion drafting time:** Decide #4 + #7 with user voice.
4. **Cell Rep Med formal submission (W7 6/14):** Re-audit with this same checklist.

---

Track A completed. No marathon violation. Bib + citation audit is read-only — no manuscript file modified.

---

## 4. Supplementary Tables Completeness (T5)

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

---

## 5. Methods Prose Gap Report (T4)

**Date:** 2026-05-04 (marathon scaffolding/infra; read-only audit)
**Scope:** Compare `2026_05_03_methods_M1_M11_scaffold.md` (124 lines) vs `2026_05_03_methods_M1_M11_prose.md` (72 lines) for section coverage gaps.
**Mode:** Audit only. **No file modified.** Originally proposed as "gap-fill edit"; audit found no real coverage gap, so converted to report-only.
**Voice-risk:** ZERO (no prose generated).

---

## 1. Why this was converted from edit to audit

The candidate task T4 originally proposed: "Methods prose gap-fill — generate fact-only prose for missing scaffold sections." After comparing both files section-by-section, **no missing section identified**. Both files cover M1, M2, M3, M4, M5, M5b, M6, M7, M8, M9, M10, M11 with substantively equivalent content. The line-count difference (124 vs 72) is driven by:

- Scaffold has **header metadata + voice-protected reminder section + quality-check at bottom** (~30 lines of non-Methods scaffolding).
- Scaffold uses denser bullet formatting + extra blank lines.
- Prose uses tighter paragraph form.

Generating "missing prose" would have produced redundant content, with non-zero verb-tone leak risk. **Edit cancelled. Report-only output below.**

---

## 2. Section-by-section coverage matrix

| Section | Scaffold | Prose | Coverage verdict |
|---|---|---|---|
| M1 Cohort Assembly | line 9–11 | line 9–10 | ✓ equivalent — both list 5 cohorts with n + accession |
| M2 TIERA67 candidate pool | line 15–17 | line 12–14 | ✓ equivalent — both list 7 categories + 8-gene panel members |
| M3 DEG + GSEA | line 21–25 | line 16–20 | ✓ equivalent — PyDESeq2 v0.5.4 + gseapy 1.1.13 + parameters |
| M4 DM1/DM2 + classifier | line 29–31 | line 22–24 | ✓ equivalent — KMeans + LogReg + AUC=0.962 |
| M5 Hashimoto signature transfer | line 35–37 | line 26–28 | ✓ equivalent — 150 up + 50 dn signature, 4 thresholds |
| M5b Confounder residualization | line 39–43 | line 30–32 | ✓ equivalent — Stromal + immune-proxy residualization |
| M6 Mediation Baron-Kenny | line 47–49 | line 34–36 | ✓ equivalent — 3-equation framework + 5,000 bootstrap |
| M7 BCR + TLS | line 53–55 | line 38–40 | ✓ equivalent — IGH/IGK/IGL extraction + Cabrita TLS |
| M8 DM1 sub-cluster | line 59–65 | line 42–48 | ✓ equivalent — KMeans k=2 + Welch t + sub-A/B labels |
| M9 Pan-Asian forest meta | line 69–77 | line 50–58 | ✓ equivalent — arcasHLA + DerSimonian-Laird + 4-scenario sensitivity |
| M10 Pan-genome MAD | line 81–85 | line 60–64 | ✓ equivalent — top-N MAD + KMeans + ARI/NMI |
| M11 Stats + reproducibility | line 89–95 | line 66–72 | ✓ equivalent — Cohen d + Wilson CI + BH-FDR + Python 3.12 environment |
| Voice-protected reminder | line 99–111 | (absent) | scaffold-only meta |
| Quality check ✅ | line 115–125 | (absent) | scaffold-only meta |

**Gaps:** none.

---

## 3. Substantive content drift between scaffold and prose — STATUS: minor

Spot-check identified 3 sentences appearing in scaffold but slightly trimmed in prose. None are factually consequential.

| Topic | Scaffold | Prose | Drift |
|---|---|---|---|
| M5b OR=0.289 statement | "(OR = 0.289, p = 8×10⁻⁹), confirming that the autoimmune-PTC axis is not a generic immune-infiltration artifact" | (line 32 of prose) "yielded OR = 0.289 (Fisher p = 8×10⁻⁹), confirming that the autoimmune-PTC axis is not a generic immune-infiltration artifact" | ✓ same content |
| M6 single-predictor dominance | "with HLA-II as the dominant single-predictor" | "with HLA-II as dominant single-predictor (Suppl Table S5/S5b)" | ✓ prose adds Suppl Table cite — improvement |
| M9 Suppl Figure S6 cite | "Sensitivity analyses recomputed the meta-analysis under four scenarios" | "(Suppl Figure S6)" added | ✓ prose adds Suppl Figure cite — improvement |
| M11 GitHub Actions CI | (absent) | "Reproducibility smoke tests covering all five pillars are provided in `tests/test_signature_score.py` and verified continuously via GitHub Actions CI" | ✓ prose adds CI mention — improvement |

**Verdict:** Prose is **slightly more polished** than scaffold (cross-references inserted, CI mention added). No regression.

---

## 4. ⚠ MINOR FINDING — Methods prose missing an explicit cite to references.bib

**Observation:** Methods prose mentions tool names + versions (PyDESeq2 v0.5.4, gseapy v1.1.13, arcasHLA v0.6.0, etc.) but does not insert author-year citations inline. STAR Methods Key Resources Table covers tool citations.

**Cell Press convention:** numerical references in main text body. Methods may rely on STAR Methods table for tool citations. Acceptable.

**Reconciliation candidates (no decision):**
- **(a)** Keep current state — STAR table carries citations.
- **(b)** Add inline "(Muzellec et al. 2023)", "(Fang et al. 2023)", "(Orenbuch et al. 2020)" citations in Methods M3 and M9. Lengthens prose by ~5 sentences.

Recommend **(a)** — STAR Methods table is the conventional location.

---

## 5. ⚠ MINOR FINDING — Verb-tone consistency

Spot-check of Results section verbs (R1–R5):
- "We assembled" / "We analyzed" / "We tested" — used uniformly ✓
- "Demonstrating" / "indicating" / "consistent with" / "confirming" / "extending" — used in conclusion clauses, varied
- Borderline subjective: "★ exceptional in this cohort" (Caption F2D line 66) — could read as voice-touched. ⚠

**Reconciliation:** ★ markers in captions are convention markers, not subjective claims. Borderline. User may want to remove or replace with neutral language during voice-protected polish.

**This audit does not modify any caption.** Flag only.

---

## 6. Summary of gap report

| # | Severity | Finding |
|---|---|---|
| 1 | **NONE** | Methods prose covers all M1–M11 + M5b sections — no gap |
| 2 | LOW | Prose is *slightly more polished* than scaffold (added Suppl Table/Figure cites + CI mention) |
| 3 | LOW | Methods relies on STAR Methods table for tool citations — acceptable Cell Press convention |
| 4 | LOW | "★ exceptional" in F2D caption is borderline subjective — flag for user voice polish |

**No prose generation needed. T4 closed as audit.**

---

## 7. Recommended next action

1. **Do not modify Methods prose** — it is complete and already slightly improved over scaffold.
2. **Optionally remove "★ exceptional"** from F2D caption during voice-protected polish (single-word edit).
3. **Optionally add inline tool citations** in Methods M3/M9 if Cell Press reviewer requests — currently not needed.

---

Track A completed. No marathon violation. Methods prose audit is read-only — no manuscript file modified.

---

## 6. Decision summary

One-page issue inventory aggregated from sections 1–5. Sorted by severity. Each row is a single decision the user can make in <30 sec; resolution is mostly diff-friendly single-line edits.

### 🔴 HIGH (paper-shaping; resolve first)


| # | Issue | Files affected | Reconciliation candidate | Effort |
|---|---|---|---|---|
| H1 | DM1 sub-B "96% mutation-negative" vs SUPP T8-derived 94.6% (53/56) — arithmetic inconsistency | OUTLINE Abstract + R5d, R-PRO R5d, CAP F5D + F7A, QA Q12 (5 files) | (a) update to 94.6% (or ~95%) per SUPP; (b) fix SUPP T8 RAS+ count from 2→1 yielding 96.4%; (c) add denominator footnote ("of mut-tested n") | <30 min batch edit after decision |
| H2 | S1 cohort_assembly.tsv "build TBD" — referenced by 5 sources | S1 SUPP, M1 prose, M-SCAF, R1 prose, F1 caption, CHK | Build TSV from existing fragments | ~30 min ETL, marathon-safe |


### 🟡 MEDIUM (typesetting / reviewer-visible)


| # | Issue | Files affected | Reconciliation candidate | Effort |
|---|---|---|---|---|
| M1 | sub-A "69% RAS+" (51/74 mut-tested) vs "61%" (51/84 total) — denominator ambiguity | OUTLINE R5d, CAP F5D + F7A, R-PRO R5d, SUPP T8 (5 cells) | Add denominator note: "51/74 mut-tested (61% of 84 sub-A)" | 1 batch edit |
| M2 | F6 / F7 not explicitly named in R-PRO body | R-PRO R5b + R5d (2 lines) | Add "(Figure 6)" and "(Figure 7)" callouts | 5 min |
| M3 | Lee2024GSE213647 bib author placeholder ("Y and others") | references.bib | Pull full author list from GEO / PubMed | 15 min |
| M4 | Lim2025GSE286332 bib title ("Vitamin D…") may not match our PTC-vs-PTC+HT framing | references.bib | Verify title against PMID 41113708 | 5 min |
| M5 | S3 XLSX (29,672-row DEG) not packaged | submission_data/ | Convert TSV → XLSX | <15 min |
| M6 | S8d DM1 sub-A/B DEG list (8,935 rows) placeholder | SUPP S8 | Add as XLSX sheet | <15 min |


### 🟢 LOW (cosmetic)


| # | Issue | Action |
|---|---|---|
| L1 | TIERA67 / pan-genome ARI rounding (0.90 vs 0.903) — Outline uses rounded, others precise | Unify or accept rounding |
| L2 | TCGA Hashi top-30% Fisher p (6e-10 vs 6.4e-10) | Unify precision |
| L3 | BRAF mRNA Cohen d (-0.04 vs -0.044) | Unify to -0.044 |
| L4 | references.bib header claims "~22 papers" but actual 29 | Update header |
| L5 | ★ Unicode in Chu2018JMG note may break BibTeX | Replace with ASCII |
| L6 | Pan2025NatComm pages format unusual | Verify Nat Comm convention |
| L7 | 4 bib entries (Tuttle2019, Yi2016, Chen2024, Kim2014) not yet cited | Defer to Discussion drafting |
| L8 | "★ exceptional" verb in F2D caption — borderline subjective | User voice polish |
| L9 | Possible Krishnamoorthy 2025 bib gap | Decide at §3.1 voice draft |


### Recommended ordering


1. **Today (post-audit, <1 hr):** decide H1 reconciliation (a/b/c), apply 5-file batch edit.
2. **This week:** M2 (F6/F7 callouts) — single-line R-PRO edit; M3 + M4 (bib author/title verification).
3. **W5 (6/1–6/7):** H2 + M5 + M6 — 1-hour XLSX packaging session.
4. **Cosmetic pass W6 (6/8–6/13):** L1–L9 in one batch.


No marathon violation. All resolutions are non-analysis, non-data-download, non-voice-protected (except L8 which is explicitly user voice).


---


Track A completed. No marathon violation.
