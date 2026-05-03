# Paper 1 — Numerical Consistency Audit (T1)

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
