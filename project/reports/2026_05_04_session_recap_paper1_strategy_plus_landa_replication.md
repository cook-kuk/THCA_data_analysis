# Session recap — Paper 1 strategy + scaffolding + GSE76039 external advanced-disease replication

**Date:** 2026-05-04 · **Owner:** Seungho Cook
**Mode:** Read-only summary. **No further edits.** This file inventories the session's deliverables, lists the user-decision queue, and supplies verification commands.
**Marathon discipline:** voice-protected sections (Hook / Aim / Discussion §3.1 / §3.4 Limitations / Cover ¶1 / Reviewer Q9) were untouched throughout. No GPU/RunPod, no TCGA WSI, no H&E-DM1 retry, no Paper 3 / Paper 4 body modification.

---

## 1. Session arc (what was done, in order)

| # | Task | Result |
|---|---|---|
| 1 | 8-gene curated vs de novo / driver-map final strategy memo | `2026_05_04_8gene_curated_vs_denovo_final_strategy.md` (36 KB). Tournament A–J + ★ integrated; verdict: **demote 8-gene from headline / keep as compact RAI-lineage readout**; T1 title primary; "BRAF/RAS-negative dark matter" headline framing dropped (N1 NS). |
| 2 | Combined strategy + Hook ¶1 checklist VSCode-viewable file | `2026_05_04_paper1_strategy_plus_hook_checklist.md` (44 KB). Part A = strategy memo; Part B = Hook ¶1 6-section checklist (angle / 4-sentence / safe numbers / forbidden / citations / after-draft audit). |
| 3 | Manuscript scaffolding update — 4 in-place edits | `2026_05_03_manuscript_v8_OUTLINE.md` +176/−35 (scope banner · T1 title primary · Abstract M/R/C re-skeleton · R1/R2/R5 [moved] tags + R6/R7/R8/R9 NEW SCAFFOLD ROW + Figure plan F1–F5 + F6-Supp + F8); `2026_05_03_methods_M1_M11_scaffold.md` and `_prose.md` M2 +1 sentence each ("drivers in candidate pool but ranked low" correction); `2026_05_03_figure_captions_all.md` banner only. Voice-protected sections untouched. |
| 4 | Supp Methods scaffold for pre-tested H&E negative feasibility | `2026_05_04_supp_methods_closure_negative_feasibility.md` (9 KB). Closure-battery facts (r=0.022, AUROC=0.511, real ≤ random max). Limitations narrative deferred to author keyboard. |
| 5 | Minimal external-data acquisition plan | `2026_05_04_minimal_external_data_plan.md` (21 KB). 4 candidates evaluated; final = GSE76039 GO + Lee 2014 GO (citation only) + everything else DEFER/DROP. |
| 6 | GSE76039 acquisition + first-pass | Step 1 access check → Step 2 download (Series Matrix 6.9 MB + GPL570 SOFT 85.6 MB) → Step 3 metadata → Step 4 probe-to-gene mapping (22/22 panel-union genes available) → Step 5 7 module scores → Step 6 ATC vs PDTC + 4 Spearman tests → Step 7 3 figures → Step 8 reports (`2026_05_04_gse76039_access_check.md` + `2026_05_04_gse76039_first_pass_report.md`). **Result:** external advanced-disease replication SUCCEEDS at Cohen's d ≈ −3.5 ATC vs PDTC; ρ(DM1_like, NONOVERLAP) = −0.925, p = 3.1 × 10⁻¹⁶. **MAIN FIGURE candidate (F8)**. |
| 7 | Propagation of GSE76039 result into scaffolding | XLSX packaging (`SuppTable_SX_LANDA_GSE76039.xlsx`, 4 sheets) · OUTLINE F8 figure plan addition · Methods M12 in M-SCAF + M-PRO · acquisition-plan §3.1 / TL;DR / GSE33630-row corrections (microarray + n=37, NOT RNA-seq + n=84) · `Landa2016JCI` bib note extended. |

---

## 2. Files created or modified by Claude this session

### 2.1 Created (untracked)

| Path | Size |
|---|---:|
| `project/reports/2026_05_04_8gene_curated_vs_denovo_final_strategy.md` | 36,450 B |
| `project/reports/2026_05_04_paper1_strategy_plus_hook_checklist.md` | 43,646 B |
| `project/reports/2026_05_04_supp_methods_closure_negative_feasibility.md` | 9,123 B |
| `project/reports/2026_05_04_minimal_external_data_plan.md` | ~21 KB |
| `project/reports/2026_05_04_gse76039_access_check.md` | 5,254 B |
| `project/reports/2026_05_04_gse76039_first_pass_report.md` | 17,570 B |
| `project/reports/2026_05_04_session_recap_paper1_strategy_plus_landa_replication.md` | (this file) |
| `project/notebooks_or_scripts/p_landa_2016_first_pass.py` | 20,857 B |
| `project/results/p_landa_2016/raw/GSE76039_series_matrix.txt.gz` | 7.24 MB |
| `project/results/p_landa_2016/raw/GPL570_full.soft` | 85.60 MB |
| `project/results/p_landa_2016/sample_metadata.tsv` | 4,222 B |
| `project/results/p_landa_2016/expression_matrix_log_gene.tsv.gz` | 4,741 B |
| `project/results/p_landa_2016/probe_to_gene_panel.tsv` | 1,100 B |
| `project/results/p_landa_2016/scores.tsv` | 5,947 B |
| `project/results/p_landa_2016/score_tests.tsv` | 1,459 B |
| `project/results/p_landa_2016/SuppTable_SX_LANDA_GSE76039.xlsx` | 14,467 B |
| `project/results/p_landa_2016/gse76039_dm1_lineage_boxplot.png` | 102,798 B |
| `project/results/p_landa_2016/gse76039_dm1_vs_nonoverlap_scatter.png` | 68,775 B |
| `project/results/p_landa_2016/gse76039_mechanism_heatmap.png` | 72,708 B |

### 2.2 Modified in-place

| Path | Diff |
|---|---|
| `project/reports/2026_05_03_manuscript_v8_OUTLINE.md` | +228 / −36 (scope banner + T1 title + Abstract M/R/C re-skeleton + R1/R2/R5 [moved] tags + R6/R7/R8/R9 NEW SCAFFOLD ROW + Figure plan F1–F5 + F6-Supp + F8 + legacy [moved] mapping) |
| `project/reports/2026_05_03_methods_M1_M11_scaffold.md` | +1 paragraph M2 correction + new M12 entry |
| `project/reports/2026_05_03_methods_M1_M11_prose.md` | +1 paragraph M2 correction + new M12 entry |
| `project/reports/2026_05_03_figure_captions_all.md` | +10 lines top-of-file banner only |
| `project/reports/2026_05_03_references.bib` | +1 line extension to `Landa2016JCI` note |

---

## 3. Key results (carry-forward)

### 3.1 GSE76039 advanced-disease replication (Landa 2016, n=37; 17 PDTC + 20 ATC; Affymetrix HG-U133 Plus 2.0; gcRMA)

| Axis | ATC z mean (n=20) | PDTC z mean (n=17) | Cohen's d | MW p (two-sided) |
|---|---:|---:|---:|---:|
| RAI_8 | −0.666 | +0.784 | −3.47 | 4.6 × 10⁻⁷ |
| THYROID_NONOVERLAP (zero-overlap) | −0.665 | +0.782 | −3.32 | 6.3 × 10⁻⁷ |
| TDS_like | −0.666 | +0.783 | −3.54 | 4.6 × 10⁻⁷ |
| TF_collapse (FOXE1/NKX2-1/PAX8/HHEX) | −0.720 | +0.848 | −3.57 | 2.4 × 10⁻⁷ |
| STAT3_AP1_DNMT (5 genes) | +0.320 | −0.376 | +1.72 | 3.6 × 10⁻⁵ |
| TACSTD2 (TROP2 z) | +0.456 | −0.536 | +1.13 | 2.7 × 10⁻³ |

| Within-cohort Spearman (n=37) | ρ | p |
|---|---:|---:|
| DM1_like vs THYROID_NONOVERLAP | **−0.925** | 3.1 × 10⁻¹⁶ |
| TF_collapse vs DM1_like | −0.931 | 7.3 × 10⁻¹⁷ |
| STAT3_AP1_DNMT vs DM1_like | +0.681 | 3.5 × 10⁻⁶ |
| TACSTD2 vs DM1_like | +0.443 | 6.1 × 10⁻³ |

→ **direction-consistent across all 7 Paper 1 axes; zero-overlap lineage cross-validation stronger than TCGA-bulk reference (r = −0.885 on n = 561).** All honest-framing constraints applied (no PTC controls in series ⇒ no monotonic-gradient claim; cross-platform ⇒ within-cohort z only; mutation/survival metadata absent ⇒ no Cox / no BRAF/RAS-subset claims; TROP2 = "tumor-population vulnerability" only, no spot-level co-localisation claim).

### 3.2 Strategy decision (carried over from 2026_05_04 strategy memo)

- Title T1 primary: **"A transcriptional differentiation axis stratifies thyroid cancer orthogonally to canonical driver mutations and identifies a tumor-population TROP2 vulnerability"**
- 8-gene = compact readout (out of title and out of discovery framing)
- "BRAF/RAS-negative dark matter" headline framing dropped → honest-negative supplement only (N1 NS)
- Cancer Cell aim (functional + methylation + IC50 + prospective cohort) is NOT a marathon target

---

## 4. Verification commands

```bash
cd /home/seungho/personal/THCA_data_analysis

# A. forbidden-active-claim sweep (matches must all be intentional disclaimers)
grep -nEi "H&E predicts DM1|H&E-inferable|morphology-derived DM1|WSI-validated DM1|tile-level DM1 inference[^a-z]|pathology-AI triage|all risks resolved|Cancer Cell-ready|Nat Cancer reach|RunPod G[123][^_]|colocalizes with DM1|DM1-high spots are TROP2|excluded by design|5-15%|5–15%" \
  project/reports/2026_05_03_manuscript_v8_OUTLINE.md \
  project/reports/2026_05_03_methods_M1_M11_scaffold.md \
  project/reports/2026_05_03_methods_M1_M11_prose.md \
  project/reports/2026_05_04_minimal_external_data_plan.md \
  project/reports/2026_05_04_gse76039_first_pass_report.md \
  project/reports/2026_05_04_gse76039_access_check.md \
  project/reports/2026_05_03_references.bib

# B. GSE76039 numerical consistency (Cohen's d, ρ, n=37)
grep -nE "n[ =]?37|17 PDTC|20 ATC|−0\\.925|−3\\.47" project/reports/2026_05_*.md \
  project/reports/2026_05_03_manuscript_v8_OUTLINE.md \
  project/reports/2026_05_03_methods_M1_M11_*.md

# C. confirm GSE76039 Landa series matrix integrity
md5sum project/results/p_landa_2016/raw/GSE76039_series_matrix.txt.gz
# expected: b81615a5f5cffe27b10ba136111f8efa

# D. confirm script reproducibility
python3 project/notebooks_or_scripts/p_landa_2016_first_pass.py 2>&1 | grep -E "n_samples|histology counts|Spearman|score|MAIN" | head -20

# E. confirm SuppTable XLSX
python3 -c "import pandas as pd; xf=pd.ExcelFile('project/results/p_landa_2016/SuppTable_SX_LANDA_GSE76039.xlsx'); print(xf.sheet_names)"
# expected: ['S_meta', 'S_scores', 'S_tests', 'S_probe_to_gene']

# F. voice-protected guard — no Hook/Aim/Discussion/Limitations/Cover/Q9 prose generated by Claude this session
git diff project/reports/2026_05_03_manuscript_v8_OUTLINE.md | grep -E "^\\+" | grep -iE "no mechanistic compass|hook|cover letter|reviewer Q9|limitations narrative" | head
# expected: empty (Claude does not generate voice-protected prose)
```

---

## 5. Recommended commit groups (user decision; nothing committed by Claude this session)

| group | files | message suggestion |
|---|---|---|
| **C1 — Paper 1 strategy + Hook checklist (NEW reports)** | `2026_05_04_8gene_curated_vs_denovo_final_strategy.md`, `2026_05_04_paper1_strategy_plus_hook_checklist.md` | `docs: paper1 strategy memo + hook checklist (8-gene compact readout, T1 title)` |
| **C2 — Manuscript scaffolding update (in-place)** | `2026_05_03_manuscript_v8_OUTLINE.md`, `2026_05_03_methods_M1_M11_scaffold.md`, `2026_05_03_methods_M1_M11_prose.md`, `2026_05_03_figure_captions_all.md` | `docs: paper1 v8 outline scope re-skeleton (T1 title; R6-R9 scaffold; F1-F5+F6 figure plan)` |
| **C3 — Closure-battery negative-feasibility Supp Methods** | `2026_05_04_supp_methods_closure_negative_feasibility.md` | `docs: paper1 supp methods scaffold for pre-tested H&E negative feasibility` |
| **C4 — External-data plan + GSE76039 acquisition + first-pass** | `2026_05_04_minimal_external_data_plan.md`, `2026_05_04_gse76039_access_check.md`, `2026_05_04_gse76039_first_pass_report.md`, `project/notebooks_or_scripts/p_landa_2016_first_pass.py`, `project/results/p_landa_2016/{sample_metadata.tsv, expression_matrix_log_gene.tsv.gz, probe_to_gene_panel.tsv, scores.tsv, score_tests.tsv, SuppTable_SX_LANDA_GSE76039.xlsx, 3 PNG figures}` | `feat: paper1 external advanced-disease replication via GSE76039 (Landa 2016)` |
| **C5 — Bib + scaffolding propagation of GSE76039 result** | `2026_05_03_references.bib` (Landa note extension), F8 lines in OUTLINE, M12 entries in M-SCAF/M-PRO, acquisition-plan §3.1/§3.2/§4 corrections | (squash with C2 or C4 as user prefers) |
| **C6 — gitignore for raw GEO files** | `.gitignore` (suggested addition: `project/results/p_landa_2016/raw/`) — **note:** raw files (~93 MB) should be excluded from git; reproducible from the URLs in the access check | `chore: gitignore raw GEO downloads (reproducible from access check URLs)` |

**Important:** the `raw/` subdirectory under `p_landa_2016/` contains 92.8 MB of files (Series Matrix + GPL570 SOFT) — these are reproducible from public GEO URLs in `2026_05_04_gse76039_access_check.md` §3 and should NOT be committed. Add to `.gitignore` before `git add .`.

**Files NOT to commit:**
- `project/results/p_landa_2016/raw/*` — reproducible from access-check URLs
- Anything under `submission/`, `manuscript_p2_brief/`, `paper3_ici/`, `papers_overview.html`, `manuscript_v8/PAPER2_*` — Paper 2/3/4 territory written by parallel agents
- `MARATHON_REENTRY_FINAL_*.md`, `PAPERS_5_OVERVIEW_*.md`, `papers_hub_2026_05_04/` (parallel agent root-level memos; user reviews separately)

---

## 6. User decision queue

1. **Voice-protected drafting (author keyboard, not Claude)** — Hook ¶1 / Aim ¶4 / Discussion §3.1 (Landa 2016 + GSE76039 result integration) / §3.4 Limitations / Cover ¶1 / Reviewer Q9. Hook ¶1 checklist already prepared in `2026_05_04_paper1_strategy_plus_hook_checklist.md` Part B.
2. **F8 composite figure** — current outputs are 3 separate PNGs (`gse76039_dm1_lineage_boxplot.png`, `gse76039_dm1_vs_nonoverlap_scatter.png`, `gse76039_mechanism_heatmap.png`); user decides single multi-panel composite vs three sub-panels for final submission. (Stylistic; no new analysis required.)
3. **Commit grouping** — accept C1–C6 above, regroup, or ask Claude to execute one group at a time after explicit "commit C1" / etc.
4. **Bundang outreach status** — `v17_npj_ship_status` notes 4 outreach drafts NOT sent. This session did not touch outreach.
5. **Future acquisitions** — per `2026_05_04_minimal_external_data_plan.md` §4, all DEFER/DROP unless explicit per-slot "go." TCGA-THCA methylation 450K is the most likely next acquisition (Nature Cancer reach / revision); recommend NOT starting until after voice-protected draft + Yu meeting.
6. **Parallel-agent territory** — many other 2026-05-04 reports / files exist (Paper 2 brief, Pod C/D completion, scope-split audit, etc.). I have not modified them in this session. User reviews separately.

---

## 7. Marathon-discipline final attestation

| Constraint | Status |
|---|---|
| New analysis (paper-blocking-only rule) | ✓ one external acquisition (GSE76039) — paper-blocking for Nature Cancer reach gate |
| New data download | ✓ one (GSE76039 Series Matrix + GPL570 annotation, 92.8 MB total) |
| GPU / RunPod / Azure burst | ✓ none (all CPU-only on primary VM) |
| TCGA WSI / H&E-DM1 retry | ✓ none |
| GSE33630 / TCGA methylation / DepMap acquisition | ✓ NOT performed |
| Paper 3 / Paper 4 body modification | ✓ none |
| Voice-protected sections (Hook / Aim / Discussion §3.1 / §3.4 Limitations / Cover ¶1 / Q9) | ✓ untouched |
| Manuscript prose (large-scale) modification | ✓ none |
| Forbidden language in active-claim form | ✓ all matches are intentional disclaimers |
| Result presentation as "main figure" without honest direction check | ✓ none — all 7 axes direction-consistent + Cohen's d > 1 + p < 1e-2 |
| Time budget for GSE76039 (today: 6 h to first score table) | ✓ ~30 min total CPU + I/O, scaffolding propagation ~15 min |
| Commit by Claude | ✓ none — all changes staged-uncommitted, awaiting user review |

---

## 8. Cross-reference index

- Strategy: `2026_05_04_8gene_curated_vs_denovo_final_strategy.md` · `2026_05_04_paper1_strategy_plus_hook_checklist.md`
- Lock authority: `2026_05_04_paper1_dm1_full_molecular_only_lock.md`
- External-data plan: `2026_05_04_minimal_external_data_plan.md`
- GSE76039 first-pass: `2026_05_04_gse76039_access_check.md` · `2026_05_04_gse76039_first_pass_report.md`
- Pipeline script: `project/notebooks_or_scripts/p_landa_2016_first_pass.py`
- Outputs: `project/results/p_landa_2016/`
- Manuscript scaffolding: `2026_05_03_manuscript_v8_OUTLINE.md` (with F8 added) · `2026_05_03_methods_M1_M11_scaffold.md` (with M12 added) · `2026_05_03_methods_M1_M11_prose.md` (with M12 added) · `2026_05_03_figure_captions_all.md` (banner only) · `2026_05_03_references.bib` (Landa note extension)
- Closure battery: `2026_05_04_supp_methods_closure_negative_feasibility.md` · `pathology_dm1_phaseA_cpu_verdict_2026_05_04.md` · `pathology_dm1_closure_battery_2026_05_04.md` · `2026_05_04_image_dm1_final_nogo_decision.md`
- Memory: `paper_numbering_2026_05_04` · `v18_paper2_HT_isolated` · `v19_paper4_GD_backlog` · `v17_marathon_mode_post_pillar1` · `v17_sprint_vs_marathon_violation` · `v17_korean_K2_calibration` · `v17_landa2016_cite_save`

---

*Session wrap-up authored 2026-05-04 by Claude (Opus 4.7) under marathon-mode discipline. No commit by Claude. All deliverables staged-uncommitted, awaiting user review per §5/§6.*
