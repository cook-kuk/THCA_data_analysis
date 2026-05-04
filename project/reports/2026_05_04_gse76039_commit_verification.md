# GSE76039 — commit-safe verification (STEPS 1–5)

**Date:** 2026-05-04 · **Owner:** Seungho Cook
**Mode:** Read-only verification + `.gitignore` cleanup. **No commit by Claude.** Awaiting per-group user approval.
**Marathon discipline:** voice-protected sections (Hook / Aim / Discussion §3.1 / §3.4 Limitations / Cover ¶1 / Q9) untouched. No new analysis, no new download, no GPU/RunPod, no GSE33630 / TCGA-methylation acquisition, no Paper 3/4 body modification, no H&E-DM1 retry, no TCGA WSI.
**Verdict:** Commit-safe. Raw files (~89 MB) excluded by new `.gitignore` rule; processed files (~292 KB) commit-ready.

---

## 1. Raw files excluded from git (STEP 2)

`.gitignore` updated (single 3-line append):

```
# Landa 2016 GSE76039 raw GEO downloads (reproducible from public GEO URLs in
# project/reports/2026_05_04_gse76039_access_check.md §2-§3; ~89 MB)
project/results/p_landa_2016/raw/
```

Verification:

```
$ git check-ignore -v project/results/p_landa_2016/raw/GSE76039_series_matrix.txt.gz
.gitignore:74:project/results/p_landa_2016/raw/  project/results/p_landa_2016/raw/GSE76039_series_matrix.txt.gz
$ git check-ignore -v project/results/p_landa_2016/raw/GPL570_full.soft
.gitignore:74:project/results/p_landa_2016/raw/  project/results/p_landa_2016/raw/GPL570_full.soft
$ git check-ignore project/results/p_landa_2016/scores.tsv  # processed file
(no output — NOT ignored, OK to commit)
```

Excluded raw files (reproducible from URLs in `2026_05_04_gse76039_access_check.md` §2/§3):

| File | Size | md5 (recorded for reproducibility) | URL (also in access-check) |
|---|---:|---|---|
| `GSE76039_series_matrix.txt.gz` | 7.0 MB | `b81615a5f5cffe27b10ba136111f8efa` | `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE76nnn/GSE76039/matrix/GSE76039_series_matrix.txt.gz` |
| `GPL570_full.soft` | 82 MB | `ab5b044c4795d307e9201cda9a0e6a4f` | `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?targ=self&form=text&view=full&acc=GPL570` |

**Total excluded:** ~89 MB (no committed payload).

---

## 2. Processed files to commit (STEP 2)

| Path | Size | Type |
|---|---:|---|
| `project/results/p_landa_2016/sample_metadata.tsv` | 4.2 KB | TSV |
| `project/results/p_landa_2016/expression_matrix_log_gene.tsv.gz` | 4.7 KB | gzip TSV (panel-only gene matrix) |
| `project/results/p_landa_2016/probe_to_gene_panel.tsv` | 1.1 KB | TSV |
| `project/results/p_landa_2016/scores.tsv` | 5.9 KB | TSV |
| `project/results/p_landa_2016/score_tests.tsv` | 1.5 KB | TSV |
| `project/results/p_landa_2016/SuppTable_SX_LANDA_GSE76039.xlsx` | 14.5 KB | XLSX (4 sheets: S_meta, S_scores, S_tests, S_probe_to_gene) |
| `project/results/p_landa_2016/gse76039_dm1_lineage_boxplot.png` | 101 KB | PNG (4-panel ATC vs PDTC) |
| `project/results/p_landa_2016/gse76039_dm1_vs_nonoverlap_scatter.png` | 68 KB | PNG (DM1_like vs NONOVERLAP scatter) |
| `project/results/p_landa_2016/gse76039_mechanism_heatmap.png` | 72 KB | PNG (10-gene mechanism heatmap) |

**Total committed (results):** 292 KB.

**Reproducibility script:** `project/notebooks_or_scripts/p_landa_2016_first_pass.py` (20.4 KB).

**Reports:** `project/reports/2026_05_04_gse76039_access_check.md` (5.3 KB) + `project/reports/2026_05_04_gse76039_first_pass_report.md` (17.6 KB) + `project/reports/2026_05_04_minimal_external_data_plan.md` (~21 KB).

---

## 3. Numerical consistency check (STEP 4 §3)

Verified across all GSE76039-touching files (access check, first-pass report, acquisition plan, OUTLINE F8, methods M12 in M-SCAF + M-PRO, session recap, references.bib `Landa2016JCI` note):

| Quantity | Value | All files agree? |
|---|---|---|
| Cohort total n | 37 | ✓ |
| PDTC count | 17 | ✓ |
| ATC count | 20 | ✓ |
| PTC controls in series | 0 | ✓ (explicit "no PTC controls") |
| Normal samples in series | 0 | ✓ |
| Platform | GPL570 / Affymetrix HG-U133 Plus 2.0 | ✓ |
| Pre-processing | gcRMA log₂ | ✓ |
| RAI_8 ATC vs PDTC Cohen's d | −3.47 | ✓ |
| RAI_8 MW p (two-sided) | 4.6 × 10⁻⁷ | ✓ |
| THYROID_NONOVERLAP Cohen's d | −3.32 | ✓ |
| THYROID_NONOVERLAP MW p | 6.3 × 10⁻⁷ | ✓ |
| TF_collapse Cohen's d (largest) | −3.57 | ✓ |
| STAT3_AP1_DNMT Cohen's d | +1.72 | ✓ |
| TACSTD2 Cohen's d (smallest) | +1.13 | ✓ |
| ρ(DM1_like, NONOVERLAP) on n=37 | −0.925 | ✓ |
| p for that ρ | 3.1 × 10⁻¹⁶ | ✓ |
| ρ(TF_collapse, DM1_like) | −0.931 | ✓ |
| ρ(STAT3_AP1_DNMT, DM1_like) | +0.681 | ✓ |
| ρ(TACSTD2, DM1_like) | +0.443 | ✓ |
| Reference (TCGA bulk r DM1 vs NONOVERLAP, n=561) | −0.885 | ✓ |
| Earlier "RNA-seq (n=84)" claim | corrected to microarray n=37 in all places | ✓ (TL;DR table, §3.1 banner, §3.2 GSE33630 row, §4 don't-acquire row, §recap) |

**Status:** All numerical claims internally consistent. No drift.

---

## 4. Overclaim sweep (STEP 3 + STEP 4 §4)

User-listed forbidden phrases swept across nine files (`access_check.md`, `first_pass_report.md`, `minimal_external_data_plan.md`, `session_recap`, `OUTLINE`, `M-SCAF`, `M-PRO`, `references.bib`, `p_landa_2016_first_pass.py`):

| Phrase | Active-claim matches | Disclaimer/guard matches | Status |
|---|---:|---:|---|
| "proves progression" | **0** | 0 | ✓ clean |
| "PTC to PDTC to ATC" / "PTC → PDTC → ATC" | **0** | 6 (all framed as "do NOT claim ... monotonic 3-stage gradient" / "limits ladder claim") | ✓ clean |
| "monotonic" | **0** | 5 (all guard form: "NOT ... monotonic" / "monotonic-gradient claims are not made") | ✓ clean |
| "Landa proves" | **0** | 0 | ✓ clean |
| "external survival" | **0** | 1 (immediately disclaimed: "GSE76039 partially closes this gate — molecular validation, no survival") | ✓ clean |
| "BRAF/RAS subset" / "BRAF-/RAS- subset" | **0** for GSE76039 | 5 (all about TCGA-internal N1 honest negative or explicit "cannot test ... in this cohort" guard) | ✓ clean |
| "fusion" | **0** | 3 (all DROP / disclaimer: MSK-IMPACT DROP, "Fusion-map figure (no fusion calls available)") | ✓ clean |
| "methylation" | **0** | 8 (all "absent" / "DROP" / "DEFER" / Cancer-Cell-aim future requirement) | ✓ clean |
| "proves" / "establishes" | **0** | 2 (both guard: "supports ... NOT proves" / "approves" was a false-positive from "approves" ≠ "proves") | ✓ clean |
| "all risks resolved" | **0** | 0 | ✓ clean |
| "Cancer Cell-ready" | **0** | 0 | ✓ clean |
| "Nat Cancer reach" / "Nature Cancer reach" | **0 in manuscript-scaffolding files (OUTLINE, M-SCAF, M-PRO)** | 8 in internal planning reports (strategy memo cross-references, "what would be required for Nature Cancer reach"); these are internal-only and must NOT leak into submission text | ✓ clean for manuscript scope; planning-internal use OK |
| "H&E predicts DM1" / "H&E-inferable" / "morphology-derived DM1" / "WSI-validated DM1" | **0** | 2 (forbidden-language list; explicit "in any active-claim form") | ✓ clean |
| "colocalizes with DM1" / "DM1-high spots are TROP2" | **0** | 2 (explicit guards, including TROP2 framing in first-pass report §6.3) | ✓ clean |

**Sweep verdict:** zero active-claim violations. All flagged phrases occur only as disclaimers, guards, exclusion lists, or explicit forbidden-wording references. The internal planning files use "Nature Cancer reach" only in internal venue-ladder discussion (strategy memo §12.6 / acquisition plan §3.3); none of this language is in the manuscript-scaffolding files (OUTLINE / M-SCAF / M-PRO / references.bib).

---

## 5. Recommended placement of GSE76039 result

| Option | Rationale | Recommendation |
|---|---|---|
| **A. Main Figure F8 (3-panel composite)** | Effect sizes are extreme (Cohen's d ≈ −3.5 ATC vs PDTC); zero-overlap lineage cross-validation ρ = −0.925 is stronger than TCGA-bulk reference; mechanism panel direction-consistent across TF / STAT3-DNMT / TROP2; closes strategy-memo §12.6 Nature-Cancer-reach gate (partially: molecular only, no survival). Available source PNGs: `gse76039_dm1_lineage_boxplot.png`, `gse76039_dm1_vs_nonoverlap_scatter.png`, `gse76039_mechanism_heatmap.png`. | **GO if main-figure space allows** (Paper 1 currently has F1–F5 + F6-Supp per OUTLINE); F8 is a clean addition that does not cannibalise other panels. |
| B. Supplementary figure | Conservative fallback if author/Yu prefer keeping main-figure count at 5–6. The same three panels move to Suppl Figure SX_landa. | Acceptable fallback. Author/Yu decision. |
| C. Internal-only / drop | NOT recommended. Result is direction-consistent, effect-size-large, externally validated, manuscript-safe. Burying it would forfeit the only external advanced-disease replication available without re-acquisition. | NOT recommended. |

**My recommendation:** **A** (Main Figure F8) — composite of three current PNGs, with caption framing per §6 below. If space-constrained, **B** (Supplementary). The OUTLINE F8 entry already accommodates either choice.

---

## 6. Exact safe wording (STEP 4 §6)

**For Methods M12** (already in `2026_05_03_methods_M1_M11_scaffold.md` and `_prose.md`): use the existing language; numerical consistency confirmed §3 above.

**For Results section (whenever author drafts it — voice-protected; this is a wording reference, not generated prose):**

| Concept | Approved wording (exact) |
|---|---|
| Result framing | "**independent advanced-disease replication** of the DM1_like / RAI-lineage axis in the Landa 2016 GSE76039 cohort (n = 37; 17 PDTC + 20 ATC; Affymetrix HG-U133 Plus 2.0)" |
| Within-cohort effect | "**distinguishes ATC from PDTC** at extreme effect size (RAI_8 Cohen's d = −3.47, MW p = 4.6 × 10⁻⁷; THYROID_NONOVERLAP Cohen's d = −3.32, MW p = 6.3 × 10⁻⁷)" |
| Lineage cross-validation | "**recapitulates the non-overlapping thyroid-lineage anti-correlation** observed in TCGA bulk (Spearman ρ = −0.925 in GSE76039, n = 37, p = 3.1 × 10⁻¹⁶, vs r = −0.885 in TCGA, n = 561)" |
| Mechanism replication | "consistent with the lineage-TF collapse + STAT3 / DNMT activation framework (TF_collapse Cohen's d = −3.57; STAT3_AP1_DNMT Cohen's d = +1.72; both p < 1 × 10⁻⁴)" |
| TROP2 replication | "tumor-population TROP2 / TACSTD2 elevation in the more dedifferentiated (ATC) end of advanced disease (Cohen's d = +1.13, p = 2.7 × 10⁻³)" |
| Cross-platform note | "within-cohort z-score normalisation only; the TCGA-trained absolute-form classifier was not transferred across the microarray–RNA-seq platform boundary" |
| Overall verdict | "**supportive external advanced-disease replication; not primary proof of stage progression**" |

**For Figure caption (F8):** use the panel descriptions from `2026_05_03_manuscript_v8_OUTLINE.md` Figures §F8 (already present) — they are fact-only and consistent with the wording above.

**For Discussion §3.1 (voice-protected, author keyboard):** the result becomes a single number to integrate into the existing Landa 2016 framing — no Claude prose generated. Suggested integration spot: after the existing sentence introducing Landa 2016 PDTC/ATC dedifferentiation precedent, add a result-tying clause along the lines of "and we replicate the dedifferentiation-axis signal on Landa's transcriptomic subset (GSE76039 n = 37) at extreme effect size" — **author keyboard** for actual phrasing.

---

## 7. Forbidden wording (must NOT appear in manuscript text or figure captions)

| Phrase | Why forbidden |
|---|---|
| "GSE76039 proves PTC → PDTC → ATC progression" / "monotonic gradient" / "stage ladder" | No PTC controls in series; only ATC vs PDTC contrast within cohort |
| "Landa cohort confirms our survival association" / "external survival validation" / "external Cox" | No survival metadata in GEO record |
| "Replicates BRAF/RAS-negative dark-matter result" / "validates Cox in BRAF-/RAS- subset" | No mutation status in GEO record; N1 BRAF-/RAS- subset stays as TCGA-internal honest negative |
| "DM1-high spots are TROP2-high in advanced disease" / "spatial co-localisation in Landa cohort" | GSE76039 is bulk microarray; no spatial information; Q3 spot-level NEG remains the TCGA-spatial finding |
| "Cancer-Cell-ready" / "Nature Cancer reach achieved" | Internal venue-ladder language; never submission-facing |
| "GSE76039 RNA-seq" / "GSE76039 n=84" | Both incorrect (microarray, n=37). Earlier draft of acquisition plan had this error; corrected post-acquisition |
| "Fusion calls in Landa cohort" / "RET / NTRK / ALK fusion replication" | Not in GEO record; out of scope per strategy memo §8 |
| "Methylation-supported TF→DNMT mechanism" | Methylation not assayed in GSE76039; mechanism remains "consistent with," not "demonstrates" |
| "H&E-inferable from Landa data" / "image triage" | Closure battery NO-GO applies; GSE76039 is bulk transcriptomics |
| "Functionally validated in Landa cohort" | No functional perturbation; future work |

---

## 8. Commit groups proposed (STEP 5) — **awaiting user "commit Cx" approval; no commits by Claude**

| Group | Files | Size | Suggested message |
|---|---|---:|---|
| **C2 — manuscript scaffolding (Paper 1 v8 outline scope re-skeleton, GSE76039 propagated)** | `2026_05_03_manuscript_v8_OUTLINE.md` (M) · `2026_05_03_methods_M1_M11_scaffold.md` (M) · `2026_05_03_methods_M1_M11_prose.md` (M) · `2026_05_03_figure_captions_all.md` (M) · `2026_05_03_references.bib` (M, Landa note extended) | +228/−36 OUTLINE + ~+30 across others | `docs: paper1 v8 scope re-skeleton (T1 title; R6-R9 scaffold; F1-F5 + F8 figure plan; M2 driver-pool correction; M12 GSE76039)` |
| **C3 — closure-battery negative-feasibility Supp Methods scaffold** | `2026_05_04_supp_methods_closure_negative_feasibility.md` (NEW) | 9 KB | `docs: paper1 supp methods scaffold for pre-tested H&E negative feasibility` |
| **C4a — GSE76039 reports + script + .gitignore** | `2026_05_04_gse76039_access_check.md` (NEW) · `2026_05_04_gse76039_first_pass_report.md` (NEW) · `2026_05_04_minimal_external_data_plan.md` (NEW) · `2026_05_04_gse76039_commit_verification.md` (NEW, this file) · `2026_05_04_session_recap_paper1_strategy_plus_landa_replication.md` (NEW) · `project/notebooks_or_scripts/p_landa_2016_first_pass.py` (NEW) · `.gitignore` (M, +3 lines for raw/) | ~80 KB + 21 KB script | `feat: paper1 external advanced-disease replication via GSE76039 (Landa 2016) — reports + script + .gitignore` |
| **C4b — GSE76039 processed results + figures + supp table** | `project/results/p_landa_2016/sample_metadata.tsv` (NEW) · `expression_matrix_log_gene.tsv.gz` · `probe_to_gene_panel.tsv` · `scores.tsv` · `score_tests.tsv` · `SuppTable_SX_LANDA_GSE76039.xlsx` · `gse76039_dm1_lineage_boxplot.png` · `gse76039_dm1_vs_nonoverlap_scatter.png` · `gse76039_mechanism_heatmap.png` | 292 KB | `feat: paper1 GSE76039 processed scores + supp table + figures` |
| **C1 — strategy memo + Hook checklist (NEW reports)** | `2026_05_04_8gene_curated_vs_denovo_final_strategy.md` (NEW) · `2026_05_04_paper1_strategy_plus_hook_checklist.md` (NEW) | 80 KB | `docs: paper1 strategy memo + hook checklist (8-gene compact readout, T1 title, dark-matter headline framing dropped)` |

### Order recommendation

The cleanest commit order is **C1 → C2 → C3 → C4a → C4b** (strategy first, then scaffolding, then Supp Methods, then GSE76039 reports + script + .gitignore, then the GSE76039 result artifacts). This way the GSE76039 commits land on top of an OUTLINE that already references F8 / M12.

**Alternative:** if user prefers a single landed PR per concern, group as `(C1 + C2 + C3)` then `(C4a + C4b)` — also fine.

### Files explicitly NOT in any of the above groups

- Anything `M` (modified) under `project/manuscript_p2_brief/`, `project/submission/`, `project/three_papers_index.*`, `project/manuscript_v8/0[0-9]_*.md`, `project/manuscript_v8/PAPER2_*`, `project/manuscript_v8/_*.md`, `project/results/terminology_correction_2026_05_04/*` → **Paper 2 / Paper 4 territory written by parallel agents this session;** Claude did not modify them and recommends user reviews / commits them separately under their own scope.
- Anything `??` (untracked) at repo root (`MARATHON_REENTRY_FINAL_*.md`, `PAPERS_5_OVERVIEW_*.md`, `POD_C*.md`, `GIGATIME_STRATEGY_*.md`, `papers_hub_2026_05_04/`) → parallel-agent root-level memos; user reviews separately.
- Anything `??` under `project/manuscript_v8/{10,14}_*.md` → parallel-agent compiled drafts; user reviews separately.
- Anything `??` under `project/notebooks_or_scripts/{poll_and_recover_pod_D.sh, post_pod_*.py}` → image-DM1 closure / Pod B/C/D postprocess scripts (B-class deprecated angle, per `POD_CLEANUP_AND_PAPER1_LOCK_2026_05_04.md`); user decides whether to commit historical or leave untracked.
- `project/results/p_landa_2016/raw/*` → ignored by `.gitignore` (this commit-verification §1).

---

## 9. Verification commands the user can run

```bash
cd /home/seungho/personal/THCA_data_analysis

# A. confirm raw files ignored, processed files NOT ignored
git check-ignore project/results/p_landa_2016/raw/GSE76039_series_matrix.txt.gz
git check-ignore project/results/p_landa_2016/raw/GPL570_full.soft
# expected: each prints the .gitignore rule line
git check-ignore project/results/p_landa_2016/scores.tsv \
                  project/results/p_landa_2016/SuppTable_SX_LANDA_GSE76039.xlsx
# expected: empty (NOT ignored)

# B. forbidden-language sweep (must return only intentional disclaimers)
grep -niE "proves progression|Landa proves|external survival|BRAF/RAS subset|fusion[^-]|methylation" \
  project/reports/2026_05_04_gse76039_*.md \
  project/reports/2026_05_04_minimal_external_data_plan.md \
  project/reports/2026_05_03_manuscript_v8_OUTLINE.md \
  project/reports/2026_05_03_methods_M1_M11_*.md \
  project/reports/2026_05_03_references.bib

# C. numerical consistency
grep -nE "n[ =]?37|17 PDTC|20 ATC|−0\\.925|−3\\.47" \
  project/reports/2026_05_*.md \
  project/reports/2026_05_03_manuscript_v8_OUTLINE.md \
  project/reports/2026_05_03_methods_M1_M11_*.md

# D. md5 check on ignored raw files
md5sum project/results/p_landa_2016/raw/GSE76039_series_matrix.txt.gz
# expected: b81615a5f5cffe27b10ba136111f8efa
md5sum project/results/p_landa_2016/raw/GPL570_full.soft
# expected: ab5b044c4795d307e9201cda9a0e6a4f

# E. reproducibility smoke (re-running the script gives the same numbers)
python3 project/notebooks_or_scripts/p_landa_2016_first_pass.py 2>&1 | tail -20
# expected: Spearman -0.925, ATC vs PDTC d ≈ -3.47, etc.
```

---

## 10. Marathon-discipline final attestation

| Constraint | Status |
|---|---|
| New analysis | ✓ none in this verification turn |
| New data download | ✓ none |
| GPU / RunPod / Azure burst | ✓ none |
| TCGA WSI / H&E-DM1 retry | ✓ none |
| GSE33630 / TCGA methylation / DepMap | ✓ NOT acquired |
| Paper 3 / Paper 4 body modification | ✓ none |
| Voice-protected sections (Hook / Aim / Discussion §3.1 / §3.4 Limitations / Cover ¶1 / Q9) | ✓ untouched |
| Forbidden-active-claim violations | ✓ zero (sweep §4) |
| Numerical drift | ✓ none (consistency §3) |
| Raw GEO files committed to git | ✓ NOT (excluded by `.gitignore` §1) |
| Commit by Claude | ✓ none — all changes staged-uncommitted; awaiting per-group user "commit Cx" |

---

## 11. Cross-reference index

- Strategy: `2026_05_04_8gene_curated_vs_denovo_final_strategy.md` · `2026_05_04_paper1_strategy_plus_hook_checklist.md`
- Lock authority: `2026_05_04_paper1_dm1_full_molecular_only_lock.md`
- External-data plan: `2026_05_04_minimal_external_data_plan.md`
- GSE76039 access check: `2026_05_04_gse76039_access_check.md`
- GSE76039 first-pass: `2026_05_04_gse76039_first_pass_report.md`
- This commit verification: `2026_05_04_gse76039_commit_verification.md`
- Session recap: `2026_05_04_session_recap_paper1_strategy_plus_landa_replication.md`
- Pipeline script: `project/notebooks_or_scripts/p_landa_2016_first_pass.py`
- Outputs: `project/results/p_landa_2016/{sample_metadata.tsv, expression_matrix_log_gene.tsv.gz, probe_to_gene_panel.tsv, scores.tsv, score_tests.tsv, SuppTable_SX_LANDA_GSE76039.xlsx, 3 PNG figures}`
- Manuscript scaffolding: `2026_05_03_manuscript_v8_OUTLINE.md` (F8) · `2026_05_03_methods_M1_M11_scaffold.md` (M12) · `2026_05_03_methods_M1_M11_prose.md` (M12) · `2026_05_03_figure_captions_all.md` (banner) · `2026_05_03_references.bib` (Landa note)
- Closure battery: `2026_05_04_supp_methods_closure_negative_feasibility.md` · `pathology_dm1_phaseA_cpu_verdict_2026_05_04.md` · `pathology_dm1_closure_battery_2026_05_04.md` · `2026_05_04_image_dm1_final_nogo_decision.md`
- Memory: `paper_numbering_2026_05_04` · `v18_paper2_HT_isolated` · `v19_paper4_GD_backlog` · `v17_marathon_mode_post_pillar1` · `v17_sprint_vs_marathon_violation` · `v17_korean_K2_calibration` · `v17_landa2016_cite_save`

---

*Verification authored 2026-05-04 by Claude (Opus 4.7) under marathon-mode discipline. Read-only audit + `.gitignore` cleanup only. No commit. Awaiting user per-group "commit Cx" approval.*
