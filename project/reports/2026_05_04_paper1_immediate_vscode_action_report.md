# Paper 1 — Immediate VS Code Action Report
**Date:** 2026-05-04
**Mode:** ship-prep only — no new analysis, no new data, no Hook prose by Claude
**Owner:** Seungho Cook
**Scope window:** today's local VS Code session

---

## 1. Current Paper 1 identity (canonical)

- **Paper 1 = active main paper.**
- **Topic:** driver-orthogonal transcriptional differentiation axis in PTC.
- **8-gene readout:** demoted from title/discovery framing → re-cast as a *compact RAI-lineage silencing readout*. It is an instrument, not the discovery.
- **Forbidden re-framings during marathon:**
  - "8-gene signature paper"
  - "fusion-driven, epigenetically silenced"
  - "promoter hypermethylation drives 8-gene"
  - any title that puts 8-gene first
- **Voice-protected sections (Claude must not draft prose):** Hook, Aim, Discussion §3.1, Limitations, Cover letter Para 1, Reviewer Q9.

---

## 2. GSE76039 (Landa 2016) — verification summary

Cohort: GSE76039, primary thyroid carcinomas, Affymetrix HG-U133 Plus 2.0 (GPL570), gcRMA-processed.

**Sample size (verified from `sample_metadata.tsv`):** n = 37 total
- ATC (anaplastic): 20
- PDTC (poorly differentiated): 17
- PTC: 0 (no PTC controls in series — ladder framing not available from this cohort alone)

**Probe → gene panel (verified from `probe_to_gene_panel.tsv`):** 21 genes mapped via best-probe selection across the RAI_8, THYROID_NONOVERLAP, TF_collapse, STAT3_AP1_DNMT_DNA-methylation, and TACSTD2 (TROP2) modules; alias `TITF1 → NKX2-1` honored.

**Key score results (verified from `score_tests.tsv`, n=37, ATC vs PDTC except Spearman):**

| Test | Cohen d | MW p (two-sided) | direction |
|---|---|---|---|
| RAI_8 score, ATC vs PDTC | **−3.469** | 4.56e−07 | RAI lineage suppressed in ATC |
| THYROID_NONOVERLAP score, ATC vs PDTC | **−3.317** | 6.27e−07 | non-overlap module also suppressed in ATC |
| TDS_like score, ATC vs PDTC | −3.538 | 4.56e−07 | combined 16-gene RAI/thyroid axis suppressed in ATC |
| TF_collapse score (FOXE1/NKX2-1/PAX8/HHEX), ATC vs PDTC | −3.566 | 2.40e−07 | thyroid TF program suppressed in ATC |
| STAT3_AP1_DNMT score, ATC vs PDTC | +1.720 | 3.64e−05 | inflammation/AP-1/DNMT axis up in ATC |
| TACSTD2 (TROP2) z, ATC vs PDTC | +1.131 | 2.68e−03 | TROP2 elevated at tumor (bulk) level in ATC |

**Key within-cohort Spearman correlations (n=37):**

| Pair | ρ | p |
|---|---|---|
| **DM1_like vs THYROID_NONOVERLAP** | **−0.9246** | **3.08e−16** |
| TF_collapse vs DM1_like | −0.9308 | 7.27e−17 |
| STAT3_AP1_DNMT vs DM1_like | +0.6809 | 3.52e−06 |
| TACSTD2 vs DM1_like | +0.4426 | 6.08e−03 |

**Interpretation that IS allowed:**
- "External advanced-disease replication" of the driver-orthogonal differentiation axis: in an independent ATC/PDTC cohort, the same RAI-lineage silencing readout separates groups with very large effect, and it reproduces the within-cohort orthogonality between the RAI_8 readout and the THYROID_NONOVERLAP module (ρ ≈ −0.92).
- Tumor-level (bulk) TROP2 elevation in advanced disease is consistent with TROP2 as a *tumor-population* vulnerability signal.

**Interpretation that is FORBIDDEN (do not write):**
- "Proves PTC → PDTC → ATC progression"
- "Monotonic 3-stage progression gradient"
- "DM1-high spots are TROP2-high in advanced disease" (GSE76039 is bulk microarray; no spatial information)
- Any survival, mutation, BRAF/RAS-status, age, or stage-ordinal claim from this cohort (the GEO record does not contain those fields)

**Cross-platform discipline applied:** TCGA-trained absolute-form classifier was **not** transferred to GPL570 microarray. Only within-cohort z-score normalisation was used, in accordance with `v17_korean_K2_calibration` (TCGA absolute-form → non-RNA-seq is wrong-direction).

---

## 3. H&E → DM1 status

- **Closure battery 2026-05-04 verdict D = NO-GO.** Discarded.
- **No re-entry during the 2026-05-04 → 2026-06-13 manuscript marathon.**
- Forbidden phrases (must not appear as active claims): "H&E-inferable molecular subtype", "H&E predicts DM1", "morphology-derived DM1", "pathology-AI triage of DM1", "tile-level DM1 inference", "WSI-validated DM1".
- Sweep (§6) confirms these phrases appear **only inside deprecation guards / forbidden-language lists**, never as active claims.

---

## 4. Files ready for commit (this session)

Grouped in §[Commit groups proposed] of this session — see prompt response. The candidate set is:

- `project/reports/2026_05_04_8gene_curated_vs_denovo_final_strategy.md`
- `project/reports/2026_05_04_paper1_strategy_plus_hook_checklist.md`
- `project/reports/2026_05_03_manuscript_v8_OUTLINE.md`
- `project/reports/2026_05_03_methods_M1_M11_scaffold.md`
- `project/reports/2026_05_03_methods_M1_M11_prose.md`
- `project/reports/2026_05_03_figure_captions_all.md`
- `project/reports/2026_05_04_supp_methods_closure_negative_feasibility.md`
- `project/reports/2026_05_04_minimal_external_data_plan.md`
- `project/reports/2026_05_04_gse76039_access_check.md`
- `project/reports/2026_05_04_gse76039_first_pass_report.md`
- `project/notebooks_or_scripts/p_landa_2016_first_pass.py`
- `project/results/p_landa_2016/sample_metadata.tsv`
- `project/results/p_landa_2016/expression_matrix_log_gene.tsv.gz`
- `project/results/p_landa_2016/probe_to_gene_panel.tsv`
- `project/results/p_landa_2016/scores.tsv`
- `project/results/p_landa_2016/score_tests.tsv`
- `project/results/p_landa_2016/SuppTable_SX_LANDA_GSE76039.xlsx`
- `project/results/p_landa_2016/gse76039_dm1_lineage_boxplot.png`
- `project/results/p_landa_2016/gse76039_dm1_vs_nonoverlap_scatter.png`
- `project/results/p_landa_2016/gse76039_mechanism_heatmap.png`
- `project/reports/2026_05_04_paper1_immediate_vscode_action_report.md` (this report)
- `project/reports/2026_05_04_hook_workspace_READY.md`

---

## 5. Files explicitly EXCLUDED from commit

- `project/results/p_landa_2016/raw/` — ~89 MB raw GEO downloads (Series Matrix + GPL570 SOFT). Already covered by `.gitignore` line 74. Reproducible from the URLs in `2026_05_04_gse76039_access_check.md`.
- WSI files, embeddings, raw tiles, RunPod runtime files.
- Paper 2 / Paper 3 / Paper 4 files (out of scope this session).
- Voice-protected prose drafts (none should be created by Claude this session anyway).

---

## 6. Overclaim sweep — results (PASS)

Phrases swept across active scaffold + GSE76039 report + final strategy doc:

```
proves progression               → 0 hits
PTC to PDTC to ATC               → 0 hits
H&E predicts DM1                 → 1 hit  (forbidden-language list — false positive)
WSI-validated                    → 1 hit  (forbidden-language list — false positive)
spatially colocalizes            → 1 hit  (forbidden-language list — false positive)
Nat Cancer reach                 → 1 hit  (forbidden-language list — false positive)
fusion-driven, epigenetically silenced → 0 hits
promoter hypermethylation        → 0 hits

monotonic                        → 3 hits — ALL guards
  - reports/2026_05_04_gse76039_first_pass_report.md:173 (explicit "do not claim ... monotonic")
  - reports/2026_05_03_manuscript_v8_OUTLINE.md:301 ("NOT as 'PTC → PDTC → ATC monotonic gradient'")
  - reports/2026_05_03_methods_M1_M11_prose.md:70 ("monotonic-gradient claims are not made from this cohort alone")

H&E-inferable                    → 4 hits — ALL in forbidden-language lists
DM1-high spots are TROP2-high    → 4 hits — ALL in forbidden / Q3 NEG / manuscript-safe-framing context
all risks resolved               → 2 hits — both inside "→ false" or forbidden lists
Cancer Cell-ready                → 3 hits — all in forbidden / venue-language guard lists
```

**Verdict:** every phrase hit is either zero-hit or appears strictly inside a deprecation guard / forbidden-language list / explicit negation. **No active claim contains a forbidden overclaim phrase.** Sweep PASS.

---

## 7. Next action

1. User decides commit strategy (see commit-group proposal in the session response).
2. Then: **user opens `04_intro_1_1_hook.md` and writes Hook ¶1 by hand**, using the fact-only workspace at `project/reports/2026_05_04_hook_workspace_READY.md`.
3. Claude will only audit Hook ¶1 after the user's draft (factual accuracy, overclaim, scope contamination, citation risk) — Claude will not draft prose.

---

**End of report — Claude did not draft Hook/Aim/Discussion/Limitations/Cover/Q9 prose in this session.**
