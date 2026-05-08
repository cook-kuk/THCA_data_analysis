# Paper 1 — Comprehensive Audit MASTER

**Date:** 2026-05-07
**Scope:** Paper 1 one-page comprehensive audit + synthesis + web packaging
**Mode:** audit + synthesis only — voice-protected manuscript prose 작성 NOT
**URL:** http://40.82.129.113:8012/manuscript_v8/p1_onepage_audit.html

---

## 1. Inventory summary

- 451 result files across `project/results/` related to Paper 1 (DM1 / RAI / external / GPL570 / 76039 / TCGA / score / TF / NONOVERLAP / TIERA / driver / feature / validation)
- 60+ Paper 1 reports under `project/reports/`
- Inputs read: `2026_05_04_paper1_final_strategy_after_external_validation.md`, `2026_05_03_manuscript_v8_OUTLINE.md`, `2026_05_03_methods_M1_M11_scaffold.md`, `2026_05_03_figure_captions_all.md`, `2026_05_04_hook_workspace_READY.md`, `2026_05_04_external_expression_FULL_RESULTS.md`
- Web server: PID 1690659 (THCA-Secure), docroot = `/home/seungho/personal/THCA_data_analysis/project/`, port 8012

---

## 2. Dataset appropriateness audit

(See `2026_05_06_paper1_dataset_suitability_registry.tsv` for 21 dataset × 19-column full table.)

**Main 11 (discovery + external + mechanism support):** TCGA-THCA, TCGA HM450, cBioPortal SV/TERT, MSK-IMPACT (Landa 2016), GSE76039, GSE33630, GSE65144, GSE213647 (Lee 2024), K2/PRJEB11591, GSE184362 (Pu 2021), GSE193581 (Lu 2023).

**Supplementary 6:** GSE29265 (paired-N), GSE53157 (PDTC sensitivity, n=5 underpowered), GSE241184 (Phase 1 sc), GSE250521 (spatial supp; sample-mean stage trend MARGINAL), Pozdeyev 2018 supp table, TIERA67 candidate pool table.

**Drop 4:** GSE126698 (Series Matrix metadata-only; raw=SRA forbidden), GSE60542 (lymphoid confounding), GSE286332 (Paper 2 territory), GSE151179 (Paper 3 territory).

**Misleading-if-overinterpreted:** MSK-IMPACT + GSE76039 100% advanced — SA6 bias panel + Methods M1 disclosure mandatory. K2 mini-index inflate 4.9-12.5× — within-sample-centered profile only. GSE53157 PDTC arm n=5 → direction-consistent only. GSE250521 sample-mean stage trend p>0.6 → supportive supp only.

---

## 3. Gene-set hierarchy audit

(See `2026_05_06_paper1_gene_set_hierarchy.tsv` for 17-entity table.)

**Final hierarchy (top → bottom):**
- Pan-genome top-5000 MAD (data-driven, ARI 0.92) — discovery layer
- TIERA67 (literature-curated 67-gene; 7 sub-categories; ARI 0.90) — discovery layer
- TDS-16 / TDS_core (canonical Yoo 2016 16-gene; AUC 0.975) — comparator
- **RAI_8 (8-gene main deployable readout; AUC 0.962, ARI 0.49 by design)**
- THYROID_NONOVERLAP (8 zero-overlap orthogonal lineage genes; cross-panel ρ +0.44)
- TF_collapse (4 TF backbone: FOXE1, NKX2-1, PAX8, HHEX)
- STAT3_AP1_DNMT mechanism arm (5 genes)
- DM1_like_score (= − RAI_8) + DM1 cluster label (KMeans k=2)
- DM2 / DM1 sub-A / sub-B / not_DM partition labels

**Drop / demote:** TROP2 main claim, 10-gene q5_thyroid_tf_network ad-hoc, "55-gene driver-excluded" superseded.

---

## 4. Figure logic audit

(See `2026_05_06_paper1_figure_logic_registry.tsv` for 25-figure registry.)

**Recommended main figure plan (current strategy):**
- F1. Clinical workflow current vs future
- F2. Compact 8-gene readout + axis robustness (ARI ladder)
- F3. External GPL570 4-cohort + Korean (DM1_like vs NONOVERLAP scatter grid + forest)
- F4. Aggressive cohort + survival meta (PFI primary; OS secondary; bias panel SA6)
- F5. Mechanism support (HM450 promoter methylation + TF collapse)

**Supp:** TIERA67 ARI ladder, single-cell 3-cohort, spatial GSE250521 (marginal), MSK-IMPACT bias panel SA6, panel-size sensitivity scan ★ NEW.

**Drop:** 10-gene q5 (TRPS1 ad-hoc), H&E/WSI pathology, TROP2 main claim, GSE126698 metadata-only, GSE60542 lymphoid.

---

## 5. Clinical workflow and translational interpretation

**1. Current workflow:** surgery → empirical RAI → 6-12 months delay → Tg trend → recognition of failure → systemic / molecular review.

**2. Proposed hypothesis (NOT validated):** surgery → compact RNA readout (8-gene RT-qPCR) → DM1 high-risk flag → earlier review of systemic / molecular options.

**3. Important guardrails:**
- NOT a validated treatment selection tool
- NOT replacing RAI today
- NOT prospectively validated
- NOT a direct clinical recommendation

**4. Why this matters clinically:**
- reduced delay (months potential)
- triage logic grounded in canonical thyroid biology
- low-cost panel possibility (RT-qPCR / NanoString; deployable in any pathology lab)
- biologically interpretable — every gene is daily-clinical-recognized molecule

---

## 6. Claim boundary + reviewer-risk audit

(See `2026_05_06_paper1_claim_boundary_and_reviewer_risk.tsv` for 21 topic × 8-column table.)

**Top 10 reviewer risks:**
1. Predicts RAI response — VERY HIGH (forbidden claim)
2. Clinical utility overclaim — VERY HIGH (forbidden)
3. TROP2 demotion — VERY HIGH (forbidden in title/main)
4. H&E pathology branch — VERY HIGH (forbidden; closure complete)
5. Just dedifferentiation relabeled — HIGH
6. OS swap / age confound — HIGH
7. Stage doesn't separate dramatically — HIGH (orthogonality)
8. Sub-A/B versioning gap — HIGH (caveat explicit)
9. Spatial signal weak — HIGH (supportive supp only)
10. External cohorts aggressive-biased — MEDIUM (SA6 disclosure)

**Allowed wording:** compact RAI-lineage readout / direction-consistent lineage silencing / candidate triage scaffold / consistent with / supports / PFI primary / East-Asian generalizability.

**Forbidden wording:** validated RAI predictor / treatment selection / clinical utility proven / TROP2 vulnerability (title) / H&E inferable / demonstrates causality.

---

## 7. Gene panel combination scan ★ NEW analysis

(Output: `project/manuscript_v8/assets/p1_onepage_audit/panel_combos.tsv` + `fig_panel_combos.png`)

**Top findings:**
- TIERA TDS-extended 12: AUC 0.917 (best)
- Lineage TF + Effector minimal 6 (FOXE1+NKX2-1+PAX8+TG+TPO+DIO1): AUC 0.893 + ARI 0.410 — **outperforms canonical RAI_8**
- Minimal 4 (FOXE1+NKX2-1+TG+TPO): AUC 0.872 — 4 genes match RAI_8
- NONOVERLAP_8 alone: AUC 0.877 — orthogonal panel matches RAI_8 → §4.3 defense reinforced
- Mechanism arm 5: AUC 0.171 (= 1 − 0.829) — direction reverse (DM1 ↑) confirms orthogonal mechanism axis

**Recommended action:** Add panel scan as supplementary figure → reviewer Q "why 8?" → strong answer that 8 is canonical compromise; 6/4 alternatives shown for completeness.

---

## 8. Output inventory

| File | Purpose |
|------|---------|
| `2026_05_06_paper1_inventory_files.txt` | full result-file inventory |
| `2026_05_06_paper1_dataset_suitability_registry.tsv` | 21 datasets × 19 columns |
| `2026_05_06_paper1_gene_set_hierarchy.tsv` | 17 entities |
| `2026_05_06_paper1_figure_logic_registry.tsv` | 25 figures |
| `2026_05_06_paper1_claim_boundary_and_reviewer_risk.tsv` | 21 topics |
| `2026_05_06_paper1_comprehensive_audit_master.md` | this file |
| `project/manuscript_v8/p1_onepage_audit.html` | one-page web report |
| `project/manuscript_v8/assets/p1_onepage_audit/` | 19 figures + scripts + panel_combos.tsv |

---

## 9. Forbidden-action audit (sweep total)

- ✅ No new GEO accessions queried
- ✅ No raw FASTQ / CEL / WES download
- ✅ No GPU / RunPod
- ✅ No H&E / WSI / pathology retry
- ✅ No TROP2 main claim restoration
- ✅ No Paper 2/3/4 territory work
- ✅ No voice-protected prose written
- ✅ No commit performed
