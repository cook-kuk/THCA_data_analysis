# Paper 1 — DM1_FULL bundle molecular-only LOCK

**Date:** 2026-05-04 · **Owner:** Seungho Cook
**Authority:** This memo + sibling `2026_05_04_paper1_dm1_full_conflict_audit.md` (claim-by-claim source verification by parallel auditor) + `2026_05_04_image_dm1_final_nogo_decision.md` (closure battery decision) + `2026_05_04_paper2_post_image_dm1_nogo_status.md` + `2026_05_04_marathon_state_after_image_dm1_nogo.md`.

**Scope of this memo**: lock the manuscript-safe interpretation envelope for `PAPER1_DM1_FULL_2026_05_04.md`. Defines what may be cited, what is deprecated, and the allowed / forbidden language. **No new analysis, no manuscript prose modifications, no RunPod/GPU dispatch, no new download.**

---

## 1. Executive verdict

**Paper 1 push continues — molecular only.**

- Molecular DM1 / RAI axis is **strengthened** by the multi-layer bundle audited in `2026_05_04_paper1_dm1_full_conflict_audit.md` §2 (16/16 source TSVs/PNGs verified).
- H&E → DM1 image angle remains **dropped** (closure battery 2026-05-04, verdict D = full 폐기). Re-entry blocked under the 5 conditions in `2026_05_04_image_dm1_final_nogo_decision.md` §5.
- RunPod G1 / G2 / G3 dispatch and TCGA WSI Phase C download remain **blocked**.
- Paper 1 manuscript should rely **only** on molecular / bulk RNA / spatial transcriptomics / TF-network / TROP2 / drop-one-out evidence for its claims.

`PAPER1_DM1_FULL_2026_05_04.md` is therefore an **internal molecular analysis bundle** — not a manuscript draft, not a submission doc — and its §9 "GPU Phase B 패키지" + §11 "RunPod 에서 G1+G2+G3 실행" are **superseded** by the closure battery NO-GO. The conflict auditor (parallel agent) has already added inline `[MANUSCRIPT-SAFE candidate]` / `[INTERNAL-ONLY]` / `[DEPRECATED]` per-section tags and a top-of-file banner to the bundle. **No further banner required from this lock.**

---

## 2. Manuscript-safe claims (cite-allowed in Paper 1)

Sourced from `2026_05_04_paper1_dm1_full_conflict_audit.md` §2 (each row source-verified):

| # | Claim | Source TSV | n | Paper 1 placement |
|---|---|---|---|---|
| 1 | TCGA bulk DM1 vs THYROID_NONOVERLAP r = −0.885 (p = 5.1×10⁻¹⁸⁸) | `s_tcga_thca_scored.tsv` | 561 | **Main** validation panel |
| 2 | TCGA DM1 cluster vs DM1_like score (DM1=0.00, DM2=−0.66, MW p≪0.001) | `s_tcga_thca_scored.tsv` | 157 | **Main** anchor |
| 3 | PFI dichotomized HR = 2.04 [1.15–3.61] p=0.015; DFI multivariate HR=1.41 p=0.025 (BRAF/RAS/age/stage adjusted) | `o_multivariate_cox.tsv`, `o_multivariate_cox.png` | 461 | **Main** KM + Cox forest |
| 4 | TF activity collapse — FOXE1 −1.21 / NKX2-1 −0.65; activation STAT3 +1.55 / FOSL1 +0.98 / DNMT1 +0.74 | `m_tf_activity_diff.tsv`, `m_master_regulator_volcano.png` | 261 | **Main** mechanism panel |
| 5 | TROP2 / TACSTD2 Δmean +4.33, FDR=1.2×10⁻³² (DE in DM1-high vs DM1-low) | `t_DM1_DE_genes.tsv`, `t_drug_target_volcano.png` | 261 | **Main** volcano + drug overlay (tumor-level only) |
| 6 | Tumor vs normal TROP2/DM1 window — DM1 p=1.8×10⁻²⁰, TROP2 1.4× tumor>normal | `n2_thca_tumor_vs_normal.tsv`, `n2_tumor_vs_normal.png` | 517 | **Main** sub-panel |
| 7 | Hallmark GSEA — IL6_JAK_STAT3 +8.0 (p=6×10⁻¹⁵), OXPHOS −5.5 (p=1×10⁻⁷), EMT, IFN-γ, TNF-α/NF-κB up | `q1_gsea_hallmark.tsv`, `q1_gsea_hallmark.png` | 261 | **Main** mechanism (or Supp if §4 panel busy) |
| 8 | TF network coordination — FOXE1-PAX8 r=+0.55, PAX8-DIO1 r=+0.64, HHEX-PAX8 r=+0.54 | `q5_thyroid_tf_network.tsv`, `q5_thyroid_tf_network.png` | 261 | **Supp** |
| 9 | Bootstrap 95% CI on 12-slide sample-mean: [−0.997, −0.916] | `d3_bootstrap_ci.tsv` | 12 | **Supp** |
| 10 | 8-gene drop-one-out r ∈ [−0.987, −0.981] (Δ < 1%) | `a1_drop_one_out.tsv`, `a1_drop_one_out_summary.tsv`, `a1_drop_one_out.png` | 28 | **Supp** referee Q (panel robustness) |
| 11 | Moran's I per slide mean 0.36, 26/28 perm p < 0.05 | `a2_morans_i.tsv`, `d1_permutation_morans.tsv`, `a2_morans_i.png` | 28 | **Supp** spatial structure |
| 12 | Bivariate Moran (DM1 × NONOVERLAP) all 12 negative, p ≤ 0.01 | `d2_bivariate_morans.tsv` | 12 | **Supp** spatial anti-correlation |
| 13 | Margin-distance gradient — 23/28 negative ρ | `d5_margin_distance.tsv` | 28 | **Supp** |
| 14 | Pan-cancer THCA r=−0.892 outlier (33 cancer types) | `p_pancancer_scored.tsv`, `p_pancancer_summary.tsv`, `p_pancancer.png` | 33 | **Supp** thyroid-specific framing |
| 15 | Dark-matter subset (BRAF-neg ∩ RAS-neg, n=156) PFI HR=1.20, p=0.80 (NS) | `n1_dark_matter_subset.tsv`, `n1_dark_matter.png` | 156 | **Supp** + reviewer Q (honest negative) |
| 16 | Harmony 28-slide UMAP (ST QC backbone) | `c2_harmony_28slides.h5ad`, `c2_harmony_umap.png` | 28 | **Supp** Methods |

These are the **citation envelope** for Paper 1's molecular/spatial story. Anything outside this envelope is **not** a Paper 1 claim under the current marathon.

---

## 3. Excluded / deprecated claims (must NOT appear in Paper 1 manuscript-facing text)

| Claim | Reason | Authority |
|---|---|---|
| H&E → DM1 prediction at tile level | NO-GO: r=0.022, AUROC=0.511, real ≤ random max | `pathology_dm1_phaseA_cpu_verdict_2026_05_04.md`, `pathology_dm1_closure_battery_2026_05_04.md` |
| TCGA WSI external validation of DM1 | Phase C never executed, blocked under marathon | `2026_05_04_image_dm1_final_nogo_decision.md` §3 |
| RunPod G1 / G2 / G3 pipeline | Conflicts with closure battery; same data + pipeline already failed | conflict_audit §3 |
| Image-DM1 IP / patent claim | Color-artifact unpatentable; ATC classifier prior art dense | `2026_05_04_image_dm1_final_nogo_decision.md` §3 |
| "All risks resolved" / "TDS overlap, inflammation artifact 모두 해소" | N1 dark-matter subset NS (HR=1.20, p=0.80); Q3 TROP2 spot-level NEG (ρ=−0.016); H&E NO-GO. ≥3 risks remain. | conflict_audit §2 (subtractive honesty) |
| "Cancer Cell 도전권 20-30%" / "Nat Cancer 50-65% reach" / venue-probability framing | Speculative internal planning; never submission-facing | conflict_audit §1 |
| "DM1-high spots are TROP2-high" / spot-level colocalization | Q3 cross-slide mean ρ(TROP2, DM1_resid) = −0.016. Sample-level (N2) tumor>normal works; spot-level does NOT. | conflict_audit §2 (Q3) |
| "H&E-inferable molecular subtype" | Closure battery NO-GO | as above |

---

## 4. Required cautious language (use these exact phrasings or close paraphrases)

| context | preferred phrasing |
|---|---|
| epidemiology / outcome association | "DM1_like score **supports** independent prognostic association" — NOT "proves"/"establishes" |
| TF-mechanism | "**consistent with** a 4-step framework (TF collapse → DNMT silencing → STAT3/AP-1 activation → TROP2 re-expression)" — NOT "demonstrates the framework" |
| TROP2 therapeutic claim | "**tumor-level therapeutic vulnerability**" or "tumor-population vulnerability" — NOT "DM1-axis-specific TROP2 colocalization" |
| spatial transcriptomics | "**supportive spatial validation** in 28 ST slides" — NOT "spatial validation alone resolves" |
| wet-lab status | "**functional validation remains required** (TROP2 IHC, sacituzumab IC50, organoid)" — must appear in Limitations |
| H&E negative result (if cited) | "We **pre-tested** H&E-based prediction of the depth-residualized DM1 axis at GSE250521 hires resolution and found **no learnable signal beyond random panels**. Image-based triage of the molecular subtype is therefore not pursued in this work." |

---

## 5. Forbidden language (Paper 1 manuscript text — STOP if seen)

- `proves` / `establishes` (use `supports` / `consistent with`)
- `all risks resolved` / `모든 risk 해소` (false — see §3)
- `Cancer Cell-ready` / `Cancer Cell 도전권` / `Nat Cancer reach` / venue-probability percentages
- `H&E-inferable` / `H&E predicts DM1` / `pathology-AI triage of DM1`
- `TROP2 spatially colocalizes with DM1-high spots` (Q3 negative)
- `DM1-high spots are TROP2-high` (Q3 negative)
- `tile-level DM1 inference` / `WSI-validated DM1`
- `RunPod G1/G2/G3` / `image-DM1 Phase B/C` (in any active-claim form)
- `morphology-derived DM1`

If a draft contains any of these, replace per §4 or delete.

---

## 6. Figure recommendation (Paper 1)

### Main figure candidates (cite-ready under §2 envelope)

| panel | content | source rows in §2 |
|---|---|---|
| A | TCGA DM1 cluster vs DM1_like score (anchor) | #2 |
| B | TCGA DM1 vs THYROID_NONOVERLAP scatter (validation) | #1 |
| C | KM + Cox forest (PFI / DFI multivariate) | #3 |
| D | TF mechanism (FOXE1↓ NKX2-1↓ STAT3↑ FOSL1↑ DNMT1↑) volcano | #4 |
| E | TROP2 / TACSTD2 DE volcano + drug overlay | #5 |
| F | Tumor vs normal TROP2/DM1 window | #6 |

Verify source PNGs/TSVs exist before final selection (conflict_audit §2 confirmed all 16 source files present).

### Supplementary candidates

- ST clustering (Harmony 28-slide UMAP) — #16
- 8-gene drop-one-out robustness — #10
- Pan-cancer THCA-outlier specificity — #14
- Hallmark GSEA mechanism cross-check — #7
- TF network coordination (FOXE1-PAX8 / PAX8-DIO1) — #8
- Bootstrap 95% CI 12-slide sample-mean — #9
- Moran's I + bivariate Moran spatial structure — #11, #12
- Margin-distance gradient — #13
- N1 dark-matter underpowered (honest-negative supp) — #15
- Closure battery negative-feasibility (1-page Methods supplement) — closure battery files

### Internal-only (do NOT include in submission package)

- §9 GPU Phase B / RunPod G1-G2-G3 package recommendation — superseded
- §11 "RunPod 에서 G1+G2+G3 실행" next-step block — superseded
- venue-probability framing (Cancer Cell 20-30%, Nat Cancer 50-65%) — speculative
- H&E-DM1 / pathology-AI / WSI-projection wording in any active-claim form

---

## 7. Final lock

**Paper 1 marathon proceeds molecular-only.** No retry of image-DM1, no RunPod / TCGA-WSI dispatch, no foundation-model attempts during the 5/4–6/13 window.

Re-entry to image-DM1 is permitted **only** under the 5 conditions of `2026_05_04_image_dm1_final_nogo_decision.md` §5 (full-res WSI + external molecular labels + post-marathon explicit decision + foundation model access acquired + pre-registered hypothesis). None are met currently.

`PAPER1_DM1_FULL_2026_05_04.md` is preserved on disk as an **internal molecular analysis bundle**. The conflict auditor's inline tags (`[MANUSCRIPT-SAFE candidate]` / `[INTERNAL-ONLY]` / `[DEPRECATED]`) provide per-section guidance — that is sufficient. No additional Claude-side banner editing required.

---

## 8. Cross-references

- `2026_05_04_paper1_dm1_full_conflict_audit.md` — claim-by-claim source verification (auditor agent)
- `2026_05_04_image_dm1_final_nogo_decision.md` — closure decision authority
- `2026_05_04_paper2_post_image_dm1_nogo_status.md` — Paper 2 scope cleanup
- `2026_05_04_marathon_state_after_image_dm1_nogo.md` — marathon snapshot
- `CLOSURE_BATTERY_2026_05_04.md` — top-level closure summary
- `pathology_dm1_phaseA_cpu_verdict_2026_05_04.md` — Phase A NO-GO baseline
- `pathology_dm1_closure_battery_2026_05_04.md` — full forensic battery
- `STATUS_PAPER1_2026_05_02.md` — Paper 1 status (existing)
- Memory: `v17_marathon_mode_post_pillar1`, `v17_sprint_vs_marathon_violation`, `v18_paper2_HT_isolated`, `v17_npj_ship_status`, `v17_dark_matter_pivot_2026_04_29`
