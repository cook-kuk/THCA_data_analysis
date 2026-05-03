# Paper 1 — `PAPER1_DM1_FULL_2026_05_04.md` × Closure Battery Conflict Audit

**Author:** Seungho Cook
**Date:** 2026-05-04 (post-closure recheck, marathon mode)
**Scope:** Reconcile new claims in `PAPER1_DM1_FULL_2026_05_04.md` (top-level, 295 lines) against the closure-battery NO-GO and three sibling memos. Lock what Paper 1 may actually cite. **No analysis run. No manuscript prose modified. No RunPod / GPU dispatch. No new download.**

**Inputs:**
- `PAPER1_DM1_FULL_2026_05_04.md` — full sprint bundle, includes a §9 RunPod G1/G2/G3 recommendation
- `CLOSURE_BATTERY_2026_05_04.md` — top-level Korean NO-GO (option D 폐기)
- `PHASE_A_NOGO_REPORT_2026_05_04.md` — Phase A 16-fold LOSO baseline NO-GO
- `project/reports/2026_05_04_image_dm1_final_nogo_decision.md` — decision authority memo
- `project/reports/2026_05_04_paper2_post_image_dm1_nogo_status.md` — Paper 2 scope cleanup
- `project/reports/2026_05_04_post_closure_recheck.md` — 9-point bundle integrity recheck
- `SITUATION_BACKGROUND_ACTIVITY_2026_05_04.md` — RunPod dispatch state
- Per-claim source TSVs in `project_external_st/results/extra/` (12 files verified to exist)

---

## 1. Executive verdict

**Paper 1 = molecular DM1 axis, full GO. Paper 2-territory image-DM1 = full STOP.**

`PAPER1_DM1_FULL` mixes two scope worlds:

| Sections | Scope (per `paper_numbering_2026_05_04` memory) | Verdict |
|---|---|---|
| §0–§8, §10–§11 (10 sprint layers + venue + sprint score) | **Paper 1** = DM1 molecular dark matter | ✓ all source TSVs exist; cite-ready with caveats below |
| §9 *"GPU Phase B 패키지 (RunPod 용)"* + §11 *"RunPod 에서 G1+G2+G3 실행"* | **Paper 2** = H&E→DM1 pathology projection — already CLOSED **D 폐기** 2026-05-04 | ✗ **conflicts** with closure battery + 3 sibling memos. **Discard from this bundle.** |
| Bundled venue claims (Cancer Cell 20-30%, Nat Cancer 50-65% reach) | speculative | internal-memo only, never submission-facing |

The bundle is therefore **not internally consistent** — its narrative escalates Paper 1 to a Cancer-Cell push *contingent on* G1/G2/G3 hits that the closure battery has already rejected at hard NO-GO. The molecular core is independently strong and does not need that escalation. Lock it on its own merits, peel off the pathology recommendation.

**Action items:** preserve §0–§8/§10/§11 as a Paper 1 internal scoping memo (not a submission doc) with venue language softened; mark §9 + RunPod G1/G2/G3 dispatch as **superseded by closure battery NO-GO**; manuscript marathon resumes unchanged.

---

## 2. What is genuinely strengthened (Paper 1 molecular core)

Per-claim verification: each row checks (a) source TSV/PNG exists, (b) sample size, (c) whether the claim is novel relative to the existing 2026-05-03 Paper 1 manuscript scaffold, (d) Paper 1 placement.

| Claim | Source file (verified ✓) | n | Novel vs. 5/3 scaffold? | Paper 1 placement |
|---|---|---|---|---|
| TCGA bulk DM1 vs THYROID_NONOVERLAP r = −0.885 (p = 5.1e-188) | `s_tcga_thca_scored.tsv`, `s_tcga_4panel.png` | 561 | ★ already framework-anchor | **Main figure candidate** (validation panel) |
| Bootstrap 95% CI on 12-slide sample-mean: [−0.997, −0.916] | `d3_bootstrap_ci.tsv` | 12 | extends Layer 1 | **Supplementary** |
| PFI dichotomized HR=2.04 [1.15–3.61] p=0.015; DFI multivar HR=1.41 p=0.025 | `o_multivariate_cox.tsv`, `o_multivariate_cox.png` | 461 | ★ supports clinical claim | **Main figure** (KM + Cox forest) |
| 8-gene drop-one-out r ∈ [−0.987, −0.981] (Δ < 1%) | `a1_drop_one_out.tsv`, `a1_drop_one_out_summary.tsv`, `a1_drop_one_out.png` | 28 slides | ★ panel robustness | **Supplementary** (referee Q) |
| Moran's I per slide mean 0.36, 26/28 perm p<0.05 | `a2_morans_i.tsv`, `d1_permutation_morans.tsv`, `a2_morans_i.png` | 28 | ★ spatial structure | **Supplementary** |
| Bivariate Moran (DM1 × NONOVERLAP) all 12 negative, p ≤ 0.01 | `d2_bivariate_morans.tsv` | 12 | ★ spatially anti-correlated | **Supplementary** |
| Margin distance gradient — 23/28 negative ρ | `d5_margin_distance.tsv` | 28 | enriches spatial story | **Supplementary** |
| TF activity collapse FOXE1 −1.21 / NKX2-1 −0.65; activation STAT3 +1.55 / FOSL1 +0.98 / DNMT1 +0.74 | `m_tf_activity_diff.tsv`, `m_master_regulator_volcano.png` | 261 | ★ 4-step mechanism | **Main figure candidate** (mechanism panel) |
| TROP2/TACSTD2 Δmean +4.33, FDR=1.2e-32 | `t_DM1_DE_genes.tsv`, `t_drug_target_volcano.png` | 261 | ★ druggable target | **Main figure** (volcano + drug overlay) |
| Tumor-vs-normal: DM1 p=1.8e-20, TROP2 1.4× tumor>normal | `n2_thca_tumor_vs_normal.tsv`, `n2_tumor_vs_normal.png` | 517 | ★ therapeutic window | **Main figure** (sub-panel) |
| Pan-cancer THCA r=−0.892 outlier | `p_pancancer_scored.tsv`, `p_pancancer_summary.tsv`, `p_pancancer.png` | 33 cancers | ★ thyroid-specific framing | **Supplementary** |
| Hallmark GSEA — IL6_JAK_STAT3 +8.0 (p=6e-15), OXPHOS −5.5 (p=1e-7), EMT, IFN-γ, TNF-α/NF-κB up | `q1_gsea_hallmark.tsv`, `q1_gsea_hallmark.png` | 261 | ★ independent validation of M | **Supplementary** (could be Main if mechanism is main figure) |
| TF network coordination — FOXE1-PAX8 r=+0.55, PAX8-DIO1 r=+0.64, HHEX-PAX8 r=+0.54 | `q5_thyroid_tf_network.tsv`, `q5_thyroid_tf_network.png` | 261 | ★ coordinated collapse | **Supplementary** |
| Dark matter subset (BRAF-neg ∩ RAS-neg, n=156) PFI HR=1.20 NS | `n1_dark_matter_subset.tsv`, `n1_dark_matter.png` | 156 | ★ honest negative | **Supplementary** + reviewer Q (honesty asset) |
| TCGA DM1 cluster vs DM1_like score: DM1=0.00, DM2=−0.66, MW p<<0.001 | `s_tcga_thca_scored.tsv` | 157 | ★ consilience with manuscript dark-matter axis | **Main figure** (anchor panel; already in v8 outline?) |
| Cell-type × DM1 quartile distribution | `a3_celltype_x_dm1.tsv`, `a3_DM1_quartile_x_celltype.tsv`, `a3_per_spot_celltype.tsv.gz`, `a3_dm1_quartile_celltype.png` | 28 | tangential — useful for ST methods | **Supplementary** |
| Harmony 28-slide UMAP (c2) | `c2_harmony_28slides.h5ad`, `c2_harmony_umap.png` | 28 slides | ST QC backbone | **Supplementary** (Methods) |

All 12 result TSVs and the listed PNGs exist on disk (`ls` verified). No claim depends on a missing source.

**Subtractive honesty (already in PAPER1_DM1_FULL):**

- **N1 dark matter subset is NS** (HR=1.20, p=0.80, n=156). Treat as positive epistemic asset — the audit has already pre-empted the obvious referee Q.
- **Q3 TROP2 spatial colocalization is NEGATIVE** (mean ρ=−0.016 across 28 slides; sample-level still works). The bundle's reframe — "tumor=TROP2-high overall, DM1-axis=dedifferentiation; sacituzumab population is tumor (vs normal), not DM1-high specifically" — is the right narrative; Cancer-Cell-tier "DM1-high spots = TROP2-high spots" is **not** claimable.

These two negatives are the bundle's most credible feature, and they are independent of any pathology / RunPod work.

---

## 3. What conflicts with the closure battery (and must be peeled off)

`PAPER1_DM1_FULL` §9 *"GPU Phase B 패키지"* + §11 *"다음 step → RunPod 에서 G1+G2+G3 실행"* directly contradict the closure battery.

| `PAPER1_DM1_FULL` claim/recommendation | Closure-battery counter-evidence | Verdict |
|---|---|---|
| §9 — package 4.35 GB ready, run `g1_embed.py --model resnet50` then ridge LOSO for DM1 prediction | `PHASE_A_NOGO_REPORT_2026_05_04.md` §4: same pipeline, same data (3,200 GSE250521 tiles, 16-fold LOSO), DM1_resid Ridge **r=0.022, AUROC=0.511**. Real ≤ random max. | **Already done. Already failed. Re-running ≡ violating spec.** |
| §9 — UNI foundation model retry (gated, ~60 min) | `2026_05_04_image_dm1_final_nogo_decision.md` §2.5 + §1 option E: published UNI/CONCH/Virchow2 gain over ResNet50 = +0.05–0.15 Spearman. Gap to GO = +0.28. Even max-gain lands at ~0.17, below BORDERLINE 0.20. | **Closed.** Re-entry blocked under §5 of decision memo (5 conditions, none currently met). |
| §9 — TCGA H&E (n=500) "digital pathology Cancer Cell finisher" 1–2 day | `2026_05_04_image_dm1_final_nogo_decision.md` §3: *"TCGA WSI Phase C (~500 GB download): never executed, never to execute under current marathon"*. `2026_05_04_paper2_post_image_dm1_nogo_status.md` §1: TCGA WSI Phase C marked **deprecated**. | **Closed.** No download. |
| §9 — CellDART G3 (DM1-high spots = which cell type) 2-3 h | Not directly NO-GO'd by closure battery, but it's a Paper-2-territory pathology-spot-deconvolution sprint. Per `2026_05_04_paper2_post_image_dm1_nogo_status` and `v17_marathon_mode_post_pillar1` memory: *"Pillar 1 STRONG = 분석 끝. D9+ sprint 없음. 새 분석은 paper-blocking 만"*. CellDART deconvolution does not gate any Paper 1 claim already in scaffold. | **Out of marathon scope.** Defer post-bioRxiv. |
| §11 — venue language: *"G1 LOSO Spearman > 0.4 hit 시 Nat Cancer 비율 70%+"* | G1 already at r=0.022 with full bundle. No path from 0.022 → 0.40 via foundation model alone. Conditional venue claim is therefore moot. | **Drop the conditional.** |
| §10 — *"PID 168753 별도 session 이 embeddings_resnet50_448.npz 로 G1 의 일부 이미 진행 중"* | Closure battery already ran 448 px (see `PHASE_A_NOGO_REPORT` §6 / closure §C): **r=0.017, slightly worse than 224's 0.022**. The "PID 168753" 448 work is the same line of inquiry, already concluded NO-GO. | **PID stale.** No "complementary" finding to expect. |
| Sprint score table §8 row reading "Q3 TROP2 spatial colocalization: ❌ NEGATIVE — sample-level still works" | Consistent with closure battery findings (spot-level H&E ↔ molecular labels mismatched, dominant cause #4 in closure). | **No conflict** — Paper 1 keeps the sample-level claim; do not upgrade to spot-level. |

**Summary of the conflict:** PAPER1_DM1_FULL was authored as if the H&E→DM1 leg were still under active investigation. The closure battery + final NO-GO decision memo + Paper 2 scope cleanup memo, all dated **the same day** (2026-05-04), establish that the H&E→DM1 angle is officially DROPPED. §9 of the bundle is therefore obsolete on arrival.

**Background activity check** (`SITUATION_BACKGROUND_ACTIVITY_2026_05_04.md`): Pod B (WSI) was killed before this audit started; Pods C/D (AlphaFold + K2 STAR) remain user-managed but are unrelated to image-DM1; the dispatch loop polled but did not dispatch. **No live RunPod compute is currently running on the H&E-DM1 angle.** Lock holds without further intervention.

---

## 4. Paper 1 usable claims — A/B/C tiering

Three tiers: **A = main-figure candidate**, **B = supplementary candidate**, **C = internal/do-not-use**.

### Tier A — Main figure candidates (Paper 1)

| # | Claim | Source | Why main |
|---|---|---|---|
| A1 | DM1 vs DM2 cluster axis (DM1 cluster mean=0.00, DM2 cluster mean=−0.66, MW p<<0.001) | `s_tcga_thca_scored.tsv` | Anchor — connects expression score to manuscript's mutation-based dark-matter axis (consilience) |
| A2 | TCGA bulk Pearson r=−0.885 between DM1_like and THYROID_NONOVERLAP (n=561) | `s_tcga_thca_scored.tsv`, `s_tcga_4panel.png` | Cross-validation between two independent gene signatures of the same axis |
| A3 | PFI dichotomized HR=2.04 (top30 vs bot30); DFI multivar HR=1.41 | `o_multivariate_cox.tsv`, `o_multivariate_cox.png` | Clinical relevance; survives BRAF/RAS/age/stage adjustment for DFI |
| A4 | TF activity 4-step framework (FOXE1/NKX2-1 collapse → DNMT1/3B up → STAT3/AP-1 up) | `m_tf_activity_diff.tsv`, `m_master_regulator_volcano.png` | Mechanism — single panel can show the 4 layers |
| A5 | TROP2 (TACSTD2) Δmean+4.33, FDR=1.2e-32 + sacituzumab govitecan FDA-approved | `t_DM1_DE_genes.tsv`, `t_drug_target_volcano.png` | Therapeutic axis (paired with A6 below for tumor-vs-normal therapeutic window) |
| A6 | Tumor-vs-normal therapeutic window: TROP2 normal RPKM 8.33 → tumor 11.76 → metastatic 13.12; tumor-vs-normal p=1.3e-13 | `n2_thca_tumor_vs_normal.tsv`, `n2_tumor_vs_normal.png` | Sacituzumab off-target safety argument |

If Paper 1 already has a different main-figure plan (per `2026_05_03_manuscript_v8_OUTLINE.md`), A1 + A4 + A5 + A6 should be considered as **augmentations** for revision rounds, not slot displacements. Author keyboard decides which become main.

### Tier B — Supplementary candidates

| # | Claim | Source |
|---|---|---|
| B1 | 12-slide bootstrap CI [−0.997, −0.916] | `d3_bootstrap_ci.tsv` |
| B2 | Drop-one-out 8-gene panel robustness (Δ < 1%) | `a1_drop_one_out_summary.tsv` |
| B3 | Moran's I per slide (mean 0.36, 26/28 perm p<0.05) + stage gradient (PT/PTC/LPTC 0.39–0.47, ATC 0.13) | `a2_morans_i.tsv`, `d1_permutation_morans.tsv` |
| B4 | Bivariate Moran (DM1 × NONOVERLAP), all 12 negative | `d2_bivariate_morans.tsv` |
| B5 | Margin-distance gradient (23/28 negative ρ) | `d5_margin_distance.tsv` |
| B6 | Pan-cancer specificity (THCA r=−0.892 outlier across 33 TCGA cancers) | `p_pancancer_summary.tsv`, `p_pancancer.png` |
| B7 | Hallmark GSEA (IL6_JAK_STAT3, OXPHOS, EMT, IFN-γ, TNF-α/NF-κB) | `q1_gsea_hallmark.tsv` |
| B8 | TF network coordination (FOXE1/PAX8/DIO1/HHEX pairwise r) | `q5_thyroid_tf_network.tsv` |
| B9 | Dark-matter subset NS (n=156, HR=1.20) — **honest negative** | `n1_dark_matter_subset.tsv` |
| B10 | Q3 TROP2 spatial spot-level NOT colocalized — **honest negative** | `q3_trop2_spatial.tsv` |
| B11 | Cell-type × DM1 quartile (CellDART/decon outputs) | `a3_*` |
| B12 | Harmony 28-slide UMAP (Methods Suppl.) | `c2_harmony_umap.png` |

### Tier C — Internal only / do-not-use in submission

| # | Item | Reason |
|---|---|---|
| C1 | All §9 RunPod G1/G2/G3 recommendations | Closure-battery NO-GO, see §3 above |
| C2 | TCGA H&E (n=500) "digital pathology Cancer Cell finisher" wording | Phase C deprecated 2026-05-04 |
| C3 | UNI/CONCH/Virchow2 retry plan | Foundation-model gain insufficient (gap 15×); blocked by 5-condition gate |
| C4 | "G1 LOSO Spearman > 0.4 hit 시 Nat Cancer 비율 70%+" conditional venue claim | G1 already 0.022 — premise impossible |
| C5 | Wet-lab roadmap (TROP2 IHC / cell-line IC50 / organoid / xenograft) §4 | Marathon mode = manuscript only; wet-lab is post-bioRxiv |
| C6 | Methylation correlate / STAT3 target enrichment / paired scRNA deconvolution / mutational signature COSMIC SBS / additional thyroid cohorts (GSE33630, GSE76039) §4 | Not paper-blocking; defer per `v17_marathon_mode_post_pillar1` |

---

## 5. Paper 1 excluded claims (explicit)

**Excluded from Paper 1 manuscript** until a separate post-marathon authorization gates each one:

1. **H&E → DM1 image-axis prediction** at any tile size, any embedder, any cohort. Bound only as a *negative-feasibility* statement in Methods Supplementary (per `2026_05_04_paper2_post_image_dm1_nogo_status.md` §5 recommendation): one sentence acknowledging the closure battery, no metric advertised as a positive result.
2. **TCGA WSI external validation** of DM1 program (image-based). Not in scope.
3. **Image-DM1 IP / patent claim.** Closure battery deemed unpatentable (color-artifact); ATC classifier prior art is dense.
4. **RunPod G1 (ResNet50 LOSO)** — already at NO-GO. **G1 (UNI)** — gain insufficient. **G2 (TCGA H&E)** — cohort never downloaded. **G3 (CellDART)** — out of marathon scope.
5. **"Digital pathology" / "pathology AI projection" / "tile-level DM1 inference" / "morphology-derived DM1" / "H&E foundation model triage"** as positive claim wording (per `2026_05_04_paper2_post_image_dm1_nogo_status.md` §3 forbidden-words extension; same applies to Paper 1 because the closure is cohort-level NOT paper-level).
6. **Paper 3 (ICI Track A) and Paper 4 (GD backlog)** — not touched.

These remain blocked **until the user issues a new post-marathon command explicitly re-opening one item** under all 5 re-entry conditions of `2026_05_04_image_dm1_final_nogo_decision.md` §5.

---

## 6. Figure / supplementary recommendations

This is a recommendation only. Author keyboard decides during voice-protected drafting (Hook ¶1, Aim ¶4, etc., per `v17_sprint_vs_marathon_violation`).

### Recommended Paper 1 main figure additions / consolidations

Treat these as *additions* to the existing v8 outline figure plan; do not displace the manuscript's existing anchor figures unless author decides:

- **Figure (mechanism / 4-step):** A4 TF volcano + small inset of FOXE1/NKX2-1 vs STAT3/AP-1/DNMT directionality (`m_master_regulator_volcano.png`).
- **Figure (clinical / Cox):** A3 Cox forest + KM curve for dichotomized DM1 (top30 vs bot30) (`o_multivariate_cox.png`, `s_tcga_km_dm1.png`).
- **Figure (therapeutic):** A5 + A6 paired panel — TROP2 Δmean volcano + tumor-vs-normal RPKM violin (`t_drug_target_volcano.png`, `n2_tumor_vs_normal.png`).

### Recommended supplementary figures

- B3 + B4 spatial Moran panel (Moran's I bar by stage + bivariate Moran)
- B6 pan-cancer THCA outlier panel
- B7 Hallmark GSEA bar
- **Already on disk and signed-off:** `SuppFig_X1_GSE250521_spatial_PT_to_ATC.{png,pdf}`, `SuppFig_X2_DM1_thyroid_nonoverlap_crossvalidation.{png,pdf}`, `SuppFig_X3_GSE248205_autoimmune_negative_control.{png,pdf}` (existing freeze; do not regenerate)

### Methods supplementary

- Closure-battery one-paragraph negative-feasibility note (per Paper 2 cleanup memo §5; equivalent text fits Paper 1 Methods Supplementary if author wants to pre-empt referee "have you tried H&E?" questions).

### Suppl tables

- B9 dark-matter subset NS table (honest negative)
- B10 Q3 spot-level TROP2 colocalization table (honest negative)
- B2 drop-one-out summary table

---

## 7. Venue language — required rewrite

`PAPER1_DM1_FULL` §0 and §3 contain probability claims:

| Original line | Status | Replacement (internal-memo only) |
|---|---|---|
| *"Sci Rep base → Cell Rep Med 90%, Nat Cancer 50-65%, Cancer Cell 20-30% 도전권"* | speculative | *"Submission strategy (internal): Sci Rep is the floor; JCI Insight / Cell Rep Med is the realistic ceiling for a computational-only manuscript without wet-lab validation. Reach venues (Nat Cancer, Cancer Cell) require functional validation that is outside the marathon scope."* |
| *"Sci Rep 100% / JCI Insight 99% / Cell Rep Med 90% / Nat Cancer 50-65% / Cancer Cell 20-30%"* (§3 table) | speculative | Same as above; do not use percentages. |
| *"G1 LOSO Spearman > 0.4 hit 시 Nat Cancer 비율 70%+"* (§11) | doubly speculative — premise impossible per closure | **Delete.** |
| *"Cell Rep Med 92% comfortable / Nat Cancer 55-65% reach / Cancer Cell 15-25%"* (§11) | speculative | **Delete.** |

**Submission-facing documents (cover letter, response-to-reviewers, abstract, brief, README, slide deck):** **never** include venue-probability percentages. Allowed wording: target venue named (singular), no probability number, no comparison ranking.

**Internal-only documents (this audit, planning memos, status notes):** rough qualitative bands acceptable (*"floor / realistic ceiling / reach"*), still no numeric percentages.

---

## 8. Final lock (post-audit state)

### 8.1 Continue (no change)

- **Molecular DM1 axis** as Paper 1 spine: cross-validation (TCGA bulk + ST), clinical (PFI/DFI Cox), mechanism (TF / GSEA / TF network), therapeutic (TROP2 + tumor-vs-normal) — all source files on disk, all claims independent of any pathology/image work.
- **Marathon mode 5/4–6/13:** manuscript writing only; new analyses paper-blocking only (per `v17_marathon_mode_post_pillar1`).
- **Voice-protected sections** owned by author keyboard (Hook ¶1, Aim ¶4, Discussion §3.1 with Landa 2016 cite save, §3.4 Limitations, Cover ¶1, Reviewer Q9). Per `v17_sprint_vs_marathon_violation`.
- **Existing scaffold** (`2026_05_03_manuscript_v8_OUTLINE.md` and siblings) and existing supplementary freeze (`project/supplementary/spatial_freeze_2026_05_03/`) remain authoritative for Paper 1 layout.
- **Closure of HM issues** per `2026_05_04_paper1_HM_closure_report.md` and §1 of post-closure recheck (H1/H2/M1/M2/M5/M6 closed; M3/M4 deferred bib lookups; L1–L9 cosmetic batch deferred to W6).

### 8.2 Dropped / blocked (no exception without explicit user authorization)

- **H&E → DM1 pathology model** (any embedder, any tile size, any cohort) — closure-battery final NO-GO.
- **TCGA WSI external validation** — Phase C deprecated; ~500 GB download blocked.
- **Image-DM1 IP / patent claim** — color-artifact unpatentability.
- **RunPod G1 / G2 / G3** — all three blocked. The package at `/data/spatial/phaseB_gpu_pkg.tar.gz` should not be re-dispatched. PID 168753 (separate-session embeddings_resnet50_448) line of inquiry already concluded NO-GO; no "complementary" result will materialize.
- **Foundation model UNI/CONCH/Virchow2 retry** — gated by 5 re-entry conditions of decision memo §5; none currently met.
- **Wet-lab roadmap** — post-bioRxiv only.
- **Paper 3 / Paper 4** — not touched (Paper 3 chmod 444 freeze preserved).

### 8.3 Action items (no compute)

1. **`PAPER1_DM1_FULL_2026_05_04.md` annotation (optional, not done by this audit):** if the user wants this document to remain on the top-level for reference, prepend a one-paragraph header noting that §9, §11 *"다음 step"* RunPod block, and venue percentages are **superseded by 2026-05-04 closure battery / final NO-GO / this audit**, and that the live Paper 1 plan is `2026_05_03_manuscript_v8_OUTLINE.md` + scaffold siblings. This is a single-file annotation (no manuscript prose touched). User authorization required.
2. **Submission-facing documents:** sweep for percentage-style venue claims (`grep -rIn -E "Cell Rep Med [0-9]+%|Nat Cancer [0-9]+%|Cancer Cell [0-9]+%" project/manuscript_p1_*` etc.). If found, delete or replace per §7. Defer to author keyboard if any sit in voice-protected sections.
3. **Manuscript marathon resumes:** `bib-m3m4` (M3/M4 bib lookups, ~20 min, single-file edit) or `voice-hook` (author keyboard). Per `2026_05_04_post_closure_recheck.md` §4. **No new analysis. No RunPod. No download.**

### 8.4 Audit attestation

- **No file modified** except this audit report (`project/reports/2026_05_04_paper1_dm1_full_conflict_audit.md`).
- **No analysis run.**
- **No RunPod / GPU dispatch.**
- **No new download.**
- **Paper 3 / Paper 4 not touched.**
- **Voice-protected manuscript sections not drafted.**
- All 12 result TSVs cited in §2 verified to exist on disk (`ls` check).
- Closure-battery NO-GO (`CLOSURE_BATTERY_2026_05_04.md` option D 폐기 + `2026_05_04_image_dm1_final_nogo_decision.md` final D adopted) is treated as authoritative; this audit does not relitigate it.

---

*Audit closed. Bundle reconciled. Paper 1 molecular axis = continue. Image-DM1 = stays dropped. RunPod G1/G2/G3 = stays blocked. Manuscript marathon resumes.*
