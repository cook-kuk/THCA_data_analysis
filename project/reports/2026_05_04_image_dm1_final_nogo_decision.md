# Final NO-GO decision — H&E → DM1 image angle

**Date:** 2026-05-04 · **Owner:** Seungho Cook
**Status:** **CLOSED** — image-DM1 / pathology-AI angle officially dropped from active work.
**Authority:** This memo + `CLOSURE_BATTERY_2026_05_04.md` + `project/reports/pathology_dm1_closure_battery_2026_05_04.md` + `project/reports/pathology_dm1_phaseA_cpu_verdict_2026_05_04.md`.

---

## 1. Final decision

**D — Full 폐기.**

Original hypothesis (Paper 1's depth-residualized 8-gene DM1_like_score_resid is predictable from spot-aligned H&E tiles) is **rejected** under all tested model classes, target preprocessings, tile sizes, and aggregation strategies.

Decision-criteria walkthrough (per spec):

| option | gate | result |
|---|---|---|
| A. keep original image-DM1 | DM1_resid Spearman ≥ 0.20 OR AUROC ≥ 0.65 + beats controls + not color/depth artifact | **✗** best 0.080 (simple RGB on resid), all settings well below |
| B. pivot to raw thyroid-lineage morphology | raw Spearman ≥ 0.30 OR AUROC ≥ 0.70 + resid weak + controls weak | **✗** raw best 0.224 / AUROC 0.702 from RGB color stats, but driven by stain darkness ↔ sequencing-depth artifact (per spec, **not a valid scientific claim**). ResNet50 raw alone is 0.061. |
| C. pivot to gross histology / stage | stage or ATC strong, molecular weak | **△** ATC AUROC 0.689 (BORDERLINE), stage ordinal Spearman −0.047. ATC anaplastic morphology classification has dense prior art and is not a novel angle. |
| **D. full 폐기** | all molecular r < 0.15 AND AUROC < 0.60 + stage weak + controls similar | **✓** all conditions met |
| E. foundation model retry | ResNet raw or stage shows signal, OR DM1_resid r ≥ 0.10 | **✗** raw 0.061, resid 0.022 — spec says "otherwise no" |

→ **D adopted. E withheld** (future-only, not during current 5/4–6/13 manuscript marathon).

## 2. Why — five evidence pillars

1. **DM1_resid LOSO 224 px ResNet50 ImageNet** — pooled Spearman r = **0.022**, AUROC = **0.511**. At chance.
2. **Real signal indistinguishable from random.** 30 random 8-gene panels (resid mode): mean r = +0.002, max = +0.031. Real DM1_resid (0.022) sits below random max. 50 random panels (raw mode): mean r = +0.019, max = +0.062 vs real RAI_8 raw 0.067 — also at random max.
3. **Raw signal is stain ↔ depth artifact.** 9-dim RGB mean/std/entropy baseline beats 2048-dim ResNet50 features on raw DM1 (0.224 vs 0.061). Only signal that exists is H&E color depth, which residualization (correctly) removes. Residualization is **not** the failure mode — it removes a non-biological confound.
4. **All rescue paths failed.** 5 residualization variants (B1–B5: log_counts only, log_ngenes only, no resid, rank-norm, epi-top-50%) all r ≤ 0.060. Tile size 448 px gave Spearman 0.017 (slightly **worse** than 224's 0.022) — gate (+0.05) not met, 672 cancelled. Slide-level top25%-mean DM1_resid = +0.209 (n=16, marginal, not significant).
5. **Foundation model expected gain insufficient.** Published UNI/CONCH/Virchow2 over ResNet50 ImageNet typically delivers +0.05–0.15 Spearman on tile-level molecular tasks. Gap to GO threshold (Spearman ≥ 0.30) is +0.28. Even if foundation gains the maximum +0.15, the result lands at ~0.17 — still below the BORDERLINE 0.20 floor.

Dominant failure cause among five hypotheses tested:

| hypothesis | ranking |
|---|---|
| (1) ResNet50 ImageNet too weak | partial (9-dim RGB outperforms on raw) |
| (2) 224 px tile too small | **REJECTED** (448 worse, not better) |
| (3) depth-residualization removes morphology-visible signal | partial **but justified** — removed signal was stain artifact, not biology |
| (4) **spot-level ST molecular labels fundamentally mismatched with tile-level H&E** | **dominant** (stage ordinal r = −0.05, all continuous molecular targets r ≤ 0.08, only ATC gross histology works) |
| (5) hires PNG resolution / single cohort | partial (fullres + multi-cohort would not flip the verdict given (4)) |

## 3. What is dropped

| item | status |
|---|---|
| image-DM1 triage as Paper 2 active pillar | **dropped** |
| pathology projection wording in Paper 2 manuscript / brief | mark **deprecated** (separate status memo) |
| TCGA WSI Phase C (~500 GB download) | **never executed**, never to execute under current marathon |
| multi-cohort LODO Phase B (GSE230424 + GSE248205) | **dropped** as image-DM1 path |
| RunPod / Azure GPU spend on image-DM1 retries | **stopped** ($30 RunPod balance preserved for unrelated future work) |
| HF UNI/CONCH/Virchow2 access acceptance — pathology embedder | **withheld** (no application during marathon) |
| Paper 2 / IP claim based on H&E predictability of DM1 | **dropped** |

## 4. What remains (unaffected)

| item | status |
|---|---|
| Paper 1 DM1 molecular subtype (8-gene driver-excluded sub-stratifier) | **active**, unaffected |
| Paper 1 ST supplementary (GSE250521 spatial validation as supplementary figure) | **active**, unaffected |
| Paper 1 bioRxiv 6/13 target | **active** |
| Paper 2 HT-overlap PTC molecular-only scope (Pillar I v2 + Pillar II/III TBD) | **active**, unaffected by image-DM1 closure |
| Paper 2 GSE286332 PTC vs PTC+HT signature + TCGA Hashimoto-like generalization | **active**, unaffected |
| Paper 2 BCR clonal + TLS (D5-P6) + DM1 sub-B NBNR (D6-P7) findings | **active**, unaffected |
| Paper 3 Korean GD HLA / pan-Asian | **frozen** (per v19 gating; no touch) |
| Marathon 5/4–6/13 manuscript writing schedule | **active**, image-DM1 removal frees ~4–6 weeks of would-be analysis time |

## 5. Allowed future re-entry conditions

Re-opening this angle is **only** permitted when ALL of the following are true:

1. **Full-resolution scanner WSI available** for ≥1 thyroid cohort with paired bulk RNA-seq (e.g., Bundang FFPE prospective cohort once delivered) — Visium hires PNG (0.17× downsample) is insufficient.
2. **External molecular labels** (e.g., bulk 8-gene RNA score on the same WSI patients) for honest external validation, not the tile-level ST surrogate.
3. **Explicit post-marathon decision** — the current 5/4–6/13 6-week sprint is for manuscript writing only. Re-entry must wait until at minimum Paper 1 bioRxiv submission (target 2026-06-13) is in.
4. **Foundation model access** (UNI / CONCH / Virchow2) acquired and validated on a non-thyroid public benchmark first, to ensure the embedder is not the bottleneck before attempting our task.
5. **Pre-registered hypothesis** distinguishing what would constitute a positive result vs another color-artifact false positive.

If any of (1)–(5) is false, **do not retry**.

## 6. Cross-references

- `CLOSURE_BATTERY_2026_05_04.md` (top-level Korean summary)
- `project/reports/pathology_dm1_closure_battery_2026_05_04.md` (long-form English [1]–[10] sections)
- `project/reports/pathology_dm1_phaseA_cpu_verdict_2026_05_04.md` (Phase A 224-only NO-GO, pre-closure)
- `PHASE_A_NOGO_REPORT_2026_05_04.md` (top-level Phase A summary)
- `PATHOLOGY_DM1_FEASIBILITY_2026_05_03.md` (initial feasibility plan)
- `project/results/03_pathology_poc/closure_battery_metrics.tsv` (26 rows, A/B/D/E/F numbers)
- `project/results/03_pathology_poc/negative_controls_summary.tsv` (resid mode controls)
- `project/results/03_pathology_poc/negative_controls_raw_summary.tsv` (raw mode controls)
- `project/results/03_pathology_poc/loso_metrics_resnet50.tsv` (224 + 448 LOSO)
- `project/reports/2026_05_04_paper2_post_image_dm1_nogo_status.md` (Paper 2 scope cleanup)
- `project/reports/2026_05_04_marathon_state_after_image_dm1_nogo.md` (marathon snapshot)

---

*Decision is final under current marathon (5/4–6/13). No further compute, no further pod spin-up, no further cohort downloads on this angle. Re-open only via the 5 conditions in §5.*
