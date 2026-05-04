# Paper 2A — Re-entry protocol v0 (future-work design doc)

**Date:** 2026-05-04 · **Owner:** Seungho Cook
**Status:** v0 design — Claude-drafted scaffolding for a future image-DM1 re-entry attempt under the 5 conditions of `2026_05_04_image_dm1_final_nogo_decision.md` §5. **Not for execution during current marathon (5/4–6/13).** Reference document for post-marathon planning, grant proposals, and Bundang FFPE pilot scoping.

**Authority:** `2026_05_04_image_dm1_final_nogo_decision.md` §5 specifies the 5 re-entry conditions. This document defines what each condition would look like quantitatively, using GigaTIME (Valanarasu et al., 2026 Cell) as the multimodal-AI-baseline reference.

**No execution. No new analysis. No download. No GPU. No RunPod.** This is a planning document only.

---

## 1. The 5 re-entry conditions (verbatim from decision memo)

1. **Full-resolution scanner WSI available** for ≥ 1 thyroid cohort with paired bulk RNA-seq
2. **External molecular labels** (e.g., bulk 8-gene RNA score on the same WSI patients)
3. **Explicit post-marathon decision** (after Paper 1 bioRxiv submission, target 2026-06-13)
4. **Foundation model access** (UNI / CONCH / Virchow2) acquired and validated on a non-thyroid public benchmark first
5. **Pre-registered hypothesis** distinguishing what would constitute a positive result vs another color-artifact false positive

Current status: **0 / 5 met.**

---

## 2. Quantitative spec for each condition

### Condition 1 — Full-resolution scanner WSI

| metric | minimum | preferred | rationale |
|---|---|---|---|
| Patients (paired H&E + bulk RNA) | 100 | 1,000+ | 100 = LOSO statistical power for r ≥ 0.20 detection at α = 0.05; 1,000+ approaches GigaTIME per-cancer power |
| Resolution | 0.5 µm/pixel | 0.25 µm/pixel | nuclear morphometry boundary; current Visium hires ≈ 1.0 µm/pixel insufficient |
| Scanner | clinical-grade (Hamamatsu, Aperio, Leica) | any | scanner artifact stability matters for foundation models |
| Tissue blocks | FFPE primary tumor | + matched normal + recurrence | FFPE = pathology workflow standard; matched normal supports tumor-vs-normal validation |
| Cohort source | Bundang FFPE prospective (current outreach stage, n target ≥ 50 first) | + Korean multi-center retrospective | Bundang n > 50 is one of the 4 Paper 3 entry conditions; alignment opportunity |

Estimated time-to-meet: 6–18 months pending Bundang outreach delivery.

### Condition 2 — External molecular labels

| metric | minimum | preferred |
|---|---|---|
| RNA modality | bulk 8-gene RT-qPCR | bulk RNA-seq + ST validation subset |
| Label scope | DM1_like_score per case | + RAI_8 + TDS + other axes |
| Per-cohort RNA-H&E pairing | 100% of WSI cases | + temporal (pre-/post-treatment) |
| Independent cohort | ≥ 1 (TCGA-THCA already paired RNA + WSI but at frozen-section quality) | ≥ 2 (TCGA + Bundang + ideally another international) |

Note: TCGA-THCA already provides paired bulk RNA + WSI for ≈ 500 patients. The blocker has been download (~500 GB) explicitly forbidden under marathon mode (`2026_05_04_image_dm1_final_nogo_decision.md` §3). Post-marathon, a stratified subsample (n ≈ 50–100) is the realistic first cohort.

### Condition 3 — Post-marathon decision

| event | gating action |
|---|---|
| Paper 1 bioRxiv submission (target 2026-06-13) | unlocks reconsideration of image-DM1 |
| Paper 2B (Hashimoto-overlap molecular) submission | further unlocks scope discussion |
| Yu advisor review of Paper 1 outcome + funding alignment | required green light |
| Funding source for compute (~$1k–$10k/year) | required |

Estimated time-to-meet: late 2026 at earliest, contingent on Paper 1 publication trajectory.

### Condition 4 — Foundation model access

| asset | status | path to access |
|---|---|---|
| UNI (MahmoodLab) | gated, not requested | HuggingFace request, 1–7 day approval |
| CONCH (MahmoodLab) | gated, not requested | HuggingFace request, 1–7 day approval |
| Virchow2 (paige-ai) | gated, not requested | HuggingFace request + commercial-vs-academic check |
| Prov-GigaPath (Microsoft) | gated, not requested | HuggingFace request |
| GigaTIME (Microsoft) | not located publicly | wait for release (typical 6–12 months post-publication) |
| H-optimus-0 (bioptimus) | publicly available | direct HuggingFace |
| Hibou-B (HistAI) | publicly available | direct HuggingFace |

Validation requirement (per condition #4 spec): each candidate must be validated on a **non-thyroid public benchmark** (e.g., CAMELYON16/17 metastasis detection, TCGA-BRCA molecular subtype, BACH classification) before attempting our task. This prevents "embedder bottleneck" misattribution of negative results.

Estimated time-to-meet: 1–2 months for non-gated; 1–8 weeks for gated approval; pre-validation 1–2 weeks of GPU work.

### Condition 5 — Pre-registered hypothesis

The pre-registered hypothesis must distinguish:

| outcome | interpretation | gate |
|---|---|---|
| Spearman r ≥ 0.30 + AUROC ≥ 0.70 + controls weak + foundation > simple-RGB on residualized target | Real biological signal at full-res with foundation model — GO image-DM1 angle |
| Spearman r ∈ [0.20, 0.30) | BORDERLINE — extend to multi-cohort LODO before claiming |
| Spearman r ≥ 0.30 BUT simple-RGB also ≥ 0.20 on residualized | Color-artifact suspect — NO-GO unless stain normalization removes it |
| Spearman r < 0.20 OR AUROC < 0.65 | Confirmed NO-GO — close angle definitively |

Pre-registration document (template): would specify exact dataset, model, training protocol, evaluation metrics, gates, and stopping rules **before** any embedding or training. Filed publicly (e.g., OSF, GitHub release) before execution.

Estimated time-to-meet: 1 day to write, 1 week to formalize.

---

## 3. Pipeline architecture (post-conditions-met)

### Phase R0 — Pre-validation (per condition #4, no thyroid data)

```
inputs: CAMELYON16 (or chosen benchmark) WSIs + ground truth
compute: ~50 GPU-hours per foundation model
output: validated_embedders.json with per-model AUC vs published benchmark
gate: each embedder must reach within 2% of published benchmark AUC
```

### Phase R1 — Foundation embedding on existing closure-battery tiles

```
inputs: existing GSE250521 3,200 tiles 224 px (already on disk)
compute: ~5 GPU-hours per validated foundation model
output: embeddings_<model>_224.npz
gate: foundation Spearman r > ResNet50 r + 0.05 (improvement check); else stop
```

### Phase R2 — Foundation embedding on full-res Bundang WSI

```
inputs: Bundang FFPE n ≥ 100, paired RNA + 8-gene score
compute: ~50 GPU-hours per slide × n slides ≈ 5,000 GPU-hours
output: per-slide WSI tiles + per-tile embeddings
gate: pre-registered Spearman ≥ 0.30 OR AUROC ≥ 0.70
```

### Phase R3 — RNA-paired multimodal training (GigaTIME-class)

```
inputs: paired (H&E tile, RNA score) at scale (≥ 1k patients ideal, 100+ minimum)
compute: ~1,000+ GPU-hours
output: thyroid-specific multimodal foundation model
gate: external validation on TCGA-THCA WSI subset
```

### Phase R4 — External validation + claim lock

```
inputs: TCGA-THCA WSI ≥ 50, paired bulk RNA
compute: ~100 GPU-hours
output: external Spearman + AUROC with bootstrap CI
gate: ≥ 0.25 Spearman OR ≥ 0.65 AUROC + controls weak + biological coherence
```

### Estimated total cost (post-conditions, ballpark)

| line | wall time | compute cost | other |
|---|---|---|---|
| R0 pre-validation | 1 month | $1,000 | embedder access fees free for academic |
| R1 GSE250521 retry | 1 week | $200 | none |
| R2 Bundang full-res | 6 months | $5,000–$15,000 | scanner time, IRB, slides |
| R3 multimodal training | 6 months | $20,000–$50,000 | model training infrastructure |
| R4 external validation | 2 months | $1,000 | TCGA WSI download |
| **Total** | **12–18 months** | **$30,000–$70,000** | + Bundang cohort delivery + grant funding |

This is a **funded research program**, not a marathon-mode side activity.

---

## 4. Decision tree for re-entry triggering

```
After Paper 1 bioRxiv (2026-06-13)
  ↓
Funding secured (NIH R01 / Korean equivalent / industry partnership)?
  ├── No → defer indefinitely, retain as future-work in Discussion of subsequent papers
  └── Yes → Bundang FFPE n ≥ 50 delivered?
        ├── No → outreach acceleration; meanwhile Phase R0 with public benchmarks
        └── Yes → Yu advisor green light + pre-registration filed → execute R0–R4 sequentially
```

---

## 5. Comparison to GigaTIME's resource scale (reality check)

| dimension | GigaTIME (achieved) | Our R3 minimum (planned) | Ratio |
|---|---|---|---|
| Paired training cells | 40 × 10⁶ | 100 patients × ~10⁵ tiles/patient ≈ 10⁷ | ~4× lower |
| Patients | 14,256 | 100–1,000 | 14×–140× lower |
| Cancer types | 24 | 1 (thyroid) | thyroid-specialized |
| Hospitals | 51 | 1–2 | 25×–50× lower |
| Compute infrastructure | Microsoft Research / Providence Health | TBD academic | unknown |
| Time to publish | ~3 years (Providence cohort 2022 onset) | ~3 years estimate | comparable |

We do not need to match GigaTIME's 14k-patient scale — thyroid-specialized depth at 100–1,000 patients with high-quality paired RNA labels can deliver a thyroid molecular axis predictor that GigaTIME's general approach cannot. The framing is "specialized vs general", not "competing".

---

## 6. What this document does NOT do

- ❌ Authorize any execution during marathon (5/4–6/13)
- ❌ Commit to specific funding sources, cohorts, or timelines
- ❌ Replace or modify `2026_05_04_image_dm1_final_nogo_decision.md` §5 (5 conditions remain canonical)
- ❌ Justify any RunPod / Azure GPU spend during marathon
- ❌ Claim foundation models will work — only specifies what would constitute a fair test

## 7. What this document IS for

- ✅ Reference for grant proposals and Bundang FFPE pilot scoping
- ✅ Discussion §3 "Future direction" 1-sentence cite anchor (voice-protected, author keyboard)
- ✅ Reviewer Q response anchor (cross-ref `2026_05_04_paper2a_reviewer_q_gigatime_scaffold.md` Q1, Q2, Q7)
- ✅ Yu advisor briefing material for venue + future-work decisions
- ✅ Internal record of what re-entry would actually require, so the marathon-mode "stop" decision is durable

---

## Cross-references

- `2026_05_04_image_dm1_final_nogo_decision.md` §5 (5 re-entry conditions canonical)
- `2026_05_04_paper2a_gigatime_methods_comparison_table.md` (S1 comparison table)
- `2026_05_04_paper2a_reviewer_q_gigatime_scaffold.md` (S2 reviewer Q)
- `2026_05_04_gigatime_thca_relevance_scan.md` (S3 scan)
- `2026_05_04_paper1_dm1_full_molecular_only_lock.md` (manuscript-safe envelope)
- `2026_05_04_marathon_state_after_image_dm1_nogo.md` (per-paper marathon snapshot)
