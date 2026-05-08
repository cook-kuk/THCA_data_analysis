# Sprint plan — Foundation-model + CLAM Image-DM1 v2 (full sprint)

**Date:** 2026-05-07
**Owner:** Seungho Cook
**Sprint name:** `p2_image_dm1_v2_foundation_clam`
**Mode:** **Full sprint executed** (Phase 0 + 1 + 2 + 3)
**Goal:** Either (a) strengthen Paper 1 spatial supp with rigorous pathology layer, or (b) launch Paper 2 image-DM1 sprint, depending on Phase 2 kill-switch.

---

## 0. NO-GO closure unlock — explicit record

The 2026-05-04 closure (`pathology_dm1_phaseA_cpu_verdict_2026_05_04.md` + `image_dm1_final_nogo_decision.md`) is hereby **explicitly re-opened** by user authorization 2026-05-07 with the following justification:

| Aspect | Closure 2026-05-04 (NO-GO) | This sprint v2 (re-opened) |
|--------|----------------------------|----------------------------|
| Architecture | ResNet50 (ImageNet-pretrained, generic) | **UNI / CONCH** (Mahmood lab; pretrained on millions of WSI from 100+ tissue types) |
| Pooling | Tile-mean / simple aggregation | **CLAM-style attention MIL** (Lu 2021 Nat BME) |
| Cell-level | None | **SAM segmentation** per spot |
| Foundation | none | UNI = Vision Transformer Huge / CONCH = vision-language multimodal |
| Embedding | 2048-d ResNet50 | 1024-d UNI or 512-d CONCH |
| Training | tile-level direct | bag-level attention MIL |

→ This is **architecturally distinct** from closure attempt. Closure ResNet50 baseline is retained as **comparison anchor** (ablation control), not as the main model.

**Audit trail:** This file + memory entry `image_dm1_final_nogo_decision` should be updated with v2-unlock pointer after Phase 3 results.

---

## 1. Sprint scope

**Goal options (choose at Phase 3):**
- (a) Paper 1 spatial supplement strengthening — descriptive + foundation-model embedding overlay on GSE250521
- (b) Paper 2 image-DM1 launching — TCGA WSI CLAM model with held-out validation
- (c) Both — PASS at Phase 2 kill-switch unlocks (b) automatically

**Out-of-scope (preserved boundaries):**
- ❌ Voice-protected manuscript prose (Hook/Aim/Disc §3.1/Limitations/Cover/Q9)
- ❌ Paper 3 (ICI) / Paper 4 (Korean GD HLA) territory
- ❌ Bundang Korean FFPE outreach (parallel; not blocking this sprint)
- ❌ TROP2 main claim restoration
- ❌ ICI / immunotherapy / HLA work

---

## 2. Phase plan + kill switches

### Phase 0 — Setup (now, no GPU, no cost)
- [x] Sprint workspace created: `project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/`
- [ ] This SPRINT_PLAN.md
- [ ] Inventory existing assets in detail
- [ ] Phase 1 dispatch package (Python scripts ready for RunPod)
- [ ] Phase 2 dispatch package
- [ ] Phase 3 integration template
- [ ] Cost cap acknowledged

### Phase 1 — GSE250521 deep-dive (1× A100, 1–2 hr, ~$3-5)
**Input:** existing 3.4 GB H&E tiles + per-spot RAI_8/DM1 scores
**Steps:**
1. UNI foundation model embedding (1024-d) on each tile (~30k tiles × 16 slides)
2. SAM cell segmentation per spot (Visium 55 µm spots)
3. Per-cell DM1 score density mapping (interpolate spot-level score to cell-level)
4. Within-slide UNI-embedding ↔ DM1 spatial correlation
5. Output: spatial figure + correlation report

**Kill switch (decision at end of Phase 1):**
- **PASS:** within-slide UNI-embedding 1st PC ↔ DM1 score Spearman ρ > 0.3 in ≥ 50% of slides → proceed to Phase 2
- **STOP:** ρ ≤ 0.2 in ≥ 75% slides → fundamental image-DM1 signal absent at foundation-model resolution; abort Phase 2; document as honest negative for Paper 1 supp

### Phase 2 — TCGA-THCA WSI CLAM (1-2× A100, 4-8 hr, ~$15-25)
**Input:** TCGA-THCA diagnostic WSI selective subset (60-100 slides; balanced DM1/DM2/not_DM)
**Steps:**
1. GDC public download (NOT controlled — diagnostic WSI is open access)
2. Tile extraction at 20× / 256-px (CLAM standard)
3. UNI / CONCH foundation embedding per tile
4. CLAM attention-MIL training: DM1 vs DM2 binary classification
5. LOSO + held-out 20% slide-level validation
6. Comparison with closure ResNet50 baseline (ablation)

**Kill switch (decision at end of Phase 2):**
- **PASS:** held-out AUC > 0.70 (closure ResNet50 baseline ~0.55; CLAM-foundation aim ≥ 0.70) → unlock Paper 2 launch path
- **MARGINAL:** 0.60 < AUC ≤ 0.70 → Paper 1 supp only (limited utility)
- **FAIL:** AUC ≤ 0.60 → reaffirm closure; document v2 negative; STOP

### Phase 3 — Results integration (no GPU, 1-2 hr)
**Outputs depending on Phase 2 result:**
- PASS:
  - Paper 1 supp updated with strong spatial + WSI evidence
  - Paper 2 first-pass figure draft (image_dm1_v2 panel set)
  - audit page Future-validation board updated
  - SPRINT_RESULT.md PASS report
  - memory update: `v19_paper2_image_dm1_v2_PASS`
- MARGINAL:
  - Paper 1 spatial supp updated only
  - Paper 2 launch decision deferred
- FAIL:
  - Paper 1 supp gets honest negative figure
  - closure re-confirmed with v2 architecture caveat
  - memory update: `v19_paper2_image_dm1_v2_NEGATIVE`

---

## 3. Cost cap

| Item | Estimate |
|------|----------|
| Phase 1 RunPod (1× A100, 2 hr) | $4 |
| Phase 2 RunPod (1× A100, 8 hr) or (4× A100, 2 hr) | $15-25 |
| TCGA WSI download bandwidth | ~$0 (Azure/AWS within egress limits) |
| Foundation model weights (UNI/CONCH) | $0 (free download) |
| **Total cap** | **$50** (sufficient cushion) |

**Hard cost stop:** if RunPod cumulative bill > $40 before Phase 3 completion, abort and report status.

---

## 4. Out-of-scope safeguards (recap)

| Boundary | Status |
|----------|--------|
| Voice-protected manuscript prose | NOT touched in this sprint |
| Hook / Aim / Discussion / Limitations / Cover / Q9 | NOT written |
| Paper 3 / 4 / 9 territory | NOT touched |
| TROP2 main claim | NOT restored |
| RAI predictor / treatment selection / clinical utility claim | NOT made |
| Bundang Korean FFPE outreach | parallel; not blocking |
| Closure ResNet50 retry (same architecture) | NOT — using foundation model + CLAM |

---

## 5. Files / outputs

```
project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/
├── SPRINT_PLAN.md                       # this file
├── scripts/                             # locally developed
│   ├── phase1_uni_embed.py             # UNI on GSE250521
│   ├── phase1_sam_segmentation.py      # SAM cell-level
│   ├── phase1_overlay_correlation.py   # DM1 score ↔ UNI PC
│   ├── phase2_gdc_download.py          # TCGA WSI download
│   ├── phase2_clam_train.py            # CLAM attention-MIL
│   └── phase2_clam_eval.py             # held-out + LOSO eval
├── phase1_gse250521/
│   ├── uni_embeddings.npz
│   ├── sam_segmentations/
│   ├── correlation_per_slide.tsv
│   ├── spatial_overlay_figs/
│   └── PHASE1_REPORT.md
├── phase2_tcga_clam/
│   ├── slide_manifest.tsv
│   ├── tiles_per_slide/
│   ├── clam_train_log.txt
│   ├── clam_results.tsv
│   └── PHASE2_REPORT.md
├── phase3_integration/
│   ├── SPRINT_RESULT.md                # PASS/MARGINAL/FAIL verdict
│   ├── paper1_supp_update_proposal.md
│   ├── paper2_launch_proposal.md (if PASS)
│   └── memory_update_proposal.md
├── reports/
│   └── (mirrored to project/reports/ at end)
└── runpod_dispatch/
    ├── DISPATCH_README.md              # SSH steps
    ├── env_setup.sh                    # CUDA + UNI + CLAM install
    └── run_all.sh                      # one-shot orchestrator
```

---

## 6. Next steps (now)

1. ~~SPRINT_PLAN.md~~ ✅
2. Phase 0 inventory + scripts (next)
3. Dispatch package (RunPod commands ready to copy-paste)
4. User SSH's RunPod + runs `run_all.sh` (or phase by phase)
5. Phase 3 integration locally

---

**Sprint authorized: 2026-05-07.** **NO-GO closure re-opened with foundation-model + CLAM architecture justification.**
