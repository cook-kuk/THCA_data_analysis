# Paper 2 status — post image-DM1 NO-GO scope clean-up

**Date:** 2026-05-04 · **Owner:** Seungho Cook
**Trigger:** Closure battery final NO-GO (`2026_05_04_image_dm1_final_nogo_decision.md`).
**Scope of this memo:** scope/wording cleanup only. **No edits to voice-protected manuscript prose.**

---

## 1. What changes in Paper 2

| element | before | after |
|---|---|---|
| H&E → DM1 projection as Paper 2 active pillar | listed in some early planning (initial Phase A feasibility doc) | **dropped** |
| Pathology AI / image-DM1 triage as headline contribution | tentative | **dropped** |
| TCGA WSI Phase C external validation | gated, never run | **deprecated** — do not reference as planned |
| Multi-cohort LODO (GSE230424 + GSE248205) for image-DM1 | Phase B plan | **deprecated** for image-DM1 purpose. The two datasets remain available for *molecular-only* use under Paper 2's existing HT-overlap scope. |

## 2. What stays unchanged in Paper 2

Paper 2 = **Hashimoto-overlap PTC mechanism, molecular-only.** Per `STATUS_PAPER2_2026_05_02.md` and `v18_paper2_HT_isolated` memory.

| pillar / asset | status |
|---|---|
| Pillar I v2 forest (Korean PTC pool n=874 vs **AFND South Korea baseline**) | **active**, committed |
| Terminology audit (Task B, 5/4) | **done** |
| GSE286332 PTC vs PTC+HT (n=18) signature: 10,380 DEGs, HLA-II d=+3.65, IFN-γ FDR=2e-4 | **active** |
| TCGA Hashimoto-like generalization (GSE286332 sig → DM2 enrichment OR up to 5×, p=6e-10) | **active** |
| BCR clonal + TLS (D5-P6): TLS d=+1.96, IGHV clonality d>0.5, AICDA up | **active** |
| DM1 sub-B = NBNR cluster (D6-P7, n=56, 96% mutation-neg, 4× Hashimoto-like rate) | **active** |
| Pillar II / Pillar III scaffolding (TBD) | **active**, marathon target |
| Cross-paper reciprocal Discussion line (PTC ≠ GD) | **active** |

## 3. Forbidden in Paper 2 (additions to existing list)

Add the following to the Paper 2 forbidden-wording list (existing list in `STATUS_PAPER2_2026_05_02.md` covers GD-specific terms and conflations; this memo extends with image-DM1 deprecations):

- `H&E-DM1 projection` / `H&E DM1 prediction` / `image-DM1 triage`
- `pathology AI projection of 8-gene` / `morphology-derived DM1`
- `TCGA WSI external validation of DM1 program` (image-based)
- `tile-level DM1 inference`
- `H&E foundation model triage` (in active-claim form — only allowed as a stated negative-feasibility result)

If any of the above appear in current draft text, drafts, or briefs, replace with one of:
- (deletion is preferred when sentence is salvageable without it)
- "DM1 molecular subtype (RNA-seq-based, see Paper 1)"
- "We pre-tested H&E-based prediction of the depth-residualized DM1 axis at GSE250521 hires resolution and found no learnable signal beyond random panels (see closure battery, Methods Supplementary X)."

## 4. Files to audit (for image-DM1 wording)

**Author should grep / review these on next marathon block** (status memo only, no edits made by Claude here):

```
project/manuscript_p2_brief/                  # brief HTML/PDF, README, COHORT_ACCESS_GUIDE
project/reports/2026_05_03_paper2_pillar1_for_web_claude.md
project/STATUS_PAPER2_2026_05_02.md
project/results/p2_pillar1_forest_v2/PILLAR1_FOREST_V2_SUMMARY.md
PATHOLOGY_DM1_FEASIBILITY_2026_05_03.md       # initial plan — already superseded by closure
```

**Recommended audit command** (run by author at next marathon block):

```bash
grep -rIn -E "image[- ]?DM1|H&?E.{0,3}DM1|tile.level DM1|pathology[- ]AI projection|WSI external validation of DM1" \
    project/manuscript_p2_brief project/reports project/STATUS_PAPER*.md
```

For each hit: either delete, or convert to negative-feasibility framing per §3.

## 5. Methods supplementary slot

Suggested home for the closure battery as a **pre-registered failed hypothesis**:

- Paper 1 Methods Supplementary, sub-section: *"H&E predictability of the depth-residualized DM1 axis (negative feasibility)"* — strengthens the molecular-only argument by explicitly bounding what H&E can and cannot recover at this resolution.
- Paper 2 Discussion, one-line acknowledgement: *"We pre-tested but did not find a tile-level H&E correlate of the depth-residualized DM1 axis at hires Visium resolution; H&E-based triage of the molecular subtype is therefore not pursued in this work (see closure battery)."*

These are author-executed insertions during marathon writing. Claude does not draft voice-protected prose.

## 6. Allowed future re-entry to image-DM1 (mirrors decision memo §5)

Re-opening permitted **only** when:

1. Full-resolution scanner WSI for ≥1 cohort + paired bulk RNA-seq
2. External molecular labels (not ST surrogate)
3. Post-marathon explicit decision (after Paper 1 bioRxiv 6/13)
4. Foundation model access acquired + validated on non-thyroid benchmark
5. Pre-registered hypothesis for re-test

Until all five hold: no retry.

## 7. Marathon impact

- **Time freed:** ~4–6 weeks of would-be analysis (multi-cohort ingestion + TCGA WSI download + foundation model embedding + LODO/external validation)
- **Re-deploy to:** Pillar II/III scaffolding for Paper 2; Paper 1 voice-protected sections (author keyboard); cross-paper Discussion harmonization
- **Risk reduced:** removes 4–6 week dependency on uncertain GPU/HF access and 500 GB TCGA WSI download

## 8. Cross-references

- `2026_05_04_image_dm1_final_nogo_decision.md` (decision authority)
- `2026_05_04_marathon_state_after_image_dm1_nogo.md` (marathon-level snapshot)
- `CLOSURE_BATTERY_2026_05_04.md`
- `project/reports/pathology_dm1_closure_battery_2026_05_04.md`
- `STATUS_PAPER2_2026_05_02.md` (existing Paper 2 status — still authoritative for scope/forbidden words; this memo extends, does not replace)
