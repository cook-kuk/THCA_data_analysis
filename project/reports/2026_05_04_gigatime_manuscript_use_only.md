# GigaTIME — manuscript-use-only consolidation (Paper 1 utility)

**Date:** 2026-05-04 · **Owner:** Seungho Cook
**Scope:** Compress the GigaTIME (Valanarasu et al., 2026 Cell 189[2]:386, DOI 10.1016/j.cell.2025.11.016) bundle (commit `df54927`) into manuscript-use-only utility. **Methods supplement insert + cover letter author bullets + reviewer Q1–Q8 one-pager + forbidden-claims list.** Nothing else carries into Paper 1 manuscript or submission package.

**Position lock:** GigaTIME = **task-class boundary / future-work comparator only.** NOT competitor. NOT THCA-specific claim. NOT main-text material.

**Voice rule:** Hook ¶1 / Aim ¶4 / Discussion §3.1 / Limitations §3.4 / Cover Para 1 / Reviewer Q9 — author keyboard, NOT scaffolded here.

---

## A. Methods supplement insert (≤ 250 words, drop-in for `07_star_methods.md`)

```markdown
### Comparison to recent multimodal-AI baselines for H&E and tumor biology

To position our pre-registered negative result against the current
state of multimodal AI in tumor pathology, we note that recent
population-scale work (GigaTIME; Valanarasu et al., 2026, Cell
189[2]:386) trains a cross-modal H&E → virtual multiplex-IF
translator on ~40 × 10⁶ paired cells across 14,256 patients and
24 cancer types, recovering 1,234 protein–biomarker–staging–survival
associations and validating on 10,200 TCGA patients. That work
addresses a different task class (image-to-image protein-channel
translation given paired multimodal training at scale) than ours
(continuous regression of a depth-residualized 8-gene transcriptional
axis from frozen ImageNet feature transfer on 3,200 spot-aligned
hires Visium tiles from 16 thyroid slides). The two studies are
complementary rather than competing: GigaTIME demonstrates that H&E
carries substantial signal for the tumor immune microenvironment when
paired with molecular ground-truth at the order of 10⁷ cells, whereas
our pre-registered test bounds what is recoverable from H&E alone
for a thyroid transcriptional axis at our scale and resolution.
Bridging the two — RNA-paired H&E training on a thyroid cohort
approaching ≥ 10² patients with full-resolution scanner WSI and
matched 8-gene RNA labels — is identified as the principal future-
work direction for any image-based triage of the molecular subtype,
and is conditioned on the five re-entry criteria specified separately.
Cancer-type and protein-panel composition of GigaTIME's 24 cancers /
21 proteins are not enumerated in publicly accessible metadata as of
this writing; we therefore do not claim direct THCA-specific
comparison.
```

(247 words. Cite as `\citep{Valanarasu2026}`. BibTeX entry already in `project/manuscript_v8/03_intro_references.bib` via commit `2ffe929`.)

---

## B. Cover letter Para 2–3 — author bullet menu (NOT prose)

**Voice rule:** Para 1 = author keyboard. Para 2–3 = author writes prose using bullets below. **No Claude prose.**

### Para 2 — molecular axis contribution (bullet menu)

- 8-gene RAI / DM1 transcriptional axis, driver-excluded curated pool
- multi-cohort RNA + ST validation (TCGA n = 561, 28 ST slides, 12-slide bootstrap)
- independent prognostic effect after BRAF/RAS/age/stage adjustment (DFI multivariate HR = 1.41, p = 0.025; n = 461)
- 4-step mechanism framework: lineage TF collapse → DNMT activation → STAT3/AP-1 → TROP2 re-expression
- TROP2 = tumor-population therapeutic vulnerability (NOT spot-level DM1 colocalization; Q3 NEG ρ = −0.016)
- subtractive honesty embedded: N1 dark-matter NS, Q3 spot-level NEG, H&E NO-GO acknowledged
- functional validation explicitly noted as required (TROP2 IHC, sacituzumab IC50, organoid)

### Para 3 — orthogonal complementarity to multimodal AI frontier (bullet menu)

- recent multimodal AI (GigaTIME; Valanarasu et al., 2026 Cell 189[2]:386) translates H&E → virtual mIF protein channels at population scale (40 × 10⁶ paired cells, 14,256 patients, 24 cancer types)
- different task class (image-to-image protein translation) at different scale (40 × 10⁶ paired vs our 3.2 × 10³ tiles)
- our 8-gene RNA axis = molecular ground-truth axis that any thyroid-specific multimodal training would target
- pre-registered negative feasibility (closure battery) bounds H&E-alone recoverability at our scale + resolution
- five re-entry conditions specify when image-based triage could be productively re-opened (full-res WSI + paired molecular labels + foundation-model access + post-marathon decision + pre-registered hypothesis)
- positioning: complementary lanes, not competing claim

### Forbidden in Para 2 / Para 3

- ❌ venue-probability percentages (any form)
- ❌ "Cancer Cell-ready" / "Nat Cancer reach" / aspiration tier labels
- ❌ "GigaTIME does X for thyroid" / "we replicate GigaTIME for THCA"
- ❌ "image-DM1 / pathology-AI Phase B/C" in any active form
- ❌ "DM1-high spots are TROP2-high" / spot-level TROP2 colocalization
- ❌ "proves" / "establishes" / "all risks resolved"
- ❌ Hook-style rhetorical openers (Para 1 carries voice)

---

## C. Reviewer Q1–Q8 one-page (Q9 = voice-protected, EXCLUDED here)

| Q | likely topic | answer key (≤ 2 sentences) |
|---|---|---|
| **Q1** | Why frozen ResNet50 ImageNet, not UNI/CONCH/Virchow2? | Pre-registered gate (Spearman ≥ 0.30 GO; ≥ 0.20 BORDERLINE). Published foundation-model gain over ResNet50 on tile-level molecular tasks ≈ +0.05 to +0.15; even maximal gain lands at ≈ 0.17, below BORDERLINE. Foundation retry retained as re-entry condition #4 (access + non-thyroid pre-validation required). |
| **Q2** | Why didn't you do GigaTIME-style multimodal AI? | Different task class (image-to-image protein translation vs continuous transcriptional regression) at different scale (40 × 10⁶ paired cells vs our 3,200 tiles). Our 8-gene RNA axis supplies the molecular ground-truth a future thyroid-specific multimodal training would target; the studies are complementary, not competing. |
| **Q3** | Was GSE250521 hires Visium resolution sufficient? | Hires PNG is downsampled by tissue_hires_scalef ≈ 0.56 from full-res; our 224 px tile ≈ 110 µm in tissue space — adequate for tissue architecture, insufficient for nuclear morphometry. Resolution explicitly marked as re-entry condition #1 (full-resolution scanner WSI required). |
| **Q4** | Could slide-level aggregation rescue the spot-level negative? | We tested mean / top-25%-mean / high-risk fraction; only top-25%-mean reached marginal r = +0.209 at n = 16 slides — too few for significance and consistent with per-fold variance. Reported as exploratory; not used to flip GO/NO-GO verdict. |
| **Q5** | Were other tile sizes tested? | 448 px gave Spearman = 0.017, slightly worse than 224 px (0.022). Pre-registered Δr ≥ +0.05 improvement gate not met; 672 px embedding cancelled per pre-registration after tile extraction (saving ≈ 25 min CPU). Larger context does not bridge the floor. |
| **Q6** | Why does a 9-dim RGB baseline beat 2,048-dim ResNet50 on raw target? | RGB features (mean, SD, entropy per channel) capture tile-level color statistics — a stain-darkness ↔ sequencing-depth artifact (independently confirmed by DM1_raw ~ log_counts ρ = −0.45). On the depth-residualized target, both fall to near-zero (RGB 0.080, ResNet50 0.022), validating the residualization design rather than refuting it. |
| **Q7** | Could you have used GigaTIME's released model directly? | GigaTIME model weights and code repository have not been located in publicly accessible search as of submission (GitHub microsoft/GigaTIME → 404; HuggingFace search → none). Microsoft Research typically releases code 6–12 months post-publication; the same closure-battery pipeline can be re-applied if/when assets become available. |
| **Q8** | Does ATC vs non-ATC AUROC = 0.689 contradict the H&E negative? | No. ATC's gross anaplastic morphology is pathologically distinct and well within ImageNet-feature-transfer reach; this captures gross histology, not within-tumor molecular gradients. Continuous molecular targets (stage ordinal Spearman = −0.047; epithelial / proliferation r ≤ 0.03; depth-residualized DM1 r = 0.022) all fail; the two findings are mutually consistent. |

**Q9 — voice-protected (author keyboard).** Reference: `2026_05_04_voice_hook_fact_brief.md` §C dark-matter content. Closure-battery negative is evidence of pre-registered honest hypothesis testing.

---

## D. GigaTIME-related forbidden claims list (additions to existing forbidden envelope)

In addition to the existing manuscript-safe envelope (`2026_05_04_paper1_dm1_full_molecular_only_lock.md` §3, §5), the following GigaTIME-specific claims are **explicitly forbidden** in any Paper 1 submission-facing text (manuscript, supplement, cover letter, response to reviewers):

| forbidden claim | reason |
|---|---|
| "GigaTIME includes thyroid cancer" / "we cite GigaTIME's THCA result" | THCA inclusion in their 24 cancers NOT confirmed in public metadata |
| "GigaTIME's panel includes TROP2 / FOXE1 / NKX2-1 / NIS / [our axis genes]" | 21-protein panel composition NOT confirmed in public metadata |
| "GigaTIME outperforms ResNet50 on transcriptional regression" | different task; comparison undefined |
| "We replicate / extend / reproduce GigaTIME for thyroid" | task class + scale + cohort all different |
| "GigaTIME model / weights are publicly available" | NOT located as of 2026-05-04 search |
| "Our 8-gene axis = the GigaTIME-equivalent for thyroid" | overstates; different methodology and intent |
| "GigaTIME validates our H&E negative result" | NOT a direct test; we cite as task-boundary context only |
| "Our work is in the GigaTIME class / lane" | lane distinction (RNA-axis vs image-translation) is the framing |
| Citing specific GigaTIME AUC / accuracy / per-cancer numbers | not located in available metadata; use only abstract numbers (40M cells, 14,256 patients, 24 cancers, 1,234 associations) |
| Using GigaTIME publication as venue-probability argument in submission text | venue arguments live in internal Yu-advisor memo only, NEVER submission-facing |
| "Image-DM1 angle is supported by GigaTIME's success" | image-DM1 = D 폐기; GigaTIME = different task; cannot transitively rescue |
| Recommending GigaTIME-style retry without 5/5 re-entry conditions met | conditions are canonical; do not bypass in any text |

---

## E. Where each utility lands (target file map)

| utility | target | who acts |
|---|---|---|
| §A 250-word Methods supp insert | `project/manuscript_v8/07_star_methods.md` (new sub-section under STAR Methods supplementary) | author copies + light edit |
| §B Para 2 / Para 3 bullet menu | `project/manuscript_v8/08_cover_letter.md` (Para 2 + Para 3 author writes prose; Para 1 voice-protected) | author writes prose using bullets |
| §C Q1–Q8 one-page | `project/manuscript_v8/09_reviewer_qa.md` (Q9 author keyboard) | author copies + refines |
| §D forbidden claims | reference card; consulted whenever GigaTIME wording appears in any draft | author + Claude both check |

---

## F. Cross-references (single source of truth)

- Bib: `project/manuscript_v8/03_intro_references.bib` `@Valanarasu2026` (commit `2ffe929`)
- Decision: `project/reports/2026_05_04_image_dm1_final_nogo_decision.md` §5 (5 re-entry conditions)
- Lock: `project/reports/2026_05_04_paper1_dm1_full_molecular_only_lock.md` (manuscript-safe envelope)
- Bundle (full long-form): commit `df54927`, 6 docs in `project/reports/2026_05_04_*gigatime*.md` and `project/reports/2026_05_04_paper2a_*.md`
- Hook fact brief: `project/reports/2026_05_04_voice_hook_fact_brief.md` (uncommitted, author-ready)

---

GigaTIME positioned as future-work comparator only. Return to voice-hook.
