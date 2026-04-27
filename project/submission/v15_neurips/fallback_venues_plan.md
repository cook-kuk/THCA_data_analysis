# v15 — Fallback Venue Plan (if NeurIPS 2026 rejects)

_Generated 2026-04-27. Sequencing for re-submission targets if
NeurIPS 2026 main track rejects in September 2026._

## Context

NeurIPS 2026 main-track acceptance rate is historically ~25 %.
Self-review (`self_review_simulation.md`) predicts ~30 % accept
probability for v15 specifically — borderline. **70 % chance of
rejection**, in which case the paper needs a fallback target.

This document plans the next 3 venues in priority order.

---

## Tier 1 — TMLR (Transactions on ML Research)

**Target:** open-access ML journal. No deadlines (rolling
submission). Decision in 8–12 weeks. Publishes accepted papers
immediately on acceptance with no camera-ready window.

**Why TMLR is the right next stop:**
- No deadline pressure → submit any time after NeurIPS rejection
  email lands
- Accepts theory + empirical content equally (DIAL fits)
- Published-by-acceptance means the paper is citable as soon as
  it lands; useful for academic record
- TMLR has explicit policy that NeurIPS/ICML/ICLR rejections are
  not held against the paper
- Standard rebuttal format

**Adaptation effort:** ~10 hours
- Use NeurIPS reviewer feedback to strengthen weakest points
  (probably empirical: add Camelyon17 if R2-style critique was
  the cause of rejection)
- Re-format from NeurIPS double-column (or single-col + 9 page)
  to TMLR style (LMCS-format LaTeX class, no page limit)
- Re-write abstract to journal style (slightly less compressed)
- Move appendix A proof into the main body (no page limit)

**Submission timing:** 1-2 weeks after NeurIPS rejection (mid-
September 2026) → first decision early-mid December 2026. If
revisions requested, second-round decision ~February 2027.

---

## Tier 2 — ICLR 2027 (Jan/Feb 2027 deadline expected)

**Target:** ICLR main track. Deadline historically late September /
early October for ICLR. ICLR 2027 abstract deadline ~Sep 2026
(coincident with NeurIPS decision date), full paper ~Oct 2026.

**Why ICLR is a strong fallback:**
- Comparable prestige to NeurIPS, slightly broader theory acceptance
- ICLR open-review process means rebuttals are public; advantageous
  for honest-disclosure papers like v15 (the v5.2 retraction story
  reads better when fully public)
- ICLR 2027 likely accepts NeurIPS-formatted PDFs with minor edits

**Adaptation effort:** ~8 hours (less than TMLR because format is
similar)
- Port `\documentclass{neurips_2026}` → `\documentclass{iclr2027}`
- Add ICLR-specific reproducibility statement (similar to NeurIPS
  checklist; usually 90 % overlap)
- Address top 1-2 reviewer concerns from NeurIPS rebuttal that
  predicted Tier 1 rejection

**Submission timing:** if NeurIPS rejects mid-September 2026,
~10 days to adapt → ICLR submit by Oct 2026. ICLR decisions are
usually early Jan 2027.

**Tradeoff vs TMLR:** ICLR is venue-prestige equivalent but cycle
is 4 months. TMLR is faster (8 weeks) but less venue-prestige.
**Recommendation:** if NeurIPS reviewer feedback is constructive
(R1 + R3 positive, only R2 critical), submit ICLR. If feedback
suggests deeper rewrites (R1 + R2 both critical), submit TMLR
with the rewrites.

---

## Tier 3 — AISTATS 2027 (Oct 2026 deadline expected)

**Target:** Statistical-ML focused conference; theory-leaning.
AISTATS 2027 deadline historically October 2026. Stats-tilted
reviewer pool likely receptive to Theorem 1.

**Why AISTATS:**
- Best fit for "calibration paper with closed-form theorem" framing
- Reviewer pool has higher tolerance for synthetic-only empirics
- Lower-prestige than NeurIPS but legitimate venue
- Same paper, lighter empirical demand

**Adaptation effort:** ~12 hours
- Reformat to AISTATS style (slightly different from NeurIPS)
- Strengthen Theorem 1 to the level AISTATS reviewers expect
  (potentially adding the deferred Conjecture 1 proof if possible)
- Tighten the empirical section (AISTATS sometimes prefers fewer,
  more rigorous experiments)

**Submission timing:** if Tier 1 (TMLR) takes too long or returns
rejection, AISTATS in Oct 2026 → decisions Jan-Feb 2027.

---

## Tier 4 — Workshop tracks (NeurIPS / ICML workshops)

**Target:** AdvML, OOD-Robustness, ML for Health workshops at
NeurIPS 2026 (workshops happen alongside the main conference) or
ICML 2027 workshops.

**Why workshops:**
- Higher acceptance rate (40-60 %)
- 4-page extended-abstract format → trim main paper, retain
  experiments
- No camera-ready costs
- Useful as "publication before journal submission" if Tier 1-3
  all reject

**Adaptation effort:** 6 hours (trim to 4 pages + restructure for
workshop format)

**Tradeoff:** workshop publication has lower prestige than main
conference; may be perceived as "lesser publication" on academic
CV. Use only if Tiers 1–3 fully reject.

---

## Decision tree

```
NeurIPS 2026 reject (Sep 2026)
│
├── if reviewer concerns are LIGHT (~minor revisions) →
│       ICLR 2027 (Oct 2026 deadline; similar format, ~8 h adapt)
│
├── if reviewer concerns are MEDIUM (need real-data benchmark) →
│       implement Camelyon17 + DomainNet (CR-2 from camera-ready plan,
│       ~18 h work) → TMLR (Tier 1, late Sep 2026; rolling)
│
├── if reviewer concerns are HEAVY (Theorem 1 attacked) →
│       AISTATS 2027 (Tier 3, Oct 2026; stats-tilted reviewers)
│
└── if all the above run out of cycle (mid-2027) →
        NeurIPS / ICML workshop (Tier 4, 4-page version)
```

---

## Worst-case timeline

If everything rejects in sequence:

| Submit | Decision | Outcome |
|---|---|---|
| NeurIPS 2026 main | Sep 2026 | reject |
| ICLR 2027 | Jan 2027 | reject |
| TMLR (Tier 1, in parallel) | Mar 2027 | accept (likely, 8-12 wk cycle) |

**Worst plausible: TMLR accept by March 2027.** v15 NEVER LOSES the
ability to be published; the only question is venue prestige.

---

## What to keep stable across venues

The following do not change regardless of venue:

- Theorem 1 statement and proof
- 5 main experiments + their numerical claims
- v5.2 retraction disclosure
- Anonymous code repository (if anonymised)
- 52 references (may grow ~5 across cycles as new related work appears)

Only adapt: format / page limit / reproducibility-checklist style.

---

## Recommendation summary

1. **Submit NeurIPS 2026 main track on May 4-6, 2026.** This sprint
   has prepared everything. ~30% accept probability.
2. **If reject, default to ICLR 2027** (Oct 2026 deadline; ~8 hours
   to adapt) unless reviewer feedback specifically points to
   real-data weakness, in which case implement CR-2 from
   `camera_ready_prep_plan.md` first and submit TMLR.
3. **AISTATS 2027** is the conservative theory-tilted backup if
   ICLR also rejects (Oct 2026 deadline overlaps with ICLR cycle —
   pick one; don't double-submit).
4. **Workshops are last resort only.** Better to publish in TMLR
   than to settle for a workshop track.

The paper will be published. The only question is when and where.
