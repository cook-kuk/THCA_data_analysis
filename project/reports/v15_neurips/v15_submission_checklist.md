# v15 NeurIPS 2026 Submission Checklist

_Generated 2026-04-27. Single-source pre-submission checklist.
Companion to `v15_status_audit.md` (completeness measurement) and
`v15_missing_sections.md` (gap drill-down)._

## 0. NeurIPS 2026 deadline (verify before relying on this!)

The NeurIPS 2026 official call had not yet been published online at
the time of this audit; deadlines below are inferred from the
NeurIPS 2024 / 2025 schedules and should be **verified against the
official call** at https://neurips.cc/Conferences/2026 before
locking in any submission plan.

| Milestone                 | NeurIPS 2024 actual | NeurIPS 2025 actual | NeurIPS 2026 likely |
|---------------------------|--------------------:|--------------------:|--------------------:|
| Abstract registration     |          May 15 2024 |          May 11 2025 |    **mid-May 2026** |
| Full paper deadline       |          May 22 2024 |          May 15 2025 |  **mid-late May 2026** |
| Supplementary deadline    |       May 29 2024   |          May 22 2025 |        late May 2026 |
| Author response (rebuttal)|          Jul–Aug 2024 |        Jul–Aug 2025  |        Jul–Aug 2026 |
| Decisions                 |          Sep 2024   |          Sep 2025   |             Sep 2026 |
| Camera-ready              |        ~Oct 2024   |          ~Oct 2025  |             ~Oct 2026 |

**Action: confirm dates by checking** `https://neurips.cc/` **and**
`https://openreview.net/group?id=NeurIPS.cc/2026` **on first
opportunity.** This document assumes the canonical mid-May 2026
deadline.

## 1. Workdays-to-submission burndown

Assuming mid-May 2026 deadline and "today is 2026-04-27":

| Phase                                  | Workdays | Calendar |
|----------------------------------------|---------:|----------|
| Proof verification (10 TODOs)          |  0.5–1   | by Apr 30 |
| Bibliography expand 12 → 40+           |    0.25  | by May 1 |
| Broader-impact rewrite + ethics + code-availability statements | 0.25 | by May 1 |
| Reproducibility Checklist answers      |   0.25   | by May 2 |
| NeurIPS template port + appendix inline + figure PDF export | 0.5  | by May 4 |
| End-to-end re-read + reviewer simulation | 0.5     | by May 6 |
| **Buffer for unexpected reviewer-style holes** | 1–2 | May 6 – May 14 |
| **Submission**                         |    —     | by May 15 (target) |

Total focused writing: ~3 working days. Buffer to mid-May: ~7 working
days. Schedule is tight but realistic if Track 1 (paper-1 v5p2
methodology rewrite — see `v5p2_morning_followup_prompt.md`) is run
in parallel on different days.

## 2. Pre-submission gates (MUST pass)

### A. Mathematical content

- [ ] **Theorem 2 proof completed** — every TODO marker in
  `theorem2_proof.tex` either resolved or downgraded.
  Fallback if proof doesn't go through: rename Theorem → Proposition,
  restrict to Gaussian-only, defer general case to future work.
  See `v15_status_audit.md` §2.
- [ ] **Theorem 2 proof inlined into main PDF appendix** — currently
  the paper §A is a pointer; NeurIPS submission must be one PDF.
  Do this *after* the TODOs are resolved.
- [ ] **§5.2 ROC-AUC caption disambiguation** — paper says 0.78,
  task2 checkpoint records both 0.78 (broad def) and 0.62 (narrow
  def). Rewrite caption to make the definition explicit.

### B. NeurIPS submission requirements

- [ ] **NeurIPS 2026 .sty installed** — fetch
  `neurips_2026.sty` from `https://neurips.cc/Conferences/2026/CallForPapers`
  when published. Until then, develop in `neurips_2025.sty` (the 2026
  version is typically a backwards-compatible minor revision).
- [ ] **9-page main body limit** — current draft is well under
  (505 lines of LaTeX is roughly 8 typeset pages); NeurIPS 2025
  allowed 9 + unlimited appendix.
- [ ] **All figures as PDF/PNG inline** (not Plotly HTML iframes).
  Five figures need export: `theorem2_decomposition`,
  `stress_test_phase_diagram`, `tta_benchmark`, `scaling_plot`,
  `cross_domain_heatmap`. Use plotly's
  `fig.write_image(path, engine='kaleido')`.
- [ ] **Bibliography in NeurIPS bib style** — NeurIPS uses an
  unsorted-numerical bibliography; the current `\bibliographystyle{plain}`
  is acceptable for review but the camera-ready uses a custom style.
- [ ] **Anonymisation** — author block must be removed from the
  submitted PDF (use `[anonymous]` or NeurIPS's hide-author macro).
  Email address `kukshomr@gmail.com` currently in title block →
  remove for review version.
- [ ] **Anonymous code repository** — push v15 scripts to a fresh
  anonymous GitHub or use https://anonymous.4open.science for the
  review version.

### C. Required statements

- [ ] **Abstract ≤ 250 words** — currently 201 words ✅ (passes).
- [ ] **Broader Impact statement** — current paragraph is too brief
  for 2024+ NeurIPS expectations. Use template in
  `v15_missing_sections.md` §A.
- [ ] **Ethics statement** — 3-line minimal version in
  `v15_missing_sections.md` §E (synthetic data, no human subjects).
- [ ] **Code & data availability statement** — drop-in template in
  `v15_missing_sections.md` §C.
- [ ] **Reproducibility Checklist** — 14 questions, ~10 already
  answerable from existing checkpoints, 4 need minor additions
  (hyperparameter table, hardware spec, LICENSE file). See
  `v15_missing_sections.md` §F.

## 3. v17 self-citation triage — verdict: clean

**Question:** Should v15 cite v17, or v17 cite v15, or neither?

**Findings:**

- **v15 → v17 citation: NO.** The v15 paper does not mention v17,
  TROP2, sacituzumab, TACSTD2, DM1/DM2, or dark-matter clusters
  anywhere. This is correct: v15 is a theoretical / methodology
  paper about DIAL; v17 is a clinical / biological paper about a
  specific subtype axis. They share no scientific overlap.
- **v17 → v15 citation: YES (already present in v17 dashboard).**
  The v17 reports cite "DIAL" extensively: see
  `reports/v17/index.html`, `reports/v17/v17_summary.md`,
  `reports/v17p2/index.html`. The v17 dashboard navigation also
  surfaces a "v15 NeurIPS" link. v17 treats DIAL as established
  prior work and reports DIAL audit numbers ("DIA-AUC 0.969–1.000")
  on its own dark-matter sub-cohorts.

**Direction:** v17 → v15. No circular self-cite. Clean.

**v17 stale-number caveat:** the v17 reports were authored
*before* the v5.2 retraction (which only landed 2026-04-25). v17's
DIAL audit numbers come from the v5.1 leaky-pipeline run on the
dark-matter subset. The v17 paper-2 (TROP2 standalone) draft
already had its v5.1 leak paragraph rewritten on 2026-04-25/27;
the v17 *index dashboard* still shows leaky DIAL numbers and would
benefit from a v5.2 footnote when v17 itself ships. **This is not
a blocker for v15** — v15 is the source of truth for DIAL methodology
and stands on its own.

**Recommended citation language for v17 → v15** (drop-in for any
v17 paper that wants to cite the DIAL methodology):

> "DIAL (Direction-Invariant AUC Leakage), a post-hoc diagnostic for
> subspace-aligned conditional shift under linear batch correction
> [Cook 2026 NeurIPS]; we apply DIAL to evaluate dark-matter
> sub-cohort transfer to GSE27155 / GSE76039."

## 4. Pre-submission audit (run last, single-pass)

- [ ] Spell-check + grammar pass
- [ ] All `\TODO` macros removed or replaced with completed text
- [ ] All `\cite{X}` keys exist in `\thebibliography`
- [ ] All `\ref{X}` labels resolve
- [ ] All `\input{X}` files render without compilation error
- [ ] Camera-ready compile clean (no overfull hbox warnings if
  possible; allowed but ugly)
- [ ] Page count fits NeurIPS limit (9 main + unlimited appendix)
- [ ] Anonymisation scrubbed (author name, email, GitHub URLs,
  acknowledgements that reveal identity)
- [ ] Final PDF metadata stripped of author info (LaTeX leaks
  username via `\pdfauthor` defaults)

## 5. Post-submission to-do (after submit, before reviews)

- [ ] Rebuttal preparation: anticipate reviewer questions on
  - Theorem 2 Gaussian assumption tightness
  - Synthetic-only empirics; "where is the real-data?"
  - DIAL classifier-dependence; "is this really a 'metric' or a
    classifier-specific number?"
  - Foundation-model scaling negative result; "did you try
    pre-trained models X, Y, Z?"
- [ ] If accepted: camera-ready additions
  - De-anonymise author block
  - Add Acknowledgements (PI 유 교수님, v5.2 audit team)
  - Switch anonymous GitHub link to public
  - Optional: extend §5.5 cross-domain to one real-data benchmark
    (e.g., Camelyon17 or DomainNet) if reviewers ask

## 6. Suggested submission day-of plan

When the actual deadline date arrives (~May 15 2026):

1. Final compile + visual diff vs previous accepted version
2. Submit via OpenReview ≥ 4 hours before hard deadline
3. Verify PDF rendering on OpenReview (some `\hyperref` settings
   break silently)
4. Submit anonymous code link as supplementary material
5. Confirm submission ID and email confirmation
6. **Don't refresh the OpenReview tab obsessively.** Reviewers
   come ~3 weeks later.

## 7. Risk register (honest)

| Risk | Likelihood | Mitigation |
|------|-----------|------------|
| Theorem 2 proof step (c) doesn't go through | Medium | Downgrade to Gaussian-only Proposition (1 h) |
| Reviewer flags v5.2 retraction as "honesty problem" | Low | Disclosure is already explicit; methodology contribution stands |
| Reviewer wants real-data empirics | Medium | Promise camera-ready addition; v5.2 leak example IS real data, framed correctly |
| NeurIPS 2026 .sty has incompatible breaking changes | Low | Develop in 2025 .sty as fallback; port at last minute |
| Out-of-deadline because of v17p35 / paper-2 sprint conflict | Medium-High | Schedule explicitly: v15 paper writing on Mon/Wed/Fri, v17/p2 work on Tue/Thu |

---

**Decision gate (re-read before each work session):**
the v15 paper is **submission-realistic for NeurIPS 2026 mid-May**
provided one focused proof-verification day, one focused writing day,
and the v17p35 dark-matter sprint does not consume the entire
remaining calendar. If forced to choose between the two, paper-1
v5p2 rewrite has higher publication probability than v15 NeurIPS
(NeurIPS acceptance rate ~25 %, Bioinformatics ~50 %); but v15
NeurIPS has higher *visibility* impact if accepted. Both are worth
doing. The realistic plan is: **finish v15 in the first half of
May, paper-1 / paper-2 in the second half**.
