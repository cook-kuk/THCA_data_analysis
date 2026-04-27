# v15 — 9-Day Submission Countdown

_From 2026-04-27 (today) to 2026-05-06 (full paper deadline AOE)._

> Critical dates (verified 2026-04-27 at neurips.cc):
> - **Phase 1 abstract registration: May 4, 2026 AOE** (7 days)
> - **Phase 2 full paper + supplementary: May 6, 2026 AOE** (9 days)
> - All authors must have an OpenReview profile by Phase 1 deadline.

## Day-by-day

### **Day 0 — 2026-04-27 (Mon)** — TODAY

- ✅ Sprint complete (24 tasks).
- ✅ All 5 scripts re-verified end-to-end.
- ✅ PDF compiles clean, 15 pages, anonymisation passes.
- 🟡 *You*: confirm your OpenReview profile is active. Visit
  https://openreview.net/signup or the existing profile page. **2 min**.

### Day 1 — 2026-04-28 (Tue) — buffer / personal review day

- 🟡 *You*: read the PDF cover-to-cover once. **30 min**.
  Path: `submission/v15_neurips/v15_neurips_SUBMIT_FINAL.pdf`
- 🟡 *You*: skim `rebuttal_prep_v1.md` to confirm the 5 anticipated
  probes feel right. **15 min**.
- ⚠️ Do NOT touch the LaTeX yet — first read should be reviewer-eyes.

### Day 2 — 2026-04-29 (Wed) — anonymous repo

- 🟡 *You*: upload `anonymous_code.zip` to https://anonymous.4open.science.
  Get the URL of the form `https://anonymous.4open.science/r/v15-dial-XXXXXX`.
  **10 min**.
- 🟡 *You*: replace placeholder URL in
  `reports/v15_neurips/submit/v15_missing_drafts_v2.tex` line ~58.
- 🟡 *You*: recompile:
  ```bash
  cd reports/v15_neurips/submit
  pdflatex v15_neurips_SUBMIT.tex
  bibtex v15_neurips_SUBMIT
  pdflatex v15_neurips_SUBMIT.tex
  pdflatex v15_neurips_SUBMIT.tex
  cp v15_neurips_SUBMIT.pdf ../../../submission/v15_neurips/v15_neurips_SUBMIT_FINAL.pdf
  ```
- 🟡 *You*: re-run anonymity grep:
  ```bash
  python3 -c "
  from pdfminer.high_level import extract_text
  txt = extract_text('submission/v15_neurips/v15_neurips_SUBMIT_FINAL.pdf')
  for m in ['Seungho','kukshomr','Cornerstone','SNUBH']:
      assert m not in txt, f'leak: {m}'
  print('anonymisation: clean')
  "
  ```
  **15 min total**.

### Day 3 — 2026-04-30 (Thu) — final read + spelling

- 🟡 *You*: full read-through with grammar/spelling tools. Use
  Grammarly or LaTeX-aware spell-check. **45 min**.
- 🟡 *You*: verify all 5 figures display correctly in the PDF
  (titles, axis labels, legend). **10 min**.

### Day 4 — 2026-05-01 (Fri) — buffer / fix any issues

- ⚪ Reserved for any issues surfaced during Day 3 read.
- ⚪ If no issues: rest day; do NOT keep tweaking.

### Day 5 — 2026-05-02 (Sat) — weekend / OPTIONAL final pass

- ⚪ Weekend rest. Do NOT submit yet — Phase 1 deadline is May 4.
- 🟢 If time permits: read the rebuttal_prep_v1 deeper, mentally
  preparing for the 5 probes that will land in August.

### Day 6 — 2026-05-03 (Sun) — pre-submit dry-run

- 🟡 *You*: do a "dry run" of the OpenReview submission flow:
  log in, navigate to NeurIPS 2026 submission portal, verify the
  form fields exist, check upload limits. **15 min**. Don't submit.
- 🟡 *You*: verify deadline timezone. AOE = Anywhere on Earth =
  UTC-12, so "May 4 AOE" effectively gives until 2026-05-05 11:59
  UTC for the abstract registration.

### **Day 7 — 2026-05-04 (Mon) — PHASE 1 ABSTRACT DEADLINE** ⚠️

- 🟢 *You*: register the abstract on OpenReview. Use text from
  `submission/v15_neurips/openreview_form_fields.md`:
  - Title (single line)
  - Abstract (250 words; copy-paste)
  - Keywords (10 keywords)
  - Authors / co-authors registered with OpenReview profiles
  - Primary subject area: "Domain adaptation / transfer learning /
    out-of-distribution generalization"
  - Conflict of interest: standard
- 🟢 **Submit before 11:59 UTC on May 5** (= AOE end of May 4).
  Do NOT push to deadline-day.
- 🟢 Take screenshot of submission ID + confirmation email.
- ⚠️ **Once abstract is submitted, no new abstracts can be added —
  Phase 2 must use this exact title.**

### Day 8 — 2026-05-05 (Tue) — final compile + submission prep

- 🟡 *You*: final pdflatex pass on the abstract-locked title.
- 🟡 *You*: re-verify the anonymous repo URL still points to
  the right snapshot. **15 min**.

### **Day 9 — 2026-05-06 (Wed) — PHASE 2 PAPER DEADLINE** 🎯

- 🟢 *You*: upload Phase 2 to OpenReview against the existing
  Phase 1 abstract record:
  - `v15_neurips_SUBMIT_FINAL.pdf` (main)
  - `anonymous_code.zip` (supplementary)
- 🟢 **Submit before 11:59 UTC on May 7** (= AOE end of May 6).
- 🟢 Take screenshot of final submission ID + confirmation email.
- 🟢 Verify PDF rendering on OpenReview after upload.

---

## Time budget summary

| Day | Activity | Time |
|---|---|---:|
| 0 (today) | OpenReview profile | 2 min |
| 1 | First read-through | 45 min |
| 2 | Anon repo upload + recompile | 25 min |
| 3 | Spell/grammar pass + figure check | 55 min |
| 4 | Buffer | 0–60 min |
| 5 | Weekend | 0 min |
| 6 | Submission dry-run | 15 min |
| 7 | **Phase 1 submit** | 15 min |
| 8 | Final prep | 15 min |
| 9 | **Phase 2 submit** | 15 min |
| **Total** | | **~3 hours** |

The 60-minute total in `v15_NEURIPS_SUBMIT_DUMP.md` was the
mechanical time. Adding the personal-review and timezone-buffer
work brings it to ~3 hours over 9 days, which is realistic
alongside other paper / project work.

## Risk: do NOT do these things

- ⛔ Do NOT push to deadline-day. AOE is forgiving but OpenReview
  servers are slow at peak times.
- ⛔ Do NOT modify the abstract or title between Phase 1 and Phase 2.
- ⛔ Do NOT submit a different version of the PDF to a different
  field; OpenReview submissions must be self-consistent.
- ⛔ Do NOT skip the OpenReview profile pre-check; new accounts
  cannot submit the same day they're created.
- ⛔ Do NOT include de-anonymising hints (acknowledgements with
  named PIs, GitHub URLs to your personal repo, email addresses).

## Backup plan

If OpenReview is down on submission day (rare but happens), email
papers@neurips.cc with timestamped PDF as evidence of intent.
This is the standard NeurIPS recovery path.
