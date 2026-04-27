# Current State & Open Decisions — For Outside Opinion

_Generated 2026-04-27. Self-contained problem dump for web Claude
(claude.ai). PhD-level thyroid-cancer ML project, 2 papers in parallel
preparation, deadline pressure._

---

## TL;DR (read first)

I have **2 papers within ~9 days of deadline**, and I just caught
two material problems that I want a second opinion on before I
submit:

1. **v17 (npj Precision Oncology, ready to submit this week)** —
   the TERT prognostic claim collapses under stage adjustment.
   Univariate HR = 6.31 (p = 3×10⁻⁴) → multivariate HR = 1.88
   (p = 0.29). I've reframed it as "molecular handle for Stage
   III/IV" but I'm not 100 % sure that's strong enough.
2. **v15 (NeurIPS 2026, abstract May 4 / paper May 6)** —
   methodology paper on a metric called DIAL. Self-review predicts
   mean reviewer score 6.0/10 (borderline accept, ~30 % accept
   probability). Should I submit main track or workshop?

I want a fresh outside read because I've been in the loop for
6 months and might be missing something.

---

## Project context (1 paragraph)

PhD-level thyroid cancer transcriptomic subtyping → drug
actionability → paper. Two papers landing simultaneously:

- **v17 npj** = clinical / biology paper. 8-gene RAI-responsiveness
  panel that outperforms BRAF V600E; DM1/DM2 transcriptomic axis;
  hot/cold immune separation; TERT 4-group survival stratification
  (the reframed bit).
- **v15 NeurIPS** = methodology paper. DIAL = "Direction-Invariant
  AUC Leakage" diagnostic for cross-cohort batch correction. Theorem
  + 5 synthetic experiments + 1 honestly-disclosed retraction
  worked-example.

These are different papers on different cohorts, but they share an
author and a 6-month sprint cycle.

---

## What just happened (yesterday + today)

**v8.1 → v5.2 retraction (resolved 2 days ago).**
Earlier in the project, v5.1 reported THCA "4 of 5 classifiers flip
post-ComBat" with DIAL = 0.494. A self-audit found this was a
leakage artefact: ComBat was fit on the full pooled matrix
*before* the LODO split. Under proper LODO ComBat (fit-on-train-
only), all 5 THCA classifiers move to `true_biology` (DIAL = 0).
v5.1 paper retracted; v15 NeurIPS paper now uses the retraction as
a *worked example* of what DIAL detects (this is the cleanest
honest-disclosure framing we could find).

**v17 honest audit (just landed today).**
TERT promoter mutations were recovered from cBioPortal and a
4-group survival stratification (BRAF only / RAS only / TERT⁺ /
triple-neg) was reported with logrank p = 4.92×10⁻⁶. A pre-
submission honest audit identified three robustness gaps:

- Stage and age confounding (TERT⁺ patients are median age 64.6
  vs 46 elsewhere; 61 % Stage III/IV vs 21–32 %)
- Small event count (16 total events / 6 in TERT⁺)
- Single-event sensitivity not tested

I ran the full audit:

| Test | Result |
|---|---|
| TERT univariate Cox HR | **6.31** (95 % CI 2.34–17.04, p = 3×10⁻⁴) |
| TERT + age Cox | 2.54 (CI 0.87–7.40, p = 0.087) |
| TERT + stage Cox | 2.92 (CI 0.87–9.80, p = 0.082) |
| **TERT + stage + age + sex (full multivariate)** | **1.88 (CI 0.58–6.09, p = 0.29)** |
| 1000-iter bootstrap median p | 2.6×10⁻⁵ (95 % CI 3.2×10⁻¹⁵ – 0.21) |
| 36-patient leave-one-out worst-case p | 7.5×10⁻⁴ (all 36 LOO p < 0.001) |

**Reading:** the univariate signal is real and LOO-stable
(no single event drives it). But after multivariate adjustment for
stage + age + sex, TERT does not retain stage-and-age-independent
prognostic value. Age is the dominant predictor (concordance 0.954
in the full model, but TERT itself is no longer the discriminating
variable).

I rewrote the manuscript R8 paragraph honestly:
- Title was "TERT defines an independent aggressive axis"
- Now: "TERT-promoter status as a molecular handle for advanced
  age and stage"
- Both univariate and multivariate Cox HRs are explicitly reported
- Bootstrap and LOO sensitivity are added as robustness panels
- Limitations item 11 explicitly notes the stage-confounding
  caveat and frames Korean-cohort cross-validation as the next
  test

---

## Problems by severity

### 🔴 Critical — decision-blocking

#### Problem 1: Is the v17 Scenario-B reframe publishable as-is?

The 8-gene RAI panel (the *primary* contribution) is unaffected by
this audit — AUC 0.954 vs 0.822 for BRAF-V600E-only, ΔAUC = +0.132,
external-validation AUC 0.974 on GSE76039 PDTC+ATC. That part is
solid.

The TERT 4-group (the secondary contribution that I was excited
about) is now framed honestly as a Stage-III/IV-enrichment marker
rather than a prognostic axis. The univariate stats are real but
the multivariate stats say "stage and age dominate".

**My read**: this is publishable; reviewers respect honest
disclosure and the bootstrap+LOO robustness checks address the
small-n concern. But I am not 100 % sure. **Three options:**

- **A. Submit v17 this week** with Scenario-B reframe (the work
  I just did); accept that npj acceptance probability drops from
  ~75 % to ~65–70 %.
- **B. Demote TERT 4-group to Supplementary**; lead with
  8-gene RAI panel only as a single-contribution paper.
- **C. Wait 4–8 weeks** for Korean cohort data (active outreach
  pending) which would let us re-test multivariate HR on an
  independent cohort. ~62 % BRAF-mutation rate in Korean PTC vs
  ~45 % in TCGA — substantively different population.

I think A is right. But the "wait for Korean cohort" option is
genuinely tempting because the multivariate finding might be a
TCGA-specific artefact (TCGA is exceptional-prognosis, only 16
events out of 504), and a Korean cohort might show TERT retains
stage-independent significance. Question for outside opinion:
does the bird-in-the-hand A beat the bird-in-the-bush C?

#### Problem 2: v15 NeurIPS borderline (P ~30 %); main track or workshop?

NeurIPS 2026 deadlines: **abstract May 4 AOE, paper May 6 AOE**
(confirmed at neurips.cc/Conferences/2026; I had earlier
extrapolated mid-May which was 7–10 days too generous).

I ran a 3-reviewer self-simulation with NeurIPS rubric (soundness,
presentation, contribution, overall):

| Reviewer | Profile | Soundness | Presentation | Contribution | Overall |
|---|---|---:|---:|---:|---:|
| R1 | theorist | 3 | 3 | 3 | 6 |
| R2 | empiricist | 2 | 3 | 3 | 5 |
| R3 | practitioner | 3 | 4 | 4 | 7 |
| **mean** | | **2.7** | **3.3** | **3.3** | **6.0** |

NeurIPS 2025 main-track acceptance rate ~25 %. Predicted accept
P ~30 %. The biggest risk is R2 (empiricist) caring about the
synthetic-only empirics. The paper has 5 synthetic experiments
+ 1 worked example of a real-data leak retraction; no real-data
benchmark like Camelyon17 / DomainNet.

**Three options:**

- **A. Submit main track May 4-6.** ~30 % accept; if rejected,
  fallback to ICLR 2027 (Oct 2026) or TMLR (rolling, ~50 % accept,
  no deadline pressure).
- **B. Submit workshop track instead** (NeurIPS has Workshop on
  Distribution Shifts, OOD, etc.). ~50 % accept; lower prestige;
  4-page extended abstract format.
- **C. Skip NeurIPS, go straight to TMLR** (rolling, ~50 % accept,
  reasonable prestige for ML, longer cycle).

I lean A (the upside justifies the 70 % rejection risk because
TMLR is a viable fallback and doesn't expire). But the
"submit-anyway" reasoning could be motivated reasoning. Question:
am I anchoring on prestige?

### 🟡 Material — worth discussing

#### Problem 3: v17 cover letter still claims old TERT framing

The current `cover_letter_v2.md` was written 2 days ago with the
"independent prognostic axis" framing. With the v3 reframe, the
cover letter needs corresponding updates. ~30 min of work but not
yet done.

#### Problem 4: v17p35 dashboard pages may carry pre-v5.2 numbers

A separate v17 dashboard hub (`reports/v17/`, `v17p2/`, `v17p3/`,
`v17p35/`) was authored before the v5.2 retraction landed. Those
dashboard pages cite "DIA-AUC 0.969–1.000" on dark-matter sub-
cohort transfer to GSE27155 / GSE76039 — those numbers were
computed under the v5.1 leaky-ComBat protocol that was retracted.

I added v5.2 retraction banners to 16 documents but the v17
dashboards may still show leaky numbers. If a reviewer Googles the
v17 paper they might find the dashboard and notice the
inconsistency.

#### Problem 5: Korean cohort timing (active outreach but uncertain)

The v17 paper Limitations now explicitly cite Korean / Asian
cohort cross-validation as the next step. I have an active
outreach to 분당서울대 (Bundang Seoul National University Hospital,
where my PI 유 교수님 has connections). IRB approval is pending;
realistic data delivery is **2–6 months**. If accepted, a Korean
cohort with substantially higher BRAF rate (~62 % vs ~45 % in
TCGA) is exactly the right test.

But I have no leverage to accelerate IRB. So this is something that
will land *after* v17 npj review, not before.

### 🟢 Operational — just need execution

- v15 anonymous code repo upload (10 min user task)
- v15 PDF recompile after URL replacement (5 min)
- v17 cover letter v3 update (30 min)
- 유 교수님 (PI) author confirm (5 min)
- Both papers under review simultaneously = August rebuttal load
  (~10 hours each, paged across 1–2 weeks each)

---

## What I want web Claude to weigh in on

I am NOT asking for help on the technical work (statistics, code,
LaTeX). All of that is already done and verified.

I AM asking for an **outside read** on three strategic questions:

### Question 1 (most important)

> Given that v17's TERT 4-group claim collapsed under multivariate
> adjustment but the 8-gene RAI panel is unaffected, **is the
> Scenario-B reframe ("TERT⁺ as a molecular handle for Stage III/IV
> identification") actually publishable in npj Precision Oncology?**
> Or am I rationalising a weakened story?
>
> The honest version: univariate HR 6.31 (p=3e-4), multivariate
> HR 1.88 (p=0.29), bootstrap median p 2.6e-5, LOO worst-case
> p 7.5e-4.
>
> Three valid options: (A) submit reframed, (B) demote TERT to
> supplementary and lead with 8-gene panel only, (C) wait for
> Korean cohort.

### Question 2

> NeurIPS 2026 main track acceptance ~25 %, my predicted accept ~30 %.
> The paper has a Theorem with explicit Gaussian-shared-Σ assumption,
> 5 synthetic experiments, and 1 honestly-disclosed retraction
> worked-example (the v5.2 leak). NO real-data benchmark.
>
> **Should I submit NeurIPS main track May 4-6 (P ~30%, fallback
> is TMLR rolling), or skip directly to TMLR (P ~50%, no deadline
> pressure, lower prestige)?** Workshop is the third option.
>
> Rebuttal-ready answers exist for the 5 anticipated reviewer
> probes. The honest disclosure of the v5.2 retraction reads as a
> strength to me, but a reviewer could read it as "this paper's
> motivation was wrong; weak premise". Am I being honest about
> that risk?

### Question 3

> Two papers under review simultaneously means August rebuttal
> phase coincides with v17 npj's first revision and v15 NeurIPS's
> author response. **Is this a reasonable load for a single PhD
> author, or am I setting myself up for a poor rebuttal on at
> least one paper?**
>
> Constraints: I have a PI (유 교수님) who can review but not
> author-respond on my behalf. I have an active 분당서울대 outreach
> that may need follow-up in August. I have a v17p35 dashboard
> that needs cleaning before public release. I have an Anthropic
> application I'd like to send by end of summer.

---

## Self-honesty notes

Things I'm aware I might be wrong about:

1. **I might be over-trusting LOO sensitivity as evidence of
   robustness.** All 36 LOO iterations p < 0.001 sounds strong,
   but LOO is sensitive only to *single*-patient removal, not to
   the choice of cohort or to the underlying multivariate
   confounding. The multivariate Cox is the more important test,
   and that one says HR = 1.88 (p = 0.29).

2. **Self-review reviewer scoring is probably biased upward.**
   I can predict R1 (theorist) and R3 (practitioner) accurately,
   but R2 (empiricist) is the one I tend to under-weight, and
   exactly the one most likely to reject for synthetic-only.

3. **The retraction story might read worse than I think.** I keep
   framing it as "honest disclosure → reviewer respect", but a
   tired reviewer might just read "the original motivation was
   wrong" and downvote.

4. **PhD-level workload assessment is hard from inside.** If a
   senior researcher told me "submit one, hold the other for 3
   months", I would probably listen. I'm asking that question
   genuinely.

---

## Numbers / files (for verifiability if needed)

- v17 paper: `submission/npj/manuscript_v3.md` (4 452 words)
- v17 robustness analysis: `results/v17_actual/A1{A,B,C}_*.{tsv,json}`
- v17 robustness dashboard: `reports/html/pages/v17_npj_robustness.html`
- v15 paper: `submission/v15_neurips/v15_neurips_SUBMIT_FINAL.pdf`
  (15 pages, 356 KB)
- v15 self-review: `submission/v15_neurips/self_review_simulation.md`
- v15 deadline: confirmed May 4 / May 6 AOE
- v5.2 retraction summary: `reports/v5p2/v5p2_critical_assessment.md`

---

## What I am NOT asking

- Help with the technical statistics — verified
- Help with the LaTeX / formatting — done
- Help with the proofs — closed in v15 audit
- Help with code reproducibility — verified end-to-end run
- Help with figure rendering — 5/5 figures generated cleanly
- Strategic re-design of the papers — too late, deadline is in 9 days

I am ONLY asking for an outside opinion on the three strategic
questions above (mainly Q1: should v17 ship as-is with the
TERT reframe, or wait, or restructure).

---

## Single-line summary for web Claude

I have two papers within 9 days of deadline. The bigger one
(v17 npj) just got its secondary contribution honestly downgraded
to a "molecular handle for Stage III/IV" framing after a
multivariate Cox audit. I want an outside read on whether to ship
v17 this week with the reframe, wait for Korean cohort data, or
restructure. And whether v15 NeurIPS at predicted P ~30 % is worth
the deadline-rush risk vs going to TMLR (rolling, P ~50 %,
no deadline).
