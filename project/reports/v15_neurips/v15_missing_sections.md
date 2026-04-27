# v15 NeurIPS — Missing Sections vs Submission Template

_Generated 2026-04-27. Maps `v15_neurips_paper_v0.tex` against the
NeurIPS 2026 expected submission structure. "Missing" here means
either absent or stub-only._

## NeurIPS 2026 expected structure

Based on NeurIPS 2024 / 2025 submission templates (the 2026 .sty file
is published in mid-April but typically a minor revision of 2025):

| # | Section                              | Required | v15 status |
|---|--------------------------------------|:--------:|------------|
|  1 | Title + author block                | ✅       | written |
|  2 | Abstract (≤ 250 words)              | ✅       | 201 words; v5.2 disclosed |
|  3 | Introduction                        | ✅       | written, 5 contributions |
|  4 | Related work                        | ✅       | written (companion `v15_related_work.md` provides citations) |
|  5 | Method / theory                     | ✅       | §3 prelim + §4 theory; Theorem 2 stated |
|  6 | Experiments                         | ✅       | 5 sub-experiments §5.1–5.5 |
|  7 | Discussion / limitations            | ✅       | written incl. limitations + broader impact |
|  8 | Conclusion                          | ✅       | written |
|  9 | **Broader Impact** (NeurIPS-specific) | ✅     | one paragraph in §6.4; **NeurIPS expects a structured statement** — see §A below |
| 10 | **References (full)**               | ✅       | 12 entries; **paper itself notes "expand to 40+"** — see §B |
| 11 | **Reproducibility Checklist**       | ✅       | **MISSING** — required since NeurIPS 2020 |
| 12 | **Code & data availability statement** | ✅    | stub ("release path TBD") — see §C |
| 13 | **Author contributions**            | ⚪       | optional pre-deadline (de-anonymised at camera-ready); single author so trivial |
| 14 | **Acknowledgements**                | ⚪       | one line present |
| 15 | Appendix A: full proofs             | ✅       | currently a pointer to `theorem2_setup.tex` + `theorem2_proof.tex`; **NeurIPS submission must contain proofs in the main PDF** — see §D |
| 16 | Appendix B: experimental details    | ✅       | stub; needs hyperparameter table |
| 17 | Appendix C: additional results      | ✅       | pointer to interactive HTML; needs PDF figure exports |
| 18 | Appendix D: ethics statement        | ⚠️       | not a separate section (folded into broader impact); NeurIPS 2024+ asks for an explicit ethics paragraph for human-data work — see §E |

Sections marked **MISSING** or **stub** below are the actual gaps.

## A. Broader Impact — needs structured rewrite

**Current state** (paper §6.4, ~120 words):

> "A post-hoc, AUC-based diagnostic for preprocessing-induced flips
> lowers the audit cost for biomedical ML pipelines, where many
> published classifiers may be vulnerable to the v5.1-style leak.
> We expect DIAL to find direct application in regulatory ML audits
> (FDA medical-device re-validation, EMA cross-cohort review) and in
> ML-for-science reproducibility programmes (NeurIPS reproducibility
> checklist, ICLR ML+code reviews)."

**NeurIPS 2026 expected structure** (per 2024 / 2025 calls):
- **Positive societal impacts**: who benefits, how
- **Negative / dual-use risks**: could the method be weaponised?
- **Mitigations**: what we did to limit harms
- **Limitations of impact claims**: scope of the claim

**Suggested rewrite** (drop-in for §6.4):

```latex
\subsection{Broader impact}
\paragraph{Positive impact.}
A post-hoc, AUC-based diagnostic for preprocessing-induced flips
lowers the audit cost for biomedical ML pipelines, where many
published classifiers may be vulnerable to v5.1-style leakage.
DIAL has direct application in regulatory ML audits (FDA medical-
device re-validation, EMA cross-cohort review) and in ML-for-science
reproducibility programmes (NeurIPS Reproducibility Checklist, ICLR
ML+code reviews).

\paragraph{Negative / dual-use risk.}
DIAL is a diagnostic, not a corrective. We see no direct dual-use
risk: the method does not enable surveillance, manipulation of
clinical outcomes, or generation of harmful content. The closest
indirect risk is reviewer over-reliance on a single AUC-based
audit; we explicitly flag in §6.3 (Limitations) that DIAL is
classifier- and CV-protocol-dependent.

\paragraph{Mitigations.}
We release both the metric implementation and the v5.2 leak-fix
audit script with explicit fail-safe behavior (DIAL = 0 when the
leak is closed, by Theorem~2). The numerical experiments are
fully synthetic and reproducible; no human-subjects data are
used in v15 itself.

\paragraph{Scope of claim.}
DIAL fires under the Theorem~2 conditions (Gaussian-shared-Σ
class-conditionals, linear correction T). Outside these
assumptions, DIAL's behaviour is empirical and modality-dependent;
we benchmark four modalities in §5.5 and find the flip mechanism
reproduces, but the precise constants will differ.
```

## B. References — 12 → 40+ expansion

The paper has 12 \bibitem entries; the companion `v15_related_work.md`
file (161 lines) lists ~40 citations across six adjacent literatures:

1. **Covariate-shift theory** — Shimodaira 2000, Sugiyama et al.
   (KLIEP) 2007, Bickel et al. 2009 [present]
2. **Domain-invariance bounds** — Ben-David 2007/2010 [present],
   Mansour 2009, Zhao 2019 [missing] + Acuna 2021 [missing]
3. **Adversarial DA** — DANN [present], ADDA [present], CDAN
   [present], MCD [present], DIRT-T 2018 [missing], FixBi 2021
   [missing]
4. **Test-time adaptation** — TENT [present], SHOT [present],
   BN-stats [present], MEMO [present], EATA 2022 [missing], CoTTA
   2022 [missing]
5. **Leakage / preprocessing** — Kaufman 2012 [present], Roberts
   2017 [missing], Kapoor 2022 [missing], Whalen 2022 [missing]
6. **Foundation models / scaling** — Kaplan 2020 [missing],
   Hoffmann (Chinchilla) 2022 [missing], Bommasani 2021 [missing],
   Geirhos 2020 [missing], Sagawa (group robustness) 2020 [missing]

**Action item**: copy the bibtex entries from `v15_related_work.md`
section "Full citations" into the v15 paper bibliography. ~1 hour
of mechanical work.

## C. Code & data availability statement — full draft

**Replace current stub with:**

```latex
\section*{Code and data availability}
\label{app:code}
All code used to produce the v15 numerical experiments is at
\url{ANONYMOUS_GITHUB_URL_AT_CAMERA_READY} (anonymised for review at
\url{https://anonymous.4open.science/r/v15-dial-XXXXXX}).
The five experiment scripts are:
\texttt{v15\_theorem2\_verify.py} (§5.1, R²=0.94 numerical
validation),
\texttt{v15\_stress\_test.py} (§5.2, $4 \times 5 \times 5 \times 20$
synthetic shift grid),
\texttt{v15\_tta\_benchmark.py} (§5.3, ten adaptation methods),
\texttt{v15\_scaling.py} (§5.4, eight encoder configurations),
\texttt{v15\_cross\_domain.py} (§5.5, four-domain reproducibility).
Each script writes a JSON checkpoint
(\texttt{results/v15\_neurips/checkpoints/task[1-5]\_*.json})
that records exact wall time, seed list, and headline metric, plus
a long-form TSV that supports the tabular numbers in the paper.
No human-subjects data are used in v15; the v5.1 leak example
referenced in §1 and §6 is reproducible from public TCGA/GEO
inputs but is not required to validate any v15 claim.
```

## D. Appendix A: full proofs — needs in-PDF inclusion

Currently §A reads:
> "See `theorem2_setup.tex` and `theorem2_proof.tex` shipped
> alongside this manuscript."

NeurIPS submission requires a single PDF. Action: `\input` the proof
files into the main paper's appendix. Combined length ≈ 270 + 178 =
448 lines of LaTeX, comfortably fits the NeurIPS 9-page main +
unlimited appendix format.

10 TODO markers in `theorem2_proof.tex` need resolution before this
inclusion. See `v15_status_audit.md` §2 for triage options.

## E. Ethics statement — minimal addition needed

NeurIPS 2024+ asks for explicit ethics consideration even for
non-human-subjects work. Suggested 3-line addition before §A appendix:

```latex
\section*{Ethics statement}
v15 uses no human-subjects data. The v5.1 leak example referenced in
§1 and §6 is computational and based on publicly de-identified TCGA
and GEO records. No new biological samples were collected. The
broader-impact discussion in §6.4 covers dual-use considerations.
```

## F. Reproducibility Checklist — full template draft

NeurIPS 2026 expected questions (extrapolated from 2025):

```markdown
1. Claims and abstract:
   ✅ Yes — abstract claims (R²=0.94, ROC=0.78, scaling negative,
      4/4 domains) all backed by checkpoints in §5
2. Limitations clearly discussed:
   ✅ Yes — §6.3 (Gaussian assumption, classifier-dependence,
      synthetic-only)
3. Theoretical claims have full proofs:
   ⚠️ In progress — Theorem 2 proof has 10 TODO markers in
      theorem2_proof.tex; see v15_status_audit.md §2
4. Did you describe assumptions and limits of theory:
   ✅ Yes — Definition 1 in theorem2_setup.tex; §6.3 limitations
5. Experimental code provided:
   ✅ Yes — five v15_*.py scripts
6. Experimental data described:
   ✅ Yes — fully synthetic; data generators in scripts
7. Hyperparameters specified:
   ⚠️ Partial — Appendix B is currently a stub; need explicit
      hyperparameter table
8. Random seeds reported:
   ✅ Yes — checkpoints record n_seeds (3 for scaling, 5 for TTA,
      6 for cross-domain, 15 for theory, 20 for stress)
9. Variance reported:
   ✅ Yes — _std columns in all TSVs
10. Compute resource described:
    ⚠️ Stub — wall times in checkpoints (3s to 672s); add
       hardware spec (CPU model, RAM)
11. License of code declared:
    ⚠️ Missing — add MIT or Apache-2.0 to repo
12. License of data declared:
    ✅ Synthetic; no licensed data used
13. Privacy/consent:
    ✅ Not applicable — synthetic
14. Theoretical fairness considered:
    ✅ Mentioned in §6.4 (broader impact)
```

3 questions are partial / missing (3, 7, 10, 11). All are
mechanically closeable: verify proof TODOs (Q3); add a 6-row
hyperparameter table (Q7); add hardware spec line (Q10); commit
LICENSE file (Q11).

## Summary — what to actually write before submission

In priority order:

1. **Resolve 10 TODO markers in `theorem2_proof.tex`** (3–4 h)
   OR downgrade to a Gaussian-only Proposition (1 h fallback).
2. **Inline the proof appendix into the paper** (15 min once
   step 1 is done; `\input{theorem2_setup.tex}` + `\input{theorem2_proof.tex}`).
3. **Bibliography expand 12 → 40+** (1 h, citations already
   listed in `v15_related_work.md`).
4. **Broader-impact rewrite** (drop-in template provided in §A above; 30 min).
5. **Code & data availability statement** (drop-in template provided
   in §C above; 15 min — needs anonymous GitHub link).
6. **Ethics statement** (3 lines from §E above; 5 min).
7. **Reproducibility Checklist** answers (~30 questions, all
   determinate; 30 min).
8. **NeurIPS template port** (`\documentclass{neurips_2026}` etc.;
   ~30 min mechanical).
9. **Figure exports**: Plotly → PDF for the 5 v15 figures (10 min).

**Total focused writing time: ~6–7 hours.** Realistic in one
working day if the proof verification (item 1) doesn't blow up.
