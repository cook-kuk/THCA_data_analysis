# v15 NeurIPS — Final Submission Audit

_Generated 2026-04-27 by S3-A._

## Compile status

| Check                              | Result |
|------------------------------------|--------|
| `pdflatex` clean (no errors)       | ✅ pass |
| `bibtex` resolved all citations    | ✅ pass |
| `\TODO` macros remaining           | ✅ 0 |
| `\ref{}` undefined                 | ✅ 0 |
| `\cite{}` undefined                | ✅ 0 |
| Page count                         | 12 pages (main ~8 + appendix ~4); within NeurIPS 2025 10-page main limit |
| PDF size                           | 345 KB |

## Anonymisation

| Check                                  | Result |
|----------------------------------------|--------|
| Author block in title                  | ✅ "Anonymous / Anonymous Affiliation / anonymous@anonymous.org" |
| Email `kukshomr@gmail.com` in PDF      | ✅ not found via `strings` grep |
| "Seungho Cook" in PDF                  | ✅ not found |
| "Cornerstone Partners" in PDF          | ✅ not found |
| Anonymous code repo path placeholder   | ✅ `https://anonymous.4open.science/r/v15-dial-XXXXXX` |
| Acknowledgements removed for review    | ✅ none in the SUBMIT.tex |

## Content gates

| Item                                  | Result |
|---------------------------------------|--------|
| Theorem 2 status                      | ✅ **STAYS as Theorem** (8/8 TODOs resolved per `theorem2_proof_audit.md`) |
| Bibliography entries                  | ✅ 52 (target ≥ 40) |
| Missing sections (broader impact, ethics, code, repro) | ✅ all 4 inlined via `v15_missing_drafts_v2.tex` |
| 5 figures inlined as PDF              | ✅ all 5 (`figures/figure[1-5].pdf`) |
| Hyperparameters table (Repro Q7)      | ✅ Appendix `app:hparams` |
| Reproducibility checklist (14 Q)      | ✅ all 14 answered |
| Ethics statement                      | ✅ `\section*{Ethics statement}` |
| Code & data availability statement    | ✅ `\section*{Code and data availability}` |
| LICENSE in code repo                  | ✅ MIT, 1 KB |
| Anonymous code zip                    | ✅ 69 KB; 5 scripts + 5 checkpoints + 9 TSVs |

## Numerical consistency

| Claim                              | Source                                | Verified |
|------------------------------------|---------------------------------------|----------|
| Theorem 2 R² = 0.94                | `task1_theorem2.json` r²=0.940        | ✅ |
| β_parallel / β_cov ≈ 33×           | 0.0493 / 0.0015 = 32.9×               | ✅ paper says 33× |
| ROC = 0.78 (broader def)           | `task2_stress.json` roc_auc_any=0.78  | ✅ caption disambiguated to "any flip-positive shift type" |
| ROC = 0.62 (narrow def)            | `task2_stress.json` roc_auc_subspace=0.62 | ✅ also reported |
| Concept mag=1.0 DIAL=0.393         | stress_test_summary.tsv                | ✅ |
| DANN-lite target AUC 0.39          | tta methods_comparison.tsv 0.392      | ✅ |
| Scaling slope r²=0.19              | task4_scaling.json r²=0.193           | ✅ "near zero" |
| 4/4 domains flip                   | task5_cross_domain.json all True      | ✅ |
| 8/8 druggable retained (paper-2 numbers) | n/a — v15 is methodology, no biomarker claims | ✅ correctly excluded |

## Reviewer simulation (5 anticipated probes)

### 1. "Theorem 2 Gaussian-shared-Σ assumption — is this tight?"

**Anticipated complaint:** the assumption is restrictive; many real
distributions are non-Gaussian. **Defence:** §6 Limitations explicitly
flags this. Conjecture 1 in Appendix A indicates the IB-identity
extension. Real-data §5.5 4-modality validation shows the *flip
mechanism* reproduces beyond Gaussian (vision, NLP, clinical, biology),
even if the precise constants differ.
**Strength: medium-strong.** Honest scoping is the standard NeurIPS-
acceptable framing.

### 2. "All empirics synthetic — where's the real-data?"

**Anticipated complaint:** modern NeurIPS expects at least one real
dataset. **Defence:** §5.5 4-domain validation; the leak retraction in
§1 references a real biomedical pipeline that DIAL successfully
flagged. Camera-ready can extend to Camelyon17 or DomainNet
(committed in §6 Future Work).
**Strength: medium.** Reviewer may push for camera-ready commitment.

### 3. "DIAL is classifier-dependent — is it really a 'metric'?"

**Anticipated complaint:** DIAL = f(classifier family, CV protocol);
not a property of the data alone. **Defence:** §6 Limitations lists
this. Rebuttal: DIAL is a *protocol-aware* diagnostic, like AUC
itself. Theorem 2 establishes it is invariant up to monotone
score transformations, which is the relevant invariance class.
**Strength: medium-strong.**

### 4. "Foundation-model scaling negative result — did you try
pre-trained X/Y/Z?"

**Anticipated complaint:** the 8 encoders are MLPs / random projections /
PCA; reviewer may want CLIP, DINOv2, or a domain-pretrained encoder.
**Defence:** §5.4 explicit caveat that "shift-aware pre-training"
(not raw scale) is the protective regime. Not claiming foundation
models are useless; claiming raw parameter count alone doesn't
help. Camera-ready can add 1–2 pretrained-encoder rows.
**Strength: medium.** Honest caveat may be enough.

### 5. "v5.1 retraction — is this a honesty problem for the paper?"

**Anticipated complaint:** original v5.1 motivation is retracted;
does the paper still stand? **Defence:** Abstract, §1, §5 all
explicitly disclose. The methodological contribution
(Theorem 2 + leak detector) is *strengthened* by the retraction —
the retraction is precisely the worked example showing DIAL
caught a real leak.
**Strength: strong.** Disclosure = honesty; reviewers respect it.

## Submit-ready verdict

**🟢 GO.** The v15 NeurIPS PDF compiles cleanly, anonymisation passes,
all 9 audit gates pass, and the 5 anticipated reviewer probes have
defensible answers. The single outstanding item is for the author to
verify: (a) confirm NeurIPS 2026 deadline at https://neurips.cc/, and
(b) upload the anonymous code repo to https://anonymous.4open.science
and replace the placeholder URL in the code-availability statement
before submission.

**Estimated submission readiness: 95 %.** Remaining 5 % is the two
manual items above (≤ 30 minutes of work).
