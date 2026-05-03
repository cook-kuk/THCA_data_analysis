# Paper 3 ICI — DIAL Audit Plan

**Working title:** HLA loss, neoantigen architecture, and immune ecotypes define ICI vulnerability in molecularly dark thyroid cancer
**Track:** A — design only. No DIAL execution. No bootstrap, no permutation, no Cox fitting.
**Author:** Seungho Cook
**Status date:** 2026-05-04

**Claim guard:** DIAL audit is the *gating step* that decides which pan-cancer ICI signatures are allowed to enter Paper 3's integrated readiness score in thyroid. Without DIAL, the integrated score has no validity argument.

---

## 1. Context — DIAL framework lineage

- DIAL = Direction-Invariance ALgorithm. Originated in `v52_lodo_finding`: under proper LODO ComBat, all 5 THCA classifiers returned DIAL=0.000 / true biology, falsifying the v5.1 DIAL=0.494 flip artifact.
- Codified into `v18_agentic_research framework` (5 patterns + 6 components).
- Applied in Paper 3 to **pan-cancer ICI signatures** rather than to THCA classifiers — same algebra, different inputs.

**Definition for Paper 3 use.** A signature S is *direction-invariant* across context set C if its sign of association with the chosen endpoint (ICI response or T-cell-inflamed phenotype) is preserved across all leave-one-out subsets of C, with effect size > minimum. Failure modes:
- **Flip** — sign reverses in ≥1 leave-one-out fold.
- **Collapse** — effect size shrinks to noise (|t| < 1) in ≥1 fold.
- **Tissue-flip** — sign holds within trained tissue (e.g., melanoma) but reverses in held-out tissue (here: thyroid).

Only signatures that pass DIAL on the pan-cancer reference set are eligible to enter Paper 3's integrated score in thyroid.

---

## 2. Reference set & endpoints

**Pan-cancer ICI reference (per `paper3_ici_dataset_registry.md` §4):**
- IMvigor210 (UC, atezolizumab) — Tier 1
- Hugo GSE78220 (melanoma, anti-PD-1) — Tier 1
- Riaz GSE91061 (melanoma, nivolumab; pre-treatment only) — Tier 1
- Gide melanoma (PRJEB23709) — Tier 2
- Liu melanoma (dbGaP) — Tier 2 if access granted
- Kim gastric (PRJEB25780) — Tier 3
- Cho NSCLC — Tier 3

**Endpoints:**
- Primary — RECIST objective response (CR/PR vs SD/PD).
- Secondary — PFS (Cox HR).
- Tertiary — T-cell-inflamed phenotype (TIS-positive vs negative) as a *biological* anchor independent of clinical response.

**Held-out tissue = thyroid:**
- TCGA-THCA aggressive subset (PDTC/ATC + Stage IV PTC + dedifferentiated) as proxy for "ICI-eligible thyroid".
- GSE76039 (PDTC + ATC) as second held-out cohort.
- No clinical ICI response endpoint in thyroid — endpoint in held-out is *T-cell-inflamed phenotype* (TIS) as proxy. Direction-invariance from response → TIS is checked first within reference cohorts; only signatures invariant under both are claimed.

---

## 3. Audit pipeline (defined; NOT executed)

```
For each signature S in pan-cancer registry:
    For each cohort c in reference set:
        Compute S on c (ssGSEA / IMPRES rules / Cytolytic GM)
        Fit logistic regression: response ~ S + clinical_covariates
        Record beta_S, se, p

    Direction-invariance check 1 (within-reference, response endpoint):
        sign_consistent = all(sign(beta_S, c) == majority_sign across c)
        effect_size_min  = min(|beta_S, c| > threshold)
        DIAL_response    = sign_consistent AND effect_size_min

    Direction-invariance check 2 (within-reference, TIS endpoint):
        Repeat with TIS_18 score as endpoint
        DIAL_TIS = sign_consistent AND effect_size_min

    Tissue-transfer check (held-out thyroid, TIS endpoint):
        Compute S on TCGA-THCA aggressive + GSE76039
        Fit linear regression: TIS ~ S
        DIAL_thyroid_TIS = (sign(beta_S, thyroid) == majority sign in reference)

    DIAL verdict:
        PASS if DIAL_response AND DIAL_TIS AND DIAL_thyroid_TIS
        FLIP if any check has reversed sign
        COLLAPSE if effect size below threshold
        AMBIGUOUS otherwise

    Record verdict in paper3_ici_DIAL_results.tsv
```

Bootstrap — 1000 resamples per cohort to attach CI to beta_S. Permutation — 1000 label-shuffled fits to attach empirical p.

---

## 4. Signatures audited

From `paper3_ici_signature_registry.md` (binding list):

| Signature ID | DIAL prior risk | Why this risk |
|---|---|---|
| `IFNG_AYERS6` | Low | Robust across tissues in published literature. |
| `TIS_18` | Low | T-cell-inflamed phenotype well-established. |
| `CYTOLYTIC` | Low | 2-gene; very stable. |
| `MHC1_CORE` | Low | Mechanistic — direction predictable. |
| `MHC2_CORE` | **Medium** | Tumor-intrinsic Class II expression in thyroid (HT-overlap signal in Paper 2) may flip thyroid direction. |
| `TLS_CABRITA9` | Low | Established in melanoma; expected to hold in thyroid. |
| `TIDE_LIKE` | **High** | Composite includes M2 + Treg + MDSC; some components flip in thyroid. |
| `IMPRES` | **High** | 15-pair rule trained against melanoma; pairwise rules may not transfer. |
| `EXHAUSTION_INDEX` | **High** | Exhaustion can co-occur with non-inflamed in thyroid; semantic flip likely. |
| `M2_TAM` | **Medium** | Thyroid-specific macrophage biology may invert "M2 = bad" rule. |
| `MDSC_LIKE` | **Medium** | Confounded with neutrophil signal in thyroid. |
| `TREG_CORE` | Medium | Treg in thyroid can mark inflamed-but-controlled state. |
| `EFFECTOR_T` | Low | Stable. |
| `CXCL13_AXIS` | Low | Tumor-reactive proxy stable. |

---

## 5. Statistical thresholds

- Sign-consistency required across ≥6 of 7 reference cohorts (allow 1 outlier).
- Effect-size threshold |beta_S| > 0.2 (logistic, standardized predictor).
- Empirical permutation p < 0.05.
- Tissue-transfer check |beta_S, thyroid TIS regression| ≥ 0.1 with same sign as reference majority.
- Multiple-testing correction — Benjamini-Hochberg across the signature set; report q-values.

---

## 6. Pre-specified failure handling

- **All composite ICI scores fail (TIDE_LIKE, IMPRES, EXHAUSTION_INDEX) →** report this as the headline DIAL finding ("pan-cancer composites do not transfer to thyroid"), use only the surviving signatures (likely IFNG_AYERS6, TIS_18, CYTOLYTIC, MHC1_CORE, TLS_CABRITA9, EFFECTOR_T, CXCL13_AXIS) in Module E. This is a publication-strengthening result, not a setback.
- **MHC2_CORE flips in thyroid →** investigate whether thyrocyte-intrinsic Class II expression (Hashimoto-like) is the driver, citing Paper 2. Report as biology-driven flip with mechanistic interpretation.
- **All signatures pass with same direction →** weakest publication outcome (no thyroid-specific story). Pivot prose toward "thyroid behaves canonically" framing; lean harder on Module C HLA LOH + neoantigen architecture for the headline.

---

## 7. Thyroid-adjusted readiness score

Once DIAL audit closes, define:

```
ICI_VULN_DIAL = sum over signatures S where DIAL(S) == PASS of
                (sign_majority_S * w_S * z_S)
              − sum over Module C terms (HLA loss, neoantigen presented)
```

Weights for surviving signatures = equal (default) and DIAL-effect-size-weighted (sensitivity).

This score is reported at PTC vs PDTC vs ATC, BRAF / RAS / fusion / dark-matter, and TDS-tertile strata.

---

## 8. What this audit does *not* do

- **Does not** establish that thyroid patients with high `ICI_VULN_DIAL` will respond to ICI. That requires thyroid ICI-treated raw RNA-seq (which we do not have). Frame as readiness *prioritization*, not response prediction.
- **Does not** retrofit thresholds after seeing thyroid results — DIAL verdict locks before thyroid scoring.
- **Does not** replace TIDE proper — we audit TIDE-like composite and report TIDE proper (via TIDE server) as supplementary at Track B.

---

## 9. Outputs (filenames to be produced in Track B Module D)

- `project/results/paper3_ici/dial/per_cohort_betas.tsv`
- `project/results/paper3_ici/dial/dial_verdict.tsv` (per-signature PASS/FLIP/COLLAPSE/AMBIGUOUS)
- `project/results/paper3_ici/dial/sign_consistency_heatmap.pdf`
- `project/results/paper3_ici/dial/thyroid_transfer_betas.tsv`
- `project/results/paper3_ici/dial/audit_report.md`

---

## 10. Sequence within 12-week plan

DIAL audit runs Wk 9, after bulk ecotype (Wk 3–4), scRNA atlas (Wk 5–6), and HLA/neoantigen (Wk 7–8). Pan-cancer ICI cohorts must be acquired by Wk 7. If a cohort lags, downgrade DIAL to remaining cohorts and document.

---

Track A completed. No marathon violation.
