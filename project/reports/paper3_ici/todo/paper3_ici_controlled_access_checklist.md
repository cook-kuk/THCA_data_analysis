# Paper 3 ICI — Controlled-Access Application Checklist

**Status:** Application work permitted ONLY in spare cycles that do not displace Paper 1/2 marathon work. Long-lead applications (4–8 weeks) may begin in advance because submission ≠ compute. This file is editable; check items off as you submit.
**Author:** Seungho Cook
**Created:** 2026-05-04

**Marathon-compatibility rule:** Submitting a dbGaP / EGA application is permitted as a *background administrative task* only if it takes <30 min of focused effort per session. Compute / data download / analysis remains blocked. If an application requires significant time (IRB letter, data-use agreement drafting), defer to Track B kickoff.

---

## A. TCGA / GDC controlled-access (paired tumor-normal BAM)

**Why:** LOHHLA needs paired BAM, not just MAF. MC3 MAF (open-access) suffices for somatic-MAF-based neoantigen prediction.

- [ ] Decide: do we actually need paired BAM, or is MC3 MAF + arcasHLA-on-RNA-BAM sufficient?
  - If MC3 MAF + RNA BAM sufficient → **no controlled application required**, mark this section closed.
  - If LOHHLA on paired WES BAM is desired → continue.
- [ ] Confirm institutional eRA Commons / dbGaP eligibility (PI / signing official identified)
- [ ] Draft Data Access Request — Research Use Statement aligned with Paper 3 working title
- [ ] IRB / ethics letter (if institutional)
- [ ] Local Data-Use Certification signature
- [ ] Submit dbGaP application phs000178 (TCGA general)
- [ ] Wait for DAC decision (~2–6 weeks)
- [ ] Receive approval; configure GDC token

**Estimated lead time:** 4–8 weeks. Start at Track B Wk -8 (8 weeks before Wk 1) if LOHHLA is desired.

---

## B. EGA controlled-access — Yoo SK 2019 Korean ATC (if EGA-hosted)

**Why:** Korean ATC / advanced DTC cohort for Asian dedifferentiation generalization (dataset registry §1 Tier 1).

- [ ] Verify whether Yoo 2019 is EGA-controlled or open (check supp + EGA portal — defer to Track B kickoff)
- [ ] If EGA-controlled:
  - [ ] Identify DAC contact
  - [ ] Draft Data Access Agreement
  - [ ] Submit application
  - [ ] Wait for DAC decision (~4–6 weeks typical)
- [ ] If denied: drop to Tier 2; substitute with TCGA aggressive subset; document in `paper3_ici_dataset_registry_VERIFIED.tsv` at Track B Wk 1.

**Estimated lead time:** 4–6 weeks if controlled.

---

## C. dbGaP — Liu et al. melanoma (phs000452) — DIAL audit Tier 2

**Why:** Pan-cancer ICI reference cohort for direction-invariance audit (registry §4).

- [ ] Verify accession current (phs000452 or successor)
- [ ] Confirm whether prior dbGaP approval covers melanoma study
- [ ] If new application required:
  - [ ] Same TCGA-style flow (Section A)
  - [ ] Submit
  - [ ] Wait
- [ ] If denied or timeline > 8 weeks: drop to Tier 3; retain Tier 1 melanoma + UC anchors (IMvigor210, Hugo, Riaz) for DIAL.

---

## D. Other potentially-controlled cohorts (verify-then-decide)

- [ ] **Cho NSCLC ICI** — pin accession at Track B kickoff; controlled status unknown
- [ ] **Han 2024 JCI Insight scRNA** — pin accession at Track B kickoff; access mode unknown
- [ ] **Ning / Liao / Zheng spatial cohorts** — pin accessions at Track B kickoff

For each: verify access mode at Track B Wk 1 (NOT during marathon). Application submission deferred to Track B kickoff unless explicitly cleared.

---

## E. IMvigor210 — open-access R-package

**Why:** DIAL Tier 1 anchor (UC, atezolizumab).

- [ ] Confirm `IMvigor210CoreBiologies` R-package is still distributable (Roche / Mariathasan 2018)
- [ ] No application needed — open-access bundle
- [ ] Reserve R environment / renv lockfile for Wk 1

---

## F. Application status board (single line per app)

| Cohort | Status | Submitted date | Decision date | Notes |
|---|---|---|---|---|
| TCGA paired BAM | not_yet_submitted | — | — | Decide need first (Section A) |
| Yoo 2019 Korean ATC | verify_first | — | — | EGA status unknown |
| Liu melanoma | not_yet_submitted | — | — | Tier 2; degrade-tolerable |
| Cho NSCLC | verify_first | — | — | |
| Han 2024 scRNA | verify_first | — | — | |
| Ning / Liao / Zheng spatial | verify_first | — | — | |

Update this table as items move; do NOT delete rows.

---

## G. Marathon-displacement guard

- [ ] Each application session ≤ 30 min focused effort
- [ ] No application work during scheduled Paper 1 writing blocks
- [ ] If any application requires Paper 1/2 attention to draft (e.g., need to cite Paper 1 in research-use statement), defer until Paper 1 bioRxiv submission
- [ ] Track time spent on applications; if cumulative > 4 hr/week, halt and re-decide

---

Paper 3 Track A frozen. Return to Paper 1/2 marathon.
