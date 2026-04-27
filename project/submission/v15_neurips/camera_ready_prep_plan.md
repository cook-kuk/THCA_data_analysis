# v15 — Camera-Ready Prep Plan (if accepted)

_Generated 2026-04-27. Activates only if NeurIPS 2026 sends an
acceptance decision in September 2026. Pre-committed work to
deliver by camera-ready deadline (~October 2026)._

## Binding commitments from rebuttal_prep_v1.md

The following promises were made in the anticipated rebuttal
responses. If the paper is accepted, these become camera-ready
deliverables:

### CR-1. New §6.4 "How to report DIAL in practice"

**Source:** Probe 3 (R2 reviewer concern about classifier-dependence).

**Content (≤ 200 words + 1 small box):**
- DIAL must be reported with the (classifier, CV protocol, feature
  pre-processing) triple, analogous to AUC reporting requirements
- A 5-row reporting checklist for biomedical pipelines:
  1. classifier family + hyperparameters
  2. cross-validation protocol (must be LODO for cross-cohort claims)
  3. batch-correction operator + fit protocol (train-only vs pooled)
  4. n_seeds + variance reporting
  5. DIAL ranges + interpretation thresholds
- Drop-in template for FDA / EMA audit logs (CSV schema)

**Effort:** 4 hours.

### CR-2. New §5.6 "Real-data benchmark"

**Source:** Probe 2 (R2 reviewer concern about synthetic-only empirics).

**Content:**
- Camelyon17 (medical imaging, n ≈ 400, 5 hospital sources):
  apply DIAL to the LR pipeline that already exists in the WILDS
  benchmark (Koh et al. 2021). Expected: confirmable flip on at
  least one cross-hospital pairing.
- DomainNet (n ≈ 600k, 6 domains): apply DIAL to standard ResNet-18
  feature extractor + linear classifier. Expected: positive flip
  on the painting → real domain pair (highest baseline shift).
- Single new figure (3 panels): Camelyon17 by-hospital DIAL,
  DomainNet pairwise DIAL heatmap, summary table.

**Effort:** 18 hours (12 GPU hours + 6 hours analysis / writing).
Hardware: borrow 1 NVIDIA L4 from a colleague or use Lambda Labs
for ~$8 (cheaper than effort-equivalent CPU).

### CR-3. Extended §5.4 with foundation-model rows

**Source:** Probe 4 (R3 reviewer concern about foundation-model coverage).

**Content:**
- Add 3 rows to Table 4: DINOv2-small (frozen), CLIP-ViT-B/32
  (frozen text encoder applied to NLP-shaped data), BiomedCLIP
  (frozen, on biology-shaped data)
- Updated narrative: "shift-aware pre-training partially protects
  on the BiomedCLIP row (DIAL drops to X) but not on DINOv2 frozen
  (DIAL stays Y) — confirming that domain-aligned pre-training is
  the necessary protective property, not raw scale or generic
  vision objectives."

**Effort:** 6 hours (small synthetic n; CPU-feasible with HuggingFace
encoder API).

### CR-4. New Appendix E "Scope of applicability"

**Source:** Probe 1 (R1 reviewer concern about Gaussian-Σ tightness).

**Content:**
- Explicit table distinguishing **proven** regimes (Gaussian
  shared-Σ + linear T) from **observed** regimes (the 4 modalities
  of §5.5 + the new §5.6 real-data rows)
- One non-Gaussian counter-example: heavy-tailed conditional
  (Cauchy class-conditionals) where DIAL is observed positive but
  the explicit form of Φ is unknown — frames the gap honestly
- Lists three open problems for future work

**Effort:** 6 hours (1 small synthetic experiment + writing).

### CR-5. Re-anonymisation / de-anonymisation

- De-anonymise author block + email
- Add Acknowledgements (PI 유 교수님, v5.2 audit team, anonymous
  reviewers for camera-ready feedback)
- Switch anonymous GitHub link to public release URL
- Re-verify PDF metadata (set `\pdfauthor{Seungho Cook}` etc.)

**Effort:** 1 hour.

---

## Total camera-ready effort

| Item | Hours | Calendar |
|---|---:|---|
| CR-1 reporting protocol | 4 | week 1 |
| CR-2 real-data benchmark | 18 | weeks 1-3 (overlapping with GPU access) |
| CR-3 FM rows | 6 | week 2 |
| CR-4 scope appendix | 6 | week 3 |
| CR-5 de-anonymisation | 1 | week 4 |
| QA / final compile | 4 | week 4 |
| **Total** | **39 h** | 4-week buffer |

NeurIPS gives ~3-4 weeks for camera-ready. **39 hours fits within a
3-week sprint at ~13 hr/week** — manageable alongside other work.

---

## Hard constraints to monitor

1. **Page limit** — NeurIPS camera-ready often expands to 10 pages
   (vs 9 main for submission). Verify NeurIPS 2026 specifies 9 vs
   10 main + unlimited appendix.
2. **Anonymity check at OpenReview** — the de-anonymised version
   must pass OpenReview's automated PDF-metadata strip.
3. **Code release timing** — public GitHub URL must exist before
   camera-ready upload (~1 day before).
4. **Conflict-of-interest forms** — NeurIPS asks accepted authors
   for an updated COI form 1 week before camera-ready.
5. **DOI / ArXiv** — recommended to also post a non-anonymous
   ArXiv preprint when camera-ready is accepted, for citation
   propagation.

---

## What does NOT change for camera-ready

- Theorem 1 statement and proof (Appendix A is final)
- The 5 v15 main experiments (§5.1–§5.5) — only §5.6 is added
- Bibliography (52 entries — may add ~3-5 if new comparisons cite
  new works, but stable)
- Limitations narrative — already honest, no need to soften
- v5.1 retraction disclosure — must remain explicit; this is a
  feature

---

## Risk: not accepted

If NeurIPS 2026 returns a reject decision (~70% probability per
self_review_simulation.md), the camera-ready effort plan above
becomes the **revision plan** for the next venue. See
`fallback_venues_plan.md` for the sequencing.

The 39-hour upgrade is venue-portable: it strengthens the paper
for **any** subsequent submission (ICML 2027, TMLR, AISTATS 2027).
None of the work is wasted regardless of NeurIPS outcome.
