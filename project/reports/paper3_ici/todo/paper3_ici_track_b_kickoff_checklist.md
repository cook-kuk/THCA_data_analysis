# Paper 3 ICI — Track B Kickoff Verification Checklist

**Status:** Track B BLOCKED. This file is editable; check items off only when actually executed at Track B kickoff. Track A design bundle (`../PAPER3_ICI_TRACK_A_BUNDLE.md`) is frozen and read-only.
**Author:** Seungho Cook
**Created:** 2026-05-04 (frozen design date)

**Purpose:** Single-pass verification that all entry conditions are satisfied before any Paper 3 compute begins. If any unchecked item remains, Track B does not start.

---

## A. Hard entry gates (G1–G6 from go/no-go verdict §2)

- [ ] **G1** — Paper 1 bioRxiv submission confirmed (DOI received)
- [ ] **G2** — Paper 2 Task A/B/C closure confirmed (per `v18_paper2_HT_isolated.md`)
- [ ] **G3** — TCGA dbGaP / GDC controlled access decided (granted, denied, or confirmed-not-needed for the chosen scope)
- [ ] **G4** — At least 3 of 6 scRNA cohorts confirmed accessible
- [ ] **G5** — At least 3 of 5 Tier 1/2 pan-cancer ICI cohorts confirmed accessible
- [ ] **G6** — Explicit user command "Track B 시작" recorded (not "고고" / "다 해줘", per `v17_sprint_vs_marathon_violation.md`)

If any G1–G6 unchecked: **STOP**. Do not proceed to Section B.

---

## B. Design-bundle re-read (binding before any compute)

- [ ] Re-read `../PAPER3_ICI_TRACK_A_BUNDLE.md` Section 1 (dataset registry) end-to-end
- [ ] Re-read Section 2 (signature registry)
- [ ] Re-read Section 3 (scRNA atlas plan)
- [ ] Re-read Section 4 (HLA & neoantigen feasibility)
- [ ] Re-read Section 5 (DIAL audit plan)
- [ ] Re-read Section 6 (figure plan)
- [ ] Re-read Section 7 (go/no-go verdict)
- [ ] Re-read Section 8 (12-week execution plan)
- [ ] Re-read memory `v19_paper3_ici_track_a.md`
- [ ] Confirm "ICI vulnerability / readiness / immunogenomic prioritization" framing — verify no prose drift toward "ICI response predictor"

---

## C. Cross-paper boundary sanity (one-pass before any compute)

- [ ] Paper 1 (DM1 molecular) — confirm Paper 3 will NOT duplicate DM1 cluster definition; will reuse Paper 1 driver / fusion / TDS calls only
- [ ] Paper 2 (H&E-DM1 / TCGA-validation) — confirm Paper 3 will NOT duplicate PTC+HT TLS / AICDA / BCR analyses; will reference at most as a comparator subgroup
- [ ] Paper 4 (GD HLA backlog) — confirm Paper 3 will NOT introduce cookHLA SNP imputation; will not introduce GD/Graves' disease analyses
- [ ] Symmetric boundary lines drafted for Paper 3 Discussion (referencing P1/P2/P4)

---

## D. Compute environment readiness

- [ ] Azure burst credits available (`v17_arcasHLA_korean_k2` $4.80 pattern reference)
- [ ] GPU instance availability for scVI training (1–2 GPU-day budget)
- [ ] CPU instance availability for LOHHLA / Mutect2 (200–400 CPU-hr)
- [ ] NetMHCpan / NetMHCIIpan licenses installed (academic — DTU)
- [ ] OptiType / Polysolver / xHLA / arcasHLA installed and version-pinned
- [ ] LOHHLA installed and version-pinned
- [ ] scvi-tools / scanpy / Seurat installed and version-pinned
- [ ] celltypist / ProjecTILs reference atlases downloaded
- [ ] Output directory `project/results/paper3_ici/` created with subtree from registry / atlas / hla / neo / dial / integrated

---

## E. Internal source-of-truth handshake

- [ ] Paper 1 ETL driver/fusion calls accessible read-only
- [ ] TCGA-THCA harmonized RNA-seq matrix accessible read-only
- [ ] PRJEB11591 arcasHLA results accessible read-only (per `v17_arcasHLA_korean_k2.md`) — DO NOT re-run
- [ ] GSE286332 (PTC+HT) accessible read-only — Paper 2 boundary discipline (no re-analysis here)
- [ ] Paper 1 spatial pipeline outputs accessible read-only

---

## F. Statistical / reporting standards lock

- [ ] Significance: q < 0.05 BH-corrected
- [ ] Effect size minimum threshold confirmed (|beta| > 0.2 standardized for DIAL)
- [ ] Bootstrap × 1000 standard for CI
- [ ] Permutation × 1000 standard for empirical p
- [ ] Score sign convention locked (no silent sign-flips)
- [ ] Cohort z-score before pooling (no platform-mean leakage)
- [ ] Multiple-testing pre-registered

---

## G. Stop conditions during Track B (kill switches K1–K6 from verdict §4)

Re-confirm understanding (no checkbox — these are continuous):
- K1 — never claim "ICI response predictor" without thyroid ICI raw RNA-seq
- K2 — Paper 1 ship slip → halt
- K3 — DIAL zero PASS → pivot to feasibility-only paper, re-decide
- K4 — Paper 1/2/4 cross-paper contamination → halt + rewrite
- K5 — driver-call disagreement with Paper 1 ETL → halt + reconcile
- K6 — user pivot priority → halt + re-decide

---

## H. Final pre-Wk-1 sign-off

- [ ] All A–F items checked
- [ ] User explicit "Track B 시작" recorded (date + time)
- [ ] This checklist file moved to `archive/track_b_kickoff_checklist_<YYYYMMDD>.md`
- [ ] Wk 1 of `paper3_ici_12week_execution_plan.md` begins

---

Paper 3 Track A frozen. Return to Paper 1/2 marathon.
