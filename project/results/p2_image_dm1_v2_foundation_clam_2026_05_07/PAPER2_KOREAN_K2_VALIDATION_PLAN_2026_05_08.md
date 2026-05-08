# Paper 2 — Korean K2 External Validation Plan

- Document date: 2026-05-08
- Owner: Seungho Cook
- Scope: External validation strategy for Paper 2 H&E → DM1 image classifier (TCGA-trained) on independent Korean cohort
- Status: PLAN — execution pending H&E availability verification (Section 3)
- Related artifacts: `PAPER2_OUTLINE_2026_05_08.md`, `SPRINT_PLAN.md` (this directory)

---

## 1. Aim

- Validate the H&E → DM1 image classifier (TCGA-trained UNI/ViT-L + CLAM, best-fold weights) on an independent Korean PTC cohort
- Quantify generalization across staining/scanner/population shift (TCGA US multi-center vs Korean single/few-center)
- Pre-specified primary metrics (locked before inference):
  - Slide-level AUC with 95% bootstrap CI (1000 reps, stratified)
  - Calibration: Brier score, calibration intercept, calibration slope, reliability curve
  - Decision-curve analysis: net benefit across thresholds 0.1–0.9
  - Sensitivity, specificity, PPV, NPV at TCGA Youden-optimal threshold (frozen pre-inference)
- Secondary: attention-map qualitative review for biological plausibility (lymphoid aggregates, tumor borders) on Korean tissue

---

## 2. Cohort: Korean K2 (PRJEB11591)

- n = 260 PTC samples, RNA-seq + clinical metadata
- Source: Yoo 2016, SNU-GMI (Seoul National University Gangbuk Genomic Medicine Institute)
- Public access: ENA accession PRJEB11591
- Already replicated locally:
  - arcasHLA RNA-seq HLA imputation (n=260, see memory `v17_arcasHLA_korean_k2`)
  - DPB1*05:01 56% replication of Kim 2014 finding
- Critical disambiguation (memory `v17_K2_vs_bundang_distinction`):
  - K2 / KOREAN_K2 = PRJEB11591 Yoo 2016 SNU-GMI public dataset
  - Bundang SNUH = outreach-stage, NO data in hand — do not confuse in plan or paper text

---

## 3. Critical question: Does Korean K2 have FFPE H&E available?

- Status: UNKNOWN as of 2026-05-08 — verification is the gating action
- Expected outcome: most likely H&E NOT bundled with PRJEB11591 (RNA-seq study, not pathology study)

### 3.1 Plan A — ENA / PRJEB11591 supplementary check
- Inspect PRJEB11591 ENA project page for any image / pathology related files or linked archives
- Check Yoo 2016 paper supplements for slide availability statements
- Check NCBI BioStudies, Figshare, Zenodo for companion deposits keyed on study accession
- Deliverable (2026-05-09 EOD): `K2_HE_AVAILABILITY_REPORT.md` with one of three verdicts: AVAILABLE / PARTIAL / NOT_AVAILABLE

### 3.2 Plan B — Corresponding author contact
- Target: Yoo SK et al. corresponding author (email TBD from paper masthead)
- Ask: (a) FFPE block / scanned WSI availability for the PRJEB11591 cohort, (b) sharing terms (DTA, IRB), (c) any subset already digitized
- Cc: Yu lab via existing channel
- Deliverable: contact log entry with response or non-response by 2026-05-15

### 3.3 Plan C — Bundang SNUH collaboration request
- Already in outreach; longer timeline (institutional approval, DTA, IRB)
- NOT counted as primary external cohort for original Paper 2 submission
- Hold for revision-stage validation (post-submission)

### 3.4 Decision tree
- Plan A AVAILABLE → execute Section 4 in full; K2 becomes primary external cohort
- Plan A PARTIAL (e.g., subset of n) → execute Section 4 on subset; report as "partial external validation, n=k of 260"
- Plan A NOT_AVAILABLE + Plan B positive within 7 days → schedule WSI transfer, execute Section 4 on delivered subset
- Plan A NOT_AVAILABLE + Plan B null → execute Section 5 fallback; submit Paper 2 with TCGA-only and pre-registered K2/Bundang validation as revision deliverable

---

## 4. If H&E available — execution plan

### 4.1 Pre-processing
- WSI download → tile extraction using identical Phase 2 pipeline (reproducibility lock)
- Tile spec: 256×256 px @ 20× (matches TCGA training)
- UNI/ViT-L feature extraction with frozen weights (no re-training, no fine-tune)
- Quality control:
  - Tissue mask area distribution (flag slides with <5% tissue)
  - Tile count per slide (flag <50 or >5000)
  - Stain intensity histogram vs TCGA reference (flag >2σ outliers)
  - Optional Macenko/Vahadane stain normalization (pre-specify ON or OFF before inference; recommend ON given known TCGA→Korean staining drift)
- Deliverable: `K2_tiles_manifest.csv` (slide_id, n_tiles, tissue_pct, stain_QC_flag)

### 4.2 Label generation
- DM1 / DM2 sub-cluster labels derived from K2 RNA-seq
- CRITICAL: use within-sample-centered 8-gene profile, NOT absolute-form TCGA-trained LogReg (memory `v17_korean_k2_calibration` — TPM inflates 10–100× and absolute classifier is wrong-direction)
- Sub-clustering pipeline replicates Paper 1 K2 cluster definition (link exact script + commit hash before run)
- Cross-tab K2-derived DM1 vs Korean ethnicity-specific HLA features (AFND South Korea allele frequency pool)
- Deliverable: `K2_DM1_labels.csv` (slide_id, dm1_label, dm1_prob_rna, hla_features) with version-locked centering script

### 4.3 Inference
- Load TCGA-trained CLAM model (best-fold weights, frozen, SHA-256 logged)
- Per-slide DM1 probability + attention map
- Frozen pre-inference artifacts:
  - TCGA Youden-optimal threshold τ*
  - Decision-curve threshold grid 0.1–0.9 in 0.05 steps
- Compute on K2:
  - AUC + 95% CI (bootstrap 1000 stratified)
  - Sensitivity, specificity, PPV, NPV at τ*
  - Brier, calibration intercept, calibration slope, reliability curve (10 bins)
  - Decision-curve net benefit vs treat-all / treat-none
- Deliverable: `K2_inference_metrics.json` + figures/tables

### 4.4 Robustness checks (pre-specified subsets)
- Hashimoto-overlap status (PTC-HT vs PTC-only) — Korean PTC-HT subset enriched per GSE286332 priors
- Stage T1 vs T2 vs T3–T4
- Molecular driver: BRAF vs RAS vs driver-negative (NBNR)
- Report subset AUC + CI; flag any subset with n < 20 as exploratory only
- DM1 sub-B (mutation-negative) slice: explicit cross-link to K2 NBNR finding (memory `v17_D6P7_dm1_subB_NBNR`)

---

## 5. If H&E NOT available — fallback

- Plan B (already covered in 3.2): GSE286332 Korean PTC vs PTC+HT (n=18) — too small for AUC validation; usable only as qualitative attention-map sanity check IF FFPE H&E exists
- Plan C: Bundang SNUH outreach — Paper 2 revision-stage supplemental study, NOT submission deliverable
- Honest reporting in original Paper 2:
  - Limitations section: "External validation pending; submission cohort is TCGA-only"
  - Pre-register K2/Bundang validation protocol publicly (OSF or GitHub-tagged release) before submission to lock revision response credibility
- Alternative external surfaces to explore in parallel (low priority):
  - CPTAC thyroid pathology slide set (US, but independent of TCGA training set)
  - GTEx thyroid normal tiles as negative-control distribution check

---

## 6. Timeline

- Week 1 (2026-05-08 → 2026-05-14): Section 3 H&E availability verification (Plan A + Plan B initiated)
- Week 2 (2026-05-15 → 2026-05-21): WSI download + Section 4.1 pre-processing on whatever subset is available
- Week 3 (2026-05-22 → 2026-05-28): Section 4.2 label generation + Section 4.3 inference
- Week 4 (2026-05-29 → 2026-06-04): Section 4.4 robustness + figure/table assembly
- Week 5 (2026-06-05 → 2026-06-11): Write external validation supplement (or revision-ready document if H&E unavailable)
- Marathon deadline 2026-06-13: original Paper 2 submission with TCGA-only primary + K2 status statement matched to actual outcome

---

## 7. Risks

- H&E unavailable in K2 — most likely; mitigated by Plan B/C and pre-registered revision protocol
- Staining / scanner batch effect (TCGA vs Korean protocols) — mitigated by stain normalization in 4.1; if AUC drops materially without normalization but recovers with it, report both
- Small Hashimoto-overlap subset in K2 → underpowered for sub-stratification — flag as exploratory, do not over-claim
- Korean WSI access permissions (institutional approval, DTA, IRB) — start in parallel with Plan A even before availability confirmed, to avoid serial blocking
- Label noise from RNA-derived DM1 — mitigated by within-sample-centered profile (calibration memory) + sensitivity analysis around label confidence cutoff
- Time pressure vs marathon 6/13 — fallback (Section 5) keeps submission on schedule even if validation slips

---

## 8. Success criteria (pre-specified)

- Primary: AUC ≥ 0.70 with 95% CI lower bound > 0.55 (above ResNet50 closure baseline established in Phase 2)
- Calibration slope in [0.8, 1.2]; intercept |b| < 0.2
- Decision curve: net benefit > both treat-all and treat-none across threshold range 0.3–0.7
- Robustness: no single pre-specified subset shows AUC < 0.60 with adequate n; if it does, flagged transparently and discussed
- Failure handling: if AUC < 0.65, do NOT hide — report, attribute to identifiable shift (stain, scanner, population), and propose targeted fix (domain adaptation) as future work

---

## 9. Reporting structure (paper)

- Section position: Paper 2 §3.4 "External validation (Korean K2)" if H&E secured by submission; otherwise Supplement Section S4 stub with pre-registered protocol
- Tables:
  - Table E1: external cohort characteristics (n, age, sex, stage, driver, HT status)
  - Table E2: TCGA vs K2 AUC / sens / spec / PPV / NPV side-by-side
  - Table E3: calibration metrics (Brier, intercept, slope)
- Figures:
  - Figure E1: external ROC curve with TCGA reference overlay
  - Figure E2: calibration / reliability plot
  - Figure E3: attention heatmap exemplars (DM1 high-prob, low-prob, FN, FP)
  - Figure E4 (optional): decision-curve net benefit
- Code + frozen weights + threshold τ* released alongside paper for reproducibility

---

## 10. Open questions requiring user decision

- Q1: Approve direct contact to Yoo 2016 corresponding author (Plan B email) — yes/no, and whether to cc Yu lab
- Q2: Stain normalization default ON or OFF for primary inference (recommend ON; needs lock before inference run)
- Q3: If H&E unavailable, accept "TCGA-only + pre-registered revision protocol" framing for submission, or hold submission past 6/13 deadline
