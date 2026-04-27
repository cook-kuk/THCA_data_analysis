# v17p35 Sprint — Situation Report for Web Claude

_Written 2026-04-27 09:00 KST. Self-contained. Read this cold and you should be able to answer "what's the right next sprint?" without further context fetching._

---

## 0 · Meta

- **User**: Seungho Cook (kukshomr@gmail.com). PI = 유 교수님. Project = Thyroid cancer transcriptomic subtyping → drug actionability → paper.
- **Project root**: `/opt/thyroid-dash/project/`
- **Two AI agents collaborating on the same repo**: (a) Claude Code CLI ("나" — the agent writing this doc), (b) Codex CLI ("친구" — running asynchronously, started ~2026-04-26 22:00 KST). They write to the same `notebooks_or_scripts/` and `results/`. No merge conflicts so far because they hit different filenames.
- **Goal venue cascade (honest)**: npj Precision Oncology primary (target ~70%), Genome Medicine secondary (~40%), Nature Communications reach (<15%). Bioinformatics + JCO Precision Oncology as fallback split paper.
- **Current paper draft scope**: 7 main figures, ~5,000 words, `reports/v17p3/v17p3_paper_outline_FINAL.md` (currently 22-line stub — see §6 below).
- **Today's situation in one sentence**: The codex friend completed Tier 1 (FIX1–FIX5) of the v17p35 12-task sprint, then stopped at 08:50:30 KST. Tier 2 (AMP1–AMP5), Tier 3 (SYNTH1–SYNTH2), Tier 4 (DRAFT1–DRAFT4) — 11 tasks — are entirely unstarted. Friend is asking which AMP to do first.

---

## 1 · Project context (one paragraph each)

### 1.1 Biology
PTC (papillary thyroid carcinoma) is canonically split into BRAF-like vs RAS-like via BRS52 expression signature (Chakravarty 2011). v17 sub-clustering of the **driver-negative** PTC fraction discovered a finer transcriptomic axis we now call **DM1 / DM2** (Dark Matter 1 / 2). DM1 is MAPK-active, hot-tumour-immune, IFN-γ/TNF-α/Inflammatory enriched, ATC-trajectory-proximal. DM2 is OxPhos-enriched, differentiation-marker-high (TPO/DIO1/SLC5A8), RAI-uptake-high. Critically, DM1/DM2 are now shown to be **largely orthogonal to BRAF/RAS mutation status** (v17p3 A2): 281 BRAF-mutant tumours → 99% DM1, 54 RAS-mutant → 72% DM2, but 17 outliers (2 BRAF/DM2-like, 15 RAS/DM1-like) violate the canonical map.

### 1.2 v17 phase history
- **v17 Phase 1** (`results/v17/`, 65 files): Initial dark-matter discovery + driver landscape + trajectory + DIAL audit.
- **v17 Phase 2** (`results/v17p2/`, 48 files): Robustness (consensus matrix, bootstrap), proxy GSEA, clinical, external proxy.
- **v17 Phase 3** (`results/v17p3/`, 67 files): A1–A6 advanced analyses + F1–F3 fixes (external recovery, proper GSEA, alternative endpoints).
- **v17p3.5 / v17p35** (this sprint, current): Fix 4 broken tasks, amplify with 5 new angles, synthesise master figures, draft paper. **Spec is in the user's most recent message; orchestrator is at `notebooks_or_scripts/v17p35_orchestrator.py`.**

### 1.3 Honest novelty position (from prior-art check 2026-04-25)
Already established literature that constrains our claims:
- **Liu 2018 (PMID 31949805)** — TROP2 ↔ BRAF V600E in PTC, IHC.
- **Bychkov 2018 (PMID 29228520)** — TROP2 prognostic in PTC.
- **Kalfert 2024 (PMID 38696857)** — BRAF × TACSTD2 mRNA in PTC + paired LNM.
- **Nieto-Jiménez 2023 (PMID 37740463)** — sacituzumab govitecan flagged for thyroid as niche indication.
- **Dum 2022 (PMID 35477165)** — TROP2 TMA n=18,563 with PTC nodal-metastasis link.
- **Grothey 2021 (PMID 33836264)** — BRAF V600E mCRC + irinotecan = resistance polarity (counter-polarity citation; thyroid finding is opposite).

What we genuinely contribute (after priority reductions):
1. **DM1/DM2 axis itself** — the BRAF/RAS-orthogonal transcriptomic split with a new biology readout (immune state + RAI gradient + TF circuit).
2. **PRISM-level SN-38 BRAF-selective killing in thyroid** — opposite polarity to CRC. Novel pharmacology observation.
3. **Cell-line corroboration of BRAF×TROP2 in CCLE thyroid lines** (top 3 TACSTD2 lines all V600E).
4. **BRAF sub-stratification of TROP2 ADC trials** (NCT06235216 / NCT07521670 are BRAF-blind).
5. **Two-gate biomarker (BRAF-like classifier + TACSTD2 z-score)** — no prior paper proposes this combination.

Verdict files: `reports/v14_novelty_verdict.md` + `reports/v14_novelty_verdict_addendum.md`.

---

## 2 · v17p35 sprint — task spec

User's prompt defines 14+ tasks across 4 tiers with explicit DAG: Tier 1 (parallel) → Tier 2 (parallel after T1) → Tier 3 (synth after T2) → Tier 4 (draft after T3).

| Tier | Task | Spec headline | Paper figure target |
|---|---|---|---|
| 1 | FIX1 | Drug response (PRISM × CCLE × DM1/DM2) — PERCEPTION-style | Figure 7 (Genome Med 결정타) |
| 1 | FIX2 | CNV/Aneuploidy (Taylor 2018) + Thorsson 2018 immune subtype | Figure 5 보강 |
| 1 | FIX3 | Cox multivariable + alternative endpoints (PFI/DFI/LN+/Mets/ATA) | Figure 2 |
| 1 | FIX4 | Proxy vs proper GSEA agreement | Supplementary |
| 1 | FIX5 | scRNA per-patient + immune cell type breakdown + per-cell DM signature | Figure 5 |
| 2 | AMP1 | BRAF/DM2-like 2명 + RAS/DM1-like 15명 outlier deep dive | Figure 4 (paper core, BRAF/RAS 직교성 결정타) |
| 2 | AMP2 | A6 PDTC RAI score reframe (perfect separation, opposite direction) | Figure 3 (paper hook) |
| 2 | AMP3 | F2 + A1 + A5 immune integrated "Hot/Cold landscape" + TIDE | Figure 5 (paper highlight) |
| 2 | AMP4 | TF circuit (PAX8/NKX2-1/FOXE1) + 8-gene RAI prediction model | Figure 6 (npj 결정타) |
| 2 | AMP5 | Pan-cancer DM signature transfer to LUAD/COAD/LGG/SKCM | Supplementary |
| 3 | SYNTH1 | Master figure composer (7 main + 15 supplementary) | Submission ready |
| 3 | SYNTH2 | Korean HTML gallery v3 (Pretendard, Plotly, mobile responsive) | Internal review |
| 4 | DRAFT1 | Paper outline FINAL + Title × 5 + Abstract × 3 versions | Manuscript |
| 4 | DRAFT2 | Cover letter × 3 venue (Bioinformatics, npj, Genome Med) | Submission |
| 4 | DRAFT3 | Reviewer defense (15-20 attacks × rebuttal paragraphs) | Pre-submission |
| 4 | DRAFT4 | External review dump v3 | Internal review |

---

## 3 · What the codex friend actually delivered (audit)

### 3.1 Friend's timeline (orchestrator log)
```
~2026-04-26 22:00 KST   friend starts working (created v17p3_*.py, v17p2_*.py)
2026-04-27 02:02:51     v17p3 A3 done {n_drugs_fdr<0.1: 0, best_dm1_drug: null}    ← EMPTY
2026-04-27 02:03        v17p3 paper_outline_FINAL.md, response_to_reviewers_v2.md, dump_v3.md written
2026-04-27 08:21:16     AUDIT_full_state.md written (friend's own meta-audit)
2026-04-27 08:43        v17p35 sprint started — FIX4 first
2026-04-27 08:49–08:50  FIX1, FIX2, FIX3, FIX5 all completed
2026-04-27 08:50:30     friend stops. No further activity.
```

### 3.2 Friend's v17p35 deliverables — Tier 1 only

**Scripts written** (8):
- `v17p35_common.py` (15 lines, just paths)
- `v17p35_orchestrator.py` (41 lines — TIER1 list contains only `FIX4`; FIX1/2/3/5 were called separately by hand)
- `v17p35_FIX1_drug_response.py`
- `v17p35_FIX2_genomic_immune.py`
- `v17p35_FIX3_clinical_full.py`
- `v17p35_FIX4_gsea_compare.py`
- `v17p35_FIX5_scrna_full.py`

**Tier 2/3/4 scripts**: none. AMP1–AMP5, SYNTH1–SYNTH2, DRAFT1–DRAFT4 do not exist.

**Tables produced** (22): per-task TSVs in `results/v17p35/tables/`.

**Figures produced** (15): per-task HTMLs in `results/v17p35/figs/`.

**Dashboard**: `reports/v17p35/index.html` exists but is a **56-line placeholder** (2,779 bytes), not the SYNTH-2 spec deliverable.

### 3.3 Substance-level audit of each FIX

| FIX | Status | Substance | Hole |
|---|---|---|---|
| **FIX1** | ✅ Output exists, but flawed | 2 drugs FDR<0.1 (clonazepam, baicalein, propofol, luteolin top of list) | **`n_dm2_like_cell_lines = 0`** — all 13 CCLE thyroid lines classified DM1-like → no DM2 comparison group → "DM1-selective" definition vacuous. Top drugs are GABA-receptor modulators which are biology-implausible for thyroid actionability — likely false positives from BRAF-like cell-line GABRB2 expression Tanimoto matching. Cannot anchor Genome Med Figure 7 with this. |
| **FIX2** | ⚠ Half-done | CNV: `cnv_nonzero_frac_dm2_minus_dm1 = 0.12` (real). Aneuploidy violin figure exists. | **Thorsson immune subtype: `thorsson_status = MISSING`**. Friend did not fall back to cBioPortal `thca_tcga_pan_can_atlas_2018` clinical attributes (which contains the immune subtype). AMP1 + AMP3 both depend on this. |
| **FIX3** | ✅ Real | 4 endpoints p<0.05 (best = `rai_score_recalc`); Cox multivariable HR=0.82. KM panel + forest + endpoint volcano figures exist. | Cox HR 0.82 is not significant in either direction — should be framed honestly as "cluster does not independently predict survival after adjusting for stage/age" rather than as an endpoint hit. |
| **FIX4** | ✅ Real | proxy/proper agreement = 1.0 (5 pathways overlap, all sign-concordant) | None significant. |
| **FIX5** | ⚠ Half-done | 7 patients, 8 immune cell types, single-cell DM signature exists. | **`mutation_status_values = ['NA']`** — GEO GSE184362 patient-level metadata fetch failed. Per-patient correlation between dominance vs BRAF status / age / histology is therefore impossible. Can only show dominance heterogeneity. |

### 3.4 Friend's pre-v17p35 deliverables (v17p3 paper outline + reviewer response)
- `reports/v17p3/v17p3_paper_outline_FINAL.md` — **22 lines, stub.** 5 title candidates + 3 abstract angles + 7 figure list. **No actual paragraphs**: no Results, no Discussion, no Methods, no Limitations.
- `reports/v17p3/v17p3_response_to_reviewers_v2.md` — **9 lines, routing table.** 7 reviewer concerns each pointing to a single figure file. Not a defense doc.
- `reports/v17p3/v17p3_external_review_dump_v3.md` — 43 lines (light).
- `reports/v17p3/RAW_DUMP_for_external_review.md` — 2,072 lines (substantial; this is real).

### 3.5 Friend's pre-v17p3 substantive analyses (these are real)
v17p3 tables — 67 files, all with substance:
- A1 scRNA — 7 patients × dominance ratio 0.76–1.00, 6 mixed.
- A2 BRAF/RAS orthogonality — full cohort DM score table.
- A3 drug response — **EMPTY** (1-byte tumor_predicted_response.tsv; friend papered over this in paper outline by citing the empty file as Figure 7 evidence).
- A4 pancancer transfer — 5-cancer applicability table.
- A5 immune evasion — PD-L1/IFN-γ/HLA panel real.
- A6 PDTC/RAI/TF circuit — RAI score 8.34–11.65 (PDTC) vs 2.89–7.62 (ATC), AUC 0.012 = **perfect separation, opposite direction** (this is the AMP2 hook).
- F1 external 5-cohort — 2/5 cohorts achieve DIA_AUC > 0.85.
- F2 proper GSEA — 51 hallmark pathways with NES + FDR (real).
- F3 alternative endpoints — 4 endpoints p<0.05.

---

## 4 · Critical holes (must fix before AMP)

### 4.1 FIX1 threshold — DM2-like = 0 cell lines

```
FIX1_celline_dm_scores.tsv (13 lines):
  All 13 lines: dm_like = "DM1_like"
  No DM2-like cell lines → "DM1-selective drug" undefined
```

**Hypothesised root cause**: DM1 score uses MAPK-active gene set with high baseline in established thyroid cell lines (which all derive from aggressive variants). Threshold needs to be relative (top-50 percentile vs bottom-50) rather than absolute, OR DM2-like reference signature needs different weighting.

**Fix recommendation**: Either (a) re-rank cell lines by DM1/DM2 ratio and split at median (forcing 6 vs 7), or (b) add a BRAF/RAS-context negative control set (FTC-133, TT, normal-thyroid-line if available) and recompute z-scores within thyroid context only.

### 4.2 FIX2 Thorsson — backup needed

Friend's spec listed three sources (Thorsson 2018 Immunity, iAtlas GitHub, supplementary Table S1). All apparently failed.

**Backup**: cBioPortal study `thca_tcga_pan_can_atlas_2018` has `Subtype_Immune_Model_Based` clinical attribute. URL: `https://www.cbioportal.org/api/studies/thca_tcga_pan_can_atlas_2018/clinical-data?clinicalDataType=PATIENT&attributeId=SUBTYPE`. This gives Thorsson C1–C6 immune subtype per patient.

### 4.3 FIX5 GSE184362 metadata — backup needed

Friend's GEO metadata fetch for patient-level BRAF status etc. returned all NA.

**Backup**: GEO Series Matrix file at `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184362/matrix/GSE184362_series_matrix.txt.gz` has per-sample characteristics. If still empty, try ArrayExpress mirror or fall back to publication supplementary.

---

## 5 · Friend's proposed AMP order vs my recommendation

### 5.1 Friend's plan

```
AMP2 → AMP3 → AMP1 (only 3 of 5 AMP listed)
```

Friend defensible logic: AMP2 = cleanest win (A6 data complete); AMP3 = data ready (F2 + A1 + A5); AMP1 last because FIX2 Thorsson hole weakens outlier characterisation.

**AMP4 and AMP5 not mentioned in friend's plan — this is the gap.**

### 5.2 My recommended order

```
0a. (preflight) FIX1 threshold rerun → 6 DM2-like cell lines forced
0b. (preflight) FIX2 Thorsson cBioPortal backup → Thorsson subtype filled
1.  AMP2 (RAI gradient trajectory reframe)         ← Figure 3 hook
2.  AMP4 (TF circuit + RAI 8-gene decision tool)   ← Figure 6 npj 결정타 (friend's plan missing!)
3.  AMP3 (Hot/Cold immune integrated landscape)    ← Figure 5
4.  AMP1 (Outlier 2 BRAF/DM2 + 15 RAS/DM1 deep)    ← Figure 4 paper core
5.  AMP5 (Pan-cancer transfer)                      ← Supplementary, last
```

**Why AMP4 cannot be skipped**: User's prompt explicitly says "AMP4: Decision tool (Figure 6, npj 결정타)". Skipping it drops npj P(accept) from ~70% to ~55% per the user's own venue model. Friend's plan would ship without the npj hook.

### 5.3 Order-dependency analysis

| AMP | Depends on | Blocked by current holes? |
|---|---|---|
| AMP2 | A6 PDTC data (already complete in v17p3) | No |
| AMP4 | A6 TF circuit (complete) + 8-gene panel definition | No |
| AMP3 | F2 GSEA (complete) + A1 scRNA (complete) + A5 immune (complete) + FIX2 Thorsson | **Yes — needs FIX2 backup first** |
| AMP1 | A2 outlier list (complete) + FIX2 Thorsson + FIX5 patient metadata | **Yes — needs FIX2 backup first** |
| AMP5 | A4 pancancer transfer (complete) | No |

So actual blocking dependency: **FIX2 Thorsson backup must precede AMP3 + AMP1.** AMP2, AMP4, AMP5 are independent.

---

## 6 · Tier 3/4 — the real ship gate

### 6.1 SYNTH-1 (master figures) — 0%

Friend has 15 individual fix-level figures. Spec requires **7 main paper figures + 15 supplementary**, each composed from multiple panels. None composed. This is a 60–90 minute scripting job (Plotly subplot composition) but absolutely required for submission.

### 6.2 SYNTH-2 (Korean dashboard) — 0%

`reports/v17p35/index.html` is a 56-line placeholder. Spec requires Pretendard / Noto Sans KR, 7 figure gallery + 15 supp, interactive table view, methods in Korean, decision tool live, Lighthouse score 90+. This is a 60–120 minute build.

### 6.3 DRAFT-1 (paper outline FINAL) — 10%

Currently 22-line stub at `reports/v17p3/v17p3_paper_outline_FINAL.md`. Spec requires:
- Title × 5 (have)
- Abstract × 3 versions (have angles, no actual abstracts)
- Significance statement
- Figure 1–7 caption Korean + English
- Results 7 paragraph skeleton in English
- Discussion 5 paragraphs
- Methods (all v17 phases integrated)
- Limitations (honest)
- References (extracted from v14_priorart + new search)

This is the largest remaining piece — probably 4–6 hours of writing if all analysis is final.

### 6.4 DRAFT-3 (reviewer defense) — 5%

Currently 9-line routing table. Spec requires 15–20 reviewer attacks × rebuttal paragraphs.

**The user (Seungho) has explicitly flagged that the existing v17p3 outline + response are stubs the friend papered over.** Friend's "_FINAL" naming is misleading.

---

## 7 · Honest venue probability (post-Tier-1, before any fix)

| Venue | IF | P(accept) — current state | P(accept) — after my recommended sprint |
|---|---|---:|---:|
| **npj Precision Oncology** | ~6.3 | ~50% — A3/FIX1 hollow, AMP4 missing, paper outline stub | ~70% — if AMP4 + FIX1 threshold + DRAFT-1 real + SYNTH-1 done |
| **Genome Medicine** | ~9.0 | ~25% — drug actionability weak | ~40% — if FIX1 produces credible DM-selective drug list AND AMP3 hot/cold is figure-grade |
| **Bioinformatics** (methods paper fallback) | ~5.4 | ~50% | ~60% with v5.2 ComBat clean |
| **JCO Precision Oncology** (clinical letter, separate submission) | ~6.3 | ~45% (already independently set up via v13 standalone draft) | ~50% with AMP4 decision tool added as Figure 1B |

If we ship **without** AMP4 and **without** fixing FIX1 threshold: npj drops to ~50%, Genome Med to ~25%. If user submits the friend's current state today as-is, **reviewer will open `A3_top_drugs_dm1_selective.tsv` (cited as Figure 7 evidence) and see a header-only file. Instant reject.**

---

## 8 · Concrete "what to do next" for web Claude to validate

### 8.1 If web Claude agrees with my plan, the next sprint message to the codex friend should be:

```
Friend, AMP plan adjustment:

PRE-AMP: 두 hole 먼저 막자.
  preflight-1: FIX1 DM1/DM2 threshold rerun — median split or relative ranking
               → forced 6 DM2-like cell lines so Mann-Whitney has comparison group
  preflight-2: FIX2 Thorsson backup via cBioPortal
               URL: https://www.cbioportal.org/api/studies/thca_tcga_pan_can_atlas_2018/
                    clinical-data?clinicalDataType=PATIENT&attributeId=SUBTYPE
               → fill Thorsson C1-C6 per TCGA-THCA sample

AMP order:
  1. AMP2 trajectory reframe (PDTC RAI 8.34-11.65 vs ATC 2.89-7.62, AUC 0.012 → reframe)
  2. AMP4 RAI decision tool (8 gene LogReg + ROC + PDTC validation)  ← npj 결정타, 친구 plan에 빠짐
  3. AMP3 immune Hot/Cold synthesis (F2 + A1 + A5 + FIX2 Thorsson)
  4. AMP1 outlier deep (2 BRAF/DM2-like + 15 RAS/DM1-like)
  5. AMP5 pancancer (supplementary)

그 다음 무조건 SYNTH-1 (7 main + 15 supp figure 합성) + DRAFT-1 (real paper outline 작성) + DRAFT-3 (15-20 reviewer defense). 
Tier 3/4 안 채우면 submit-ready 아님.
```

### 8.2 Open questions for web Claude (please answer)

1. **Is my AMP order correct?** Specifically: should AMP4 (RAI decision tool) really be #2, or can it move later?
2. **Is FIX1 threshold fix the right approach?** Or should we accept "all cell lines DM1-like" and reframe FIX1 as "BRAF-like cell-line drug response" (no DM split at cell-line level)? This would lose Figure 7 but gain methodological cleanliness.
3. **Is the npj 70% target realistic given the 5 prior-art papers eroding novelty?** Specifically Liu 2018 + Kalfert 2024 + Nieto-Jiménez 2023 already establish BRAF↔TROP2 + SG-for-thyroid. We're left with DM1/DM2 axis + SN-38 polarity + decision tool. Is that enough?
4. **Should we skip Genome Medicine and go straight to npj + JCO-PO split?** Genome Med raises drug-actionability bar; JCO-PO has lower bar. Two submissions might be lower combined risk.
5. **Should v5.2 ComBat LODO fix be done before or after this sprint?** It's flagged as Bioinformatics fallback blocker. Currently un-touched. Doing it now means npj sprint pauses; doing it after means Bioinformatics fallback stays open.

### 8.3 Resource constraints to flag for web Claude

- **No GPU available.** scRNA, deep learning models will be slow.
- **Codex friend works asynchronously** — they may pick up tasks I queue, or they may not. Coordination is via shared filesystem and the user routing decisions between us.
- **AlphaFold3 access blocked** (gated). AF2 only.
- **PRISM 19Q4 only** — newer releases may have CBD/dronabinol but a scheduled agent (2026-05-23) is checking.
- **User has explicitly set scope boundaries** — no wet-lab; no new experiments; in-silico only this sprint.
- **Time budget for this sprint as user wrote it**: 1–2 hours of compute, but realistically with stops/retries = ~half a working day to fully ship Tier 2/3/4.

---

## 9 · File inventory (for reference)

### v17p35 outputs as of 2026-04-27 09:00 KST

```
results/v17p35/
├── orchestrator.log                          # 169 bytes initial → expanded after FIX1-5
├── tables/
│   ├── FIX1_celline_dm_scores.tsv            # 14 lines (13 cell lines, all DM1)
│   ├── FIX1_top_drugs_dm1_selective.tsv      # 24 lines
│   ├── FIX1_top_drugs_dm2_selective.tsv      # 24 lines
│   ├── FIX1_tumor_predicted_response.tsv
│   ├── FIX1_moa_enrichment.tsv
│   ├── FIX1_summary.json                     # n_drugs_fdr<0.1=2, n_dm2_like=0
│   ├── FIX2_aneuploidy_per_sample.tsv
│   ├── FIX2_thorsson_subtype.tsv             # MISSING data
│   ├── FIX2_genomic_summary_per_dm.tsv
│   ├── FIX2_summary.json                     # thorsson_status=MISSING
│   ├── FIX3_cdr_clinical_full.tsv
│   ├── FIX3_alternative_endpoints_full.tsv
│   ├── FIX3_cox_multivariable.tsv
│   ├── FIX3_km_logrank_per_endpoint.tsv
│   ├── FIX3_summary.json                     # cox_HR=0.82, n_endpoints<0.05=4
│   ├── FIX4_proxy_vs_proper_full.tsv
│   ├── FIX4_agreement_summary.tsv
│   ├── FIX4_summary.json                     # agreement=1.0
│   ├── FIX5_per_patient_full.tsv
│   ├── FIX5_immune_cell_type_breakdown.tsv
│   ├── FIX5_per_cell_dm_signature.tsv
│   └── FIX5_summary.json                     # mutation_status=['NA']
└── figs/
    ├── FIX1_drug_volcano.html
    ├── FIX1_moa_pathway_radar.html
    ├── FIX1_top_drugs_heatmap.html
    ├── FIX1_tumor_prediction_dist.html
    ├── FIX2_aneuploidy_violin.html
    ├── FIX2_cnv_genome_view.html
    ├── FIX2_thorsson_distribution.html
    ├── FIX3_cox_forest.html
    ├── FIX3_endpoint_volcano.html
    ├── FIX3_km_panel.html
    ├── FIX4_proxy_proper_scatter.html
    ├── FIX4_top_pathway_concordance.html
    ├── FIX5_immune_composition_stacked.html
    ├── FIX5_per_patient_dm_landscape.html
    └── FIX5_single_cell_dm_score_umap.html

reports/v17p35/
└── index.html                                # 56-line placeholder, NOT the SYNTH-2 deliverable

notebooks_or_scripts/
├── v17p35_common.py                          # 15 lines
├── v17p35_orchestrator.py                    # 41 lines, TIER1=[FIX4] only
├── v17p35_FIX1_drug_response.py
├── v17p35_FIX2_genomic_immune.py
├── v17p35_FIX3_clinical_full.py
├── v17p35_FIX4_gsea_compare.py
├── v17p35_FIX5_scrna_full.py
└── (no AMP, SYNTH, DRAFT scripts exist)
```

### Pre-existing v17p3 substance (real, anchors the paper)

```
results/v17p3/tables/  (67 files)
  A1 scRNA — 7 patient dominance, real
  A2 BRAF/RAS orthogonality — full cohort DM score, real
  A3 drug response — EMPTY (header-only, papered over in outline)
  A4 pancancer — 5 cancer transfer
  A5 immune evasion — real
  A6 PDTC/RAI/TF — real
  F1 external 5-cohort — 2/5 robust
  F2 proper GSEA — 51 hallmark pathways
  F3 alternative endpoints — 4 p<0.05

reports/v17p3/
  v17p3_paper_outline_FINAL.md          22-line stub
  v17p3_response_to_reviewers_v2.md     9-line routing table
  v17p3_external_review_dump_v3.md      43 lines (light)
  RAW_DUMP_for_external_review.md       2,072 lines (substantial)
```

### Pre-existing v14 prior-art context (anchors honest framing)

```
reports/v14_novelty_verdict.md            (90 lines, written 2026-04-25)
reports/v14_novelty_verdict_addendum.md   (74 lines, written 2026-04-25)
results/v14_ccle/v14_priorart/
  all_queries_results.tsv (62 PubMed records, 8 queries)
  high_relevance_papers.md (10 HIGH-relevance papers)
  trial_status.md (NCT06235216, NCT07521670)
```

---

## 10 · One-paragraph executive summary

The codex friend completed **5 of 14+ tasks** in the v17p35 sprint — all Tier 1 (FIX) — and stopped at 08:50 KST. Three of those tasks have **substantive output** (FIX3, FIX4, FIX1 partially), two have **fetch-failure stubs** (FIX2 Thorsson missing, FIX5 patient metadata NA), and **all 11 Tier 2/3/4 tasks are unstarted** (5 AMP, 2 SYNTH, 4 DRAFT). Friend now proposes AMP2 → AMP3 → AMP1; this is **75% correct but skips AMP4** which the user's own prompt named the "npj 결정타." My recommended sprint inserts two preflight fixes (FIX1 threshold + FIX2 Thorsson cBioPortal backup), reorders to AMP2 → AMP4 → AMP3 → AMP1 → AMP5, and then mandates SYNTH-1 + DRAFT-1 + DRAFT-3 to actually reach submit-ready. Without those Tier 3/4 deliverables the v17p3 paper outline stays a 22-line stub and the v17p3 reviewer response stays a 9-line routing table, both of which the friend wrote but did not actually finish. Honest venue probability post-full-sprint: **npj ~70%, Genome Med ~40%**, conditional on FIX1 threshold producing a credible DM-selective drug list and DRAFT-1 having real paragraphs. Submitting the friend's current state today would fail at peer-review-stage one when a reviewer opens `A3_top_drugs_dm1_selective.tsv` cited as Figure 7 evidence and finds a header-only file.

---

_End of situation report. Filename: `reports/v17p35/situation_for_web_claude_20260427.md`. Self-contained. Forward to web Claude as a single attachment for verdict._
