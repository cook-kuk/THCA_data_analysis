# Paper 1 v8 outline — scope split audit (2026-05-04)

**Author:** Seungho Cook (audit by Claude, read-only)
**Date:** 2026-05-04 (post image-DM1 NO-GO, post Paper 2 HT-isolated reframe)
**Mode:** read-only scope audit. **No file modified. No analysis run. No download.**
**Voice-risk:** ZERO. Voice-protected sections (Hook / Aim / Disc §3.1 / §3.4 Limitations / Cover ¶1 / Q9) untouched.

**Authority sources:**
- `project/reports/2026_05_03_manuscript_v8_OUTLINE.md` (object of audit)
- `project/reports/2026_05_04_paper1_dm1_full_molecular_only_lock.md` (Paper 1 manuscript-safe envelope, §2 16 claims)
- `project/reports/2026_05_04_paper1_dm1_full_conflict_audit.md` (claim-by-claim source verification)
- `project/reports/2026_05_04_paper2_post_image_dm1_nogo_status.md` (Paper 2B HT-only scope)
- `PAPERS_5_OVERVIEW_2026_05_04.md` (canonical 5-paper map)
- Memory: `paper_numbering_2026_05_04`, `v18_paper2_HT_isolated`, `v19_paper3_ici_track_a`, `v19_paper4_GD_backlog`, `v17_marathon_mode_post_pillar1`, `v17_sprint_vs_marathon_violation`

---

## 1. Executive verdict

**v8 outline must be SPLIT before any further manuscript writing.**

v8 was authored **2026-05-03**, *before* (a) the 2026-05-04 Yu-meeting decision that reframed Paper 2 to HT-only with Pillar I baseline = AFND South Korea (and Chu 2018 GD → Paper 4), and (b) the 2026-05-04 closure battery that retired the image-DM1 angle. As a result, the current v8 fuses Paper 1, Paper 2B, and Paper 4 content under a single "5-Pillar" frame.

Three independent audits now point the same direction:

| Audit | Finding |
|---|---|
| `paper1_dm1_full_molecular_only_lock` §2 (16 claims envelope) | Paper 1's manuscript-safe claim set is molecular / spatial / TF / Cox / TROP2-tumor-level; **the HLA-forest, GSE286332-DEG, HLA-II-mediation, BCR/TLS items are not in this envelope.** |
| `v18_paper2_HT_isolated` (memory anchor) | Paper 2 = HT-overlap PTC molecular-only; Pillar I v2 = Korean PTC vs Korean baseline; GD/Chu 2018 → Paper 4 |
| `paper2_post_image_dm1_nogo_status` §2 | Paper 2B "active" pillar list explicitly includes Pillar I v2 forest, GSE286332, TCGA HM-like, BCR/TLS, DM1 sub-B = NBNR — **all items currently inside v8 R1/R2/R5a/R5b/R5c.** |

In other words, ~60% of v8 Results body and ~70% of v8 figures are **not Paper 1 content** under the current canonical numbering.

**Recommendation:** Author should split v8 outline into a Paper-1-only outline (molecular axis, ship target 6/13 bioRxiv) and a Paper-2B outline (HT-overlap mechanism, Pillar II/III to be defined by Yu input). Paper 4 reserve items go to backlog. **Do not start writing prose on v8 as-is** — the prose would have to be torn back out 2 weeks later when reviewer asks "why is HLA in a molecular-axis paper?"

---

## 2. Why the current v8 outline is a problem (concretely)

### 2.1 Reframe sequence v8 missed

| Date | Event | Effect |
|---|---|---|
| **2026-04-30** Yu meeting | Marathon mode declared (5/4–6/13); Pillar 1 strong = analysis end | v8 was supposed to be Paper 1 manuscript scaffold — but was assembled before Paper 2 split decision |
| **2026-05-03** | v8 outline written | snapshot of "5-Pillar Cancer Paper" frame at this time |
| **2026-05-04 morning** Yu meeting | Paper 2 reframed → HT-only; Pillar I v2 = Korean PTC vs **Korean baseline** (AFND), not Chu 2018 GD | v8 R1 / Pillar 1 forest no longer Paper 1 substrate |
| **2026-05-04** | Phase A NO-GO + closure battery NO-GO | image-DM1 dropped (already absent from v8 prose, but v8 still referenced "5-Pillar" frame that included autoimmune mechanism in Paper 1) |
| **2026-05-04** | Canonical numbering memo | Paper 1 = DM1 molecular dark matter; Paper 2B = HT-overlap; Paper 3 frozen; Paper 4 = Korean GD HLA |
| **2026-05-04** | `paper1_dm1_full_molecular_only_lock` | Paper 1 envelope locked to 16 molecular claims; HLA-forest / HT-overlap / BCR-TLS not in it |

→ v8 is **two days stale** on scope.

### 2.2 Specific contamination pattern

| v8 element | Per-2026-05-03 frame | Per-2026-05-04 canonical |
|---|---|---|
| Title candidate (1) "autoimmune-driven dedifferentiation" | OK in 5-pillar | hybrid Paper 1 + Paper 2B |
| Title candidate (2) "autoimmune-PTC overlap" | OK | Paper 2B |
| Title candidate (3) "Pan-Asian transcriptional and HLA framework" | OK | hybrid Paper 1 + Paper 4 |
| Abstract Methods "K2 + Lee 2024 + GSE286332 + Chu 2018" | OK | mixes 4 papers' substrates |
| Abstract Results (1) Pan-Asian HLA / DPB1\*05:01 / Chu GD | OK | **Paper 4** |
| Abstract Results (2) GSE286332 PTC+HT | OK | **Paper 2B** |
| Abstract Results (3) BRAF/RAS mRNA neutrality | OK | **Paper 1** |
| Abstract Results (4) Pan-genome ARI | OK | **Paper 1** |
| Abstract Results (5) HLA-II mediation / TLS / BCR / Hashimoto-like | OK | **Paper 2B** |
| Intro ¶3 "Korean / Asian-specific HLA architecture" + DPB1\*05:01 + Chu 2018 | OK | **Paper 4** intro |
| Methods M3 (PyDESeq2 + GSEA on GSE286332) | OK | **Paper 2B** |
| Methods M5 (TCGA Hashimoto-like signature transfer) | OK | **Paper 2B** |
| Methods M6 (mediation analysis) | OK | **Paper 2B** |
| Methods M7 (BCR repertoire + TLS) | OK | **Paper 2B** |
| Methods M9 (Pan-Asian HLA forest, **vs Chu GD**) | OK | **Paper 4** (per v18_paper2_HT_isolated decision) |
| Results R1 Pillar 1 (Korean PTC vs Chu GD forest) | OK | **Paper 4** |
| Results R2 Pillar 2 (GSE286332 PTC+HT DEG) | OK | **Paper 2B** |
| Results R3 Pillar 3 (driver mRNA neutrality) | OK | **Paper 1** |
| Results R4 Pillar 4 (pan-genome ARI) | OK | **Paper 1** |
| Results R5a Pillar 5a (HLA-II mediation) | OK | **Paper 2B** |
| Results R5b Pillar 5b (TCGA Hashimoto-like generalization) | OK | **Paper 2B** |
| Results R5c Pillar 5c (BCR clonal + TLS) | OK | **Paper 2B** |
| Results R5d (DM1 sub-B = NBNR cluster, mutation breakdown) | OK | **Paper 1 main**, with optional cross-paper bridge sentence to Paper 2B |
| Discussion D2 (ATA 2015 + DM1/DM2) | OK | **Paper 1** |
| Discussion D3 (Pan-Asian autoimmune-thyroid continuum) | OK | **Paper 4** |
| Figure F1 (Pan-Asian HLA forest) | OK | **Paper 4** |
| Figure F2 (GSE286332 multi-panel) | OK | **Paper 2B** |
| Figure F3 (driver mRNA neutrality) | OK | **Paper 1** |
| Figure F4 (pan-genome ARI) | OK | **Paper 1** |
| Figure F5 (autoimmune mechanism 4-panel) | OK | **Paper 2B** |
| Figure F6 (TCGA Hashimoto-like distribution) | OK | **Paper 2B** |
| Figure F7 (DM1 sub-A vs sub-B + K2 NBNR) | OK | **Paper 1 main**, with cross-paper bridge cell |
| Suppl SF2 (Mediation Baron-Kenny) | OK | **Paper 2B** |
| Suppl SF4 (GSE286332 BCR repertoire) | OK | **Paper 2B** |
| Suppl SF6 (4-scenario HLA forest sensitivity) | OK | **Paper 4** |

### 2.3 Paper 1's actual molecular envelope is mostly absent from v8

Per `paper1_dm1_full_molecular_only_lock` §2, the 16 manuscript-safe Paper 1 claims include items that **are not in v8 R1–R5 at all**:

| Claim (lock §2) | Status in v8 |
|---|---|
| #1 TCGA bulk DM1 vs THYROID_NONOVERLAP r=−0.885 (n=561) | absent from v8 R |
| #2 TCGA cluster anchor (DM1=0.00, DM2=−0.66, MW p≪0.001) | only in §1 layer-9 of DM1_FULL bundle, **absent from v8** |
| #3 PFI HR=2.04 dichotomized + DFI HR=1.41 multivariate | absent from v8 |
| #4 TF mechanism — FOXE1/NKX2-1↓ + STAT3/FOSL1↑ + DNMT1↑ | absent from v8 |
| #5 TROP2 +4.33 FDR=1.2×10⁻³² + tumor-level | absent from v8 |
| #6 Tumor vs normal TROP2/DM1 (n=517) | absent from v8 |
| #7 GSEA Hallmark IL6_JAK_STAT3 / OXPHOS down | absent from v8 |
| #8 TF coordination network (FOXE1-PAX8 etc.) | absent from v8 |
| #9 Bootstrap 95% CI 12-slide | absent from v8 |
| #10 Drop-one-out 8-gene robustness | absent from v8 |
| #11 Moran's I (per slide) 26/28 | absent from v8 |
| #12 Bivariate Moran (DM1×NONOVERLAP) | absent from v8 |
| #13 Margin-distance gradient | absent from v8 |
| #14 Pan-cancer THCA outlier (33 types) | absent from v8 |
| #15 N1 dark-matter NS (HR=1.20) | partially — v8 §2.N1 limitation only |
| #16 Harmony 28-slide UMAP | absent from v8 |

→ v8 R1–R5 contains **almost zero** of Paper 1's actual manuscript-safe claim envelope. It is currently a **Paper 2B + Paper 4 scaffold** with three Paper 1 results inserted (R3 driver mRNA, R4 pan-genome ARI, R5d sub-B mutation).

---

## 3. Section-by-section classification

### 3.1 Working title candidates

| # | Title | Verdict |
|---|---|---|
| 1 | "...autoimmune-driven dedifferentiation, independent of canonical driver mutations" | **Paper 2B** (autoimmune-driven dedifferentiation is HT mechanism) |
| 2 | "DM1/DM2 — eight-gene transcriptional axis reveals autoimmune-PTC overlap as a distinct molecular subtype" | **Paper 2B-leaning hybrid** |
| 3 | "Beyond BRAF/RAS/TERT: Pan-Asian transcriptional and HLA framework" | **Paper 1 + Paper 4 hybrid** |

→ **None of the three is a Paper-1-only title.** Paper 1 needs a new title that anchors on (a) molecular DM1/DM2 axis, (b) recurrence-prognostic, (c) driver-independent. Voice-protected — Author keyboard. Suggested seed angles (**not draft**, just framing seeds for the author):
- "RAI-lineage transcriptional axis defines recurrence-prognostic dark matter beyond BRAF/RAS/TERT in PTC"
- "An eight-gene differentiation axis stratifies PTC recurrence orthogonally to driver mutations"
- "DM1/DM2 — a driver-independent molecular subtype of papillary thyroid carcinoma with TF-network and tumor-level therapeutic vulnerability"

### 3.2 Abstract

| Block | Current content | Per-canonical paper |
|---|---|---|
| Background | DTC recurrence + ATA gap | **Paper 1** OK |
| Methods (TCGA + Lee 2024) | molecular cohort | **Paper 1** OK |
| Methods (GSE286332 + Chu 2018) | autoimmune cohorts | **Paper 2B** + **Paper 4** |
| Result (1) Pan-Asian HLA / Chu GD | DPB1\*05:01 forest | **Paper 4** |
| Result (2) GSE286332 PTC+HT DEG | HT-overlap mechanism | **Paper 2B** |
| Result (3) Driver mRNA neutrality | BRAF/RAS transcript | **Paper 1** |
| Result (4) Pan-genome ARI | cluster robustness | **Paper 1** |
| Result (5) Autoimmune-PTC mechanism | HLA-II mediation, BCR/TLS, Hashimoto-like, sub-B NBNR | **Paper 2B** (sub-B mutation breakdown is Paper 1) |
| Conclusion | hybrid framing | needs Paper-1-only re-write |

→ Abstract needs **complete re-draft** for Paper 1. Paper 1 should drop Result (1), (2), (5a/b/c) and add the molecular-envelope claims (TCGA cross-validation, Cox, TF mechanism, TROP2 tumor-level). **Voice-protected when written.**

### 3.3 Introduction ¶1–¶4

| Block | Verdict |
|---|---|
| ¶1 Hook (epi + ATA gap + Landa 2016) | **Paper 1** OK; voice-protected |
| ¶2 TCGA / Yoo / Krishnamoorthy-Landa / Pan-Ge | **Paper 1** OK |
| ¶3 Korean/Asian HLA + DPB1\*05:01 + Chu 2018 + cookHLA | **Paper 4** (move) — Paper 1 does not need this |
| ¶4 Aim + 5-pillar overview | **rewrite** — Paper 1 needs new aim that anchors on molecular DM1/DM2 axis only; voice-protected |

### 3.4 Methods M1–M11

| M | Topic | Verdict |
|---|---|---|
| M1 Cohort | TCGA/Lee Paper 1; GSE286332 Paper 2B; Chu 2018 Paper 4; K2 paper-flexible (Paper 1 if molecular-only, Paper 2B if HT/autoimmune emphasis) | **split** |
| M2 TIERA67 candidate pool | **Paper 1** |
| M3 PyDESeq2 + GSEA on GSE286332 | **Paper 2B** |
| M4 DM1/DM2 classifier (centered profile) | **Paper 1** |
| M5 TCGA Hashimoto-like signature transfer | **Paper 2B** |
| M6 Mediation (Baron-Kenny) | **Paper 2B** |
| M7 BCR repertoire + TLS | **Paper 2B** |
| M8 DM1 sub-cluster (KMeans k=2) | **Paper 1 main**; cross-paper bridge OK |
| M9 Pan-Asian HLA forest meta (vs Chu GD) | **Paper 4** (Paper 2B Pillar I v2 = Korean baseline, not Chu) |
| M10 Pan-genome MAD selection | **Paper 1** |
| M11 Statistics + reproducibility | shared (all papers) |

Paper 1 Methods after split = M1 (subset: TCGA + Lee + K2) + M2 + M4 + M8 + M10 + M11 + (NEW) M-spatial (Moran's I + bivariate Moran + margin distance) + (NEW) M-TF-inference (ULM CollecTRI) + (NEW) M-Cox (PFI/DFI multivariate).

### 3.5 Results R1–R5

| R | Topic | Verdict |
|---|---|---|
| R1 Pillar 1 Korean PTC vs Chu GD forest | **Paper 4** |
| R2 Pillar 2 GSE286332 PTC+HT DEG / GSEA / 8-gene / HLA-II | **Paper 2B** |
| R3 Pillar 3 Driver mRNA neutrality | **Paper 1** |
| R4 Pillar 4 Pan-genome cluster robustness | **Paper 1** |
| R5a Pillar 5a HLA-II mediation (Baron-Kenny) | **Paper 2B** |
| R5b Pillar 5b TCGA Hashimoto-like generalization | **Paper 2B** |
| R5c Pillar 5c BCR clonal + TLS | **Paper 2B** |
| R5d Pillar 5d DM1 sub-B = NBNR (mutation breakdown) | **Paper 1 main** + cross-paper bridge sentence to Paper 2B (4× HM-like rate) |

### 3.6 Discussion D1–D4

| D | Topic | Verdict |
|---|---|---|
| D1 Five-pillar synthesis | **rewrite** for Paper 1 (drop autoimmune pillars); voice-protected |
| D2 ATA 2015 + 8-gene complementarity | **Paper 1** OK; voice-protected |
| D3 Pan-Asian autoimmune-thyroid continuum | **Paper 4** (move) |
| D4 Limitations | **Paper 1** OK (but seed list needs update with 6 risks from `PAPER1_DM1_FULL_2026_05_04.md` §0); voice-protected |

### 3.7 Figures F1–F7

| F | Title | Verdict |
|---|---|---|
| F1 Pan-Asian HLA forest meta | **Paper 4** |
| F2 GSE286332 PTC+HT multi-panel | **Paper 2B** |
| F3 Driver mRNA neutrality | **Paper 1** |
| F4 Pan-genome ARI + 8-gene clinical interpretability | **Paper 1** |
| F5 Autoimmune-PTC mechanism 4-panel (P_DM1 mediation + TCGA HM-like + BCR/TLS + sub-B) | **Paper 2B** main; sub-B sub-panel can carry to Paper 1 |
| F6 TCGA Hashimoto-like distribution + cross-cohort | **Paper 2B** |
| F7 DM1 sub-A vs sub-B + K2 NBNR clinical phenotype | **Paper 1 main**; cross-paper bridge cell to Paper 2B |

### 3.8 Suppl figures + tables

| Asset | Verdict |
|---|---|
| SF1 TIERA67 pan-genome rank ladder | **Paper 1** |
| SF2 Mediation Baron-Kenny detail | **Paper 2B** |
| SF3 Korean GSE213647 replication detail | **split** — molecular replication = Paper 1 supp; HT-like rate replication = Paper 2B supp |
| SF4 GSE286332 BCR repertoire detail | **Paper 2B** |
| SF5 ATA 2015 risk tier × DM cluster cross-tab | **Paper 1** |
| SF6 4-scenario HLA forest sensitivity | **Paper 4** |
| Suppl Table 1 (per-allele HLA full) | **Paper 4** |
| Suppl Table 2 (top 200 PTC-vs-PTC+HT DEGs) | **Paper 2B** |

---

## 4. Paper 1 KEEP list (final)

**Cohorts:** TCGA-THCA n=500 (DM1 140 / DM2 360); Lee 2024 GSE213647 n=632 (Korean molecular replication); K2 PRJEB11591 n=260 (HLA imputation backbone — but allele-level analysis itself goes to Paper 2B/Paper 4). Plus GSE250521 28-slide ST + 12 external Visium for spatial validation.

**Methods (Paper 1 envelope):**
- M1 (subset) Cohort assembly: TCGA + Lee 2024 + K2 (molecular only, no HLA forest)
- M2 TIERA67 candidate pool
- M4 DM1/DM2 centered-profile classifier (TCGA AUC 0.962)
- M8 DM1 sub-cluster (KMeans k=2 → sub-A/sub-B, mutation breakdown)
- M10 Pan-genome MAD selection (ARI 0.918 vs TIERA67 0.903)
- M11 Statistics + reproducibility
- **NEW M-Cox** (PFI dichotomized HR 2.04; DFI multivariate HR 1.41)
- **NEW M-TF-inference** (ULM CollecTRI 1185 TFs, n=261 high vs low)
- **NEW M-spatial** (Moran's I + bivariate Moran + margin distance + Harmony UMAP)
- **NEW M-pan-cancer** (33 TCGA types within-cohort r)
- **NEW M-tumor-vs-normal** (n=517: 49 normal + 460 primary + 8 met)
- **NEW M-drop-one-out** (8 RAI gene panel)

**Results (proposed Paper 1 R1–R5):**
- R1 (NEW) Cohort + DM1/DM2 cluster definition + DM1_like score anchor (TCGA cluster vs score, MW p≪0.001) + cross-validation (TCGA bulk r=−0.885 + spatial replicate r=−0.99)
- R2 Driver mRNA neutrality (former v8 R3)
- R3 Pan-genome cluster robustness (former v8 R4)
- R4 (NEW) Mechanism — TF collapse + GSEA Hallmark + TF coordination + 4-step framework
- R5 Clinical outcome (Cox PFI/DFI) + tumor-level TROP2 vulnerability (n=517) + DM1 sub-B = NBNR (former v8 R5d, with cross-paper bridge sentence)

**Figures (proposed Paper 1 F1–F5):**
- F1 TCGA cluster anchor + cross-validation scatter
- F2 Driver mRNA neutrality (former v8 F3)
- F3 Pan-genome ARI + TIERA67 ranking (former v8 F4)
- F4 TF mechanism volcano + GSEA + TF coordination
- F5 KM Cox + TROP2 tumor-level + DM1 sub-A/B mutation breakdown (former v8 F7 reorganized)

**Suppl (proposed Paper 1 SF1–SF8):**
- SF1 Harmony 28-slide UMAP (ST QC backbone)
- SF2 8-gene drop-one-out robustness
- SF3 Pan-cancer THCA outlier
- SF4 Moran's I + bivariate Moran + margin distance
- SF5 Bootstrap 95% CI 12-slide sample-mean
- SF6 TIERA67 pan-genome rank ladder (former v8 SF1)
- SF7 Korean GSE213647 *molecular* replication (DM1_like score r in Korean cohort)
- SF8 N1 dark-matter underpowered (HR=1.20, honest negative) + closure-battery negative-feasibility 1-pager

**Discussion (proposed Paper 1 D1–D4):**
- D1 Synthesis (rewrite — molecular axis only); voice-protected
- D2 ATA 2015 + 8-gene complementarity (kept from v8); voice-protected
- D3 (NEW, replaces v8 D3) Cross-paper boundary statement: paper2B→cross-link, paper4→cross-link, paper3→deferred
- D4 Limitations (kept; seed = `PAPER1_DM1_FULL_2026_05_04.md` §0 6 risks + §2 N1 + Q3 negative); voice-protected

**Cross-paper bridge sentences (allowed; no full-narrative cross-over):**
- "DM1 sub-B (n=56, 53/56 mutation-negative) is enriched for a Hashimoto-like rate that is examined as the HT-overlap mechanism in a companion paper (Cook et al., in preparation)."
- "An autoimmune-PTC overlap mechanism for the DM2 sub-population is dissected in a companion paper."

---

## 5. Paper 2B MOVE list

**From v8 → into Paper 2B parking lot (HT-overlap mechanism manuscript, Pillar I/II/III):**

| v8 element | Paper 2B placement |
|---|---|
| Methods M3 PyDESeq2 + GSEA on GSE286332 | Paper 2B Methods |
| Methods M5 TCGA Hashimoto-like signature transfer | Paper 2B Methods |
| Methods M6 Mediation (Baron-Kenny + bootstrap) | Paper 2B Methods |
| Methods M7 BCR repertoire + TLS | Paper 2B Methods |
| Results R2 GSE286332 PTC vs PTC+HT (10,380 DEGs, GSEA, HLA-II d=+3.65, 8-gene d=−1.60) | Paper 2B Results (likely main result) |
| Results R5a HLA-II 140% mediation | Paper 2B Results |
| Results R5b TCGA Hashimoto-like generalization (OR up to 5×, p=6.4e-10) | Paper 2B Results |
| Results R5c BCR clonal + TLS (TLS d=+1.96; IGHV clonality; AICDA) | Paper 2B Results |
| Figure F2 GSE286332 multi-panel | Paper 2B figure |
| Figure F5 Autoimmune mechanism 4-panel | Paper 2B figure |
| Figure F6 TCGA Hashimoto-like distribution + cross-cohort | Paper 2B figure |
| Suppl SF2 Mediation Baron-Kenny detail | Paper 2B supp |
| Suppl SF4 GSE286332 BCR repertoire detail | Paper 2B supp |
| Suppl Table 2 (top 200 DEGs) | Paper 2B supp |

**Active-not-in-v8 Paper 2B items (from `paper2_post_image_dm1_nogo_status` §2, complete the parking lot):**

| Item | Notes |
|---|---|
| Pillar I v2 forest (Korean PTC pool n=874 vs **AFND South Korea baseline**) | NOT vs Chu 2018 (which is Paper 4). v2 substrate per `v18_paper2_HT_isolated`. |
| DM1 sub-B = NBNR cross-paper bridge | sub-B mutation breakdown stays in Paper 1; HM-like rate annotation can also be referenced from Paper 2B narrative |
| Pillar II / Pillar III definitions | TBD per Yu professor input; not in v8 |
| Cross-paper reciprocal Discussion line (PTC ≠ GD) | already noted in Paper 2 status memo |

→ Paper 2B is **not yet ready for outline** (Pillar II/III missing). Next Yu-meeting input required before Paper 2B outline can be drafted.

---

## 6. Paper 4 MOVE list

**From v8 → Paper 4 (Korean GD HLA / Pan-Asian backlog):**

| v8 element | Paper 4 placement |
|---|---|
| Intro ¶3 Korean/Asian-specific HLA + DPB1\*05:01 + Chu 2018 + cookHLA | Paper 4 Intro |
| Methods M9 Pan-Asian HLA forest meta (vs Chu GD) | Paper 4 Methods |
| Results R1 Pan-Asian HLA cohort n=874 (3-arm forest with Chu 2018) | Paper 4 Results |
| Discussion D3 Pan-Asian autoimmune-thyroid susceptibility continuum | Paper 4 Discussion |
| Figure F1 Pan-Asian HLA forest plot | Paper 4 figure |
| Suppl SF6 4-scenario HLA forest sensitivity | Paper 4 supp |
| Suppl Table 1 per-allele HLA full | Paper 4 supp |
| Cross-paper bridge: TSAb / TSI / hyperthyroidism / ophthalmopathy | Paper 4 (NOT in any current paper, future scope) |

→ Paper 4 is **backlog-gated**; do not draft now. Per `v19_paper4_GD_backlog`, work begins after Paper 1·2·3 settle. Move v8's HLA-forest content here as parking material.

---

## 7. DROP / INTERNAL list

**DROP from any manuscript draft (not just v8):**

| Item | Why |
|---|---|
| H&E → DM1 prediction at tile level | Closure battery NO-GO 2026-05-04, verdict D = full 폐기 |
| RunPod G1 / G2 / G3 sprints | superseded by closure battery |
| TCGA WSI external validation of DM1 | image-axis dropped + WSI download forbidden in marathon |
| Phase B / Phase C image-DM1 plans | all paths sealed |
| Image-DM1 IP / patent claim | color artifact unpatentable; ATC classifier prior art dense |
| "TROP2 spot-level colocalization with DM1-high spots" | Q3 ρ = −0.016 negative |
| "DM1-high spots are TROP2-high" | Q3 negative |
| "All risks resolved" / "모든 risk 해소" | ≥3 risks remain (N1 NS / Q3 negative / wet-lab missing); explicitly forbidden by lock §5 |
| "Cancer Cell 도전권 20-30%" / "Nat Cancer 50-65% reach" | venue probability = internal-only |

**INTERNAL-only (preserve on disk, never in submission):**

| Item | Where |
|---|---|
| `PAPER1_DM1_FULL_2026_05_04.md` §0 venue paragraph | internal status |
| `PAPER1_DM1_FULL_2026_05_04.md` §3 venue probability table | internal decision aid |
| `PAPER1_DM1_FULL_2026_05_04.md` §4 wet-lab + computational roadmap | marathon-외 영역 |
| `PAPER1_DM1_FULL_2026_05_04.md` §9 GPU Phase B package | DEPRECATED |
| `PAPER1_DM1_FULL_2026_05_04.md` §10 parallel VM overlap | DEPRECATED |
| `PAPER1_DM1_FULL_2026_05_04.md` §11 final recommend | DEPRECATED |
| `2026_05_03_marathon_W1_decision_brief.md` venue percentages | internal |

---

## 8. Proposed Paper 1 revised outline (post-split)

> Outline only. Voice-protected sections (Hook, Aim, Disc §3.1, §3.4 Limitations, Cover ¶1, Q9) **not drafted here** — author keyboard.

### 8.1 Title (author voice; seeds only)

- "RAI-lineage transcriptional axis defines recurrence-prognostic dark matter beyond BRAF/RAS/TERT in PTC"
- "An eight-gene differentiation axis stratifies PTC recurrence orthogonally to driver mutations"
- "DM1/DM2 — a driver-independent molecular subtype of PTC with TF-network and tumor-level therapeutic vulnerability"

### 8.2 Abstract skeleton (~200–250 words; author voice)

- **Background.** PTC recurrence (15–35%) lacks a transcriptional axis on top of ATA 2015 anatomic risk.
- **Methods.** TCGA-THCA n=500 dynamic-risk-stratified DM1/DM2 by 8-gene RAI panel; centered-profile classifier (CV AUC 0.962); applied to Korean Lee 2024 (GSE213647 n=632) for replication; 28-slide GSE250521 ST + 12 external Visium for spatial validation; ULM TF inference (CollecTRI 1185 TFs); Cox PFI/DFI on TCGA n=461 with BRAF/RAS/age/stage adjustment.
- **Results.** (1) DM1/DM2 cluster anchor: TCGA bulk DM1_like vs THYROID_NONOVERLAP r = −0.885 (n=561, p=5.1×10⁻¹⁸⁸); cluster vs score MW p≪0.001 (n=157). (2) Driver mRNA neutrality: BRAF transcript Cohen d = −0.044 (p=0.567); driver-anchor cluster ARI = 0. (3) Pan-genome cluster robustness: TIERA67 ARI = 0.903 ≈ unrestricted top-5000 MAD ARI = 0.918. (4) Mechanism: 4-step framework — FOXE1/NKX2-1 collapse → DNMT activation → STAT3/AP-1 activation → TROP2 re-expression; cross-checked by GSEA Hallmark (IL6_JAK_STAT3 +8.0, OXPHOS −5.5) and TF coordination (FOXE1-PAX8 r=+0.55). (5) Outcome + therapeutic vulnerability: PFI dichotomized HR = 2.04 (p=0.015); DFI multivariate HR = 1.41 (p=0.025); TROP2/TACSTD2 +4.33 FDR=1.2×10⁻³² (tumor-level); tumor 1.4× normal; DM1 sub-B (n=56, 53/56 mutation-negative) = NBNR cluster.
- **Conclusion.** A transcriptional differentiation axis stratifies PTC orthogonally to driver mutations and supports a recurrence-prognostic, mechanism-anchored, tumor-level therapeutic readout. (HT-overlap mechanism + HLA-forest pursued in companion papers.)
- **Keywords.** PTC, transcriptional axis, RAI lineage, DM1/DM2, BRAF/RAS, TERT, TF network, TROP2, recurrence.

### 8.3 Results R1–R5

| R | Title | Source claims (lock §2) |
|---|---|---|
| R1 | Cohort + DM1/DM2 cluster definition + cross-validation | #1 #2 + spatial replicate (#9 Bootstrap CI) |
| R2 | Driver mRNA neutrality (former v8 R3) | (existing source rows in `2026_05_03_results_R1_R5_prose.md`) |
| R3 | Pan-genome cluster robustness (former v8 R4) | as above |
| R4 | TF-driven 4-step mechanism + GSEA + TF coordination | #4 #7 #8 |
| R5 | Clinical outcome + tumor-level TROP2 + DM1 sub-B = NBNR | #3 #5 #6 #15 + sub-B from R5d |

### 8.4 Figures F1–F5

| F | Content | Lock §6 panels |
|---|---|---|
| F1 | TCGA DM1/DM2 cluster anchor + DM1_like vs NONOVERLAP scatter | A + B |
| F2 | Driver mRNA neutrality | (former v8 F3) |
| F3 | Pan-genome ARI + TIERA67 ranking | (former v8 F4) |
| F4 | TF mechanism volcano + GSEA + TF coordination panel | D |
| F5 | KM Cox + TROP2 tumor-level + DM1 sub-A/B mutation breakdown | C + E + F + (sub-B sub-panel from former v8 F7) |

### 8.5 Suppl

- SF1 Harmony 28-slide UMAP
- SF2 8-gene drop-one-out
- SF3 Pan-cancer THCA outlier
- SF4 Moran's I + bivariate + margin distance
- SF5 Bootstrap 12-slide CI
- SF6 TIERA67 pan-genome rank ladder
- SF7 Korean GSE213647 molecular replication
- SF8 N1 dark-matter NS + closure-battery 1-page negative-feasibility

---

## 9. Proposed Paper 2B parking lot

> Move target. **Do not draft prose now** — Pillar II/III still TBD per Yu input.

### 9.1 Items moved from v8 (already analyzed; ready as parking material)

- v8 Methods M3 / M5 / M6 / M7
- v8 Results R2 (GSE286332 PTC vs PTC+HT)
- v8 Results R5a / R5b / R5c (mediation / TCGA HM-like / BCR-TLS)
- v8 Figure F2 / F5 / F6
- v8 Suppl SF2 / SF4
- v8 Suppl Table 2

### 9.2 Pillar I v2 (Paper 2B substrate)

- Korean PTC pool n=874 vs **AFND South Korea baseline** (NOT Chu 2018; Chu → Paper 4)
- Substrate: K2 PRJEB11591 n=235 valid HLA + Lee 2024 GSE213647 n=630 valid HLA + GSE286332 PTC arm n=9
- Per `v18_paper2_HT_isolated` + `2026_05_03_paper2_pillar1_for_web_claude.md` (existing scaffold)

### 9.3 Pillar II / Pillar III (TBD)

- Definition pending Yu professor input
- Marathon-permitted scaffolding only (no analysis); manuscript drafting blocked until Pillar II/III defined

### 9.4 DM1 sub-B = NBNR cross-paper bridge

- Paper 1 carries the mutation breakdown (sub-B 53/56 = 94.6% mut-neg)
- Paper 2B carries the 4× HM-like rate annotation as cross-link
- Both papers can cite the same `n1_dark_matter_subset.tsv` source

### 9.5 Forbidden in Paper 2B (already documented in `paper2_post_image_dm1_nogo_status` §3)

H&E-DM1 / pathology AI / TCGA WSI / image-DM1 / morphology-derived DM1 — never appear in Paper 2B prose.

---

## 10. Immediate next action

### Should current v8 outline be split before further writing?

**YES — split before any further writing.**

Concretely:

1. **Freeze current v8 as historical artifact.** v8 is the 2026-05-03 5-pillar snapshot. Preserve unchanged on disk (already chmod-implied via M-status). It can be re-read for source data inventory but should NOT receive new prose.
2. **Author drafts a Paper-1-only v9 outline** (post-split) using §8 of this audit as scaffold. Scope = molecular DM1/DM2 axis only, per `paper1_dm1_full_molecular_only_lock` §2 envelope. **Voice-protected sections (Title, Hook ¶1, Aim ¶4, Disc §3.1, §3.4 Limitations, Cover ¶1) author-keyboard during W1–W2.**
3. **Author parks Paper 2B candidate material** into a separate scaffold file (e.g., `project/manuscript_p2_brief/v9_outline_paper2B_parking.md`) using §9 of this audit. **Do not draft prose** until Yu-meeting input defines Pillar II/III.
4. **Author moves Paper 4 candidate material** (HLA forest + DPB1\*05:01 + Chu 2018) into Paper 4 backlog (e.g., `project/reports/paper4_backlog/v8_extracted_HLA_material.md`). Backlog only; do not work during marathon.
5. **Update `STATUS_PAPER1_2026_05_02.md` (or successor)** to reflect v8→v9 split decision.
6. **Resume marathon** after split: voice-hook (Hook ¶1) author-keyboard, then Aim ¶4, then Disc §3.1 (Landa 2016 cite save).

### What Claude does NOT do

- ❌ Modify v8 outline prose
- ❌ Draft v9 outline (only seed angles in §8 above)
- ❌ Run any new analysis or download
- ❌ Touch Paper 3 / Paper 4 design bundles
- ❌ Author voice-protected sections (Title, Hook, Aim, Disc, Limitations, Cover, Q9)
- ❌ Use RunPod API or kill RunPod pods
- ❌ Touch H&E-DM1 / image-axis residuals

### Recommended next user-facing action

`v9-outline-split` — Author opens `project/reports/2026_05_03_manuscript_v8_OUTLINE.md` in read-only and a fresh new file `project/reports/2026_05_05_manuscript_v9_paper1_only_outline.md` for editing, then drafts a Paper-1-only outline using §8 of this audit as scaffold. Split parking files for Paper 2B and Paper 4 in parallel.

---

## 11. Marathon attestation

| Constraint | Status |
|---|---|
| Voice-protected sections | ✓ untouched |
| Paper 3 design bundle (chmod 444) | ✓ untouched |
| Paper 3 Track B | ✓ NOT started |
| Paper 4 backlog | ✓ untouched (only referenced as parking target) |
| New analysis launched | ✓ none |
| New data download | ✓ none |
| WebFetch | ✓ 0회 |
| H&E-DM1 retry | ✓ NOT attempted |
| TCGA WSI download | ✓ NOT attempted |
| RunPod API | ✓ NOT used |
| `.runpod_pods.json` modified | ✓ NOT modified |
| Manuscript prose modified | ✓ NOT modified (audit only) |
| v8 outline file edited | ✓ NOT edited (read-only audit) |

---

Audit closed. v8 outline = stale 5-pillar snapshot; needs split into Paper 1 (molecular only) + Paper 2B parking + Paper 4 backlog before any prose drafting.
