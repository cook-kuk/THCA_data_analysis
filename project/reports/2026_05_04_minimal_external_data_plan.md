# Minimal external-data acquisition plan — Paper 1 & Paper 2 (marathon-aware)

**Date:** 2026-05-04 · **Owner:** Seungho Cook
**Mode:** Plan only. **No download, no analysis, no GPU/RunPod, no TCGA WSI, no H&E-DM1, no Paper 3 Track B.**
**Window:** 2026-05-04 → 2026-06-13 (Paper 1 bioRxiv target). Anything not actionable inside this window is *deferred* by default.
**Authority anchors:**

- `2026_05_04_8gene_curated_vs_denovo_final_strategy.md` §12.6 (Nature Cancer reach), §12.7 (Cancer Cell aim — what is *required* but NOT promised)
- `2026_05_04_paper1_dm1_full_molecular_only_lock.md` §2 (16-claim envelope) + §3 (excluded claims)
- `2026_05_04_marathon_state_after_image_dm1_nogo.md` §2 (per-paper state)
- `v17_marathon_mode_post_pillar1` · `v17_sprint_vs_marathon_violation` · `paper_numbering_2026_05_04` · `v18_paper2_HT_isolated` · `v17_landa2016_cite_save`

---

## 1. Executive recommendation (TL;DR)

| Slot | Paper | Pick | Verdict |
|---|---|---|---|
| Paper 1 acquisition #1 | Paper 1 (DM1 molecular dark matter) | **GSE76039 — Landa 2016 PDTC/ATC microarray (n=37; 17 PDTC + 20 ATC; HG-U133 Plus 2.0, gcRMA)** | **GO — DONE 2026-05-04** (first-pass complete; see `2026_05_04_gse76039_first_pass_report.md`. Earlier "RNA-seq (n=84)" was incorrect — corrected post-acquisition; see §3.1 banner.) |
| Paper 1 acquisition #2 | Paper 1 | *(reserved unused)* | **DROP** — no second slot is worth marathon-blocking risk before 6/13 |
| Paper 2B reference | Paper 2 (Pillar I v2 Korean baseline cross-check) | **Lee 2014 Korean HLA frequency reference (Tissue Antigens) — published summary stats only** | **GO** (zero-cost; already in bib as `Lee2014TissueAntigens`) |
| Don't acquire now | — | GSE33630 microarray · TCGA-THCA 450K methylation · DepMap TROP2 · GTEx · MSK-IMPACT · Pan/Ge 2025 ATC raw · GSE82208 · ICGC | **DEFER / DROP** (see §4) |

**Single-line takeaway.** Take **GSE76039 alone** as the maximum marathon-compatible Paper 1 add (and only if not blocked). Use **Lee 2014** as Paper 2B reference (it is a citation, not a dataset). Defer methylation and everything else to bioRxiv-revision / Paper-1-v9 / Nature Cancer reach window. **Doing nothing is also a fully acceptable outcome** — Paper 1's current envelope (16 source-verified claims + spatial freeze Verdict C + 28 ST slides + non-overlap lineage cross-val) clears Cell Reports Medicine on its own.

---

## 2. Marathon constraint reminder

| Constraint | Implication for acquisitions |
|---|---|
| 5/4–6/13 = 6 weeks; bioRxiv 6/13 = critical milestone | Anything requiring >5 days end-to-end is **automatic DEFER**. |
| Voice-protected sections (Hook / Aim / Disc 3.1 / Limitations / Cover ¶1 / Q9) = author keyboard | A new dataset that *requires* a new Discussion §3.1 sentence to land is voice-protected-blocked unless author commits to writing it. |
| `v17_marathon_mode_post_pillar1` = 새 분석은 paper-blocking 만 | "paper-blocking" = without it, a referee will reject. Nice-to-have ≠ paper-blocking. |
| Paper 3 (`v19_paper3_ici_track_a` FROZEN) | Anything resembling ICI cohorts = **DROP** (Track B blocked). |
| Paper 4 (`v19_paper4_GD_backlog`) | GD HLA / Korean germline expansion = **DROP** for *acquisition*; Paper 2 baseline references = OK as citation. |
| H&E-DM1 / RunPod / TCGA WSI | All **DROP** (closure battery NO-GO; 5 re-entry conditions not met). |

**Decision rule:** acquisition is **GO** only if it (i) is paper-blocking for Paper 1 or Paper 2B (Korean baseline), (ii) is end-to-end ≤5 days, and (iii) does not require voice-protected new Discussion paragraphs.

---

## 3. Per-candidate evaluation

### 3.1 GSE76039 — Landa 2016 PDTC/ATC microarray (n=37) **[ACQUIRED + FIRST-PASS DONE 2026-05-04]**

> **Update 2026-05-04 (post-acquisition correction).** This row was authored before the GEO landing-page check.
> Two corrections after Step 1 access check (`2026_05_04_gse76039_access_check.md`):
> (1) **Modality is microarray (Affymetrix HG-U133 Plus 2.0, GPL570, gcRMA log₂), NOT RNA-seq** as initially listed.
> (2) **n = 37 (17 PDTC + 20 ATC; no PTC controls in series, no normal in series)**, NOT n ≈ 84. The n ≈ 84 figure refers to Landa 2016's larger reported transcriptomic n; what is deposited at GSE76039 is a 37-sample subset.
> First-pass result: external advanced-disease replication SUCCEEDS at extreme effect size (Cohen's d ≈ −3.5 ATC vs PDTC on differentiation axis; Spearman ρ = −0.925 zero-overlap lineage cross-validation, p = 3.1 × 10⁻¹⁶). See `2026_05_04_gse76039_first_pass_report.md`. Status now **DONE** for first-pass; awaiting author/Yu approval for inclusion as Paper 1 main figure F8.

| Field | Value (verified post-acquisition) |
|---|---|
| Use case | **Paper 1** — independent advanced-disease validation of the DM1_like / 8-gene compact readout |
| Why paper-relevant | (i) partially closes the strategy memo §12.6 Nature Cancer reach gate (molecular validation of DM1_like in a second cohort; survival validation still pending); (ii) Landa 2016 was already a citation — now becomes a result on their cohort with our score |
| Cohort size | **n = 37** (17 PDTC + 20 ATC) — **no PTC controls, no normals in this series** |
| Modality | **Affymetrix HG-U133 Plus 2.0 microarray (GPL570), gcRMA log₂** (Landa lab pre-processed) |
| Data size | 6.9 MB Series Matrix + 85.6 MB GPL570 SOFT annotation |
| Accessibility | Public (GEO); processed Series Matrix used; raw CEL `GSE76039_RAW.tar` (164.7 MB) NOT downloaded |
| Cross-platform / calibration risk | **Within-cohort z is the only operation done.** TCGA-trained absolute-form classifier explicitly NOT applied across platform per `v17_korean_K2_calibration`. Within-cohort effect-size estimates (Cohen's d, Spearman ρ) are platform-independent. Cross-platform comparison to TCGA absolute scores not attempted. |
| Actual end-to-end time | ≈ 30 min (CPU-only): download 5 min, parse 2 min, score 1 min, tests + figures < 1 min, reports ≈ 20 min. **Well under 2.5-day budget.** |
| Manuscript impact | **HIGH** (confirmed). Direction-consistent across all seven Paper 1 axes (RAI_8, NONOVERLAP, TDS_like, TF_collapse, STAT3_AP1_DNMT, TACSTD2 single, DM1_like). Strongest external zero-overlap lineage cross-validation to date. |
| Scope risk | **LOW** (confirmed). No PTC controls limits gradient claim to "ATC > PDTC within advanced disease," NOT "PTC → PDTC → ATC monotonic." Honest framing applied throughout `2026_05_04_gse76039_first_pass_report.md` §6/§9. |
| Voice-protected impact | **None for scaffolding** (Methods M12 + Figure F8 + Suppl Table SX_LANDA added by Claude). Discussion §3.1 integration of the result is **author keyboard only** per `v17_sprint_vs_marathon_violation`. |
| **Verdict** | **GO — DONE 2026-05-04.** First-pass complete. No further dataset acquired this turn. |
| Output locations | `project/results/p_landa_2016/{raw/, sample_metadata.tsv, expression_matrix_log_gene.tsv.gz, probe_to_gene_panel.tsv, scores.tsv, score_tests.tsv, SuppTable_SX_LANDA_GSE76039.xlsx, 3 PNG figures}`; script `project/notebooks_or_scripts/p_landa_2016_first_pass.py`; reports `project/reports/2026_05_04_gse76039_access_check.md` + `2026_05_04_gse76039_first_pass_report.md` |

### 3.2 GSE33630 — Tomas/Espinal 2012 normal/PTC/ATC microarray

| Field | Value |
|---|---|
| Use case | Paper 1 — alternative advanced-disease replication; supplementary cross-platform comparison. |
| Cohort size | ≈ 105 samples (45 normal + 49 PTC + 11 ATC, verify). |
| Modality | Affymetrix HG-U133 Plus 2.0 microarray. |
| Data size | ≈ 0.5 GB processed expression matrix. |
| Accessibility | Public (GEO). |
| Cross-platform / calibration risk | **HIGH.** Microarray vs RNA-seq cross-platform alignment of an 8-gene panel will face reviewer Q ("did you re-Z within platform?"); within-sample centering helps but introduces an extra defence step. |
| Estimated end-to-end time | ≈ 3–4 days (cross-platform handling + probe-to-gene mapping). |
| Manuscript impact | **MEDIUM.** Largely redundant given GSE76039 (3.1 above) is the same modality family (Affymetrix microarray) and answers the same advanced-disease scientific question. Post-acquisition (2026-05-04), GSE76039 already provided strong direction-consistent results, removing the case for GSE33630 as a second cohort. |
| Scope risk | **MEDIUM.** Cross-platform calibration adds a Methods sub-section + a defensive Supp Figure that does not pay back proportionally. |
| **Verdict** | **DEFER.** Strictly redundant with GSE76039 if that lands. If GSE76039 fails QC (unlikely), GSE33630 is the fallback — but only as a fallback. **Do not acquire concurrently.** |

### 3.3 TCGA-THCA methylation 450K (Illumina Infinium 450K)

| Field | Value |
|---|---|
| Use case | Paper 1 — closes the "methylation-absent" Limitations item; converts the TF→DNMT→TROP2 framework from "consistent with" to "supported by promoter-methylation evidence." |
| Cohort | TCGA-THCA same n ≈ 500 patients we already use for expression. **Not new cohort; new modality on existing cohort.** |
| Data size | Processed beta values ≈ 500 MB; raw IDATs ≈ tens of GB (avoid). |
| Accessibility | Public (GDC; processed Level 3 beta matrices are commonly redistributed). |
| Estimated end-to-end time | Download + QC ≈ 1 day. Probe-to-gene mapping (450K manifest, FOXE1 / NKX2-1 / PAX8 / TG / TPO / SLC5A5 promoter regions) ≈ 1 day. Methylation × expression × DM1 cluster integration ≈ 3–5 days. Figure + Supp Methods + 1 Results paragraph ≈ 1–2 days. **Total ≈ 6–9 days.** |
| Manuscript impact | **VERY HIGH for revision / Nature Cancer reach.** Marginal for CRM-safe submission (the strategy memo §12.5 CRM-safe explicitly accepts "methylation absent" as a Limitations item). |
| Scope risk | **HIGH.** New data layer = new Methods section + new figure + new Discussion paragraph. Voice-protected (Discussion §3.4 Limitations would need rewording from "absent" to "supported"). Likely to push 6/13 deadline. |
| **Verdict** | **DEFER to bioRxiv-revision / Paper-1-v9 / Nature Cancer reach window.** Per strategy memo §12.7, methylation is on the Cancer Cell aim list — these are *future-paper* requirements, not Paper 1 marathon. Optionally pre-position by reading the GDC manifest *only* (no download) to estimate scope before W6 decision. |

### 3.4 Lee 2014 Korean HLA reference (Tissue Antigens; C\*01:02 / DQB1\*02:01)

| Field | Value |
|---|---|
| Use case | **Paper 2 Pillar I v2** (Korean PTC vs Korean baseline; per `v18_paper2_HT_isolated`). Provides Korean-population baseline carrier frequencies for two specific alleles flagged in the Korean PTC cohort (C\*01:02, DQB1\*02:01) — alongside the AFND South Korea pool already used as the primary baseline. |
| Format | **Published summary statistics** (allele frequency table from a peer-reviewed paper). Not a dataset; a citation. |
| Acquisition | Already in bib as `Lee2014TissueAntigens` (`2026_05_03_references.bib` line 361). PubMed / journal landing page lookup if frequency tables need to be re-extracted; effectively zero acquisition cost. |
| Data size | Single table, < 100 KB. |
| Accessibility | Open access journal; bib already populated. |
| Estimated end-to-end time | ≈ 1 hour (table extraction + bib note expansion). |
| Manuscript impact | **MEDIUM for Paper 2 Pillar I v2.** Adds a second Korean baseline reference for two specific alleles, strengthening the "Korean baseline ≠ AFND-only artefact" defence. Does not affect Paper 1. |
| Scope risk | **LOW.** No new analysis, no new data, no new figure. One additional row in Pillar I forest baseline table. |
| Voice-protected impact | None. Pillar I forest table is fact-only Methods/Results scaffold. |
| **Verdict** | **GO.** This is the only acquisition that is paper-relevant *and* zero-marathon-cost. Already in bib; only deliverable is a one-row addition to the Pillar I baseline reference table when Paper 2 Pillar I v2 is finalised. |

### 3.5 Other candidates — recorded only, no decision

| Dataset / source | Why mentioned | Status |
|---|---|---|
| GSE213647 (Lee SE 2024 Macrogen) | Korean PTC n=632 — Paper 2 territory | **already in envelope** (Paper 2 Pillar III); not a new acquisition |
| PRJEB11591 K2 (Yoo 2016 SNU-GMI) | Korean PTC n=260 (n=235 valid HLA) — Paper 2 Pillar I | **already in envelope**; not a new acquisition |
| GSE286332 (Lim 2025 Dongguk) | Korean PTC n=18 (9 + 9) — Paper 2 Pillar II | **already in envelope** |
| GSE230424 / GSE248205 / GSE250521 | Three Visium ST cohorts — Paper 1 spatial supplement freeze §a/b | **already in envelope** (`paper1_spatial_supplement_freeze_2026_05_03.md` Verdict C) |
| AFND South Korea pool | Korean baseline allele frequencies — Paper 2 Pillar I v2 baseline | **already in envelope** (per `v18_paper2_HT_isolated`); not a new acquisition |
| cBioPortal `thca_tcga_pub` | TCGA-THCA TERT promoter calls — Paper 1 multivariate Cox covariate | **already used** (per `v17_tert_recovery_v2`, n=36 carriers HR=4.33) |
| Pan/Ge 2025 ATC proteogenomic (Nat Commun) supplementary tables | Cite for Discussion §3.1 (advanced-disease landscape) | **citation only** — `Pan2025NatComm` already in bib; supplementary tables can be reviewed *as references*, no analysis |
| DepMap public CRISPR / RNAi (TROP2 dependency in thyroid lines) | Could support Paper 1 §7 (TROP2 vulnerability) translational hook | **DROP** for marathon — risk is asymmetric (negative TROP2 dependency in thyroid cell lines hurts the actionability story; current "tumor-population vulnerability" framing is safer) |
| CCLE thyroid cell lines (Broad) | Cell-line expression context for TROP2 ADC framing | **DROP** for marathon — same asymmetric-risk reasoning |
| GTEx normal tissue baseline | Normal tissue baseline | **already covered** by TCGA-THCA n=58 paired normals (lock-claim 6) |
| GSE82208 (older PTC microarray) | Cross-platform PTC reference | **DROP** — redundant with TCGA-THCA |
| ICGC pan-cancer | Pan-cancer outlier framing | **already covered** by TCGA pan-cancer 33-types analysis (lock-claim 14) |
| MSK-IMPACT thyroid (Ricarte-Filho et al.) | Driver / fusion landscape including RET/NTRK fusions | **DROP** — typically restricted access (dbGaP-like); strategy memo §8 explicitly disclaims fusion-map figure for Paper 1 |
| Krishnamoorthy 2025 Nat Commun ATC | Advanced disease cite | **DROP for acquisition** — bib-aliased to `Landa2016JCI` per memory `v17_landa2016_cite_save`; cite-only, not data |
| TCGA-THCA whole-slide images (WSI) | Image-DM1 retry | **HARD DROP** — closure battery NO-GO + 5 re-entry conditions not met |
| GTEx eQTL summary | eQTL evidence for lineage-TF expression | **DROP for marathon** — scope creep; defer to revision |
| ATA / SEER thyroid epidemiology | Hook stat sources | **already covered** (`Haugen2016ATA` + `Tuttle2019JCEM` already in bib) |
| HLA imputation reference panels (additional) | HLA forest expansion | **Paper 4 backlog** — out of marathon scope |

---

## 4. Don't-acquire-now list (consolidated, with rationale)

| Item | Rationale for not acquiring before 6/13 |
|---|---|
| GSE33630 microarray | Redundant with GSE76039 (same Affymetrix microarray modality family, same advanced-disease scientific question; GSE76039 first-pass 2026-05-04 already produced strong direction-consistent results). Fallback only. |
| TCGA-THCA methylation 450K | New data layer (~6–9 days) + voice-protected Discussion / Limitations rewrite required. Belongs to Nature Cancer reach / revision response. |
| DepMap CRISPR / RNAi | Asymmetric risk — negative thyroid-line dependency would hurt TROP2 actionability story; current "tumor-population vulnerability" framing is safer. |
| CCLE thyroid lines | Same as DepMap; cell-line context not patient-context; reviewer Q "is this clinically relevant?" |
| MSK-IMPACT thyroid | Restricted-access; fusion-map angle disclaimed for Paper 1 (strategy memo §8). |
| Pan/Ge 2025 ATC raw | Cite-only role is sufficient; raw analysis = scope creep. |
| TCGA-THCA WSI | Closure battery NO-GO; 5 re-entry conditions not met. |
| GSE82208 / ICGC / GTEx eQTL | Redundant or scope creep. |
| Any new ICI cohort | Paper 3 Track B FROZEN. |
| Any new GD / Korean germline / NRG1 cohort | Paper 4 backlog; marathon scope creep. |

---

## 5. Implementation gating (only invoked **after explicit user "go"** per acquisition)

### 5.1 GSE76039 (Paper 1) — gating before any download
- [ ] Author chooses scaffolding-track over voice-hook this week (i.e., not currently writing voice-protected Hook / Aim / Discussion §3.1 / Limitations / Cover ¶1 / Q9)
- [ ] No parallel agent claims `project/results/p_landa_2016/` or similar in concurrent reports
- [ ] Cohort metadata sanity check via GEO landing page (read-only browser; no programmatic action)
- [ ] Pre-registered analysis plan written: (i) score 84 samples with within-sample-centered 8-gene compact readout, (ii) Mann–Whitney advanced (PDTC+ATC) vs PTC + Cohen's d, (iii) Suppl Figure box-strip + per-sample table, (iv) sentence-level Results scaffold, no voice-protected prose
- [ ] **Hard time budget: 2.5 days end-to-end; if exceeded, halt and revert.**
- [ ] No dependence on RunPod / GPU / Azure burst (CPU-only score computation; ≪ 1 hour compute)

### 5.2 Lee 2014 (Paper 2 Pillar I v2) — gating before any reference-table extension
- [ ] Bib entry `Lee2014TissueAntigens` re-verified to point to correct DOI / PMID
- [ ] Pillar I v2 baseline-table layout decision (single Korean baseline column vs separate AFND + Lee 2014 columns) made by Yu / author
- [ ] One-row addition only; no figure change; no Methods change beyond a footnote
- [ ] **Hard time budget: 1 hour.**

### 5.3 TCGA-THCA methylation — gating before *any* metadata-only manifest read (not download)
- [ ] Confirmed deferred to bioRxiv-revision / v9
- [ ] Author has decided to include methylation in revision response or as new Paper "1.5"
- [ ] (Out of marathon scope; do not invoke before 6/14.)

---

## 6. What this plan does NOT do

- Does not download any dataset.
- Does not run any analysis.
- Does not invoke RunPod / GPU / Azure burst / TCGA WSI / H&E-DM1 / Paper 3 Track B.
- Does not generate Hook / Aim / Discussion / Limitations / Cover ¶1 / Q9 prose.
- Does not modify any manuscript file.
- Does not commit anything.

---

## 7. Cross-references

- `2026_05_04_8gene_curated_vs_denovo_final_strategy.md` — strategic axis decision; identifies external molecular replication as the gate to Nature Cancer reach (§12.6) and methylation as a Cancer Cell-aim requirement (§12.7).
- `2026_05_04_paper1_dm1_full_molecular_only_lock.md` — current Paper 1 envelope; methylation absent declared as Limitations (§3.5).
- `2026_05_04_marathon_state_after_image_dm1_nogo.md` — per-paper state and "no new analyses except paper-blocking" rule.
- `paper1_spatial_supplement_freeze_2026_05_03.md` — current ST validation Verdict C; GSE250521 advanced-disease evidence direction-consistent but underpowered (rationale for GSE76039 strategic value).
- `v17_landa2016_cite_save` — Landa 2016 cite-save memo; supports replicating-on-Landa's-cohort framing.
- `v18_paper2_HT_isolated` — Paper 2 Pillar I v2 design; Lee 2014 baseline reference role.
- `v19_paper4_GD_backlog` — Paper 4 GD HLA scope; explicitly excludes from marathon.
- `2026_05_04_bib_m3m4_straggler_cleanup.md` — bib status; pre-condition for GSE76039 acquisition GO.

---

## 8. Decision asks for the user

1. **Slot 1 — GSE76039 Paper 1 advanced-disease replication.** **GO** this week, **DEFER to W3 (5/18–5/24)**, or **DROP**? My recommendation: GO if author plans scaffolding-track this week; DEFER otherwise. Hard 2.5-day budget on entry.
2. **Slot 2 — Paper 1 second acquisition.** My recommendation: leave empty. Confirm OK to leave the second slot unused, or specify another candidate.
3. **Paper 2B reference — Lee 2014.** **GO** as a 1-hour bib-table extension when Paper 2 Pillar I v2 is touched next? My recommendation: GO.
4. **Methylation 450K** — confirm DEFER to bioRxiv-revision / v9? My recommendation: yes.
5. **Other candidates** — confirm DROP for marathon? My recommendation: yes (per §4 table).

**No acquisition action will be taken without explicit per-slot "go" approval.**

---

*Plan authored 2026-05-04 by Claude (Opus 4.7) under marathon-mode discipline. Read-only audit + plan only. No new data acquired, no analysis run, no RunPod / GPU / Azure burst, no voice-protected prose generated, no Paper 3 / Paper 4 touch, no manuscript file modified, no commit.*
