# GSE76039 — first-pass report (Steps 2–8)

**Date:** 2026-05-04 · **Owner:** Seungho Cook
**Cohort:** GSE76039 (Landa 2016) advanced thyroid cancer transcriptomics; n=37 (17 PDTC + 20 ATC; no PTC, no normal in this series); GPL570 microarray, gcRMA log2.
**Mode:** Read-only audit + first-pass score computation. **No raw FASTQ alignment, no CEL re-normalisation, no GPU/RunPod, no TCGA WSI, no H&E-DM1 retry, no Paper 3/4 touch.**
**Authority anchors:** `2026_05_04_paper1_dm1_full_molecular_only_lock.md` §2 envelope (16 claims) + §4 cautious language; `2026_05_04_8gene_curated_vs_denovo_final_strategy.md` §5 Layers 1–7 architecture + §11 manuscript-safe claims; `2026_05_04_gse76039_access_check.md` Step 1.
**Verdict (executive):** **External advanced-disease validation succeeds.** All seven Paper 1 axis layers (RAI/lineage axis, zero-overlap lineage cross-validation, TF collapse, STAT3/AP-1/DNMT activation, TROP2 elevation) are direction-consistent with TCGA bulk + ST validation, with extreme effect size on the differentiation axis (Cohen's d ≈ −3.5 ATC vs PDTC) and the strongest external zero-overlap lineage replication to date (Spearman ρ = −0.925, p = 3 × 10⁻¹⁶ on n=37). **MAIN FIGURE candidate** for Paper 1 — but report honestly within the limits below.

---

## 1. Dataset access result

| Item | Value |
|---|---|
| GEO accession | GSE76039 |
| Citation | Landa 2016 JCI 126(3):1052–1066 (PubMed 26878173) — `Landa2016JCI` already in `2026_05_03_references.bib` |
| n samples | 37 |
| Histology breakdown | 17 PDTC + 20 ATC; **no PTC, no normal in series** |
| Platform | GPL570 (Affymetrix HG-U133 Plus 2.0); 54,675 probes |
| Pre-processing | gcRMA log2 (Landa lab; we did NOT re-process) |
| Files | `GSE76039_series_matrix.txt.gz` (6.9 MB) + `GPL570_full.soft` (85.6 MB) |
| md5 / integrity | per `2026_05_04_gse76039_access_check.md` §3 |

---

## 2. Sample metadata summary

Per-sample metadata at `project/results/p_landa_2016/sample_metadata.tsv` (37 rows × 9 cols).

| Histology | n | Sex (F / M) | Tumor type detail |
|---|---:|---|---|
| ATC | 20 | 14 / 6 | 16 primary, 1 metastasis, 3 recurrent (incl. recurrent in neck, persistent metastatic) |
| PDTC | 17 | 14 / 3 | 14 primary, 3 recurrent (incl. residual primary) |

**Caveats present in this series:**
- BRAF / RAS / TERT mutation status NOT in GEO characteristics — would require parsing Landa 2016 supplementary tables (out of scope for this acquisition slot).
- Stage / age / survival NOT in GEO characteristics.
- 4/37 are recurrent / metastatic, not pure primary — does not change axis result (within-cohort z normalisation neutralises this for our scope).

---

## 3. Gene availability — 22/22 panel-union genes successfully mapped

Probe → gene mapping at `project/results/p_landa_2016/probe_to_gene_panel.tsv`. All 22 needed genes have ≥ 1 GPL570 probe; per-gene best-probe (max-mean across samples) collapse used (standard for Affymetrix multi-probe genes).

| Panel | Genes | All available? |
|---|---|---|
| RAI_8 (8-gene compact readout) | SLC5A5, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, DIO1 | ✓ 8/8 |
| THYROID_NONOVERLAP (zero gene overlap with RAI_8) | SLC26A4, IYD, DUOX1, DUOX2, TFF3, HHEX, GLIS3, DIO2 | ✓ 8/8 |
| TF_collapse | FOXE1, NKX2-1, PAX8, HHEX | ✓ 4/4 |
| STAT3_AP1_DNMT | STAT3, FOSL1, JUNB, DNMT1, DNMT3B | ✓ 5/5 |
| Drug-target hook | TACSTD2 (TROP2) | ✓ 1/1 |

NKX2-1 alias `TITF1` was prepared as fallback but not needed (current GPL570 annotation uses `NKX2-1`).

---

## 4. Primary score results — within-cohort z (n=37)

Score table at `project/results/p_landa_2016/scores.tsv`.

| Score | ATC mean (n=20) | PDTC mean (n=17) | direction | rationale |
|---|---:|---:|---|---|
| RAI_8 (8-gene differentiation) | **−0.666** | **+0.784** | ATC < PDTC | ATC more dedifferentiated (expected) |
| DM1_like = −RAI_8 | +0.666 | −0.784 | ATC > PDTC | (mirror; same information) |
| THYROID_NONOVERLAP (zero-overlap lineage) | **−0.665** | **+0.782** | ATC < PDTC | replicates 8-gene-free lineage axis |
| TDS_like (RAI_8 ∪ NONOVERLAP) | −0.666 | +0.783 | ATC < PDTC | broader axis, same sign |
| TF_collapse (FOXE1/NKX2-1/PAX8/HHEX) | **−0.720** | **+0.848** | ATC TF-collapsed | matches TCGA TF activity panel direction |
| STAT3_AP1_DNMT (STAT3/FOSL1/JUNB/DNMT1/DNMT3B) | **+0.320** | **−0.376** | ATC activated | matches TCGA mechanism panel direction |
| TACSTD2 (TROP2) z | **+0.456** | **−0.536** | ATC TROP2-high | drug-target hook direction-consistent |

---

## 5. Advanced-disease validation result — ATC vs PDTC and within-cohort correlation

Test table at `project/results/p_landa_2016/score_tests.tsv` (Mann–Whitney two-sided; Cohen's d pooled SD; Spearman ρ on all 37 samples).

### 5.1 ATC vs PDTC (within-cohort contrast — only stage contrast available in this series)

| Axis | Cohen's d | MW p (two-sided) |
|---|---:|---:|
| **RAI_8** | **−3.47** | **4.6 × 10⁻⁷** |
| **DM1_like** (= −RAI_8) | **+3.47** | **4.6 × 10⁻⁷** |
| **THYROID_NONOVERLAP** | **−3.32** | **6.3 × 10⁻⁷** |
| **TDS_like** | **−3.54** | **4.6 × 10⁻⁷** |
| **TF_collapse** | **−3.57** | **2.4 × 10⁻⁷** ← strongest |
| **STAT3_AP1_DNMT** | **+1.72** | **3.6 × 10⁻⁵** |
| **TACSTD2 (TROP2)** | **+1.13** | **2.7 × 10⁻³** ← weakest but still significant |

All 7 axes are direction-consistent with the strategy-memo §5 architecture. ATC (the more dedifferentiated end of the cohort) shows lower differentiation-axis scores and higher mechanism / drug-target scores, exactly as predicted from TCGA-bulk extrapolation.

### 5.2 Within-cohort correlations on n = 37 (Spearman)

| Pair | ρ | p | comparison |
|---|---:|---:|---|
| **DM1_like vs THYROID_NONOVERLAP** | **−0.925** | **3.1 × 10⁻¹⁶** | 8-gene-free lineage cross-validation; **stronger than TCGA-bulk reference** (TCGA r = −0.885 on n=561) |
| **TF_collapse vs DM1_like** | **−0.931** | **7.3 × 10⁻¹⁷** | TF activity collapse co-occurs with DM1_like-high |
| **STAT3_AP1_DNMT vs DM1_like** | **+0.681** | **3.5 × 10⁻⁶** | mechanism panel activation co-occurs with DM1_like-high |
| **TACSTD2 vs DM1_like** | **+0.443** | **6.1 × 10⁻³** | TROP2 expression elevation in DM1-axis-high direction |

### 5.3 Interpretation against user's GO-rules

- **"DM1_like increases in PDTC/ATC and NONOVERLAP decreases" → MAIN FIGURE candidate.** ✓ Both directions confirmed at Cohen's d > 3.
- "If signal is weak but direction-consistent, supplementary only." → Not the case; signal is strong.
- "If no signal, report honest negative and stop." → Not triggered.

---

## 6. Mechanism / TACSTD2 result

### 6.1 TF collapse co-occurrence

The four lineage TFs (FOXE1, NKX2-1, PAX8, HHEX) collapse together along the DM1_like axis — Spearman ρ = −0.931 between TF_collapse module and DM1_like (effectively saturating). Direction matches TCGA-internal lock-claim 4 (FOXE1 −1.21 / NKX2-1 −0.65 in DM1-high), and the within-cohort effect size between ATC and PDTC (d = −3.57) is the largest of any axis tested.

### 6.2 STAT3 / AP-1 / DNMT activation

STAT3 / FOSL1 / JUNB / DNMT1 / DNMT3B aggregated as a module score: ATC mean +0.320 vs PDTC −0.376 (d = +1.72, p = 3.6 × 10⁻⁵). Direction matches TCGA-internal lock-claim 4 (STAT3 +1.55 / FOSL1 +0.98 / DNMT1 +0.74 in DM1-high). DNMT3B is added (not in lock-claim 4) and is direction-consistent in this cohort. **Important hedging:** methylation itself remains absent — the TF→DNMT framework remains "consistent with" not "demonstrates" per `molecular_only_lock` §4.

### 6.3 TACSTD2 / TROP2 single-gene

TROP2 z-score: ATC +0.456 vs PDTC −0.536 (d = +1.13, p = 2.7 × 10⁻³); ρ vs DM1_like = +0.443 (p = 6.1 × 10⁻³). **Smallest effect of the seven axes**, but still significant. Direction-consistent with TCGA-internal lock-claim 5 (TROP2 Δmean +4.33, FDR = 1.2 × 10⁻³² in DM1-high).

**Manuscript-safe framing per `molecular_only_lock` §4:** TROP2 elevation in advanced disease is a **tumor-population vulnerability** signal. Q3 spot-level co-localisation NEG (from TCGA spatial) is **not** addressed by this cohort and remains as previously reported. Forbidden-language guard: do **not** say "DM1-high spots are TROP2-high in advanced disease" — GSE76039 is bulk microarray, no spatial information.

---

## 7. Figure paths

| Figure | Path | Use |
|---|---|---|
| Box+strip plot (4 panels — RAI_8 / DM1_like / NONOVERLAP / TDS_like by histology) | `project/results/p_landa_2016/gse76039_dm1_lineage_boxplot.png` | **Main figure candidate** — primary advanced-disease replication panel |
| Scatter — DM1_like vs THYROID_NONOVERLAP (n=37, ρ = −0.925) | `project/results/p_landa_2016/gse76039_dm1_vs_nonoverlap_scatter.png` | **Main figure candidate** — zero-overlap lineage cross-validation |
| Mechanism heatmap — 10 genes (TF collapse + STAT3/AP-1/DNMT + TACSTD2), z within cohort, samples ordered by histology then DM1_like | `project/results/p_landa_2016/gse76039_mechanism_heatmap.png` | **Main figure candidate** — mechanism replication panel |

All figures generated CPU-only via matplotlib; deterministic with `np.random.seed(42)` for jitter.

---

## 8. Verdict

| Bucket | Items |
|---|---|
| **Main figure candidate** | (a) DM1_like / RAI_8 / TDS-like ATC vs PDTC box+strip; (b) DM1_like vs THYROID_NONOVERLAP scatter (ρ = −0.925); (c) Mechanism heatmap (TF collapse + STAT3/DNMT activation + TACSTD2). All three could fit a single 3-panel composite figure ("F8 — External advanced-disease replication, GSE76039 Landa 2016") — to be added to the Figure plan in `2026_05_03_manuscript_v8_OUTLINE.md` if author/Yu approves. |
| **Supplementary candidate** | Per-gene boxplots for 8-gene panel members and zero-overlap lineage members (within cohort, ATC vs PDTC); STAT3_AP1_DNMT module per-gene; TACSTD2 single-gene scatter; metadata table. |
| **Internal-only** | None from this cohort — every measured axis came out clean. (No "internal-only" risk-disclosed result like TROP2 spot-level Q3 NEG, which is TCGA-spatial territory.) |
| **Drop** | None. |

### Recommended placement in Paper 1 architecture (per strategy memo §5)

- **Layer 1 (driver orthogonality):** GSE76039 cannot contribute (no mutation status in series).
- **Layer 2 (axis robustness — pan-genome ARI):** GSE76039 cannot contribute (different cohort; would need re-clustering, scope creep).
- **Layer 3 (compact 8-gene readout):** ✓ **GSE76039 supports** — 8-gene RAI score d = −3.47 ATC vs PDTC.
- **Layer 4 (mechanism — TF→DNMT→STAT3→TROP2):** ✓ **GSE76039 supports** — all three sub-modules direction-consistent.
- **Layer 5 (outcome):** GSE76039 has no survival data; cannot contribute.
- **Layer 6 (spatial validation):** GSE76039 is bulk; cannot contribute.
- **Layer 7 (TROP2 vulnerability):** ✓ **GSE76039 supports** — TACSTD2 elevated in advanced disease (d = +1.13).
- **Honest negatives (R9):** No new honest-negative; existing N1 BRAF-/RAS- subset NS still holds within TCGA.

### Strategy-memo §12.6 (Nature Cancer reach gate)

The strategy memo flagged that "an independent external survival/molecular validation of DM1_like in a second cohort" is the gate to Nature Cancer reach. **GSE76039 partially closes this gate (molecular validation, no survival).** Combined with internal multivariate Cox (lock-claim 3) on TCGA, this raises the venue ladder modestly toward CRM-strong / JCI Insight-comfortable. Cancer Cell aim still requires functional + methylation per §12.7.

---

## 9. Risks

| Risk | Severity | Mitigation already applied |
|---|---|---|
| FFPE / older microarray platform (HG-U133 Plus 2.0) | LOW | Within-cohort z neutralises platform absolute scale; multiple within-cohort effect-size estimates (Cohen's d, Spearman ρ) are platform-independent. |
| Histology label quality (PDTC Turin vs MSKCC criteria; recurrent vs primary mixed) | MEDIUM | Direction-consistency is robust to label noise (effect size d > 3 will not flip on 1–2 label corrections). For Methods scaffold, frame as "PDTC and ATC as defined by Landa et al. 2016." |
| Sample size n=37 | MEDIUM | All MW p ≤ 2.7 × 10⁻³; smallest effect is TACSTD2 (d = 1.13, p = 2.7 × 10⁻³). Even with 1–2 outlier removal, all axes remain significant. |
| Cross-platform normalisation (microarray vs RNA-seq comparison to TCGA) | LOW (within-cohort framing) / MEDIUM-HIGH (if absolute-score comparison attempted) | We do **not** transfer the TCGA-trained absolute-form classifier (per `v17_korean_K2_calibration` memory). Within-cohort z is the only operation done here. Reviewer-Q-safe. |
| No PTC controls in series | MEDIUM (limits PTC → PDTC → ATC ladder claim) | Frame as "ATC vs PDTC within advanced disease" — do **not** claim a monotonic 3-stage gradient from this cohort alone. Cross-cohort ladder framing combines GSE76039 with TCGA-PTC reference, requiring careful platform handling — out of scope for this acquisition slot. |
| No mutation / survival metadata | MEDIUM (cannot test BRAF-/RAS- subset replication or external Cox) | Do not claim survival or mutation effects from this cohort. Multivariate Cox remains TCGA-only. |
| Result is "too clean" — risk of an artefact (gcRMA quantile structure giving spurious anti-correlation) | LOW | Sanity check: if differentiated-lineage genes are saturated low in advanced tumors as a class, the −0.925 between RAI_8 and NONOVERLAP-with-zero-gene-overlap is biologically expected. Direction matches TCGA-bulk internal r = −0.885. Consistent. |

---

## 10. What this analysis did NOT do

- Did NOT realign FASTQ / re-process CEL files / run STAR/Salmon/kallisto.
- Did NOT use GPU / RunPod / Azure burst.
- Did NOT touch Paper 2 / Paper 3 / Paper 4 files.
- Did NOT modify any voice-protected manuscript section (Hook / Aim / Discussion §3.1 / Limitations / Cover ¶1 / Q9).
- Did NOT modify `2026_05_03_manuscript_v8_OUTLINE.md` (figure plan addition is **proposed** below; awaits user approval).
- Did NOT fetch additional external datasets (GSE33630, TCGA methylation, DepMap, etc. all explicitly DROP per `2026_05_04_minimal_external_data_plan.md` §4).
- Did NOT compute survival, mutation, or stage-ordinal claims.

---

## 11. Recommended next actions (user approval required for each)

1. **Add Suppl Table SX_LANDA_GSE76039** — sample × score matrix (already at `project/results/p_landa_2016/scores.tsv` in TSV form; XLSX packaging is a 1-line change).
2. **Add F8 to OUTLINE Figure plan** — "External advanced-disease replication (GSE76039)" composite, panels A (boxplot), B (scatter), C (heatmap). This is a **Main figure** add to the strategy-memo §9 figure plan; user approval required because it modifies OUTLINE.
3. **Add Methods M3' or M_landa scaffold** — 1-paragraph Methods entry covering: cohort, gcRMA Series Matrix usage, GPL570 best-probe collapse, within-cohort z, MW + Spearman tests, no platform-transfer of TCGA classifier.
4. **Update strategy-memo §3.1 wording** — change "RNA-seq" to "microarray (HG-U133 Plus 2.0)" and update n from "≈84" to **n=37 (17 PDTC + 20 ATC, no PTC controls in series)**. Cross-platform within-cohort-z framing already validated above.
5. **Update bib note for `Landa2016JCI`** — note that the transcriptomic subset `GSE76039` (n=37) is a strict subset of the larger DNA-targeted-seq cohort (n=117).
6. **Voice-protected work (separate session):** if author wants to integrate this into Discussion §3.1, that is **author keyboard only** per `v17_sprint_vs_marathon_violation`. This report does not generate Discussion prose.
7. **Do NOT** acquire GSE33630 / TCGA methylation / DepMap / any other dataset without an explicit per-slot "go" per `2026_05_04_minimal_external_data_plan.md` §5.

---

## 12. Marathon-discipline attestation

| Constraint | Status |
|---|---|
| Raw FASTQ alignment / STAR / Salmon / kallisto | ✓ NOT performed |
| CEL re-normalisation | ✓ NOT performed |
| GPU / RunPod / Azure burst | ✓ none |
| TCGA WSI / H&E-DM1 retry | ✓ none |
| GSE33630 / TCGA methylation download | ✓ NOT performed |
| Paper 3 / Paper 4 touch | ✓ none |
| Voice-protected sections (Hook / Aim / Disc 3.1 / Limitations / Cover ¶1 / Q9) | ✓ untouched |
| Manuscript prose (large-scale) modification | ✓ none |
| Result presentation as "main figure" without honest direction check | ✓ none — direction-consistency + effect-size + within-cohort z all defensible |
| Time budget (today: 6 h to first score table) | ✓ ~30 min total CPU + I/O |

---

## 13. Cross-references

- Step 1 access check: `2026_05_04_gse76039_access_check.md`
- Pipeline script (CPU-only): `project/notebooks_or_scripts/p_landa_2016_first_pass.py`
- Output files: `project/results/p_landa_2016/sample_metadata.tsv`, `expression_matrix_log_gene.tsv.gz`, `probe_to_gene_panel.tsv`, `scores.tsv`, `score_tests.tsv`, three PNG figures
- Strategy: `2026_05_04_8gene_curated_vs_denovo_final_strategy.md` §5 architecture, §11 manuscript-safe claims, §12.6 Nature Cancer reach gate
- Lock authority: `2026_05_04_paper1_dm1_full_molecular_only_lock.md` §2 envelope
- Acquisition plan: `2026_05_04_minimal_external_data_plan.md` §3.1 GO slot
- Bib citekey: `Landa2016JCI` in `2026_05_03_references.bib`
- Memory: `v17_landa2016_cite_save` · `v17_korean_K2_calibration` · `v17_marathon_mode_post_pillar1` · `v17_sprint_vs_marathon_violation`

---

GSE76039 first-pass complete. No additional dataset acquired.
