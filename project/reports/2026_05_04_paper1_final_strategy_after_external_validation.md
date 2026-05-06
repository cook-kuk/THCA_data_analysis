# Paper 1 final strategy after GPL570 external expression validation

**Date:** 2026-05-04
**Scope:** Paper 1 only. Strategy/scaffold integration after external expression validation.
**Canonical inputs:** `2026_05_04_external_expression_FULL_RESULTS.md`, `2026_05_04_gpl570_validation_FULL_VIEW.md`, `2026_05_04_8gene_curated_vs_denovo_final_strategy.md`, manuscript v8 scaffold reports.
**Discipline:** No new datasets, no new GEO search, no new analysis, no raw CEL/FASTQ processing, no manuscript voice-protected prose.

---

## 1. Executive verdict

**Lineage axis strengthened.** The GPL570 external expression validation sweep is a hit for the driver-orthogonal thyroid-lineage differentiation axis. GSE33630, GSE29265, GSE65144, and GSE53157 reproduce the expected RAI_8 / THYROID_NONOVERLAP / TDS_like / TF_collapse direction, and the zero-overlap THYROID_NONOVERLAP anti-correlation with DM1_like is replicated across cohorts.

**TROP2 demoted.** TACSTD2/TROP2 bulk-microarray replication is failed or heterogeneous in the GPL570 validation pack. TROP2 must not remain in the title or main claim. It can remain as exploratory/supplement/internal and as a hypothesis-generating tumor-level observation from TCGA, but not as a general external-validation headline.

**No more datasets.** The external expression validation slot is closed for the marathon. GSE76039 plus the GPL570 pack is sufficient for the current claim envelope. Do not add GEO datasets, DepMap/CCLE/PRISM, methylation, WSI, H&E-DM1, or RunPod work.

---

## 2. New main claim

**Primary claim:** “A driver-orthogonal thyroid-lineage differentiation axis stratifies thyroid cancer beyond canonical driver mutations.”

This is the Paper 1 identity. The novelty is not a new 8-gene discovery panel and not a TROP2-targetability claim. The novelty is the driver-orthogonal, lineage-silencing transcriptional axis with external expression replication.

---

## 3. Title ranking

1. **A thyroid-lineage differentiation axis stratifies thyroid cancer beyond canonical driver mutations**
2. **A driver-orthogonal thyroid-lineage axis recapitulates advanced-disease dedifferentiation in thyroid cancer**
3. **A transcriptional differentiation axis defines lineage silencing in thyroid cancer**

**Deprecated old title:** “A transcriptional differentiation axis stratifies thyroid cancer orthogonally to canonical driver mutations and identifies a tumor-population TROP2 vulnerability”

Reason for deprecation: the lineage axis externally replicates; TACSTD2/TROP2 bulk-microarray replication does not.

---

## 4. Role of the 8-gene readout

- The 8-gene panel remains a compact RAI-lineage readout: SLC5A5, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, DIO1.
- It is not the title and not the discovery premise.
- It operationalizes the axis for scoring and interpretation.
- Its defensibility comes from canonical thyroid/RAI biology, retention of most TDS-16 discriminative signal, drop-one-out robustness, and zero-overlap validation using THYROID_NONOVERLAP.

---

## 5. External expression validation

**Anchor:** GSE76039 (Landa 2016; GPL570; PDTC/ATC advanced-disease anchor).

**GPL570 pack:** GSE33630, GSE29265, GSE65144, GSE53157.

**Core result:** DM1_like and THYROID_NONOVERLAP show replicated zero-overlap anti-correlation across cohorts, supporting a true lineage-silencing axis rather than an 8-gene panel artifact.

**Direction-consistent modules:** RAI_8, THYROID_NONOVERLAP, TDS_like, TF_collapse, and STAT3_AP1_DNMT are direction-consistent across testable advanced-disease contrasts. GSE53157 is supportive/sensitivity only because the PDTC arm is underpowered.

**Placement:** main figure candidate or strong supplementary figure. Recommended: make the external expression validation panel a main-figure candidate if space permits; otherwise make it the first strong supplement and cite it directly in Results.

---

## 6. Mechanism layer

Mechanism remains supportive, not causal proof:

- TF_collapse supports lineage transcription-factor loss.
- STAT3/AP-1/DNMT directionality supports the proposed mechanism arm.
- Use “consistent with” and “supports”; do not use “demonstrates” for TF-to-DNMT or methylation mechanisms.
- TACSTD2/TROP2 is no longer part of the main mechanism headline. It can remain as exploratory tumor-level biology or internal/supplementary follow-up.

---

## 7. TROP2 demotion note

TROP2/TACSTD2 status after external validation:

- **Not title.**
- **Not main claim.**
- **Not general external validation.**
- **Not “DM1-high spots are TROP2-high.”**
- **Not “TROP2-targetable” as a headline.**
- May remain as exploratory/supplement/internal: TCGA bulk signal, tumor-vs-normal signal, and possible future ADC hypothesis.
- GPL570 bulk-microarray result is heterogeneous/negative, so any downstream text must explicitly avoid implying external replication of TROP2.

---

## 8. Safe wording

- “external expression validation”
- “direction-consistent lineage silencing”
- “advanced-disease replication”
- “zero-overlap module validation”
- “driver-orthogonal thyroid-lineage differentiation axis”
- “compact RAI-lineage readout”
- “supportive mechanism arm”

## 9. Forbidden wording

- “TROP2-targetable” as headline
- “TROP2 vulnerability” in the title or main claim
- “progression proven”
- “clinical validation”
- “survival validation”
- “fusion validation”
- “DM1-high spots are TROP2-high”
- “H&E-inferable”
- “all risks resolved”
- “Cancer Cell-ready” or venue-reach claims

---

## 10. Figure placement decision

Recommended Paper 1 figure strategy:

- **F1:** Driver landscape orthogonality.
- **F2:** Axis robustness and compact 8-gene readout.
- **F3 or F4:** External expression validation across GSE76039 plus GPL570 pack.
- **Mechanism:** keep TF_collapse and STAT3/AP1/DNMT support; TROP2 is supplement/internal only.
- **Spatial:** supportive supplement unless main-figure space opens.

External validation figure panels:

- Boxplots: RAI_8, DM1_like, THYROID_NONOVERLAP by histology/cohort.
- Scatter: DM1_like vs THYROID_NONOVERLAP within cohort.
- Forest: direction consistency across dataset × contrast × score.
- QC: gene-panel coverage heatmap in supplement or methods supplement.

---

## 11. Commit proposal only

**E1 strategy integration**
```
project/reports/2026_05_04_paper1_final_strategy_after_external_validation.md
```

**E2 manuscript scaffold integration**
```
project/reports/2026_05_03_manuscript_v8_OUTLINE.md
project/reports/2026_05_03_methods_M1_M11_scaffold.md
project/reports/2026_05_03_methods_M1_M11_prose.md
project/reports/2026_05_03_figure_captions_all.md
```

**E3 external validation report + scripts**
```
project/reports/2026_05_04_external_expression_FULL_RESULTS.md
project/results/p_external_expression_validation/scripts/run_external_validation.py
project/results/p_external_expression_validation/scripts/make_figures.py
```

**E4 processed outputs + figures**
```
project/results/p_external_expression_validation/sample_metadata.tsv
project/results/p_external_expression_validation/external_gene_coverage.tsv
project/results/p_external_expression_validation/external_score_tests.tsv
project/results/p_external_expression_validation/external_direction_consistency.tsv
project/results/p_external_expression_validation/external_spearman.tsv
project/results/p_external_expression_validation/external_sample_scores.tsv.gz
project/results/p_external_expression_validation/*_expression_gene_log.tsv.gz
project/results/p_external_expression_validation/*.png
```

**E5 gitignore**
```
.gitignore
```

**Exclude:** `project/results/p_external_expression_validation/raw/`, Paper 2/3/4 files, voice-protected prose, WSI/tiles/embeddings/RunPod runtime.

---

## 12. Final operating instruction

External validation is integrated as a lineage-axis win, not a TROP2 win. Return to the user-written Hook after scaffold updates.
