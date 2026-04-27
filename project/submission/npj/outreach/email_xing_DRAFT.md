# DRAFT — Mingzhao Xing outreach (English, formal academic)

⚠ **DRAFT — do not send as-is. Read-through and replace placeholders before sending.**

**To**: mxing1@jhmi.edu  _(verify current email; Xing may now be at SUSTech — check ORCID 0000-0003-0010-0265)_
**CC**: [Yu Kyungho email TBD]
**Subject**: Collaboration inquiry — TERT/BRAF transcriptomic axis extension (TCGA-THCA, npj submission preparation)

---

Dear Professor Xing,

I write to share computational work that builds directly on your foundational framework of BRAF V600E and TERT promoter mutations in papillary thyroid cancer (Liu R, Bishop J, Zhu G, et al. _JAMA Oncology_ 2017; Xing M, Liu R, Liu X, et al. _JCO_ 2014).

I am Seungho Cook, an independent computational researcher based in Seoul, working in partnership with [Yu Kyungho — full name TBD] at [Affiliation TBD]. We have re-analysed 513 TCGA-THCA primary tumours and recovered 36 TERT-promoter-mutated patients via the cBioPortal `thca_tcga_pub` mirror (Sanger-validated by TCGA Cell 2014). Three findings, in order of relevance to your work:

1. **Four-group stratification (BRAF / RAS / TERT⁺ / triple-negative)** yields logrank p = 4.92×10⁻⁶ for overall survival in TCGA-THCA, confirming the Liu–Xing genetic-duet framework in this cohort. Four-group p = 3.78×10⁻⁵.

2. **An unsupervised transcriptomic axis (DM1/DM2)** captures continuous differentiation/immune state, statistically distinct from mutation-based stratification (Spearman ρ = 0.49 with the canonical BRAF/RAS dichotomy). The TERT axis and the DM1/DM2 axis are independent (Fisher p = 0.65 within the DM-clustered subcohort, underpowered at n = 5 TERT⁺ overlap), consistent with two complementary stratification layers.

3. **An 8-gene RAI-responsiveness panel** (NIS, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, DIO1) trained on this axis achieves 5-fold CV AUC 0.954, outperforming BRAF V600E status alone (AUC 0.822, **ΔAUC = +0.132**). External validation on GSE76039 (PDTC + ATC) yields direction-correct AUC 0.974.

We are preparing for _npj Precision Oncology_ submission within the next 1–2 weeks. I would like to:

(a) **Cite your work as the foundational genetic framework** that motivated our 4-group analysis (Xing 2014, Liu 2017).

(b) **Inquire whether your 1,051-patient JHU cohort (Liu 2017) or the 507-patient cohort (Xing 2014) has expression data available**, which would substantially strengthen external validation of our 8-gene panel.

(c) **Explore co-author or correspondence collaboration** if your group is interested in the transcriptomic extension of the four-genotype framework.

I would be honoured to send the current manuscript for your review. Please let me know if there is interest, and I am happy to schedule a brief video call at your convenience (Seoul time, but flexible).

With great respect for your contributions to this field,

Seungho Cook
Independent Researcher (with part-time PhD affiliation, Seoul National University Graduate School of Convergence Science and Technology)
Seoul, Republic of Korea
Email: kukshomr@gmail.com

Co-investigator: [Yu Kyungho — full name TBD], [Affiliation TBD]

---

**Internal note (delete before sending)**:
- Replace `[Yu Kyungho — full name TBD]` and `[Affiliation TBD]` with confirmed values.
- Verify Xing's current institutional email (he split time between Hopkins and SUSTech historically).
- Consider attaching the manuscript v2 PDF if Xing prefers to review before responding.
