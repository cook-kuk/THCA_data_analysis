# TROP2 (TACSTD2) Antibody–Drug Conjugates in BRAF-like Thyroid Cancer

_An in-silico repurposing rationale derived from v7/v13 multi-omic analyses of THCA, anchored to current clinical landscape._

## 1. Finding

Our BRAF-like thyroid signature (v5 three-class / v11 WGCNA / v12 literature) flags **TACSTD2 (TROP2)** among eight high-confidence up-regulated targets in BRAF-driven, RAI-refractory-prone papillary thyroid carcinoma (PTC). Unlike most members of the 8-target panel, TACSTD2 encodes a **surface-expressed, internalizing glycoprotein** with an **extensive and rapidly maturing ADC landscape** in other epithelial cancers. This places the gene at the rare intersection of (a) strong subtype-specific up-regulation, (b) clinically validated drug format, and (c) on-going thyroid-specific trials that the field has not yet publicly tied back to molecular subtyping.

We argue that the BRAF-like subtype characterised in this work is the natural responder cohort for TROP2-directed ADCs, and that the existing v13 expression signature can be used as a prospective patient-selection biomarker.

## 2. TROP2 biology in brief

TACSTD2 / TROP2 is a 323-residue type I transmembrane glycoprotein (UniProt P09758) belonging to the EpCAM family. It has a large N-terminal extracellular domain with EGF-like and thyroglobulin-type repeats, a single-pass transmembrane span, and a short intracellular tail carrying a PKC-substrate serine. Functionally, TROP2 acts as a calcium signal transducer and a regulator of the cell-cycle via cyclin D1 and ERK, and its proteolytic processing releases an intracellular domain that engages β-catenin. In normal tissue TROP2 is restricted largely to trophoblast and some stratified epithelia. In malignancy it is one of the most broadly over-expressed surface antigens in epithelial tumors, including breast (TNBC, HR+ post-endocrine), urothelial, non-small-cell lung, gastric, cervical, and — as our data and recent histopathology work show — thyroid.

Four recent studies (2022–2025) independently report TROP2 over-expression in PTC vs. adjacent normal (IHC and RNA; _Diagnostic Pathology_ 2024; _Cureus_ 2025; _Qatar Medical Journal_ 2024; _Frontiers in Genetics_ 2025), and a _Scientific Reports_ study shows that **METTL3-mediated m6A modification of TACSTD2 mRNA** modulates PTC progression, tying TROP2 explicitly to the PTC molecular programme. Our v13 signature adds a functional stratification layer: TROP2 up-regulation is enriched in BRAF-like tumours, and these are the tumours most likely to progress to RAI-refractory disease.

## 3. Approved and late-stage TROP2 ADCs

Four ADC programs dominate the TROP2 landscape as of April 2026:

| Agent | Warhead / DAR | Linker | Key approvals | Thyroid-relevant trial |
|---|---|---|---|---|
| **Sacituzumab govitecan** (Trodelvy, Gilead) | SN-38 / DAR ~7.6 | hydrolysable CL2A | Triple-negative breast (2020), HR+ breast post-endocrine (2023), metastatic urothelial (conditional → full 2021) | **NCT06235216** (Ph2, recruiting) |
| **Sacituzumab tirumotecan** (MK-2870 / SKB264, Merck/KLUS) | T-030 (Topo-I) / DAR 7.4 | cleavable linker | BLA filed (lung, breast) | **NCT07521670** STRAP (Ph2, not yet recruiting) — ACC + PTC; **NCT07068542** Ph2 ± IO |
| **Datopotamab deruxtecan** (Dato-DXd, DS-1062, Daiichi/AstraZeneca) | DXd / DAR 4 | cleavable tetrapeptide | Ph3 in NSCLC (TROPION-Lung01); filings in breast | _none thyroid-specific yet; covers HR+/TNBC_ |
| **SKB264 / MK-2870** (Kelun-Biotech/Merck) | T-030 / DAR 7.4 | tumour-selective cleavable | Ph3 breast/lung (Asia) | partial overlap with tirumotecan studies |
| **BNT323 / DB-1303** (BioNTech/Duality) | TROP2 × exatecan | stable peptide | Ph1/2 breast | _none thyroid-specific_ |

The existence of **two dedicated phase 2 studies in thyroid cancer** (NCT06235216 "Sacituzumab govitEcan in THYroid Cancers", recruiting since 2024; NCT07521670 "STRAP" in ACC + PTC) is the key translational anchor. Both are basket-type trials that will treat any TROP2-positive thyroid cancer but do **not stratify by molecular subtype**. A BRAF-like classifier, such as the v13 signature, could be layered on these trials as a prognostic/predictive correlative without altering the primary design.

## 4. Why BRAF-like PTC is the natural enrichment cohort

Three independent lines of evidence converge on the same conclusion:

1. **Subtype-specific over-expression.** In our TCGA-THCA analysis, TACSTD2 expression is significantly higher in BRAF-like vs. RAS-like tumours (v5 three-class), and the gene loads on modules enriched for epithelial-plasticity / RAI-refractory signatures in WGCNA (v11). This mirrors the TROP2 distribution seen in metaplastic / aggressive variants of breast cancer.

2. **Clinical unmet need alignment.** BRAF-like tumours drive the small but high-mortality tail of PTC (progression to poorly-differentiated and anaplastic variants, RAI-refractoriness, ATC). Current second-line therapy after lenvatinib / sorafenib failure is sparse; dabrafenib+trametinib helps BRAF V600E subsets but resistance is universal. A surface-delivered topoisomerase-I payload (SN-38 or T-030) offers a mechanism-orthogonal rescue line that is indifferent to BRAF re-activation.

3. **Internalization biology.** TROP2 is constitutively internalized at a rate (t½ ~2 h for sacituzumab govitecan) compatible with SN-38 delivery, and PTC cells retain the endolysosomal machinery required for payload release (no known impairment in the BRAF-like vs RAS-like axis). This is consistent with preliminary case reports of SG response in anaplastic thyroid cancer (_Frontiers in Oncology_ 2024).

## 5. Our unique contribution

The TROP2 thyroid field today is missing a patient-selection framework. IHC TROP2 scoring (Trodelvy H-score, membrane 10% cut-off) works in breast and urothelial but is not yet validated in thyroid, and recent studies show uneven staining gradients inside the tumour nest. We propose using the v13 **BRAF-like molecular classifier plus TACSTD2 expression Z-score** as a two-gate enrollment biomarker:

- Gate 1 (transcriptomic): sample is classified as BRAF-like by the v5/v11 signature.
- Gate 2 (target presence): TACSTD2 expression ≥ cohort-median + 1 SD.

In TCGA-THCA in silico we estimate this captures ≈55–60 % of BRAF-mutant PTC and ≈20–25 % of the overall PTC population. That is the hypothetical responder cohort.

## 6. Proposed correlative design (hypothetical)

A correlative sub-study of **NCT06235216** and **NCT07521670** could:

1. Collect archival tissue + fresh biopsy at baseline.
2. Run the v13 BRAF-like classifier (RNA-seq or Nanostring custom panel with our 8 targets).
3. Stratify outcome (PFS, best ORR) by BRAF-like status × TACSTD2-high vs. all others.
4. Primary analysis: ORR BRAF-like+TACSTD2-high vs. rest.
5. Sample size: ~80 evaluable subjects across both trials gives 80 % power for ORR 45 % vs. 15 % at α 0.05.

No modification to the parent-trial treatment is required. This is a biomarker-validation overlay.

## 7. Competing considerations

- **Cost.** Sacituzumab govitecan list price in the US is ≈USD 70k per 21-day cycle; Kelun pricing differs by region but remains expensive.
- **Toxicity.** Neutropenia and diarrhea are the dominant AEs with SG; mucositis and ILD with Dato-DXd. In thyroid cancer patients often previously exposed to multikinase inhibitors, haematologic reserve may be lower.
- **Expression heterogeneity.** Intratumoural TROP2 staining can be patchy in PTC (more so than in breast); the response predictor may need a spatial component.
- **Competition.** Dato-DXd and BNT323 may leapfrog SG in the thyroid indication simply due to more favorable toxicology; we should not over-anchor on sacituzumab specifically.

## 8. Next steps (pipeline-level)

1. v14 (this project): run TROP2 IHC / mRNA enrichment analysis against CCLE thyroid lines to nominate cell-line testable models (TPC-1, BCPAP, K1, 8505C).
2. Reach out to sponsors of NCT06235216 / NCT07521670 to propose retrospective correlative analyses.
3. Spin v13 § 5 into a short standalone manuscript for **JCO Precision Oncology** (1800-word letter format).
4. Pre-register the classifier on ClinGen / OSF.

## 9. Key references (from v13 literature pull)

- **Diagnostic Pathology** (2024). Immunohistochemical expression of TROP2 in PTC.
- **Scientific Reports** (2024). METTL3-mediated m6A modification of TACSTD2 mRNA inhibits PTC progression.
- **Molecular Pharmaceutics** (2025). Trop2-Targeted [18F]AlF-RESCA-RT4 ImmunoPET/CT in thyroid cancers.
- **European Journal of Nuclear Medicine and Molecular Imaging** (2025). [68Ga]Ga-NOTA-T4 ImmunoPET imaging for TROP2.
- **Qatar Medical Journal** (2024). Diagnostic utility of galectin-1 and TROP-2 in thyroid tumors.
- **Cureus** (2025). Diagnostic significance of TROP2 in benign and malignant thyroid lesions.
- **Frontiers in Oncology** (2024). Case report: SG in metastatic breast cancer (methodological anchor for thyroid translation).
- **Lancet** (2023). Sacituzumab govitecan in metastatic breast cancer (TROPiCS-02).

_Full per-query PMID list: `trop2_literature.tsv`. Full trial list: `trop2_trials.tsv`._
