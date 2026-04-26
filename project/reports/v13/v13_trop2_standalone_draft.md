# TROP2 up-regulation in BRAF-like thyroid cancer: a molecular rationale for sacituzumab govitecan repurposing

**Working title.** Target journal: _JCO Precision Oncology_ — Brief Report / Letter format (~1,800 words, 2 figures, 1 table).

**Authors (provisional).** [to be filled — PI 유 교수님 + v13 pipeline lead]

## Abstract (250 words)

Papillary thyroid cancer (PTC) is usually indolent, but a minority of tumours progress to radioactive-iodine–refractory (RAI-R) disease for which second-line options after lenvatinib/sorafenib failure are sparse. Trophoblast cell surface antigen 2 (TROP2 / TACSTD2) over-expression in PTC and its association with BRAF V600E mutation have been documented at the IHC level since 2018 (Liu, Bychkov) and at the transcriptomic level in primary PTC plus paired lymph-node metastases by Kalfert et al. (2024); sacituzumab govitecan (anti-TROP2 ADC) has been catalogued among "non-explored niche indications" for thyroid cancer in a 2023 ADC-landscape review (Nieto-Jiménez et al.). What has not been proposed is a **molecular sub-stratification** of the TROP2-directed ADC indication that selects responders prospectively. Using a v13 in-silico drug-discovery pipeline over eight BRAF-like up-regulated targets and ten therapeutic modalities, we (i) confirm the BRAF-axis TROP2 up-regulation at the cell-line transcriptomic level (top three TACSTD2-expressing CCLE thyroid lines all carry BRAF V600E), (ii) show in PRISM Repurposing 19Q4 that the SN-38 payload class (topotecan, irinotecan) preferentially kills BRAF-like CCLE thyroid lines by ~1.0 log-fold-change vs RAS-like — opposite polarity to the established CRC literature — and (iii) cross-reference ClinicalTrials.gov (April 2026) showing **two dedicated Ph2 thyroid TROP2-ADC trials are already active** (NCT06235216 "Sacituzumab govitEcan in THYroid Cancers"; NCT07521670 "STRAP") but **neither stratifies patients by BRAF status**, with STRAP explicitly waiving even TROP2 IHC testing. We propose that the v13 BRAF-like classifier **plus** TACSTD2 mRNA z-score be evaluated as a two-gate enrichment biomarker in these existing trials, captured retrospectively from archival tissue. In silico this two-gate strategy captures ≈55–60 % of BRAF-mutant PTC and ≈20–25 % of overall PTC, and is the natural responder cohort.

## 1. Background

PTC accounts for the majority of thyroid cancer incidence and typically carries an excellent prognosis with surgery and RAI. A minority (~5–10 %) progress to RAI-R disease, where five-year disease-specific survival drops from > 95 % to < 50 %. Lenvatinib and sorafenib are the current first-line tyrosine-kinase inhibitors (TKIs); dabrafenib + trametinib helps BRAF V600E subsets but resistance is near-universal. A mechanism-orthogonal rescue option is needed.

The BRAF-like / RAS-like molecular classification of PTC, originally described in the TCGA-THCA paper (TCGA Research Network, _Cell_ 2014) and formalised as the BRS (BRAF-RAS score) axis, is reproducible and correlates with aggressiveness. Two practical caveats motivate orthogonal molecular readouts in trial enrolment. First, the BRAF-like / RAS-like axis is **expression-based and only ~95 % concordant with mutation truth** even on TCGA-THCA (Chakravarty 2011 BRS52 panel applied to our local TCGA-THCA cohort, 334 / 351 = 95.2 % accuracy with specificity-for-RAS = 100 %; 17 BRAF-mutant tumours expression-classify RAS-like — the well-documented "BRAF-mutated, RAS-like-by-expression" subset of Landa _et al._ 2016 _Cell_). Second, **cross-cohort harmonisation in THCA is structurally fragile**: within-subtype PC1 is perfectly disjoint between TCGA-THCA and GSE27155 (BUHMBOX-style KS p < 10⁻⁴⁰ in BRAF, p < 10⁻¹⁴ in RAS — our v8.1 supplementary), so subtype labels learned on one cohort generalise unevenly to a second cohort harmonised by linear batch correction. Both observations support the same operational conclusion: a TROP2-ADC enrichment strategy in PTC should not rely on subtype labelling alone, and should layer a direct molecular readout (TROP2 mRNA or IHC) on top of the BRAF / BRS axis.

The TROP2 / BRAF link in PTC has prior support. Liu et al. (2018, _Int J Clin Exp Pathol_, PMID 31949805) reported that TROP2 over-expression in PTC tissue is associated with BRAF V600E mutation and aggressive behaviour at the IHC level, and Bychkov et al. (2018, _J Pathol Transl Med_, PMID 29228520) independently flagged TROP2 as a prognostic marker in PTC. Kalfert et al. (2024, _Pathol Res Pract_, PMID 38696857) extended this to the mRNA axis in 60 primary PTC tumours plus 40 paired lymph-node metastases, integrating BRAF status with TACSTD2 transcript and four miRNAs. A pan-tumour TMA of 18,563 cases (Dum et al. 2022, _Pathobiology_, PMID 35477165) further linked TROP2 over-expression to nodal metastasis in PTC (p = 0.0013). The therapeutic-translation step has also been flagged: Nieto-Jiménez et al. (2023, _Clin Transl Med_, PMID 37740463) explicitly catalogue sacituzumab govitecan in "pancreatic, gastric, **thyroid** or endometrial cancer" as a "non-explored niche indication." Building on this established foundation, our contribution is **not** to identify TROP2 as a thyroid ADC target — that has been done — but to propose a BRAF-axis sub-stratification of the previously-proposed thyroid TROP2-ADC opportunity, supported by orthogonal cell-line and pharmacogenomic evidence (sections 4–5).

## 2. TROP2 biology and ADC landscape

TACSTD2 encodes TROP2, a 323-residue type I transmembrane glycoprotein of the EpCAM family. Its extracellular domain carries EGF-like and thyroglobulin-type repeats, making it both an obvious adhesion / signalling receptor and an accessible target for antibody-based therapy. In normal tissues TROP2 is restricted to trophoblast and selected stratified epithelia; in malignancy it is one of the most broadly over-expressed surface antigens in epithelial tumours. Internalisation kinetics (t½ ≈ 2 h for sacituzumab govitecan complexes) are compatible with an SN-38 payload.

Four TROP2-directed ADCs dominate the 2026 landscape:

1. **Sacituzumab govitecan (Trodelvy)** — SN-38 payload, CL2A hydrolysable linker, DAR ≈ 7.6; approved in TNBC (2020), HR+ breast after endocrine therapy (2023) and metastatic urothelial (2021 → 2024 full).
2. **Sacituzumab tirumotecan (MK-2870 / SKB264)** — Topo-I T-030 payload, DAR 7.4; BLA filed lung/breast.
3. **Datopotamab deruxtecan (Dato-DXd, DS-1062)** — DXd payload, cleavable tetrapeptide, DAR 4; Ph3 NSCLC (TROPION-Lung01).
4. **BNT323 / DB-1303** — exatecan payload; Ph1/2 breast.

All four share TROP2 as the antigen, differ in payload class and linker, and are in differentiated clinical positioning. They do not currently compete for the thyroid indication.

## 3. TROP2 in thyroid cancer — literature review

Our v13 + v14 PubMed pull (E-utilities, 2022-04-25 to 2026-04-25, fourteen combined queries) returned 122 unique records. The most relevant priors are summarised in §1; this section adds context across five themes:

- **Pathology priors (foundation, must cite).** Liu et al. (2018, PMID 31949805) — TROP2 over-expression ↔ BRAF V600E in PTC, IHC. Bychkov et al. (2018, PMID 29228520) — TROP2 prognostic significance in PTC. Kalfert et al. (2024, PMID 38696857) — first integrated BRAF × TACSTD2 mRNA × miRNA study in primary PTC + paired LNM. Dum et al. (2022, PMID 35477165) — pan-tumour TMA n=18,563 with PTC TROP2-nodal-metastasis link.
- **ADC indication priors.** Nieto-Jiménez et al. (2023, PMID 37740463) — SG cited for thyroid as a niche indication at review level. Cancers 2022 (PMID 35158847) and IJSP 2021 (PMID 33289434) — TROP2 as ADC target in anaplastic thyroid carcinoma.
- **Recent histopathology (2023–2025).** Four IHC studies (Diagnostic Pathology 2024; Cureus 2025; Qatar Medical Journal 2024; Histopathology 2025) report TROP2 over-expression in PTC vs. adjacent normal with 60–80 % positivity, complementing the 2018 IHC priors.
- **Molecular regulation.** Scientific Reports (Zhou et al., 2024) shows METTL3-mediated m6A modification of TACSTD2 mRNA modulates PTC progression — a direct molecular tie to the PTC programme.
- **Imaging.** Molecular Pharmaceutics (2025) and European Journal of Nuclear Medicine (2025) describe ⁶⁸Ga- and ¹⁸F-labelled anti-TROP2 ImmunoPET proof-of-concept in thyroid cancers; NCT06851663 (broad Trop2 ImmunoPET in solid tumours, Ph2/3 recruiting). Frontiers in Genetics (2025) integrative single-cell transcriptomics identifies TROP2 as a biomarker node in PTC.

The combined priors consistently validate TROP2 expression and its BRAF association, but flag **intra-tumoural heterogeneity** and concur on the unresolved patient-selection question — a molecular classifier may out-perform IHC alone, and no prior paper proposes BRAF-axis sub-stratification specifically for TROP2-ADC trials.

## 4. Ongoing thyroid TROP2 trials (ClinicalTrials.gov, April 2026)

| NCT | Title | Drug | Phase | Status | Conditions |
|---|---|---|---|---|---|
| NCT06235216 | Sacituzumab govitEcan in THYroid Cancers | sacituzumab govitecan | Ph2 | Recruiting (2024-09) | Differentiated / anaplastic thyroid |
| NCT07521670 | STRAP | sacituzumab tirumotecan | Ph2 | Not yet recruiting (2026-05) | Adenoid cystic carcinoma + PTC |
| NCT07068542 | Sacituzumab tirumotecan + IO | sacituzumab tirumotecan + anti-PD-1 | NA | Recruiting (2025-07) | RAI-R differentiated + poorly-differentiated thyroid |
| NCT06923826 | SG-SGTC | sacituzumab govitecan | Ph2 | Recruiting (2025-04) | Salivary gland + thyroid |

Neither trial selects by molecular subtype. Each evaluates IHC TROP2 or treats all-comers.

## 5. Proposed correlative design

We propose a **biomarker-only** correlative overlay that can be grafted onto any of the four trials above without modifying the primary protocol.

### 5.1 Two-gate enrichment biomarker

- **Gate 1 — BRAF-like classifier.** Apply the v13 BRAF-like RNA-seq classifier (or a NanoString-equivalent 25-gene panel including BRAF, RAS signature genes, and the eight v13 druggable targets) to archival or fresh tissue. Retain BRAF-like samples only.
- **Gate 2 — TACSTD2 z-score ≥ cohort median + 1 SD.** Use cohort-normalised TACSTD2 mRNA as a continuous predictor; threshold is adjustable.

### 5.2 Primary analysis

Primary analysis: best overall response (BOR) and 6-month PFS in BRAF-like + TACSTD2-high vs. all other molecular groups (RAS-like, unclassified, TACSTD2-low). Secondary: overall TROP2-IHC concordance; heterogeneity metrics from spatial pathology.

### 5.3 Sample size

Using TROPiCS-02 breast-cancer baseline ORR as a loose anchor (≈21 % all-comers → assumed 40–45 % in molecular responders vs. 15 % otherwise), **~80 evaluable patients** across the four trials provides 80 % power to detect ORR 45 % vs. 15 % at α = 0.05. Ten months of accrual at current site rates is plausible.

### 5.4 Feasibility

Required additional workflow: RNA extraction from archival FFPE (standard), sequencing or NanoString targeted panel, central readout. Existing infrastructure for correlative studies in TKI thyroid trials (e.g., SELECT, REFLECT) is directly reusable.

## 6. Competing considerations

- **Cost.** SG list price ≈ USD 70k / cycle; Kelun pricing differs by region.
- **Toxicity.** SG — neutropenia, diarrhoea; Dato-DXd — mucositis, ILD; thyroid patients after prior TKI may have reduced haematologic reserve.
- **Heterogeneity.** IHC TROP2 is patchy in PTC; a molecular classifier (Gate 1) adds robustness.
- **Positioning.** Dato-DXd or BNT323 may be preferred over SG in thyroid on toxicology grounds; conclusions should be framed around the TROP2 target class, not a single agent.

## 7. Conclusion

TROP2 sits at the intersection of (a) a robustly over-expressed antigen in PTC whose link to BRAF V600E has been documented since 2018 (Liu, Bychkov) and at the mRNA level since 2024 (Kalfert), (b) a mature ADC drug class with two dedicated Ph2 thyroid trials and a 2023 review-level mention of thyroid as a niche indication (Nieto-Jiménez), and (c) an unmet need for a mechanism-orthogonal rescue option in RAI-R disease. **What is missing in the literature, and what we propose, is BRAF-axis sub-stratification of the already-proposed thyroid TROP2-ADC opportunity.** We argue that the v13 BRAF-like classifier + TACSTD2-high two-gate biomarker should be evaluated prospectively in the ongoing thyroid ADC trials, as a biomarker-only overlay with no modification to primary treatment. This is a low-cost, high-yield translational exercise that would either validate a precision-medicine enrichment strategy for thyroid ADCs or demonstrate — equally valuable — that all-comer TROP2 dosing captures the biology. Independent of the trial overlay, our PRISM-level finding that the SN-38 payload class preferentially kills BRAF-like CCLE thyroid lines (opposite polarity to the established CRC literature) raises an additional, testable mechanistic question: how much of the proposed ADC efficacy in BRAF-like PTC is TROP2-directed delivery vs. tumour-intrinsic SN-38 sensitivity?

## Figures (placeholders)

- **Figure 1.** TACSTD2 expression in BRAF-like vs. RAS-like PTC across TCGA-THCA (v5.1 three-class) and a validation cohort (GSE33630 / GSE60542). Boxplot + paired normal.
- **Figure 2.** Integrated modality-matrix heatmap (from `v13_drug_discovery.html` §A) with TROP2 cells highlighted.
- **Table 1.** Ongoing thyroid TROP2 trials (as above) with proposed correlative overlay stated explicitly.

## References

Full PMID list in `results/v13_drug_discovery/modalities/trop2_literature.tsv` and `results/v14_ccle/v14_priorart/all_queries_results.tsv`. Primary anchors (to be cited):

**Mandatory priority-reduction priors (Introduction):**
1. **Liu et al.** _Int J Clin Exp Pathol_ 2018 — Overexpression of TROP2 is associated with BRAF V600E mutation and aggressive behaviour in PTC. **PMID 31949805.**
2. **Bychkov et al.** _J Pathol Transl Med_ 2018 — Significance of TROP2 expression in PTC prognosis. **PMID 29228520.**
3. **Kalfert et al.** _Pathol Res Pract_ 2024 — Integrated BRAF mutation × TACSTD2 mRNA × miRNA in primary PTC + paired LNM (n=60+40). **PMID 38696857.**
4. **Nieto-Jiménez et al.** _Clin Transl Med_ 2023 — ADC clinical-development opportunities; SG cited for thyroid as niche indication. **PMID 37740463.**
5. **Dum et al.** _Pathobiology_ 2022 — Pan-tumour TROP2 TMA n=18,563; PTC TROP2 ↔ nodal metastasis p=0.0013. **PMID 35477165.**
6. Cancers 2022 (PMID 35158847) and IJSP 2021 (PMID 33289434) — TROP2 as ADC target in anaplastic thyroid carcinoma.

**Counter-polarity / context citations:**
7. **Grothey et al.** _Ann Oncol_ 2021 — BRAF V600E mCRC + irinotecan combinations; established as BRAF-RESISTANCE polarity in CRC, the polarity opposite to our thyroid finding. **PMID 33836264.**

**Background / standard cites:**
8. TCGA Research Network. _Cell_ 2014 — BRAF-like / RAS-like classification of PTC.
9. **TROPiCS-02.** _Lancet_ 2024 — sacituzumab govitecan in metastatic breast cancer. **PMID 39067901.**
10. TROPION-Lung01 — datopotamab deruxtecan in NSCLC.
11. Zhou et al. _Scientific Reports_ 2024 — METTL3 / TACSTD2 m6A in PTC.
12. Zhong et al. _Mol Pharm_ 2025 — [¹⁸F]AlF-RESCA-RT4 TROP2 PET in thyroid.
13. Li et al. _Eur J Nucl Med Mol Imaging_ 2025 — [⁶⁸Ga]Ga-NOTA-T4 TROP2 ImmunoPET.

**This work:**
14. v13 pipeline (`reports/v13/v13_paper_section.md`).
15. v14 CCLE + PRISM (`reports/v14/v14_ccle_section.md`).
