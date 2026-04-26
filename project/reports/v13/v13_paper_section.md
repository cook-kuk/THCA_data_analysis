# v13 — Integrated multi-modality therapeutic landscape for BRAF-like thyroid cancer

_Proposed as Results §5.X in the Bioinformatics submission, to follow the v12 literature-validation section. ~1,100 words; tables/figures can be exported from `reports/html/pages/v13_drug_discovery.html`._

## 5.X Integrated multi-modality therapeutic landscape

### 5.X.1 Motivation

The v7 drug-repurposing atlas enumerated 83 drug–target pairs across our eight BRAF-like druggable targets but was restricted to small-molecule repurposing and to evidence already in DrugCentral / ChEMBL / OpenTargets. In this section we widen the search to **ten therapeutic modalities** — small molecule (SM), monoclonal antibody (mAb), antibody–drug conjugate (ADC), bispecific antibody, siRNA, PROTAC, CAR-T, radioligand, mRNA replacement, and ASO — and integrate the evidence with a **structure-informed virtual screen** and a **standardised ADMET filter**. The goal is to provide a practitioner-oriented "go-now / develop-next / novel-chemistry" triage for each target.

### 5.X.2 Pipeline overview

The v13 pipeline (Figure 5.X.1) consists of seven chained, idempotent Python scripts:

1. **Structure acquisition.** For each target UniProt ID we query RCSB PDB and fall back to AlphaFold DB v6. Crystal fragments (< 100 residues in the first model) are rejected in favour of the AF2 full-chain model. Five targets (CYP1B1, LDLR, TACSTD2, GABRB2, PTPRE) yielded crystal structures; three (TMPRSS4, PLEKHA6, B3GNT3) used AF2, all with mean pLDDT > 71.
2. **Pocket detection.** A lightweight α-sphere implementation (scipy KD-tree, alpha-sphere radius 3 Å, burial cutoff 8 Å) enumerates candidate pockets; known crystallographic ligand sites are extracted from HETATM records and merged. Each pocket receives a pocket-ligandability-proxy (PLB) score from hydrophobic fraction, aromatic fraction and α-sphere count.
3. **Druggability tiers.** Best-pocket PLB is combined with UniProt-derived subcellular location to assign each target a tier and a preferred-modality list.
4. **Virtual screen.** 2,443 FDA-approved small molecules are pulled from ChEMBL (molecule.json, max_phase=4, small-molecule filter). For each target we curate literature-anchored reference binders (e.g., α-naphthoflavone + flavonoids for CYP1B1; diazepam / etomidate / propofol for GABRB2; camostat / nafamostat for TMPRSS4; SN-38 / irinotecan for TACSTD2) and compute a Morgan/ECFP4 Tanimoto similarity for the full library.
5. **ADMET filter.** For the top-50 per-target hits (400 rows) we compute RDKit-native ADMET proxies: Lipinski / Veber rule passes, QED, ESOL logS, BBB-likelihood, hERG-risk heuristic, PAINS + Brenk substructural alerts and a small panel of hepatotox flags. A composite ADMET favourability is reported.
6. **Composite prioritisation.** score = 0.25·binding_norm + 0.20·admet_fav + 0.15·maturity + 0.15·off_patent + 0.10·thyroid_pubmed_norm + 0.10·modality_fit + 0.05·repurpose_bias. Tier A requires score ≥ 0.70 AND max_phase = 4 AND no red flag.
7. **Modality matrix.** 8 × 10 cells each carry an evidence class (approved / clinical / preclinical / hypothesis / none), an agent count and a representative agent; curated from biology plus live ClinicalTrials.gov v2 queries for named agents.

### 5.X.3 Results

**Pocket / druggability.** Five of eight targets score Tier-High (PLB ≥ 0.70): CYP1B1 (0.85, crystal ligand), GABRB2 (0.89, GABA site), TACSTD2 (0.83, crystal ligand), LDLR (0.70, ligand-adjacent pocket) and PLEKHA6 (0.81, hydrophobic scaffold pocket). TMPRSS4 (0.51) and B3GNT3 (0.58) are Tier-Medium; PTPRE (0.39) is Tier-Low, consistent with the known historical difficulty of phosphatase small-molecule campaigns.

**Virtual screen.** Of 400 per-target hits, the screen recapitulates the expected on-target chemistry wherever reference binders exist — diazepam / propofol / etomidate dominate GABRB2; rosuvastatin / simvastatin / ezetimibe dominate LDLR; resveratrol / luteolin / cannabidiol dominate CYP1B1. The screen also surfaces **topotecan as the closest FDA-approved surrogate for the SN-38 payload of sacituzumab govitecan** (max Tanimoto 0.69 against TACSTD2 references), and **CRBN-binders (thalidomide / lenalidomide / pomalidomide)** as the top match to PLEKHA6 reference chemistry — not as direct ligands but as obvious PROTAC-warhead starting points (Section 5.X.4).

**ADMET.** 216/400 hits carry at least one red flag (PAINS, Brenk, hERG > 0.6, or hepatotox substructure), leaving **184 clean** candidates. GABRB2 and PLEKHA6 have the highest mean ADMET favourability (0.77); TACSTD2 and TMPRSS4 the lowest (0.72–0.74), reflecting the lipophilic / larger-MW nature of their reference chemistries.

**Composite tiers.** Of the 400 hits, 59 are Tier-A, 283 Tier-B, 56 Tier-C and 2 Tier-D. Notable Tier-A leaders: propofol (GABRB2, 0.963), simvastatin (LDLR, 0.936), etomidate (GABRB2, 0.922), dronabinol (CYP1B1, 0.764), indomethacin (PTPRE, 0.726), topotecan (TACSTD2, 0.696 — borderline Tier A/B).

**Modality matrix (highlights).** The 80-cell matrix contains 9 approved cells, 11 clinical cells, 24 preclinical, 9 hypothesis and 27 none (Figure 5.X.2). Four cells span the approved–clinical divide and define our priority programmes:

- **TACSTD2 × ADC** — approved sacituzumab govitecan (breast, urothelial), clinical sacituzumab tirumotecan and datopotamab deruxtecan; **two active Ph2 thyroid trials** (NCT06235216, NCT07521670).
- **LDLR × siRNA** — approved inclisiran (PCSK9 siRNA, raises LDLR indirectly); clinical lepodisiran.
- **GABRB2 × radioligand** — [¹¹C]flumazenil PET approved (benzodiazepine site) — diagnostic not therapeutic.
- **CYP1B1 × SM** — resveratrol and quercetin in multiple Ph2/Ph3 studies; dronabinol from the Tier-A screen.

### 5.X.4 Novel-modality hypotheses

For targets where no approved / clinical small molecule exists we identify a preferred emerging modality from the matrix + screen output:

- **PLEKHA6 (PH-domain scaffold) → PROTAC.** IMiD analogues (thalidomide, lenalidomide, pomalidomide) occupy the top three similarity-screen hits. Their role is as _warhead chemistry_ rather than direct binders — providing cereblon-recruiting fragments for a de-novo PLEKHA6 degrader campaign. Designed PROTAC chemistry is a natural v15 (REINVENT4) output.
- **PTPRE (receptor PTP) → PROTAC / siRNA.** Low PLB (0.39) plus non-selective pan-PTP chemistry argue against a classical SM campaign. Phosphatase PROTACs (cf. SHP2 degraders) are a growing space; siRNA is already preclinical.
- **B3GNT3 (Golgi glycosyltransferase) → siRNA / ASO.** The *Cell* 2018 B3GNT3-knockdown / anti-PD-1 synergy result provides strong biological rationale; the intra-Golgi topology makes antibody modalities unsuitable.
- **TMPRSS4 (plasma-membrane serine protease) → siRNA + targeted mAb.** Nafamostat cross-reactivity is non-selective; a specific mAb (Kim 2020) validates the surface epitope.

### 5.X.5 Flagship finding — BRAF sub-stratification of TROP2 ADC

The **TACSTD2 × ADC** cell is the most actionable output of the pipeline, but the underlying biology is not novel: the BRAF V600E ↔ TROP2 link in PTC has been documented since 2018 (Liu, _Int J Clin Exp Pathol_, PMID 31949805; Bychkov, _J Pathol Transl Med_, PMID 29228520) and extended to the mRNA axis in primary PTC + paired LNM by Kalfert et al. (2024, _Pathol Res Pract_, PMID 38696857). Sacituzumab govitecan (SG) for thyroid was previously catalogued among "non-explored niche indications" by Nieto-Jiménez et al. (2023, _Clin Transl Med_, PMID 37740463). Our contribution is therefore not the discovery of TROP2 as a thyroid ADC target — that has been done — but a **molecular sub-stratification** of the previously-proposed thyroid ADC opportunity: two dedicated thyroid Ph2 trials (NCT06235216 "Sacituzumab govitEcan in THYroid Cancers", recruiting; NCT07521670 "STRAP", adenoid-cystic + PTC) are already under way, but **neither stratifies patients by BRAF status, and STRAP explicitly waives even TROP2 IHC testing**. We propose a correlative overlay using the v13 BRAF-like classifier plus TACSTD2 mRNA z-score as a two-gate enrichment biomarker. A full rationale and proposed correlative design — including a separate, novel pharmacogenomic finding (PRISM-level BRAF-selective SN-38 sensitivity in CCLE thyroid lines, opposite polarity to the established CRC literature) — are in `v13_trop2_standalone_draft.md` and `reports/v14/v14_ccle_section.md`, with the trial-overlay framing intended for a short-form _JCO Precision Oncology_ letter.

### 5.X.6 Limitations

- **Similarity-based screen.** Docking (GNINA / Vina) was not run due to environment constraints; Tanimoto surrogates over-weight chemotypes near the reference set and recover known binders well but can miss novel scaffolds. A GNINA re-run over the Tier-B shortlist is scheduled for v14.
- **ADMET proxies.** Heuristic hERG and hepatotox scoring are not substitute for deep-learning ADMET (ADMET-AI, Chemprop). Proportional agreement on flag-vs-no-flag ranking is expected; absolute numbers should not be interpreted as clinical predictions.
- **AlphaFold2 fallbacks.** Three targets rely on AF2 (TMPRSS4, PLEKHA6, B3GNT3). AlphaFold3-quality models via 유 교수님 access are noted for a v14 structure-refinement pass.

### 5.X.7 Data availability

All raw tables and the interactive dashboard are published at `/pages/v13_drug_discovery.html`, with a pipeline-level reproducibility footnote in the same page's §G. Per-target pocket TSVs, per-target top-50 screen, 400-row prioritised candidate table and the modality-matrix JSON are under `results/v13_drug_discovery/` and mirrored to `reports/html/data/v13/` for direct browser download.

_Files for paper production: `reports/v13/v13_paper_section.md` (this file); figure export: `reports/html/pages/v13_drug_discovery.html` (PNG capture of modality heatmap + screen scatter)._
