# High-IF inspired THCA paper strategy - 2026-05-08

User instruction: do not chase low-impact thyroid papers; inspect very high-impact journal logic and derive non-copycat THCA directions.

## What I looked at

I used these papers as design references, not as things to copy.

| High-impact template | Journal/date | What matters for us | URL |
|---|---|---|---|
| Non-invasive profiling of the tumour microenvironment with spatial ecotypes | Nature, 2026 | The current top-tier language is "spatial ecotypes": multicellular ecosystems, conserved across cancers, linked to clinical state and even recoverable from another modality. | https://www.nature.com/articles/s41586-026-10452-4 |
| Tumour evolution and microenvironment interactions in 2D and 3D space | Nature, 2024 | The biology is not just spatial clustering. It is clonal/subclonal evolution plus local immune/stromal neighbourhoods, ideally across sections or 3D-like structure. | https://www.nature.com/articles/s41586-024-08087-4 |
| Inferring allele-specific copy number aberrations and tumor phylogeography from spatial transcriptomics | Nature Methods, 2024 | Spatial RNA can be pushed toward copy-number clone maps and phylogeography. This directly upgrades our T00 CNV finding. | https://www.nature.com/articles/s41592-024-02438-9 |
| A multimodal whole-slide foundation model for pathology | Nature Medicine, 2025 | High-impact pathology AI now values vision-language alignment, rare disease retrieval, zero/few-shot, and patient/slide-level clinical utility. | https://www.nature.com/articles/s41591-025-03982-3 |
| AI-enabled virtual spatial proteomics from histopathology | Nature Medicine, 2026 | The strong version of H&E AI is not "predict one marker"; it generates/interprets virtual spatial molecular layers and validates clinical outcome links. | https://www.nature.com/articles/s41591-025-04060-4 |
| A visual-omics foundation model to bridge histopathology with spatial transcriptomics | Nature Methods, 2025 | H&E and spatial transcriptomics can be aligned by retrieval/embedding, but the compelling claim is evidence-grounded molecular retrieval. | https://www.nature.com/articles/s41592-025-02707-1 |
| Novae: a graph-based foundation model for spatial transcriptomics data | Nature Methods, 2025 | Spatial domains should be modeled as graph/context units, not independent spots. Good inspiration for thyroid spatial domain hierarchy. | https://www.nature.com/articles/s41592-025-02899-6 |
| Conserved spatial subtypes and cellular neighborhoods of CAFs | Cancer Cell, 2025 | CAF/ECM is high-impact only when framed as conserved spatial subtype/neighborhood biology, not a generic stromal score. | https://doi.org/10.1016/j.ccell.2025.03.004 |

## High-IF pattern extracted

Top-tier papers are currently built around one of five moves:

1. Define a new unit of tumor biology: spatial ecotype, clone territory, microregion, multicellular module.
2. Show that the unit is not a marker artifact: validated across modalities, datasets, or perturbation/clinical context.
3. Link malignant-intrinsic state to ecosystem state: genome/CNV/lineage state organizes immune, stromal, and differentiation neighborhoods.
4. Use multimodal AI as evidence retrieval, not decoration: every prediction should map back to observed tissue/spatial data.
5. Provide a clear clinical boundary: what can be inferred non-invasively/from H&E, and what remains molecularly hidden.

## Best non-copycat direction

### H1. Copy-number-defined residual thyroid cancers organize spatial molecular ecosystems

This is the strongest direction.

Proposed title:

> Copy-number states organize spatial molecular ecosystems in driver-negative thyroid cancer

Core claim:

> True driver-negative thyroid cancers are not a leftover bin. A recurrent arm-level CNV class, enriched in DM2/FVPTC-like tumors, defines a residual malignant state that organizes thyroid differentiation, TROP2, CAF/ECM, TLS/immune, and RAI-related spatial ecosystems.

Why this is not copying:

- Nature HTAN/CalicoST papers give the design grammar: spatial tumor evolution plus local ecosystem.
- Our biological object is thyroid-specific: true driver-negative THCA, Class6-DM2, 7p/7q/12q/16p/2p/2q/16q arm-CNV signature.
- Our pilot result is already strong: cBioPortal ARMDRIVER signature, DM2 35.4% vs DM1 1.3%, OR 41.1, p=1.6e-7.
- Robustness already survived cPTC, FVPTC, purity high/low, Stage I, and Ras-like strata.

High-IF version of the paper:

1. Bulk genome layer: define CNV-residual class in TCGA/cBio/GDC.
2. Spatial layer: map TROP2, CAF/ECM, TLS, RAI/DM axes across Visium/GeoMx.
3. Spatial-CNV bridge: infer arm-level expression-CNV or CalicoST-like clone territories in spatial slides.
4. Multimodal layer: show which spatial ecosystem states are recoverable by morphology + clinical text, and which are not.
5. Clinical boundary: driver-negative does not mean genomically quiet; a subset is CNV-ecosystem organized.

Current grade: A-.

Upgrade gate to A/A+:

- Run spatial-CNV clone inference on Visium, even if first-pass is expression-CNV by chromosome arm.
- Add at least one orthogonal validation: GeoMx compartment, IHC/IF for TROP2 + CAF/TLS, or independent ST cohort.
- Keep RAI-refractory claims out unless treated outcome data exist.

## Other high-IF derived ideas

| Rank | Topic | High-IF hook | Local support | Risk | Verdict |
|---:|---|---|---|---|---|
| 1 | H1 CNV-residual spatial ecosystem | Genome evolution + spatial ecotype + thyroid residual class | cBio ARMDRIVER p=1.6e-7; stratified robustness; spatial axes exist | Needs spatial-CNV bridge | Flagship |
| 2 | H2 Spatial phylogeography of driver-negative thyroid cancer | CalicoST-style CNA clone geography from Visium | Arm-CNV signature is strong; Visium data available | Requires method execution and careful validation | Do next |
| 3 | H3 Morphology-visible vs context-recoverable vs molecular-hidden thyroid programs | TITAN/Haiku/HEX-style multimodal boundary map | TROP2 fusion retrieval rho 0.532/AUROC 0.825; DM/RAI remains weak | Needs real WSI embeddings for top-tier | Companion |
| 4 | H4 Thyroid spatial ecotypes from bulk-to-spatial transfer | Nature spatial ecotype framing applied to THCA | TROP2, CAF/ECM, TLS, immune, RAI modules already scored | Needs robust ecotype definitions, not hand-picked scores | Strong if integrated |
| 5 | H5 CAF/ECM-instructed thyroid lineage suppression | Cancer Cell CAF neighborhood grammar | DCN/COL1A2/SFRP2/COL1A1 spatial ligands; CAF/M1_M2 vs RAI negative | Stromal abundance confounding | Mechanism layer |
| 6 | H6 TROP2 as a spatial niche, not a bulk biomarker | Virtual spatial proteomics / spatial marker rediscovery | TROP2 Moran d=4.99; retrieval AUROC 0.825 | Standalone TROP2 is crowded and ADC-overclaim-prone | Sub-aim |
| 7 | H7 Immune-rich protective thyroid ecotype | Spatial ecotype / immune geography | DSS M1_M2 HR 0.347 p=0.005; CD8 HR 0.277 p=0.011 | Hashimoto vs antitumor immunity separation | Sub-aim |
| 8 | H8 ATC transition as ecosystem rewiring | Cancer progression rewires multicellular modules | ATC myeloid/checkpoint signals; GATCI advanced-CNA context | n=4 spatial ATC | Exploratory |
| 9 | H9 DICER1/DGCR8 compartment progression route | Rare genotype + compartment-resolved progression | GeoMx PDTC vs microPTC lineage/HLA shifts | genotype-histology confounding | Secondary paper |
| 10 | H10 NIS/SLC5A5 non-methylation escape | Mechanism exception inside RAI-lineage silencing | SLC5A5 methylation null while TPO/TG/TSHR/PAX8 strong | Needs histone/TF/miRNA data | Mechanistic follow-up |
| 11 | H11 TERT-only microclass | Residual risk class | TERT log-rank strong in recovered data | TERT-only n too small; known biology | Do not lead |
| 12 | H12 GLS/metabolic vulnerability | Functional dependency angle | Some pan-cancer dependency signal | thyroid model coverage weak; wet-lab required | Hold |

## What not to do

Do not lead with these:

- "TROP2 ADC in thyroid cancer" as the main paper. Too crowded and too translational without functional validation.
- "H&E predicts RAI/DM." Our closure tests say this is weak/negative.
- "TERT-only triple-negative class." Interesting but too small.
- "HLA cancer association." Prior audit marked this as contaminated/claim-risky.
- "GLS synthetic lethality." Needs direct thyroid wet-lab or stronger model coverage.
- "Spatial ligand-receptor in ATC" alone. Sample-limited and too easy to look like a descriptive atlas.

## Recommended paper architecture

Figure 1. Define the true-driver-negative residual cohort and Class6-DM2 CNV architecture.

Figure 2. Validate the CNV signature in cBio/GDC and across histology, purity, stage, and BRAF/RAS-like strata.

Figure 3. Map spatial ecosystems: TROP2 niche, CAF/ECM, TLS/immune, RAI/DM, thyroid differentiation.

Figure 4. Spatial-CNV bridge: chromosome-arm expression-CNV or CalicoST-like clone territories overlaid with ecosystems.

Figure 5. Multimodal boundary: TROP2 is context-recoverable by morphology + metadata retrieval; DM/RAI remains molecular-hidden.

Figure 6. Clinical interpretation and validation plan: driver-negative thyroid cancer contains a CNV-organized residual state, not a homogeneous quiet group.

## Final call

The paper with the best high-IF ceiling is not a new standalone marker story. It is:

> **A spatial-genomic ecosystem paper for true driver-negative thyroid cancer.**

Use T00 as the spine. Use T02/T04/T07/T12 as ecosystem and multimodal layers. Use T10 as the honest negative boundary. This is the only direction from the current evidence that resembles the ambition of Nature/Cancer Cell/Nature Medicine logic without copying any single paper.

