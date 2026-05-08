# Haiku arXiv 2605.00925 review and THCA multimodal direction

Date: 2026-05-08

## 1. Paper reviewed

**Cui et al., "Linking spatial biology and clinical histology via Haiku"**

- arXiv: https://arxiv.org/abs/2605.00925
- submitted: 2026-04-30
- code: https://github.com/zhihuanglab/Haiku
- model/demo data: gated Hugging Face repos (`zhihuanglab/Haiku`, `zhihuanglab/Haiku-demo-data`)

## 2. What Haiku actually contributes

Haiku is not just H&E biomarker prediction. It is a tri-modal contrastive framework aligning:

- H&E histology patches
- multiplexed immunofluorescence / CODEX spatial proteomics patches
- structured clinical metadata and patch-level text descriptions

Key scale and results reported by the authors:

- 26.7M paired spatial proteomics patches from 3,218 tissue sections, 1,606 patients, 11 organ types.
- Full mIF corpus: 7,600 slices; held-out evaluation includes 336 paired and 198 mIF-only slices.
- Patient-level train/test split.
- Cross-modal retrieval: H&E-to-mIF Recall@50 0.611, mIF-to-H&E 0.604, Text-to-mIF 0.169.
- Fusion retrieval adds metadata text to H&E and improves biomarker inference; mean PCC 0.718 across 52 biomarkers.
- Survival prediction C-index 0.737.
- Counterfactual analysis changes only clinical metadata while fixing H&E morphology, then inspects shifts in retrieved mIF patches.

The strongest conceptual point for us is **metadata-conditioned retrieval**, not the raw model scale. Haiku says: morphology alone is incomplete, clinical context provides a semantic prior, and the predicted molecular state should be grounded in retrieved real spatial measurements rather than generated pixels.

## 3. Limitations relevant to our use

- Haiku is arXiv-only as of 2026-05-08.
- The released model/data are gated and large: model about 3.2 GB, demo data about 3.5 GB, GPU recommended.
- Counterfactual examples are explicitly hypothesis-generating and single-patient proof-of-concepts.
- It uses mIF/CODEX proteomics, while our strongest public thyroid layer is Visium spatial RNA plus GeoMx ROI. Direct replication is not the right goal.

## 4. THCA-specific multimodal pilot run

I ran a Haiku-lite pilot using our existing spot-level THCA spatial/morphology table.

Input:

- `project/results/03_pathology_poc/spark_st_joint_per_spot.tsv.gz`
- `project/results/spatial_full_2026_05_06/spatial_per_spot_TROP2.tsv.gz`
- 3,200 Visium-anchored spots from 16 slides
- modalities: local morphology/cell-composition features + stage-as-clinical-text proxy + spatial RNA targets
- validation: leave-one-slide-out ridge models

Output:

- `project/results/high_impact_topic_pilots_2026_05_08/multimodal_haiku_lite/multimodal_haiku_lite_report.md`
- `project/results/high_impact_topic_pilots_2026_05_08/multimodal_haiku_lite/multimodal_haiku_lite_metrics.tsv`
- `project/results/high_impact_topic_pilots_2026_05_08/multimodal_haiku_lite/multimodal_haiku_lite_delta.tsv`
- `project/notebooks_or_scripts/multimodal_haiku_lite_pilot_2026_05_08.py`

Headline results:

| target | morphology-only rho | clinical-text-only rho | morphology+text rho | morphology+text AUROC q75 | interpretation |
|---|---:|---:|---:|---:|---|
| TROP2 | -0.044 | 0.372 | **0.564** | **0.844** | strongest Haiku-like signal |
| DM1_like | 0.189 | -0.104 | 0.237 | 0.693 | weak/moderate; still not H&E-recoverable |
| RAI_8 | 0.189 | -0.104 | 0.237 | 0.568 | weak/moderate |
| TDS_like | 0.163 | -0.105 | 0.196 | 0.551 | weak |
| Epithelial | 0.160 | -0.018 | 0.181 | 0.535 | weak |
| CAF_ECM | 0.104 | -0.051 | 0.135 | 0.607 | weak/moderate q75 |
| EMT | 0.037 | -0.044 | 0.038 | 0.473 | negative |

Readout:

- TROP2 is the only axis where the Haiku-like fusion idea clearly works in this pilot.
- Morphology alone does not recover TROP2, but clinical context plus morphology does.
- DM1/RAI stay weak, matching our earlier closure battery. This is useful: it means the multimodal story should **not** claim routine H&E can solve the differentiation axis.

## 5. New multimodal paper topic

**T12. Thyroid-Haiku-lite: clinical-text conditioned retrieval of spatial molecular niches in thyroid cancer**

Core thesis:

> In thyroid cancer, spatial TROP2 is morphology-invisible but clinically/contextually recoverable, whereas the DM/RAI lineage axis remains only weakly recoverable from morphology and metadata, defining a boundary between image-visible and molecularly hidden thyroid cancer programs.

Best version:

- Build a thyroid-specific retrieval atlas, not a large foundation model.
- Query modalities: H&E/spot morphology, clinical text metadata, CNV/driver context where available.
- Reference modalities: Visium spatial RNA niches, GeoMx ROI programs, TCGA/cBio/GATCI molecular states.
- Use evidence-grounded retrieval: every inferred niche should map back to an observed Visium/GeoMx region.
- Counterfactual prompts: keep morphology fixed, change metadata text such as PTC -> LPTC -> ATC, HT-overlap yes/no, driver-negative CNV class yes/no, and inspect retrieved spatial programs.

Current grade: **B+**

Upgrade gate to A-:

- Replace simple ridge pilot with actual retrieval/contrastive embedding, even if small.
- Add full-resolution H&E or real WSI embeddings for the 16 GSE250521 slides.
- Use TROP2, TLS, DM1/RAI, CAF/ECM as multi-target retrieval axes.
- Keep all claims evidence-grounded; no synthetic molecular maps.

## 6. Updated topic ranking after Haiku/multimodal check

| rank | topic | grade now | why |
|---:|---|---|---|
| 1 | T00 CNV-defined true-driver-negative residual class | A- | cBioPortal ARMDRIVER_CN independently supports Class6-DM2 enrichment; strongest biology flagship. |
| 2 | T02 TROP2 tumor-specific spatial niche | A- | strongest spatial signal; Haiku-lite pilot makes it multimodal-relevant. |
| 3 | T12 Thyroid-Haiku-lite multimodal retrieval atlas | B+ | interesting and timely, but needs actual retrieval/contrastive step before A-. |
| 4 | T01 DICER1/DGCR8 GeoMx progression route | B+ | good GeoMx signal, but genotype/histology confounding. |
| 5 | T04 CAF/ECM ligand remodeling | B+ | useful as T12 reference axis or T00 mechanism layer. |
| 6 | T07 immune/TLS protective paradox | B+ | clinically interesting, needs external survival validation. |
| 7 | T03 Hashimoto-overlap TLS niche | B | good spatial organization, n=4 limitation. |
| 8 | T06 NIS/SLC5A5 non-methylation silencing | B | mechanistic follow-up, not standalone yet. |
| 9 | T08 DICER1/EIF1AX rare-driver residual class | B- | n too small unless pooled externally. |
| 10 | T05 checkpoint ligand-receptor routes | B- | sample-limited and JCI Insight overlap risk. |
| 11 | T11 ATC spatial coherence collapse | B | improved by GATCI CNA context, still n=4 spatial. |
| 12 | T10 morphology-invisible DM/RAI axis | negative-use-only | valuable as a boundary/negative-control module inside T12. |
| 13 | T09 TERT-only triple-negative microclass | C | too small. |

## 7. Recommendation

The best high-IF strategy is **not** to abandon T00 for multimodal. The cleanest plan is:

1. Lead with **T00 CNV-defined true-driver-negative residual class** as the biology flagship.
2. Develop **T02 + T12** as the high-novelty spatial/multimodal companion: TROP2 niche plus Haiku-inspired evidence-grounded retrieval.
3. Use T10 as an honest negative: DM/RAI is not sufficiently visible from routine morphology, so direct spatial/molecular measurement remains necessary.

If we want one "fun" multimodal paper direction, choose **T12**, but pitch it as a thyroid-specific evidence-grounded retrieval atlas rather than a foundation model.
