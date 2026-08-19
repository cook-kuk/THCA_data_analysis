#!/usr/bin/env python3
"""Build a Paper 2 reference leverage package.

This is scaffold/figure/table work only. It does not write protected manuscript
voice sections.
"""
from __future__ import annotations

import html
import json
import shutil
from pathlib import Path

import cv2
import numpy as np
import pandas as pd


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
P2 = ROOT / "project/results/p2_image_dm1_v2_foundation_clam_2026_05_07"
SUPP = P2 / "analysis_supp"
OUT = SUPP / "paper2_reference_leverage_2026_05_10"
HUB = ROOT / "project/papers_hub_2026_05_04"
LIVE = Path("/var/www/papers/papers_hub_2026_05_04")
ASSET = "paper2_reference_leverage"
ASSET_LOCAL = HUB / "assets" / ASSET
ASSET_LIVE = LIVE / "assets" / ASSET
PAGE_LOCAL = HUB / "paper2_reference_leverage_pack.html"
PAGE_LIVE = LIVE / "paper2_reference_leverage_pack.html"

BG = (13, 17, 23)
PANEL = (23, 30, 41)
PANEL2 = (31, 40, 54)
LINE = (68, 79, 94)
INK = (235, 240, 246)
MUTED = (156, 169, 188)
GOLD = (94, 178, 214)
TEAL = (196, 213, 99)
BLUE = (255, 166, 87)
RED = (114, 123, 255)
GREEN = (157, 211, 53)
PURPLE = (219, 139, 240)
WHITE = (245, 246, 248)

CLAIM_LAYERS = [
    ("he_molecular", "H&E -> molecular state"),
    ("he_spatial_rna", "H&E -> spatial RNA"),
    ("virtual_spatial", "virtual spatial omics"),
    ("foundation", "foundation pathology"),
    ("weak_supervision", "WSI weak supervision"),
    ("ffpe_spatial", "FFPE/spatial tissue"),
    ("thyroid_biology", "thyroid/RAI biology"),
    ("guardrail", "benchmark/guardrail"),
]


def read_json(path: Path) -> dict:
    with path.open() as handle:
        return json.load(handle)


def metrics() -> dict:
    return {
        "uni": read_json(SUPP / "audit_uni_loto/UNI_LOTO_SUMMARY.json"),
        "g250_stage": read_json(
            SUPP / "path2space_stage_generalization_controls_2026_05_09/GSE250521_STAGE_GENERALIZATION_SUMMARY.json"
        ),
        "g250_spec": read_json(
            SUPP / "path2space_gse250521_random_module_specificity_2026_05_09/GSE250521_RANDOM_MODULE_SPECIFICITY_SUMMARY.json"
        ),
        "g230": read_json(SUPP / "gse230424_pathology_thyroid_axis_2026_05_09/GSE230424_PATHOLOGY_THYROID_SUMMARY.json"),
        "g230_resid": read_json(
            SUPP
            / "gse230424_pathology_thyroid_axis_2026_05_09/residual_target_controls/GSE230424_RESIDUAL_TARGET_SUMMARY.json"
        ),
        "decile": read_json(SUPP / "path2space_decile_dose_response_2026_05_09/PATH2SPACE_DECILE_DOSE_RESPONSE_SUMMARY.json"),
        "hotspot": read_json(SUPP / "path2space_hotspot_concordance_2026_05_09/PATH2SPACE_HOTSPOT_CONCORDANCE_SUMMARY.json"),
        "aitd": read_json(SUPP / "gse248205_pathology_aitd_axis_2026_05_09/GSE248205_PATHOLOGY_AITD_SUMMARY.json"),
    }


def put(img, text, xy, scale=0.55, color=INK, thick=1, width=None, gap=8):
    x, y = xy
    font = cv2.FONT_HERSHEY_SIMPLEX
    lines = []
    for raw in str(text).split("\n"):
        if width is None:
            lines.append(raw)
            continue
        line = ""
        for word in raw.split():
            trial = word if not line else f"{line} {word}"
            if cv2.getTextSize(trial, font, scale, thick)[0][0] <= width or not line:
                line = trial
            else:
                lines.append(line)
                line = word
        lines.append(line)
    step = cv2.getTextSize("Ag", font, scale, thick)[0][1] + gap
    for i, line in enumerate(lines):
        cv2.putText(img, line, (x, y + i * step), font, scale, color, thick, cv2.LINE_AA)
    return y + max(len(lines), 1) * step


def box(img, xyxy, color=PANEL, outline=LINE):
    x0, y0, x1, y1 = xyxy
    cv2.rectangle(img, (x0, y0), (x1, y1), color, -1)
    cv2.rectangle(img, (x0, y0), (x1, y1), outline, 1, cv2.LINE_AA)


def title(img, kicker, heading, sub):
    put(img, kicker.upper(), (70, 78), 0.58, GOLD, 2)
    put(img, heading, (70, 142), 1.1, INK, 2, 1700)
    put(img, sub, (70, 205), 0.56, MUTED, 1, 1560)


def split_tags(tags: str) -> set[str]:
    return {tag.strip() for tag in tags.split(",") if tag.strip()}


def make_reference_table() -> pd.DataFrame:
    rows = [
        [
            "R01",
            "H&E-to-spatial transcriptomics",
            "Cell 2026",
            "AI-predicted spatial transcriptomics unlocks breast cancer biomarkers from pathology",
            "10.1016/j.cell.2026.04.023",
            "https://www.cell.com/cell/abstract/S0092-8674(26)00458-7",
            "he_spatial_rna,virtual_spatial,guardrail",
            "This is the closest high-impact precedent: H&E-predicted spatial transcriptomics enabling biomarker discovery.",
            "Foreground in Fig 1 and presub; state that Paper 2 applies the logic to thyroid DM1/RAI.",
            "Breast cancer precedent, not thyroid; do not claim same model or same biomarker class.",
            "critical",
            "local PDF plus web DOI check",
        ],
        [
            "R02",
            "H&E-to-spatial transcriptomics",
            "Nat Biomed Eng 2020",
            "Integrating spatial gene expression and breast tumour morphology via deep learning",
            "10.1038/s41551-020-0578-x",
            "https://www.nature.com/articles/s41551-020-0578-x",
            "he_spatial_rna,virtual_spatial,guardrail",
            "ST-Net established local gene-expression prediction from H&E as a serious biomedical-engineering problem.",
            "Use as the classic ST prediction precedent.",
            "Breast cancer and paired ST training; not a thyroid validation.",
            "critical",
            "publisher",
        ],
        [
            "R03",
            "H&E-to-transcriptomics",
            "Nat Commun 2020",
            "A deep learning model to predict RNA-Seq expression of tumours from whole slide images",
            "10.1038/s41467-020-17678-4",
            "https://www.nature.com/articles/s41467-020-17678-4",
            "he_molecular,virtual_spatial,weak_supervision",
            "HE2RNA supports the broad claim that WSI morphology contains transcriptomic signal.",
            "Use for the general H&E-to-RNA premise and spatialization precedent.",
            "Bulk RNA labels, not true thyroid spatial RNA.",
            "critical",
            "publisher",
        ],
        [
            "R04",
            "Benchmark/guardrail",
            "Nat Commun 2025",
            "Benchmarking the translational potential of spatial gene expression prediction from histology",
            "10.1038/s41467-025-56618-y",
            "https://www.nature.com/articles/s41467-025-56618-y",
            "he_spatial_rna,guardrail",
            "Provides editor-facing precedent that the field now judges H&E-to-spatial expression by benchmarking and translational guardrails.",
            "Use to justify honest caveats, held-out testing, and random-module controls.",
            "Benchmark paper; does not validate DM1 biology.",
            "critical",
            "publisher",
        ],
        [
            "R05",
            "Robust H&E-to-ST",
            "Nat Commun 2026",
            "Robust and interpretable prediction of gene markers and cell types from spatial transcriptomics data",
            "10.1038/s41467-026-68487-0",
            "https://www.nature.com/articles/s41467-026-68487-0",
            "he_spatial_rna,foundation,guardrail",
            "STimage makes robustness, uncertainty, interpretability, and multi-platform validation part of the current field standard.",
            "Use to explain why Paper 2 co-locates caveats with results.",
            "Do not imply Paper 2 performs STimage-level uncertainty modeling.",
            "high",
            "publisher",
        ],
        [
            "R06",
            "Single-cell H&E-to-expression",
            "Nat Commun 2026",
            "sCellST predicts single-cell gene expression from H&E images",
            "10.1038/s41467-025-67965-1",
            "https://www.nature.com/articles/s41467-025-67965-1",
            "he_spatial_rna,virtual_spatial,guardrail",
            "Moves the field from spot-scale prediction toward cell-level molecular interpretation from H&E.",
            "Use as a frontier reference; Paper 2 stays spot/module-level.",
            "Not thyroid; not used by our pipeline.",
            "high",
            "publisher",
        ],
        [
            "R07",
            "ST-pathology foundation model",
            "npj Digit Med 2025",
            "STPath: a generative foundation model for integrating spatial transcriptomics and whole-slide images",
            "10.1038/s41746-025-02020-3",
            "https://www.nature.com/articles/s41746-025-02020-3",
            "he_spatial_rna,foundation,virtual_spatial,guardrail",
            "Shows that ST+WSI foundation modeling is now an explicit digital-medicine direction.",
            "Use to position Paper 2 as a disease-focused biomarker application rather than a new universal model.",
            "Do not claim our model has STPath-scale pretraining or organ breadth.",
            "high",
            "publisher",
        ],
        [
            "R08",
            "Super-resolution ST + histology",
            "Nat Biotechnol 2024",
            "Inferring super-resolution tissue architecture by integrating spatial transcriptomics with histology",
            "10.1038/s41587-023-02019-9",
            "https://www.nature.com/articles/s41587-023-02019-9",
            "he_spatial_rna,virtual_spatial",
            "iStar supports the spatial-resolution argument: histology can improve molecular maps beyond spot resolution.",
            "Use as adjacent method precedent.",
            "It integrates measured ST; Paper 2 infers thyroid DM1 spatial state from H&E features.",
            "medium",
            "publisher",
        ],
        [
            "R09",
            "Multi-tissue H&E-to-ST",
            "Nat Commun 2024",
            "Self-supervised learning for characterising histomorphological diversity and spatial RNA expression prediction across 23 human tissue types",
            "10.1038/s41467-024-50317-w",
            "https://www.nature.com/articles/s41467-024-50317-w",
            "he_spatial_rna,foundation,guardrail",
            "Supports the general tissue-diversity framing for self-supervised histology representations.",
            "Use to keep H&E framing broad without claiming pan-disease validation.",
            "General tissue paper; not thyroid-specific.",
            "medium",
            "publisher",
        ],
        [
            "R10",
            "Virtual spatial proteomics",
            "Nat Med 2026",
            "AI-enabled virtual spatial proteomics from histopathology for interpretable biomarker discovery in lung cancer",
            "10.1038/s41591-025-04060-4",
            "https://www.nature.com/articles/s41591-025-04060-4",
            "he_molecular,virtual_spatial,foundation,guardrail",
            "HEX is a major clinical precedent that H&E can generate virtual spatial molecular maps with outcome relevance.",
            "Use as the strongest translational analog; Paper 2 is RNA/thyroid, not proteomics/lung.",
            "Do not claim equal clinical validation or immunotherapy response prediction.",
            "critical",
            "publisher",
        ],
        [
            "R11",
            "H&E mutation prediction",
            "Nat Med 2018",
            "Classification and mutation prediction from non-small cell lung cancer histopathology images using deep learning",
            "10.1038/s41591-018-0177-5",
            "https://www.nature.com/articles/s41591-018-0177-5",
            "he_molecular,weak_supervision,guardrail",
            "Classic proof that H&E can predict molecular alterations from routine slides.",
            "Use as historical root of image-to-molecular biomarkers.",
            "Mutation prediction is not transcriptomic/spatial prediction.",
            "high",
            "publisher",
        ],
        [
            "R12",
            "Pan-cancer molecular pathology",
            "Nat Cancer 2020",
            "Pan-cancer image-based detection of clinically actionable genetic alterations",
            "10.1038/s43018-020-0087-6",
            "https://www.nature.com/articles/s43018-020-0087-6",
            "he_molecular,weak_supervision,guardrail",
            "Shows that routine FFPE H&E can reflect a broad range of molecular alterations across cancers.",
            "Use for generality of molecular information in histology.",
            "Do not imply Paper 2 is pan-cancer.",
            "high",
            "publisher",
        ],
        [
            "R13",
            "Pan-cancer pathomics",
            "Nat Cancer 2020",
            "Pan-cancer computational histopathology reveals mutations, tumor composition and prognosis",
            "10.1038/s43018-020-0085-8",
            "https://www.nature.com/articles/s43018-020-0085-8",
            "he_molecular,virtual_spatial,weak_supervision",
            "Supports the idea that histology features correlate with mutations, tumor composition, expression, and prognosis.",
            "Use to frame tumor composition and spatially resolved morphology.",
            "Pan-cancer association, not thyroid DM1 validation.",
            "high",
            "publisher",
        ],
        [
            "R14",
            "Weakly supervised WSI",
            "Nat Biomed Eng 2021",
            "Data-efficient and weakly supervised computational pathology on whole-slide images",
            "10.1038/s41551-020-00682-w",
            "https://www.nature.com/articles/s41551-020-00682-w",
            "foundation,weak_supervision,guardrail",
            "CLAM legitimizes attention/MIL-style WSI learning with slide-level labels and localization.",
            "Use in Methods/architecture justification.",
            "Method precedent, not biological validation.",
            "critical",
            "publisher",
        ],
        [
            "R15",
            "Clinical-grade WSI AI",
            "Nat Med 2019",
            "Clinical-grade computational pathology using weakly supervised deep learning on whole slide images",
            "10.1038/s41591-019-0508-1",
            "https://www.nature.com/articles/s41591-019-0508-1",
            "weak_supervision,guardrail",
            "Shows clinical-grade weakly supervised WSI AI at scale.",
            "Use sparingly for WSI AI validation ethos and scale expectations.",
            "Not molecular prediction.",
            "medium",
            "publisher",
        ],
        [
            "R16",
            "Pathology foundation model",
            "Nat Med 2024",
            "Towards a general-purpose foundation model for computational pathology",
            "10.1038/s41591-024-02857-3",
            "https://www.nature.com/articles/s41591-024-02857-3",
            "foundation,guardrail",
            "UNI is the core model precedent and directly matches Paper 2's feature backbone.",
            "Foreground in Methods, Fig 1, and data-use table.",
            "Do not claim we invented UNI or trained a foundation model.",
            "critical",
            "publisher",
        ],
        [
            "R17",
            "Vision-language pathology FM",
            "Nat Med 2024",
            "A visual-language foundation model for computational pathology",
            "10.1038/s41591-024-02856-4",
            "https://www.nature.com/articles/s41591-024-02856-4",
            "foundation",
            "CONCH supports pathology foundation-model maturation and transferable representations.",
            "Use as context in the foundation-model paragraph.",
            "Not used in current Paper 2 results.",
            "medium",
            "publisher",
        ],
        [
            "R18",
            "Whole-slide foundation model",
            "Nature 2024",
            "A whole-slide foundation model for digital pathology from real-world data",
            "10.1038/s41586-024-07441-w",
            "https://www.nature.com/articles/s41586-024-07441-w",
            "foundation,weak_supervision,guardrail",
            "Prov-GigaPath supports whole-slide context and real-world pretraining as field direction.",
            "Use as field context, especially if discussing external robustness.",
            "Not our backbone unless explicitly tested.",
            "high",
            "publisher",
        ],
        [
            "R19",
            "Clinical-grade pathology FM",
            "Nat Med 2024",
            "A foundation model for clinical-grade computational pathology and rare cancers detection",
            "10.1038/s41591-024-03141-0",
            "https://www.nature.com/articles/s41591-024-03141-0",
            "foundation,guardrail",
            "Virchow supports large-scale pathology FM claims and rare-cancer utility.",
            "Use in impact paragraph only as field context.",
            "Not our model; do not imply rare-cancer validation.",
            "medium",
            "publisher",
        ],
        [
            "R20",
            "Spatial multi-omics",
            "Cell 2020",
            "High-Spatial-Resolution Multi-Omics Sequencing via Deterministic Barcoding in Tissue",
            "10.1016/j.cell.2020.10.026",
            "https://www.sciencedirect.com/science/article/pii/S0092867420313908",
            "ffpe_spatial,virtual_spatial",
            "DBiT-seq is the wet-lab spatial multi-omics root behind the user's 2026 epi-Patho-DBiT precedent.",
            "Use as spatial-omics lineage reference.",
            "Experimental spatial omics, not computational H&E inference.",
            "medium",
            "publisher",
        ],
        [
            "R21",
            "FFPE spatial epigenomics",
            "Nat Commun 2026",
            "Spatially decoding genotype-associated epigenetic landscapes in human lymphoma FFPE tissues via epi-Patho-DBiT",
            "10.1038/s41467-026-71576-9",
            "https://www.nature.com/articles/s41467-026-71576-9",
            "ffpe_spatial,virtual_spatial,guardrail",
            "The user-provided Nature Communications precedent: pathology/FFPE tissue can unlock spatial molecular disease biology.",
            "Use as adjacent NC bridge: Paper 2 is the computational H&E front door.",
            "Lymphoma epigenomics, not H&E-to-RNA AI and not thyroid.",
            "critical",
            "publisher",
        ],
        [
            "R22",
            "Thyroid genomic foundation",
            "Cell 2014",
            "Integrated Genomic Characterization of Papillary Thyroid Carcinoma",
            "10.1016/j.cell.2014.09.050",
            "https://gdc.cancer.gov/about-data/publications/thca_2014",
            "thyroid_biology,he_molecular",
            "TCGA-THCA is the core public thyroid reference and provides the genomics/pathology context.",
            "Use in thyroid biology and dataset provenance.",
            "Does not establish H&E-to-DM1 prediction by itself.",
            "critical",
            "GDC/publisher record",
        ],
        [
            "R23",
            "Advanced thyroid biology",
            "J Clin Invest 2016",
            "Genomic and transcriptomic hallmarks of poorly differentiated and anaplastic thyroid cancers",
            "10.1172/JCI85271",
            "https://www.jci.org/articles/view/85271",
            "thyroid_biology,guardrail",
            "Landa et al. anchors advanced thyroid dedifferentiation and transcriptomic context.",
            "Use as thyroid/RAI biology reference and boundary for advanced disease.",
            "Advanced PDTC/ATC, not H&E spatial prediction.",
            "critical",
            "publisher",
        ],
        [
            "R24",
            "Breast histo-genomics review",
            "npj Breast Cancer 2023",
            "Biological insights and novel biomarker discovery through deep learning approaches in breast cancer histopathology",
            "10.1038/s41523-023-00518-1",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC10079681/",
            "he_molecular,he_spatial_rna",
            "Useful review language for histo-genomics and biomarker discovery from H&E.",
            "Use only as background; prioritize primary papers in main claims.",
            "Review article, not primary validation.",
            "low",
            "PMC",
        ],
    ]
    columns = [
        "ref_id",
        "category",
        "venue_year",
        "title",
        "doi",
        "url",
        "claim_tags",
        "what_it_strengthens",
        "best_use",
        "boundary",
        "priority",
        "source_kind",
    ]
    return pd.DataFrame(rows, columns=columns)


def make_derived_tables(refs: pd.DataFrame, m: dict) -> dict[str, pd.DataFrame]:
    layers = []
    for tag, label in CLAIM_LAYERS:
        hits = refs[refs["claim_tags"].map(lambda tags: tag in split_tags(tags))]
        layers.append(
            [
                tag,
                label,
                ", ".join(hits["ref_id"].tolist()),
                len(hits),
                {
                    "he_molecular": "Routine H&E contains recoverable molecular signal.",
                    "he_spatial_rna": "The spatial transcriptome can be approximated from H&E under matched-validation designs.",
                    "virtual_spatial": "Virtual spatial molecular maps are now a high-impact translational direction.",
                    "foundation": "Foundation pathology features are a legitimate backbone for transferable morphology.",
                    "weak_supervision": "Slide-level and weakly supervised WSI learning is accepted when leakage is controlled.",
                    "ffpe_spatial": "Pathology/FFPE sections are valid substrates for spatial molecular discovery.",
                    "thyroid_biology": "DM1/RAI framing rests on established thyroid differentiation biology.",
                    "guardrail": "Benchmarking, OOD validation, uncertainty, and caveat transparency are field expectations.",
                }[tag],
                {
                    "he_molecular": "Fig 1 premise; Table 1 reference block.",
                    "he_spatial_rna": "Fig 3 spatial validation; presub novelty paragraph scaffold.",
                    "virtual_spatial": "Title/subtitle framing and graphical abstract.",
                    "foundation": "Methods and model schematic.",
                    "weak_supervision": "Methods, leakage-control subsection, reviewer Q.",
                    "ffpe_spatial": "Introduction factual background; Nature Communications bridge.",
                    "thyroid_biology": "Dataset rationale and biological target definition.",
                    "guardrail": "Limitations/caveat boxes, not voice-protected prose.",
                }[tag],
            ]
        )
    claim_map = pd.DataFrame(
        layers,
        columns=["claim_tag", "claim_layer", "supporting_refs", "n_refs", "claim_utility", "paper2_placement"],
    )
    journal_map = pd.DataFrame(
        [
            [
                "Cell Reports Medicine",
                "current best route",
                "R01, R02, R03, R04, R10, R14, R16, R22, R23",
                "Translate H&E-to-spatial/molecular biomarker recovery into a thyroid cancer use case.",
                "Use current evidence; no need to pretend clinical deployment.",
            ],
            [
                "Nature Communications",
                "plausible with K2/Bundang validation",
                "R01, R02, R03, R04, R05, R10, R16, R21, R22, R23",
                "Multidisciplinary cross-modal biomarker recovery from routine pathology.",
                "Needs external H&E validation to avoid looking like a single-cohort application.",
            ],
            [
                "The Lancet Digital Health",
                "presub moonshot",
                "R01, R04, R05, R10, R14, R15, R16, R18",
                "Digital pathology triage into molecular/spatial biomarker states.",
                "Needs clinical workflow, external validation, and decision threshold clarity.",
            ],
            [
                "Nature Medicine",
                "not direct submit now",
                "R10, R11, R15, R16, R18, R19",
                "Only if K2/Bundang plus clinical utility/outcome layer emerges.",
                "Current Paper 2 lacks prospective clinical-impact evidence.",
            ],
            [
                "Nature Cancer",
                "biology presub only",
                "R12, R13, R21, R22, R23",
                "Cancer-state biology through image-recovered molecular/spatial signals.",
                "Needs stronger mechanism or therapy-response biology beyond prediction.",
            ],
        ],
        columns=["journal", "recommended_route", "refs_to_foreground", "impact_frame", "blocking_boundary"],
    )
    figure_map = pd.DataFrame(
        [
            ["Fig 1", "General architecture", "R01, R02, R03, R10, R16, R21", "H&E front door to latent molecular/spatial state"],
            ["Fig 2", "TCGA H&E image-DM1 classifier", "R11, R12, R13, R14, R16", f"UNI LOTO AUC {m['uni']['pooled_overall_auc']:.3f}"],
            [
                "Fig 3",
                "GSE250521 spatial validation",
                "R01, R02, R04, R05, R06, R07, R08",
                f"slide-centered rho {m['g250_stage']['dm1_existing_slide_centered_rho']:.3f}; PTC rho {m['g250_stage']['dm1_existing_ptc_rho']:.3f}",
            ],
            [
                "Fig 4",
                "GSE230424 external thyroid spatial support",
                "R04, R05, R10, R16, R23",
                f"H&E rho {m['g230']['top_he_sample_centered_rho']:.3f}; residual-target rho {m['g230_resid']['dm1_coord_qc_residual_target_he_rho']:.3f}",
            ],
            [
                "Fig 5",
                "Robustness/caveat moat",
                "R04, R05, R07, R14, R15",
                f"random-module residual p {m['g250_spec']['actual_residual_empirical_p']:.4f}; AITD rho {m['aitd']['ap_tls_he_sample_centered_rho']:.3f}",
            ],
            ["Fig 6", "Journal/validation unlock", "R01, R10, R16, R18, R21", "K2/Bundang H&E external validation is the ceiling-raiser"],
            ["Table 1", "Data and model provenance", "R14, R16, R22, R23", "Make every cohort and model dependency explicit"],
            ["Table 2", "Reference-to-claim matrix", "R01-R24", "Shows that the larger claim is anchored, not inflated"],
        ],
        columns=["paper2_item", "role", "refs_to_cite", "caption_or_table_use"],
    )
    overclaim = pd.DataFrame(
        [
            [
                "H&E can recover latent molecular/spatial biomarker states.",
                "R01, R02, R03, R10, R16",
                "H&E is a direct spatial-omics assay.",
                "Say prediction/recovery/triage; reserve measurement for wet-lab spatial assays.",
            ],
            [
                "Paper 2 validates this architecture in thyroid DM1/RAI.",
                "R22, R23 plus Paper 2 results",
                "The architecture is pan-disease validated.",
                "Call thyroid the validated exemplar; mention generality as a design frame.",
            ],
            [
                "Foundation features make the approach modern and reusable.",
                "R16, R17, R18, R19",
                "Paper 2 introduces a new foundation model.",
                "State UNI/CLAM were used as existing building blocks.",
            ],
            [
                "Spatial/FFPE pathology tissue is an NC-level discovery substrate.",
                "R20, R21",
                "epi-Patho-DBiT proves H&E-to-RNA prediction.",
                "Use as adjacent wet-lab precedent only.",
            ],
            [
                "The field expects benchmarking and caveat honesty.",
                "R04, R05, R07",
                "Modest correlations are failure.",
                "Present spatial pattern recovery, residual controls, and negative controls together.",
            ],
        ],
        columns=["safe_claim", "supporting_refs", "forbidden_overclaim", "wording_rule"],
    )
    return {"references": refs, "claim_map": claim_map, "journal_map": journal_map, "figure_map": figure_map, "overclaim": overclaim}


def draw_architecture_map(tables: dict[str, pd.DataFrame], m: dict) -> Path:
    img = np.full((1800, 2400, 3), BG, dtype=np.uint8)
    title(
        img,
        "Reference architecture",
        "Build Paper 2 On Four High-Impact Reference Pillars",
        "The impact move is not more citations. It is assigning each precedent to a claim layer and boundary.",
    )
    center = (1120, 900)
    box(img, (790, 700, 1450, 1040), PANEL2, GOLD)
    put(img, "Paper 2", (855, 780), 0.95, GOLD, 2)
    put(img, "Routine H&E -> DM1/RAI-linked spatial RNA state", (855, 850), 0.62, INK, 2, 520)
    put(
        img,
        f"UNI AUC {m['uni']['pooled_overall_auc']:.3f} | GSE250521 rho {m['g250_stage']['dm1_existing_slide_centered_rho']:.3f} | GSE230424 rho {m['g230']['top_he_sample_centered_rho']:.3f}",
        (855, 955),
        0.48,
        MUTED,
        1,
        520,
    )
    pillars = [
        ("H&E -> RNA/ST", "R01 R02 R03 R04 R05 R06 R07 R08 R09", (90, 370, 720, 735), TEAL),
        ("Virtual spatial omics", "R01 R08 R10 R20 R21", (90, 930, 720, 1295), BLUE),
        ("Foundation WSI AI", "R14 R15 R16 R17 R18 R19", (1520, 370, 2320, 735), GREEN),
        ("Thyroid biology", "R22 R23 + Paper 2 data", (1520, 930, 2320, 1295), GOLD),
        ("Benchmark guardrails", "R04 R05 R07 R14 R15 R21", (610, 1375, 1770, 1645), RED),
    ]
    for head, refs, xyxy, color in pillars:
        box(img, xyxy, PANEL, color)
        x0, y0, x1, _ = xyxy
        put(img, head, (x0 + 28, y0 + 60), 0.72, color, 2, x1 - x0 - 56)
        put(img, refs, (x0 + 28, y0 + 120), 0.54, INK, 2, x1 - x0 - 56)
        put(img, "Use as support, not as overclaim permission.", (x0 + 28, y0 + 190), 0.45, MUTED, 1, x1 - x0 - 56)
        cv2.arrowedLine(img, ((x0 + x1) // 2, (y0 + xyxy[3]) // 2), center, color, 3, cv2.LINE_AA, tipLength=0.04)
    box(img, (95, 1680, 2305, 1745), (42, 34, 22), GOLD)
    put(img, "Manuscript rule: every large claim must cite a precedent AND name its boundary.", (130, 1722), 0.62, GOLD, 2)
    path = OUT / "F01_reference_architecture_map.png"
    cv2.imwrite(str(path), img)
    return path


def draw_precedent_matrix(tables: dict[str, pd.DataFrame]) -> Path:
    refs = tables["references"]
    img = np.full((2350, 2600, 3), BG, dtype=np.uint8)
    title(
        img,
        "Precedent-to-claim matrix",
        "Which Reference Supports Which Paper 2 Claim Layer",
        "Gold cells are direct support; empty cells should not be used for that claim.",
    )
    x0, y0 = 70, 330
    id_w, title_w, col_w, row_h = 82, 690, 175, 72
    heads = ["ID", "Reference"] + [label for _, label in CLAIM_LAYERS]
    widths = [id_w, title_w] + [col_w] * len(CLAIM_LAYERS)
    x = x0
    for head, w in zip(heads, widths):
        box(img, (x, y0, x + w, y0 + 84), (42, 34, 22), GOLD)
        put(img, head, (x + 8, y0 + 35), 0.36, GOLD, 2, w - 16, 5)
        x += w
    colors = [TEAL, TEAL, BLUE, GREEN, GREEN, PURPLE, GOLD, RED]
    for i, row in refs.iterrows():
        yy = y0 + 84 + i * row_h
        tags = split_tags(row.claim_tags)
        base = PANEL if i % 2 == 0 else PANEL2
        box(img, (x0, yy, x0 + id_w, yy + row_h), base)
        put(img, row.ref_id, (x0 + 13, yy + 45), 0.44, GOLD if row.priority == "critical" else INK, 2)
        box(img, (x0 + id_w, yy, x0 + id_w + title_w, yy + row_h), base)
        put(img, f"{row.venue_year} | {row.title}", (x0 + id_w + 10, yy + 31), 0.32, INK, 1, title_w - 20, 5)
        x = x0 + id_w + title_w
        for j, (tag, _) in enumerate(CLAIM_LAYERS):
            fill = (38, 47, 60)
            outline = LINE
            box(img, (x, yy, x + col_w, yy + row_h), fill, outline)
            if tag in tags:
                cv2.circle(img, (x + col_w // 2, yy + row_h // 2), 15, colors[j], -1, cv2.LINE_AA)
            x += col_w
    path = OUT / "F02_precedent_to_claim_matrix.png"
    cv2.imwrite(str(path), img)
    return path


def draw_journal_ladder(tables: dict[str, pd.DataFrame]) -> Path:
    img = np.full((1700, 2300, 3), BG, dtype=np.uint8)
    title(
        img,
        "Journal leverage ladder",
        "Use Different Reference Sets For Different Editorial Doors",
        "This separates the realistic route from the presub moonshot route.",
    )
    colors = [GREEN, TEAL, GOLD, BLUE, RED]
    for i, row in tables["journal_map"].iterrows():
        yy = 330 + i * 245
        color = colors[i]
        box(img, (85, yy, 2215, yy + 185), PANEL if i % 2 == 0 else PANEL2, color)
        cv2.rectangle(img, (85, yy), (93, yy + 185), color, -1)
        put(img, row.journal, (125, yy + 52), 0.7, color, 2, 400)
        put(img, row.recommended_route.upper(), (125, yy + 105), 0.42, MUTED, 1, 400)
        put(img, row.refs_to_foreground, (575, yy + 55), 0.52, INK, 2, 490)
        put(img, row.impact_frame, (1110, yy + 48), 0.48, INK, 1, 510)
        put(img, "Boundary: " + row.blocking_boundary, (1640, yy + 48), 0.43, RED if i >= 2 else MUTED, 1, 520)
    box(img, (85, 1575, 2215, 1640), (42, 34, 22), GOLD)
    put(img, "Decision: current evidence is CRM-ready; K2/Bundang H&E is the Nat Commun/Lancet DH unlock.", (120, 1618), 0.58, GOLD, 2)
    path = OUT / "F03_journal_reference_ladder.png"
    cv2.imwrite(str(path), img)
    return path


def draw_boundary_map(tables: dict[str, pd.DataFrame]) -> Path:
    img = np.full((1700, 2300, 3), BG, dtype=np.uint8)
    title(
        img,
        "Citation boundary map",
        "How To Sound Bigger Without Overclaiming",
        "Each safe claim has a matching forbidden claim and a wording rule.",
    )
    y0 = 330
    for i, row in tables["overclaim"].iterrows():
        yy = y0 + i * 245
        box(img, (90, yy, 2210, yy + 190), PANEL if i % 2 == 0 else PANEL2)
        put(img, "SAFE", (125, yy + 45), 0.5, GREEN, 2)
        put(img, row.safe_claim, (125, yy + 92), 0.47, INK, 1, 520)
        put(img, "REFS", (700, yy + 45), 0.5, GOLD, 2)
        put(img, row.supporting_refs, (700, yy + 92), 0.52, GOLD, 2, 340)
        put(img, "DO NOT CLAIM", (1090, yy + 45), 0.5, RED, 2)
        put(img, row.forbidden_overclaim, (1090, yy + 92), 0.47, RED, 1, 420)
        put(img, "WORDING RULE", (1570, yy + 45), 0.5, BLUE, 2)
        put(img, row.wording_rule, (1570, yy + 92), 0.47, INK, 1, 560)
    path = OUT / "F04_citation_boundary_map.png"
    cv2.imwrite(str(path), img)
    return path


def draw_timeline(tables: dict[str, pd.DataFrame]) -> Path:
    img = np.full((1550, 2300, 3), BG, dtype=np.uint8)
    title(
        img,
        "Field timeline",
        "The Literature Has Moved Toward Paper 2's Bigger Frame",
        "2018-2026: mutation prediction -> H&E-to-RNA -> foundation models -> virtual spatial omics.",
    )
    events = [
        ("2018", "R11 Coudray", "H&E predicts mutations", TEAL),
        ("2019", "R15 Campanella", "clinical-grade weak supervision", GREEN),
        ("2020", "R02 ST-Net / R03 HE2RNA", "H&E-to-spatial RNA and RNA-seq", TEAL),
        ("2020", "R12/R13/R20", "pan-cancer molecular pathomics + DBiT", BLUE),
        ("2021", "R14 CLAM", "interpretable MIL for WSI", GREEN),
        ("2024", "R16-R19", "foundation pathology models mature", GOLD),
        ("2024", "R08/R09", "histology-ST integration broadens", TEAL),
        ("2025", "R04/R07", "benchmarking + STPath foundation direction", RED),
        ("2026", "R01/R05/R06/R10", "Cell/Nat Med/Nat Commun virtual spatial molecular pathology", GOLD),
        ("2026", "R21", "FFPE spatial epigenomics NC precedent", PURPLE),
    ]
    x_axis = 260
    cv2.line(img, (x_axis, 320), (x_axis, 1380), LINE, 4, cv2.LINE_AA)
    for i, (year, ref, note, color) in enumerate(events):
        yy = 350 + i * 112
        cv2.circle(img, (x_axis, yy), 18, color, -1, cv2.LINE_AA)
        put(img, year, (75, yy + 10), 0.58, color, 2)
        box(img, (315, yy - 45, 2125, yy + 48), PANEL if i % 2 == 0 else PANEL2, color)
        put(img, ref, (345, yy - 6), 0.58, color, 2, 430)
        put(img, note, (820, yy - 6), 0.52, INK, 1, 1160)
    box(img, (90, 1450, 2210, 1545), (42, 34, 22), GOLD)
    put(
        img,
        "Paper 2 should be framed as a disease-specific validation inside the 2026 virtual spatial molecular pathology wave.",
        (125, 1510),
        0.6,
        GOLD,
        2,
        1980,
    )
    path = OUT / "F05_reference_timeline_2018_2026.png"
    cv2.imwrite(str(path), img)
    return path


def draw_figure_citation_plan(tables: dict[str, pd.DataFrame]) -> Path:
    img = np.full((1800, 2400, 3), BG, dtype=np.uint8)
    title(
        img,
        "Figure citation plan",
        "Where Each Reference Set Should Enter Paper 2",
        "This makes the manuscript look designed, not citation-stuffed.",
    )
    rows = tables["figure_map"]
    y0 = 320
    widths = [170, 430, 520, 1060]
    heads = ["Item", "Role", "Refs", "Use"]
    x = 70
    for head, w in zip(heads, widths):
        box(img, (x, y0, x + w, y0 + 72), (42, 34, 22), GOLD)
        put(img, head, (x + 12, y0 + 45), 0.52, GOLD, 2)
        x += w
    for i, row in rows.iterrows():
        yy = y0 + 72 + i * 155
        vals = [row.paper2_item, row.role, row.refs_to_cite, row.caption_or_table_use]
        x = 70
        for j, (w, val) in enumerate(zip(widths, vals)):
            box(img, (x, yy, x + w, yy + 155), PANEL if i % 2 == 0 else PANEL2)
            color = GOLD if j == 0 else TEAL if j == 2 else INK
            put(img, val, (x + 12, yy + 45), 0.43, color, 1, w - 24)
            x += w
    path = OUT / "F06_figure_citation_plan.png"
    cv2.imwrite(str(path), img)
    return path


def write_outputs(tables: dict[str, pd.DataFrame], figs: list[Path], m: dict) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for name, df in tables.items():
        df.to_csv(OUT / f"paper2_reference_{name}.tsv", sep="\t", index=False)
    summary = {
        "created": "2026-05-10",
        "purpose": "Paper 2 reference leverage pack for high-impact framing",
        "reference_count": int(len(tables["references"])),
        "critical_reference_ids": tables["references"].query("priority == 'critical'")["ref_id"].tolist(),
        "headline_frame": "Routine H&E is a general front door to latent molecular/spatial biomarker states; thyroid DM1/RAI is the validated exemplar.",
        "current_best_route": "Cell Reports Medicine",
        "ceiling_unlock": "K2/Bundang external H&E validation",
        "headline_metrics": {
            "uni_loto_auc": m["uni"]["pooled_overall_auc"],
            "gse250521_slide_centered_rho": m["g250_stage"]["dm1_existing_slide_centered_rho"],
            "gse230424_he_rho": m["g230"]["top_he_sample_centered_rho"],
            "gse230424_residual_target_rho": m["g230_resid"]["dm1_coord_qc_residual_target_he_rho"],
            "gse250521_random_module_residual_p": m["g250_spec"]["actual_residual_empirical_p"],
            "gse248205_aitd_rho": m["aitd"]["ap_tls_he_sample_centered_rho"],
        },
        "figures": [fig.name for fig in figs],
        "voice_boundary": "Figure/table/scaffold only; no protected manuscript voice sections were generated.",
    }
    (OUT / "PAPER2_REFERENCE_LEVERAGE_PACK.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    md = [
        "# Paper 2 reference leverage pack",
        "",
        "Purpose: convert high-impact precedents into a claim-by-claim scaffold for Paper 2.",
        "",
        f"References curated: {len(tables['references'])}",
        "",
        "Headline framing:",
        "- Routine H&E is a general front door to latent molecular/spatial biomarker states.",
        "- Thyroid DM1/RAI is the validated exemplar, not a pan-disease proof.",
        "- K2/Bundang external H&E validation is the ceiling unlock.",
        "",
        "Figures:",
    ]
    for fig in figs:
        md.append(f"- `{fig.name}`")
    md += ["", "Voice boundary: scaffold/table/figure package only."]
    (OUT / "SUMMARY.md").write_text("\n".join(md) + "\n", encoding="utf-8")


def publish(figs: list[Path]) -> None:
    ASSET_LOCAL.mkdir(parents=True, exist_ok=True)
    ASSET_LIVE.mkdir(parents=True, exist_ok=True)
    for fig in figs:
        shutil.copy2(fig, ASSET_LOCAL / fig.name)
        shutil.copy2(fig, ASSET_LIVE / fig.name)


def html_table(df: pd.DataFrame) -> str:
    out = [
        "<table><thead><tr>"
        + "".join(f"<th>{html.escape(str(c))}</th>" for c in df.columns)
        + "</tr></thead><tbody>"
    ]
    for _, row in df.iterrows():
        cells = []
        for c in df.columns:
            value = html.escape(str(row[c]))
            if c == "url":
                value = f'<a href="{value}" target="_blank" rel="noopener">{value}</a>'
            cells.append(f"<td>{value}</td>")
        out.append("<tr>" + "".join(cells) + "</tr>")
    out.append("</tbody></table>")
    return "\n".join(out)


def write_page(tables: dict[str, pd.DataFrame], figs: list[Path], m: dict) -> None:
    fig_buttons = "\n".join(
        f'<button class="fig-option" data-src="assets/{ASSET}/{fig.name}" data-title="{html.escape(fig.stem)}">'
        f'<img src="assets/{ASSET}/{fig.name}" alt=""><span>{html.escape(fig.stem.replace("_", " "))}</span></button>'
        for fig in figs
    )
    html_doc = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Paper 2 Reference Leverage Pack</title>
<style>
:root {{ --bg:#0d1117; --panel:#151b23; --panel2:#1c2634; --ink:#e6edf3; --muted:#98a6ba; --line:#303846; --gold:#d6b25e; --blue:#57a6ff; --green:#9dd335; --red:#ff7b72; }}
* {{ box-sizing:border-box; }}
body {{ margin:0; background:var(--bg); color:var(--ink); font-family:"JetBrains Mono", ui-monospace, Menlo, monospace; line-height:1.55; }}
a {{ color:var(--blue); text-decoration:none; overflow-wrap:anywhere; }}
a:hover {{ text-decoration:underline; }}
.hero {{ min-height:64vh; display:grid; align-items:end; padding:58px 5vw 42px; border-bottom:1px solid var(--line); background:linear-gradient(180deg, rgba(13,17,23,.1), rgba(13,17,23,.97)), url("assets/{ASSET}/{figs[0].name}") center / cover no-repeat; }}
.kicker {{ color:var(--gold); text-transform:uppercase; font-size:13px; font-weight:900; }}
h1 {{ margin:10px 0 14px; max-width:1180px; font-family:"Cormorant Garamond", Georgia, serif; font-size:clamp(42px,7vw,88px); line-height:.96; letter-spacing:0; }}
.lead {{ max-width:1120px; color:#d8dee8; font-size:18px; }}
.stats {{ display:grid; grid-template-columns:repeat(6,minmax(0,1fr)); gap:10px; max-width:1320px; margin-top:24px; }}
.stat {{ border:1px solid rgba(214,178,94,.45); background:rgba(13,17,23,.84); padding:14px; min-height:92px; }}
.stat .v {{ color:var(--gold); font-size:23px; font-weight:900; }}
.stat .l {{ color:var(--muted); font-size:12px; margin-top:5px; }}
.layout {{ display:grid; grid-template-columns:280px minmax(0,1fr); gap:32px; max-width:1580px; margin:0 auto; padding:34px 24px 84px; }}
nav {{ position:sticky; top:0; align-self:start; max-height:100vh; overflow:auto; padding:18px 0; }}
nav a {{ display:block; padding:9px 12px; border-left:2px solid transparent; color:var(--muted); font-size:13px; }}
nav a:hover {{ color:var(--ink); border-left-color:var(--gold); text-decoration:none; }}
section {{ border-top:1px solid var(--line); padding:34px 0; }}
section:first-child {{ border-top:0; padding-top:0; }}
h2 {{ margin:0 0 18px; font-family:"Cormorant Garamond", Georgia, serif; font-size:34px; letter-spacing:0; }}
h3 {{ margin:0 0 10px; font-size:17px; }}
.num {{ color:var(--gold); margin-right:10px; }}
.grid {{ display:grid; grid-template-columns:repeat(3,minmax(0,1fr)); gap:16px; }}
.card {{ border:1px solid var(--line); background:var(--panel); border-radius:8px; padding:18px; }}
.muted {{ color:var(--muted); }}
table {{ width:100%; border-collapse:collapse; margin:16px 0; font-size:12px; }}
th,td {{ border-bottom:1px solid var(--line); padding:10px 8px; text-align:left; vertical-align:top; }}
th {{ color:var(--gold); background:rgba(214,178,94,.06); position:sticky; top:0; }}
td:nth-child(1) {{ color:var(--gold); font-weight:800; }}
.viewer {{ display:grid; grid-template-columns:370px minmax(0,1fr); gap:16px; }}
.fig-list {{ display:grid; gap:10px; align-content:start; max-height:780px; overflow:auto; padding-right:6px; }}
.fig-option {{ border:1px solid var(--line); background:var(--panel); color:var(--ink); border-radius:8px; padding:10px; display:grid; grid-template-columns:96px minmax(0,1fr); gap:10px; text-align:left; cursor:pointer; font:inherit; }}
.fig-option.active,.fig-option:hover {{ border-color:var(--gold); background:#202938; }}
.fig-option img {{ width:96px; height:66px; object-fit:cover; border-radius:4px; display:block; }}
.stage {{ border:1px solid var(--line); background:#090b10; border-radius:8px; overflow:hidden; }}
.toolbar {{ display:flex; gap:8px; align-items:center; padding:12px; background:var(--panel); border-bottom:1px solid var(--line); }}
.toolbar strong {{ flex:1 1 auto; }}
button.control {{ border:1px solid var(--line); background:var(--panel2); color:var(--ink); border-radius:6px; padding:8px 10px; font:inherit; font-size:12px; cursor:pointer; }}
button.control:hover {{ border-color:var(--gold); color:var(--gold); }}
.canvas {{ height:min(78vh,860px); overflow:auto; padding:16px; background:#090b10; }}
.canvas img {{ display:block; width:100%; min-width:640px; max-width:none; height:auto; margin:0 auto; }}
@media(max-width:980px) {{ .layout,.viewer {{ grid-template-columns:1fr; }} nav {{ position:static; max-height:none; display:grid; grid-template-columns:repeat(2,minmax(0,1fr)); }} .stats,.grid {{ grid-template-columns:1fr; }} }}
</style>
</head>
<body>
<header class="hero">
  <div>
    <div class="kicker">Paper 2 · reference leverage pack · 2026-05-10</div>
    <h1>More References, But Organized As Impact Architecture</h1>
    <p class="lead">A claim-by-claim literature map for raising Paper 2 from a thyroid image model into an H&E-to-molecular/spatial biomarker recovery manuscript.</p>
    <div class="stats">
      <div class="stat"><div class="v">{len(tables['references'])}</div><div class="l">curated precedents</div></div>
      <div class="stat"><div class="v">{m['uni']['pooled_overall_auc']:.3f}</div><div class="l">UNI LOTO AUC</div></div>
      <div class="stat"><div class="v">{m['g250_stage']['dm1_existing_slide_centered_rho']:.3f}</div><div class="l">GSE250521 rho</div></div>
      <div class="stat"><div class="v">{m['g230']['top_he_sample_centered_rho']:.3f}</div><div class="l">GSE230424 rho</div></div>
      <div class="stat"><div class="v">{m['hotspot']['gse230424_top10_lift']:.1f}x</div><div class="l">external hotspot lift</div></div>
      <div class="stat"><div class="v">K2</div><div class="l">ceiling unlock</div></div>
    </div>
  </div>
</header>
<div class="layout">
<nav>
  <a href="#verdict">01 Verdict</a>
  <a href="#figures">02 Figures</a>
  <a href="#references">03 References</a>
  <a href="#claim-map">04 Claim Map</a>
  <a href="#journal-map">05 Journal Map</a>
  <a href="#figure-map">06 Figure Plan</a>
  <a href="#boundaries">07 Boundaries</a>
  <a href="#sources">08 Local Sources</a>
</nav>
<main>
<section id="verdict">
  <h2><span class="num">01</span>Reference Leverage Verdict</h2>
  <div class="grid">
    <div class="card"><h3>Impact frame</h3><p class="muted">Routine H&E is the general front door; thyroid DM1/RAI is the validated exemplar.</p></div>
    <div class="card"><h3>Strongest new analogs</h3><p class="muted">Cell 2026 H&E-to-spatial transcriptomics and Nat Med 2026 HEX make the virtual spatial molecular pathology frame timely.</p></div>
    <div class="card"><h3>Boundary</h3><p class="muted">References raise plausibility; they do not replace K2/Bundang external H&E validation.</p></div>
  </div>
</section>
<section id="figures">
  <h2><span class="num">02</span>Reference Figure Viewer</h2>
  <div class="viewer">
    <div class="fig-list" id="figList">{fig_buttons}</div>
    <div class="stage">
      <div class="toolbar"><strong id="figTitle">Figure</strong><button class="control" id="zoomOut">Zoom out</button><button class="control" id="zoomReset">Reset</button><button class="control" id="zoomIn">Zoom in</button></div>
      <div class="canvas" id="canvas"><img id="mainFig" src="assets/{ASSET}/{figs[0].name}" alt="Paper 2 reference figure"></div>
    </div>
  </div>
</section>
<section id="references"><h2><span class="num">03</span>Curated Precedents</h2>{html_table(tables['references'])}</section>
<section id="claim-map"><h2><span class="num">04</span>Claim-Layer Map</h2>{html_table(tables['claim_map'])}</section>
<section id="journal-map"><h2><span class="num">05</span>Journal Reference Map</h2>{html_table(tables['journal_map'])}</section>
<section id="figure-map"><h2><span class="num">06</span>Paper 2 Figure/Table Citation Plan</h2>{html_table(tables['figure_map'])}</section>
<section id="boundaries"><h2><span class="num">07</span>Overclaim Guardrails</h2>{html_table(tables['overclaim'])}</section>
<section id="sources">
  <h2><span class="num">08</span>Local Sources</h2>
  <table><tbody>
    <tr><td>Output directory</td><td><code>{OUT.relative_to(ROOT)}</code></td></tr>
    <tr><td>Builder script</td><td><code>project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/scripts/build_paper2_reference_leverage_pack.py</code></td></tr>
    <tr><td>User PDF</td><td><code>AI-predicted spatial transcriptomics unlocks breast cancer biomarkers from pathology.pdf</code></td></tr>
    <tr><td>Moonshot pack</td><td><a href="paper2_moonshot_editorial_pack.html">paper2_moonshot_editorial_pack.html</a></td></tr>
    <tr><td>Submission overview</td><td><a href="paper2_submission_overview.html">paper2_submission_overview.html</a></td></tr>
  </tbody></table>
</section>
</main>
</div>
<script>
const buttons = [...document.querySelectorAll('.fig-option')];
const img = document.getElementById('mainFig');
const title = document.getElementById('figTitle');
const canvas = document.getElementById('canvas');
let zoom = 1;
function setZoom(next) {{ zoom = Math.min(3.5, Math.max(.65, next)); img.style.width = `${{zoom * 100}}%`; }}
buttons.forEach((button, index) => {{
  button.addEventListener('click', () => {{
    buttons.forEach(b => b.classList.remove('active'));
    button.classList.add('active');
    img.src = button.dataset.src;
    title.textContent = button.dataset.title.replaceAll('_', ' ');
    canvas.scrollTo({{ top: 0, left: 0 }});
    setZoom(1);
  }});
  if (index === 0) button.classList.add('active');
}});
document.getElementById('zoomOut').addEventListener('click', () => setZoom(zoom - .2));
document.getElementById('zoomReset').addEventListener('click', () => {{ setZoom(1); canvas.scrollTo({{ top: 0, left: 0 }}); }});
document.getElementById('zoomIn').addEventListener('click', () => setZoom(zoom + .2));
title.textContent = buttons[0].dataset.title.replaceAll('_', ' ');
</script>
</body>
</html>
"""
    PAGE_LOCAL.write_text(html_doc, encoding="utf-8")
    shutil.copy2(PAGE_LOCAL, PAGE_LIVE)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    m = metrics()
    refs = make_reference_table()
    tables = make_derived_tables(refs, m)
    figs = [
        draw_architecture_map(tables, m),
        draw_precedent_matrix(tables),
        draw_journal_ladder(tables),
        draw_boundary_map(tables),
        draw_timeline(tables),
        draw_figure_citation_plan(tables),
    ]
    write_outputs(tables, figs, m)
    publish(figs)
    write_page(tables, figs, m)
    print(json.dumps({"out": str(OUT), "page": str(PAGE_LOCAL), "figures": [p.name for p in figs]}, indent=2))


if __name__ == "__main__":
    main()
