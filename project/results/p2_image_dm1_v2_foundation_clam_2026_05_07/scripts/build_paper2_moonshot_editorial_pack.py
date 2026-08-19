#!/usr/bin/env python3
"""Build a Paper 2 moonshot editorial package."""
from __future__ import annotations

import json
import shutil
from pathlib import Path

import cv2
import numpy as np
import pandas as pd


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
P2 = ROOT / "project/results/p2_image_dm1_v2_foundation_clam_2026_05_07"
SUPP = P2 / "analysis_supp"
OUT = SUPP / "paper2_moonshot_editorial_pack_2026_05_10"
HUB = ROOT / "project/papers_hub_2026_05_04"
LIVE = Path("/var/www/papers/papers_hub_2026_05_04")
ASSET = "paper2_moonshot_editorial_pack"
ASSET_LOCAL = HUB / "assets" / ASSET
ASSET_LIVE = LIVE / "assets" / ASSET
PAGE_LOCAL = HUB / "paper2_moonshot_editorial_pack.html"
PAGE_LIVE = LIVE / "paper2_moonshot_editorial_pack.html"

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
WHITE = (245, 246, 248)


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
    put(img, heading, (70, 142), 1.12, INK, 2, 1600)
    put(img, sub, (70, 205), 0.56, MUTED, 1, 1520)


def fit_image(path: Path, size: tuple[int, int]) -> np.ndarray:
    tw, th = size
    src = cv2.imread(str(path))
    if src is None:
        canvas = np.full((th, tw, 3), PANEL2, dtype=np.uint8)
        put(canvas, f"missing {path.name}", (20, th // 2), 0.45, RED, 1, tw - 40)
        return canvas
    h, w = src.shape[:2]
    scale = min(tw / max(w, 1), th / max(h, 1))
    nw, nh = max(1, int(w * scale)), max(1, int(h * scale))
    resized = cv2.resize(src, (nw, nh), interpolation=cv2.INTER_AREA)
    canvas = np.full((th, tw, 3), WHITE, dtype=np.uint8)
    x0 = (tw - nw) // 2
    y0 = (th - nh) // 2
    canvas[y0 : y0 + nh, x0 : x0 + nw] = resized
    return canvas


def make_tables(m: dict) -> dict[str, pd.DataFrame]:
    claims = pd.DataFrame(
        [
            [
                "Current-safe",
                "Cell Reports Medicine / Communications Medicine",
                "Routine H&E image analysis recovers a DM1/RAI-linked spatial RNA state in thyroid cancer, as proof-of-principle for an H&E-to-molecular-state architecture.",
                "TCGA UNI LOTO AUC 0.852; GSE250521 rho 0.440; GSE230424 rho 0.644.",
                "Deployment-ready clinical tool; pan-disease validated platform; causal MAPK/RAI mechanism.",
            ],
            [
                "NC-tier",
                "Nature Communications",
                "Routine histopathology contains recoverable spatial transcriptomic biomarker structure; thyroid DM1/RAI is the validated exemplar.",
                "Current evidence plus K2/Bundang H&E external validation or concrete revision-stage validation.",
                "General pan-cancer validation without additional cohorts; prospective clinical decision support.",
            ],
            [
                "Moonshot",
                "The Lancet Digital Health",
                "H&E can serve as a general digital front-door assay for latent molecular/spatial biomarker states, validated here in thyroid cancer.",
                "K2/Bundang external H&E; clinical triage scenario; implementation relevance.",
                "Practice-changing claim without clinical workflow validation.",
            ],
            [
                "Cancer-biology reach",
                "Nature Cancer",
                "H&E recovers a cancer-state readout linked to DM1/RAI biology and immune overlap in thyroid cancer.",
                "Paper 1 public/preprint plus mechanism/therapy relevance.",
                "New causal cancer mechanism from image model alone.",
            ],
        ],
        columns=["claim_level", "target", "allowable_claim", "minimum_evidence", "forbidden_overclaim"],
    )
    editorial = pd.DataFrame(
        [
            ["Novelty", 9, "General H&E front-door assay framing: cross-modal H&E -> spatial RNA biomarker recovery."],
            ["Validation", 6, "TCGA + GSE250521 + GSE230424 strong; K2/Bundang external H&E would raise this."],
            ["Clinical pull", 5, "Triage logic exists, but prospective/clinical utility is not yet demonstrated."],
            ["Visual clarity", 9, "Separated 81-figure browser plus submission overview/impact dossier."],
            ["Caveat honesty", 9, "QC/smoothness, AITD no-go, and residual limits are co-located."],
            ["Moonshot readiness", 6, "Enough for presub; not enough for direct LDH/Nature Medicine full submission."],
        ],
        columns=["dimension", "score_10", "editor_note"],
    )
    attachment = pd.DataFrame(
        [
            ["A", "Submission overview", "paper2_submission_overview.html", "Figure/table/data-use provenance"],
            ["B", "Impact dossier", "paper2_impact_upgrade.html", "Evidence moat + validation unlock"],
            ["C", "High-impact journal strategy", "paper2_high_impact_journal_strategy.html", "Venue ladder and decision tree"],
            ["D", "Separated figure browser", "paper2_cv2_visual_summary.html", "81 source figures with zoom/fullscreen"],
            ["E", "K2/Bundang validation plan", "PAPER2_KOREAN_K2_VALIDATION_PLAN_2026_05_08.md", "Highest-impact unlock"],
        ],
        columns=["slot", "attachment", "path_or_url", "purpose"],
    )
    validation = pd.DataFrame(
        [
            ["1", "Secure H&E availability", "K2/Bundang", "highest", "Determines whether NC/LDH route is realistic"],
            ["2", "Run locked inference", "same model / no retuning", "highest", "Avoids leakage and reviewer criticism"],
            ["3", "Report external AUC/rho", "DM1 labels or proxy score", "highest", "Main missing validation"],
            ["4", "Add calibration/decision threshold", "external cohort", "medium", "Needed for clinical/digital-health framing"],
            ["5", "Update Fig 6 and Table 1/2", "submission pack", "high", "Turns roadmap into evidence"],
        ],
        columns=["step", "task", "data", "priority", "why_it_raises_impact"],
    )
    presub = pd.DataFrame(
        [
            ["Thesis", "Routine H&E can act as a general digital front door to latent molecular/spatial biomarker states; this study validates the architecture in thyroid DM1/RAI biology.", "Fact slot only; author writes final sentence."],
            ["Evidence", f"UNI LOTO AUC {m['uni']['pooled_overall_auc']:.3f}; GSE250521 rho {m['g250_stage']['dm1_existing_slide_centered_rho']:.3f}; GSE230424 rho {m['g230']['top_he_sample_centered_rho']:.3f}; hotspot lift {m['hotspot']['gse230424_top10_lift']:.1f}x.", "Do not omit QC/smoothness caveat."],
            ["Boundary", f"GSE230424 residual-target rho {m['g230_resid']['dm1_coord_qc_residual_target_he_rho']:.3f}; GSE248205 AITD control rho {m['aitd']['ap_tls_he_sample_centered_rho']:.3f}.", "This honesty increases editor trust."],
            ["Ask", "Does this fit as digital pathology biomarker recovery / cross-modal translational medicine?", "Presub inquiry, not final cover letter."],
        ],
        columns=["component", "fact_scaffold", "guardrail"],
    )
    generalization = pd.DataFrame(
        [
            [
                "Input assay",
                "Routine H&E",
                "Standard histopathology stain across cancer and inflammatory pathology workflows",
                "Validated here using thyroid H&E cohorts",
                "A pan-disease validation claim",
            ],
            [
                "Model layer",
                "Foundation pathology features + MIL",
                "Domain-adaptable morphology representation rather than a thyroid-only hand-crafted feature set",
                "Applied here to DM1/RAI state recovery",
                "A new foundation model invention",
            ],
            [
                "Target layer",
                "Latent molecular/spatial biomarker state",
                "Target can be changed: RAI/DM1, immune niche, therapy vulnerability, fibrosis, hypoxia, etc.",
                "Current target is DM1/RAI-linked thyroid state",
                "All targets validated today",
            ],
            [
                "Validation layer",
                "Spatial transcriptomics + external H&E",
                "Same validation pattern can be reused in other diseases with matched H&E/spatial RNA",
                "GSE250521 and GSE230424 support the pattern",
                "Prospective clinical utility",
            ],
            [
                "Clinical layer",
                "Reflex-test triage",
                "H&E-negative patients avoid reflex testing; H&E-positive patients move to targeted molecular confirmation",
                "Current status is hypothesis/strategy",
                "Deployment-ready threshold",
            ],
        ],
        columns=["layer", "paper2_unit", "generalizable_design", "current_evidence", "not_claimed"],
    )
    precedent = pd.DataFrame(
        [
            [
                "Nature Communications 2026",
                "Li et al., epi-Patho-DBiT",
                "Spatially decoding genotype-associated epigenetic landscapes in human lymphoma FFPE tissues via epi-Patho-DBiT",
                "Shows that archival/pathology FFPE tissue can be a substrate for spatial molecular discovery, not only histologic diagnosis.",
                "Use as adjacent NC precedent for pathology tissue -> spatial omics biology; not as an AI benchmark.",
                "Paper 2 is computational H&E-to-spatial-RNA prediction; Li et al. is wet-lab spatial epigenomics from FFPE.",
                "https://www.nature.com/articles/s41467-026-71576-9",
            ],
            [
                "What it strengthens",
                "General H&E/FFPE framing",
                "Pathology sections can encode spatial molecular states and disease-genotype biology.",
                "Raises the credibility of a generalized H&E front-door assay framing.",
                "Mention in Introduction/Discussion factual background slots, not protected voice prose here.",
                "Do not imply they predicted RNA from H&E or studied thyroid cancer.",
                "https://doi.org/10.1038/s41467-026-71576-9",
            ],
            [
                "Paper 2 contrast",
                "Computational vs experimental spatial omics",
                "Our route uses routine H&E images plus foundation morphology features to predict a spatial RNA biomarker state.",
                "Makes the manuscript complementary: cheaper/front-door triage vs direct spatial profiling.",
                "Use in Figure 1/claim ladder as a precedent box.",
                "Do not claim equivalence to direct ATAC/CUT&Tag measurement.",
                "https://www.nature.com/articles/s41467-026-71576-9",
            ],
        ],
        columns=["precedent", "short_name", "article_title", "why_it_matters", "how_to_use", "boundary", "source_url"],
    )
    return {
        "claims": claims,
        "editorial": editorial,
        "attachment": attachment,
        "validation": validation,
        "presub": presub,
        "generalization": generalization,
        "precedent": precedent,
    }


def draw_he_generalization_frame(tables: dict[str, pd.DataFrame]) -> Path:
    img = np.full((1700, 2200, 3), BG, dtype=np.uint8)
    title(
        img,
        "General H&E framing",
        "H&E As A General Front-Door Assay",
        "Make the paper read as a reusable cross-modal architecture, while keeping the evidence anchored to thyroid DM1/RAI.",
    )
    rows = tables["generalization"]
    x0, y0 = 70, 320
    cols = [260, 310, 560, 420, 410]
    heads = ["Layer", "Paper 2 unit", "Generalizable design", "Current evidence", "Not claimed"]
    x = x0
    for w, head in zip(cols, heads):
        box(img, (x, y0, x + w, y0 + 62), (42, 34, 22), GOLD)
        put(img, head, (x + 12, y0 + 40), 0.46, GOLD, 2)
        x += w
    for i, row in rows.iterrows():
        yy = y0 + 62 + i * 195
        vals = [row.layer, row.paper2_unit, row.generalizable_design, row.current_evidence, row.not_claimed]
        x = x0
        for j, (w, val) in enumerate(zip(cols, vals)):
            color = RED if j == 4 else INK if j in {0, 1, 3} else MUTED
            box(img, (x, yy, x + w, yy + 195), PANEL if i % 2 == 0 else PANEL2)
            put(img, val, (x + 12, yy + 38), 0.43, color, 1, w - 24)
            x += w
    box(img, (70, 1425, 2030, 1605), PANEL2)
    put(img, "Editorial sentence skeleton", (105, 1480), 0.78, GOLD, 2)
    put(
        img,
        "Routine H&E is the universal input; thyroid DM1/RAI is the validated exemplar; spatial transcriptomics is the molecular readout; K2/Bundang is the external unlock.",
        (105, 1540),
        0.56,
        INK,
        1,
        1760,
    )
    path = OUT / "F00_general_he_front_door_frame.png"
    cv2.imwrite(str(path), img)
    return path


def draw_graphical_abstract(m: dict) -> Path:
    img = np.full((1700, 2200, 3), BG, dtype=np.uint8)
    title(
        img,
        "Moonshot graphical abstract",
        "Routine H&E -> General Molecular/Spatial State Recovery",
        "The visual pitch: H&E is the universal assay; thyroid DM1/RAI is the validated proof-of-principle; caveats define the current boundary.",
    )
    nodes = [
        ("Routine H&E", "universal pathology input", (210, 505), GOLD),
        ("Foundation morphology", "UNI / ViT-L + MIL", (620, 505), GREEN),
        ("Target state", "molecular / spatial RNA", (1030, 505), TEAL),
        ("Validated exemplar", f"thyroid DM1/RAI rho {m['g250_stage']['dm1_existing_slide_centered_rho']:.3f}", (1440, 505), BLUE),
        ("External unlock", "K2/Bundang H&E", (1850, 505), RED),
    ]
    for head, note, (cx, cy), color in nodes:
        box(img, (cx - 170, cy - 105, cx + 170, cy + 105), PANEL2, color)
        put(img, head, (cx - 135, cy - 35), 0.58, color, 2, 270)
        put(img, note, (cx - 135, cy + 32), 0.43, INK, 1, 270)
    for i in range(len(nodes) - 1):
        a = (nodes[i][2][0] + 170, nodes[i][2][1])
        b = (nodes[i + 1][2][0] - 170, nodes[i + 1][2][1])
        cv2.arrowedLine(img, a, b, GOLD, 3, cv2.LINE_AA, tipLength=0.08)

    cards = [
        ("Image anchor", f"UNI LOTO AUC {m['uni']['pooled_overall_auc']:.3f}", GREEN),
        ("Spatial recovery", f"PTC rho {m['g250_stage']['dm1_existing_ptc_rho']:.3f} / LPTC {m['g250_stage']['dm1_existing_lptc_rho']:.3f}", TEAL),
        ("External support", f"Residual-target rho {m['g230_resid']['dm1_coord_qc_residual_target_he_rho']:.3f}", BLUE),
        ("Hotspot localization", f"{m['hotspot']['gse250521_top10_lift']:.1f}x / {m['hotspot']['gse230424_top10_lift']:.1f}x", GOLD),
        ("No-overclaim lock", f"AITD control rho {m['aitd']['ap_tls_he_sample_centered_rho']:.3f}", RED),
    ]
    y = 900
    for i, (head, val, color) in enumerate(cards):
        x = 70 + i * 400
        box(img, (x, y, x + 360, y + 160), PANEL)
        cv2.rectangle(img, (x, y), (x + 360, y + 8), color, -1)
        put(img, head.upper(), (x + 20, y + 42), 0.43, MUTED, 1)
        put(img, val, (x + 20, y + 98), 0.72, color, 2, 310)
        put(img, "editor-visible metric", (x + 20, y + 132), 0.38, INK, 1)

    box(img, (70, 1260, 2030, 1540), PANEL2)
    put(img, "Moonshot thesis", (105, 1325), 0.9, GOLD, 2)
    put(
        img,
        "H&E is not the disease-specific novelty. The novelty is using a universal pathology input to recover a latent molecular/spatial biomarker state. This paper validates the architecture in thyroid DM1/RAI biology and explicitly marks what remains unvalidated.",
        (105, 1400),
        0.62,
        INK,
        1,
        1780,
    )
    path = OUT / "F01_moonshot_graphical_abstract.png"
    cv2.imwrite(str(path), img)
    return path


def draw_claim_ladder(tables: dict[str, pd.DataFrame]) -> Path:
    img = np.full((1600, 2100, 3), BG, dtype=np.uint8)
    title(img, "Claim ladder", "Use The Strongest Claim Each Journal Can Safely Bear", "This is how to go bigger without overclaiming.")
    claims = tables["claims"]
    y0 = 320
    for i, row in claims.iterrows():
        yy = y0 + i * 285
        color = [GREEN, TEAL, GOLD, BLUE][i]
        box(img, (70, yy, 2020, yy + 220), PANEL if i % 2 == 0 else PANEL2)
        cv2.rectangle(img, (70, yy), (78, yy + 220), color, -1)
        put(img, row.claim_level, (110, yy + 48), 0.68, color, 2, 280)
        put(img, row.target, (410, yy + 48), 0.58, INK, 2, 420)
        put(img, row.allowable_claim, (850, yy + 45), 0.52, INK, 1, 560)
        put(img, "Evidence: " + row.minimum_evidence, (850, yy + 118), 0.44, MUTED, 1, 560)
        put(img, "Do not claim: " + row.forbidden_overclaim, (1450, yy + 55), 0.46, RED, 1, 500)
    path = OUT / "F02_moonshot_claim_ladder.png"
    cv2.imwrite(str(path), img)
    return path


def draw_editor_scorecard(tables: dict[str, pd.DataFrame]) -> Path:
    img = np.full((1650, 2100, 3), BG, dtype=np.uint8)
    title(img, "Editor scorecard", "What Looks Strong, What Still Blocks The Moonshot", "A compact editorial readiness map for advisor triage.")
    y0 = 330
    for i, row in tables["editorial"].iterrows():
        yy = y0 + i * 150
        score = int(row.score_10)
        color = GREEN if score >= 8 else GOLD if score >= 6 else RED
        box(img, (70, yy, 2020, yy + 108), PANEL if i % 2 == 0 else PANEL2)
        put(img, row.dimension, (105, yy + 50), 0.68, INK, 2, 380)
        cv2.rectangle(img, (530, yy + 35), (1030, yy + 65), (43, 50, 62), -1)
        cv2.rectangle(img, (530, yy + 35), (530 + int(500 * score / 10), yy + 65), color, -1)
        cv2.rectangle(img, (530, yy + 35), (1030, yy + 65), LINE, 1, cv2.LINE_AA)
        put(img, f"{score}/10", (1065, yy + 62), 0.54, color, 2)
        put(img, row.editor_note, (1190, yy + 43), 0.48, MUTED, 1, 740)
    box(img, (70, 1280, 2020, 1540), PANEL2)
    put(img, "Upgrade rule", (105, 1340), 0.9, GOLD, 2)
    put(img, "K2/Bundang H&E validation is the one move that improves validation, clinical pull, Nature Communications readiness, and Lancet Digital Health presub credibility at the same time.", (105, 1415), 0.62, INK, 1, 1780)
    path = OUT / "F03_editor_scorecard.png"
    cv2.imwrite(str(path), img)
    return path


def draw_attachment_stack(tables: dict[str, pd.DataFrame]) -> Path:
    img = np.full((1500, 2100, 3), BG, dtype=np.uint8)
    title(img, "Presub attachment stack", "What The Editor Should See First", "A high-impact presub should be short, visual, and traceable.")
    for i, row in tables["attachment"].iterrows():
        yy = 330 + i * 210
        color = [GOLD, GREEN, TEAL, BLUE, RED][i]
        box(img, (120, yy, 1980, yy + 150), PANEL if i % 2 == 0 else PANEL2)
        cv2.circle(img, (175, yy + 75), 34, color, -1, cv2.LINE_AA)
        put(img, row.slot, (164, yy + 86), 0.72, BG, 2)
        put(img, row.attachment, (240, yy + 58), 0.66, color, 2, 470)
        put(img, row.path_or_url, (760, yy + 48), 0.44, MUTED, 1, 500)
        put(img, row.purpose, (1320, yy + 58), 0.52, INK, 1, 560)
    path = OUT / "F04_presub_attachment_stack.png"
    cv2.imwrite(str(path), img)
    return path


def draw_k2_unlock(tables: dict[str, pd.DataFrame]) -> Path:
    img = np.full((1500, 2100, 3), BG, dtype=np.uint8)
    title(img, "K2/Bundang unlock board", "The One Dataset That Raises The Ceiling", "This is the practical action board for moving from strong CRM to plausible NC/LDH presub.")
    for i, row in tables["validation"].iterrows():
        yy = 335 + i * 210
        color = GREEN if row.priority == "highest" else GOLD
        box(img, (90, yy, 2010, yy + 148), PANEL if i % 2 == 0 else PANEL2)
        put(img, row.step, (125, yy + 62), 0.72, color, 2)
        put(img, row.task, (230, yy + 55), 0.66, INK, 2, 450)
        put(img, row.data, (730, yy + 55), 0.54, color, 2, 310)
        put(img, row.why_it_raises_impact, (1090, yy + 54), 0.52, MUTED, 1, 780)
    path = OUT / "F05_k2_bundang_unlock_board.png"
    cv2.imwrite(str(path), img)
    return path


def draw_nc_precedent_bridge(tables: dict[str, pd.DataFrame]) -> Path:
    img = np.full((1700, 2200, 3), BG, dtype=np.uint8)
    title(
        img,
        "Nature Communications precedent bridge",
        "Pathology Tissue -> Spatial Molecular Biology Is Already NC-Level",
        "Use Li et al. 2026 as an adjacent precedent: FFPE/pathology tissue can unlock spatial molecular states; Paper 2 makes H&E the computational front door.",
    )
    left = (80, 340, 985, 890)
    right = (1115, 340, 2020, 890)
    box(img, left, PANEL2, GOLD)
    box(img, right, PANEL2, TEAL)
    put(img, "Li et al. Nat Commun 2026", (120, 405), 0.78, GOLD, 2)
    put(img, "epi-Patho-DBiT on FFPE lymphoma", (120, 465), 0.58, INK, 2, 780)
    put(
        img,
        "Direct spatial epigenomics: spatial-FFPE-ATAC and spatial-FFPE-CUT&Tag decode chromatin and histone-modification landscapes from archival pathology tissue.",
        (120, 535),
        0.52,
        MUTED,
        1,
        780,
    )
    put(img, "High-impact precedent: pathology tissue can be a spatial molecular discovery substrate.", (120, 725), 0.54, GOLD, 2, 780)

    put(img, "Paper 2", (1155, 405), 0.78, TEAL, 2)
    put(img, "H&E-to-spatial-RNA biomarker recovery", (1155, 465), 0.58, INK, 2, 780)
    put(
        img,
        "Computational front door: routine H&E plus foundation morphology features recover a DM1/RAI-linked spatial RNA state in thyroid cancer.",
        (1155, 535),
        0.52,
        MUTED,
        1,
        780,
    )
    put(img, "Complementary claim: cheap triage into molecular/spatial states, not direct epigenome measurement.", (1155, 725), 0.54, TEAL, 2, 780)

    cv2.arrowedLine(img, (985, 615), (1115, 615), GOLD, 4, cv2.LINE_AA, tipLength=0.12)
    put(img, "extends the logic", (985, 575), 0.44, GOLD, 1)

    rows = [
        ("Shared logic", "Pathology material can reveal spatial molecular organization."),
        ("Paper 2 upgrade", "H&E becomes a general, low-cost entry point into molecular/spatial state prediction."),
        ("Boundary", "Do not claim direct measurement, pan-disease validation, or equivalence to ATAC/CUT&Tag."),
        ("Where to use", "Figure 1 precedent box; Introduction factual background; Discussion comparison slot."),
    ]
    y = 1040
    for i, (head, body) in enumerate(rows):
        yy = y + i * 135
        color = [GREEN, TEAL, RED, GOLD][i]
        box(img, (120, yy, 1980, yy + 95), PANEL if i % 2 == 0 else PANEL2)
        put(img, head, (155, yy + 55), 0.62, color, 2, 320)
        put(img, body, (520, yy + 55), 0.52, INK if i != 2 else RED, 1, 1320)

    path = OUT / "F06_natcom_2026_precedent_bridge.png"
    cv2.imwrite(str(path), img)
    return path


def write_outputs(tables: dict[str, pd.DataFrame], figs: list[Path], m: dict) -> None:
    for name, df in tables.items():
        df.to_csv(OUT / f"paper2_moonshot_{name}.tsv", sep="\t", index=False)
    summary = {
        "created": "2026-05-10",
        "purpose": "editor-facing moonshot package for Paper 2",
        "moonshot_thesis": "Routine H&E can act as a general digital front-door assay for latent molecular/spatial biomarker states; thyroid DM1/RAI is the validated exemplar.",
        "highest_impact_unlock": "K2/Bundang external H&E validation",
        "nature_communications_precedent": {
            "article": "Li et al. 2026, Spatially decoding genotype-associated epigenetic landscapes in human lymphoma FFPE tissues via epi-Patho-DBiT",
            "doi": "10.1038/s41467-026-71576-9",
            "use": "Adjacent precedent for pathology/FFPE tissue as a spatial molecular discovery substrate.",
            "boundary": "Not an H&E-to-RNA AI benchmark and not thyroid cancer.",
        },
        "headline_metrics": {
            "uni_loto_auc": m["uni"]["pooled_overall_auc"],
            "gse250521_slide_centered_rho": m["g250_stage"]["dm1_existing_slide_centered_rho"],
            "gse230424_he_rho": m["g230"]["top_he_sample_centered_rho"],
            "gse230424_residual_target_rho": m["g230_resid"]["dm1_coord_qc_residual_target_he_rho"],
            "gse230424_hotspot_lift": m["hotspot"]["gse230424_top10_lift"],
        },
        "figures": [p.name for p in figs],
    }
    (OUT / "PAPER2_MOONSHOT_EDITORIAL_PACK.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    md = [
        "# Paper 2 moonshot editorial pack",
        "",
        "Purpose: make Paper 2 read as a general H&E-to-molecular-state recovery architecture, not a small thyroid-only classifier.",
        "",
        "Highest-impact unlock: K2/Bundang external H&E validation.",
        "",
        "Figures:",
    ]
    for fig in figs:
        md.append(f"- `{fig.name}`")
    md += ["", "Voice boundary: this is figure/table/scaffold work only, not final cover-letter prose."]
    (OUT / "SUMMARY.md").write_text("\n".join(md) + "\n", encoding="utf-8")


def publish(figs: list[Path]) -> None:
    ASSET_LOCAL.mkdir(parents=True, exist_ok=True)
    ASSET_LIVE.mkdir(parents=True, exist_ok=True)
    for fig in figs:
        shutil.copy2(fig, ASSET_LOCAL / fig.name)
        shutil.copy2(fig, ASSET_LIVE / fig.name)


def html_table(df: pd.DataFrame) -> str:
    out = ["<table><thead><tr>" + "".join(f"<th>{c}</th>" for c in df.columns) + "</tr></thead><tbody>"]
    for _, row in df.iterrows():
        out.append("<tr>" + "".join(f"<td>{row[c]}</td>" for c in df.columns) + "</tr>")
    out.append("</tbody></table>")
    return "\n".join(out)


def write_page(tables: dict[str, pd.DataFrame], figs: list[Path], m: dict) -> None:
    fig_buttons = "\n".join(
        f'<button class="fig-option" data-src="assets/{ASSET}/{fig.name}" data-title="{fig.stem}"><img src="assets/{ASSET}/{fig.name}" alt=""><span>{fig.stem.replace("_", " ")}</span></button>'
        for fig in figs
    )
    html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Paper 2 Moonshot Editorial Pack</title>
<style>
:root {{ --bg:#0d1117; --panel:#151b23; --panel2:#1c2634; --ink:#e6edf3; --muted:#98a6ba; --line:#303846; --gold:#d6b25e; --blue:#57a6ff; }}
* {{ box-sizing:border-box; }}
body {{ margin:0; background:var(--bg); color:var(--ink); font-family:"JetBrains Mono", ui-monospace, Menlo, monospace; line-height:1.55; }}
a {{ color:var(--blue); text-decoration:none; }}
a:hover {{ text-decoration:underline; }}
.hero {{ min-height:64vh; display:grid; align-items:end; padding:58px 5vw 42px; border-bottom:1px solid var(--line); background:linear-gradient(180deg, rgba(13,17,23,.12), rgba(13,17,23,.96)), url("assets/{ASSET}/{figs[0].name}") center / cover no-repeat; }}
.kicker {{ color:var(--gold); text-transform:uppercase; font-size:13px; font-weight:900; }}
h1 {{ margin:10px 0 14px; max-width:1180px; font-family:"Cormorant Garamond", Georgia, serif; font-size:clamp(42px,7vw,88px); line-height:.96; letter-spacing:0; }}
.lead {{ max-width:1060px; color:#d8dee8; font-size:18px; }}
.stats {{ display:grid; grid-template-columns:repeat(5,minmax(0,1fr)); gap:10px; max-width:1200px; margin-top:24px; }}
.stat {{ border:1px solid rgba(214,178,94,.45); background:rgba(13,17,23,.84); padding:14px; min-height:92px; }}
.stat .v {{ color:var(--gold); font-size:23px; font-weight:900; }}
.stat .l {{ color:var(--muted); font-size:12px; margin-top:5px; }}
.layout {{ display:grid; grid-template-columns:280px minmax(0,1fr); gap:32px; max-width:1540px; margin:0 auto; padding:34px 24px 84px; }}
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
table {{ width:100%; border-collapse:collapse; margin:16px 0; font-size:13px; }}
th,td {{ border-bottom:1px solid var(--line); padding:10px 8px; text-align:left; vertical-align:top; }}
th {{ color:var(--gold); background:rgba(214,178,94,.06); }}
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
    <div class="kicker">Paper 2 · moonshot editorial pack · 2026-05-10</div>
    <h1>Make H&E The General Front Door, With Thyroid As The Proof</h1>
    <p class="lead">An editor-facing package: general H&E framing, graphical abstract, claim ladder, editorial scorecard, presub attachment stack, and K2/Bundang unlock board.</p>
    <div class="stats">
      <div class="stat"><div class="v">{m['uni']['pooled_overall_auc']:.3f}</div><div class="l">UNI LOTO AUC</div></div>
      <div class="stat"><div class="v">{m['g250_stage']['dm1_existing_slide_centered_rho']:.3f}</div><div class="l">GSE250521 rho</div></div>
      <div class="stat"><div class="v">{m['g230']['top_he_sample_centered_rho']:.3f}</div><div class="l">GSE230424 rho</div></div>
      <div class="stat"><div class="v">{m['hotspot']['gse230424_top10_lift']:.1f}x</div><div class="l">external hotspot lift</div></div>
      <div class="stat"><div class="v">K2</div><div class="l">single biggest unlock</div></div>
    </div>
  </div>
</header>
<div class="layout">
<nav>
  <a href="#verdict">01 Verdict</a>
  <a href="#figures">02 Figures</a>
  <a href="#generalization">03 H&amp;E Generalization</a>
  <a href="#precedent">04 NC Precedent</a>
  <a href="#claims">05 Claim Ladder</a>
  <a href="#editorial">06 Editor Scorecard</a>
  <a href="#attachments">07 Attachments</a>
  <a href="#validation">08 K2 Unlock</a>
  <a href="#presub">09 Presub Scaffold</a>
  <a href="#sources">10 Sources</a>
</nav>
<main>
<section id="verdict">
  <h2><span class="num">01</span>Moonshot Verdict</h2>
  <div class="grid">
    <div class="card"><h3>Best Big Claim</h3><p class="muted">Routine H&E can be a general front-door assay for latent molecular/spatial biomarker states, validated here in thyroid DM1/RAI.</p></div>
    <div class="card"><h3>Biggest Unlock</h3><p class="muted">K2/Bundang H&E external validation raises Nat Commun and Lancet Digital Health credibility simultaneously.</p></div>
    <div class="card"><h3>Guardrail</h3><p class="muted">Do not claim deployment-ready clinical decision support or causal MAPK/RAI mechanism from image data alone.</p></div>
  </div>
</section>
<section id="figures">
  <h2><span class="num">02</span>Moonshot Figure Viewer</h2>
  <div class="viewer">
    <div class="fig-list" id="figList">{fig_buttons}</div>
    <div class="stage">
      <div class="toolbar"><strong id="figTitle">Figure</strong><button class="control" id="zoomOut">Zoom out</button><button class="control" id="zoomReset">Reset</button><button class="control" id="zoomIn">Zoom in</button></div>
      <div class="canvas" id="canvas"><img id="mainFig" src="assets/{ASSET}/{figs[0].name}" alt="Paper 2 moonshot figure"></div>
    </div>
  </div>
</section>
<section id="generalization"><h2><span class="num">03</span>H&amp;E Generalization Frame</h2>{html_table(tables['generalization'])}</section>
<section id="precedent"><h2><span class="num">04</span>Nature Communications Precedent</h2>{html_table(tables['precedent'])}</section>
<section id="claims"><h2><span class="num">05</span>Claim Ladder</h2>{html_table(tables['claims'])}</section>
<section id="editorial"><h2><span class="num">06</span>Editor Scorecard</h2>{html_table(tables['editorial'])}</section>
<section id="attachments"><h2><span class="num">07</span>Presub Attachment Stack</h2>{html_table(tables['attachment'])}</section>
<section id="validation"><h2><span class="num">08</span>K2/Bundang Validation Unlock</h2>{html_table(tables['validation'])}</section>
<section id="presub"><h2><span class="num">09</span>Presubmission Fact Scaffold</h2>{html_table(tables['presub'])}</section>
<section id="sources">
  <h2><span class="num">10</span>Local Sources</h2>
  <table><tbody>
    <tr><td>Output directory</td><td><code>{OUT.relative_to(ROOT)}</code></td></tr>
    <tr><td>Builder script</td><td><code>project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/scripts/build_paper2_moonshot_editorial_pack.py</code></td></tr>
    <tr><td>Submission overview</td><td><a href="paper2_submission_overview.html">paper2_submission_overview.html</a></td></tr>
    <tr><td>High-impact journal route</td><td><a href="paper2_high_impact_journal_strategy.html">paper2_high_impact_journal_strategy.html</a></td></tr>
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
    PAGE_LOCAL.write_text(html, encoding="utf-8")
    shutil.copy2(PAGE_LOCAL, PAGE_LIVE)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    m = metrics()
    tables = make_tables(m)
    figs = [
        draw_he_generalization_frame(tables),
        draw_graphical_abstract(m),
        draw_claim_ladder(tables),
        draw_editor_scorecard(tables),
        draw_attachment_stack(tables),
        draw_k2_unlock(tables),
        draw_nc_precedent_bridge(tables),
    ]
    write_outputs(tables, figs, m)
    publish(figs)
    write_page(tables, figs, m)
    print(json.dumps({"out": str(OUT), "page": str(PAGE_LOCAL), "figures": [p.name for p in figs]}, indent=2))


if __name__ == "__main__":
    main()
