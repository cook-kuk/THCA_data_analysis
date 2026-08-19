#!/usr/bin/env python3
"""Build Paper 2 submission overview, figure roadmap, and data-use dossier."""
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
OUT = SUPP / "paper2_submission_overview_2026_05_10"
HUB = ROOT / "project/papers_hub_2026_05_04"
LIVE = Path("/var/www/papers/papers_hub_2026_05_04")
ASSET = "paper2_submission_overview"
ASSET_LOCAL = HUB / "assets" / ASSET
ASSET_LIVE = LIVE / "assets" / ASSET
PAGE_LOCAL = HUB / "paper2_submission_overview.html"
PAGE_LIVE = LIVE / "paper2_submission_overview.html"

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
        "g250": read_json(SUPP / "path2space_inspired_reanalysis_2026_05_09/PATH2SPACE_REANALYSIS_SUMMARY.json"),
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
        "downsample": read_json(
            SUPP / "path2space_downsample_robustness_2026_05_09/PATH2SPACE_DOWNSAMPLE_ROBUSTNESS_SUMMARY.json"
        ),
        "hotspot": read_json(SUPP / "path2space_hotspot_concordance_2026_05_09/PATH2SPACE_HOTSPOT_CONCORDANCE_SUMMARY.json"),
        "aitd": read_json(SUPP / "gse248205_pathology_aitd_axis_2026_05_09/GSE248205_PATHOLOGY_AITD_SUMMARY.json"),
    }


def put(
    img: np.ndarray,
    text: str,
    xy: tuple[int, int],
    scale: float = 0.55,
    color: tuple[int, int, int] = INK,
    thick: int = 1,
    width: int | None = None,
    gap: int = 8,
) -> int:
    x, y = xy
    font = cv2.FONT_HERSHEY_SIMPLEX
    lines: list[str] = []
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


def box(img: np.ndarray, xyxy: tuple[int, int, int, int], color: tuple[int, int, int] = PANEL, outline: tuple[int, int, int] = LINE) -> None:
    x0, y0, x1, y1 = xyxy
    cv2.rectangle(img, (x0, y0), (x1, y1), color, -1)
    cv2.rectangle(img, (x0, y0), (x1, y1), outline, 1, cv2.LINE_AA)


def title(img: np.ndarray, kicker: str, heading: str, sub: str) -> None:
    put(img, kicker.upper(), (70, 78), 0.58, GOLD, 2)
    put(img, heading, (70, 142), 1.12, INK, 2, 1600)
    put(img, sub, (70, 205), 0.56, MUTED, 1, 1500)


def fit_image(path: Path, size: tuple[int, int]) -> np.ndarray:
    tw, th = size
    src = cv2.imread(str(path))
    if src is None:
        canvas = np.full((th, tw, 3), PANEL2, dtype=np.uint8)
        put(canvas, f"missing\n{path.name}", (20, th // 2), 0.46, RED, 1, tw - 40)
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
    datasets = pd.DataFrame(
        [
            [
                "TCGA-THCA H&E WSI",
                "GDC / TCGA",
                "54 UNI LOTO slides; 59 ViT-L/CLAM final classification slides",
                "Primary image-DM1 classifier and attention heatmaps",
                "Main Fig 2; Main Fig 5; Table 1; Table S1-S2",
                "slide_manifest.tsv; clam_per_slide_predictions.tsv; UNI_LOTO_SUMMARY.json",
            ],
            [
                "TCGA-THCA RNA labels",
                "Paper 1 DM1/DM2 8-gene labels",
                "DM1/DM2 label source for classifier",
                "Ground-truth molecular label for H&E model",
                "Main Fig 2; Table 1; Table S1",
                "phase2_tcga_clam/slide_manifest.tsv",
            ],
            [
                "GSE250521",
                "GEO / Lu 2023 thyroid Visium + H&E",
                f"{m['g250']['n_slides']} slides; {m['g250']['n_spots_with_embeddings']} spots",
                "Spatial transcriptomic recovery and stage generalization",
                "Main Fig 3; Main Fig 4; Table 1; Table S3-S5",
                "path2space_spot_predictions.tsv.gz; PATH2SPACE_REANALYSIS_SUMMARY.json",
            ],
            [
                "GSE230424",
                "GEO / Zhang 2024 thyroid Visium + raw H&E",
                f"{m['g230']['n_samples']} samples; {m['g230']['n_spots']} spots",
                "External thyroid H&E-to-DM1/RAI support plus residual/QC caveat",
                "Main Fig 4; Main Fig 5; Table 1; Table S4-S6",
                "gse230424_pathology_predictions.tsv.gz; GSE230424_PATHOLOGY_THYROID_SUMMARY.json",
            ],
            [
                "GSE248205",
                "GEO / AITD spatial cohort",
                f"{m['aitd']['n_samples']} samples; {m['aitd']['n_spots']} spots",
                "Negative-control/caveat lock for local H&E immune-axis overclaim",
                "Main Fig 5; Table S6; reviewer response table",
                "GSE248205_PATHOLOGY_AITD_SUMMARY.json",
            ],
            [
                "Korean K2 / Bundang H&E",
                "PRJEB11591 / future institutional",
                "Not yet used as H&E validation",
                "Highest-impact revision unlock and Nat Commun ceiling raiser",
                "Main Fig 6 roadmap; Table 4; presubmission plan",
                "PAPER2_KOREAN_K2_VALIDATION_PLAN_2026_05_08.md",
            ],
        ],
        columns=["dataset", "source", "sample_size", "paper_role", "where_used", "key_files"],
    )

    figures = pd.DataFrame(
        [
            [
                "Fig 1",
                "Whole-study overview and clinical translation path",
                "H&E -> foundation model -> DM1/RAI spatial state -> reflex molecular testing",
                "TCGA, GSE250521, GSE230424, GSE248205",
                "Conceptual workflow plus actual cohort boxes",
                "Lead figure for Cell Reports Medicine / Communications Medicine",
            ],
            [
                "Fig 2",
                "Image-DM1 classifier performance",
                f"UNI LOTO AUC {m['uni']['pooled_overall_auc']:.3f}; RAS-like AUC {m['uni']['pooled_ras_like_auc']:.3f}",
                "TCGA-THCA H&E WSI + DM1/DM2 labels",
                "ROC, LOTO bars, attention heatmap examples",
                "Core model performance figure",
            ],
            [
                "Fig 3",
                "GSE250521 spatial RNA recovery",
                f"slide-centered rho {m['g250_stage']['dm1_existing_slide_centered_rho']:.3f}; PTC rho {m['g250_stage']['dm1_existing_ptc_rho']:.3f}; LPTC rho {m['g250_stage']['dm1_existing_lptc_rho']:.3f}",
                "GSE250521 Visium + H&E",
                "16 separated slide overlays, LOSO/stage-centered summary",
                "Turns classifier into cross-modal biomarker recovery",
            ],
            [
                "Fig 4",
                "External thyroid validation in GSE230424",
                f"H&E rho {m['g230']['top_he_sample_centered_rho']:.3f}; residual-target rho {m['g230_resid']['dm1_coord_qc_residual_target_he_rho']:.3f}",
                "GSE230424 raw H&E + Visium",
                "H&E overlay, residual model, sample context",
                "Independent thyroid support; caveat co-located",
            ],
            [
                "Fig 5",
                "Robustness, hotspot, dose-response, and caveat moat",
                f"hotspot lift {m['hotspot']['gse250521_top10_lift']:.1f}x/{m['hotspot']['gse230424_top10_lift']:.1f}x; decile rho {m['decile']['gse250521_raw_decile_spearman']:.3f}/{m['decile']['gse230424_raw_decile_spearman']:.3f}",
                "GSE250521, GSE230424, GSE248205",
                "Hotspot, decile, downsample, random-module, AITD no-go",
                "Reviewer-proofing figure",
            ],
            [
                "Fig 6",
                "Impact and validation unlock map",
                "Cell Reports Medicine primary; K2/Bundang H&E unlock for Nat Commun attempt",
                "All cohorts plus future validation plan",
                "Venue ladder, risk register, external validation roadmap",
                "Submission strategy and editor-facing close",
            ],
        ],
        columns=["figure", "title", "headline_message", "data_used", "panels", "submission_use"],
    )

    tables = pd.DataFrame(
        [
            ["Table 1", "Cohort and data inventory", "Every cohort, modality, n, role, and use location", "Required transparency table"],
            ["Table 2", "Model and validation performance", "AUC/rho/hotspot/dose-response metrics by analysis layer", "Main quantitative summary"],
            ["Table 3", "Reviewer caveat and mitigation matrix", "QC, smoothness, small-N, calibration, external validation gaps", "Preempts reviewer attacks"],
            ["Table 4", "Journal/validation strategy", "Venue ladder and K2/Bundang external validation unlock", "Submission planning table"],
            ["Table S1", "TCGA slide manifest", "file_id, labels, fold, prediction, tissue-site variables", "Reproducibility"],
            ["Table S2", "Per-slide classifier predictions", "OOF predictions, fold summary, calibration", "Model audit"],
            ["Table S3", "GSE250521 spot-level predictions", "spot-level observed/predicted DM1/RAI, slide/stage annotations", "Spatial reproducibility"],
            ["Table S4", "GSE230424 residual/QC models", "H&E, QC-only, coord+QC, residual-target metrics", "Caveat transparency"],
            ["Table S5", "Hotspot/dose-response/downsample controls", "Top-decile lift, decile gradients, 200-spot robustness", "Robustness audit"],
            ["Table S6", "Negative controls", "GSE248205 AITD, smoothness-adjusted failures, no-go statements", "Boundary lock"],
        ],
        columns=["table", "title", "contents", "purpose"],
    )

    flow = pd.DataFrame(
        [
            ["Step 1", "Define molecular target", "Paper 1 DM1/DM2 8-gene state", "TCGA RNA labels", "Classifier label"],
            ["Step 2", "Train pathology model", "H&E WSI tile features + MIL", "TCGA H&E", "AUC/attention"],
            ["Step 3", "Project into spatial RNA", "Path2Space-inspired LOSO prediction", "GSE250521", "rho, stage controls"],
            ["Step 4", "External thyroid validation", "Raw H&E + Visium support", "GSE230424", "rho, residual, hotspots"],
            ["Step 5", "Guard against overclaim", "Negative/caveat controls", "GSE248205 + smoothness/random modules", "claim boundary"],
            ["Step 6", "Raise ceiling", "Korean H&E validation plan", "K2/Bundang", "revision unlock"],
        ],
        columns=["step", "action", "analysis", "data", "output"],
    )

    captions = pd.DataFrame(
        [
            ["Fig 1", "Overview", "Schematic of the Paper 2 cross-modal biomarker-recovery workflow and cohort roles.", "Conceptual; no numerical claim beyond cohort roles."],
            ["Fig 2", "Classifier", "Out-of-fold H&E image-DM1 performance in TCGA-THCA with attention examples.", "Model performance, not deployment-ready calibration."],
            ["Fig 3", "Spatial recovery", "GSE250521 slide-level and stage-controlled spatial RNA recovery of DM1/RAI state.", "Strong in PTC/LPTC; ATC weaker."],
            ["Fig 4", "External validation", "GSE230424 H&E-to-DM1/RAI validation on independent thyroid Visium H&E.", "Support with QC/tissue-density caveat."],
            ["Fig 5", "Robustness/caveat", "Hotspot, decile, downsample, random-module, smoothness, and AITD negative-control checks.", "Bounded support; not causal mechanism proof."],
            ["Fig 6", "Impact roadmap", "Venue ladder and K2/Bundang H&E external validation unlock map.", "Submission strategy; not data result."],
        ],
        columns=["figure", "short_label", "caption_scaffold", "claim_boundary"],
    )
    return {"datasets": datasets, "figures": figures, "tables": tables, "flow": flow, "captions": captions}


def draw_flow_figure(tables: dict[str, pd.DataFrame], m: dict) -> Path:
    img = np.full((1320, 2100, 3), BG, dtype=np.uint8)
    title(img, "Paper 2 submission overview", "One Story: H&E Recovers A Spatial DM1/RAI Biomarker State", "Use this as the master figure-flow map for manuscript organization and advisor review.")
    steps = tables["flow"]
    y = 320
    x0 = 80
    w = 315
    gap = 32
    for i, row in steps.iterrows():
        x = x0 + (i % 3) * (w + gap)
        yy = y + (i // 3) * 320
        color = [GOLD, GREEN, TEAL, BLUE, GOLD, RED][i]
        box(img, (x, yy, x + w, yy + 245), PANEL if i % 2 == 0 else PANEL2)
        cv2.rectangle(img, (x, yy), (x + w, yy + 8), color, -1)
        put(img, row.step, (x + 18, yy + 42), 0.48, color, 2)
        put(img, row.action, (x + 18, yy + 84), 0.62, INK, 2, w - 36)
        put(img, row.analysis, (x + 18, yy + 132), 0.45, MUTED, 1, w - 36)
        put(img, f"Data: {row.data}", (x + 18, yy + 188), 0.43, INK, 1, w - 36)
        put(img, f"Output: {row.output}", (x + 18, yy + 222), 0.43, color, 1, w - 36)
        if i < len(steps) - 1 and i % 3 != 2:
            cv2.arrowedLine(img, (x + w + 4, yy + 122), (x + w + gap - 8, yy + 122), color, 2, cv2.LINE_AA, tipLength=0.28)
    y2 = 1000
    put(img, "Headline numbers to keep visible", (80, y2), 0.8, GOLD, 2)
    cards = [
        ("Classifier", f"AUC {m['uni']['pooled_overall_auc']:.3f}", "UNI leave-one-TSS-out"),
        ("Spatial", f"rho {m['g250_stage']['dm1_existing_slide_centered_rho']:.3f}", "GSE250521 slide-centered"),
        ("External", f"rho {m['g230']['top_he_sample_centered_rho']:.3f}", "GSE230424 H&E"),
        ("Hotspot", f"{m['hotspot']['gse230424_top10_lift']:.1f}x", "GSE230424 top-decile lift"),
    ]
    for i, (head, val, note) in enumerate(cards):
        x = 80 + i * 500
        box(img, (x, y2 + 65, x + 455, y2 + 225), PANEL2)
        put(img, head.upper(), (x + 22, y2 + 103), 0.46, MUTED, 1)
        put(img, val, (x + 22, y2 + 160), 1.0, GOLD if i != 2 else GREEN, 2)
        put(img, note, (x + 22, y2 + 198), 0.43, INK, 1, 400)
    path = OUT / "F00_paper2_whole_story_overview.png"
    cv2.imwrite(str(path), img)
    return path


def draw_data_use_map(tables: dict[str, pd.DataFrame]) -> Path:
    img = np.full((1450, 2200, 3), BG, dtype=np.uint8)
    title(img, "Dataset provenance", "Which Data Was Used Where", "Every major dataset is mapped to figure, table, and claim use so the submission package stays traceable.")
    cols = [330, 310, 430, 450, 520]
    heads = ["Dataset", "Sample size", "Paper role", "Where used", "Key files"]
    x0, y0 = 55, 310
    x = x0
    for w, h in zip(cols, heads):
        box(img, (x, y0, x + w, y0 + 62), (42, 34, 22), GOLD)
        put(img, h, (x + 12, y0 + 40), 0.48, GOLD, 2)
        x += w
    for r, row in tables["datasets"].iterrows():
        yy = y0 + 62 + r * 168
        x = x0
        vals = [row.dataset, row.sample_size, row.paper_role, row.where_used, row.key_files]
        for c, (w, val) in enumerate(zip(cols, vals)):
            box(img, (x, yy, x + w, yy + 168), PANEL if r % 2 == 0 else PANEL2)
            put(img, val, (x + 12, yy + 35), 0.43, INK if c != 4 else MUTED, 1, w - 24)
            x += w
    path = OUT / "F01_paper2_data_use_map.png"
    cv2.imwrite(str(path), img)
    return path


def draw_roadmap(tables: dict[str, pd.DataFrame]) -> Path:
    img = np.full((1320, 2200, 3), BG, dtype=np.uint8)
    title(img, "Figure and table roadmap", "Main Figures, Supplement Figures, And Tables", "A submission-ready map of what each item should contain and why it exists.")
    figures = tables["figures"]
    y = 300
    for i, row in figures.iterrows():
        yy = y + i * 154
        color = GREEN if i in {1, 2, 3} else GOLD if i in {0, 4} else BLUE
        box(img, (70, yy, 2030, yy + 118), PANEL if i % 2 == 0 else PANEL2)
        cv2.rectangle(img, (70, yy), (78, yy + 118), color, -1)
        put(img, row.figure, (100, yy + 42), 0.68, color, 2)
        put(img, row.title, (230, yy + 42), 0.62, INK, 2, 430)
        put(img, row.headline_message, (720, yy + 35), 0.46, MUTED, 1, 620)
        put(img, row.submission_use, (1415, yy + 35), 0.46, INK, 1, 520)
    put(img, "Tables", (70, 1285), 0.85, GOLD, 2)
    for i, row in tables["tables"].head(6).iterrows():
        yy = 1360 + i * 95
        box(img, (70, yy, 2030, yy + 72), PANEL2 if i % 2 == 0 else PANEL)
        put(img, row.table, (100, yy + 43), 0.56, GOLD, 2)
        put(img, row.title, (250, yy + 43), 0.52, INK, 1, 430)
        put(img, row.contents, (740, yy + 36), 0.44, MUTED, 1, 760)
        put(img, row.purpose, (1570, yy + 36), 0.44, INK, 1, 380)
    path = OUT / "F02_paper2_figure_table_roadmap.png"
    cv2.imwrite(str(path), img)
    return path


def draw_representative_bundle() -> Path:
    img = np.full((1500, 2200, 3), BG, dtype=np.uint8)
    title(img, "Representative figure bundle", "Suggested Main-Figure Visual Flow", "Representative existing assets are arranged in the order an editor should understand the paper.")
    reps = [
        ("Fig 1 overview", HUB / "assets/paper2_impact_upgrade/F01_paper2_impact_ladder.png"),
        ("Fig 2 classifier", HUB / "assets/paper2_individual_figures/analysis_supp__audit_uni_loto__fig_H_uni_loto.png"),
        ("Fig 3 GSE250521 spatial", HUB / "assets/paper2_individual_figures/figures__fig10_spatial_grid_16.png"),
        ("Fig 4 GSE230424 external", HUB / "assets/paper2_cv2_summary/fig02_gse230424_cv2_he_overlay_mosaic.png"),
        ("Fig 5 robustness", HUB / "assets/paper2_impact_upgrade/F02_paper2_evidence_moat.png"),
        ("Fig 6 unlock", HUB / "assets/paper2_impact_upgrade/F04_paper2_validation_unlock_map.png"),
    ]
    for i, (label, path) in enumerate(reps):
        x = 70 + (i % 2) * 1035
        y = 315 + (i // 2) * 370
        box(img, (x, y, x + 975, y + 318), PANEL)
        thumb = fit_image(path, (440, 250))
        img[y + 48 : y + 298, x + 20 : x + 460] = thumb
        put(img, label, (x + 490, y + 75), 0.68, GOLD, 2, 430)
        put(img, f"Asset: {path.name}", (x + 490, y + 130), 0.42, MUTED, 1, 430)
        put(
            img,
            "Use as the representative visual anchor for this stage of the manuscript. Keep full source figure inspectable in the separated browser.",
            (x + 490, y + 178),
            0.45,
            INK,
            1,
            410,
        )
    path = OUT / "F03_paper2_representative_figure_bundle.png"
    cv2.imwrite(str(path), img)
    return path


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


def write_reports(tables: dict[str, pd.DataFrame], figs: list[Path], m: dict) -> None:
    for name, df in tables.items():
        df.to_csv(OUT / f"paper2_submission_{name}.tsv", sep="\t", index=False)
    summary = {
        "created": "2026-05-10",
        "purpose": "Paper 2 high-impact submission figure/table/data-use overview",
        "headline": {
            "uni_loto_auc": m["uni"]["pooled_overall_auc"],
            "gse250521_slide_centered_rho": m["g250_stage"]["dm1_existing_slide_centered_rho"],
            "gse230424_he_rho": m["g230"]["top_he_sample_centered_rho"],
            "gse230424_residual_target_rho": m["g230_resid"]["dm1_coord_qc_residual_target_he_rho"],
            "hotspot_lift_gse250521": m["hotspot"]["gse250521_top10_lift"],
            "hotspot_lift_gse230424": m["hotspot"]["gse230424_top10_lift"],
        },
        "figures": [p.name for p in figs],
        "tables": [f"paper2_submission_{name}.tsv" for name in tables],
    }
    (OUT / "PAPER2_SUBMISSION_OVERVIEW_SUMMARY.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    md = [
        "# Paper 2 submission overview pack",
        "",
        "This is scaffolding only: figure/table/data-use roadmap and caption skeletons. It does not write protected manuscript sections.",
        "",
        "## Representative main figure flow",
    ]
    for _, row in tables["figures"].iterrows():
        md.append(f"- {row.figure}: {row.title} -- {row.headline_message}")
    md += ["", "## Data-use inventory"]
    for _, row in tables["datasets"].iterrows():
        md.append(f"- {row.dataset}: {row.paper_role}; used in {row.where_used}.")
    md += ["", "## Files"]
    for fig in figs:
        md.append(f"- `{fig.name}`")
    (OUT / "SUMMARY.md").write_text("\n".join(md) + "\n", encoding="utf-8")


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
<title>Paper 2 Submission Overview</title>
<style>
:root {{ --bg:#0d1117; --panel:#151b23; --panel2:#1c2634; --ink:#e6edf3; --muted:#98a6ba; --line:#303846; --gold:#d6b25e; --teal:#63d5c4; --blue:#57a6ff; --red:#ff7b72; }}
* {{ box-sizing:border-box; }}
body {{ margin:0; background:var(--bg); color:var(--ink); font-family:"JetBrains Mono", ui-monospace, Menlo, monospace; line-height:1.55; }}
a {{ color:var(--blue); text-decoration:none; }}
a:hover {{ text-decoration:underline; }}
.hero {{ min-height:62vh; display:grid; align-items:end; padding:58px 5vw 42px; border-bottom:1px solid var(--line); background:linear-gradient(180deg, rgba(13,17,23,.16), rgba(13,17,23,.97)), url("assets/{ASSET}/{figs[0].name}") center / cover no-repeat; }}
.kicker {{ color:var(--gold); text-transform:uppercase; font-size:13px; font-weight:900; }}
h1 {{ margin:10px 0 14px; max-width:1180px; font-family:"Cormorant Garamond", Georgia, serif; font-size:clamp(42px,7vw,88px); line-height:.96; letter-spacing:0; }}
.lead {{ max-width:1060px; color:#d8dee8; font-size:18px; }}
.stats {{ display:grid; grid-template-columns:repeat(6,minmax(0,1fr)); gap:10px; max-width:1240px; margin-top:24px; }}
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
.grid {{ display:grid; grid-template-columns:repeat(2,minmax(0,1fr)); gap:16px; }}
.card {{ border:1px solid var(--line); background:var(--panel); border-radius:8px; padding:18px; }}
.muted {{ color:var(--muted); }}
table {{ width:100%; border-collapse:collapse; margin:16px 0; font-size:13px; }}
th,td {{ border-bottom:1px solid var(--line); padding:10px 8px; text-align:left; vertical-align:top; }}
th {{ color:var(--gold); background:rgba(214,178,94,.06); }}
code {{ background:var(--panel2); border:1px solid var(--line); padding:2px 5px; border-radius:4px; color:#d8dee8; }}
.viewer {{ display:grid; grid-template-columns:370px minmax(0,1fr); gap:16px; }}
.fig-list {{ display:grid; gap:10px; align-content:start; max-height:760px; overflow:auto; padding-right:6px; }}
.fig-option {{ border:1px solid var(--line); background:var(--panel); color:var(--ink); border-radius:8px; padding:10px; display:grid; grid-template-columns:96px minmax(0,1fr); gap:10px; text-align:left; cursor:pointer; font:inherit; }}
.fig-option.active,.fig-option:hover {{ border-color:var(--gold); background:#202938; }}
.fig-option img {{ width:96px; height:66px; object-fit:cover; border-radius:4px; display:block; }}
.stage {{ border:1px solid var(--line); background:#090b10; border-radius:8px; overflow:hidden; }}
.toolbar {{ display:flex; gap:8px; align-items:center; padding:12px; background:var(--panel); border-bottom:1px solid var(--line); }}
.toolbar strong {{ flex:1 1 auto; }}
button.control {{ border:1px solid var(--line); background:var(--panel2); color:var(--ink); border-radius:6px; padding:8px 10px; font:inherit; font-size:12px; cursor:pointer; }}
button.control:hover {{ border-color:var(--gold); color:var(--gold); }}
.canvas {{ height:min(78vh,840px); overflow:auto; padding:16px; background:#090b10; }}
.canvas img {{ display:block; width:100%; min-width:640px; max-width:none; height:auto; margin:0 auto; }}
@media(max-width:980px) {{ .layout,.viewer {{ grid-template-columns:1fr; }} nav {{ position:static; max-height:none; display:grid; grid-template-columns:repeat(2,minmax(0,1fr)); }} .stats {{ grid-template-columns:repeat(2,minmax(0,1fr)); }} .grid {{ grid-template-columns:1fr; }} }}
</style>
</head>
<body>
<header class="hero">
  <div>
    <div class="kicker">Paper 2 · submission overview · 2026-05-10</div>
    <h1>Representative Figures, Tables, Flow, And Data Provenance</h1>
    <p class="lead">A high-impact submission scaffold: main figure flow, data-use map, table plan, caption skeletons, and exact links from each dataset to each figure/table/claim.</p>
    <div class="stats">
      <div class="stat"><div class="v">6</div><div class="l">main figure plan</div></div>
      <div class="stat"><div class="v">10</div><div class="l">table roadmap entries</div></div>
      <div class="stat"><div class="v">{m['uni']['pooled_overall_auc']:.3f}</div><div class="l">UNI LOTO AUC</div></div>
      <div class="stat"><div class="v">{m['g250_stage']['dm1_existing_slide_centered_rho']:.3f}</div><div class="l">GSE250521 rho</div></div>
      <div class="stat"><div class="v">{m['g230']['top_he_sample_centered_rho']:.3f}</div><div class="l">GSE230424 rho</div></div>
      <div class="stat"><div class="v">K2</div><div class="l">highest unlock</div></div>
    </div>
  </div>
</header>
<div class="layout">
<nav>
  <a href="#verdict">01 Verdict</a>
  <a href="#viewer">02 Overview Figures</a>
  <a href="#flow">03 Flow Steps</a>
  <a href="#figures">04 Figure Plan</a>
  <a href="#tables">05 Table Plan</a>
  <a href="#datasets">06 Data Use</a>
  <a href="#captions">07 Caption Skeletons</a>
  <a href="#sources">08 Sources</a>
</nav>
<main>
<section id="verdict">
  <h2><span class="num">01</span>Submission Package Verdict</h2>
  <div class="grid">
    <div class="card"><h3>Representative Story</h3><p class="muted">Lead with cross-modal biomarker recovery: H&E predicts a DM1/RAI-linked spatial RNA state, validated in TCGA, GSE250521, and GSE230424, with caveats visible.</p></div>
    <div class="card"><h3>Figure Logic</h3><p class="muted">Fig 1 explains the whole system; Fig 2 proves the classifier; Fig 3-4 prove spatial/external transfer; Fig 5 handles robustness/caveats; Fig 6 shows the high-impact validation unlock.</p></div>
  </div>
</section>
<section id="viewer">
  <h2><span class="num">02</span>Overview Figure Viewer</h2>
  <div class="viewer">
    <div class="fig-list" id="figList">{fig_buttons}</div>
    <div class="stage">
      <div class="toolbar"><strong id="figTitle">Figure</strong><button class="control" id="zoomOut">Zoom out</button><button class="control" id="zoomReset">Reset</button><button class="control" id="zoomIn">Zoom in</button></div>
      <div class="canvas" id="canvas"><img id="mainFig" src="assets/{ASSET}/{figs[0].name}" alt="Paper 2 overview figure"></div>
    </div>
  </div>
</section>
<section id="flow"><h2><span class="num">03</span>Whole-Story Flow</h2>{html_table(tables['flow'])}</section>
<section id="figures"><h2><span class="num">04</span>Main Figure Plan</h2>{html_table(tables['figures'])}</section>
<section id="tables"><h2><span class="num">05</span>Table Plan</h2>{html_table(tables['tables'])}</section>
<section id="datasets"><h2><span class="num">06</span>Data-Use Inventory</h2>{html_table(tables['datasets'])}</section>
<section id="captions"><h2><span class="num">07</span>Caption Skeletons And Claim Boundaries</h2>{html_table(tables['captions'])}</section>
<section id="sources">
  <h2><span class="num">08</span>Output Sources</h2>
  <table><tbody>
    <tr><td>Output directory</td><td><code>{OUT.relative_to(ROOT)}</code></td></tr>
    <tr><td>Builder script</td><td><code>project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/scripts/build_paper2_submission_overview_pack.py</code></td></tr>
    <tr><td>Separated figure browser</td><td><a href="paper2_cv2_visual_summary.html">paper2_cv2_visual_summary.html</a></td></tr>
    <tr><td>Impact dossier</td><td><a href="paper2_impact_upgrade.html">paper2_impact_upgrade.html</a></td></tr>
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
        draw_flow_figure(tables, m),
        draw_data_use_map(tables),
        draw_roadmap(tables),
        draw_representative_bundle(),
    ]
    write_reports(tables, figs, m)
    publish(figs)
    write_page(tables, figs, m)
    print(json.dumps({"out": str(OUT), "page": str(PAGE_LOCAL), "figures": [p.name for p in figs]}, indent=2))


if __name__ == "__main__":
    main()
