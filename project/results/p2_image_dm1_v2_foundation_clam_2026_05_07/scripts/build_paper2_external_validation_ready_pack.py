#!/usr/bin/env python3
"""Build a Paper 2 external-validation-ready packet.

Scaffold/infra/figure/table work only. This does not generate protected
manuscript voice sections or final collaborator email prose.
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
OUT = SUPP / "paper2_external_validation_ready_2026_05_10"
HUB = ROOT / "project/papers_hub_2026_05_04"
LIVE = Path("/var/www/papers/papers_hub_2026_05_04")
ASSET = "paper2_external_validation_ready"
ASSET_LOCAL = HUB / "assets" / ASSET
ASSET_LIVE = LIVE / "assets" / ASSET
PAGE_LOCAL = HUB / "paper2_external_validation_ready_pack.html"
PAGE_LIVE = LIVE / "paper2_external_validation_ready_pack.html"

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
    put(img, heading, (70, 142), 1.08, INK, 2, 1650)
    put(img, sub, (70, 205), 0.56, MUTED, 1, 1550)


def make_tables(m: dict) -> dict[str, pd.DataFrame]:
    request = pd.DataFrame(
        [
            ["Required", "De-identified H&E WSI", "SVS/NDPI/MRX/SCN/TIFF or pyramidal TIFF", "At least one diagnostic tumor block per case", "No PHI in filename or metadata"],
            ["Required", "Case/slide manifest", "TSV template provided", "case_id, slide_id, block_id, histology, scanner, stain, QC flags", "Use coded IDs only"],
            ["Required", "Primary molecular label", "Binary DM1-like status or continuous RAI/DM1 proxy", "RNA score, 8-gene/RAI score, or validated local proxy", "State missingness and assay platform"],
            ["Strongly preferred", "Driver/genotype labels", "BRAF, RAS, RET/NTRK/ALK fusion, TERT", "Use for subgroup and confounder audit", "Not required for primary image inference"],
            ["Preferred", "Pathology annotations", "Tumor ROI, lymphoid aggregate/TLS, necrosis, fibrosis, artifact", "Can be coarse polygon or slide-level", "Reader blinded to model output if possible"],
            ["Optional", "Matched RNA or spatial RNA", "bulk RNA-seq, NanoString, panel, Visium, GeoMx", "Enables continuous rho and residual controls", "Do not delay H&E-only validation"],
        ],
        columns=["priority", "data_item", "format", "minimum_content", "guardrail"],
    )
    manifest = pd.DataFrame(
        [
            ["case_id", "string", "required", "De-identified case key; no MRN/name/date of birth", "P2EXT_0001"],
            ["slide_id", "string", "required", "De-identified slide key; unique per WSI", "P2EXT_0001_A"],
            ["wsi_filename", "string", "required", "Actual de-identified image filename", "P2EXT_0001_A.svs"],
            ["block_id", "string", "preferred", "Coded FFPE block ID", "B01"],
            ["cohort_site", "string", "required", "K2, Bundang, or other site label", "Bundang"],
            ["scanner", "string", "required", "Scanner model if known", "Aperio_AT2"],
            ["magnification", "string", "required", "20x/40x or objective", "40x"],
            ["mpp", "float", "preferred", "Microns per pixel", "0.25"],
            ["stain_protocol", "string", "preferred", "H&E protocol/batch if known", "HE_batch_01"],
            ["histology", "string", "required", "PTC/FVPTC/FTC/PDTC/ATC/etc.", "PTC"],
            ["tumor_percent", "float", "preferred", "Estimated tumor cellularity percent", "65"],
            ["hashimoto_status", "string", "preferred", "present/absent/unknown", "present"],
            ["braf_status", "string", "preferred", "V600E/wildtype/unknown", "V600E"],
            ["ras_status", "string", "preferred", "mutant/wildtype/unknown", "wildtype"],
            ["tert_status", "string", "preferred", "mutant/wildtype/unknown", "unknown"],
            ["rai8_score", "float", "preferred", "Continuous RAI_8 or local equivalent", "-0.72"],
            ["dm1_label", "string", "preferred", "positive/negative/unknown", "positive"],
            ["rna_available", "string", "required", "yes/no", "yes"],
            ["pathology_qc_flag", "string", "required", "pass/fail/review", "pass"],
            ["exclude_reason", "string", "required_if_fail", "Reason if QC fail", ""],
        ],
        columns=["field", "type", "status", "definition", "example"],
    )
    sap = pd.DataFrame(
        [
            ["Freeze", "Record code hash, model files, feature extractor, preprocessing, and thresholds before external data are opened", "No retuning after cohort receipt"],
            ["QC", "Reject blank, out-of-focus, non-tumor, severe artifact, wrong stain, duplicate, or PHI-containing files", "Keep exclusion table"],
            ["Inference", "Run existing Paper 2 H&E model or UNI feature pipeline on all QC-pass slides", "No site-specific normalization beyond pre-specified stain/QC handling"],
            ["Primary binary endpoint", "External AUROC for DM1-like label if binary label exists", "Report CI and exact n"],
            ["Primary continuous endpoint", "Spearman rho with RAI_8/DM1 continuous score if RNA/proxy exists", "Report raw and covariate-adjusted variants"],
            ["Secondary hotspot endpoint", "Top-decile image-DM1 hotspot enrichment against pathology-reader regions", "Exploratory unless pre-registered"],
            ["Confounder audit", "Scanner/site/stain/histology/BRAF/RAS/TERT/HT strata where available", "Show failed strata instead of hiding them"],
            ["Decision", "Assign green/yellow/red route by locked thresholds", "No threshold movement after looking at data"],
        ],
        columns=["module", "locked_action", "guardrail"],
    )
    gates = pd.DataFrame(
        [
            ["Green", "n >= 80 QC-pass slides and AUROC >= 0.75 or rho >= 0.30", "Calibration not collapsed; site/stain residual remains directionally positive", "Nat Commun route becomes realistic; Lancet DH presub possible with workflow framing"],
            ["Yellow", "n >= 40 QC-pass slides and AUROC 0.65-0.75 or rho 0.20-0.30", "Signal survives at least one confounder audit", "Use as external support for Cell Reports Medicine / Communications Medicine"],
            ["Red", "AUROC < 0.65 and rho < 0.20 or signal reverses", "No post hoc rescue as main result", "Keep Paper 2 current route; report as limitation or hold for more data"],
            ["No-go", "QC-pass n < 25 or labels unusable", "Do not interpret performance", "Use only as feasibility/data-acquisition note"],
        ],
        columns=["decision", "minimum_performance", "required_context", "journal_action"],
    )
    reader = pd.DataFrame(
        [
            ["slide_id", "string", "Slide key from manifest", "required"],
            ["reader_id", "string", "Blinded pathology reader code", "required"],
            ["tumor_roi_available", "yes/no", "Tumor ROI marked or not", "required"],
            ["lymphoid_aggregate_tls", "0/1/2/3", "Absent/focal/moderate/extensive", "preferred"],
            ["fibrosis_sclerosis", "0/1/2/3", "Absent/focal/moderate/extensive", "preferred"],
            ["follicular_colloid_pattern", "0/1/2/3", "Absent/focal/moderate/extensive", "preferred"],
            ["solid_trabecular_growth", "0/1/2/3", "Absent/focal/moderate/extensive", "preferred"],
            ["necrosis_anaplasia", "0/1/2/3", "Absent/focal/moderate/extensive", "preferred"],
            ["artifact_or_low_quality", "0/1/2/3", "None/mild/moderate/severe", "required"],
            ["model_hotspot_overlap_rating", "0/1/2/3", "After model unblinding: none/weak/moderate/strong", "secondary"],
        ],
        columns=["field", "type", "definition", "status"],
    )
    deliverables = pd.DataFrame(
        [
            ["D0", "Receive de-identified manifest + WSI file list", "Completeness check; PHI audit", "Data receipt log"],
            ["D1", "Run WSI QC + model inference", "QC-pass count; tile/slide predictions", "QC table + prediction TSV"],
            ["D2", "Primary endpoint analysis", "AUROC/rho with confidence intervals", "External validation figure"],
            ["D3", "Confounder and calibration audit", "Scanner/stain/site/histology/driver strata", "Risk buy-down table"],
            ["D4", "Pathologist hotspot audit", "Reader overlap + morphology categories", "Hotspot plausibility panel"],
            ["D5", "Journal decision", "Green/yellow/red route assignment", "Updated Fig 6 / presub packet"],
        ],
        columns=["day", "task", "readout", "deliverable"],
    )
    routes = pd.DataFrame(
        [
            ["Nature Communications", "External H&E green result", "K2/Bundang locked validation plus current spatial RNA evidence", "Do not submit NC on references alone"],
            ["The Lancet Digital Health", "External H&E green + threshold/workflow", "Reflex-test triage scenario, calibration, reader audit", "Presub first only"],
            ["Cell Reports Medicine", "Current evidence or yellow external result", f"UNI AUC {m['uni']['pooled_overall_auc']:.3f}; spatial rho {m['g250_stage']['dm1_existing_slide_centered_rho']:.3f}/{m['g230']['top_he_sample_centered_rho']:.3f}", "Most realistic current route"],
            ["Communications Medicine", "Current evidence + reference/validation-ready packet", "Clinical translational framing with caveat moat", "Backup route"],
        ],
        columns=["journal", "trigger", "evidence_package", "boundary"],
    )
    return {
        "request": request,
        "manifest": manifest,
        "sap": sap,
        "gates": gates,
        "reader": reader,
        "deliverables": deliverables,
        "routes": routes,
    }


def draw_request_packet(tables: dict[str, pd.DataFrame]) -> Path:
    img = np.full((1700, 2300, 3), BG, dtype=np.uint8)
    title(img, "External validation request packet", "What We Need From K2/Bundang To Raise The Ceiling", "This converts the max-impact idea into a concrete data handoff.")
    rows = tables["request"]
    y = 330
    for i, row in rows.iterrows():
        yy = y + i * 195
        color = GOLD if row.priority == "Required" else TEAL if row.priority == "Strongly preferred" else BLUE
        box(img, (90, yy, 2210, yy + 145), PANEL if i % 2 == 0 else PANEL2, color)
        put(img, row.priority, (125, yy + 52), 0.5, color, 2, 300)
        put(img, row.data_item, (460, yy + 52), 0.56, INK, 2, 420)
        put(img, row.format, (930, yy + 45), 0.43, MUTED, 1, 410)
        put(img, row.minimum_content, (1385, yy + 45), 0.43, INK, 1, 430)
        put(img, row.guardrail, (1845, yy + 45), 0.41, RED if "PHI" in row.guardrail else GOLD, 1, 330)
    path = OUT / "F01_external_validation_request_packet.png"
    cv2.imwrite(str(path), img)
    return path


def draw_manifest_schema(tables: dict[str, pd.DataFrame]) -> Path:
    img = np.full((1800, 2400, 3), BG, dtype=np.uint8)
    title(img, "Sample manifest schema", "The File That Prevents Reviewer Chaos", "Every WSI needs enough metadata to audit scanner, stain, histology, driver, RNA, and QC effects.")
    fields = tables["manifest"]
    col_w = [310, 230, 230, 900, 380]
    heads = ["Field", "Type", "Status", "Definition", "Example"]
    y0 = 315
    x = 70
    for head, w in zip(heads, col_w):
        box(img, (x, y0, x + w, y0 + 62), (42, 34, 22), GOLD)
        put(img, head, (x + 12, y0 + 40), 0.43, GOLD, 2)
        x += w
    for i, row in fields.iterrows():
        yy = y0 + 62 + i * 70
        vals = [row.field, row.type, row.status, row.definition, row.example]
        x = 70
        for j, (w, val) in enumerate(zip(col_w, vals)):
            box(img, (x, yy, x + w, yy + 70), PANEL if i % 2 == 0 else PANEL2)
            color = GOLD if j == 0 else RED if val == "required" else INK
            put(img, val, (x + 10, yy + 39), 0.34, color, 1, w - 20, 4)
            x += w
    path = OUT / "F02_sample_manifest_schema.png"
    cv2.imwrite(str(path), img)
    return path


def draw_locked_workflow(tables: dict[str, pd.DataFrame]) -> Path:
    img = np.full((1550, 2300, 3), BG, dtype=np.uint8)
    title(img, "Locked statistical analysis plan", "No Retuning After External Data Arrive", "This is the reviewer-safe validation workflow.")
    rows = tables["sap"]
    for i, row in rows.iterrows():
        x = 100 + (i % 4) * 535
        y = 330 + (i // 4) * 480
        color = [GOLD, TEAL, GREEN, BLUE, BLUE, PURPLE, RED, GOLD][i]
        box(img, (x, y, x + 455, y + 350), PANEL if i % 2 == 0 else PANEL2, color)
        cv2.circle(img, (x + 45, y + 48), 28, color, -1, cv2.LINE_AA)
        put(img, str(i + 1), (x + 35, y + 59), 0.62, BG, 2)
        put(img, row.module, (x + 90, y + 58), 0.58, color, 2, 330)
        put(img, row.locked_action, (x + 28, y + 145), 0.44, INK, 1, 390)
        put(img, row.guardrail, (x + 28, y + 280), 0.4, RED if i < 3 else MUTED, 1, 390)
        if i < len(rows) - 1 and i % 4 != 3:
            cv2.arrowedLine(img, (x + 455, y + 170), (x + 535, y + 170), color, 3, cv2.LINE_AA, tipLength=0.12)
    box(img, (110, 1330, 2190, 1445), (42, 34, 22), GOLD)
    put(img, "Core rule: report the failed states. The credibility gain comes from locked analysis, not from perfect results.", (145, 1400), 0.58, GOLD, 2, 1900)
    path = OUT / "F03_locked_analysis_workflow.png"
    cv2.imwrite(str(path), img)
    return path


def draw_go_nogo(tables: dict[str, pd.DataFrame]) -> Path:
    img = np.full((1600, 2300, 3), BG, dtype=np.uint8)
    title(img, "Go/no-go thresholds", "Make The Journal Decision Before Seeing The Data", "Green raises the ceiling; yellow supports current route; red is reported honestly.")
    colors = {"Green": GREEN, "Yellow": GOLD, "Red": RED, "No-go": PURPLE}
    for i, row in tables["gates"].iterrows():
        yy = 340 + i * 285
        color = colors[row.decision]
        box(img, (110, yy, 2190, yy + 210), PANEL if i % 2 == 0 else PANEL2, color)
        put(img, row.decision, (150, yy + 70), 0.82, color, 2, 220)
        put(img, row.minimum_performance, (430, yy + 52), 0.5, INK, 1, 520)
        put(img, row.required_context, (1010, yy + 52), 0.48, MUTED, 1, 440)
        put(img, row.journal_action, (1510, yy + 52), 0.5, color, 2, 560)
    path = OUT / "F04_go_nogo_threshold_board.png"
    cv2.imwrite(str(path), img)
    return path


def draw_reader_audit(tables: dict[str, pd.DataFrame]) -> Path:
    img = np.full((1700, 2300, 3), BG, dtype=np.uint8)
    title(img, "Pathologist hotspot audit", "Turn Heatmaps Into Clinically Legible Morphology", "A blinded reader sheet reduces black-box concern without claiming clinical utility.")
    rows = tables["reader"]
    for i, row in rows.iterrows():
        yy = 330 + i * 120
        box(img, (100, yy, 2200, yy + 82), PANEL if i % 2 == 0 else PANEL2)
        put(img, row.field, (135, yy + 50), 0.48, GOLD if i < 3 else TEAL, 2, 430)
        put(img, row.type, (610, yy + 50), 0.42, MUTED, 1, 280)
        put(img, row.definition, (930, yy + 50), 0.44, INK, 1, 740)
        put(img, row.status, (1710, yy + 50), 0.44, RED if row.status == "required" else GOLD, 2, 330)
    path = OUT / "F05_pathologist_hotspot_audit_sheet.png"
    cv2.imwrite(str(path), img)
    return path


def draw_timeline(tables: dict[str, pd.DataFrame]) -> Path:
    img = np.full((1550, 2300, 3), BG, dtype=np.uint8)
    title(img, "Five-day execution board", "What Happens Once The Slides Arrive", "This makes the collaborator request concrete and fast.")
    rows = tables["deliverables"]
    for i, row in rows.iterrows():
        yy = 330 + i * 185
        color = [GOLD, TEAL, GREEN, BLUE, PURPLE, RED][i]
        box(img, (130, yy, 2170, yy + 125), PANEL if i % 2 == 0 else PANEL2, color)
        put(img, row.day, (170, yy + 65), 0.78, color, 2)
        put(img, row.task, (330, yy + 55), 0.58, INK, 2, 560)
        put(img, row.readout, (940, yy + 50), 0.46, MUTED, 1, 520)
        put(img, row.deliverable, (1510, yy + 55), 0.5, color, 2, 530)
    box(img, (130, 1520, 2170, 1630), (42, 34, 22), GOLD)
    put(img, "If the data are clean, the validation figure can be submission-ready within days, not weeks.", (165, 1588), 0.58, GOLD, 2, 1880)
    path = OUT / "F06_five_day_execution_board.png"
    cv2.imwrite(str(path), img)
    return path


def write_outputs(tables: dict[str, pd.DataFrame], figs: list[Path], m: dict) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for name, df in tables.items():
        df.to_csv(OUT / f"paper2_external_validation_{name}.tsv", sep="\t", index=False)
    collaborator_md = [
        "# Paper 2 external validation packet",
        "",
        "Purpose: define a de-identified H&E validation handoff for Paper 2.",
        "",
        "Required files:",
        "- De-identified diagnostic H&E WSI files.",
        "- Completed sample manifest template.",
        "- Primary molecular label or continuous RAI/DM1 proxy if available.",
        "",
        "Locked analysis boundary:",
        "- No model retuning after external data receipt.",
        "- All QC failures and missing labels are reported.",
        "- Green/yellow/red journal decision follows the pre-specified thresholds.",
        "",
        "Voice boundary: this is not collaborator email prose.",
    ]
    (OUT / "README_FOR_COLLABORATOR_HANDOFF.md").write_text("\n".join(collaborator_md) + "\n", encoding="utf-8")
    summary = {
        "created": "2026-05-10",
        "purpose": "external validation ready handoff for Paper 2",
        "top_unlock": "K2/Bundang de-identified H&E WSI plus locked manifest",
        "required_files": [
            "H&E WSI files",
            "paper2_external_validation_manifest.tsv",
            "primary molecular label or continuous RAI/DM1 proxy if available",
        ],
        "headline_metrics_context": {
            "uni_loto_auc": m["uni"]["pooled_overall_auc"],
            "gse250521_slide_centered_rho": m["g250_stage"]["dm1_existing_slide_centered_rho"],
            "gse230424_he_rho": m["g230"]["top_he_sample_centered_rho"],
            "gse230424_residual_target_rho": m["g230_resid"]["dm1_coord_qc_residual_target_he_rho"],
        },
        "figures": [fig.name for fig in figs],
        "voice_boundary": "scaffold/infra/table/figure only",
    }
    (OUT / "PAPER2_EXTERNAL_VALIDATION_READY_PACK.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    md = [
        "# Paper 2 external-validation-ready pack",
        "",
        "Top purpose: make K2/Bundang H&E validation executable immediately after data receipt.",
        "",
        "Figures:",
    ]
    for fig in figs:
        md.append(f"- `{fig.name}`")
    md += ["", "Tables:", *[f"- `paper2_external_validation_{name}.tsv`" for name in tables], "", "Voice boundary: no protected manuscript prose generated."]
    (OUT / "SUMMARY.md").write_text("\n".join(md) + "\n", encoding="utf-8")


def publish(figs: list[Path]) -> None:
    ASSET_LOCAL.mkdir(parents=True, exist_ok=True)
    ASSET_LIVE.mkdir(parents=True, exist_ok=True)
    for fig in figs:
        shutil.copy2(fig, ASSET_LOCAL / fig.name)
        shutil.copy2(fig, ASSET_LIVE / fig.name)
    for extra in OUT.glob("paper2_external_validation_*.tsv"):
        shutil.copy2(extra, ASSET_LOCAL / extra.name)
        shutil.copy2(extra, ASSET_LIVE / extra.name)
    for extra in [OUT / "README_FOR_COLLABORATOR_HANDOFF.md", OUT / "PAPER2_EXTERNAL_VALIDATION_READY_PACK.json", OUT / "SUMMARY.md"]:
        shutil.copy2(extra, ASSET_LOCAL / extra.name)
        shutil.copy2(extra, ASSET_LIVE / extra.name)


def html_table(df: pd.DataFrame) -> str:
    out = [
        "<table><thead><tr>"
        + "".join(f"<th>{html.escape(str(c))}</th>" for c in df.columns)
        + "</tr></thead><tbody>"
    ]
    for _, row in df.iterrows():
        out.append("<tr>" + "".join(f"<td>{html.escape(str(row[c]))}</td>" for c in df.columns) + "</tr>")
    out.append("</tbody></table>")
    return "\n".join(out)


def asset_link(filename: str, label: str) -> str:
    return f'<a href="assets/{ASSET}/{html.escape(filename)}" download>{html.escape(label)}</a>'


def write_page(tables: dict[str, pd.DataFrame], figs: list[Path], m: dict) -> None:
    fig_buttons = "\n".join(
        f'<button class="fig-option" data-src="assets/{ASSET}/{fig.name}" data-title="{html.escape(fig.stem)}">'
        f'<img src="assets/{ASSET}/{fig.name}" alt=""><span>{html.escape(fig.stem.replace("_", " "))}</span></button>'
        for fig in figs
    )
    download_links = " · ".join(
        [
            asset_link("paper2_external_validation_manifest.tsv", "manifest template"),
            asset_link("paper2_external_validation_reader.tsv", "reader CRF"),
            asset_link("paper2_external_validation_gates.tsv", "go/no-go gates"),
            asset_link("README_FOR_COLLABORATOR_HANDOFF.md", "handoff README"),
        ]
    )
    html_doc = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Paper 2 External Validation Ready Pack</title>
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
.stats {{ display:grid; grid-template-columns:repeat(5,minmax(0,1fr)); gap:10px; max-width:1200px; margin-top:24px; }}
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
.downloads {{ border:1px solid var(--line); background:var(--panel2); border-radius:8px; padding:14px; margin:16px 0 0; }}
table {{ width:100%; border-collapse:collapse; margin:16px 0; font-size:12px; }}
th,td {{ border-bottom:1px solid var(--line); padding:10px 8px; text-align:left; vertical-align:top; }}
th {{ color:var(--gold); background:rgba(214,178,94,.06); }}
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
    <div class="kicker">Paper 2 · external validation ready pack · 2026-05-10</div>
    <h1>Make K2/Bundang Validation Executable, Not Aspirational</h1>
    <p class="lead">A locked handoff packet: data request, sample manifest schema, statistical analysis plan, go/no-go thresholds, pathologist hotspot audit sheet, and five-day execution board.</p>
    <div class="stats">
      <div class="stat"><div class="v">WSI</div><div class="l">de-identified H&E required</div></div>
      <div class="stat"><div class="v">{m['uni']['pooled_overall_auc']:.3f}</div><div class="l">current UNI LOTO AUC</div></div>
      <div class="stat"><div class="v">{m['g250_stage']['dm1_existing_slide_centered_rho']:.3f}</div><div class="l">GSE250521 rho</div></div>
      <div class="stat"><div class="v">{m['g230']['top_he_sample_centered_rho']:.3f}</div><div class="l">GSE230424 rho</div></div>
      <div class="stat"><div class="v">NC</div><div class="l">green result ceiling</div></div>
    </div>
  </div>
</header>
<div class="layout">
<nav>
  <a href="#verdict">01 Verdict</a>
  <a href="#figures">02 Figures</a>
  <a href="#request">03 Data Request</a>
  <a href="#manifest">04 Manifest</a>
  <a href="#sap">05 Locked SAP</a>
  <a href="#gates">06 Go/No-Go</a>
  <a href="#reader">07 Reader Audit</a>
  <a href="#deliverables">08 Execution</a>
  <a href="#routes">09 Journal Routes</a>
  <a href="#sources">10 Sources</a>
</nav>
<main>
<section id="verdict">
  <h2><span class="num">01</span>Ready Packet Verdict</h2>
  <div class="grid">
    <div class="card"><h3>What changed</h3><p class="muted">The K2/Bundang ask is now operational: exact files, manifest fields, endpoints, QC rules, and decision thresholds.</p></div>
    <div class="card"><h3>Why it raises impact</h3><p class="muted">Positive locked external H&E validation is the one result that moves Paper 2 from CRM-ready to Nat Commun-plausible.</p></div>
    <div class="card"><h3>Guardrail</h3><p class="muted">No retuning after data receipt. Failed strata and QC exclusions remain visible.</p></div>
  </div>
  <div class="downloads"><strong>Downloads:</strong> {download_links}</div>
</section>
<section id="figures">
  <h2><span class="num">02</span>External Validation Figure Viewer</h2>
  <div class="viewer">
    <div class="fig-list" id="figList">{fig_buttons}</div>
    <div class="stage">
      <div class="toolbar"><strong id="figTitle">Figure</strong><button class="control" id="zoomOut">Zoom out</button><button class="control" id="zoomReset">Reset</button><button class="control" id="zoomIn">Zoom in</button></div>
      <div class="canvas" id="canvas"><img id="mainFig" src="assets/{ASSET}/{figs[0].name}" alt="Paper 2 external validation figure"></div>
    </div>
  </div>
</section>
<section id="request"><h2><span class="num">03</span>Data Request</h2>{html_table(tables['request'])}</section>
<section id="manifest"><h2><span class="num">04</span>Sample Manifest Template</h2>{html_table(tables['manifest'])}</section>
<section id="sap"><h2><span class="num">05</span>Locked Statistical Analysis Plan</h2>{html_table(tables['sap'])}</section>
<section id="gates"><h2><span class="num">06</span>Go/No-Go Thresholds</h2>{html_table(tables['gates'])}</section>
<section id="reader"><h2><span class="num">07</span>Pathologist Hotspot Audit CRF</h2>{html_table(tables['reader'])}</section>
<section id="deliverables"><h2><span class="num">08</span>Five-Day Execution Board</h2>{html_table(tables['deliverables'])}</section>
<section id="routes"><h2><span class="num">09</span>Journal Route Triggers</h2>{html_table(tables['routes'])}</section>
<section id="sources">
  <h2><span class="num">10</span>Local Sources</h2>
  <table><tbody>
    <tr><td>Output directory</td><td><code>{OUT.relative_to(ROOT)}</code></td></tr>
    <tr><td>Builder script</td><td><code>project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/scripts/build_paper2_external_validation_ready_pack.py</code></td></tr>
    <tr><td>Max-impact pack</td><td><a href="paper2_max_impact_unlock_pack.html">paper2_max_impact_unlock_pack.html</a></td></tr>
    <tr><td>Reference leverage pack</td><td><a href="paper2_reference_leverage_pack.html">paper2_reference_leverage_pack.html</a></td></tr>
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
    tables = make_tables(m)
    figs = [
        draw_request_packet(tables),
        draw_manifest_schema(tables),
        draw_locked_workflow(tables),
        draw_go_nogo(tables),
        draw_reader_audit(tables),
        draw_timeline(tables),
    ]
    write_outputs(tables, figs, m)
    publish(figs)
    write_page(tables, figs, m)
    print(json.dumps({"out": str(OUT), "page": str(PAGE_LOCAL), "figures": [p.name for p in figs]}, indent=2))


if __name__ == "__main__":
    main()
