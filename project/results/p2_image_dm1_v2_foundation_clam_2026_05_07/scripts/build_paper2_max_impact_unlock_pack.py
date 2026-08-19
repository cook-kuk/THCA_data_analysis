#!/usr/bin/env python3
"""Build a Paper 2 max-impact unlock dossier.

Scaffold/figure/table work only. This does not generate protected manuscript
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
OUT = SUPP / "paper2_max_impact_unlock_2026_05_10"
HUB = ROOT / "project/papers_hub_2026_05_04"
LIVE = Path("/var/www/papers/papers_hub_2026_05_04")
ASSET = "paper2_max_impact_unlock"
ASSET_LOCAL = HUB / "assets" / ASSET
ASSET_LIVE = LIVE / "assets" / ASSET
PAGE_LOCAL = HUB / "paper2_max_impact_unlock_pack.html"
PAGE_LIVE = LIVE / "paper2_max_impact_unlock_pack.html"

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
    put(img, heading, (70, 142), 1.08, INK, 2, 1650)
    put(img, sub, (70, 205), 0.56, MUTED, 1, 1550)


def make_tables(m: dict) -> dict[str, pd.DataFrame]:
    unlocks = pd.DataFrame(
        [
            [
                "U1",
                "K2/Bundang locked external H&E validation",
                "highest",
                "+25",
                "A real independent thyroid H&E cohort turns the story from strong internal/external spatial support into external clinical-pathology validation.",
                "Frozen model; no retuning; report AUC/rho/calibration; keep fail state visible.",
                "Nature Communications ceiling; Lancet DH presub credibility.",
            ],
            [
                "U2",
                "HEST-1k / HEST-Benchmark bridge",
                "high",
                "+14",
                "Shows the manuscript is aligned with a 2024-2026 field benchmark for histology-spatial transcriptomics and foundation-model gene-expression prediction.",
                "Do not download all >1TB now; query subset or cite benchmark; run only if compute/data budget is approved.",
                "Converts 'small thyroid AI' into field-standard benchmark-aware digital pathology.",
            ],
            [
                "U3",
                "Pathologist-reader hotspot audit",
                "high",
                "+12",
                "A blinded pathology plausibility audit makes the heatmaps clinically legible and reduces black-box concern.",
                "Reader marks morphology categories; compare with top-decile image-DM1 hotspots; no outcome claim.",
                "Improves digital-health and translational medicine fit.",
            ],
            [
                "U4",
                "Calibration + reflex-test threshold",
                "high",
                "+11",
                "Editors need a concrete clinical action: H&E-positive cases go to molecular/spatial confirmation; H&E-negative cases deprioritized.",
                "Decision-curve/calibration scaffold only unless external cohort exists.",
                "Needed for Lancet Digital Health; useful for Cell Reports Medicine.",
            ],
            [
                "U5",
                "Caveat moat as a feature",
                "medium",
                "+8",
                f"Smoothness/QC caveats are already measured: residual target rho {m['g230_resid']['dm1_coord_qc_residual_target_he_rho']:.3f}; AITD control rho {m['aitd']['ap_tls_he_sample_centered_rho']:.3f}; random-module residual p {m['g250_spec']['actual_residual_empirical_p']:.4f}.",
                "Show raw and residual results side by side; make negative control visible.",
                "Raises reviewer trust and prevents overclaim desk rejection.",
            ],
            [
                "U6",
                "Paper 1 public linkage after Paper 1 leaves",
                "medium",
                "+7",
                "Once Paper 1 is out/preprint, Paper 2 can cite the DM1/RAI biology without looking circular.",
                "Use as factual cross-reference only; do not rewrite protected discussion.",
                "Stronger biology frame for Nature Cancer/Nat Commun presub.",
            ],
            [
                "U7",
                "One-figure editor storyboard",
                "medium",
                "+6",
                "A single graphical abstract should show H&E input, foundation features, spatial RNA recovery, thyroid proof, and external unlock.",
                "Use existing separated figure viewer; avoid composite-only evidence dump.",
                "Improves editor comprehension in 30 seconds.",
            ],
        ],
        columns=["unlock_id", "unlock", "priority", "impact_delta", "why_it_raises_impact", "execution_rule", "journal_effect"],
    )
    validation = pd.DataFrame(
        [
            ["1", "Cohort lock", "K2/Bundang FFPE H&E", "Case IDs, slide IDs, scanner/stain/site metadata", "No model changes after this point"],
            ["2", "Model lock", "Existing UNI/CLAM or UNI feature pipeline", "Git hash, model weights, threshold", "No retuning on external cohort"],
            ["3", "Primary endpoint", "DM1/RAI molecular or proxy label", "AUC if binary; Spearman rho if continuous", "Report exact missingness"],
            ["4", "Secondary endpoint", "Hotspot/localization", "Top-decile enrichment and pathologist-reader overlap", "Exploratory unless pre-specified"],
            ["5", "Controls", "QC/stain/slide/site", "Residual, block, downsample, random module if ST available", "Failures shown in supplement"],
            ["6", "Decision", "Route choice", "NC/LDH presub if positive; CRM if absent/negative", "No indefinite waiting"],
        ],
        columns=["step", "module", "data", "readout", "guardrail"],
    )
    hest = pd.DataFrame(
        [
            [
                "What HEST adds",
                "A field benchmark for histology + spatial transcriptomics, including HEST-Benchmark for foundation-model gene-expression prediction.",
                "Use to show Paper 2 is benchmark-aware and not isolated from the 2024-2026 H&E-to-ST field.",
                "It is not automatically thyroid validation.",
            ],
            [
                "Current public scale",
                "HEST GitHub reports HEST-1k access with 1,276 paired spatial transcriptomics samples plus H&E-stained WSI, and HEST-Benchmark results for 25 models as of 2026-04-03.",
                "Use as a field-scale contrast to our focused thyroid validation.",
                "Full dataset is >1TB; do not casually download during manuscript sprint.",
            ],
            [
                "Fastest safe action",
                "Add a no-download benchmark-aware table and, if needed, query only metadata/subsets through HEST-Library.",
                "Improves Methods/Discussion credibility without blocking Paper 2.",
                "Do not claim HEST proves DM1/RAI.",
            ],
            [
                "Best optional experiment",
                "Run a small HEST-Benchmark-style ridge/PCA readout using our preferred embedding family on a permitted subset.",
                "Would support morphology-to-expression generality.",
                "Only if compute and data terms are acceptable.",
            ],
        ],
        columns=["item", "fact", "paper2_use", "boundary"],
    )
    risk = pd.DataFrame(
        [
            ["No independent thyroid H&E cohort", "K2/Bundang locked validation", "red -> green if positive", "Fig 6 / Table 1"],
            ["Smoothness/QC confounding", "Raw/residual/random-module/block controls already generated", "amber -> managed", "Fig 5 / caveat box"],
            ["Small disease-specific story", "Cell 2026 + HEST + HEX precedent bridge", "amber -> field-aligned", "Fig 1 / Table 2"],
            ["Black-box hotspot maps", "Blinded pathologist-reader hotspot audit", "amber -> clinical legibility", "Supplementary figure"],
            ["Clinical utility not proven", "Reflex-test threshold and calibration scaffold", "red -> presub-ready only", "Fig 6 / presub"],
            ["Foundation model novelty concern", "Position as application/validation, not new FM", "green", "Methods / references"],
        ],
        columns=["editor_risk", "impact_fix", "status_after_fix", "where_to_show"],
    )
    journal = pd.DataFrame(
        [
            ["Cell Reports Medicine", "Submit-ready if caveats remain visible", "Current evidence + reference leverage + visual dossier", "Do not over-wait for perfect K2"],
            ["Nature Communications", "Plausible after external H&E validation", "K2/Bundang positive locked result", "If no external cohort, high-risk submit/presub only"],
            ["The Lancet Digital Health", "Presub moonshot", "External H&E + calibration/reflex-test scenario", "No direct full submit without clinical utility"],
            ["Nature Medicine", "No-go now", "Would need prospective/outcome/clinical decision impact", "Keep as future route"],
            ["Nature Cancer", "Presub only if biology strengthened", "Paper 1 public + DM1/RAI mechanism linkage", "No causal mechanism from H&E alone"],
        ],
        columns=["journal", "impact_route", "minimum_unlock", "decision_rule"],
    )
    presub = pd.DataFrame(
        [
            ["Attachment 1", "One-page graphical abstract", "F01/F07 from max-impact packs", "30-second editor comprehension"],
            ["Attachment 2", "Reference-to-claim matrix", "paper2_reference_leverage_pack.html", "Shows field alignment"],
            ["Attachment 3", "Evidence moat", "paper2_impact_upgrade.html", "AUC/rho/hotspots/caveats together"],
            ["Attachment 4", "Separated figure browser", "paper2_cv2_visual_summary.html", "Editor can inspect individual panels"],
            ["Attachment 5", "External validation plan/result", "K2/Bundang locked protocol", "Determines NC/LDH ceiling"],
        ],
        columns=["slot", "attachment", "source", "why_editor_cares"],
    )
    return {"unlocks": unlocks, "validation": validation, "hest": hest, "risk": risk, "journal": journal, "presub": presub}


def draw_control_room(tables: dict[str, pd.DataFrame], m: dict) -> Path:
    img = np.full((1700, 2300, 3), BG, dtype=np.uint8)
    title(
        img,
        "Max-impact control room",
        "The Ceiling Moves Only When A Risk Disappears",
        "Current Paper 2 is strong. The next move is external validation + benchmark-aware field alignment.",
    )
    stats = [
        ("Current core", "CRM-ready", GREEN),
        ("UNI LOTO", f"AUC {m['uni']['pooled_overall_auc']:.3f}", GOLD),
        ("GSE250521", f"rho {m['g250_stage']['dm1_existing_slide_centered_rho']:.3f}", TEAL),
        ("GSE230424", f"rho {m['g230']['top_he_sample_centered_rho']:.3f}", BLUE),
        ("Residual target", f"rho {m['g230_resid']['dm1_coord_qc_residual_target_he_rho']:.3f}", PURPLE),
        ("Hotspot", f"{m['hotspot']['gse230424_top10_lift']:.1f}x", GOLD),
    ]
    for i, (head, val, color) in enumerate(stats):
        x = 80 + i * 365
        box(img, (x, 315, x + 325, 500), PANEL2, color)
        put(img, head.upper(), (x + 18, 360), 0.4, MUTED, 1)
        put(img, val, (x + 18, 430), 0.78, color, 2, 280)
    lanes = [
        ("NOW", "Cell Reports Medicine", "Strong cross-modal biomarker recovery with caveat moat", GREEN),
        ("+K2/Bundang", "Nature Communications", "Independent thyroid H&E validation removes the biggest editor objection", GOLD),
        ("+clinical threshold", "Lancet DH presub", "Reflex-test workflow makes the story digital-health relevant", BLUE),
        ("+Paper 1 public", "Nature Cancer presub", "Biology anchor becomes citable without circularity", PURPLE),
    ]
    for i, (tag, route, body, color) in enumerate(lanes):
        yy = 650 + i * 235
        box(img, (120, yy, 2180, yy + 165), PANEL if i % 2 == 0 else PANEL2, color)
        put(img, tag, (160, yy + 62), 0.68, color, 2, 300)
        put(img, route, (520, yy + 58), 0.72, INK, 2, 460)
        put(img, body, (1030, yy + 55), 0.54, MUTED, 1, 940)
    box(img, (120, 1620, 2180, 1720), (42, 34, 22), GOLD)
    put(img, "Priority order: K2/Bundang external H&E > HEST benchmark bridge > pathologist hotspot audit > calibration/reflex threshold.", (155, 1680), 0.58, GOLD, 2, 1900)
    path = OUT / "F01_max_impact_control_room.png"
    cv2.imwrite(str(path), img)
    return path


def draw_unlock_ladder(tables: dict[str, pd.DataFrame]) -> Path:
    img = np.full((1800, 2500, 3), BG, dtype=np.uint8)
    title(img, "Impact unlock ladder", "Seven Moves Ranked By Journal-Ceiling Gain", "Do U1 first if any external H&E path is real.")
    colors = [GOLD, TEAL, BLUE, BLUE, RED, PURPLE, GREEN]
    for i, row in tables["unlocks"].iterrows():
        yy = 320 + i * 195
        color = colors[i]
        box(img, (80, yy, 2380, yy + 145), PANEL if i % 2 == 0 else PANEL2, color)
        put(img, row.unlock_id, (115, yy + 60), 0.72, color, 2)
        put(img, row.unlock, (235, yy + 52), 0.58, INK, 2, 520)
        put(img, row.impact_delta, (820, yy + 64), 0.82, color, 2)
        put(img, row.why_it_raises_impact, (990, yy + 42), 0.43, MUTED, 1, 670)
        put(img, row.journal_effect, (1720, yy + 42), 0.43, color, 1, 560)
    path = OUT / "F02_impact_unlock_ladder.png"
    cv2.imwrite(str(path), img)
    return path


def draw_external_validation_protocol(tables: dict[str, pd.DataFrame]) -> Path:
    img = np.full((1600, 2300, 3), BG, dtype=np.uint8)
    title(img, "Locked external validation", "K2/Bundang Is The Single Highest-Value Experiment", "This is the protocol board for turning Paper 2 into a Nat Commun-level story.")
    y0 = 330
    for i, row in tables["validation"].iterrows():
        x = 120 + i * 350
        box(img, (x, y0, x + 300, y0 + 520), PANEL if i % 2 == 0 else PANEL2, GOLD if i < 3 else BLUE)
        cv2.circle(img, (x + 44, y0 + 48), 28, GOLD if i < 3 else BLUE, -1, cv2.LINE_AA)
        put(img, row.step, (x + 34, y0 + 59), 0.6, BG, 2)
        put(img, row.module, (x + 25, y0 + 125), 0.54, INK, 2, 245)
        put(img, row.data, (x + 25, y0 + 225), 0.43, MUTED, 1, 245)
        put(img, row.readout, (x + 25, y0 + 340), 0.43, GOLD, 1, 245)
        put(img, row.guardrail, (x + 25, y0 + 445), 0.39, RED if i < 2 else MUTED, 1, 245)
        if i < len(tables["validation"]) - 1:
            cv2.arrowedLine(img, (x + 300, y0 + 260), (x + 350, y0 + 260), GOLD, 3, cv2.LINE_AA, tipLength=0.16)
    box(img, (120, 1060, 2180, 1210), PANEL2, RED)
    put(img, "No retuning. No hidden thresholds. No missing-failure supplement. This is what converts a strong model paper into an editor-trustworthy validation paper.", (155, 1130), 0.58, INK, 2, 1900)
    path = OUT / "F03_locked_external_validation_protocol.png"
    cv2.imwrite(str(path), img)
    return path


def draw_hest_bridge(tables: dict[str, pd.DataFrame]) -> Path:
    img = np.full((1700, 2300, 3), BG, dtype=np.uint8)
    title(
        img,
        "HEST benchmark bridge",
        "Field-Scale Benchmark Awareness Without Derailing The Manuscript",
        "HEST-1k/HEST-Benchmark can raise credibility, but it should not replace thyroid validation.",
    )
    left = (90, 335, 1050, 945)
    right = (1250, 335, 2210, 945)
    box(img, left, PANEL2, TEAL)
    box(img, right, PANEL2, GOLD)
    put(img, "HEST-1k / HEST-Benchmark", (135, 420), 0.78, TEAL, 2)
    put(img, "1,276 paired ST + H&E WSI samples reported by HEST GitHub; HEST-Benchmark tests foundation models for gene-expression prediction from morphology.", (135, 505), 0.54, INK, 1, 830)
    put(img, "Use: benchmark-aware framing, optional subset experiment, field-standard comparison.", (135, 710), 0.54, TEAL, 2, 810)
    put(img, "Paper 2", (1295, 420), 0.78, GOLD, 2)
    put(img, "Focused thyroid DM1/RAI validation: TCGA H&E classifier + GSE250521/GSE230424 spatial RNA support + external H&E unlock.", (1295, 505), 0.54, INK, 1, 820)
    put(img, "Use: disease-specific translational biomarker recovery, not a generic benchmark paper.", (1295, 710), 0.54, GOLD, 2, 800)
    cv2.arrowedLine(img, (1050, 640), (1250, 640), TEAL, 4, cv2.LINE_AA, tipLength=0.12)
    rows = tables["hest"]
    for i, row in rows.iterrows():
        yy = 1070 + i * 145
        box(img, (120, yy, 2180, yy + 105), PANEL if i % 2 == 0 else PANEL2)
        put(img, row.item, (155, yy + 60), 0.56, TEAL, 2, 360)
        put(img, row.fact, (540, yy + 40), 0.42, INK, 1, 640)
        put(img, row.paper2_use, (1220, yy + 40), 0.42, GOLD, 1, 430)
        put(img, row.boundary, (1695, yy + 40), 0.42, RED, 1, 430)
    path = OUT / "F04_hest_benchmark_bridge.png"
    cv2.imwrite(str(path), img)
    return path


def draw_risk_annihilation(tables: dict[str, pd.DataFrame]) -> Path:
    img = np.full((1600, 2350, 3), BG, dtype=np.uint8)
    title(img, "Editor risk buy-down", "Turn Weaknesses Into Managed Reviewer Questions", "This is the table that should guide supplement order and presub attachments.")
    y0 = 330
    for i, row in tables["risk"].iterrows():
        yy = y0 + i * 190
        color = RED if "red" in row.status_after_fix else GOLD if "amber" in row.status_after_fix else GREEN
        box(img, (90, yy, 2260, yy + 135), PANEL if i % 2 == 0 else PANEL2, color)
        put(img, row.editor_risk, (125, yy + 48), 0.5, RED, 2, 430)
        put(img, row.impact_fix, (610, yy + 48), 0.5, INK, 1, 620)
        put(img, row.status_after_fix, (1285, yy + 55), 0.52, color, 2, 320)
        put(img, row.where_to_show, (1660, yy + 55), 0.5, GOLD, 1, 470)
    path = OUT / "F05_editor_risk_buydown.png"
    cv2.imwrite(str(path), img)
    return path


def draw_submission_packet(tables: dict[str, pd.DataFrame]) -> Path:
    img = np.full((1600, 2350, 3), BG, dtype=np.uint8)
    title(img, "Submission packet upgrade", "What To Put In Front Of The Editor", "A short, visual, traceable packet beats a long claim-heavy cover letter.")
    for i, row in tables["presub"].iterrows():
        yy = 340 + i * 220
        color = [GOLD, TEAL, BLUE, GREEN, RED][i]
        box(img, (130, yy, 2170, yy + 155), PANEL if i % 2 == 0 else PANEL2, color)
        cv2.circle(img, (190, yy + 78), 36, color, -1, cv2.LINE_AA)
        put(img, str(i + 1), (180, yy + 91), 0.75, BG, 2)
        put(img, row.attachment, (270, yy + 58), 0.62, color, 2, 470)
        put(img, row.source, (795, yy + 55), 0.45, MUTED, 1, 620)
        put(img, row.why_editor_cares, (1485, yy + 62), 0.5, INK, 1, 560)
    box(img, (130, 1535, 2170, 1640), (42, 34, 22), GOLD)
    put(img, "Presub rule: one sentence claim, one figure, one validation plan, one caveat box. Do not bury the editor in every exploratory result.", (165, 1598), 0.56, GOLD, 2, 1880)
    path = OUT / "F06_editor_presub_packet.png"
    cv2.imwrite(str(path), img)
    return path


def write_outputs(tables: dict[str, pd.DataFrame], figs: list[Path], m: dict) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for name, df in tables.items():
        df.to_csv(OUT / f"paper2_max_impact_{name}.tsv", sep="\t", index=False)
    summary = {
        "created": "2026-05-10",
        "purpose": "Paper 2 max-impact unlock plan",
        "top_unlock": "K2/Bundang locked external H&E validation",
        "second_unlock": "HEST-1k / HEST-Benchmark bridge",
        "current_route": "Cell Reports Medicine-ready",
        "ceiling_route": "Nature Communications plausible if external H&E validation is positive",
        "moonshot_route": "The Lancet Digital Health presub only after external H&E plus calibration/reflex-test framing",
        "headline_metrics": {
            "uni_loto_auc": m["uni"]["pooled_overall_auc"],
            "gse250521_slide_centered_rho": m["g250_stage"]["dm1_existing_slide_centered_rho"],
            "gse230424_he_rho": m["g230"]["top_he_sample_centered_rho"],
            "gse230424_residual_target_rho": m["g230_resid"]["dm1_coord_qc_residual_target_he_rho"],
            "gse250521_random_module_residual_p": m["g250_spec"]["actual_residual_empirical_p"],
            "aitd_negative_control_rho": m["aitd"]["ap_tls_he_sample_centered_rho"],
        },
        "figures": [fig.name for fig in figs],
        "voice_boundary": "Scaffold/figure/table work only.",
    }
    (OUT / "PAPER2_MAX_IMPACT_UNLOCK_PACK.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    md = [
        "# Paper 2 max-impact unlock pack",
        "",
        "Purpose: define the smallest set of actions that materially raises the Paper 2 journal ceiling.",
        "",
        "Top unlock: K2/Bundang locked external H&E validation.",
        "Second unlock: HEST-1k / HEST-Benchmark bridge without derailing the manuscript.",
        "",
        "Generated figures:",
    ]
    for fig in figs:
        md.append(f"- `{fig.name}`")
    md += ["", "Voice boundary: no protected manuscript prose generated."]
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
        out.append("<tr>" + "".join(f"<td>{html.escape(str(row[c]))}</td>" for c in df.columns) + "</tr>")
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
<title>Paper 2 Max-Impact Unlock Pack</title>
<style>
:root {{ --bg:#0d1117; --panel:#151b23; --panel2:#1c2634; --ink:#e6edf3; --muted:#98a6ba; --line:#303846; --gold:#d6b25e; --blue:#57a6ff; --green:#9dd335; --red:#ff7b72; }}
* {{ box-sizing:border-box; }}
body {{ margin:0; background:var(--bg); color:var(--ink); font-family:"JetBrains Mono", ui-monospace, Menlo, monospace; line-height:1.55; }}
a {{ color:var(--blue); text-decoration:none; }}
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
    <div class="kicker">Paper 2 · max-impact unlock pack · 2026-05-10</div>
    <h1>Move The Journal Ceiling By Removing The Editor's Objections</h1>
    <p class="lead">A concrete action dossier: K2/Bundang locked external validation, HEST benchmark bridge, pathologist hotspot audit, calibration/reflex-test threshold, and caveat moat.</p>
    <div class="stats">
      <div class="stat"><div class="v">U1</div><div class="l">K2/Bundang top unlock</div></div>
      <div class="stat"><div class="v">HEST</div><div class="l">benchmark bridge</div></div>
      <div class="stat"><div class="v">{m['uni']['pooled_overall_auc']:.3f}</div><div class="l">UNI LOTO AUC</div></div>
      <div class="stat"><div class="v">{m['g250_stage']['dm1_existing_slide_centered_rho']:.3f}</div><div class="l">GSE250521 rho</div></div>
      <div class="stat"><div class="v">{m['g230']['top_he_sample_centered_rho']:.3f}</div><div class="l">GSE230424 rho</div></div>
      <div class="stat"><div class="v">NC</div><div class="l">ceiling with external H&E</div></div>
    </div>
  </div>
</header>
<div class="layout">
<nav>
  <a href="#verdict">01 Verdict</a>
  <a href="#figures">02 Figures</a>
  <a href="#unlocks">03 Unlocks</a>
  <a href="#validation">04 External Validation</a>
  <a href="#hest">05 HEST Bridge</a>
  <a href="#risk">06 Risk Buy-Down</a>
  <a href="#journal">07 Journal Thresholds</a>
  <a href="#packet">08 Presub Packet</a>
  <a href="#sources">09 Sources</a>
</nav>
<main>
<section id="verdict">
  <h2><span class="num">01</span>Max-Impact Verdict</h2>
  <div class="grid">
    <div class="card"><h3>One highest-value move</h3><p class="muted">K2/Bundang locked external H&E validation. This is the only move that materially changes the Nature Communications ceiling.</p></div>
    <div class="card"><h3>Second move</h3><p class="muted">HEST-1k/HEST-Benchmark bridge. Use it to show field-standard benchmark awareness without downloading terabytes during the manuscript sprint.</p></div>
    <div class="card"><h3>Guardrail</h3><p class="muted">Do not claim clinical deployment, pan-disease validation, or direct spatial-omics measurement.</p></div>
  </div>
</section>
<section id="figures">
  <h2><span class="num">02</span>Max-Impact Figure Viewer</h2>
  <div class="viewer">
    <div class="fig-list" id="figList">{fig_buttons}</div>
    <div class="stage">
      <div class="toolbar"><strong id="figTitle">Figure</strong><button class="control" id="zoomOut">Zoom out</button><button class="control" id="zoomReset">Reset</button><button class="control" id="zoomIn">Zoom in</button></div>
      <div class="canvas" id="canvas"><img id="mainFig" src="assets/{ASSET}/{figs[0].name}" alt="Paper 2 max-impact figure"></div>
    </div>
  </div>
</section>
<section id="unlocks"><h2><span class="num">03</span>Impact Unlock Ladder</h2>{html_table(tables['unlocks'])}</section>
<section id="validation"><h2><span class="num">04</span>K2/Bundang Locked External Validation</h2>{html_table(tables['validation'])}</section>
<section id="hest"><h2><span class="num">05</span>HEST-1k / HEST-Benchmark Bridge</h2>{html_table(tables['hest'])}</section>
<section id="risk"><h2><span class="num">06</span>Editor Risk Buy-Down</h2>{html_table(tables['risk'])}</section>
<section id="journal"><h2><span class="num">07</span>Journal Thresholds</h2>{html_table(tables['journal'])}</section>
<section id="packet"><h2><span class="num">08</span>Presubmission Packet</h2>{html_table(tables['presub'])}</section>
<section id="sources">
  <h2><span class="num">09</span>Local Sources</h2>
  <table><tbody>
    <tr><td>Output directory</td><td><code>{OUT.relative_to(ROOT)}</code></td></tr>
    <tr><td>Builder script</td><td><code>project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/scripts/build_paper2_max_impact_unlock_pack.py</code></td></tr>
    <tr><td>Reference leverage pack</td><td><a href="paper2_reference_leverage_pack.html">paper2_reference_leverage_pack.html</a></td></tr>
    <tr><td>Moonshot pack</td><td><a href="paper2_moonshot_editorial_pack.html">paper2_moonshot_editorial_pack.html</a></td></tr>
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
        draw_control_room(tables, m),
        draw_unlock_ladder(tables),
        draw_external_validation_protocol(tables),
        draw_hest_bridge(tables),
        draw_risk_annihilation(tables),
        draw_submission_packet(tables),
    ]
    write_outputs(tables, figs, m)
    publish(figs)
    write_page(tables, figs, m)
    print(json.dumps({"out": str(OUT), "page": str(PAGE_LOCAL), "figures": [p.name for p in figs]}, indent=2))


if __name__ == "__main__":
    main()
