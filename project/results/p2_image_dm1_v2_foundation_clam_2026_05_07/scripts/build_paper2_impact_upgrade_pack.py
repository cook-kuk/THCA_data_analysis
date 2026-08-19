#!/usr/bin/env python3
"""Build a Paper 2 impact-upgrade pack and web dossier assets."""
from __future__ import annotations

import json
import math
import shutil
from pathlib import Path

import cv2
import numpy as np
import pandas as pd


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
P2 = ROOT / "project/results/p2_image_dm1_v2_foundation_clam_2026_05_07"
SUPP = P2 / "analysis_supp"
OUT = SUPP / "paper2_impact_upgrade_2026_05_10"
HUB = ROOT / "project/papers_hub_2026_05_04"
LIVE = Path("/var/www/papers/papers_hub_2026_05_04")
ASSET_NAME = "paper2_impact_upgrade"
ASSET_LOCAL = HUB / "assets" / ASSET_NAME
ASSET_LIVE = LIVE / "assets" / ASSET_NAME
PAGE_LOCAL = HUB / "paper2_impact_upgrade.html"
PAGE_LIVE = LIVE / "paper2_impact_upgrade.html"


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
    uni = read_json(SUPP / "audit_uni_loto/UNI_LOTO_SUMMARY.json")
    g250 = read_json(SUPP / "path2space_inspired_reanalysis_2026_05_09/PATH2SPACE_REANALYSIS_SUMMARY.json")
    g250_stage = read_json(
        SUPP / "path2space_stage_generalization_controls_2026_05_09/GSE250521_STAGE_GENERALIZATION_SUMMARY.json"
    )
    g250_spec = read_json(
        SUPP / "path2space_gse250521_random_module_specificity_2026_05_09/GSE250521_RANDOM_MODULE_SPECIFICITY_SUMMARY.json"
    )
    g230 = read_json(SUPP / "gse230424_pathology_thyroid_axis_2026_05_09/GSE230424_PATHOLOGY_THYROID_SUMMARY.json")
    g230_resid = read_json(
        SUPP
        / "gse230424_pathology_thyroid_axis_2026_05_09/residual_target_controls/GSE230424_RESIDUAL_TARGET_SUMMARY.json"
    )
    decile = read_json(SUPP / "path2space_decile_dose_response_2026_05_09/PATH2SPACE_DECILE_DOSE_RESPONSE_SUMMARY.json")
    down = read_json(SUPP / "path2space_downsample_robustness_2026_05_09/PATH2SPACE_DOWNSAMPLE_ROBUSTNESS_SUMMARY.json")
    hotspot = read_json(SUPP / "path2space_hotspot_concordance_2026_05_09/PATH2SPACE_HOTSPOT_CONCORDANCE_SUMMARY.json")
    aitd = read_json(SUPP / "gse248205_pathology_aitd_axis_2026_05_09/GSE248205_PATHOLOGY_AITD_SUMMARY.json")
    return {
        "uni": uni,
        "g250": g250,
        "g250_stage": g250_stage,
        "g250_spec": g250_spec,
        "g230": g230,
        "g230_resid": g230_resid,
        "decile": decile,
        "down": down,
        "hotspot": hotspot,
        "aitd": aitd,
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
    for raw in text.split("\n"):
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
    h = cv2.getTextSize("Ag", font, scale, thick)[0][1] + gap
    for i, line in enumerate(lines):
        cv2.putText(img, line, (x, y + i * h), font, scale, color, thick, cv2.LINE_AA)
    return y + max(len(lines), 1) * h


def box(img: np.ndarray, xyxy: tuple[int, int, int, int], color: tuple[int, int, int], outline: tuple[int, int, int] = LINE) -> None:
    x0, y0, x1, y1 = xyxy
    cv2.rectangle(img, (x0, y0), (x1, y1), color, -1)
    cv2.rectangle(img, (x0, y0), (x1, y1), outline, 1, cv2.LINE_AA)


def bar(img: np.ndarray, x: int, y: int, w: int, h: int, frac: float, color: tuple[int, int, int]) -> None:
    cv2.rectangle(img, (x, y), (x + w, y + h), (43, 50, 62), -1)
    cv2.rectangle(img, (x, y), (x + int(w * np.clip(frac, 0, 1)), y + h), color, -1)
    cv2.rectangle(img, (x, y), (x + w, y + h), LINE, 1, cv2.LINE_AA)


def section_title(img: np.ndarray, num: str, title: str, y: int) -> None:
    put(img, num, (70, y), 0.62, GOLD, 2)
    put(img, title, (132, y), 1.02, INK, 2)


def make_tables(m: dict) -> dict[str, pd.DataFrame]:
    evidence = pd.DataFrame(
        [
            ["Image classifier", "TCGA H&E / UNI LOTO", f"AUC {m['uni']['pooled_overall_auc']:.3f}", "Strong", "Core claim anchor"],
            [
                "Spatial RNA recovery",
                "GSE250521 Visium",
                f"slide-centered rho {m['g250_stage']['dm1_existing_slide_centered_rho']:.3f}",
                "Strong",
                "Cross-modal biology",
            ],
            [
                "External thyroid H&E",
                "GSE230424 Visium + raw H&E",
                f"sample-centered rho {m['g230']['top_he_sample_centered_rho']:.3f}",
                "Strong with QC caveat",
                "Independent thyroid support",
            ],
            [
                "Hotspot localization",
                "GSE250521 / GSE230424",
                f"lift {m['hotspot']['gse250521_top10_lift']:.1f}x / {m['hotspot']['gse230424_top10_lift']:.1f}x",
                "Strong",
                "Visual/editorial hook",
            ],
            [
                "Dose response",
                "Decile gradients",
                f"rho {m['decile']['gse250521_raw_decile_spearman']:.3f} / {m['decile']['gse230424_raw_decile_spearman']:.3f}",
                "Strong",
                "Reviewer robustness",
            ],
            [
                "Specificity",
                "Random modules + residuals",
                f"resid p {m['g250_spec']['actual_residual_empirical_p']:.3g}; G230 resid rho {m['g230_resid']['dm1_coord_qc_residual_target_he_rho']:.3f}",
                "Partial",
                "Keep boundary explicit",
            ],
            [
                "Negative control",
                "GSE248205 AITD",
                f"sample-centered rho {m['aitd']['ap_tls_he_sample_centered_rho']:.3f}",
                "Useful no-go",
                "Prevents overclaim",
            ],
        ],
        columns=["pillar", "dataset", "headline_metric", "strength", "impact_use"],
    )
    venue = pd.DataFrame(
        [
            ["Nature Communications", 72, "High risk", "Needs platform framing + ideally one more external WSI/Korean validation"],
            ["Cell Reports Medicine", 84, "Primary realistic high-impact", "Best match for translational pathology + biomarker recovery"],
            ["Communications Medicine", 82, "Strong backup", "Nature-family clinical/translational fit"],
            ["npj Precision Oncology", 78, "Strong backup", "Precision oncology biomarker + external validation angle"],
            ["Modern Pathology", 70, "Stable backup", "Computational/digital pathology scope fit"],
        ],
        columns=["venue", "readiness_score", "disposition", "what_would_raise_it"],
    )
    risks = pd.DataFrame(
        [
            ["Small TCGA classifier N", "Medium", "Lead with LOTO + spatial/external validation, not standalone classifier"],
            ["QC/tissue-density confounding in GSE230424", "High", "Co-locate residual and smoothness caveats with positive H&E result"],
            ["Biology vs platform novelty", "Medium", "Frame as cross-modal biomarker recovery, not a new causal thyroid mechanism"],
            ["No prospective Korean H&E yet", "High", "Make K2/Bundang the single revision-stage unlock"],
            ["Digital medicine implementation bar", "Medium", "Use presubmission for Cell Reports Medicine / Communications Medicine"],
        ],
        columns=["risk", "severity", "mitigation"],
    )
    unlock = pd.DataFrame(
        [
            ["K2/Bundang H&E external validation", "Highest", "Could move Nat Commun from long-shot to plausible"],
            ["True UNI locked rerun + release manifest", "High", "Removes foundation-model access objection"],
            ["Reader-friendly separated figure browser", "Done", "Raises editor comprehension and visual confidence"],
            ["Pre-registered revision protocol", "Medium", "Makes missing external validation less damaging"],
            ["Prospective clinical utility not available", "Boundary", "Do not sell as deployment-ready"],
        ],
        columns=["unlock", "priority", "impact"],
    )
    return {"evidence": evidence, "venue": venue, "risks": risks, "unlock": unlock}


def make_fig1(m: dict, tables: dict[str, pd.DataFrame]) -> Path:
    img = np.full((1280, 2100, 3), BG, dtype=np.uint8)
    put(img, "PAPER 2 IMPACT UPGRADE", (70, 82), 0.62, GOLD, 2)
    put(img, "Reframe as cross-modal biomarker recovery, not only a classifier", (70, 145), 1.16, INK, 2, 1450)
    put(img, "H&E image model -> spatial RNA/DM1 state -> external thyroid validation -> explicit caveat controls.", (70, 205), 0.62, MUTED, 1, 1350)

    cards = [
        ("UNI LOTO", f"AUC {m['uni']['pooled_overall_auc']:.3f}", "TCGA H&E image-DM1 anchor", GREEN),
        ("GSE250521", f"rho {m['g250_stage']['dm1_existing_slide_centered_rho']:.3f}", "slide-centered spatial RNA recovery", TEAL),
        ("GSE230424", f"rho {m['g230']['top_he_sample_centered_rho']:.3f}", "external thyroid H&E/Visium", GREEN),
        ("Hotspots", f"{m['hotspot']['gse230424_top10_lift']:.1f}x", "top-decile lift on external cohort", BLUE),
        ("Residual", f"rho {m['g230_resid']['dm1_coord_qc_residual_target_he_rho']:.3f}", "smaller signal beyond coord+QC", GOLD),
    ]
    x, y, w, h, gap = 70, 285, 370, 155, 32
    for i, (label, value, note, color) in enumerate(cards):
        x0 = x + i * (w + gap)
        box(img, (x0, y, x0 + w, y + h), PANEL)
        cv2.rectangle(img, (x0, y), (x0 + w, y + 8), color, -1)
        put(img, label.upper(), (x0 + 22, y + 42), 0.48, MUTED, 1)
        put(img, value, (x0 + 22, y + 95), 1.03, color, 2)
        put(img, note, (x0 + 22, y + 128), 0.43, INK, 1, w - 44)

    section_title(img, "01", "Venue Readiness Ladder", 535)
    venue = tables["venue"]
    y0 = 615
    for i, row in venue.iterrows():
        yy = y0 + i * 108
        color = GREEN if row.readiness_score >= 82 else GOLD if row.readiness_score >= 75 else BLUE
        box(img, (70, yy, 2025, yy + 82), PANEL)
        put(img, row.venue, (96, yy + 34), 0.68, INK, 2, 390)
        put(img, row.disposition, (515, yy + 34), 0.52, color, 2, 360)
        bar(img, 915, yy + 26, 360, 24, row.readiness_score / 100, color)
        put(img, f"{row.readiness_score}/100", (1305, yy + 46), 0.54, color, 2)
        put(img, row.what_would_raise_it, (1455, yy + 34), 0.46, MUTED, 1, 520)

    section_title(img, "02", "Impact Thesis For Editors", 1135)
    theses = [
        ("Not just AI", "The model recovers a molecular state from routine pathology and projects it into spatial transcriptomics."),
        ("Not just thyroid", "The paper tests whether histology carries latent transcriptomic biomarker structure."),
        ("Not overclaimed", "QC/smoothness caveats and AITD no-go are built into the claim boundary."),
    ]
    for i, (head, txt) in enumerate(theses):
        x0 = 70 + i * 670
        box(img, (x0, 1210, x0 + 620, 1390), PANEL2)
        put(img, head, (x0 + 24, 1255), 0.72, GOLD, 2)
        put(img, txt, (x0 + 24, 1305), 0.54, INK, 1, 555)

    path = OUT / "F01_paper2_impact_ladder.png"
    cv2.imwrite(str(path), img)
    return path


def make_fig2(m: dict, tables: dict[str, pd.DataFrame]) -> Path:
    img = np.full((1320, 2100, 3), BG, dtype=np.uint8)
    put(img, "EVIDENCE MOAT", (70, 82), 0.62, GOLD, 2)
    put(img, "Five independent reasons the story is larger than a small classifier", (70, 145), 1.08, INK, 2, 1420)

    evidence = tables["evidence"]
    x0, y0 = 70, 245
    colw = [310, 360, 370, 250, 560]
    headers = ["Pillar", "Dataset", "Headline", "Strength", "Impact use"]
    cx = x0
    for w, head in zip(colw, headers):
        box(img, (cx, y0, cx + w, y0 + 58), (42, 34, 22), GOLD)
        put(img, head, (cx + 12, y0 + 37), 0.48, GOLD, 2)
        cx += w
    for r, row in evidence.iterrows():
        yy = y0 + 58 + r * 132
        cx = x0
        for w, value in zip(colw, row.tolist()):
            box(img, (cx, yy, cx + w, yy + 132), PANEL if r % 2 == 0 else PANEL2)
            color = GREEN if value == "Strong" else GOLD if "Partial" in str(value) or "caveat" in str(value) else INK
            put(img, str(value), (cx + 12, yy + 35), 0.45, color, 1, w - 24)
            cx += w

    section_title(img, "03", "Claim Boundary That Helps, Not Hurts", 1245)
    items = [
        ("Use strongly", "H&E contains recoverable local DM1/RAI spatial signal."),
        ("Use carefully", "External thyroid support persists but QC/tissue structure is a major component."),
        ("Do not claim", "This alone proves causal MAPK/RAI mechanism or deployment-ready clinical action."),
    ]
    for i, (head, txt) in enumerate(items):
        x = 70 + i * 670
        color = GREEN if i == 0 else GOLD if i == 1 else RED
        box(img, (x, 1320, x + 620, 1510), PANEL2)
        cv2.rectangle(img, (x, 1320), (x + 620, 1328), color, -1)
        put(img, head, (x + 24, 1370), 0.72, color, 2)
        put(img, txt, (x + 24, 1425), 0.54, INK, 1, 555)

    path = OUT / "F02_paper2_evidence_moat.png"
    cv2.imwrite(str(path), img)
    return path


def make_fig3(tables: dict[str, pd.DataFrame]) -> Path:
    img = np.full((1280, 2000, 3), BG, dtype=np.uint8)
    put(img, "REVIEWER RISK REGISTER", (70, 82), 0.62, GOLD, 2)
    put(img, "Turn weaknesses into a visible control architecture", (70, 145), 1.1, INK, 2, 1320)

    risks = tables["risks"]
    y0 = 255
    for i, row in risks.iterrows():
        yy = y0 + i * 172
        severity_color = RED if row.severity == "High" else GOLD
        box(img, (70, yy, 1930, yy + 132), PANEL if i % 2 == 0 else PANEL2)
        put(img, row.risk, (105, yy + 44), 0.72, INK, 2, 520)
        put(img, row.severity, (690, yy + 44), 0.64, severity_color, 2)
        put(img, row.mitigation, (880, yy + 44), 0.55, MUTED, 1, 960)

    section_title(img, "04", "Single Highest-Impact Unlock", 1135)
    box(img, (70, 1210, 1930, 1480), PANEL2)
    put(img, "K2/Bundang H&E external validation", (110, 1275), 0.9, GOLD, 2)
    put(
        img,
        "This is the one extra dataset that changes the story from strong translational evidence to a much more credible Nature Communications attempt. Everything else is optimization.",
        (110, 1345),
        0.62,
        INK,
        1,
        1660,
    )
    path = OUT / "F03_paper2_reviewer_risk_register.png"
    cv2.imwrite(str(path), img)
    return path


def make_fig4(tables: dict[str, pd.DataFrame]) -> Path:
    img = np.full((1280, 2000, 3), BG, dtype=np.uint8)
    put(img, "NEXT VALIDATION UNLOCK MAP", (70, 82), 0.62, GOLD, 2)
    put(img, "What actually raises the ceiling from good paper to high-impact paper", (70, 145), 1.05, INK, 2, 1350)

    unlock = tables["unlock"]
    x0, y0 = 110, 285
    for i, row in unlock.iterrows():
        yy = y0 + i * 168
        color = GREEN if row.priority in {"Highest", "Done"} else GOLD if row.priority == "High" else RED if row.priority == "Boundary" else BLUE
        cv2.circle(img, (x0, yy + 44), 26, color, -1, cv2.LINE_AA)
        put(img, str(i + 1), (x0 - 9, yy + 54), 0.62, BG, 2)
        cv2.line(img, (x0, yy + 72), (x0, yy + 140), LINE, 2, cv2.LINE_AA)
        box(img, (180, yy, 1900, yy + 112), PANEL if i % 2 == 0 else PANEL2)
        put(img, row.unlock, (215, yy + 40), 0.7, INK, 2, 580)
        put(img, row.priority, (850, yy + 40), 0.58, color, 2, 180)
        put(img, row.impact, (1070, yy + 40), 0.52, MUTED, 1, 730)

    path = OUT / "F04_paper2_validation_unlock_map.png"
    cv2.imwrite(str(path), img)
    return path


def write_outputs(m: dict, tables: dict[str, pd.DataFrame], figs: list[Path]) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for name, df in tables.items():
        df.to_csv(OUT / f"paper2_impact_{name}.tsv", sep="\t", index=False)
    summary = {
        "created": "2026-05-10",
        "positioning": "cross-modal biomarker recovery from pathology into spatial transcriptomics",
        "recommended_primary": "Cell Reports Medicine presubmission inquiry",
        "high_risk_shot": "Nature Communications if K2/Bundang H&E external validation is added",
        "headline_metrics": {
            "uni_loto_auc": m["uni"]["pooled_overall_auc"],
            "gse250521_slide_centered_rho": m["g250_stage"]["dm1_existing_slide_centered_rho"],
            "gse230424_he_sample_centered_rho": m["g230"]["top_he_sample_centered_rho"],
            "gse230424_residual_target_rho": m["g230_resid"]["dm1_coord_qc_residual_target_he_rho"],
            "gse250521_hotspot_lift": m["hotspot"]["gse250521_top10_lift"],
            "gse230424_hotspot_lift": m["hotspot"]["gse230424_top10_lift"],
            "negative_control_gse248205_sample_centered_rho": m["aitd"]["ap_tls_he_sample_centered_rho"],
        },
        "figures": [p.name for p in figs],
    }
    (OUT / "PAPER2_IMPACT_UPGRADE_SUMMARY.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    md = [
        "# Paper 2 impact upgrade pack",
        "",
        "Positioning: cross-modal biomarker recovery from pathology into spatial transcriptomics.",
        "",
        "## Headline",
        f"- UNI LOTO AUC: {m['uni']['pooled_overall_auc']:.3f}.",
        f"- GSE250521 slide-centered rho: {m['g250_stage']['dm1_existing_slide_centered_rho']:.3f}.",
        f"- GSE230424 H&E sample-centered rho: {m['g230']['top_he_sample_centered_rho']:.3f}.",
        f"- Hotspot lift: {m['hotspot']['gse250521_top10_lift']:.1f}x / {m['hotspot']['gse230424_top10_lift']:.1f}x.",
        "",
        "## Recommendation",
        "- Primary: Cell Reports Medicine presubmission inquiry.",
        "- High-risk: Nature Communications only if K2/Bundang H&E external validation is secured.",
        "- Strong backup: Communications Medicine or npj Precision Oncology.",
        "",
        "## Files",
    ]
    for p in figs:
        md.append(f"- `{p.name}`")
    (OUT / "SUMMARY.md").write_text("\n".join(md) + "\n", encoding="utf-8")


def publish_assets(figs: list[Path]) -> None:
    ASSET_LOCAL.mkdir(parents=True, exist_ok=True)
    ASSET_LIVE.mkdir(parents=True, exist_ok=True)
    for fig in figs:
        shutil.copy2(fig, ASSET_LOCAL / fig.name)
        shutil.copy2(fig, ASSET_LIVE / fig.name)


def html_table(df: pd.DataFrame) -> str:
    rows = ["<table><thead><tr>" + "".join(f"<th>{c}</th>" for c in df.columns) + "</tr></thead><tbody>"]
    for _, row in df.iterrows():
        rows.append("<tr>" + "".join(f"<td>{row[c]}</td>" for c in df.columns) + "</tr>")
    rows.append("</tbody></table>")
    return "\n".join(rows)


def write_page(m: dict, tables: dict[str, pd.DataFrame], figs: list[Path]) -> None:
    fig_cards = "\n".join(
        f"""
        <button class="fig-option" data-src="assets/{ASSET_NAME}/{fig.name}" data-title="{fig.stem}">
          <img src="assets/{ASSET_NAME}/{fig.name}" alt="">
          <span>{fig.stem.replace('_', ' ')}</span>
        </button>
        """
        for fig in figs
    )
    html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Paper 2 Impact Upgrade</title>
<style>
:root {{ --bg:#0d1117; --panel:#151b23; --panel2:#1c2634; --ink:#e6edf3; --muted:#98a6ba; --line:#303846; --gold:#d6b25e; --teal:#63d5c4; --blue:#57a6ff; --red:#ff7b72; }}
* {{ box-sizing:border-box; }}
body {{ margin:0; background:var(--bg); color:var(--ink); font-family:"JetBrains Mono", ui-monospace, Menlo, monospace; line-height:1.55; }}
a {{ color:var(--blue); text-decoration:none; }}
a:hover {{ text-decoration:underline; }}
.hero {{ min-height:62vh; display:grid; align-items:end; padding:58px 5vw 42px; border-bottom:1px solid var(--line); background:linear-gradient(180deg, rgba(13,17,23,.16), rgba(13,17,23,.97)), url("assets/{ASSET_NAME}/{figs[1].name}") center / cover no-repeat; }}
.kicker {{ color:var(--gold); text-transform:uppercase; font-size:13px; font-weight:900; }}
h1 {{ margin:10px 0 14px; max-width:1180px; font-family:"Cormorant Garamond", Georgia, serif; font-size:clamp(42px,7vw,88px); line-height:.96; letter-spacing:0; }}
.lead {{ max-width:980px; color:#d8dee8; font-size:18px; }}
.stats {{ display:grid; grid-template-columns:repeat(6,minmax(0,1fr)); gap:10px; max-width:1200px; margin-top:24px; }}
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
.num {{ color:var(--gold); margin-right:10px; }}
.grid {{ display:grid; grid-template-columns:repeat(2,minmax(0,1fr)); gap:16px; }}
.card {{ border:1px solid var(--line); background:var(--panel); border-radius:8px; padding:18px; }}
.muted {{ color:var(--muted); }}
table {{ width:100%; border-collapse:collapse; margin:16px 0; font-size:13px; }}
th,td {{ border-bottom:1px solid var(--line); padding:10px 8px; text-align:left; vertical-align:top; }}
th {{ color:var(--gold); background:rgba(214,178,94,.06); }}
.viewer {{ display:grid; grid-template-columns:360px minmax(0,1fr); gap:16px; }}
.fig-list {{ display:grid; gap:10px; align-content:start; max-height:760px; overflow:auto; padding-right:6px; }}
.fig-option {{ border:1px solid var(--line); background:var(--panel); color:var(--ink); border-radius:8px; padding:10px; display:grid; grid-template-columns:92px minmax(0,1fr); gap:10px; text-align:left; cursor:pointer; font:inherit; }}
.fig-option.active,.fig-option:hover {{ border-color:var(--gold); background:#202938; }}
.fig-option img {{ width:92px; height:64px; object-fit:cover; border-radius:4px; display:block; }}
.stage {{ border:1px solid var(--line); background:#090b10; border-radius:8px; overflow:hidden; }}
.toolbar {{ display:flex; gap:8px; align-items:center; padding:12px; background:var(--panel); border-bottom:1px solid var(--line); }}
.toolbar strong {{ flex:1 1 auto; }}
button.control {{ border:1px solid var(--line); background:var(--panel2); color:var(--ink); border-radius:6px; padding:8px 10px; font:inherit; font-size:12px; cursor:pointer; }}
button.control:hover {{ border-color:var(--gold); color:var(--gold); }}
.canvas {{ height:min(78vh,840px); overflow:auto; padding:16px; background:#090b10; }}
.canvas img {{ display:block; width:100%; min-width:560px; max-width:none; height:auto; margin:0 auto; }}
code {{ background:var(--panel2); border:1px solid var(--line); padding:2px 5px; border-radius:4px; color:#d8dee8; }}
@media(max-width:980px) {{ .layout,.viewer {{ grid-template-columns:1fr; }} nav {{ position:static; max-height:none; display:grid; grid-template-columns:repeat(2,minmax(0,1fr)); }} .stats {{ grid-template-columns:repeat(2,minmax(0,1fr)); }} .grid {{ grid-template-columns:1fr; }} }}
</style>
</head>
<body>
<header class="hero">
  <div>
    <div class="kicker">Paper 2 · impact upgrade · 2026-05-10</div>
    <h1>Raise The Ceiling: From Classifier To Cross-Modal Biomarker Recovery</h1>
    <p class="lead">This pack reorganizes the evidence for high-impact submission: H&E image-DM1 classification, spatial transcriptomic recovery, external thyroid H&E validation, hotspot localization, and honest caveats in one editor-facing dossier.</p>
    <div class="stats">
      <div class="stat"><div class="v">{m['uni']['pooled_overall_auc']:.3f}</div><div class="l">UNI LOTO AUC</div></div>
      <div class="stat"><div class="v">{m['g250_stage']['dm1_existing_slide_centered_rho']:.3f}</div><div class="l">GSE250521 slide rho</div></div>
      <div class="stat"><div class="v">{m['g230']['top_he_sample_centered_rho']:.3f}</div><div class="l">GSE230424 H&E rho</div></div>
      <div class="stat"><div class="v">{m['hotspot']['gse230424_top10_lift']:.1f}x</div><div class="l">external hotspot lift</div></div>
      <div class="stat"><div class="v">81</div><div class="l">separated source figures</div></div>
      <div class="stat"><div class="v">CRM</div><div class="l">primary venue</div></div>
    </div>
  </div>
</header>
<div class="layout">
<nav>
  <a href="#verdict">01 Verdict</a>
  <a href="#figures">02 Figures</a>
  <a href="#evidence">03 Evidence</a>
  <a href="#venues">04 Venues</a>
  <a href="#risks">05 Risks</a>
  <a href="#unlocks">06 Unlocks</a>
  <a href="#sources">07 Sources</a>
</nav>
<main>
<section id="verdict">
  <h2><span class="num">01</span>Impact Verdict</h2>
  <div class="grid">
    <div class="card"><h3>Best Framing</h3><p class="muted">Cross-modal biomarker recovery: routine H&E recovers a DM1/RAI-linked spatial RNA state, validated across TCGA, GSE250521, and GSE230424 with caveats visible.</p></div>
    <div class="card"><h3>Best Submission Move</h3><p class="muted">Send Cell Reports Medicine presubmission first. Attempt Nature Communications only if K2/Bundang H&E validation can be added or promised as a concrete revision-stage unlock.</p></div>
  </div>
</section>
<section id="figures">
  <h2><span class="num">02</span>Impact Figure Viewer</h2>
  <div class="viewer">
    <div class="fig-list" id="figList">{fig_cards}</div>
    <div class="stage">
      <div class="toolbar"><strong id="figTitle">Figure</strong><button class="control" id="zoomOut">Zoom out</button><button class="control" id="zoomReset">Reset</button><button class="control" id="zoomIn">Zoom in</button></div>
      <div class="canvas" id="canvas"><img id="mainFig" src="assets/{ASSET_NAME}/{figs[0].name}" alt="Paper 2 impact figure"></div>
    </div>
  </div>
</section>
<section id="evidence"><h2><span class="num">03</span>Evidence Moat</h2>{html_table(tables['evidence'])}</section>
<section id="venues"><h2><span class="num">04</span>Venue Ladder</h2>{html_table(tables['venue'])}</section>
<section id="risks"><h2><span class="num">05</span>Reviewer Risk Register</h2>{html_table(tables['risks'])}</section>
<section id="unlocks"><h2><span class="num">06</span>Validation Unlocks</h2>{html_table(tables['unlock'])}</section>
<section id="sources">
  <h2><span class="num">07</span>Sources</h2>
  <table><tbody>
    <tr><td>Impact pack</td><td><code>{OUT.relative_to(ROOT)}</code></td></tr>
    <tr><td>Builder</td><td><code>project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/scripts/build_paper2_impact_upgrade_pack.py</code></td></tr>
    <tr><td>Separated figure browser</td><td><a href="paper2_cv2_visual_summary.html">paper2_cv2_visual_summary.html</a></td></tr>
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
function setZoom(next) {{ zoom = Math.min(3.5, Math.max(.7, next)); img.style.width = `${{zoom * 100}}%`; }}
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
    figs = [make_fig1(m, tables), make_fig2(m, tables), make_fig3(tables), make_fig4(tables)]
    write_outputs(m, tables, figs)
    publish_assets(figs)
    write_page(m, tables, figs)
    print(json.dumps({"out": str(OUT), "page": str(PAGE_LOCAL), "figures": [p.name for p in figs]}, indent=2))


if __name__ == "__main__":
    main()
