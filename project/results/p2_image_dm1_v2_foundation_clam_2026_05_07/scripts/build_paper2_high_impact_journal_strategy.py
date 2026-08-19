#!/usr/bin/env python3
"""Build Paper 2 high-impact journal routing dossier."""
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
OUT = SUPP / "paper2_high_impact_journal_strategy_2026_05_10"
HUB = ROOT / "project/papers_hub_2026_05_04"
LIVE = Path("/var/www/papers/papers_hub_2026_05_04")
ASSET = "paper2_high_impact_journal_strategy"
ASSET_LOCAL = HUB / "assets" / ASSET
ASSET_LIVE = LIVE / "assets" / ASSET
PAGE_LOCAL = HUB / "paper2_high_impact_journal_strategy.html"
PAGE_LIVE = LIVE / "paper2_high_impact_journal_strategy.html"

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


def read_json(path: Path) -> dict:
    with path.open() as handle:
        return json.load(handle)


def metrics() -> dict:
    return {
        "uni": read_json(SUPP / "audit_uni_loto/UNI_LOTO_SUMMARY.json"),
        "g250_stage": read_json(
            SUPP / "path2space_stage_generalization_controls_2026_05_09/GSE250521_STAGE_GENERALIZATION_SUMMARY.json"
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
    put(img, heading, (70, 142), 1.12, INK, 2, 1600)
    put(img, sub, (70, 205), 0.56, MUTED, 1, 1520)


def make_tables(m: dict) -> dict[str, pd.DataFrame]:
    journals = pd.DataFrame(
        [
            [
                1,
                "The Lancet Digital Health",
                "ultra-reach presub only",
                52,
                "Digital health AI journal; strongest if framed as practice-changing pathology-to-biomarker triage.",
                "Needs external clinical validation / prospective utility / implementation angle. K2/Bundang H&E is near-mandatory.",
                "Send presub inquiry only after submission-overview + K2 validation plan are attached.",
                "https://www.sciencedirect.com/journal/the-lancet-digital-health",
            ],
            [
                2,
                "Nature Medicine",
                "ultra-reach no-go now",
                35,
                "Translational/clinical medicine with direct human-health impact.",
                "Current paper lacks immediate clinical outcome impact or prospective clinical validation.",
                "Do not submit now; only revisit with prospective/clinical outcome evidence.",
                "https://www.nature.com/nm/aims",
            ],
            [
                3,
                "Nature Cancer",
                "biology reach presub only",
                45,
                "Cancer diagnostics, precision therapy, tumour heterogeneity, and systems/multi-omics scope.",
                "Needs deeper cancer-biology/mechanism or therapy-response advance beyond cross-modal prediction.",
                "Only presub if Paper 1 biology is already public and Paper 2 is framed as cancer-state recovery.",
                "https://www.nature.com/natcancer/aims",
            ],
            [
                4,
                "Nature Communications",
                "best NC-tier shot",
                72,
                "Multidisciplinary journal for important advances significant to specialists.",
                "Needs the story to read as cross-modal biomarker recovery, not a small AI classifier.",
                "Primary high-impact submission route if K2/Bundang validation is added; otherwise presub or high-risk submit.",
                "https://www.nature.com/ncomms/aims",
            ],
            [
                5,
                "Nature Biomedical Engineering",
                "engineering reach no-go now",
                42,
                "Applied biomedicine, health technology, AI, computational medicine, diagnostics.",
                "Method is not sufficiently new engineering; UNI/CLAM pipeline is applied rather than a platform invention.",
                "Do not use first unless a new validated method/platform layer is added.",
                "https://www.nature.com/natbiomedeng/aims",
            ],
            [
                6,
                "Cell Reports Medicine",
                "best practical high-impact",
                86,
                "Translational and clinical biomedical science; biomarker discovery and health technology fit.",
                "Small-N and clinical utility caveats must be visible.",
                "Submit/presub first if speed and realism matter; strongest current fit.",
                "https://www.sciencedirect.com/journal/cell-reports-medicine/about/insights",
            ],
            [
                7,
                "Communications Medicine",
                "Nature-family strong backup",
                82,
                "Clinical, translational, computational medicine, imaging, digital medicine scope.",
                "Needs medical-community relevance and strong evidence for conclusions.",
                "Good backup after Lancet DH / Nat Commun presub decisions.",
                "https://www.nature.com/commsmed/aims",
            ],
            [
                8,
                "npj Precision Oncology",
                "solid precision-oncology backup",
                78,
                "Precision oncology, biomarker studies, AI/computational oncology with external validation.",
                "Less high-impact than above, but current evidence fits well.",
                "Use after Nature/Cell-route triage or if editor feedback asks for tighter oncology venue.",
                "https://www.nature.com/npjprecisiononcology/aims",
            ],
            [
                9,
                "npj Digital Medicine",
                "not first",
                55,
                "Digital medicine and validated AI/ML models.",
                "Journal typically avoids off-the-shelf AI, purely observational, small preliminary studies.",
                "Avoid as first route unless clinical implementation/validation is strengthened.",
                "https://www.nature.com/npjdigitalmed/aims",
            ],
        ],
        columns=[
            "rank",
            "journal",
            "route",
            "readiness",
            "fit",
            "blocking_gap",
            "decision",
            "source_url",
        ],
    )
    presub = pd.DataFrame(
        [
            ["The Lancet Digital Health", "Question", "Is cross-modal H&E-to-spatial-RNA biomarker recovery sufficiently practice-changing for digital pathology?", "Attach K2/Bundang validation plan; ask before full submit."],
            ["The Lancet Digital Health", "Fact slots", f"UNI LOTO AUC {m['uni']['pooled_overall_auc']:.3f}; GSE250521 rho {m['g250_stage']['dm1_existing_slide_centered_rho']:.3f}; GSE230424 rho {m['g230']['top_he_sample_centered_rho']:.3f}; hotspot lift {m['hotspot']['gse230424_top10_lift']:.1f}x.", "Do not claim deployment-ready clinical utility."],
            ["Nature Communications", "Question", "Does the manuscript represent an important cross-modal biomarker-recovery advance for thyroid cancer specialists?", "Best route if K2/Bundang external H&E is added."],
            ["Nature Cancer", "Question", "Is there enough cancer-state biology beyond prediction?", "Only after Paper 1 public/preprint and mechanism/therapy relevance are foregrounded."],
            ["Cell Reports Medicine", "Question", "Is this a translational biomarker discovery + digital pathology article?", "Most realistic high-impact route now."],
        ],
        columns=["journal", "presub_component", "content_fact", "guardrail"],
    )
    action = pd.DataFrame(
        [
            ["Now", "Finish submission overview + journal dossier", "Done", "Use current web pages in advisor review."],
            ["Next 24-48h", "Prepare two presub packs", "Lancet DH + Nat Commun", "Fact scaffold only; author writes final voice."],
            ["This week", "Lock K2/Bundang H&E availability", "Highest", "Changes Nat Commun/Lancet DH from speculative to plausible."],
            ["If K2 yes", "Run external H&E validation", "Critical", "Add to Fig 6 and Table 1/2."],
            ["If K2 no", "Submit Cell Reports Medicine", "Pragmatic", "Do not wait indefinitely."],
        ],
        columns=["timing", "action", "priority", "reason"],
    )
    return {"journals": journals, "presub": presub, "action": action}


def draw_ladder(tables: dict[str, pd.DataFrame], m: dict) -> Path:
    img = np.full((1320, 2100, 3), BG, dtype=np.uint8)
    title(img, "High-impact routing", "Where Paper 2 Should Go, And When", "The route is not one journal. It is a staged decision: ultra-reach presub, NC-tier attempt, practical high-impact submission.")
    journals = tables["journals"]
    y0 = 300
    for i, row in journals.iterrows():
        yy = y0 + i * 112
        readiness = int(row.readiness)
        color = GREEN if readiness >= 82 else TEAL if readiness >= 72 else GOLD if readiness >= 52 else RED
        box(img, (70, yy, 2020, yy + 82), PANEL if i % 2 == 0 else PANEL2)
        put(img, f"{int(row['rank'])}", (95, yy + 51), 0.72, color, 2)
        put(img, row.journal, (150, yy + 37), 0.62, INK, 2, 390)
        put(img, row.route, (580, yy + 37), 0.48, color, 2, 310)
        cv2.rectangle(img, (940, yy + 28), (1260, yy + 52), (43, 50, 62), -1)
        cv2.rectangle(img, (940, yy + 28), (940 + int(320 * readiness / 100), yy + 52), color, -1)
        cv2.rectangle(img, (940, yy + 28), (1260, yy + 52), LINE, 1, cv2.LINE_AA)
        put(img, f"{readiness}/100", (1290, yy + 50), 0.5, color, 2)
        put(img, row.decision, (1430, yy + 34), 0.42, MUTED, 1, 520)

    box(img, (70, 1355, 2020, 1550), PANEL2)
    put(img, "Current best answer", (105, 1410), 0.82, GOLD, 2)
    put(
        img,
        "The Lancet Digital Health presub is the most ambitious play. Nature Communications is the best NC-tier route if K2/Bundang H&E external validation lands. Cell Reports Medicine is the strongest practical high-impact submission now.",
        (105, 1470),
        0.58,
        INK,
        1,
        1780,
    )
    path = OUT / "F01_high_impact_journal_ladder.png"
    cv2.imwrite(str(path), img)
    return path


def draw_decision_tree(tables: dict[str, pd.DataFrame]) -> Path:
    img = np.full((1320, 2100, 3), BG, dtype=np.uint8)
    title(img, "Submission decision tree", "Do Not Guess: Route By Validation State", "The highest-impact route depends on whether Korean/external H&E validation can be unlocked quickly.")
    nodes = [
        ("Start", "Paper 2 package ready\nclassifier + spatial + external H&E + caveats", (760, 310), GOLD),
        ("Q1", "Can K2/Bundang H&E validation be obtained now?", (760, 570), TEAL),
        ("YES", "Run external H&E validation\nthen send Nat Commun or Lancet DH presub", (330, 860), GREEN),
        ("NO", "Do not wait indefinitely\nsend Cell Reports Medicine presub/submission", (1190, 860), BLUE),
        ("Reach", "If K2 positive:\nLancet Digital Health presub\nNature Communications submission", (330, 1170), GREEN),
        ("Practical", "If no K2:\nCell Reports Medicine first\nCommunications Medicine backup", (1190, 1170), BLUE),
        ("Never first", "Avoid Nature Medicine/Nature Cancer/Nat BME as first full submission now", (760, 1500), RED),
    ]
    for _, text, (cx, cy), color in nodes:
        box(img, (cx - 300, cy - 82, cx + 300, cy + 82), PANEL2 if color != RED else (38, 27, 32), color)
        put(img, text, (cx - 265, cy - 22), 0.54, INK, 1, 530)
    arrows = [
        ((760, 392), (760, 488), GOLD),
        ((690, 640), (390, 780), GREEN),
        ((830, 640), (1130, 780), BLUE),
        ((330, 942), (330, 1088), GREEN),
        ((1190, 942), (1190, 1088), BLUE),
        ((330, 1252), (680, 1430), RED),
        ((1190, 1252), (840, 1430), RED),
    ]
    for a, b, color in arrows:
        cv2.arrowedLine(img, a, b, color, 3, cv2.LINE_AA, tipLength=0.08)
    path = OUT / "F02_submission_decision_tree.png"
    cv2.imwrite(str(path), img)
    return path


def draw_gap_map(tables: dict[str, pd.DataFrame]) -> Path:
    img = np.full((1320, 2100, 3), BG, dtype=np.uint8)
    title(img, "Gap-to-journal map", "What Raises Each Journal Route", "This converts reviewer/editor risk into a concrete next-data checklist.")
    rows = [
        ("Lancet DH", "Practice-changing clinical/digital validation", "K2/Bundang external H&E + clinical triage scenario", RED),
        ("Nature Medicine", "Immediate clinical outcome impact", "Prospective/clinical outcome data; not feasible now", RED),
        ("Nature Cancer", "Cancer biology/mechanism depth", "Tie Paper 1 biology + therapy relevance + independent validation", GOLD),
        ("Nature Comms", "Specialist-significant cross-modal advance", "K2/Bundang H&E validation + strong figure flow", GREEN),
        ("Nat Biomed Eng", "Engineering invention", "New validated method/platform beyond applied UNI/CLAM", GOLD),
        ("Cell Rep Med", "Translational biomedical fit", "Current data sufficient with transparent caveats", GREEN),
    ]
    y0 = 310
    for i, (journal, gap, unlock, color) in enumerate(rows):
        yy = y0 + i * 175
        box(img, (70, yy, 2020, yy + 128), PANEL if i % 2 == 0 else PANEL2)
        cv2.rectangle(img, (70, yy), (78, yy + 128), color, -1)
        put(img, journal, (110, yy + 50), 0.68, color, 2, 280)
        put(img, gap, (430, yy + 45), 0.55, INK, 1, 620)
        put(img, unlock, (1120, yy + 45), 0.55, MUTED, 1, 800)
    path = OUT / "F03_gap_to_journal_unlock_map.png"
    cv2.imwrite(str(path), img)
    return path


def draw_presub_pack(tables: dict[str, pd.DataFrame], m: dict) -> Path:
    img = np.full((1320, 2100, 3), BG, dtype=np.uint8)
    title(img, "Presubmission packet", "What To Send Before Full Submission", "This is a fact scaffold. The author should write the final voice.")
    cards = [
        ("1. One-sentence factual thesis", "Routine H&E can recover a DM1/RAI-linked spatial transcriptomic biomarker state in thyroid cancer."),
        ("2. Three-number evidence line", f"UNI LOTO AUC {m['uni']['pooled_overall_auc']:.3f}; GSE250521 slide-centered rho {m['g250_stage']['dm1_existing_slide_centered_rho']:.3f}; GSE230424 H&E rho {m['g230']['top_he_sample_centered_rho']:.3f}."),
        ("3. Validation statement", f"External thyroid Visium/H&E support plus residual-target rho {m['g230_resid']['dm1_coord_qc_residual_target_he_rho']:.3f}; hotspots lift {m['hotspot']['gse230424_top10_lift']:.1f}x."),
        ("4. Honest caveat", "QC/smoothness contributes to the external H&E signal; negative-control AITD cohort prevents immune-axis overclaim."),
        ("5. Editor ask", "Ask whether the cross-modal biomarker-recovery framing fits before full submission."),
        ("6. Attachment set", "Submission overview page, impact dossier page, separated figure browser, and K2/Bundang validation plan."),
    ]
    for i, (head, body) in enumerate(cards):
        x = 70 + (i % 2) * 980
        y = 305 + (i // 2) * 340
        color = [GOLD, GREEN, TEAL, BLUE, RED, GOLD][i]
        box(img, (x, y, x + 920, y + 255), PANEL2 if i % 2 else PANEL)
        cv2.rectangle(img, (x, y), (x + 920, y + 8), color, -1)
        put(img, head, (x + 26, y + 55), 0.66, color, 2, 820)
        put(img, body, (x + 26, y + 120), 0.54, INK, 1, 840)
    path = OUT / "F04_presubmission_fact_packet.png"
    cv2.imwrite(str(path), img)
    return path


def write_outputs(tables: dict[str, pd.DataFrame], figs: list[Path], m: dict) -> None:
    for name, df in tables.items():
        df.to_csv(OUT / f"paper2_high_impact_{name}.tsv", sep="\t", index=False)
    summary = {
        "created": "2026-05-10",
        "recommendation": {
            "moonshot_presub": "The Lancet Digital Health",
            "best_nc_tier": "Nature Communications with K2/Bundang H&E validation",
            "practical_high_impact": "Cell Reports Medicine",
            "backup": "Communications Medicine / npj Precision Oncology",
            "avoid_first_now": ["Nature Medicine", "Nature Cancer", "Nature Biomedical Engineering", "npj Digital Medicine"],
        },
        "headline_metrics": {
            "uni_loto_auc": m["uni"]["pooled_overall_auc"],
            "gse250521_slide_centered_rho": m["g250_stage"]["dm1_existing_slide_centered_rho"],
            "gse230424_he_rho": m["g230"]["top_he_sample_centered_rho"],
            "gse230424_residual_target_rho": m["g230_resid"]["dm1_coord_qc_residual_target_he_rho"],
            "gse230424_hotspot_lift": m["hotspot"]["gse230424_top10_lift"],
            "aitd_negative_control_rho": m["aitd"]["ap_tls_he_sample_centered_rho"],
        },
        "figures": [p.name for p in figs],
    }
    (OUT / "PAPER2_HIGH_IMPACT_JOURNAL_STRATEGY.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    md = [
        "# Paper 2 high-impact journal strategy",
        "",
        "Recommendation:",
        "- Moonshot presub: The Lancet Digital Health.",
        "- Best NC-tier route: Nature Communications after K2/Bundang H&E validation.",
        "- Practical high-impact route now: Cell Reports Medicine.",
        "- Backup: Communications Medicine / npj Precision Oncology.",
        "",
        "Protected voice note: this is a strategy and fact scaffold, not final cover-letter prose.",
        "",
        "Figures:",
    ]
    for fig in figs:
        md.append(f"- `{fig.name}`")
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
        cells = []
        for c in df.columns:
            val = row[c]
            if c == "source_url":
                cells.append(f'<td><a href="{val}">{val}</a></td>')
            else:
                cells.append(f"<td>{val}</td>")
        out.append("<tr>" + "".join(cells) + "</tr>")
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
<title>Paper 2 High-Impact Journal Strategy</title>
<style>
:root {{ --bg:#0d1117; --panel:#151b23; --panel2:#1c2634; --ink:#e6edf3; --muted:#98a6ba; --line:#303846; --gold:#d6b25e; --blue:#57a6ff; }}
* {{ box-sizing:border-box; }}
body {{ margin:0; background:var(--bg); color:var(--ink); font-family:"JetBrains Mono", ui-monospace, Menlo, monospace; line-height:1.55; }}
a {{ color:var(--blue); text-decoration:none; }}
a:hover {{ text-decoration:underline; }}
.hero {{ min-height:62vh; display:grid; align-items:end; padding:58px 5vw 42px; border-bottom:1px solid var(--line); background:linear-gradient(180deg, rgba(13,17,23,.18), rgba(13,17,23,.97)), url("assets/{ASSET}/{figs[0].name}") center / cover no-repeat; }}
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
.num {{ color:var(--gold); margin-right:10px; }}
.grid {{ display:grid; grid-template-columns:repeat(3,minmax(0,1fr)); gap:16px; }}
.card {{ border:1px solid var(--line); background:var(--panel); border-radius:8px; padding:18px; }}
.muted {{ color:var(--muted); }}
table {{ width:100%; border-collapse:collapse; margin:16px 0; font-size:13px; }}
th,td {{ border-bottom:1px solid var(--line); padding:10px 8px; text-align:left; vertical-align:top; }}
th {{ color:var(--gold); background:rgba(214,178,94,.06); }}
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
@media(max-width:980px) {{ .layout,.viewer {{ grid-template-columns:1fr; }} nav {{ position:static; max-height:none; display:grid; grid-template-columns:repeat(2,minmax(0,1fr)); }} .stats,.grid {{ grid-template-columns:1fr; }} }}
</style>
</head>
<body>
<header class="hero">
  <div>
    <div class="kicker">Paper 2 · high-impact route · 2026-05-10</div>
    <h1>Moonshot Presub, NC-Tier Attempt, Practical Landing</h1>
    <p class="lead">A staged venue strategy for Paper 2 after reframing as H&E-to-spatial-RNA biomarker recovery.</p>
    <div class="stats">
      <div class="stat"><div class="v">LDH</div><div class="l">moonshot presub</div></div>
      <div class="stat"><div class="v">NC</div><div class="l">K2-unlocked shot</div></div>
      <div class="stat"><div class="v">CRM</div><div class="l">best practical route</div></div>
      <div class="stat"><div class="v">{m['uni']['pooled_overall_auc']:.3f}</div><div class="l">UNI LOTO AUC</div></div>
      <div class="stat"><div class="v">{m['g230']['top_he_sample_centered_rho']:.3f}</div><div class="l">external H&E rho</div></div>
    </div>
  </div>
</header>
<div class="layout">
<nav>
  <a href="#verdict">01 Verdict</a>
  <a href="#figures">02 Strategy Figures</a>
  <a href="#journals">03 Journal Matrix</a>
  <a href="#presub">04 Presub Scaffold</a>
  <a href="#action">05 Action Plan</a>
  <a href="#sources">06 Sources</a>
</nav>
<main>
<section id="verdict">
  <h2><span class="num">01</span>Venue Verdict</h2>
  <div class="grid">
    <div class="card"><h3>Moonshot</h3><p class="muted">The Lancet Digital Health presubmission only. Full submission requires stronger external/clinical validation.</p></div>
    <div class="card"><h3>Best NC-Tier</h3><p class="muted">Nature Communications if K2/Bundang H&E validation is added or at least made concrete.</p></div>
    <div class="card"><h3>Best Practical</h3><p class="muted">Cell Reports Medicine is the most realistic high-impact target with the current data.</p></div>
  </div>
</section>
<section id="figures">
  <h2><span class="num">02</span>Strategy Figure Viewer</h2>
  <div class="viewer">
    <div class="fig-list" id="figList">{fig_buttons}</div>
    <div class="stage">
      <div class="toolbar"><strong id="figTitle">Figure</strong><button class="control" id="zoomOut">Zoom out</button><button class="control" id="zoomReset">Reset</button><button class="control" id="zoomIn">Zoom in</button></div>
      <div class="canvas" id="canvas"><img id="mainFig" src="assets/{ASSET}/{figs[0].name}" alt="Paper 2 strategy figure"></div>
    </div>
  </div>
</section>
<section id="journals"><h2><span class="num">03</span>Journal Matrix</h2>{html_table(tables['journals'])}</section>
<section id="presub"><h2><span class="num">04</span>Presubmission Fact Scaffold</h2>{html_table(tables['presub'])}</section>
<section id="action"><h2><span class="num">05</span>Action Plan</h2>{html_table(tables['action'])}</section>
<section id="sources">
  <h2><span class="num">06</span>Local Sources</h2>
  <table><tbody>
    <tr><td>Output directory</td><td><code>{OUT.relative_to(ROOT)}</code></td></tr>
    <tr><td>Builder script</td><td><code>project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/scripts/build_paper2_high_impact_journal_strategy.py</code></td></tr>
    <tr><td>Submission overview</td><td><a href="paper2_submission_overview.html">paper2_submission_overview.html</a></td></tr>
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
        draw_ladder(tables, m),
        draw_decision_tree(tables),
        draw_gap_map(tables),
        draw_presub_pack(tables, m),
    ]
    write_outputs(tables, figs, m)
    publish(figs)
    write_page(tables, figs, m)
    print(json.dumps({"out": str(OUT), "page": str(PAGE_LOCAL), "figures": [p.name for p in figs]}, indent=2))


if __name__ == "__main__":
    main()
