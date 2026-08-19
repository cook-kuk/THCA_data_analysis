from __future__ import annotations

import csv
import html
import shutil
from pathlib import Path


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
SRC = ROOT / "project/results/therapy_spectrum_promising_master_2026_05_11/therapy_spectrum_promising_master.tsv"
RESULT = ROOT / "project/results/car_t_bundle_2026_05_11"
HUB = ROOT / "project/papers_hub_2026_05_04"
LIVE_HUB = Path("/var/www/papers/papers_hub_2026_05_04")
ASSET_DIR = HUB / "assets/car_t_bundle_2026_05_11"
LIVE_ASSET_DIR = LIVE_HUB / "assets/car_t_bundle_2026_05_11"


CAR_T_ACCESSIONS = {
    "GSE290722",
    "GSE313971",
    "GSE308384",
    "GSE236468",
    "GSE304795",
    "GSE312384",
    "GSE292859",
    "GSE284026",
}


TOPIC_LINKS = {
    "GSE290722": "car_t_durability_2026_05_11.html",
    "GSE313971": "car_engineering_potency_2026_05_11.html",
    "GSE308384": "car_engineering_potency_2026_05_11.html",
    "GSE236468": "car_engineering_potency_2026_05_11.html",
    "GSE304795": "car_engineering_potency_2026_05_11.html",
    "GSE312384": "car_engineering_potency_2026_05_11.html",
    "GSE292859": "car_engineering_potency_2026_05_11.html",
    "GSE284026": "car_engineering_potency_2026_05_11.html",
}


def esc(s: object) -> str:
    return html.escape(str(s), quote=True)


def ensure_dirs() -> None:
    RESULT.mkdir(parents=True, exist_ok=True)
    ASSET_DIR.mkdir(parents=True, exist_ok=True)
    LIVE_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    HUB.mkdir(parents=True, exist_ok=True)
    LIVE_HUB.mkdir(parents=True, exist_ok=True)


def load_rows() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with SRC.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row["accession"] in CAR_T_ACCESSIONS:
                rows.append(row)
    rows.sort(key=lambda r: (-int(r["score_total"]), -int(r["paperability"]), r["accession"]))
    return rows


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    headers = [
        "accession",
        "title",
        "modality",
        "status",
        "score_total",
        "paperability",
        "validation",
        "bd",
        "tier",
        "class",
        "best_task",
        "why",
        "caveat",
        "source",
        "topic_link",
    ]
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=headers, delimiter="\t")
        w.writeheader()
        for row in rows:
            out = dict(row)
            out["class"] = "car-t"
            out["topic_link"] = TOPIC_LINKS.get(row["accession"], "")
            w.writerow({k: out.get(k, "") for k in headers})


def write_summary(path: Path, rows: list[dict[str, str]]) -> None:
    core = [r for r in rows if r["tier"] == "TRIPLE"]
    dual = [r for r in rows if r["tier"] == "DUAL"]
    single = [r for r in rows if r["tier"] == "SINGLE"]
    lines = [
        "# CAR-T Bundle",
        "",
        "Date: 2026-05-11",
        "",
        "This bundle merges the CAR-T-related items from the therapy-spectrum master into one sendable report.",
        "",
        "Snapshot:",
        f"- Rows: {len(rows)}",
        f"- Core / dual / single: {len(core)} / {len(dual)} / {len(single)}",
        "",
        "Included themes:",
        "- CAR-T response durability",
        "- CAR-T design / signaling / engineering",
        "- CAR-T manufacturing / potency / QC",
        "- CAR-T solid tumor support",
        "",
        "Rule:",
        "- Keep this to CAR-T only; CAR-NKT, TIL, and blinatumomab are excluded.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def table(rows: list[dict[str, str]]) -> str:
    return "\n".join(
        (
            f"<tr>"
            f"<td><a href='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={esc(r['accession'])}'>{esc(r['accession'])}</a></td>"
            f"<td>{esc(r['modality'])}</td>"
            f"<td><span class='tag {esc(r['tier']).lower()}'>{esc(r['tier'])}</span></td>"
            f"<td>{esc(r['score_total'])}</td>"
            f"<td>{esc(r['best_task'])}</td>"
            f"<td>{esc(r['why'])}</td>"
            f"<td>{esc(r['caveat'])}</td>"
            f"<td>{link}</td>"
            f"</tr>"
        )
        for r in rows
        for link in [f"<a href='{esc(r['topic_link'])}'>{esc(r['topic_link'])}</a>" if r.get("topic_link") else ""]
    )


def build_html(rows: list[dict[str, str]]) -> str:
    core = [r for r in rows if r["tier"] == "TRIPLE"]
    top = rows[:6]
    cards = "".join(
        f"<div class='card'><h3>{i+1}. {esc(r['accession'])}</h3><p><b>{esc(r['modality'])}</b></p><p>{esc(r['tier'])} | Score {esc(r['score_total'])}</p><p>{esc(r['best_task'])}</p><p class='small'>{esc(r['why'])}</p></div>"
        for i, r in enumerate(top)
    )
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8" />
<meta name="viewport" content="width=device-width,initial-scale=1.0" />
<title>CAR-T Bundle · Sendable report</title>
<style>
:root{{--bg:#061019;--panel:#0f1726;--line:#22324b;--ink:#edf3fb;--muted:#9aabc0;--gold:#ffd28a;--teal:#35d39d;--blue:#7eb6ff;--red:#ff8a6b}}
*{{box-sizing:border-box}}
html,body{{margin:0;background:var(--bg);color:var(--ink);font:14px/1.55 Inter,"Noto Sans KR",sans-serif}}
a{{color:var(--teal);text-decoration:none}} a:hover{{color:var(--gold);text-decoration:underline}}
code{{font-family:"JetBrains Mono",monospace;font-size:12px}}
h1,h2,h3{{font-family:"Cormorant Garamond","Newsreader",serif;letter-spacing:-.02em;color:#fff7dc}}
.hero{{padding:64px 34px 38px;background:radial-gradient(circle at top left,#172946 0%,#0b1220 60%,#050810 100%);border-bottom:1px solid var(--line)}}
.hero-inner{{max-width:1440px;margin:0 auto}}
.kicker{{font:700 11px "JetBrains Mono",monospace;color:var(--gold);letter-spacing:.24em;text-transform:uppercase}}
h1{{font-size:64px;line-height:.96;margin:12px 0}}
.lead{{max-width:1040px;color:#cad5e5;font:18px/1.65 "Newsreader",serif}}
.stats{{display:grid;grid-template-columns:repeat(5,1fr);gap:10px;margin-top:24px}}
.stat{{background:rgba(255,255,255,.04);border:1px solid rgba(255,210,138,.18);padding:13px 14px;border-radius:10px}}
.stat b{{display:block;color:var(--gold);font:700 28px "Cormorant Garamond",serif;line-height:1}}
.stat span{{display:block;color:var(--muted);font:700 10px "JetBrains Mono",monospace;text-transform:uppercase;letter-spacing:.08em;margin-top:4px}}
.wrap{{max-width:1440px;margin:0 auto;padding:0 34px 80px}}
section{{padding:30px 0;border-bottom:1px solid var(--line)}}
h2{{font-size:36px;margin:0 0 6px}}
h2 .num{{font:700 13px "JetBrains Mono",monospace;color:var(--gold);letter-spacing:.18em;margin-right:10px}}
.sub{{font:700 10px "JetBrains Mono",monospace;color:var(--muted);letter-spacing:.12em;text-transform:uppercase;margin-bottom:16px}}
.grid{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:14px}}
.card{{background:var(--panel);border:1px solid var(--line);border-radius:10px;padding:14px 16px}}
.card h3{{margin:0 0 8px}}
.t{{width:100%;border-collapse:collapse;margin:12px 0 16px;font-size:12px}}
.t th,.t td{{border:1px solid var(--line);padding:8px 9px;vertical-align:top}}
.t th{{background:#16243a;color:var(--gold);font:700 10px "JetBrains Mono",monospace;letter-spacing:.05em;text-transform:uppercase;text-align:left}}
.t tr:nth-child(even) td{{background:rgba(255,255,255,.02)}}
.tag{{display:inline-block;padding:2px 8px;border-radius:99px;font:700 10px "JetBrains Mono",monospace;letter-spacing:.06em;text-transform:uppercase;border:1px solid transparent}}
.tag.triple{{background:rgba(53,211,157,.12);border-color:rgba(53,211,157,.28);color:var(--teal)}}
.tag.dual{{background:rgba(126,182,255,.12);border-color:rgba(126,182,255,.28);color:var(--blue)}}
.tag.single{{background:rgba(255,210,138,.12);border-color:rgba(255,210,138,.28);color:var(--gold)}}
.box{{background:var(--panel);border:1px solid var(--line);border-left:4px solid var(--gold);border-radius:0 8px 8px 0;padding:14px 16px;margin:14px 0}}
.box.good{{border-left-color:var(--teal);background:#0a1a16}}
.box.warn{{border-left-color:var(--red);background:#1d1310}}
.small{{color:var(--muted);font-size:12px}}
@media(max-width:1100px){{.grid{{grid-template-columns:1fr 1fr}}}}
@media(max-width:980px){{.stats{{grid-template-columns:repeat(2,1fr)}}.grid{{grid-template-columns:1fr}}h1{{font-size:40px}}.wrap{{padding:0 20px 60px}}}}
</style>
</head>
<body>
<header class="hero">
  <div class="hero-inner">
    <div class="kicker">CAR-T bundle · 2026-05-11</div>
    <h1>CAR-T Only<br/><em>Sendable Bundle</em></h1>
    <p class="lead">
      This is the CAR-T-only subset pulled out of the broader therapy-spectrum work.
      It excludes CAR-NKT, TIL, blinatumomab, and non-CAR-T therapy families.
    </p>
    <div style="margin-top:18px;color:#9aabc0;font:12px 'JetBrains Mono',monospace;">
      <a href="index.html" style="color:#35d39d">Hub index</a> ·
      <a href="therapy_spectrum_promising_master_2026_05_11.html" style="color:#c79cff">promising master</a> ·
      <a href="car_t_durability_2026_05_11.html" style="color:#35d39d">CAR-T durability</a> ·
      <a href="car_engineering_potency_2026_05_11.html" style="color:#7eb6ff">CAR engineering / potency</a>
    </div>
    <div class="stats">
      <div class="stat"><b>{len(rows)}</b><span>CAR-T rows</span></div>
      <div class="stat"><b>{len(core)}</b><span>Core</span></div>
      <div class="stat"><b>{len([r for r in rows if r['tier']=='DUAL'])}</b><span>Dual</span></div>
      <div class="stat"><b>{len([r for r in rows if r['tier']=='SINGLE'])}</b><span>Single</span></div>
      <div class="stat"><b>GEO</b><span>Sendable report</span></div>
    </div>
  </div>
</header>
<div class="wrap">
<section>
  <h2><span class="num">01</span>TL;DR</h2>
  <div class="sub">What to send if you only want CAR-T</div>
  <div class="grid">{cards}</div>
  <div class="box good"><b>Best send order:</b> CAR-T durability, CAR-T engineering/potency, CAR-T signaling/construct design, then manufacturing and solid-tumor support.</div>
</section>

<section>
  <h2><span class="num">02</span>Full Table</h2>
  <div class="sub">CAR-T-related studies only</div>
  <table class="t">
    <thead><tr><th>Accession</th><th>Modality</th><th>Tier</th><th>Score</th><th>Best task</th><th>Why</th><th>Caveat</th><th>Topic page</th></tr></thead>
    <tbody>{table(rows)}</tbody>
  </table>
  <div class="box warn"><b>Boundary.</b> This bundle intentionally excludes CAR-NKT and TIL even though they are cell therapy, because you asked for CAR-T only.</div>
</section>
</div>
</body>
</html>"""


def main() -> None:
    ensure_dirs()
    rows = load_rows()
    tsv = RESULT / "car_t_bundle.tsv"
    summary = RESULT / "SUMMARY.md"
    html_path = HUB / "car_t_bundle_2026_05_11.html"

    write_tsv(tsv, rows)
    write_summary(summary, rows)
    html_path.write_text(build_html(rows), encoding="utf-8")

    shutil.copy2(tsv, ASSET_DIR / tsv.name)
    shutil.copy2(summary, ASSET_DIR / summary.name)
    shutil.copy2(tsv, LIVE_ASSET_DIR / tsv.name)
    shutil.copy2(summary, LIVE_ASSET_DIR / summary.name)
    shutil.copy2(html_path, LIVE_HUB / html_path.name)


if __name__ == "__main__":
    main()
