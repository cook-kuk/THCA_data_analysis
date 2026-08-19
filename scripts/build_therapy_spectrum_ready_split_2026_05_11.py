from __future__ import annotations

import csv
import html
import shutil
from collections import Counter
from pathlib import Path


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
SRC = ROOT / "project/results/therapy_spectrum_go_shortlist_2026_05_11/therapy_spectrum_go_shortlist.tsv"
RESULT = ROOT / "project/results/therapy_spectrum_ready_split_2026_05_11"
HUB = ROOT / "project/papers_hub_2026_05_04"
LIVE_HUB = Path("/var/www/papers/papers_hub_2026_05_04")
ASSET_DIR = HUB / "assets/therapy_spectrum_ready_split_2026_05_11"
LIVE_ASSET_DIR = LIVE_HUB / "assets/therapy_spectrum_ready_split_2026_05_11"


def esc(s: object) -> str:
    return html.escape(str(s), quote=True)


def ensure_dirs() -> None:
    RESULT.mkdir(parents=True, exist_ok=True)
    ASSET_DIR.mkdir(parents=True, exist_ok=True)
    LIVE_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    HUB.mkdir(parents=True, exist_ok=True)
    LIVE_HUB.mkdir(parents=True, exist_ok=True)


def read_rows() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with SRC.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            row["paperable"] = int(row["paperability"]) >= 5 and int(row["task_fit"]) >= 5
            row["validation_ready"] = int(row["validation"]) >= 5
            row["bd_ready"] = int(row["bd"]) >= 5
            row["flags"] = sum([row["paperable"], row["validation_ready"], row["bd_ready"]])
            if row["flags"] == 3:
                row["tier"] = "TRIPLE"
            elif row["flags"] == 2:
                row["tier"] = "DUAL"
            else:
                row["tier"] = "SINGLE"
            rows.append(row)
    rows.sort(key=lambda r: (-int(r["score_total"]), -int(r["flags"]), r["accession"]))
    return rows


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    headers = [
        "accession",
        "title",
        "modality",
        "status",
        "score_total",
        "paperability",
        "validation",
        "bd",
        "paperable",
        "validation_ready",
        "bd_ready",
        "flags",
        "tier",
        "best_task",
        "why",
        "caveat",
        "source",
    ]
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=headers, delimiter="\t")
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in headers})


def write_summary(path: Path, rows: list[dict[str, object]]) -> None:
    flags = Counter(str(row["flags"]) for row in rows)
    triple = sum(1 for row in rows if row["flags"] == 3)
    dual = sum(1 for row in rows if row["flags"] == 2)
    single = sum(1 for row in rows if row["flags"] == 1)
    lines = [
        "# Therapy Spectrum Ready Split",
        "",
        "Date: 2026-05-11",
        "",
        "This split tags GO shortlist records by whether they are paperable, validation-ready, and/or BD-ready.",
        "",
        "Snapshot:",
        f"- Rows: {len(rows)}",
        f"- Triple-ready: {triple}",
        f"- Dual-ready: {dual}",
        f"- Single-ready: {single}",
        f"- Flag counts: {', '.join(f'{k}={v}' for k, v in sorted(flags.items(), reverse=True))}",
        "",
        "Rule:",
        "- Paperable = strong mechanistic story and clear figure/decision-artifact value.",
        "- Validation-ready = strong human-relevant or orthogonal confirmation path.",
        "- BD-ready = direct translation, companion biomarker, or product design value.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def rows_table(rows: list[dict[str, object]]) -> str:
    return "\n".join(
        f"<tr>"
        f"<td><a href='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={esc(r['accession'])}'>{esc(r['accession'])}</a></td>"
        f"<td>{esc(r['modality'])}</td>"
        f"<td>{esc(r['score_total'])}</td>"
        f"<td><span class='tag {'triple' if r['flags']==3 else 'dual' if r['flags']==2 else 'single'}'>{esc(r['tier'])}</span></td>"
        f"<td>{'Y' if r['paperable'] else 'N'}</td>"
        f"<td>{'Y' if r['validation_ready'] else 'N'}</td>"
        f"<td>{'Y' if r['bd_ready'] else 'N'}</td>"
        f"<td>{esc(r['best_task'])}</td>"
        f"<td>{esc(r['why'])}</td>"
        f"</tr>"
        for r in rows
    )


def build_html(rows: list[dict[str, object]]) -> str:
    paperable = [r for r in rows if r["paperable"]]
    validation = [r for r in rows if r["validation_ready"]]
    bd = [r for r in rows if r["bd_ready"]]
    top = rows[:6]
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8" />
<meta name="viewport" content="width=device-width,initial-scale=1.0" />
<title>Therapy Spectrum Ready Split · paper / validation / BD</title>
<style>
:root{{--bg:#071019;--panel:#101a2a;--line:#22324b;--ink:#edf3fb;--muted:#9aabc0;--gold:#ffd28a;--teal:#35d39d;--blue:#7eb6ff;--red:#ff8a6b}}
*{{box-sizing:border-box}}
html,body{{margin:0;background:var(--bg);color:var(--ink);font:14px/1.55 Inter,"Noto Sans KR",sans-serif}}
a{{color:var(--teal);text-decoration:none}} a:hover{{color:var(--gold);text-decoration:underline}}
h1,h2,h3{{font-family:"Cormorant Garamond","Newsreader",serif;letter-spacing:-.02em;color:#fff7dc}}
code{{font-family:"JetBrains Mono",monospace;font-size:12px}}
.hero{{padding:64px 34px 38px;background:radial-gradient(circle at top left,#172946 0%,#0b1220 60%,#050810 100%);border-bottom:1px solid var(--line)}}
.hero-inner{{max-width:1380px;margin:0 auto}}
.kicker{{font:700 11px "JetBrains Mono",monospace;color:var(--gold);letter-spacing:.24em;text-transform:uppercase}}
h1{{font-size:60px;line-height:.98;margin:12px 0}}
.lead{{max-width:1060px;color:#cad5e5;font:18px/1.65 "Newsreader",serif}}
.stats{{display:grid;grid-template-columns:repeat(5,1fr);gap:10px;margin-top:24px}}
.stat{{background:rgba(255,255,255,.04);border:1px solid rgba(255,210,138,.18);padding:13px 14px;border-radius:10px}}
.stat b{{display:block;color:var(--gold);font:700 28px "Cormorant Garamond",serif;line-height:1}}
.stat span{{display:block;color:var(--muted);font:700 10px "JetBrains Mono",monospace;text-transform:uppercase;letter-spacing:.08em;margin-top:4px}}
.wrap{{max-width:1380px;margin:0 auto;display:grid;grid-template-columns:260px 1fr;gap:34px;padding:0 34px 80px}}
.toc{{position:sticky;top:0;align-self:start;max-height:100vh;overflow:auto;border-right:1px solid var(--line);padding:28px 18px 28px 0}}
.toc h4{{font:700 10px "JetBrains Mono",monospace;color:var(--gold);letter-spacing:.18em;text-transform:uppercase;margin:0 0 12px}}
.toc a{{display:block;color:#cbd5e5;padding:5px 0;font-size:12px}}
main{{min-width:0;padding-top:24px}}
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
@media(max-width:1100px){{.grid{{grid-template-columns:1fr 1fr}}}}
@media(max-width:980px){{.wrap{{display:block;padding:0 20px 60px}}.toc{{display:none}}.stats{{grid-template-columns:repeat(2,1fr)}}h1{{font-size:40px}}.grid{{grid-template-columns:1fr}}}}
</style>
</head>
<body>
<header class="hero">
  <div class="hero-inner">
    <div class="kicker">Ready split · 2026-05-11</div>
    <h1>Paper / Validation /<br/><em>BD Split</em></h1>
    <p class="lead">
      This page separates the GO shortlist into three practical decision layers.
      The same records can appear in multiple layers when they support more than one use case.
    </p>
    <div class="stats">
      <div class="stat"><b>{len(rows)}</b><span>GO rows</span></div>
      <div class="stat"><b>{len(paperable)}</b><span>Paperable</span></div>
      <div class="stat"><b>{len(validation)}</b><span>Validation-ready</span></div>
      <div class="stat"><b>{len(bd)}</b><span>BD-ready</span></div>
      <div class="stat"><b>{sum(1 for r in rows if r['flags']==3)}</b><span>Triple-ready</span></div>
    </div>
  </div>
</header>
<div class="wrap">
    <nav class="toc">
      <h4>Contents</h4>
      <a href="#tldr">01 TL;DR</a>
      <a href="#paper">02 Paperable</a>
      <a href="#validation">03 Validation-ready</a>
      <a href="#bd">04 BD-ready</a>
      <a href="#all">05 Full tagging table</a>
    </nav>
<main>
<section id="tldr">
  <h2><span class="num">01</span>TL;DR</h2>
  <div class="sub">The same shortlist, split into decision layers</div>
    <div style="margin:12px 0 0;color:#9aabc0;font:12px 'JetBrains Mono',monospace;">
      <a href="index.html" style="color:#35d39d">Hub index</a> ·
      <a href="therapy_spectrum_recent_registry_2026_05_11.html" style="color:#7eb6ff">recent GEO registry</a> ·
      <a href="therapy_spectrum_go_shortlist_2026_05_11.html" style="color:#ffd28a">GO shortlist</a> ·
      <a href="therapy_spectrum_promising_master_2026_05_11.html" style="color:#c79cff">promising master</a> ·
      <a href="therapy_spectrum_triple_ready_2026_05_11.html" style="color:#35d39d">triple-ready</a> ·
      <a href="therapy_spectrum_bd_ready_2026_05_11.html" style="color:#7eb6ff">BD-ready</a>
    </div>
  <div class="grid">
    <div class="card"><h3>Paperable</h3><p>{len(paperable)} rows with a strong mechanistic or figure-bearing story.</p></div>
    <div class="card"><h3>Validation-ready</h3><p>{len(validation)} rows with a clean orthogonal / human-relevant confirmation path.</p></div>
    <div class="card"><h3>BD-ready</h3><p>{len(bd)} rows with direct translation, biomarker, or product-design value.</p></div>
  </div>
  <div class="box good"><b>Immediate emphasis:</b> CAR-T response durability, CAR design / manufacturing, mRNA-LNP platform chemistry, and human antibody / biomarker studies are the cleanest triple-ready lanes.</div>
</section>

<section id="paper">
  <h2><span class="num">02</span>Paperable</h2>
  <div class="sub">Mechanism + figure value</div>
  <table class="t">
    <thead><tr><th>Accession</th><th>Modality</th><th>Score</th><th>Tier</th><th>P</th><th>V</th><th>BD</th><th>Best task</th><th>Why</th></tr></thead>
    <tbody>{rows_table(paperable[:10])}</tbody>
  </table>
</section>

<section id="validation">
  <h2><span class="num">03</span>Validation-ready</h2>
  <div class="sub">Human-relevant or orthogonally testable</div>
  <table class="t">
    <thead><tr><th>Accession</th><th>Modality</th><th>Score</th><th>Tier</th><th>P</th><th>V</th><th>BD</th><th>Best task</th><th>Why</th></tr></thead>
    <tbody>{rows_table(validation[:10])}</tbody>
  </table>
</section>

<section id="bd">
  <h2><span class="num">04</span>BD-ready</h2>
  <div class="sub">Companion biomarker / product design / decision support</div>
  <table class="t">
    <thead><tr><th>Accession</th><th>Modality</th><th>Score</th><th>Tier</th><th>P</th><th>V</th><th>BD</th><th>Best task</th><th>Why</th></tr></thead>
    <tbody>{rows_table(bd[:10])}</tbody>
  </table>
  <div class="box warn"><b>Boundary.</b> Triple-ready means strong decision value, not proof of clinical success. Treat BD-ready as a priority signal for outreach, assay packaging, or companion biomarker work.</div>
</section>

<section id="all">
  <h2><span class="num">05</span>Full Tagging Table</h2>
  <div class="sub">Complete GO shortlist with decision flags</div>
  <table class="t">
    <thead><tr><th>Accession</th><th>Modality</th><th>Score</th><th>Tier</th><th>P</th><th>V</th><th>BD</th><th>Best task</th><th>Why</th></tr></thead>
    <tbody>{rows_table(rows)}</tbody>
  </table>
</section>
</main>
</div>
</body>
</html>"""


def main() -> None:
    ensure_dirs()
    rows = read_rows()
    tsv = RESULT / "therapy_spectrum_ready_split.tsv"
    summary = RESULT / "SUMMARY.md"
    html_path = HUB / "therapy_spectrum_ready_split_2026_05_11.html"

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
