from __future__ import annotations

import csv
import html
import shutil
from pathlib import Path


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
SRC = ROOT / "project/results/therapy_spectrum_ready_split_2026_05_11/therapy_spectrum_ready_split.tsv"
RESULT = ROOT / "project/results/therapy_spectrum_priority_slices_2026_05_11"
HUB = ROOT / "project/papers_hub_2026_05_04"
LIVE_HUB = Path("/var/www/papers/papers_hub_2026_05_04")
ASSET_DIR = HUB / "assets/therapy_spectrum_priority_slices_2026_05_11"
LIVE_ASSET_DIR = LIVE_HUB / "assets/therapy_spectrum_priority_slices_2026_05_11"


def esc(s: object) -> str:
    return html.escape(str(s), quote=True)


def ensure_dirs() -> None:
    RESULT.mkdir(parents=True, exist_ok=True)
    ASSET_DIR.mkdir(parents=True, exist_ok=True)
    LIVE_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    HUB.mkdir(parents=True, exist_ok=True)
    LIVE_HUB.mkdir(parents=True, exist_ok=True)


def load_rows() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with SRC.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(row)
    rows.sort(key=lambda r: (-int(r["score_total"]), -int(r["flags"]), r["accession"]))
    return rows


def filter_rows(rows: list[dict[str, object]], mode: str) -> list[dict[str, object]]:
    if mode == "triple":
        return [r for r in rows if int(r["flags"]) == 3]
    if mode == "bd":
        return [r for r in rows if r["bd_ready"] == "True"]
    raise ValueError(mode)


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


def write_summary(path: Path, rows: list[dict[str, object]], mode: str) -> None:
    title = "Triple-ready" if mode == "triple" else "BD-ready"
    lines = [
        f"# Therapy Spectrum {title} Slice",
        "",
        "Date: 2026-05-11",
        "",
        f"This slice keeps only the GO shortlist records tagged as {title.lower()}.",
        "",
        "Snapshot:",
        f"- Rows: {len(rows)}",
        f"- Top accession: {rows[0]['accession'] if rows else 'NA'}",
        "",
        "Rule:",
        "- Use this as the shortest operational list for the next action stage.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_table(rows: list[dict[str, object]]) -> str:
    return "\n".join(
        f"<tr>"
        f"<td><a href='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={esc(r['accession'])}'>{esc(r['accession'])}</a></td>"
        f"<td>{esc(r['modality'])}</td>"
        f"<td>{esc(r['score_total'])}</td>"
        f"<td>{esc(r['best_task'])}</td>"
        f"<td>{esc(r['why'])}</td>"
        f"<td>{esc(r['caveat'])}</td>"
        f"<td><code>{esc(r['source'])}</code></td>"
        f"</tr>"
        for r in rows
    )


def build_html(rows: list[dict[str, object]], mode: str) -> str:
    title = "Triple-ready" if mode == "triple" else "BD-ready"
    kicker = "Operational slice"
    lead = (
        "The shortest list with all three flags on: paperable, validation-ready, and BD-ready."
        if mode == "triple"
        else "The shortest list with a direct business-development signal."
    )
    top = rows[:6]
    cards = "".join(
        f"<div class='card'><h3>{i+1}. {esc(r['accession'])}</h3><p><b>{esc(r['modality'])}</b></p><p>Score {esc(r['score_total'])}</p><p>{esc(r['best_task'])}</p></div>"
        for i, r in enumerate(top)
    )
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8" />
<meta name="viewport" content="width=device-width,initial-scale=1.0" />
<title>Therapy Spectrum {title} Slice</title>
<style>
:root{{--bg:#061019;--panel:#0f1726;--line:#22324b;--ink:#edf3fb;--muted:#9aabc0;--gold:#ffd28a;--teal:#35d39d;--blue:#7eb6ff;--red:#ff8a6b}}
*{{box-sizing:border-box}}
html,body{{margin:0;background:var(--bg);color:var(--ink);font:14px/1.55 Inter,"Noto Sans KR",sans-serif}}
a{{color:var(--teal);text-decoration:none}} a:hover{{color:var(--gold);text-decoration:underline}}
h1,h2,h3{{font-family:"Cormorant Garamond","Newsreader",serif;letter-spacing:-.02em;color:#fff7dc}}
code{{font-family:"JetBrains Mono",monospace;font-size:12px}}
.hero{{padding:64px 34px 38px;background:radial-gradient(circle at top left,#172946 0%,#0b1220 60%,#050810 100%);border-bottom:1px solid var(--line)}}
.hero-inner{{max-width:1280px;margin:0 auto}}
.kicker{{font:700 11px "JetBrains Mono",monospace;color:var(--gold);letter-spacing:.24em;text-transform:uppercase}}
h1{{font-size:60px;line-height:.98;margin:12px 0}}
.lead{{max-width:980px;color:#cad5e5;font:18px/1.65 "Newsreader",serif}}
.stats{{display:grid;grid-template-columns:repeat(4,1fr);gap:10px;margin-top:24px}}
.stat{{background:rgba(255,255,255,.04);border:1px solid rgba(255,210,138,.18);padding:13px 14px;border-radius:10px}}
.stat b{{display:block;color:var(--gold);font:700 28px "Cormorant Garamond",serif;line-height:1}}
.stat span{{display:block;color:var(--muted);font:700 10px "JetBrains Mono",monospace;text-transform:uppercase;letter-spacing:.08em;margin-top:4px}}
.wrap{{max-width:1280px;margin:0 auto;padding:0 34px 80px}}
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
.box{{background:var(--panel);border:1px solid var(--line);border-left:4px solid var(--gold);border-radius:0 8px 8px 0;padding:14px 16px;margin:14px 0}}
.box.good{{border-left-color:var(--teal);background:#0a1a16}}
.box.warn{{border-left-color:var(--red);background:#1d1310}}
@media(max-width:980px){{.stats{{grid-template-columns:repeat(2,1fr)}}.grid{{grid-template-columns:1fr}}h1{{font-size:40px}}.wrap{{padding:0 20px 60px}}}}
</style>
</head>
<body>
<header class="hero">
  <div class="hero-inner">
    <div class="kicker">{kicker} · 2026-05-11</div>
    <h1>Therapy Spectrum<br/><em>{title}</em></h1>
    <p class="lead">{lead}</p>
    <div style="margin-top:18px;color:#9aabc0;font:12px 'JetBrains Mono',monospace;">
      <a href="index.html" style="color:#35d39d">Hub index</a> ·
      <a href="therapy_spectrum_ready_split_2026_05_11.html" style="color:#f2c46d">ready split</a> ·
      <a href="therapy_spectrum_recent_registry_2026_05_11.html" style="color:#7eb6ff">recent registry</a> ·
      <a href="therapy_spectrum_go_shortlist_2026_05_11.html" style="color:#ffd28a">GO shortlist</a> ·
      <a href="therapy_spectrum_promising_master_2026_05_11.html" style="color:#c79cff">promising master</a>
    </div>
    <div class="stats">
      <div class="stat"><b>{len(rows)}</b><span>Rows</span></div>
      <div class="stat"><b>{rows[0]['accession'] if rows else 'NA'}</b><span>Top accession</span></div>
      <div class="stat"><b>{rows[0]['modality'] if rows else 'NA'}</b><span>Top modality</span></div>
      <div class="stat"><b>GEO</b><span>Operational slice</span></div>
    </div>
  </div>
</header>
<div class="wrap">
<section>
  <h2><span class="num">01</span>TL;DR</h2>
  <div class="sub">Immediate prioritization list</div>
  <div class="grid">{cards}</div>
  <div class="box good"><b>Use case:</b> if you only have time for one next action stage, start here.</div>
</section>
<section>
  <h2><span class="num">02</span>Full Table</h2>
  <div class="sub">Shortest operational rows</div>
  <table class="t">
    <thead><tr><th>Accession</th><th>Modality</th><th>Score</th><th>Best task</th><th>Why</th><th>Caveat</th><th>Source</th></tr></thead>
    <tbody>{build_table(rows)}</tbody>
  </table>
</section>
<section>
  <h2><span class="num">03</span>Boundary</h2>
  <div class="box warn">This slice is for prioritization only. It does not upgrade GEO into clinical evidence by itself.</div>
</section>
</div>
</body>
</html>"""


def main() -> None:
    ensure_dirs()
    rows = load_rows()
    triple = filter_rows(rows, "triple")
    bd = filter_rows(rows, "bd")

    triple_tsv = RESULT / "therapy_spectrum_triple_ready.tsv"
    bd_tsv = RESULT / "therapy_spectrum_bd_ready.tsv"
    triple_summary = RESULT / "TRIPLE_SUMMARY.md"
    bd_summary = RESULT / "BD_SUMMARY.md"
    triple_html = HUB / "therapy_spectrum_triple_ready_2026_05_11.html"
    bd_html = HUB / "therapy_spectrum_bd_ready_2026_05_11.html"

    write_tsv(triple_tsv, triple)
    write_tsv(bd_tsv, bd)
    write_summary(triple_summary, triple, "triple")
    write_summary(bd_summary, bd, "bd")
    triple_html.write_text(build_html(triple, "triple"), encoding="utf-8")
    bd_html.write_text(build_html(bd, "bd"), encoding="utf-8")

    for src in [triple_tsv, bd_tsv, triple_summary, bd_summary]:
        shutil.copy2(src, ASSET_DIR / src.name)
        shutil.copy2(src, LIVE_ASSET_DIR / src.name)
    shutil.copy2(triple_html, LIVE_HUB / triple_html.name)
    shutil.copy2(bd_html, LIVE_HUB / bd_html.name)


if __name__ == "__main__":
    main()
