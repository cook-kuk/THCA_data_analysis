from __future__ import annotations

import csv
import html
import shutil
from collections import Counter
from pathlib import Path


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
SRC = ROOT / "project/results/therapy_spectrum_ready_split_2026_05_11/therapy_spectrum_ready_split.tsv"
RESULT = ROOT / "project/results/therapy_spectrum_promising_master_2026_05_11"
HUB = ROOT / "project/papers_hub_2026_05_04"
LIVE_HUB = Path("/var/www/papers/papers_hub_2026_05_04")
ASSET_DIR = HUB / "assets/therapy_spectrum_promising_master_2026_05_11"
LIVE_ASSET_DIR = LIVE_HUB / "assets/therapy_spectrum_promising_master_2026_05_11"


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
            row["flags"] = int(row["flags"])
            row["score_total"] = int(row["score_total"])
            row["paperability"] = int(row["paperability"])
            row["validation"] = int(row["validation"])
            row["bd"] = int(row["bd"])
            row["paperable"] = row["paperable"] == "True"
            row["validation_ready"] = row["validation_ready"] == "True"
            row["bd_ready"] = row["bd_ready"] == "True"
            rows.append(row)
    rows.sort(key=lambda r: (-r["flags"], -r["score_total"], r["accession"]))
    return rows


def classify(row: dict[str, object]) -> str:
    flags = int(row["flags"])
    if flags == 3:
        return "core"
    if flags == 2 and row["bd_ready"] and row["validation_ready"]:
        return "bd-core"
    if flags == 2:
        return "dual"
    if row["validation_ready"] and not row["bd_ready"]:
        return "validation"
    if row["paperable"] and not row["validation_ready"]:
        return "paper"
    if row["bd_ready"] and not row["validation_ready"]:
        return "bd"
    return "support"


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
        "class",
        "best_task",
        "why",
        "caveat",
        "source",
    ]
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=headers, delimiter="\t")
        w.writeheader()
        for row in rows:
            out = dict(row)
            out["class"] = classify(row)
            w.writerow({k: out.get(k, "") for k in headers})


def write_summary(path: Path, rows: list[dict[str, object]]) -> None:
    cls = Counter(classify(row) for row in rows)
    lines = [
        "# Therapy Spectrum Promising Master",
        "",
        "Date: 2026-05-11",
        "",
        "This master keeps every GO-shortlist record with at least one promising flag and separates them into practical execution classes.",
        "",
        "Snapshot:",
        f"- Rows: {len(rows)}",
        f"- Class counts: {', '.join(f'{k}={v}' for k, v in sorted(cls.items()))}",
        "",
        "Rule:",
        "- Core = all three flags on.",
        "- Dual = two flags on.",
        "- Validation / paper / BD = single-axis but still useful.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def table(rows: list[dict[str, object]]) -> str:
    return "\n".join(
        f"<tr>"
        f"<td><a href='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={esc(r['accession'])}'>{esc(r['accession'])}</a></td>"
        f"<td>{esc(r['modality'])}</td>"
        f"<td>{esc(r['score_total'])}</td>"
        f"<td><span class='tag {classify(r)}'>{esc(classify(r))}</span></td>"
        f"<td>{'Y' if r['paperable'] else 'N'}</td>"
        f"<td>{'Y' if r['validation_ready'] else 'N'}</td>"
        f"<td>{'Y' if r['bd_ready'] else 'N'}</td>"
        f"<td>{esc(r['best_task'])}</td>"
        f"<td>{esc(r['why'])}</td>"
        f"<td>{esc(r['caveat'])}</td>"
        f"<td><code>{esc(r['source'])}</code></td>"
        f"</tr>"
        for r in rows
    )


def build_html(rows: list[dict[str, object]]) -> str:
    core = [r for r in rows if classify(r) == "core"]
    bd_core = [r for r in rows if classify(r) == "bd-core"]
    dual = [r for r in rows if classify(r) == "dual"]
    validation = [r for r in rows if classify(r) == "validation"]
    paper = [r for r in rows if classify(r) == "paper"]
    bd = [r for r in rows if classify(r) == "bd"]
    support = [r for r in rows if classify(r) == "support"]
    top = rows[:6]
    cards = "".join(
        f"<div class='card'><h3>{i+1}. {esc(r['accession'])}</h3><p><b>{esc(r['modality'])}</b></p><p>{esc(classify(r))} | Score {esc(r['score_total'])}</p><p>{esc(r['best_task'])}</p><p class='small'>{esc(r['why'])}</p></div>"
        for i, r in enumerate(top)
    )
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8" />
<meta name="viewport" content="width=device-width,initial-scale=1.0" />
<title>Therapy Spectrum Promising Master · All Likely Winners</title>
<style>
:root{{--bg:#061019;--panel:#0f1726;--line:#22324b;--ink:#edf3fb;--muted:#9aabc0;--gold:#ffd28a;--teal:#35d39d;--blue:#7eb6ff;--red:#ff8a6b;--violet:#c79cff}}
*{{box-sizing:border-box}}
html,body{{margin:0;background:var(--bg);color:var(--ink);font:14px/1.55 Inter,"Noto Sans KR",sans-serif}}
a{{color:var(--teal);text-decoration:none}} a:hover{{color:var(--gold);text-decoration:underline}}
code{{font-family:"JetBrains Mono",monospace;font-size:12px}}
h1,h2,h3{{font-family:"Cormorant Garamond","Newsreader",serif;letter-spacing:-.02em;color:#fff7dc}}
.hero{{padding:64px 34px 38px;background:radial-gradient(circle at top left,#172946 0%,#0b1220 60%,#050810 100%);border-bottom:1px solid var(--line)}}
.hero-inner{{max-width:1440px;margin:0 auto}}
.kicker{{font:700 11px "JetBrains Mono",monospace;color:var(--gold);letter-spacing:.24em;text-transform:uppercase}}
h1{{font-size:62px;line-height:.96;margin:12px 0}}
.lead{{max-width:1080px;color:#cad5e5;font:18px/1.65 "Newsreader",serif}}
.stats{{display:grid;grid-template-columns:repeat(6,1fr);gap:10px;margin-top:24px}}
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
.tag.core{{background:rgba(53,211,157,.12);border-color:rgba(53,211,157,.28);color:var(--teal)}}
.tag.bd-core{{background:rgba(126,182,255,.12);border-color:rgba(126,182,255,.28);color:var(--blue)}}
.tag.dual{{background:rgba(199,156,255,.12);border-color:rgba(199,156,255,.28);color:var(--violet)}}
.tag.validation{{background:rgba(255,210,138,.12);border-color:rgba(255,210,138,.28);color:var(--gold)}}
.tag.paper{{background:rgba(255,138,107,.12);border-color:rgba(255,138,107,.28);color:var(--red)}}
.tag.bd{{background:rgba(126,182,255,.12);border-color:rgba(126,182,255,.28);color:var(--blue)}}
.tag.support{{background:rgba(148,165,186,.12);border-color:rgba(148,165,186,.24);color:var(--muted)}}
.box{{background:var(--panel);border:1px solid var(--line);border-left:4px solid var(--gold);border-radius:0 8px 8px 0;padding:14px 16px;margin:14px 0}}
.box.good{{border-left-color:var(--teal);background:#0a1a16}}
.box.warn{{border-left-color:var(--red);background:#1d1310}}
.small{{color:var(--muted);font-size:12px}}
@media(max-width:1100px){{.grid{{grid-template-columns:1fr 1fr}}.stats{{grid-template-columns:repeat(3,1fr)}}}}
@media(max-width:980px){{.stats{{grid-template-columns:repeat(2,1fr)}}.grid{{grid-template-columns:1fr}}h1{{font-size:40px}}.wrap{{padding:0 20px 60px}}}}
</style>
</head>
<body>
<header class="hero">
  <div class="hero-inner">
    <div class="kicker">Promising master · 2026-05-11</div>
    <h1>All Likely Winners<br/><em>in one screen</em></h1>
    <p class="lead">
      This page keeps every GO-shortlist record that still has at least one useful execution flag.
      It is the master list for anything that could realistically become a paper, validation package, or BD pitch.
    </p>
    <div style="margin-top:18px;color:#9aabc0;font:12px 'JetBrains Mono',monospace;">
      <a href="index.html" style="color:#35d39d">Hub index</a> ·
      <a href="therapy_spectrum_ready_split_2026_05_11.html" style="color:#f2c46d">ready split</a> ·
      <a href="therapy_spectrum_triple_ready_2026_05_11.html" style="color:#35d39d">triple-ready</a> ·
      <a href="therapy_spectrum_bd_ready_2026_05_11.html" style="color:#7eb6ff">BD-ready</a>
    </div>
    <div class="stats">
      <div class="stat"><b>{len(rows)}</b><span>Promising rows</span></div>
      <div class="stat"><b>{len(core)}</b><span>Core</span></div>
      <div class="stat"><b>{len(bd_core)}</b><span>BD-core</span></div>
      <div class="stat"><b>{len(dual)}</b><span>Dual</span></div>
      <div class="stat"><b>{len(validation)}</b><span>Validation</span></div>
      <div class="stat"><b>{len(paper)+len(bd)+len(support)}</b><span>Single/support</span></div>
    </div>
  </div>
</header>
<div class="wrap">
<section>
  <h2><span class="num">01</span>TL;DR</h2>
  <div class="sub">The shortlist of likely winners</div>
  <div class="grid">{cards}</div>
  <div class="box good"><b>Best overall path:</b> start with GSE290722, GSE303153, then branch to the strongest CAR-T engineering, biomarker, and antibody nodes.</div>
</section>

<section>
  <h2><span class="num">02</span>Core</h2>
  <div class="sub">All three flags on</div>
  <table class="t">
    <thead><tr><th>Accession</th><th>Modality</th><th>Score</th><th>Class</th><th>P</th><th>V</th><th>BD</th><th>Best task</th><th>Why</th></tr></thead>
    <tbody>{table(core)}</tbody>
  </table>
</section>

<section>
  <h2><span class="num">03</span>BD-core</h2>
  <div class="sub">Commercially direct but not triple-perfect</div>
  <table class="t">
    <thead><tr><th>Accession</th><th>Modality</th><th>Score</th><th>Class</th><th>P</th><th>V</th><th>BD</th><th>Best task</th><th>Why</th></tr></thead>
    <tbody>{table(bd_core)}</tbody>
  </table>
</section>

<section>
  <h2><span class="num">04</span>Dual</h2>
  <div class="sub">Strong but needs one more axis</div>
  <table class="t">
    <thead><tr><th>Accession</th><th>Modality</th><th>Score</th><th>Class</th><th>P</th><th>V</th><th>BD</th><th>Best task</th><th>Why</th></tr></thead>
    <tbody>{table(dual)}</tbody>
  </table>
</section>

<section>
  <h2><span class="num">05</span>Single axis</h2>
  <div class="sub">Still promising, but narrower</div>
  <table class="t">
    <thead><tr><th>Accession</th><th>Modality</th><th>Score</th><th>Class</th><th>P</th><th>V</th><th>BD</th><th>Best task</th><th>Why</th></tr></thead>
    <tbody>{table(validation + paper + bd + support)}</tbody>
  </table>
  <div class="box warn"><b>Boundary.</b> GSE291443 remains reserve because it does not clear any of the three execution flags. Keep it out of the likely-winner list unless the question changes.</div>
</section>
</div>
</body>
</html>"""


def main() -> None:
    ensure_dirs()
    rows = [r for r in load_rows() if int(r["flags"]) > 0]
    tsv = RESULT / "therapy_spectrum_promising_master.tsv"
    summary = RESULT / "SUMMARY.md"
    html_path = HUB / "therapy_spectrum_promising_master_2026_05_11.html"

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
