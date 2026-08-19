from __future__ import annotations

import csv
import html
import shutil
from collections import Counter
from pathlib import Path


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
SRC = ROOT / "project/results/therapy_spectrum_recent_registry_2026_05_11/therapy_spectrum_recent_registry.tsv"
RESULT = ROOT / "project/results/therapy_spectrum_go_shortlist_2026_05_11"
HUB = ROOT / "project/papers_hub_2026_05_04"
LIVE_HUB = Path("/var/www/papers/papers_hub_2026_05_04")
ASSET_DIR = HUB / "assets/therapy_spectrum_go_shortlist_2026_05_11"
LIVE_ASSET_DIR = LIVE_HUB / "assets/therapy_spectrum_go_shortlist_2026_05_11"


def esc(s: object) -> str:
    return html.escape(str(s), quote=True)


def score(row: dict[str, object]) -> int:
    return int(row["score_total"])


def family(modality: str) -> str:
    m = modality.lower()
    if "car" in m or "til" in m or "tcr" in m or "iNKT".lower() in m:
        return "cell therapy"
    if "bispecific" in m or "engager" in m or "blinatumomab" in m or "antibody" in m:
        return "antibody / engager"
    if "parp" in m or "kinase" in m or "hormone" in m or "epigenetic" in m:
        return "small molecule / epigenetic"
    if "mrna" in m or "nucleic" in m or "vaccine" in m:
        return "nucleic acid / vaccine"
    if "adc" in m:
        return "antibody-drug conjugate"
    if "virus" in m:
        return "virus"
    return "other"


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
            if row["status"] != "GO":
                continue
            row["rank"] = ""
            row["family"] = family(row["modality"])
            rows.append(row)
    rows.sort(key=lambda r: (-int(r["score_total"]), -int(r["paperability"]), r["accession"]))
    for i, row in enumerate(rows, start=1):
        row["rank"] = i
    return rows


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    headers = [
        "rank",
        "accession",
        "title",
        "date",
        "modality",
        "family",
        "status",
        "species",
        "n",
        "score_total",
        "best_task",
        "task_fit",
        "paperability",
        "validation",
        "bd",
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
    families = Counter(row["family"] for row in rows)
    top = rows[:5]
    lines = [
        "# Therapy Spectrum GO Shortlist",
        "",
        "Date: 2026-05-11",
        "",
        "This shortlist keeps only GO records from the recent therapy-spectrum registry and ranks them by paperability + validation + BD fit.",
        "",
        "Snapshot:",
        f"- Rows: {len(rows)}",
        f"- Families: {', '.join(f'{k}={v}' for k, v in families.most_common())}",
        f"- Top five: {top[0]['accession']} ({top[0]['modality']}); {top[1]['accession']} ({top[1]['modality']}); {top[2]['accession']} ({top[2]['modality']}); {top[3]['accession']} ({top[3]['modality']}); {top[4]['accession']} ({top[4]['modality']})",
        "",
        "Rule:",
        "- Keep human or human-relevant interventions with a direct readout.",
        "- Use this as the smallest list to push into analysis, validation planning, or BD screening.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_table(rows: list[dict[str, object]]) -> str:
    trs = []
    for row in rows:
        trs.append(
            "<tr>"
            f"<td>{esc(row['rank'])}</td>"
            f"<td><a href='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={esc(row['accession'])}'>{esc(row['accession'])}</a></td>"
            f"<td>{esc(row['title'])}</td>"
            f"<td>{esc(row['date'])}</td>"
            f"<td>{esc(row['modality'])}</td>"
            f"<td>{esc(row['family'])}</td>"
            f"<td>{esc(row['species'])}</td>"
            f"<td>{esc(row['n'])}</td>"
            f"<td>{esc(row['score_total'])}</td>"
            f"<td>{esc(row['best_task'])}</td>"
            f"<td>{esc(row['why'])}</td>"
            f"<td>{esc(row['caveat'])}</td>"
            f"<td><code>{esc(row['source'])}</code></td>"
            "</tr>"
        )
    return "\n".join(trs)


def build_html(rows: list[dict[str, object]]) -> str:
    families = Counter(row["family"] for row in rows)
    top = rows[:8]
    family_cards = "".join(
        f"<div class='stat'><b>{esc(count)}</b><span>{esc(label)}</span></div>"
        for label, count in families.most_common()
    )
    top_cards = "".join(
        f"<div class='card'><h3>{i+1}. {esc(row['accession'])}</h3><p><b>{esc(row['modality'])}</b></p><p>Rank {row['rank']} | Score {esc(row['score_total'])}</p><p>{esc(row['best_task'])}</p><p class='small'>{esc(row['why'])}</p></div>"
        for i, row in enumerate(top)
    )
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8" />
<meta name="viewport" content="width=device-width,initial-scale=1.0" />
<title>Therapy Spectrum GO Shortlist · GEO Action List</title>
<style>
:root{{--bg:#071019;--panel:#0f1726;--line:#22324b;--ink:#edf3fb;--muted:#9aabc0;--gold:#ffd28a;--teal:#35d39d;--red:#ff8a6b}}
*{{box-sizing:border-box}}
html,body{{margin:0;background:var(--bg);color:var(--ink);font:14px/1.55 Inter,"Noto Sans KR",sans-serif}}
a{{color:var(--teal);text-decoration:none}} a:hover{{color:var(--gold);text-decoration:underline}}
code{{font-family:"JetBrains Mono",monospace;font-size:12px}}
h1,h2,h3{{font-family:"Cormorant Garamond","Newsreader",serif;letter-spacing:-.02em;color:#fff7dc}}
.hero{{padding:64px 34px 40px;background:radial-gradient(circle at top left,#172946 0%,#0b1220 60%,#050810 100%);border-bottom:1px solid var(--line)}}
.hero-inner{{max-width:1380px;margin:0 auto}}
.kicker{{font:700 11px "JetBrains Mono",monospace;color:var(--gold);letter-spacing:.24em;text-transform:uppercase}}
h1{{font-size:60px;line-height:.98;margin:12px 0}}
.lead{{max-width:1040px;color:#ccdae9;font:18px/1.65 "Newsreader",serif}}
.stats{{display:grid;grid-template-columns:repeat(5,1fr);gap:10px;margin-top:24px}}
.stat{{background:rgba(255,255,255,.04);border:1px solid rgba(255,210,138,.18);padding:13px 14px;border-radius:10px}}
.stat b{{display:block;color:var(--gold);font:700 28px "Cormorant Garamond",serif;line-height:1}}
.stat span{{display:block;color:var(--muted);font:700 10px "JetBrains Mono",monospace;text-transform:uppercase;letter-spacing:.08em;margin-top:4px}}
.wrap{{max-width:1380px;margin:0 auto;display:grid;grid-template-columns:260px 1fr;gap:34px;padding:0 34px 80px}}
.toc{{position:sticky;top:0;align-self:start;max-height:100vh;overflow:auto;border-right:1px solid var(--line);padding:28px 18px 28px 0}}
.toc h4{{font:700 10px "JetBrains Mono",monospace;color:var(--gold);letter-spacing:.18em;text-transform:uppercase;margin:0 0 12px}}
.toc a{{display:block;color:#cbd5e5;padding:5px 0;font-size:12px}}
main{{min-width:0;padding-top:24px}}
section{{padding:32px 0;border-bottom:1px solid var(--line)}}
h2{{font-size:36px;margin:0 0 6px}}
h2 .num{{font:700 13px "JetBrains Mono",monospace;color:var(--gold);letter-spacing:.18em;margin-right:10px}}
.sub{{font:700 10px "JetBrains Mono",monospace;color:var(--muted);letter-spacing:.12em;text-transform:uppercase;margin-bottom:16px}}
.grid{{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:14px}}
.card{{background:var(--panel);border:1px solid var(--line);border-radius:10px;padding:14px 16px}}
.card h3{{margin:0 0 8px}}
.t{{width:100%;border-collapse:collapse;margin:12px 0 16px;font-size:12px}}
.t th,.t td{{border:1px solid var(--line);padding:8px 9px;vertical-align:top}}
.t th{{background:#16243a;color:var(--gold);font:700 10px "JetBrains Mono",monospace;letter-spacing:.05em;text-transform:uppercase;text-align:left}}
.t tr:nth-child(even) td{{background:rgba(255,255,255,.02)}}
.box{{background:var(--panel);border:1px solid var(--line);border-left:4px solid var(--gold);border-radius:0 8px 8px 0;padding:14px 16px;margin:14px 0}}
.box.good{{border-left-color:var(--teal);background:#0a1a16}}
.box.warn{{border-left-color:var(--red);background:#1d1310}}
.small{{color:var(--muted);font-size:12px}}
@media(max-width:1100px){{.grid{{grid-template-columns:1fr 1fr}}}}
@media(max-width:980px){{.wrap{{display:block;padding:0 20px 60px}}.toc{{display:none}}.stats{{grid-template-columns:repeat(2,1fr)}}h1{{font-size:40px}}.grid{{grid-template-columns:1fr}}}}
</style>
</head>
<body>
<header class="hero">
  <div class="hero-inner">
    <div class="kicker">GO-only shortlist · 2026-05-11</div>
    <h1>Therapy Spectrum<br/><em>GO Shortlist</em></h1>
    <p class="lead">
      This page keeps only the records that are already GO-tier for the platform decision tree:
      human or human-relevant, intervention plus response readout, and a plausible downstream decision artifact.
      It is the smallest list to push into analysis, validation, or BD screening.
    </p>
    <div class="crumbs" style="margin-top:18px;color:#94a5ba;font:12px 'JetBrains Mono',monospace;">
      <a href="index.html" style="color:#35d39d">← Hub index</a> ·
      <a href="therapy_spectrum_recent_registry_2026_05_11.html" style="color:#7eb6ff">recent GEO registry</a> ·
      <a href="therapy_spectrum_triage_2026_05_11.html" style="color:#ffd28a">therapy-spectrum triage</a> ·
      <a href="therapy_spectrum_ready_split_2026_05_11.html" style="color:#f2c46d">paper/validation/BD split</a>
    </div>
    <div class="stats">
      <div class="stat"><b>{len(rows)}</b><span>GO records</span></div>
      <div class="stat"><b>{len(families)}</b><span>Families</span></div>
      <div class="stat"><b>{rows[0]['accession']}</b><span>Top accession</span></div>
      <div class="stat"><b>{rows[0]['modality']}</b><span>Top modality</span></div>
      <div class="stat"><b>GEO</b><span>Latest public shortlist</span></div>
    </div>
  </div>
</header>
<div class="wrap">
<nav class="toc">
  <h4>Contents</h4>
  <a href="#tldr">01 TL;DR</a>
  <a href="#top">02 Top records</a>
  <a href="#table">03 Full table</a>
  <a href="#rules">04 Rules</a>
</nav>
<main>
<section id="tldr">
  <h2><span class="num">01</span>TL;DR</h2>
  <div class="sub">High-confidence records for immediate follow-up</div>
  <div class="stats">{family_cards}</div>
  <div class="box good"><b>Best lanes:</b> CAR-T / CAR-NKT / TIL engineering, mRNA-LNP platform chemistry, bispecific dose scheduling, and PARP / resistance-to-combo logic.</div>
</section>

<section id="top">
  <h2><span class="num">02</span>Top Records</h2>
  <div class="grid">{top_cards}</div>
</section>

<section id="table">
  <h2><span class="num">03</span>Full Table</h2>
  <div class="sub">Ranked by score_total from the expanded recent registry</div>
  <table class="t">
    <thead><tr><th>Rank</th><th>Accession</th><th>Title</th><th>Date</th><th>Modality</th><th>Family</th><th>Species</th><th>n</th><th>Score</th><th>Best task</th><th>Why</th><th>Caveat</th><th>Source</th></tr></thead>
    <tbody>{build_table(rows)}</tbody>
  </table>
</section>

<section id="rules">
  <h2><span class="num">04</span>Rules</h2>
  <div class="sub">How to use the shortlist</div>
  <pre>Keep:
 - human or human-relevant interventions
 - direct readouts, not just pathway speculation
 - records that can support a decision artifact

Do not overcall:
 - mouse-only model papers as clinical proof
 - toxicity/safety papers as efficacy evidence
 - platform engineering as disease-specific validation</pre>
  <div class="box warn"><b>Boundary.</b> GEO is a prioritization source, not a clinical verdict. Use the shortlist to choose what to analyze, validate, and pitch, not to infer patient benefit on its own.</div>
</section>
</main>
</div>
</body>
</html>"""


def main() -> None:
    ensure_dirs()
    rows = load_rows()
    tsv = RESULT / "therapy_spectrum_go_shortlist.tsv"
    summary = RESULT / "SUMMARY.md"
    html_path = HUB / "therapy_spectrum_go_shortlist_2026_05_11.html"

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
