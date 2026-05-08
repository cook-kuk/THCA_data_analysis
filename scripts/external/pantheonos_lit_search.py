"""
PantheonOS paper_reporter_v2 → repurposed literature search (no PantheonOS runtime).

Replaces:
  pantheon.smart_func + pantheon.toolsets.web_browse  →  requests + crossref/biorxiv API

Usage:
  python scripts/external/pantheonos_lit_search.py \\
      --theme "DM1 thyroid molecular dark matter" \\
      --out project/reports/lit_search/dm1_dark_matter.md \\
      --max 30

Per the marathon override (v19_marathon_override_pantheonos_2026_05_08), this is
infra. It feeds Discussion related-work, NOT voice-protected sections.

Source agent prompt (verbatim adapted) lives in:
  project/external/pantheonos/paper_reporter_v2/paper_reporter_v2.py
"""
import argparse
import json
import sys
import time
from pathlib import Path
from urllib.parse import quote

import requests


CROSSREF = "https://api.crossref.org/works"
BIORXIV = "https://api.biorxiv.org/details/biorxiv"


def crossref_search(query: str, rows: int = 30) -> list[dict]:
    r = requests.get(
        CROSSREF,
        params={"query": query, "rows": rows, "select": "DOI,title,author,issued,container-title,abstract,URL"},
        headers={"User-Agent": "thca-lit-search/1.0 (mailto:kukshomr@gmail.com)"},
        timeout=30,
    )
    r.raise_for_status()
    items = r.json()["message"]["items"]
    return [
        {
            "title": (it.get("title") or [""])[0],
            "doi": it.get("DOI", ""),
            "url": it.get("URL", ""),
            "journal": (it.get("container-title") or [""])[0],
            "year": (it.get("issued", {}).get("date-parts") or [[None]])[0][0],
            "authors": ", ".join(
                f"{a.get('family','')} {a.get('given','')}".strip()
                for a in (it.get("author") or [])[:5]
            ),
            "abstract": it.get("abstract", "").replace("<jats:p>", "").replace("</jats:p>", ""),
        }
        for it in items
    ]


def write_report(theme: str, hits: list[dict], out: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"# Literature search: {theme}", "", f"Found: **{len(hits)}** papers via Crossref. Generated {time.strftime('%Y-%m-%d')}.", ""]
    for i, h in enumerate(hits, 1):
        lines += [
            f"## {i}. {h['title']}",
            f"- **Authors:** {h['authors']}",
            f"- **Journal:** {h['journal']} ({h['year']})",
            f"- **DOI:** [{h['doi']}]({h['url']})",
            "",
            f"{h['abstract'][:600]}{'...' if len(h['abstract']) > 600 else ''}",
            "",
        ]
    out.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {out} ({len(hits)} hits)")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--theme", required=True, help="Search theme/query string")
    p.add_argument("--out", required=True, help="Output markdown path")
    p.add_argument("--max", type=int, default=30)
    args = p.parse_args()

    print(f"[crossref] querying: {args.theme!r} (rows={args.max})", file=sys.stderr)
    hits = crossref_search(args.theme, rows=args.max)
    write_report(args.theme, hits, Path(args.out))


if __name__ == "__main__":
    main()
