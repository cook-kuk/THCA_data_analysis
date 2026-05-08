#!/usr/bin/env python3
"""Track 8: AFND broad fetch — country x locus, healthy-only HLA frequencies.

Hard cap: 100 HTTP GETs total. Cache every fetched HTML to cache/ so re-runs are cheap.
Boundary: NO cancer association. Healthy populations only.
"""
from __future__ import annotations
import json
import os
import re
import sys
import time
import urllib.parse
from pathlib import Path

import requests
from bs4 import BeautifulSoup

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/hla_deepdive_2026_05_08/track8_afnd_extended")
CACHE = ROOT / "cache"
CACHE.mkdir(parents=True, exist_ok=True)

# Priority countries — AFND uses these exact "hla_country" tokens.
COUNTRIES = [
    # East Asia
    "South Korea", "Japan", "China", "Taiwan", "Mongolia",
    # Southeast Asia
    "Vietnam", "Thailand", "Cambodia", "Indonesia", "Malaysia", "Philippines", "Singapore",
    # South Asia
    "India", "Pakistan", "Sri Lanka", "Nepal",
    # Reference / outgroup
    "Germany", "United Kingdom", "USA",
]
# Note: at-risk skipped because of fetch budget: Laos, Bangladesh, North Korea
# (often 0 records in AFND, will note as "not in AFND").

# 6 main loci for HLA Class I + Class II core
LOCI = ["A", "B", "C", "DRB1", "DQB1", "DPB1"]

UA = "Mozilla/5.0 (X11; Linux x86_64) THCA-research-track8/1.0 (academic; healthy HLA atlas)"

BASE = "http://www.allelefrequencies.net/hla6006a.asp"
PARAMS_FIXED = {
    "hla_locus_type": "Classical",
    "hla_show": "<",  # show all
    "hla_order": "order_1",
    "standard": "",
}

REQUEST_BUDGET = 100
SLEEP_BETWEEN = 1.5  # seconds; AFND is slow but lets be polite


def cache_path(country: str, locus: str, page: int) -> Path:
    safe_country = country.replace(" ", "_")
    return CACHE / f"{safe_country}__{locus}__p{page}.html"


def fetch_one(country: str, locus: str, page: int, *, budget: list[int]) -> Path | None:
    cp = cache_path(country, locus, page)
    if cp.exists() and cp.stat().st_size > 1000:
        return cp
    if budget[0] <= 0:
        sys.stderr.write(f"BUDGET EXHAUSTED — skipping {country}/{locus}/p{page}\n")
        return None
    params = {
        **PARAMS_FIXED,
        "hla_locus": locus,
        "hla_country": country,
        "page": str(page),
    }
    url = BASE + "?" + urllib.parse.urlencode(params)
    sys.stderr.write(f"GET ({budget[0]} left) {country}/{locus}/p{page}\n")
    try:
        r = requests.get(url, headers={"User-Agent": UA}, timeout=45)
    except Exception as exc:
        sys.stderr.write(f"  ERR {exc}\n")
        budget[0] -= 1
        return None
    budget[0] -= 1
    time.sleep(SLEEP_BETWEEN)
    if r.status_code != 200:
        sys.stderr.write(f"  HTTP {r.status_code}\n")
        return None
    cp.write_text(r.text)
    return cp


def page_count(html_path: Path) -> int:
    txt = html_path.read_text(errors="ignore")
    # patterns like "of 3"
    m = re.search(r"\bof\s*(\d+)\b", txt)
    if m:
        return int(m.group(1))
    return 1


def main():
    budget = [REQUEST_BUDGET]
    plan: list[tuple[str, str, int]] = []
    for country in COUNTRIES:
        for locus in LOCI:
            plan.append((country, locus, 1))

    # Round 1 — page 1 for every (country, locus) pair (=114 if all 19 countries × 6 loci)
    # Trim to fit half the budget.
    phase1 = plan[: REQUEST_BUDGET // 2]
    phase1_remainder = plan[REQUEST_BUDGET // 2 :]

    fetched: list[tuple[str, str, int, Path]] = []
    for country, locus, _ in phase1:
        cp = fetch_one(country, locus, 1, budget=budget)
        if cp:
            fetched.append((country, locus, 1, cp))
            n = page_count(cp)
            # Schedule remaining pages (cap to 3 to avoid runaway)
            for p in range(2, min(n, 3) + 1):
                phase1_remainder.append((country, locus, p))
        if budget[0] <= 0:
            break

    # Phase 2 — process remainder until budget exhausted, prioritizing page-1 over page-2/3
    phase1_remainder.sort(key=lambda t: t[2])  # page asc → page-1 first
    for country, locus, page in phase1_remainder:
        if budget[0] <= 0:
            break
        cp = fetch_one(country, locus, page, budget=budget)
        if cp:
            fetched.append((country, locus, page, cp))

    manifest = {
        "n_fetched": len(fetched),
        "budget_used": REQUEST_BUDGET - budget[0],
        "budget_left": budget[0],
        "items": [
            {"country": c, "locus": l, "page": p, "path": str(fp)}
            for c, l, p, fp in fetched
        ],
    }
    (ROOT / "fetch_manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"Fetched {len(fetched)} pages, budget used = {manifest['budget_used']}")


if __name__ == "__main__":
    main()
