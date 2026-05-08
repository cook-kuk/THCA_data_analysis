#!/usr/bin/env python3
"""Parse cached AFND HTML pages → long-form TSV.

Columns: country, locus, allele, population, allele_freq, sample_size, source_page.
Then merges in the existing focus-allele dump (Track 1 source) for cross-check.
"""
from __future__ import annotations
import json
import re
from pathlib import Path

import pandas as pd
from bs4 import BeautifulSoup

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/hla_deepdive_2026_05_08/track8_afnd_extended")
CACHE = ROOT / "cache"
TABLES = ROOT / "tables"
TABLES.mkdir(parents=True, exist_ok=True)

EXISTING_AFND = Path("/home/seungho/personal/THCA_data_analysis/project/manuscript_p2_brief/lit_enrich_2026_05_02/data/afnd_alleles.json")


def parse_html(path: Path, country_hint: str, locus_hint: str) -> list[dict]:
    txt = path.read_text(errors="ignore")
    soup = BeautifulSoup(txt, "html.parser")
    tables = soup.find_all("table")
    target = None
    for t in tables:
        cls = t.get("class") or []
        if "tblNormal" in cls:
            target = t
            break
    if target is None:
        return []
    rows = target.find_all("tr")
    out = []
    for r in rows[1:]:  # skip header
        cells = r.find_all(["td", "th"])
        if len(cells) < 8:
            continue
        # Real body row has 12 cells (header has 10 due to colspan):
        # [0] line, [1] allele, [2] '', [3] population, [4] %indiv, [5] allele_freq,
        # [6] '', [7] sample_size, [8] 'See' link, [9..] empty.
        line = cells[0].get_text(strip=True)
        if not re.match(r"^\d+$", line):
            continue
        allele = cells[1].get_text(strip=True)
        pop = cells[3].get_text(strip=True)
        pct_indiv = cells[4].get_text(strip=True).replace(",", "")
        afreq = cells[5].get_text(strip=True).replace(",", "")
        ssize = cells[7].get_text(strip=True).replace(",", "")
        if not allele:
            continue
        try:
            af = float(afreq) if afreq else float("nan")
        except ValueError:
            af = float("nan")
        try:
            n = int(ssize) if ssize else 0
        except ValueError:
            n = 0
        out.append({
            "country": country_hint,
            "locus": locus_hint,
            "allele": allele,
            "population": pop,
            "allele_freq": af,
            "sample_size": n,
            "pct_individuals": pct_indiv,
            "source_file": path.name,
        })
    return out


def main():
    rows: list[dict] = []
    files = sorted(CACHE.glob("*.html"))
    print(f"Parsing {len(files)} cached pages")
    for f in files:
        m = re.match(r"^(.+)__([A-Z]+\d?[A-Z]*)__p(\d+)\.html$", f.name)
        if not m:
            continue
        country = m.group(1).replace("_", " ")
        locus = m.group(2)
        rows.extend(parse_html(f, country, locus))
    df = pd.DataFrame(rows)
    print("parsed rows:", len(df), "unique pops:", df["population"].nunique() if len(df) else 0)
    if len(df):
        # Standardize allele to 2-field resolution where possible.
        def two_field(a: str) -> str:
            # Drop 3rd/4th fields; "A*02:01:01:01" → "A*02:01"
            parts = a.split(":")
            if len(parts) >= 2:
                return ":".join(parts[:2])
            return a
        df["allele_2f"] = df["allele"].apply(two_field)
        df = df[df["allele"].str.contains(":")]  # keep only 2+ field alleles
        # Map populations to broad regions
        df["region"] = df["country"].map({
            "South Korea": "East Asia",
            "Japan": "East Asia",
            "China": "East Asia",
            "Taiwan": "East Asia",
            "Mongolia": "East Asia",
            "Vietnam": "Southeast Asia",
            "Thailand": "Southeast Asia",
            "Cambodia": "Southeast Asia",
            "Indonesia": "Southeast Asia",
            "Malaysia": "Southeast Asia",
            "Philippines": "Southeast Asia",
            "Singapore": "Southeast Asia",
            "India": "South Asia",
            "Pakistan": "South Asia",
            "Sri Lanka": "South Asia",
            "Nepal": "South Asia",
            "Germany": "Outgroup_Europe",
            "United Kingdom": "Outgroup_Europe",
            "USA": "Outgroup_Mixed",
        }).fillna("Other")
    df.to_csv(TABLES / "afnd_long.tsv", sep="\t", index=False)
    print(f"wrote {TABLES/'afnd_long.tsv'}")

    # Per-country, per-locus summary
    if len(df):
        summary = (
            df.groupby(["country", "locus"], as_index=False)
              .agg(n_alleles=("allele_2f", "nunique"),
                   n_pops=("population", "nunique"),
                   total_n=("sample_size", "max"))
        )
        summary.to_csv(TABLES / "fetch_summary.tsv", sep="\t", index=False)
        print("\nFetch summary (country × locus):")
        print(summary.to_string(index=False))

    # Pull the focus alleles from existing dump for cross-check
    if EXISTING_AFND.exists():
        existing = json.loads(EXISTING_AFND.read_text())
        focus_rows = []
        for allele, by_country in existing.items():
            for ctry, items in by_country.items():
                for it in items:
                    try:
                        af = float(it["allele_freq"])
                    except (KeyError, ValueError):
                        af = float("nan")
                    try:
                        ss = int(str(it["sample_size"]).replace(",", ""))
                    except (KeyError, ValueError):
                        ss = 0
                    focus_rows.append({
                        "country": ctry,
                        "allele": allele,
                        "population": it.get("population"),
                        "allele_freq": af,
                        "sample_size": ss,
                    })
        fdf = pd.DataFrame(focus_rows)
        fdf.to_csv(TABLES / "track1_focus_alleles.tsv", sep="\t", index=False)
        print(f"existing focus alleles parsed: {len(fdf)} rows, {fdf['allele'].nunique()} alleles")


if __name__ == "__main__":
    main()
