"""Fix the Unpaywall PMC ID mismatches for Yoo2016 and Pu2021,
re-fetch the OA abstracts with the correct PMC IDs, and patch
both unpaywall_links.json and oa_abstracts.{json,md}.

Bug: v1 lit_enrich's Unpaywall fetch returned PMC4986964 for Yoo2016 (actually
Vedelek V "Testis-Specific Bb8") and PMC8523608 for Pu2021 (actually Dev SA
"Adulteration in Ayurvedic raw drugs"). Verified via Europe PMC DOI lookup:
correct PMCs are PMC4975456 (Yoo) and PMC8523550 (Pu).
"""

import json
import time
from pathlib import Path

import requests

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
DATA = ROOT / "project/manuscript_p2_brief/lit_enrich_2026_05_02/data"
MORE = ROOT / "project/manuscript_p2_brief/lit_enrich_2026_05_02/09_more_fetched"

EMAIL = "kukshomr@gmail.com"
HEADERS = {"User-Agent": f"thca-research/0.2 (mailto:{EMAIL})"}

CORRECTIONS = {
    "Yoo2016": {"old_pmc": "PMC4986964", "new_pmc": "PMC4975456"},
    "Pu2021": {"old_pmc": "PMC8523608", "new_pmc": "PMC8523550"},
}


def fetch_epmc(pmcid: str, doi: str | None = None) -> dict | None:
    pmc_norm = pmcid if pmcid.startswith("PMC") else f"PMC{pmcid}"
    url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/article/PMC/{pmc_norm}?resultType=core&format=json"
    r = requests.get(url, headers=HEADERS, timeout=30)
    if r.status_code != 200:
        return None
    j = r.json()
    result = j.get("result")
    if not result:
        return None
    authors = (result.get("authorList") or {}).get("author") or []
    return {
        "title": result.get("title"),
        "abstract": result.get("abstractText"),
        "journal": (result.get("journalInfo") or {}).get("journal", {}).get("title"),
        "year": result.get("pubYear"),
        "doi": result.get("doi"),
        "pmcid": result.get("pmcid"),
        "pmid": result.get("pmid"),
        "first_author": authors[0].get("fullName") if authors else None,
        "n_authors": len(authors) if isinstance(authors, list) else None,
    }


# 1. Fix unpaywall_links.json
print("=== Fixing unpaywall_links.json ===")
unpaywall = json.loads((DATA / "unpaywall_links.json").read_text())
for e in unpaywall:
    key = e.get("key")
    if key in CORRECTIONS:
        old = e.get("pmcid")
        new = CORRECTIONS[key]["new_pmc"]
        if old != CORRECTIONS[key]["old_pmc"]:
            print(f"  WARN: {key} pmcid is {old}, expected {CORRECTIONS[key]['old_pmc']} — skipping")
            continue
        e["pmcid"] = new
        e["europepmc_url"] = f"https://europepmc.org/article/PMC/{new}"
        e["_pmc_correction_note"] = (
            f"PMC ID corrected on 2026-05-02 from {old} (mismatched paper) to {new} (verified via DOI lookup)"
        )
        print(f"  ✓ {key}: {old} → {new}")

(DATA / "unpaywall_links.json").write_text(json.dumps(unpaywall, indent=2))

# 2. Re-fetch the corrected abstracts
print("\n=== Re-fetching corrected OA abstracts ===")
oa = json.loads((MORE / "oa_abstracts.json").read_text())
for key, fix in CORRECTIONS.items():
    print(f"  Fetching {key} with PMC {fix['new_pmc']}…")
    new_blob = fetch_epmc(fix["new_pmc"])
    if new_blob and new_blob.get("abstract"):
        first = new_blob.get("first_author") or ""
        expected_prefix = key.rstrip("0123456789")
        if expected_prefix.lower() in first.lower():
            print(f"    ✓ first_author = {first} (matches {expected_prefix})")
            oa[key] = new_blob
        else:
            print(f"    ✗ first_author = {first} STILL doesn't match {expected_prefix}!")
    else:
        print(f"    ✗ no abstract returned")
    time.sleep(1)

(MORE / "oa_abstracts.json").write_text(json.dumps(oa, indent=2))

# 3. Regenerate the oa_abstracts.md
print("\n=== Regenerating oa_abstracts.md ===")
got_abstract = sum(1 for v in oa.values() if isinstance(v, dict) and v.get("abstract"))
lines = [
    "# 09 OA abstracts (Europe PMC)",
    "",
    f"_{got_abstract}/{len(oa)} OA papers with abstract retrieved._",
    "",
    "_Note: Yoo2016 + Pu2021 PMC IDs corrected 2026-05-02 — original v1 Unpaywall fetch_",
    "_returned wrong PMCs (off-by-N collision: PMC4986964 was a Vedelek V testis paper,_",
    "_PMC8523608 was a Dev SA ayurvedic drugs paper). Corrections verified via DOI lookup._",
    "",
]
for key, blob in oa.items():
    if not isinstance(blob, dict):
        lines += [f"## {key}", "", "_(no result)_", "", "---", ""]
        continue
    lines += [
        f"## {key}",
        "",
        f"- Journal: {blob.get('journal') or '—'} ({blob.get('year') or '—'})",
        f"- Authors: {blob.get('first_author') or '—'} et al. (n={blob.get('n_authors') or '—'})",
        f"- DOI: `{blob.get('doi') or '—'}`",
        f"- PMCID: `{blob.get('pmcid') or '—'}` / PMID: `{blob.get('pmid') or '—'}`",
        "",
        f"**Title:** {blob.get('title') or '—'}",
        "",
        f"**Abstract:** {(blob.get('abstract') or '_(no abstract)_')[:2500]}",
        "",
        "---",
        "",
    ]
(MORE / "oa_abstracts.md").write_text("\n".join(lines))
print(f"  Wrote {MORE/'oa_abstracts.md'} ({(MORE/'oa_abstracts.md').stat().st_size/1024:.1f} KB)")

# 4. Re-run P4 fact-check to verify the fix
print("\n=== Re-running P4 fact-check ===")
import re
verified = json.loads((DATA / "verified_refs.json").read_text())
verified_by_key = {e.get("key"): e for e in verified if isinstance(e, dict)}
n_year, n_auth, n_both, total = 0, 0, 0, 0
for key, blob in oa.items():
    if not isinstance(blob, dict) or not blob.get("title"):
        continue
    total += 1
    expected_year = re.search(r"(\d{4})", key)
    year_match = bool(expected_year) and str(blob.get("year")) == expected_year.group(1)
    expected_author = re.match(r"([A-Z][a-z]+)", key)
    first = blob.get("first_author") or ""
    author_match = bool(expected_author) and expected_author.group(1).lower() in first.lower()
    n_year += int(year_match)
    n_auth += int(author_match)
    n_both += int(year_match and author_match)
print(f"  After fix: year {n_year}/{total}, author {n_auth}/{total}, both {n_both}/{total}")
