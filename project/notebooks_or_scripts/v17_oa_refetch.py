"""Re-fetch only OA abstracts (Europe PMC) with the PMC-prefix fix."""

import json
import time
from pathlib import Path

import requests

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
DATA = ROOT / "project/manuscript_p2_brief/lit_enrich_2026_05_02/data"
OUT = ROOT / "project/manuscript_p2_brief/lit_enrich_2026_05_02/09_more_fetched"

EMAIL = "kukshomr@gmail.com"
HEADERS = {"User-Agent": f"thca-research/0.2 (mailto:{EMAIL})"}


def fetch_epmc(pmcid: str | None, doi: str | None):
    if pmcid:
        pmc_norm = pmcid if pmcid.startswith("PMC") else f"PMC{pmcid}"
        url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/article/PMC/{pmc_norm}?resultType=core&format=json"
    elif doi:
        url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI:{doi}&resultType=core&format=json"
    else:
        return None
    try:
        r = requests.get(url, headers=HEADERS, timeout=30)
        if r.status_code != 200:
            return {"_status": r.status_code}
        j = r.json()
        # PMC core endpoint returns single 'result' dict
        if "result" in j and j["result"] and isinstance(j["result"], dict):
            result = j["result"]
        elif "resultList" in j:
            results = j.get("resultList", {}).get("result", [])
            result = results[0] if results else None
        else:
            result = None
        if not result:
            return {"_empty": True, "hitCount": j.get("hitCount")}
        # Extract authors compactly
        authors_obj = (result.get("authorList") or {}).get("author") or []
        first_author = (authors_obj[0].get("fullName") if authors_obj else None) if isinstance(authors_obj, list) else None
        return {
            "title": result.get("title"),
            "abstract": result.get("abstractText"),
            "journal": (result.get("journalInfo") or {}).get("journal", {}).get("title"),
            "year": result.get("pubYear"),
            "doi": result.get("doi"),
            "pmcid": result.get("pmcid"),
            "pmid": result.get("pmid"),
            "first_author": first_author,
            "n_authors": len(authors_obj) if isinstance(authors_obj, list) else None,
        }
    except Exception as e:
        return {"_error": str(e)}


unpaywall = json.loads((DATA / "unpaywall_links.json").read_text())
oa_entries = [e for e in unpaywall if e.get("is_oa") and (e.get("pmcid") or e.get("doi"))]

print(f"Re-fetching {len(oa_entries)} OA entries…")
oa_out = {}
for e in oa_entries:
    key = e["key"]
    res = fetch_epmc(e.get("pmcid"), e.get("doi"))
    has_abs = isinstance(res, dict) and bool(res.get("abstract"))
    oa_out[key] = res
    print(f"  {'✓' if has_abs else '✗'}  {key}  ({e.get('pmcid') or e.get('doi')})  abstract={'YES' if has_abs else 'no'}")
    time.sleep(1)

(OUT / "oa_abstracts.json").write_text(json.dumps(oa_out, indent=2))

got_abstract = sum(1 for v in oa_out.values() if isinstance(v, dict) and v.get("abstract"))
print(f"\nGot abstracts for {got_abstract}/{len(oa_out)} OA entries")

# Markdown digest
lines = ["# 09 OA abstracts (Europe PMC)", "", f"_{got_abstract}/{len(oa_out)} OA papers with abstract retrieved._", ""]
for key, blob in oa_out.items():
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
(OUT / "oa_abstracts.md").write_text("\n".join(lines))

print(f"Wrote {OUT/'oa_abstracts.json'} ({(OUT/'oa_abstracts.json').stat().st_size/1024:.1f} KB)")
print(f"Wrote {OUT/'oa_abstracts.md'} ({(OUT/'oa_abstracts.md').stat().st_size/1024:.1f} KB)")
