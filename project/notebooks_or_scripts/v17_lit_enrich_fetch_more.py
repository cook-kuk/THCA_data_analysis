"""Two parallel fetch passes:

A. Sem Scholar retry with exponential backoff (6/8 claims got 0 hits before)
B. Fetch OA abstracts via Europe PMC (10 Unpaywall PMC IDs)

Output: 09_more_fetched/ — semscholar_retry.json, oa_abstracts.json
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import pandas as pd
import requests

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
DATA = ROOT / "project/manuscript_p2_brief/lit_enrich_2026_05_02/data"
DA = ROOT / "project/manuscript_p2_brief/lit_enrich_2026_05_02/deep_analysis"
OUT = ROOT / "project/manuscript_p2_brief/lit_enrich_2026_05_02/09_more_fetched"
OUT.mkdir(parents=True, exist_ok=True)

EMAIL = "kukshomr@gmail.com"
HEADERS = {"User-Agent": f"thca-research/0.2 (mailto:{EMAIL})"}


# ----------------------------------------------------------------------------- #
# A. Sem Scholar retry with backoff
# ----------------------------------------------------------------------------- #
CLAIMS = [
    ("HLA-II PTC autoimmunity (DPB1*05:01 etc.)", "HLA DPB1 papillary thyroid Hashimoto"),
    ("Hashimoto-like signature & PTC outcomes", "Hashimoto thyroiditis papillary thyroid transcriptome"),
    ("BRAF/RAS-negative ('dark matter') PTC subtypes", "BRAF wild-type RAS wild-type papillary thyroid molecular"),
    ("BCR repertoire / TLS in thyroid cancer", "B cell receptor tertiary lymphoid thyroid"),
    ("Single-cell PTC progression", "single-cell RNA papillary thyroid carcinoma progression"),
    ("Pan-Asian HLA fine-mapping (Graves' / autoimmune thyroid)", "Pan-Asian HLA Graves thyroid fine-mapping"),
    ("8-gene / driver-excluded PTC stratification", "8-gene signature papillary thyroid stratification"),
    ("TIERA-like Korean PTC molecular cohort", "Korean papillary thyroid molecular cohort transcriptome"),
]


def sem_scholar(query: str, limit: int = 10, max_attempts: int = 5):
    url = "https://api.semanticscholar.org/graph/v1/paper/search"
    params = {
        "query": query,
        "limit": limit,
        "fields": "title,year,authors,venue,citationCount,influentialCitationCount,externalIds,abstract",
    }
    delay = 2
    for attempt in range(max_attempts):
        try:
            r = requests.get(url, params=params, headers=HEADERS, timeout=30)
            if r.status_code == 200:
                return r.json().get("data", [])
            if r.status_code == 429:
                time.sleep(delay)
                delay *= 2
                continue
            return {"_error": r.status_code, "_text": r.text[:200]}
        except Exception as e:
            time.sleep(delay)
            delay *= 2
    return {"_error": "max_attempts"}


print("[A] Sem Scholar retry with backoff…", flush=True)
ss_out = {}
for label, query in CLAIMS:
    print(f"  - {label[:60]}", flush=True)
    res = sem_scholar(query)
    ss_out[label] = res
    time.sleep(3)  # polite spacing between distinct queries
(OUT / "semscholar_retry.json").write_text(json.dumps(ss_out, indent=2))
hits_per_claim = {k: (len(v) if isinstance(v, list) else 0) for k, v in ss_out.items()}
print(f"  Sem Scholar hits per claim: {hits_per_claim}", flush=True)


# ----------------------------------------------------------------------------- #
# B. Europe PMC abstract fetch for 10 OA Unpaywall entries
# ----------------------------------------------------------------------------- #
print("\n[B] Europe PMC abstract fetch…", flush=True)
unpaywall = json.loads((DATA / "unpaywall_links.json").read_text())
oa_entries = [e for e in unpaywall if e.get("is_oa") and (e.get("pmcid") or e.get("doi"))]


def fetch_epmc(pmcid: str | None, doi: str | None):
    if pmcid:
        # IMPORTANT: Europe PMC requires the "PMC" prefix in the path component
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
        # PMC core path
        result = (
            j.get("result")
            or (j.get("resultList", {}).get("result", [{}])[0] if j.get("resultList") else None)
        )
        if not result:
            return None
        return {
            "title": result.get("title"),
            "abstract": result.get("abstractText"),
            "journal": (result.get("journalInfo") or {}).get("journal", {}).get("title"),
            "year": result.get("pubYear"),
            "doi": result.get("doi"),
            "pmcid": result.get("pmcid"),
        }
    except Exception as e:
        return {"_error": str(e)}


oa_out = {}
for e in oa_entries:
    key = e["key"]
    print(f"  - {key} ({e.get('pmcid') or e.get('doi')})", flush=True)
    res = fetch_epmc(e.get("pmcid"), e.get("doi"))
    oa_out[key] = res
    time.sleep(1)
(OUT / "oa_abstracts.json").write_text(json.dumps(oa_out, indent=2))
got_abstract = sum(1 for v in oa_out.values() if isinstance(v, dict) and v.get("abstract"))
print(f"  Got abstracts for {got_abstract}/{len(oa_out)} OA entries")

# Quick markdown digest
lines = ["# 09 OA abstracts (Europe PMC)", "", f"_{got_abstract}/{len(oa_out)} OA papers with abstract retrieved._", ""]
for key, blob in oa_out.items():
    if not isinstance(blob, dict):
        continue
    lines += [
        f"## {key}",
        "",
        f"- Journal: {blob.get('journal') or '—'} ({blob.get('year') or '—'})",
        f"- DOI: `{blob.get('doi') or '—'}`",
        f"- PMCID: `{blob.get('pmcid') or '—'}`",
        "",
        f"**Title:** {blob.get('title') or '—'}",
        "",
        f"**Abstract:** {(blob.get('abstract') or '_(no abstract)_')[:1500]}",
        "",
        "---",
        "",
    ]
(OUT / "oa_abstracts.md").write_text("\n".join(lines))
print(f"\n[done] {OUT}")
for p in sorted(OUT.glob("*")):
    print(f"  {p.name} ({p.stat().st_size/1024:.1f} KB)")
