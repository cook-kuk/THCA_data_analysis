#!/usr/bin/env python3
"""
Track 10 — Step 1: PMID / DOI resolution for 5 Korean AITD HLA papers.

Boundary: AUTOIMMUNE-ONLY (Graves' / Hashimoto / AITD); zero cancer claim.

Targets:
    Shin 2019  PLoS ONE  -- Korean pediatric AITD HLA-I/II + amino acids
    Cho  2011             -- Korean children AITD
    Park 2005             -- Korean GD HLA-DR/DQ
    Jang 2011             -- Korean GD HLA-DRB1 SBT
    Baek 2021             -- Korean class-II NGS reference

Strategy: PubMed esearch + esummary + Crossref works lookup, store top candidate
metadata + keep top-3 candidates as JSON.
"""
from __future__ import annotations
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path
import pandas as pd

OUT = Path("/home/seungho/personal/THCA_data_analysis/project/results/hla_deepdive_2026_05_08/track10_korean_lit")
TBL = OUT / "tables"
TBL.mkdir(parents=True, exist_ok=True)

NCBI = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
CROSSREF = "https://api.crossref.org/works"

# Search queries per paper (run several phrasings; pick best match)
QUERIES = {
    "Shin_2019": [
        "Shin 2019 Korean autoimmune thyroid HLA PLoS ONE",
        "Shin Korean pediatric autoimmune thyroid disease HLA",
        "HLA Korean children autoimmune thyroid 2019",
        "HLA class II amino acid Korean Graves Hashimoto pediatric",
    ],
    "Cho_2011": [
        "Cho Korean children autoimmune thyroid HLA 2011",
        "Cho Korean pediatric Graves Hashimoto HLA 2011",
        "HLA Korean children Graves' disease 2011",
    ],
    "Park_2005": [
        "Park Korean Graves disease HLA-DR DQ 2005",
        "Park Korean Graves' disease HLA DRB1 DQB1 2005",
        "HLA Korean Graves DPB1 2005 Park",
    ],
    "Jang_2011": [
        "Jang Korean Graves disease HLA-DRB1 SBT 2011",
        "HLA-DRB1 sequence-based typing Korean Graves 2011",
        "DRB1 Korean Graves thyroid 2011 Jang",
    ],
    "Baek_2021": [
        "Baek Korean HLA class II next-generation sequencing 2021",
        "HLA-DPB1 DQB1 DRB1 Korean reference NGS 2021 Baek",
        "Korean HLA-DPA1 DPB1 high-resolution NGS 2021",
    ],
}


def http_json(url: str, retries: int = 3) -> dict | None:
    for k in range(retries):
        try:
            req = urllib.request.Request(
                url,
                headers={"User-Agent": "track10-korean-lit-resolver/1.0 (mailto:kukshomr@gmail.com)"},
            )
            with urllib.request.urlopen(req, timeout=30) as r:
                return json.loads(r.read().decode("utf-8"))
        except Exception as e:
            if k == retries - 1:
                print(f"   ! HTTP fail: {url[:80]}... -> {e}")
                return None
            time.sleep(1.5 * (k + 1))


def pubmed_search(term: str, retmax: int = 5) -> list[str]:
    url = f"{NCBI}/esearch.fcgi?db=pubmed&retmode=json&retmax={retmax}&term={urllib.parse.quote(term)}"
    data = http_json(url)
    if not data:
        return []
    return data.get("esearchresult", {}).get("idlist", [])


def pubmed_summary(pmids: list[str]) -> list[dict]:
    if not pmids:
        return []
    url = f"{NCBI}/esummary.fcgi?db=pubmed&retmode=json&id={','.join(pmids)}"
    data = http_json(url)
    if not data:
        return []
    out = []
    res = data.get("result", {})
    for pmid in pmids:
        rec = res.get(pmid, {})
        if not rec:
            continue
        # Extract DOI from articleids
        doi = ""
        for aid in rec.get("articleids", []):
            if aid.get("idtype") == "doi":
                doi = aid.get("value", "")
                break
        out.append(dict(
            pmid=pmid,
            title=rec.get("title", ""),
            journal=rec.get("fulljournalname", "") or rec.get("source", ""),
            pubdate=rec.get("pubdate", ""),
            volume=rec.get("volume", ""),
            issue=rec.get("issue", ""),
            pages=rec.get("pages", ""),
            authors=";".join(a.get("name", "") for a in rec.get("authors", [])),
            doi=doi,
        ))
    return out


def crossref_search(q: str, rows: int = 5) -> list[dict]:
    url = f"{CROSSREF}?query={urllib.parse.quote(q)}&rows={rows}"
    data = http_json(url)
    if not data:
        return []
    items = data.get("message", {}).get("items", [])
    out = []
    for it in items:
        title = (it.get("title") or [""])[0]
        authors = ";".join(
            f"{a.get('family','')} {a.get('given','')}".strip()
            for a in it.get("author", [])[:6]
        )
        cont = it.get("container-title", [""])
        out.append(dict(
            doi=it.get("DOI", ""),
            title=title,
            journal=cont[0] if cont else "",
            year=(it.get("issued", {}).get("date-parts", [[None]])[0][0] or ""),
            volume=it.get("volume", ""),
            issue=it.get("issue", ""),
            page=it.get("page", ""),
            authors=authors,
            url=it.get("URL", ""),
        ))
    return out


# ----- Run -----
all_resolved = {}
records = []

for paper_id, queries in QUERIES.items():
    print(f"\n[{paper_id}] resolving ...")
    pubmed_hits = []
    seen_pmids = set()
    for q in queries:
        ids = pubmed_search(q, retmax=5)
        for pmid in ids:
            if pmid in seen_pmids:
                continue
            seen_pmids.add(pmid)
        time.sleep(0.4)
    pubmed_summary_list = pubmed_summary(list(seen_pmids))
    crossref_hits = []
    for q in queries:
        crossref_hits.extend(crossref_search(q, rows=5))
        time.sleep(0.4)
    # Dedup crossref by DOI
    seen_dois = set()
    crossref_unique = []
    for h in crossref_hits:
        if h["doi"] in seen_dois:
            continue
        seen_dois.add(h["doi"])
        crossref_unique.append(h)
    all_resolved[paper_id] = dict(
        pubmed=pubmed_summary_list,
        crossref=crossref_unique[:10],
    )
    print(f"   PubMed candidates: {len(pubmed_summary_list)}, Crossref candidates: {len(crossref_unique)}")
    if pubmed_summary_list:
        print(f"   top PubMed: PMID={pubmed_summary_list[0]['pmid']}  {pubmed_summary_list[0]['title'][:80]}")
    if crossref_unique:
        print(f"   top Crossref: DOI={crossref_unique[0]['doi']}  {crossref_unique[0]['title'][:80]}")

with open(OUT / "pmid_resolution_candidates.json", "w") as f:
    json.dump(all_resolved, f, indent=2)
print(f"\n[saved] {OUT/'pmid_resolution_candidates.json'}")

# Build a flat table for human review
rows = []
for paper_id, r in all_resolved.items():
    for src, hits in (("pubmed", r["pubmed"]), ("crossref", r["crossref"])):
        for h in hits:
            rows.append(dict(
                paper_id=paper_id,
                source=src,
                pmid=h.get("pmid", ""),
                doi=h.get("doi", ""),
                title=h.get("title", ""),
                journal=h.get("journal", ""),
                year=h.get("pubdate", h.get("year", "")),
                volume=h.get("volume", ""),
                issue=h.get("issue", ""),
                pages=h.get("pages", h.get("page", "")),
                authors=h.get("authors", ""),
                url=h.get("url", ""),
            ))
df = pd.DataFrame(rows)
df.to_csv(TBL / "T01_pmid_doi_candidates.tsv", sep="\t", index=False)
print(f"[saved] {TBL/'T01_pmid_doi_candidates.tsv'} ({len(df)} candidate rows)")
