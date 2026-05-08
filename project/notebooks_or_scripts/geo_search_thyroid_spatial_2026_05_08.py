#!/usr/bin/env python3
"""GEO search for thyroid spatial transcriptomics datasets.
Output: project/results/geo_search_2026_05_08/
Strategy: ESearch (gds db) → ESummary → filter for spatial / Visium / 10x.
"""
from __future__ import annotations
from pathlib import Path
import time, json, sys, re
import urllib.parse, urllib.request
import xml.etree.ElementTree as ET

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT/"project/results/geo_search_2026_05_08"
OUT.mkdir(exist_ok=True, parents=True)

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ALREADY_HAVE = {"GSE250521", "GSE230424", "GSE248205"}

def esearch(db, term, retmax=200):
    url = f"{EUTILS}/esearch.fcgi?db={db}&term={urllib.parse.quote(term)}&retmax={retmax}&usehistory=y"
    with urllib.request.urlopen(url, timeout=20) as r:
        text = r.read().decode()
    root = ET.fromstring(text)
    ids = [e.text for e in root.findall(".//IdList/Id")]
    count = root.findtext("Count", "0")
    return ids, int(count)

def esummary(db, ids, chunk=80):
    out = []
    for i in range(0, len(ids), chunk):
        batch = ids[i:i+chunk]
        url = f"{EUTILS}/esummary.fcgi?db={db}&id={','.join(batch)}&retmode=json"
        with urllib.request.urlopen(url, timeout=30) as r:
            data = json.loads(r.read().decode())
        for k in data.get("result", {}).get("uids", []):
            out.append(data["result"][k])
        time.sleep(0.4)  # be polite
    return out

QUERIES = [
    'thyroid AND spatial[All Fields]',
    'thyroid AND visium[All Fields]',
    'papillary thyroid AND spatial[All Fields]',
    'PTC AND visium[All Fields]',
    'thyroid carcinoma AND spatial transcriptomics[All Fields]',
    'thyroid AND 10x Genomics AND spatial[All Fields]',
    'thyroid AND spatial[Title]',
    'anaplastic thyroid AND spatial[All Fields]',
    'Hashimoto AND spatial[All Fields]',
    'Graves AND spatial[All Fields]',
    'thyroid AND CosMx[All Fields]',
    'thyroid AND Xenium[All Fields]',
    'thyroid AND GeoMx[All Fields]',
    'thyroid AND slide-seq[All Fields]',
    'thyroid AND merfish[All Fields]',
]

print("=" * 78)
print("GEO search — thyroid spatial transcriptomics")
print(f"Already have: {ALREADY_HAVE}")
print("=" * 78)

all_uids = set()
per_query = {}
for q in QUERIES:
    try:
        ids, count = esearch("gds", q, retmax=200)
        per_query[q] = (ids, count)
        all_uids.update(ids)
        print(f"  {q:60s}  hits={count:4d}  fetched={len(ids)}")
        time.sleep(0.4)
    except Exception as e:
        print(f"  {q}  ERROR: {e}")
        per_query[q] = ([], 0)

print(f"\n[summary] unique UIDs across queries: {len(all_uids)}")

# Pull summaries
print("[fetch] esummary for all UIDs ...")
records = esummary("gds", sorted(all_uids))
print(f"  fetched {len(records)} summaries")

# Save raw
with open(OUT/"geo_search_raw_summaries.json", "w") as f:
    json.dump(records, f, indent=1, default=str)

# Build TSV
import csv
rows = []
for r in records:
    acc = r.get("accession", "")
    title = r.get("title", "")
    summary = r.get("summary", "")
    gpl = r.get("gpl", "")
    n_samples = r.get("n_samples", "")
    pdat = r.get("pdat", "")
    taxon = r.get("taxon", "")
    entry_type = r.get("entrytype", "")
    suppl = r.get("suppfile", "")
    blob = " ".join([str(title).lower(), str(summary).lower(), str(suppl).lower()])
    is_spatial = any(k in blob for k in ["visium","spatial","slide-seq","slide seq","10x genomics spatial","st-seq","merfish","cosmx","xenium","geomx","spatial transcriptom","spatially-resolved","spot barcoded"])
    is_thyroid = any(k in blob for k in ["thyroid","ptc","ftc","atc","papillary thyroid","hashimoto","graves","tcga-thca","follicular thyroid"])
    novel = acc not in ALREADY_HAVE
    rows.append({
        "accession": acc,
        "is_thyroid": int(is_thyroid),
        "is_spatial": int(is_spatial),
        "novel": int(novel),
        "n_samples": n_samples,
        "taxon": taxon,
        "platform": gpl,
        "publication_date": pdat,
        "entry_type": entry_type,
        "title": title,
        "summary": summary[:380],
        "suppfile": suppl[:200],
    })

# Sort: thyroid + spatial + novel first
rows.sort(key=lambda x: (-x["is_thyroid"], -x["is_spatial"], -x["novel"], x["accession"]))

with open(OUT/"geo_search_results.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    keys = ["accession","is_thyroid","is_spatial","novel","n_samples","taxon",
            "platform","publication_date","entry_type","title","summary","suppfile"]
    w.writerow(keys)
    for r in rows:
        w.writerow([r[k] for k in keys])

# Build candidate list
cands = [r for r in rows if r["is_thyroid"] and r["is_spatial"] and r["novel"]]
print(f"\n[CANDIDATES] thyroid + spatial + novel: {len(cands)}")
for c in cands:
    print(f"  {c['accession']}  n={c['n_samples']}  [{c['platform']}]  {c['title'][:120]}")

with open(OUT/"geo_candidates.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    keys = ["accession","n_samples","platform","publication_date","entry_type","title","summary","suppfile"]
    w.writerow(keys)
    for r in cands:
        w.writerow([r[k] for k in keys])

# Also report borderline (thyroid only, no spatial detected — could be sc but adjacent)
border = [r for r in rows if r["is_thyroid"] and not r["is_spatial"] and r["novel"] and "single cell" in (r["title"]+r["summary"]).lower()]
print(f"\n[ADJACENT] thyroid + scRNA-seq (might pair w/ spatial): {len(border)}")
for c in border[:20]:
    print(f"  {c['accession']}  {c['title'][:120]}")

print(f"\n[output]\n  {OUT/'geo_search_raw_summaries.json'}\n  {OUT/'geo_search_results.tsv'}\n  {OUT/'geo_candidates.tsv'}")
