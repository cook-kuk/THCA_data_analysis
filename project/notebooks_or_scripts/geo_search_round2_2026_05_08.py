#!/usr/bin/env python3
"""Round 2 search — Europe PMC + ArrayExpress + tighter GEO queries.
Goal: find Visium / spot-level human thyroid carcinoma spatial datasets we don't already have.
"""
from pathlib import Path
import urllib.request, urllib.parse, json, time, re

OUT = Path("/home/seungho/personal/THCA_data_analysis/project/results/geo_search_2026_05_08")
ALREADY = {"GSE250521","GSE230424","GSE248205","GSE301163"}  # 301163 = GeoMx not Visium

# ---------- Europe PMC ----------
print("="*78)
print("Europe PMC — recent Visium/Slide-seq thyroid PTC/ATC/HT papers")
print("="*78)
EPMC_QUERIES = [
    '("papillary thyroid" OR PTC OR "thyroid carcinoma") AND (Visium OR "10x spatial" OR "spatial transcriptom*")',
    '("anaplastic thyroid" OR ATC) AND (Visium OR "spatial transcriptom*")',
    '"thyroid" AND ("slide-seq" OR Stereo-seq OR Xenium OR CosMx OR MERFISH) AND (carcinoma OR cancer OR PTC)',
    '"Hashimoto" AND (Visium OR "spatial transcriptom*")',
    '"Graves" AND (Visium OR "spatial transcriptom*")',
    '"thyroid" AND Visium AND FFPE',
    'thyroid AND Visium AND 2024',
    'thyroid AND Visium AND 2025',
    'thyroid AND Visium AND 2026',
]
EPMC_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query={q}&format=json&pageSize=25&resultType=core"
all_papers = {}
for q in EPMC_QUERIES:
    try:
        url = EPMC_URL.format(q=urllib.parse.quote(q))
        with urllib.request.urlopen(url, timeout=25) as r:
            data = json.loads(r.read().decode())
        results = data.get("resultList", {}).get("result", [])
        print(f"\n  Q: {q}\n   hits={len(results)}")
        for p in results:
            pmid = p.get("pmid", "") or p.get("id", "")
            title = p.get("title", "")[:140]
            authors = p.get("authorString", "")[:60]
            year = p.get("pubYear", "")
            journal = p.get("journalTitle", "")[:40]
            access = p.get("accessionList", {}).get("accession", []) if isinstance(p.get("accessionList"), dict) else []
            access_codes = [a.get("accession") for a in access] if isinstance(access, list) else []
            if pmid not in all_papers:
                all_papers[pmid] = {"title":title, "authors":authors, "year":year, "journal":journal, "access":access_codes}
            print(f"    [{year}] {journal:30s} {title[:100]}")
            if access_codes:
                for ac in access_codes:
                    print(f"        accession: {ac}")
        time.sleep(0.5)
    except Exception as e:
        print(f"  Q: {q}  ERROR: {e}")

# Save
import csv
with open(OUT/"epmc_round2.tsv","w",newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["pmid","year","journal","title","authors","accessions"])
    for pmid, p in all_papers.items():
        w.writerow([pmid, p["year"], p["journal"], p["title"], p["authors"], ";".join(p["access"])])

# ---------- ArrayExpress / BioStudies ----------
print("\n"+"="*78)
print("BioStudies / ArrayExpress — thyroid spatial")
print("="*78)
BIO_URL = "https://www.ebi.ac.uk/biostudies/api/v1/search?query={q}&pageSize=50&type=study"
BIO_QUERIES = [
    "thyroid Visium",
    "thyroid spatial transcriptomics",
    "papillary thyroid spatial",
    "PTC Visium",
    "anaplastic thyroid spatial",
]
ae_hits = []
for q in BIO_QUERIES:
    try:
        url = BIO_URL.format(q=urllib.parse.quote(q))
        with urllib.request.urlopen(url, timeout=25) as r:
            data = json.loads(r.read().decode())
        hits = data.get("hits", [])
        print(f"\n  Q: {q}  hits={len(hits)}")
        for h in hits[:30]:
            acc = h.get("accession","")
            title = h.get("title","")[:160]
            ae_hits.append({"acc":acc,"title":title,"q":q})
            print(f"    {acc}  {title[:120]}")
        time.sleep(0.4)
    except Exception as e:
        print(f"  Q: {q}  ERROR: {e}")
with open(OUT/"biostudies_round2.tsv","w",newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["accession","title","query"])
    for h in ae_hits:
        w.writerow([h["acc"], h["title"], h["q"]])

# ---------- Tighter GEO query ----------
print("\n"+"="*78)
print("GEO — tighter Visium thyroid carcinoma query")
print("="*78)
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
def esearch(db, term, retmax=200):
    url = f"{EUTILS}/esearch.fcgi?db={db}&term={urllib.parse.quote(term)}&retmax={retmax}&retmode=json"
    with urllib.request.urlopen(url, timeout=20) as r:
        data = json.loads(r.read().decode())
    return data.get("esearchresult", {}).get("idlist", []), int(data.get("esearchresult", {}).get("count", 0))

TIGHT = [
    'thyroid[All Fields] AND ("10x visium"[All Fields] OR "10x genomics visium"[All Fields])',
    '"papillary thyroid"[All Fields] AND ("10x visium"[All Fields] OR visium[All Fields])',
    '"thyroid carcinoma"[All Fields] AND visium[All Fields]',
    'thyroid AND Visium AND ("Spatial Gene Expression"[All Fields] OR "spatially resolved"[All Fields])',
    '("homo sapiens"[Organism]) AND thyroid AND visium',
]
for q in TIGHT:
    try:
        ids, count = esearch("gds", q)
        print(f"  hits={count:4d}  Q: {q}")
        time.sleep(0.4)
    except Exception as e:
        print(f"  Q: {q}  ERROR: {e}")

# ---------- ENA / SRA — PRJNA / PRJEB lookups ----------
# For datasets indexed in ENA (e.g., 10x raw FASTQs)
print("\n"+"="*78)
print("ENA — thyroid Visium runs")
print("="*78)
ENA_URL = "https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=" \
          "%28study_title%3D%22*thyroid*Visium*%22%20OR%20study_title%3D%22*thyroid*spatial*transcriptom*%22%29" \
          "&format=tsv&limit=200&fields=run_accession,study_accession,sample_accession,study_title,sample_title,instrument_platform,library_strategy"
try:
    with urllib.request.urlopen(ENA_URL, timeout=30) as r:
        ena_text = r.read().decode("utf-8", errors="ignore")
    lines = ena_text.strip().split("\n")
    print(f"  ENA rows: {len(lines)-1 if len(lines)>1 else 0}")
    for line in lines[:25]:
        print("   ", line[:200])
    with open(OUT/"ena_thyroid_visium.tsv","w") as f:
        f.write(ena_text)
except Exception as e:
    print(f"  ENA ERROR: {e}")

print(f"\n[output dir] {OUT}")
