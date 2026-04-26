#!/usr/bin/env python3
"""v14 preprint scoop scan — Europe PMC PPR (covers bioRxiv, medRxiv,
ChemRxiv, Research Square, etc.). Runs the same conceptual queries as
the prior PubMed pass but restricted to preprints, plus a temporal sweep.
Writes raw JSON + flat TSV to results/v14_priorart/.
"""
from __future__ import annotations
import json, time, urllib.parse, urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "v14_priorart"
OUT.mkdir(parents=True, exist_ok=True)

EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
UA = {"User-Agent": "thca-priorart-check/1.0 (mailto:kukshomr@gmail.com)"}

QUERIES = {
    "P1_BRAF_TROP2_thyroid":      "(BRAF AND TROP2 AND thyroid) AND SRC:PPR",
    "P2_BRAF_TACSTD2_thyroid":    "(BRAF AND TACSTD2 AND thyroid) AND SRC:PPR",
    "P3_V600E_TROP2":             "(V600E AND TROP2) AND SRC:PPR",
    "P4_thyroid_sacituzumab":     "(thyroid AND sacituzumab) AND SRC:PPR",
    "P5_thyroid_TROP2_ADC":       "(thyroid AND TROP2 AND (ADC OR \"antibody-drug conjugate\")) AND SRC:PPR",
    "P6_BRAF_irinotecan_thyroid": "(BRAF AND irinotecan AND thyroid) AND SRC:PPR",
    "P7_BRAF_topotecan_thyroid":  "(BRAF AND topotecan AND thyroid) AND SRC:PPR",
    "P8_PRISM_BRAF_thyroid":      "(PRISM AND BRAF AND thyroid) AND SRC:PPR",
    "P9_thyroid_TROP2_general":   "(thyroid AND TROP2) AND SRC:PPR",  # broad sweep
    "P10_PTC_TROP2":              "(\"papillary thyroid\" AND TROP2) AND SRC:PPR",
}

def http_get(url: str, timeout: int = 30) -> bytes:
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()

def epmc_search(query: str, page_size: int = 25) -> dict:
    qs = urllib.parse.urlencode({
        "query": query,
        "format": "json",
        "pageSize": page_size,
        "resultType": "core",
    })
    return json.loads(http_get(f"{EPMC}?{qs}"))

def main():
    print("[Europe PMC preprint scan]")
    all_results = {}
    for label, query in QUERIES.items():
        try:
            d = epmc_search(query)
            n = d.get("hitCount", 0)
            results = d.get("resultList", {}).get("result", [])
            hits = []
            for r in results:
                hits.append({
                    "id": r.get("id"),
                    "source": r.get("source"),  # PPR
                    "preprint_repo": r.get("bookOrReportDetails",{}).get("publisher") or r.get("publisher",""),
                    "year": r.get("pubYear"),
                    "title": r.get("title",""),
                    "authors": r.get("authorString","")[:200],
                    "journal_or_repo": r.get("journalTitle") or r.get("source"),
                    "abstract": (r.get("abstractText","") or "")[:1500],
                    "doi": r.get("doi"),
                    "url": r.get("fullTextUrlList",{}).get("fullTextUrl",[{}])[0].get("url") if r.get("fullTextUrlList") else None,
                })
            all_results[label] = {"query": query, "hit_count": n, "shown": len(hits), "hits": hits}
            print(f"  {label:32s}: {n} preprint hits (top {len(hits)} captured)")
        except Exception as e:
            all_results[label] = {"query": query, "error": str(e)}
            print(f"  {label:32s}: ERROR {e}")
        time.sleep(0.5)

    (OUT / "preprint_raw.json").write_text(json.dumps(all_results, indent=2, ensure_ascii=False))

    # Flat TSV
    tsv = ["query_label\tquery\tepmc_id\tsource\tyear\ttitle\tauthors\tabstract_snippet\tdoi"]
    for lbl, r in all_results.items():
        if "error" in r: continue
        for h in r["hits"]:
            ab = (h["abstract"] or "").replace("\t"," ").replace("\n"," ")[:400]
            ti = (h["title"] or "").replace("\t"," ")
            au = (h["authors"] or "").replace("\t"," ")[:120]
            tsv.append(f"{lbl}\t{r['query']}\t{h['id']}\t{h['source']}\t{h['year']}\t{ti}\t{au}\t{ab}\t{h['doi'] or ''}")
    (OUT / "preprint_scan.tsv").write_text("\n".join(tsv))
    print(f"  wrote preprint_scan.tsv ({len(tsv)-1} rows)")

    # Summary
    print("\n=== HIT COUNTS ===")
    total_hits = 0
    relevant_q = ["P1_BRAF_TROP2_thyroid","P2_BRAF_TACSTD2_thyroid","P3_V600E_TROP2",
                  "P4_thyroid_sacituzumab","P5_thyroid_TROP2_ADC","P6_BRAF_irinotecan_thyroid",
                  "P7_BRAF_topotecan_thyroid","P8_PRISM_BRAF_thyroid","P9_thyroid_TROP2_general",
                  "P10_PTC_TROP2"]
    for lbl in relevant_q:
        r = all_results.get(lbl, {})
        n = r.get("hit_count", 0) if "error" not in r else 0
        total_hits += n
        print(f"  {lbl:32s}: {n}")
    print(f"\n  Total preprint hits: {total_hits}")

if __name__ == "__main__":
    main()
