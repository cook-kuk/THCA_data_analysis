#!/usr/bin/env python3
"""v14 prior-art query: PubMed E-utilities + ClinicalTrials.gov v2.
Runs 8 PubMed queries (T1+T2) and 2 CTGov fetches (T3) in parallel,
writes raw JSON + a flat TSV summary to results/v14_priorart/.
"""
from __future__ import annotations
import json, sys, time, urllib.parse, urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "v14_priorart"
OUT.mkdir(parents=True, exist_ok=True)

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
CTGOV  = "https://clinicaltrials.gov/api/v2/studies"
UA = {"User-Agent": "thca-priorart-check/1.0 (mailto:kukshomr@gmail.com)"}

QUERIES = {
    "A_BRAF_TROP2_thyroid":     "(BRAF) AND (TROP2) AND (thyroid)",
    "B_BRAF_TACSTD2_thyroid":   "(BRAF) AND (TACSTD2) AND (thyroid)",
    "C_V600E_TROP2":            "(V600E) AND (TROP2)",
    "D_thyroid_sacituzumab":    "(thyroid) AND (sacituzumab)",
    "E_thyroid_TROP2_ADC":      "(thyroid) AND (TROP2) AND (ADC OR antibody-drug conjugate)",
    "F_BRAF_irinotecan_cancer": "(BRAF) AND (irinotecan) AND (cancer)",
    "G_BRAF_topotecan":         "(BRAF) AND (topotecan)",
    "H_SN38_BRAF_selectivity":  "(SN-38 OR SN38) AND (BRAF) AND (selectivity OR sensitivity)",
    # Task 2 extras
    "T2a_irinotecan_BRAF_CRC":  "(irinotecan) AND (BRAF) AND (colorectal)",
    "T2b_topoI_BRAF_sens":      "(topoisomerase I) AND (BRAF) AND (sensitivity)",
    "T2c_SN38_DepMap":          "(SN-38 OR SN38) AND (DepMap)",
    "T2d_PRISM_BRAF":           "(PRISM) AND (BRAF) AND (drug sensitivity)",
}

def http_get(url: str, timeout: int = 30) -> bytes:
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()

def esearch(term: str, retmax: int = 20) -> list[str]:
    qs = urllib.parse.urlencode({"db":"pubmed","term":term,"retmax":retmax,"retmode":"json"})
    data = json.loads(http_get(f"{EUTILS}/esearch.fcgi?{qs}"))
    return data.get("esearchresult", {}).get("idlist", [])

def esummary(pmids: list[str]) -> dict:
    if not pmids: return {}
    qs = urllib.parse.urlencode({"db":"pubmed","id":",".join(pmids),"retmode":"json"})
    return json.loads(http_get(f"{EUTILS}/esummary.fcgi?{qs}")).get("result", {})

def efetch_abstracts(pmids: list[str]) -> dict[str,str]:
    """Return {pmid: abstract_text} via efetch in XML, parsed minimally."""
    if not pmids: return {}
    qs = urllib.parse.urlencode({"db":"pubmed","id":",".join(pmids),"rettype":"abstract","retmode":"xml"})
    xml = http_get(f"{EUTILS}/efetch.fcgi?{qs}").decode("utf-8", errors="ignore")
    out: dict[str,str] = {}
    # crude split per <PubmedArticle> — sufficient for first-pass relevance
    import re
    blocks = re.split(r"</PubmedArticle>", xml)
    for b in blocks:
        m_pmid = re.search(r"<PMID[^>]*>(\d+)</PMID>", b)
        if not m_pmid: continue
        pmid = m_pmid.group(1)
        # collect AbstractText (may be multiple labelled sections)
        parts = re.findall(r"<AbstractText[^>]*>(.*?)</AbstractText>", b, flags=re.S)
        text = " ".join(re.sub(r"<[^>]+>","", p) for p in parts)
        out[pmid] = re.sub(r"\s+"," ", text).strip()
    return out

def run_one_query(label: str, term: str) -> dict:
    t0 = time.time()
    try:
        pmids = esearch(term, retmax=20)
        time.sleep(0.34)  # 3 req/sec polite
        summary = esummary(pmids) if pmids else {}
        time.sleep(0.34)
        abstracts = efetch_abstracts(pmids) if pmids else {}
    except Exception as e:
        return {"label": label, "term": term, "error": str(e)}
    hits = []
    for pmid in pmids:
        s = summary.get(pmid, {})
        hits.append({
            "pmid": pmid,
            "title": s.get("title",""),
            "year": (s.get("pubdate","")[:4] if s.get("pubdate") else ""),
            "journal": s.get("fulljournalname") or s.get("source",""),
            "abstract": (abstracts.get(pmid,"") or "")[:1500],
        })
    return {"label": label, "term": term, "n_hits": len(hits), "hits": hits, "elapsed_s": round(time.time()-t0,1)}

def fetch_ctgov(nct: str) -> dict:
    try:
        return json.loads(http_get(f"{CTGOV}/{nct}"))
    except Exception as e:
        return {"nct": nct, "error": str(e)}

def main():
    print("[1/3] PubMed queries (serial, polite delay)…")
    pubmed_results: dict[str, dict] = {}
    for lbl, term in QUERIES.items():
        # retry on 429
        for attempt in range(3):
            r = run_one_query(lbl, term)
            if "error" in r and "429" in r["error"]:
                wait = 3 * (attempt+1)
                print(f"  - {lbl}: 429, retry in {wait}s")
                time.sleep(wait)
                continue
            break
        pubmed_results[lbl] = r
        if "error" in r:
            print(f"  - {lbl}: ERROR {r['error']}")
        else:
            print(f"  - {lbl}: {r['n_hits']} hits ({r['elapsed_s']}s)")
        time.sleep(1.2)  # between distinct queries

    (OUT / "pubmed_raw.json").write_text(json.dumps(pubmed_results, indent=2, ensure_ascii=False))

    print("[2/3] ClinicalTrials.gov fetches…")
    nct_ids = ["NCT06235216", "NCT07521670"]
    ct_results = {nct: fetch_ctgov(nct) for nct in nct_ids}
    (OUT / "ctgov_raw.json").write_text(json.dumps(ct_results, indent=2, ensure_ascii=False))
    for nct, d in ct_results.items():
        if "error" in d:
            print(f"  - {nct}: ERROR {d['error']}")
        else:
            ps = d.get("protocolSection", {})
            status = ps.get("statusModule", {}).get("overallStatus","?")
            phase = ps.get("designModule", {}).get("phases",["?"])
            print(f"  - {nct}: {status} / phase {phase}")

    print("[3/3] Flat TSV summary…")
    tsv_lines = ["query_label\tquery_term\tpmid\tyear\tjournal\ttitle\tabstract_snippet"]
    for lbl, r in pubmed_results.items():
        if "error" in r:
            tsv_lines.append(f"{lbl}\t{r['term']}\tERROR\t\t\t{r['error']}\t")
            continue
        for h in r["hits"]:
            ab = (h["abstract"] or "").replace("\t"," ").replace("\n"," ")[:500]
            ti = (h["title"] or "").replace("\t"," ")
            jrn = (h["journal"] or "").replace("\t"," ")
            tsv_lines.append(f"{lbl}\t{r['term']}\t{h['pmid']}\t{h['year']}\t{jrn}\t{ti}\t{ab}")
    (OUT / "braf_trop2_thyroid_priorart.tsv").write_text("\n".join(tsv_lines))
    print(f"  wrote {OUT/'braf_trop2_thyroid_priorart.tsv'} ({len(tsv_lines)-1} rows)")

    # Brief stdout digest
    print("\n=== HIT COUNTS ===")
    total = 0
    for lbl, r in pubmed_results.items():
        n = r.get("n_hits", 0) if "error" not in r else 0
        total += n
        print(f"  {lbl:32s}: {n}")
    print(f"  TOTAL: {total} pubmed records")

if __name__ == "__main__":
    main()
