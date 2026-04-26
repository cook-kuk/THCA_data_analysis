#!/usr/bin/env python3
"""v16 — Tavily web search for press releases / FDA / conference talks
that PubMed misses. 8 targets + 5 drugs.
Output: results/v12_literature/news/{target_,drug_}*.json + combined TSV.
"""
import os, re, csv, json, time
from pathlib import Path
from urllib.request import Request, urlopen
from urllib.error import HTTPError

ROOT = Path("/opt/thyroid-dash/project")
RES  = ROOT / "results/v12_literature"
NEWS = RES / "news"
NEWS.mkdir(parents=True, exist_ok=True)

ENV_PATH = Path("/home/seungho/cook-forge/.env")
KEY = None
for line in ENV_PATH.read_text().splitlines():
    if line.startswith("TAVILY_API_KEY="):
        KEY = line.split("=",1)[1].strip().strip('"').strip("'")
if not KEY:
    raise SystemExit("TAVILY_API_KEY not found")

def tavily_search(query, max_results=8):
    body = json.dumps({
        "api_key": KEY,
        "query": query,
        "search_depth": "basic",
        "max_results": max_results,
        "include_answer": True,
        "topic": "general",
    }).encode()
    req = Request("https://api.tavily.com/search", data=body,
        headers={"content-type":"application/json"})
    for attempt in range(3):
        try:
            with urlopen(req, timeout=30) as r:
                return json.loads(r.read().decode())
        except HTTPError as e:
            err = e.read().decode(errors="ignore")[:200]
            if e.code in (429, 500, 502, 503): time.sleep(2 + attempt*3); continue
            raise RuntimeError(f"HTTP {e.code}: {err}")
    raise RuntimeError("Tavily retries exhausted")

TARGETS = ["CYP1B1","TACSTD2","TMPRSS4","PLEKHA6","LDLR","GABRB2","B3GNT3","PTPRE"]
DRUGS   = ["Cannabidiol","Sacituzumab govitecan","Tesevatinib","Dacomitinib","Tepotinib"]

target_queries = lambda g: f"{g} thyroid cancer FDA approval clinical trial 2024 2025 2026 pharma press release"
drug_queries   = lambda d: f"{d} thyroid cancer FDA approval phase 3 2024 2025 2026 clinical trial press release"

def first_answer(j): return (j or {}).get("answer") or ""

def main():
    t0=time.time()
    t_rows=[]; d_rows=[]
    for g in TARGETS:
        print(f"[T] {g}")
        try: j = tavily_search(target_queries(g), max_results=8)
        except Exception as e: print(f"  err {e}"); continue
        (NEWS / f"target_{g}.json").write_text(json.dumps(j, indent=2, ensure_ascii=False))
        ans = first_answer(j); results = j.get("results", [])
        # pick top-3 most relevant
        top = []
        for r in results[:5]:
            top.append({
                "title": r.get("title","")[:240],
                "url": r.get("url","")[:240],
                "domain": re.sub(r"^https?://([^/]+).*", r"\1", r.get("url","") or ""),
                "score": round(r.get("score",0), 3),
            })
        t_rows.append({"target":g, "answer":ans[:600], "n_results":len(results),
                       "top_url":top[0]["url"] if top else "",
                       "top_domain":top[0]["domain"] if top else "",
                       "top_score":top[0]["score"] if top else "",
                       "top_title":top[0]["title"] if top else "",
                       "all_top": "; ".join(f"{r['domain']}:{r['title'][:80]}" for r in top)})
        time.sleep(0.5)

    for d in DRUGS:
        print(f"[D] {d}")
        try: j = tavily_search(drug_queries(d), max_results=8)
        except Exception as e: print(f"  err {e}"); continue
        slug = re.sub(r"\W+","_",d).strip("_").lower()
        (NEWS / f"drug_{slug}.json").write_text(json.dumps(j, indent=2, ensure_ascii=False))
        ans = first_answer(j); results = j.get("results", [])
        top = [{"title":r.get("title","")[:240],"url":r.get("url","")[:240],
                "domain":re.sub(r"^https?://([^/]+).*", r"\1", r.get("url","") or ""),
                "score":round(r.get("score",0),3)} for r in results[:5]]
        d_rows.append({"drug":d, "answer":ans[:600], "n_results":len(results),
                       "top_url":top[0]["url"] if top else "",
                       "top_domain":top[0]["domain"] if top else "",
                       "top_score":top[0]["score"] if top else "",
                       "top_title":top[0]["title"] if top else "",
                       "all_top": "; ".join(f"{r['domain']}:{r['title'][:80]}" for r in top)})
        time.sleep(0.5)

    safe = lambda s: (s or "").replace("\t"," ").replace("\r"," ").replace("\n"," ¶ ")
    with (NEWS / "targets_news.tsv").open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["target","answer","n_results","top_score","top_domain","top_title","top_url","all_top"])
        for r in t_rows: w.writerow([r["target"], safe(r["answer"]), r["n_results"], r["top_score"], r["top_domain"], safe(r["top_title"]), r["top_url"], safe(r["all_top"])])
    with (NEWS / "drugs_news.tsv").open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["drug","answer","n_results","top_score","top_domain","top_title","top_url","all_top"])
        for r in d_rows: w.writerow([r["drug"], safe(r["answer"]), r["n_results"], r["top_score"], r["top_domain"], safe(r["top_title"]), r["top_url"], safe(r["all_top"])])

    idx = {"built_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
           "n_targets": len(t_rows), "n_drugs": len(d_rows),
           "elapsed_seconds": round(time.time()-t0, 1)}
    (NEWS / "v16_index.json").write_text(json.dumps(idx, indent=2, ensure_ascii=False))
    print("[done]", json.dumps(idx, indent=2))

if __name__ == "__main__": main()
