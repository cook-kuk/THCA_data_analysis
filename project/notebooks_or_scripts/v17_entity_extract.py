#!/usr/bin/env python3
"""v17 — per-abstract entity extraction via claude-haiku.
For 8 targets, take top-6 abstracts each, ask Claude to extract structured fields:
  direction · effect_size_qual · mechanism · trial_phase · cancer_context · key_finding
Output: results/v12_literature/llm/entity/{GENE}.tsv + combined v17_entity_table.tsv.
"""
import os, re, csv, json, time
from pathlib import Path
from urllib.request import Request, urlopen
from urllib.error import HTTPError

ROOT = Path("/opt/thyroid-dash/project")
RES  = ROOT / "results/v12_literature"
PUB  = RES / "pubmed_evidence"
ENT  = RES / "llm/entity"
ENT.mkdir(parents=True, exist_ok=True)

ENV_PATH = Path("/home/seungho/cook-forge/.env")
KEY = None; HAIKU = "claude-haiku-4-5-20251001"
for line in ENV_PATH.read_text().splitlines():
    if line.startswith("ANTHROPIC_API_KEY="): KEY = line.split("=",1)[1].strip().strip('"').strip("'")
    if line.startswith("ANTHROPIC_HAIKU_MODEL="):
        HAIKU = line.split("=",1)[1].strip().strip('"').strip("'")

def claude_json(system, user, max_tokens=400):
    body = json.dumps({"model":HAIKU,"max_tokens":max_tokens,"system":system,
        "messages":[{"role":"user","content":user}]}).encode()
    req = Request("https://api.anthropic.com/v1/messages", data=body,
        headers={"x-api-key":KEY,"anthropic-version":"2023-06-01","content-type":"application/json"})
    for attempt in range(3):
        try:
            with urlopen(req, timeout=60) as r:
                d = json.loads(r.read().decode())
                txt = "".join(b.get("text","") for b in d.get("content",[]))
                # extract JSON from response (may be in code block)
                m = re.search(r"\{[\s\S]+\}", txt)
                obj = json.loads(m.group(0)) if m else {"raw": txt}
                return obj, d.get("usage",{})
        except HTTPError as e:
            err = e.read().decode(errors="ignore")[:200]
            if e.code in (429,500,502,503): time.sleep(2+attempt*3); continue
            raise RuntimeError(f"HTTP {e.code}: {err}")
        except json.JSONDecodeError:
            return {"raw": txt}, d.get("usage",{}) if 'd' in dir() else {}
    raise RuntimeError("retries exhausted")

TARGETS = ["CYP1B1","TACSTD2","TMPRSS4","PLEKHA6","LDLR","GABRB2","B3GNT3","PTPRE"]

SYS = ("You extract structured entities from a single biomedical abstract.\n"
"Return ONLY a JSON object with these keys (no commentary, no markdown fence):\n"
'{"direction": "up|down|mixed|unknown",'
' "effect_size_qual": "strong|moderate|weak|unknown",'
' "mechanism": "<short string, <=80 chars, or unknown>",'
' "trial_phase": "preclinical|phase_1|phase_2|phase_3|approved|none_mentioned",'
' "cancer_context": "thyroid|breast|lung|colorectal|melanoma|glioma|pan_cancer|other|none",'
' "key_finding": "<single sentence, <=140 chars>"}\n'
"If a field is not addressable from the abstract, use 'unknown' or 'none_mentioned'. Be terse and conservative.")

def top_abstracts(gene, max_n=6):
    p = PUB / f"{gene}_pubmed_hits.tsv"
    if not p.exists(): return []
    with p.open() as f: rows = list(csv.DictReader(f, delimiter="\t"))
    def score(r):
        s = 0
        if r.get("clinical")=="1": s += 3
        if r.get("druggable")=="1": s += 2
        if r.get("direction"): s += 1
        try: s += min(int(r.get("year","0") or 0)-2014, 12)
        except: pass
        return -s
    rows.sort(key=score)
    return rows[:max_n]

def main():
    t0=time.time()
    all_rows=[]; tok_in=tok_out=ncalls=0
    for g in TARGETS:
        gene_rows=[]
        for r in top_abstracts(g):
            user = f"GENE: {g}\nABSTRACT TITLE: {r.get('title','')}\nYEAR: {r.get('year','')}"
            try: obj, u = claude_json(SYS, user, max_tokens=320)
            except Exception as e: print(f"  err {g}/{r['pmid']}: {e}"); continue
            ncalls += 1
            tok_in  += u.get("input_tokens",0) or 0
            tok_out += u.get("output_tokens",0) or 0
            row = {"target":g, "pmid":r["pmid"], "year":r.get("year",""),
                   "title":r.get("title","")[:200],
                   "direction":obj.get("direction",""),
                   "effect_size_qual":obj.get("effect_size_qual",""),
                   "mechanism":obj.get("mechanism",""),
                   "trial_phase":obj.get("trial_phase",""),
                   "cancer_context":obj.get("cancer_context",""),
                   "key_finding":obj.get("key_finding","")}
            gene_rows.append(row); all_rows.append(row)
        # per-gene tsv
        if gene_rows:
            with (ENT / f"{g}.tsv").open("w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=list(gene_rows[0].keys()), delimiter="\t")
                w.writeheader(); [w.writerow(r) for r in gene_rows]
            print(f"[{g}] {len(gene_rows)} abstracts extracted")

    if all_rows:
        with (ENT / "v17_entity_table.tsv").open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(all_rows[0].keys()), delimiter="\t")
            w.writeheader(); [w.writerow(r) for r in all_rows]

    idx = {"built_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
           "model": HAIKU, "calls": ncalls,
           "input_tokens_total": tok_in, "output_tokens_total": tok_out,
           "abstracts_extracted": len(all_rows),
           "elapsed_seconds": round(time.time()-t0, 1)}
    (ENT / "v17_index.json").write_text(json.dumps(idx, indent=2, ensure_ascii=False))
    print("[done]", json.dumps(idx, indent=2))

if __name__ == "__main__": main()
