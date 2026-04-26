#!/usr/bin/env python3
"""v15 — re-run v14 target+drug verdicts with claude-sonnet for nuance.
Output: results/v12_literature/llm/sonnet/{target_,drug_}*.md + combined TSVs.
"""
import os, re, csv, json, time
from pathlib import Path
from urllib.request import Request, urlopen
from urllib.error import HTTPError

ROOT = Path("/opt/thyroid-dash/project")
RES  = ROOT / "results/v12_literature"
PUB  = RES / "pubmed_evidence"
LLM  = RES / "llm"
SON  = LLM / "sonnet"
SON.mkdir(parents=True, exist_ok=True)

ENV_PATH = Path("/home/seungho/cook-forge/.env")
KEY = None; SONNET_MODEL = "claude-sonnet-4-5-20250929"
for line in ENV_PATH.read_text().splitlines():
    if line.startswith("ANTHROPIC_API_KEY="): KEY = line.split("=",1)[1].strip().strip('"').strip("'")
    if line.startswith("ANTHROPIC_SONNET_MODEL="):
        SONNET_MODEL = line.split("=",1)[1].strip().strip('"').strip("'")
    if line.startswith("ANTHROPIC_MODEL="):
        # often points to sonnet
        v = line.split("=",1)[1].strip().strip('"').strip("'")
        if "sonnet" in v.lower(): SONNET_MODEL = v
print(f"[v15] using model={SONNET_MODEL}")

def claude(system, user, max_tokens=600):
    body = json.dumps({"model":SONNET_MODEL,"max_tokens":max_tokens,
        "system":system,"messages":[{"role":"user","content":user}]}).encode()
    req = Request("https://api.anthropic.com/v1/messages", data=body,
        headers={"x-api-key":KEY,"anthropic-version":"2023-06-01","content-type":"application/json"})
    for attempt in range(3):
        try:
            with urlopen(req, timeout=90) as r:
                d = json.loads(r.read().decode())
                txt = "".join(b.get("text","") for b in d.get("content",[]))
                return txt.strip(), d.get("usage",{})
        except HTTPError as e:
            err = e.read().decode(errors="ignore")[:200]
            if e.code in (429,529,500,502,503): time.sleep(3+attempt*5); continue
            raise RuntimeError(f"HTTP {e.code}: {err}")
    raise RuntimeError("retries exhausted")

TARGETS = ["CYP1B1","TACSTD2","TMPRSS4","PLEKHA6","LDLR","GABRB2","B3GNT3","PTPRE"]
DRUGS   = ["Cannabidiol","Sacituzumab govitecan","Tesevatinib","Dacomitinib","Tepotinib"]

def gene_context(gene, max_n=12, max_chars=3500):
    p = PUB / f"{gene}_pubmed_hits.tsv"
    if not p.exists(): return ""
    with p.open() as f: rows = list(csv.DictReader(f, delimiter="\t"))
    def score(r):
        s=0
        if r.get("clinical")=="1": s+=4
        if r.get("druggable")=="1": s+=2
        if r.get("direction"): s+=2
        if "braf" in (r.get("query_types","") or "").lower(): s+=3
        try: s += min(int(r.get("year","0") or 0)-2014, 12)
        except: pass
        return -s
    rows.sort(key=score)
    out=[]; used=0
    for r in rows[:max_n]:
        s = f"PMID {r['pmid']} ({r.get('year','')}, {r.get('journal','')[:60]}): {r.get('title','')[:280]}"
        if used + len(s) > max_chars: break
        out.append(s); used += len(s)
    return "\n".join(out)

SYS_T = (
"You are a biomedical research synthesist preparing for a Bioinformatics paper on a thyroid-cancer "
"computational pipeline (BRAF-like vs RAS-like classification). Read the supplied PubMed abstract "
"snippets for ONE candidate druggable target and produce a publication-ready 4-paragraph evidence "
"summary with these exact labels:\n"
"BIOLOGY: 2-3 sentences on gene function and BRAF-like / MAPK / PTC-context mechanistic role. Be specific about pathway involvement.\n"
"DRUGGABILITY: 2-3 sentences naming specific small molecules / antibodies / ADCs / repurposed drugs and their development stage. Cite drug class.\n"
"CLINICAL: 1-2 sentences on existing clinical trial evidence in thyroid cancer specifically; if absent, say so explicitly.\n"
"RESEARCH GAP: 1-2 sentences on the most useful next experiment for a translational research group.\n"
"Total <=200 words. Do not invent specific PMIDs or NCT numbers. Be conservative when evidence is sparse.")

SYS_D = (
"You are a biomedical research synthesist evaluating a drug repurposing candidate for BRAF-like "
"thyroid cancer. From PubMed abstract titles, write four labelled paragraphs:\n"
"MECHANISM: 2-3 sentences on molecular target / pharmacology.\n"
"THYROID EVIDENCE: 1-2 sentences specifically on thyroid cancer literature (be honest if sparse).\n"
"CROSS-CANCER: 2 sentences on broader oncology evidence and approved indications.\n"
"REPURPOSING RATIONALE: 2 sentences on the strongest argument and the strongest counter-argument for thyroid use.\n"
"Total <=180 words. No invented citations.")

def synth_target(g):
    ctx = gene_context(g)
    if not ctx: return None
    return claude(SYS_T, f"GENE: {g}\n\nABSTRACTS:\n{ctx}\n\nWrite the four labelled paragraphs.", max_tokens=520)

def drug_titles(drug, max_n=8):
    import urllib.parse, urllib.request
    q = urllib.parse.urlencode({"db":"pubmed",
        "term":f'"{drug}"[tiab] AND cancer[tiab]',"retmode":"json","retmax":max_n,
        "tool":"thyrai-v15","email":"kukshomr@gmail.com"})
    try:
        with urllib.request.urlopen(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?{q}", timeout=20) as r:
            ids = json.loads(r.read().decode()).get("esearchresult",{}).get("idlist",[])
    except: return ""
    if not ids: return ""
    q2 = urllib.parse.urlencode({"db":"pubmed","id":",".join(ids),"retmode":"xml",
        "tool":"thyrai-v15","email":"kukshomr@gmail.com"})
    try:
        with urllib.request.urlopen(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?{q2}", timeout=30) as r:
            xml = r.read().decode()
    except: return ""
    titles = re.findall(r"<ArticleTitle[^>]*>(.*?)</ArticleTitle>", xml, re.S)
    years  = re.findall(r"<PubDate>.*?<Year>(\d{4})</Year>", xml, re.S)
    return "\n".join(f"({years[i] if i<len(years) else ''}) {re.sub(r'<[^>]+>','',t).strip()[:280]}" for i,t in enumerate(titles[:max_n]))

def synth_drug(d):
    ctx = drug_titles(d) or "(no recent abstracts retrieved)"
    return claude(SYS_D, f"DRUG: {d}\n\nABSTRACT TITLES:\n{ctx}\n\nWrite the four labelled paragraphs.", max_tokens=460)

def main():
    t0=time.time(); ts=[]; ds=[]
    for g in TARGETS:
        try: txt,u = synth_target(g)
        except Exception as e: print(f"[err] {g}: {e}"); continue
        if not txt: print(f"[skip] {g}"); continue
        (SON / f"target_{g}.md").write_text(f"# {g} (sonnet)\n\n{txt}\n", encoding="utf-8")
        ts.append({"target":g,"synthesis":txt,"input_tokens":u.get("input_tokens","?"),"output_tokens":u.get("output_tokens","?")})
        print(f"[T] {g} ✓ in={u.get('input_tokens','?')} out={u.get('output_tokens','?')}")
    for d in DRUGS:
        try: txt,u = synth_drug(d)
        except Exception as e: print(f"[err] {d}: {e}"); continue
        slug = re.sub(r"\W+","_",d).strip("_").lower()
        (SON / f"drug_{slug}.md").write_text(f"# {d} (sonnet)\n\n{txt}\n", encoding="utf-8")
        ds.append({"drug":d,"synthesis":txt,"input_tokens":u.get("input_tokens","?"),"output_tokens":u.get("output_tokens","?")})
        print(f"[D] {d} ✓ in={u.get('input_tokens','?')} out={u.get('output_tokens','?')}")

    safe = lambda s: (s or "").replace("\t"," ").replace("\r"," ").replace("\n"," ¶ ")
    with (SON / "targets_sonnet.tsv").open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t"); w.writerow(["target","synthesis","input_tokens","output_tokens"])
        for r in ts: w.writerow([r["target"], safe(r["synthesis"]), r["input_tokens"], r["output_tokens"]])
    with (SON / "drugs_sonnet.tsv").open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t"); w.writerow(["drug","synthesis","input_tokens","output_tokens"])
        for r in ds: w.writerow([r["drug"], safe(r["synthesis"]), r["input_tokens"], r["output_tokens"]])

    int_sum = lambda rs,k: sum(r[k] for r in rs if isinstance(r[k],int))
    idx = {"built_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
           "model": SONNET_MODEL,
           "n_targets": len(ts), "n_drugs": len(ds),
           "input_tokens_total": int_sum(ts+ds,"input_tokens"),
           "output_tokens_total": int_sum(ts+ds,"output_tokens"),
           "elapsed_seconds": round(time.time()-t0,1)}
    (SON / "v15_index.json").write_text(json.dumps(idx, indent=2, ensure_ascii=False))
    print("[done]", json.dumps(idx, indent=2))

if __name__ == "__main__": main()
