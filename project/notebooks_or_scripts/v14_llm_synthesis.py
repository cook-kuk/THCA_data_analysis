#!/usr/bin/env python3
"""v14 — LLM-based evidence synthesis on top of v12 literature pull.

For each of the 8 druggable targets, take the top abstracts already pulled in
v12_literature.py (results/v12_literature/pubmed_evidence/{GENE}_pubmed_hits.tsv),
construct a context, and ask Claude (haiku) for a structured 3-section verdict:
  • biology summary (BRAF-like context)
  • druggability status
  • clinical translation status

Same idea for the 5 v7 drugs.

API key sourced from cook-forge/.env (user-authorised for this run).
"""
import os, re, csv, json, time, sys
from pathlib import Path
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError

ROOT = Path("/opt/thyroid-dash/project")
RES  = ROOT / "results/v12_literature"
PUB  = RES / "pubmed_evidence"
LLM  = RES / "llm"
LLM.mkdir(parents=True, exist_ok=True)

# --- load API key from cook-forge/.env ---
ENV_PATH = Path("/home/seungho/cook-forge/.env")
KEY = None
HAIKU = "claude-haiku-4-5-20251001"
if ENV_PATH.exists():
    for line in ENV_PATH.read_text().splitlines():
        line = line.strip()
        if line.startswith("ANTHROPIC_API_KEY="):
            KEY = line.split("=", 1)[1].strip().strip('"').strip("'")
        if line.startswith("ANTHROPIC_HAIKU_MODEL="):
            HAIKU = line.split("=", 1)[1].strip().strip('"').strip("'")
if not KEY:
    print("ERROR: ANTHROPIC_API_KEY not found in", ENV_PATH); sys.exit(1)
print(f"[v14] using model={HAIKU}")

# --- minimal anthropic SDK client over urllib ---
def claude(system, user, max_tokens=400, retries=2):
    body = json.dumps({
        "model": HAIKU,
        "max_tokens": max_tokens,
        "system": system,
        "messages": [{"role":"user","content":user}],
    }).encode()
    req = Request("https://api.anthropic.com/v1/messages",
        data=body,
        headers={
            "x-api-key": KEY,
            "anthropic-version": "2023-06-01",
            "content-type": "application/json",
        })
    last = None
    for attempt in range(retries+1):
        try:
            with urlopen(req, timeout=60) as r:
                d = json.loads(r.read().decode())
                # Sum text blocks
                txt = "".join(b.get("text","") for b in d.get("content",[]))
                usage = d.get("usage",{})
                return txt.strip(), usage
        except HTTPError as e:
            err_body = e.read().decode(errors="ignore")[:300]
            last = f"HTTP {e.code}: {err_body}"
            if e.code in (429, 529, 500, 502, 503): time.sleep(2 + attempt*3)
            else: break
        except URLError as e:
            last = f"URL: {e}"; time.sleep(2 + attempt*2)
    raise RuntimeError(f"Claude call failed after retries: {last}")

# --- build context from per-gene PubMed hits ---
def gene_context(gene, max_n=8, max_chars=2400):
    p = PUB / f"{gene}_pubmed_hits.tsv"
    if not p.exists(): return ""
    rows = []
    with p.open() as f:
        r = csv.DictReader(f, delimiter="\t")
        rows = list(r)
    # Score: prioritise records that have direction or clinical or BRAF-mentioning queries
    def score(row):
        s = 0
        if row.get("direction"): s += 2
        if row.get("clinical")=="1": s += 3
        if row.get("druggable")=="1": s += 1
        if row.get("outcome")=="1": s += 1
        if "braf" in (row.get("query_types","") or "").lower(): s += 2
        try: s += min(int(row.get("year","0") or 0) - 2014, 12)
        except: pass
        return -s
    rows.sort(key=score)
    out = []
    used = 0
    for row in rows[:max_n]:
        snippet = (f"PMID {row['pmid']} ({row.get('year','')}, {row.get('journal','')[:60]}): "
                   f"{row.get('title','')[:240]}").strip()
        if used + len(snippet) > max_chars: break
        out.append(snippet)
        used += len(snippet)
    return "\n".join(out)

TARGETS = ["CYP1B1","TACSTD2","TMPRSS4","PLEKHA6","LDLR","GABRB2","B3GNT3","PTPRE"]
DRUGS   = ["Cannabidiol","Sacituzumab govitecan","Tesevatinib","Dacomitinib","Tepotinib"]

SYS_TARGET = (
"You are a careful biomedical research synthesist. The user is preparing a Bioinformatics paper "
"on a thyroid-cancer (BRAF-like vs RAS-like) computational pipeline. They have already pulled "
"PubMed abstracts for a candidate druggable target. Your job is to read the abstracts and write "
"three short paragraphs labelled exactly:\n"
"BIOLOGY: 1-3 sentences on what the gene does and what's known in BRAF-like thyroid (or related cancer) context.\n"
"DRUGGABILITY: 1-2 sentences on inhibitors / antibodies / ADCs / repurposed drugs targeting it.\n"
"CLINICAL: 1-2 sentences on clinical-trial or FDA evidence (be honest if sparse).\n"
"Total length <= 130 words. Do NOT invent citations or numbers. If evidence is missing for a section, say so."
)

SYS_DRUG = (
"You are a careful biomedical research synthesist. The user is evaluating a drug repurposing "
"candidate for BRAF-like thyroid cancer. From the supplied PubMed abstract titles, write three "
"short paragraphs labelled exactly:\n"
"MECHANISM: 1-2 sentences — what the drug does and the molecular target.\n"
"THYROID EVIDENCE: 1-2 sentences — direct thyroid-cancer literature (be honest if sparse).\n"
"CROSS-CANCER: 1-2 sentences — broader oncology evidence and any approved indication.\n"
"Total <= 120 words. No invented citations or numbers."
)

def synthesise_target(g):
    ctx = gene_context(g, max_n=10, max_chars=2400)
    if not ctx:
        return None, None
    user = f"GENE: {g}\n\nABSTRACTS:\n{ctx}\n\nWrite the three labelled paragraphs."
    return claude(SYS_TARGET, user, max_tokens=320)

def drug_context(drug, max_n=6):
    """Pull a few abstracts directly: hit pubmed via E-utilities for the drug."""
    import urllib.parse, urllib.request
    q = urllib.parse.urlencode({
        "db":"pubmed","term":f'"{drug}"[tiab] AND cancer[tiab]',
        "retmode":"json","retmax":max_n,
        "tool":"thyrai-v14","email":"kukshomr@gmail.com"
    })
    try:
        with urllib.request.urlopen(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?{q}", timeout=20) as r:
            j = json.loads(r.read().decode())
        ids = j.get("esearchresult",{}).get("idlist",[])
    except Exception:
        return ""
    if not ids: return ""
    q2 = urllib.parse.urlencode({
        "db":"pubmed","id":",".join(ids),"retmode":"xml",
        "tool":"thyrai-v14","email":"kukshomr@gmail.com"
    })
    try:
        with urllib.request.urlopen(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?{q2}", timeout=30) as r:
            xml = r.read().decode()
    except Exception:
        return ""
    titles = re.findall(r"<ArticleTitle[^>]*>(.*?)</ArticleTitle>", xml, re.S)
    years  = re.findall(r"<PubDate>.*?<Year>(\d{4})</Year>", xml, re.S)
    out = []
    for i, t in enumerate(titles[:max_n]):
        y = years[i] if i < len(years) else ""
        out.append(f"({y}) {re.sub(r'<[^>]+>','',t).strip()[:240]}")
    return "\n".join(out)

def synthesise_drug(d):
    ctx = drug_context(d)
    if not ctx:
        ctx = "(no recent PubMed abstracts retrieved)"
    user = f"DRUG: {d}\n\nABSTRACT TITLES:\n{ctx}\n\nWrite the three labelled paragraphs."
    return claude(SYS_DRUG, user, max_tokens=280)

# --- run ---
def main():
    t0 = time.time()
    target_rows = []
    for g in TARGETS:
        try:
            txt, usage = synthesise_target(g)
        except Exception as e:
            print(f"[err] {g}: {e}"); continue
        if not txt:
            print(f"[skip] {g}: no abstracts")
            continue
        (LLM / f"target_{g}.md").write_text(f"# {g} — LLM evidence synthesis\n\n{txt}\n", encoding="utf-8")
        target_rows.append({"target":g, "synthesis":txt, "input_tokens":usage.get("input_tokens","?"),
                            "output_tokens":usage.get("output_tokens","?")})
        print(f"[target] {g} ✓ in={usage.get('input_tokens','?')} out={usage.get('output_tokens','?')}")

    drug_rows = []
    for d in DRUGS:
        try:
            txt, usage = synthesise_drug(d)
        except Exception as e:
            print(f"[err] {d}: {e}"); continue
        slug = re.sub(r"\W+", "_", d).strip("_").lower()
        (LLM / f"drug_{slug}.md").write_text(f"# {d} — LLM evidence synthesis\n\n{txt}\n", encoding="utf-8")
        drug_rows.append({"drug":d, "synthesis":txt, "input_tokens":usage.get("input_tokens","?"),
                          "output_tokens":usage.get("output_tokens","?")})
        print(f"[drug]   {d} ✓ in={usage.get('input_tokens','?')} out={usage.get('output_tokens','?')}")

    # --- combined TSVs (synthesis as a single multi-line cell, escaped) ---
    def safe(s): return (s or "").replace("\t"," ").replace("\r"," ").replace("\n"," ¶ ")
    with (LLM / "targets_llm_synthesis.tsv").open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["target","synthesis","input_tokens","output_tokens"])
        for r in target_rows: w.writerow([r["target"], safe(r["synthesis"]), r["input_tokens"], r["output_tokens"]])
    with (LLM / "drugs_llm_synthesis.tsv").open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["drug","synthesis","input_tokens","output_tokens"])
        for r in drug_rows: w.writerow([r["drug"], safe(r["synthesis"]), r["input_tokens"], r["output_tokens"]])

    # --- index ---
    idx = {
        "built_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "model": HAIKU,
        "n_targets": len(target_rows),
        "n_drugs": len(drug_rows),
        "input_tokens_total": sum((r["input_tokens"] or 0) for r in target_rows+drug_rows if isinstance(r["input_tokens"], int)),
        "output_tokens_total": sum((r["output_tokens"] or 0) for r in target_rows+drug_rows if isinstance(r["output_tokens"], int)),
        "elapsed_seconds": round(time.time()-t0, 1),
    }
    (LLM / "v14_index.json").write_text(json.dumps(idx, indent=2, ensure_ascii=False))
    print("[done]", json.dumps(idx, indent=2))

if __name__ == "__main__":
    main()
