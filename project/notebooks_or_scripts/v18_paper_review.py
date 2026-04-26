#!/usr/bin/env python3
"""v18 — LLM-judged paper-readiness review of bioinformatics_additions.tex.
Uses claude-sonnet to play the role of a Bioinformatics journal reviewer
and produce a structured review (3 strengths / 3 concerns / line-edit suggestions).
Output: results/v12_literature/llm/review/v18_review.md + .json."""
import os, json, time
from pathlib import Path
from urllib.request import Request, urlopen
from urllib.error import HTTPError

ROOT = Path("/opt/thyroid-dash/project")
RES  = ROOT / "results/v12_literature"
REV  = RES / "llm/review"
REV.mkdir(parents=True, exist_ok=True)

ENV_PATH = Path("/home/seungho/cook-forge/.env")
KEY = None; SONNET = "claude-sonnet-4-5-20250929"
for line in ENV_PATH.read_text().splitlines():
    if line.startswith("ANTHROPIC_API_KEY="): KEY = line.split("=",1)[1].strip().strip('"').strip("'")
    if line.startswith("ANTHROPIC_SONNET_MODEL="):
        SONNET = line.split("=",1)[1].strip().strip('"').strip("'")
print(f"[v18] using {SONNET}")

def claude(system, user, max_tokens=1500):
    body = json.dumps({"model":SONNET,"max_tokens":max_tokens,"system":system,
        "messages":[{"role":"user","content":user}]}).encode()
    req = Request("https://api.anthropic.com/v1/messages", data=body,
        headers={"x-api-key":KEY,"anthropic-version":"2023-06-01","content-type":"application/json"})
    with urlopen(req, timeout=120) as r:
        d = json.loads(r.read().decode())
        return "".join(b.get("text","") for b in d.get("content",[])).strip(), d.get("usage",{})

tex = (ROOT / "reports/v12/v12_paper_updates/bioinformatics_additions.tex").read_text()
narrative = (ROOT / "reports/v12/v12_biology_narrative.md").read_text()[:4000]

SYS = ("You are an experienced reviewer for the journal Bioinformatics, evaluating two new "
"Results subsections (4.9 and 4.10) added to a manuscript on a thyroid-cancer "
"computational pipeline. The pipeline produces a 2,843-gene replicated biomarker slate, "
"an 8-target druggable shortlist, and a direction-invariant batch-leakage probe (DIAL). "
"You read the .tex source plus a 1,800-word companion narrative. Provide a structured review "
"with these labelled sections (use bold markdown):\n\n"
"VERDICT: <one of: accept / minor_revision / major_revision / reject> + 1 sentence.\n"
"3 STRENGTHS: numbered, one sentence each, specific.\n"
"3 CONCERNS: numbered, one sentence each, with a concrete suggested fix.\n"
"LINE-EDIT SUGGESTIONS: 3-5 bullets pointing to specific phrases or numbers that would benefit from rephrasing or qualification (cite the phrase).\n"
"MISSING TABLE/FIGURE: one suggestion for what the next supplementary should contain.\n"
"OVERALL: 2-3 sentences on the paper's contribution as currently framed.")

user = (f"### bioinformatics_additions.tex\n```latex\n{tex}\n```\n\n"
        f"### Companion narrative (excerpt, first 4000 chars)\n```\n{narrative}\n```\n\n"
        "Write the structured review now.")

def main():
    t0=time.time()
    txt, u = claude(SYS, user, max_tokens=1500)
    out = REV / "v18_review.md"
    out.write_text(f"# v18 — LLM paper-readiness review\n\n_Model: {SONNET}_  \n_Generated: "
                   f"{time.strftime('%Y-%m-%d %H:%M UTC', time.gmtime())}_\n\n{txt}\n", encoding="utf-8")
    idx = {"built_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
           "model": SONNET, "input_tokens": u.get("input_tokens","?"),
           "output_tokens": u.get("output_tokens","?"),
           "elapsed_seconds": round(time.time()-t0, 1)}
    (REV / "v18_index.json").write_text(json.dumps(idx, indent=2, ensure_ascii=False))
    print("[done]", json.dumps(idx, indent=2))
    print("---preview (first 800 chars)---")
    print(txt[:800])

if __name__ == "__main__": main()
