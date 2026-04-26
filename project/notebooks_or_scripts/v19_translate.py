#!/usr/bin/env python3
"""v19 — bilingual (KR ↔ EN) of the v12 narrative key sections.
Use claude-sonnet for high-quality translation (research-paper register).
Output: reports/v12/v12_biology_narrative.kr.md + en.md (sectioned)."""
import os, json, time, re
from pathlib import Path
from urllib.request import Request, urlopen

ROOT = Path("/opt/thyroid-dash/project")
RPT  = ROOT / "reports/v12"

ENV_PATH = Path("/home/seungho/cook-forge/.env")
KEY = None; SONNET = "claude-sonnet-4-5-20250929"
for line in ENV_PATH.read_text().splitlines():
    if line.startswith("ANTHROPIC_API_KEY="): KEY = line.split("=",1)[1].strip().strip('"').strip("'")
    if line.startswith("ANTHROPIC_SONNET_MODEL="):
        SONNET = line.split("=",1)[1].strip().strip('"').strip("'")

def claude(system, user, max_tokens=2000):
    body = json.dumps({"model":SONNET,"max_tokens":max_tokens,"system":system,
        "messages":[{"role":"user","content":user}]}).encode()
    req = Request("https://api.anthropic.com/v1/messages", data=body,
        headers={"x-api-key":KEY,"anthropic-version":"2023-06-01","content-type":"application/json"})
    with urlopen(req, timeout=120) as r:
        d = json.loads(r.read().decode())
        return "".join(b.get("text","") for b in d.get("content",[])).strip(), d.get("usage",{})

src = (RPT / "v12_biology_narrative.md").read_text()

# Pick sections 1, 11 (Introduction, Headline claims) — these are the publication-critical ones.
sections = re.split(r"^## ", src, flags=re.M)
# sections[0] = preamble (frontmatter)
# sections[1] starts with "1. Introduction"
def get_sec(n):
    for s in sections:
        if s.startswith(f"{n}."):
            return "## " + s
    return None

intro = get_sec(1) or ""
# Try multiple candidate section labels: 11. Headline / 10. / Limitations / last section
headline = (get_sec(11) or get_sec(12) or get_sec(10) or get_sec(9) or "")
if not headline.strip():
    # fall back to last labeled section
    candidates = [s for s in sections if re.match(r"^\d+\.", s)]
    if candidates:
        headline = "## " + candidates[-1]

SYS_KR = ("You are a senior biomedical research translator. Translate the supplied English research-paper "
"section into Korean (한국어), preserving:\n"
"• numerical values and units exactly\n"
"• gene symbols, drug names, and assay/method names in their original Roman form\n"
"• statistical notations like 'log10 p < -10', '25/25', 'd=+2.52' verbatim\n"
"• inline citations like (Chakravarty 2011) verbatim\n"
"• markdown structure (headings, tables, bold)\n"
"Use academic Korean register (논문체). Do not paraphrase — translate.")
SYS_EN = ("You are a senior biomedical research translator. Polish/translate the supplied Korean text into "
"English suitable for a Bioinformatics journal Results/Discussion section. Preserve all numbers, gene "
"symbols, and inline citations. Use precise scientific register; avoid casual phrasing.")

def main():
    t0=time.time(); calls=0; tin=0; tout=0
    out_kr = ["# v12 narrative (한국어 번역, sonnet)\n\n_원본 영문: `v12_biology_narrative.md`_\n"]
    for label, sec in (("Introduction", intro), ("Headline claims", headline)):
        if not sec.strip(): continue
        txt, u = claude(SYS_KR, f"### {label}\n\n{sec}", max_tokens=2000)
        calls += 1; tin += u.get("input_tokens",0) or 0; tout += u.get("output_tokens",0) or 0
        out_kr.append(f"\n---\n{txt}\n")
    (RPT / "v12_biology_narrative.kr.md").write_text("\n".join(out_kr), encoding="utf-8")

    # Also re-polish the headline-claims paragraph as a tightened English version (for abstract)
    polished, u = claude(
        ("You are a tight scientific copy-editor. Rewrite the supplied 'headline claims' bullet list "
         "into a single 110-130 word abstract paragraph for the Bioinformatics manuscript. Preserve "
         "every numeric claim verbatim. No bullet points; flowing prose. End with one sentence about "
         "prospective external validation."),
        headline, max_tokens=400)
    calls += 1; tin += u.get("input_tokens",0) or 0; tout += u.get("output_tokens",0) or 0
    (RPT / "v12_paper_updates/abstract_para.txt").write_text(polished + "\n")

    idx = {"built_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
           "model": SONNET, "calls": calls,
           "input_tokens_total": tin, "output_tokens_total": tout,
           "elapsed_seconds": round(time.time()-t0, 1),
           "outputs": ["reports/v12/v12_biology_narrative.kr.md",
                       "reports/v12/v12_paper_updates/abstract_para.txt"]}
    (RPT / "v19_translate_index.json").write_text(json.dumps(idx, indent=2, ensure_ascii=False))
    print("[done]", json.dumps(idx, indent=2))

if __name__ == "__main__": main()
