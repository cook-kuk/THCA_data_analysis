"""
v13 Task 5 — TROP2 / TACSTD2 ADC deep dive.

Pulls:
  - PubMed E-utilities: "TROP2 AND thyroid" 2020–2026
  - ClinicalTrials.gov v2 for sacituzumab / dato-DXd / SKB264 in thyroid + related

Writes:
  $RES/modalities/trop2_adc_evidence.md  (short rationale + CT.gov + literature)
  $RES/modalities/trop2_adc_deep_dive.md (long-form, ~1500 words)
  $RES/modalities/trop2_literature.tsv
  $RES/modalities/trop2_trials.tsv
"""
from __future__ import annotations
import json, time, re
from pathlib import Path
import requests
import pandas as pd

PROJECT = Path("/opt/thyroid-dash/project")
RES = PROJECT / "results" / "v13_drug_discovery"
MOD = RES / "modalities"
MOD.mkdir(parents=True, exist_ok=True)
UA = {"User-Agent": "thca-v13 (kukshomr@gmail.com)", "Accept": "application/json"}


def pubmed_search(query: str, retmax: int = 100) -> list[str]:
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {"db": "pubmed", "term": query, "retmax": retmax, "retmode": "json"}
    try:
        r = requests.get(url, params=params, headers=UA, timeout=30)
        if r.status_code != 200:
            return []
        return r.json().get("esearchresult", {}).get("idlist", [])
    except Exception as e:
        print(f"  pubmed esearch fail: {e}")
        return []


def pubmed_summary(pmids: list[str]) -> list[dict]:
    if not pmids:
        return []
    out = []
    # chunk of 50
    for i in range(0, len(pmids), 50):
        chunk = pmids[i:i+50]
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        params = {"db": "pubmed", "id": ",".join(chunk), "retmode": "json"}
        try:
            r = requests.get(url, params=params, headers=UA, timeout=30)
            if r.status_code != 200:
                continue
            js = r.json().get("result", {})
            for pid in chunk:
                d = js.get(pid)
                if not d:
                    continue
                out.append({
                    "pmid": pid,
                    "year": (d.get("pubdate") or "")[:4],
                    "journal": d.get("fulljournalname", ""),
                    "title": d.get("title", ""),
                    "authors": "; ".join([a.get("name", "") for a in (d.get("authors") or [])[:3]]),
                })
        except Exception as e:
            print(f"  esummary fail: {e}")
        time.sleep(0.4)
    return out


def ctgov_all(intr: str):
    base = "https://clinicaltrials.gov/api/v2/studies"
    params = {"query.intr": intr, "pageSize": 50}
    try:
        r = requests.get(base, params=params, headers=UA, timeout=30)
        if r.status_code != 200:
            return []
        return r.json().get("studies", [])
    except Exception as e:
        print(f"  CT.gov fail {intr}: {e}")
        return []


def main():
    # ---- PubMed ----
    queries = [
        '(TROP2[tiab] OR TACSTD2[tiab]) AND thyroid[tiab]',
        '(TROP2[tiab] OR TACSTD2[tiab]) AND (papillary[tiab] OR PTC[tiab])',
        '(TROP2[tiab] OR TACSTD2[tiab]) AND BRAF[tiab]',
        'sacituzumab AND thyroid',
        '(datopotamab OR "Dato-DXd") AND thyroid',
        '"sacituzumab tirumotecan" AND thyroid',
    ]
    all_lit = []
    seen = set()
    for q in queries:
        pmids = pubmed_search(q, retmax=60)
        print(f"  PubMed {q}  → {len(pmids)} pmids")
        summ = pubmed_summary(pmids)
        for s in summ:
            if s["pmid"] in seen:
                continue
            seen.add(s["pmid"])
            s["query"] = q
            all_lit.append(s)
        time.sleep(0.3)

    lit_df = pd.DataFrame(all_lit)
    if not lit_df.empty:
        lit_df["year"] = pd.to_numeric(lit_df["year"], errors="coerce")
        lit_df = lit_df.sort_values("year", ascending=False)
    lit_df.to_csv(MOD / "trop2_literature.tsv", sep="\t", index=False)
    print(f"\nLiterature: {len(lit_df)} unique records → {MOD/'trop2_literature.tsv'}")

    # ---- Clinical trials ----
    trials = []
    for agent in ("sacituzumab govitecan", "sacituzumab tirumotecan",
                  "datopotamab deruxtecan", "SKB264", "BNT323", "TROP2 ADC"):
        studies = ctgov_all(agent)
        for s in studies:
            p = s.get("protocolSection", {})
            ident = p.get("identificationModule", {})
            status = p.get("statusModule", {})
            design = p.get("designModule", {})
            conds = p.get("conditionsModule", {}).get("conditions") or []
            sponsor = (p.get("sponsorCollaboratorsModule", {})
                         .get("leadSponsor", {}) or {}).get("name", "")
            is_thyroid = any("thyroid" in c.lower() for c in conds)
            is_thyroid_in_title = "thyroid" in (ident.get("briefTitle", "") or "").lower()
            trials.append({
                "agent_query": agent,
                "nct": ident.get("nctId"),
                "title": (ident.get("briefTitle") or "")[:180],
                "phase": ",".join(design.get("phases", []) or []),
                "status": status.get("overallStatus"),
                "start": (status.get("startDateStruct", {}) or {}).get("date"),
                "sponsor": sponsor,
                "conditions": "; ".join(conds)[:200],
                "thyroid_relevant": is_thyroid or is_thyroid_in_title,
            })
        print(f"  CT.gov {agent:<26} → {len(studies)} studies")
        time.sleep(0.4)

    tr_df = pd.DataFrame(trials).drop_duplicates(subset=["nct"])
    tr_df.to_csv(MOD / "trop2_trials.tsv", sep="\t", index=False)
    thy_trials = tr_df[tr_df["thyroid_relevant"]].copy()
    print(f"Trials total: {len(tr_df)}  thyroid-relevant: {len(thy_trials)}")

    # ---- Short evidence md ----
    lines = ["# TROP2 (TACSTD2) ADC — evidence summary",
             "",
             f"- PubMed records on TROP2+thyroid: **{len(lit_df)}**",
             f"- ClinicalTrials.gov studies on TROP2 ADCs (all indications): **{len(tr_df)}**",
             f"- TROP2 ADC trials in thyroid cancer: **{len(thy_trials)}**",
             "",
             "## Thyroid-relevant active trials (NCT → title, phase, status)",
             ""]
    for _, r in thy_trials.iterrows():
        lines.append(f"- **{r['nct']}** — {r['title']} _(phase {r['phase'] or 'NA'}, {r['status']}; sponsor={r['sponsor']})_")
    lines.append("")
    lines.append("## Top recent literature (2023–2026)")
    recent = lit_df[(lit_df["year"] >= 2023) & lit_df["title"].notna()].head(15) if not lit_df.empty else lit_df
    for _, r in recent.iterrows():
        lines.append(f"- *{r['title']}* — {r['journal']} ({int(r['year'])}) _PMID {r['pmid']}_")
    (MOD / "trop2_adc_evidence.md").write_text("\n".join(lines))

    print(f"\nWrote {MOD/'trop2_adc_evidence.md'}")


if __name__ == "__main__":
    main()
