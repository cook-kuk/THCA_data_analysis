"""
v14 Prior-art check — PubMed E-utilities + ClinicalTrials.gov v2.

Runs 8 novelty queries, fetches titles + abstracts (first 500 chars),
and a deterministic relevance-scoring pass based on keyword matches.

Outputs:
  $RES/v14_priorart/all_queries_results.tsv
  $RES/v14_priorart/high_relevance_papers.md
  $RES/v14_priorart/trial_status.md
  $RES/v14_priorart/priorart_summary.json
"""
from __future__ import annotations
import json, time, re
from pathlib import Path
import requests
import pandas as pd
import xml.etree.ElementTree as ET

PROJECT = Path("/opt/thyroid-dash/project")
RES = PROJECT / "results" / "v14_ccle" / "v14_priorart"
RES.mkdir(parents=True, exist_ok=True)
UA = {"User-Agent": "thca-v14-priorart (kukshomr@gmail.com)"}

QUERIES = {
    "A": '(BRAF AND TROP2 AND thyroid)',
    "B": '(BRAF AND TACSTD2 AND thyroid)',
    "C": '(V600E AND TROP2)',
    "D": '(thyroid AND sacituzumab)',
    "E": '(thyroid AND TROP2 AND ADC)',
    "F": '(BRAF AND irinotecan AND cancer)',
    "G": '(BRAF AND topotecan AND cancer)',
    "H": '(SN-38 AND mutation AND selectivity)',
}

# Keywords for automatic relevance scoring — these trigger HIGH if multiple co-occur
HIGH_KEYWORDS = [
    re.compile(r"\bBRAF\b", re.I),
    re.compile(r"\b(TROP2|TACSTD2)\b", re.I),
    re.compile(r"\bthyroid\b", re.I),
    re.compile(r"\b(sacituzumab|govitecan|datopotamab|ADC|antibody[- ]drug)\b", re.I),
    re.compile(r"\b(SN-?38|irinotecan|topotecan|topoisomerase)\b", re.I),
    re.compile(r"\b(stratif|biomarker|subtype|V600E)\b", re.I),
]


def esearch(term: str, retmax: int = 20) -> list[str]:
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {"db": "pubmed", "term": term, "retmax": retmax, "retmode": "json",
              "sort": "relevance"}
    r = requests.get(url, params=params, headers=UA, timeout=30)
    if r.status_code != 200:
        return []
    return r.json().get("esearchresult", {}).get("idlist", [])


def efetch_abstracts(pmids: list[str]) -> dict[str, dict]:
    """Return dict pmid -> {title, abstract, journal, year, authors}."""
    out = {}
    if not pmids:
        return out
    # batch of 50
    for i in range(0, len(pmids), 50):
        chunk = pmids[i:i+50]
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        params = {"db": "pubmed", "id": ",".join(chunk), "rettype": "abstract", "retmode": "xml"}
        try:
            r = requests.get(url, params=params, headers=UA, timeout=45)
            if r.status_code != 200:
                continue
            root = ET.fromstring(r.content)
            for art in root.findall(".//PubmedArticle"):
                pmid = art.findtext(".//PMID")
                if not pmid:
                    continue
                title = " ".join((art.findtext(".//ArticleTitle") or "").split())
                journal = art.findtext(".//Journal/Title") or art.findtext(".//ISOAbbreviation") or ""
                year = art.findtext(".//PubDate/Year") or art.findtext(".//PubMedPubDate/Year") or ""
                abs_parts = [t.text or "" for t in art.findall(".//Abstract/AbstractText")]
                abstract = " ".join(" ".join(p.split()) for p in abs_parts if p)
                authors = []
                for au in art.findall(".//Author")[:3]:
                    last = au.findtext("LastName") or ""
                    init = au.findtext("Initials") or ""
                    if last:
                        authors.append(f"{last} {init}")
                out[pmid] = {
                    "pmid": pmid,
                    "title": title,
                    "journal": journal,
                    "year": year,
                    "authors": "; ".join(authors),
                    "abstract": abstract,
                }
        except Exception as e:
            print(f"  efetch fail: {e}")
        time.sleep(0.4)
    return out


def score_relevance(row: dict) -> str:
    text = f"{row.get('title','')} {row.get('abstract','')}"
    hits = sum(1 for k in HIGH_KEYWORDS if k.search(text))
    has_braf = bool(re.search(r"\bBRAF\b", text, re.I))
    has_trop = bool(re.search(r"\b(TROP2|TACSTD2)\b", text, re.I))
    has_thyroid = bool(re.search(r"\bthyroid\b", text, re.I))
    has_adc = bool(re.search(r"\b(sacituzumab|govitecan|datopotamab|ADC|antibody[- ]drug)\b", text, re.I))
    has_pay = bool(re.search(r"\b(SN-?38|irinotecan|topotecan|topoisomerase)\b", text, re.I))
    # HIGH = our exact thesis: BRAF + (TROP2/TACSTD2) + thyroid OR (BRAF + SN38 class + selectivity)
    if (has_braf and has_trop and has_thyroid):
        return "HIGH"
    if (has_braf and has_trop and has_adc):
        return "HIGH"
    if (has_braf and has_pay and re.search(r"\b(selectivity|biomarker|stratif|sensitiv|response)\b", text, re.I)):
        return "HIGH"
    if (has_thyroid and has_trop and has_adc):
        return "HIGH"
    if hits >= 3:
        return "MEDIUM"
    if hits >= 2:
        return "MEDIUM"
    return "LOW"


def main():
    all_rows = []
    for qid, term in QUERIES.items():
        print(f"[{qid}] {term}")
        pmids = esearch(term, retmax=20)
        print(f"    esearch → {len(pmids)} PMIDs")
        abstracts = efetch_abstracts(pmids)
        for pmid, rec in abstracts.items():
            rec["query_id"] = qid
            rec["query_term"] = term
            rec["relevance"] = score_relevance(rec)
            rec["abstract_500"] = rec["abstract"][:500]
            all_rows.append(rec)
        time.sleep(0.5)

    df = pd.DataFrame(all_rows)
    print(f"\nTotal rows: {len(df)}  unique PMIDs: {df['pmid'].nunique() if len(df) else 0}")
    if len(df):
        # keep full abstract in separate col, short in TSV
        df_out = df[["query_id", "query_term", "pmid", "year", "journal", "authors",
                      "title", "relevance", "abstract_500"]].copy()
        df_out.to_csv(RES / "all_queries_results.tsv", sep="\t", index=False)
        print(f"Wrote {RES/'all_queries_results.tsv'}")

        # per-query counts by relevance
        print("\nPer-query relevance counts:")
        print(df.groupby(["query_id", "relevance"]).size().unstack(fill_value=0))

        # Unique HIGH papers
        high = df[df["relevance"] == "HIGH"].drop_duplicates(subset=["pmid"])
        high = high.sort_values("year", ascending=False)
        print(f"\nHIGH relevance unique: {len(high)}")

        lines = ["# v14 Prior-art — HIGH-relevance papers", ""]
        lines.append(f"_Generated 2026-04-25 from {len(df)} PubMed records across 8 queries._")
        lines.append("")
        lines.append("HIGH = our exact thesis (BRAF + TROP2/TACSTD2 + thyroid) or (BRAF + SN-38-class + selectivity/biomarker).")
        lines.append("")
        for _, r in high.iterrows():
            lines.append(f"### PMID {r['pmid']} — {r['year']} · {r['journal']}")
            lines.append(f"**{r['title']}**")
            lines.append(f"_authors: {r['authors']}_  ·  queries: `{r['query_id']}` (`{r['query_term']}`)")
            lines.append("")
            lines.append(f"> {r['abstract'][:800]}...")
            lines.append("")
        (RES / "high_relevance_papers.md").write_text("\n".join(lines))
        print(f"Wrote {RES/'high_relevance_papers.md'}")

    # ---- CT.gov for NCT06235216 + NCT07521670 ----
    trials_md = ["# Trial-design check — NCT06235216 & NCT07521670",
                 "",
                 "_For each trial: status, BRAF-mention in eligibility, endpoints, sponsor._",
                 ""]
    trial_findings = {}
    for nct in ("NCT06235216", "NCT07521670"):
        try:
            r = requests.get(f"https://clinicaltrials.gov/api/v2/studies/{nct}", headers=UA, timeout=30)
            if r.status_code != 200:
                trials_md.append(f"## {nct}\n\nHTTP {r.status_code}\n")
                continue
            js = r.json()
            p = js.get("protocolSection", {})
            ident = p.get("identificationModule", {})
            status = p.get("statusModule", {})
            design = p.get("designModule", {})
            elig = p.get("eligibilityModule", {})
            outcomes = p.get("outcomesModule", {})
            sponsor = (p.get("sponsorCollaboratorsModule", {})
                         .get("leadSponsor", {}) or {}).get("name", "")
            arms = p.get("armsInterventionsModule", {})
            crit = elig.get("eligibilityCriteria", "") or ""
            crit_lower = crit.lower()
            braf_mention = bool(re.search(r"braf|v600", crit_lower))
            trop_mention = bool(re.search(r"trop[- ]?2|tacstd2", crit_lower))
            primary_out = [o.get("measure","") for o in (outcomes.get("primaryOutcomes") or [])]
            sec_out = [o.get("measure","") for o in (outcomes.get("secondaryOutcomes") or [])]
            arm_list = [a.get("label","") for a in (arms.get("armGroups") or [])]
            interv_list = [i.get("name","") for i in (arms.get("interventions") or [])]

            trial_findings[nct] = {
                "title": ident.get("briefTitle", ""),
                "status": status.get("overallStatus", ""),
                "start": (status.get("startDateStruct", {}) or {}).get("date", ""),
                "primary_completion": (status.get("primaryCompletionDateStruct", {}) or {}).get("date", ""),
                "sponsor": sponsor,
                "phase": ",".join(design.get("phases", []) or []),
                "enrollment": (design.get("enrollmentInfo", {}) or {}).get("count"),
                "braf_mentioned_in_eligibility": braf_mention,
                "trop2_mentioned_in_eligibility": trop_mention,
                "primary_outcomes": primary_out,
                "secondary_outcomes": sec_out[:3],
                "arms": arm_list,
                "interventions": interv_list,
                "inclusion_excerpt": crit[:1500],
            }
            trials_md.append(f"## {nct} — {ident.get('briefTitle','')}")
            trials_md.append("")
            trials_md.append(f"- **Status**: {status.get('overallStatus','')}")
            trials_md.append(f"- **Phase**: {','.join(design.get('phases', []) or []) or 'NA'}")
            trials_md.append(f"- **Sponsor**: {sponsor}")
            trials_md.append(f"- **Start / primary completion**: {trial_findings[nct]['start']} → {trial_findings[nct]['primary_completion']}")
            trials_md.append(f"- **Enrollment target**: {trial_findings[nct]['enrollment']}")
            trials_md.append(f"- **BRAF mentioned in eligibility?** {'**YES**' if braf_mention else 'no'}")
            trials_md.append(f"- **TROP2 mentioned in eligibility?** {'**YES**' if trop_mention else 'no'}")
            trials_md.append(f"- **Arms**: {' | '.join(arm_list) or '(not listed)'}")
            trials_md.append(f"- **Interventions**: {' | '.join(interv_list) or '(not listed)'}")
            trials_md.append(f"- **Primary outcomes**: {'; '.join(primary_out) or '(not listed)'}")
            trials_md.append(f"- **Secondary (first 3)**: {'; '.join(sec_out[:3]) or '(not listed)'}")
            trials_md.append("")
            trials_md.append(f"**Inclusion excerpt (first 1500 chars):**\n\n> {crit[:1500]}")
            trials_md.append("")
        except Exception as e:
            trials_md.append(f"## {nct}\n\nfetch failed: {e}\n")
        time.sleep(0.4)
    (RES / "trial_status.md").write_text("\n".join(trials_md))
    (RES / "priorart_summary.json").write_text(json.dumps({
        "queries": list(QUERIES.keys()),
        "n_pubmed_records": len(df) if 'df' in locals() and len(df) else 0,
        "n_unique_pmid": int(df['pmid'].nunique()) if 'df' in locals() and len(df) else 0,
        "n_high": int((df['relevance'] == 'HIGH').sum()) if 'df' in locals() and len(df) else 0,
        "trials": trial_findings,
    }, indent=2, default=str))
    print(f"\nWrote {RES/'trial_status.md'}")
    print(f"Wrote {RES/'priorart_summary.json'}")


if __name__ == "__main__":
    main()
