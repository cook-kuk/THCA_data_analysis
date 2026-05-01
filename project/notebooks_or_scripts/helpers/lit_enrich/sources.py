"""
12 free academic-data sources unified behind dataclasses.

Each client uses requests + the polite-pool email kukshomr@gmail.com.
All public methods either return a list[dict] or a dict; failures
return [] / {} so downstream tasks degrade gracefully.
"""
from __future__ import annotations

import json
import time
import urllib.parse
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from typing import Any

import requests

EMAIL = "kukshomr@gmail.com"
UA = f"THCA-lit-enrich/1.0 (mailto:{EMAIL})"
HDR = {"User-Agent": UA, "Accept": "application/json"}


def _get(url: str, params: dict | None = None, headers: dict | None = None,
         timeout: int = 30, retries: int = 2) -> requests.Response | None:
    h = {**HDR, **(headers or {})}
    for attempt in range(retries + 1):
        try:
            r = requests.get(url, params=params, headers=h, timeout=timeout)
            if r.status_code == 200:
                return r
            if r.status_code in (429, 503) and attempt < retries:
                time.sleep(2 ** attempt)
                continue
            return r
        except requests.RequestException:
            if attempt < retries:
                time.sleep(2 ** attempt)
                continue
            return None
    return None


# ─────────────────────────────────────────────────────────────────────────────
# 1. OpenAlex — 250M works, abstracts since 2022, citations graph
#    https://api.openalex.org/
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class OpenAlex:
    base: str = "https://api.openalex.org"

    def search(self, query: str, limit: int = 20, year_from: int | None = None,
               year_to: int | None = None) -> list[dict]:
        params = {"search": query, "per-page": min(limit, 50), "mailto": EMAIL}
        filters = []
        if year_from:
            filters.append(f"from_publication_date:{year_from}-01-01")
        if year_to:
            filters.append(f"to_publication_date:{year_to}-12-31")
        if filters:
            params["filter"] = ",".join(filters)
        r = _get(f"{self.base}/works", params)
        if not r or r.status_code != 200:
            return []
        return r.json().get("results", [])[:limit]

    def get_work(self, doi: str | None = None, openalex_id: str | None = None) -> dict:
        if doi:
            url = f"{self.base}/works/doi:{doi}"
        elif openalex_id:
            wid = openalex_id.split("/")[-1]
            url = f"{self.base}/works/{wid}"
        else:
            return {}
        r = _get(url, {"mailto": EMAIL})
        return r.json() if r and r.status_code == 200 else {}

    def related_works(self, openalex_id: str, limit: int = 20) -> list[dict]:
        w = self.get_work(openalex_id=openalex_id)
        rel = w.get("related_works", []) or []
        out = []
        for rid in rel[:limit]:
            doc = self.get_work(openalex_id=rid)
            if doc:
                out.append(doc)
            time.sleep(0.1)
        return out


# ─────────────────────────────────────────────────────────────────────────────
# 2. CrossRef — 145M works, DOI authority, references graph
#    https://api.crossref.org/
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class CrossRef:
    base: str = "https://api.crossref.org"

    def lookup(self, doi: str) -> dict:
        r = _get(f"{self.base}/works/{doi}", {"mailto": EMAIL})
        if not r or r.status_code != 200:
            return {}
        return r.json().get("message", {})

    def search(self, query: str, limit: int = 20, year_from: int | None = None) -> list[dict]:
        params: dict[str, Any] = {"query": query, "rows": limit, "mailto": EMAIL}
        if year_from:
            params["filter"] = f"from-pub-date:{year_from}-01-01"
        r = _get(f"{self.base}/works", params)
        if not r or r.status_code != 200:
            return []
        return r.json().get("message", {}).get("items", [])


# ─────────────────────────────────────────────────────────────────────────────
# 3. PubMed E-utilities — 35M biomedical
#    https://eutils.ncbi.nlm.nih.gov/entrez/eutils/
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class PubMed:
    base: str = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    def search(self, query: str, limit: int = 20) -> list[str]:
        params = {"db": "pubmed", "term": query, "retmax": limit,
                  "retmode": "json", "tool": "THCA-lit-enrich", "email": EMAIL}
        r = _get(f"{self.base}/esearch.fcgi", params)
        if not r or r.status_code != 200:
            return []
        return r.json().get("esearchresult", {}).get("idlist", [])

    def summary(self, pmids: list[str]) -> dict[str, dict]:
        if not pmids:
            return {}
        params = {"db": "pubmed", "id": ",".join(pmids), "retmode": "json",
                  "tool": "THCA-lit-enrich", "email": EMAIL}
        r = _get(f"{self.base}/esummary.fcgi", params)
        if not r or r.status_code != 200:
            return {}
        result = r.json().get("result", {})
        return {pid: result.get(pid, {}) for pid in pmids if pid in result}

    def search_with_summary(self, query: str, limit: int = 20) -> list[dict]:
        ids = self.search(query, limit)
        time.sleep(0.34)  # NCBI: ≤3 req/sec
        summaries = self.summary(ids)
        return [summaries[pid] for pid in ids if pid in summaries]


# ─────────────────────────────────────────────────────────────────────────────
# 4. Europe PMC — 40M, full-text search for OA
#    https://europepmc.org/RestfulWebService
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class EuropePMC:
    base: str = "https://www.ebi.ac.uk/europepmc/webservices/rest"

    def search(self, query: str, limit: int = 20, full_text_only: bool = False) -> list[dict]:
        q = query
        if full_text_only:
            q = f"({q}) AND HAS_FT:Y"
        params = {"query": q, "format": "json", "pageSize": limit, "resultType": "core"}
        r = _get(f"{self.base}/search", params)
        if not r or r.status_code != 200:
            return []
        return r.json().get("resultList", {}).get("result", [])

    def get_full_text_url(self, pmcid: str) -> str | None:
        if not pmcid:
            return None
        if not pmcid.startswith("PMC"):
            pmcid = "PMC" + pmcid
        return f"https://europepmc.org/article/PMC/{pmcid[3:]}"


# ─────────────────────────────────────────────────────────────────────────────
# 5. Semantic Scholar — 200M, influential citations
#    https://api.semanticscholar.org/graph/v1
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class SemanticScholar:
    base: str = "https://api.semanticscholar.org/graph/v1"

    def search(self, query: str, limit: int = 20) -> list[dict]:
        fields = "title,authors,year,venue,abstract,influentialCitationCount,citationCount,externalIds,openAccessPdf"
        params = {"query": query, "limit": min(limit, 100), "fields": fields}
        r = _get(f"{self.base}/paper/search", params)
        if not r or r.status_code != 200:
            return []
        return r.json().get("data", [])

    def get_paper(self, paper_id: str) -> dict:
        fields = "title,authors,year,venue,abstract,influentialCitationCount,citationCount,externalIds,references,citations"
        r = _get(f"{self.base}/paper/{paper_id}", {"fields": fields})
        return r.json() if r and r.status_code == 200 else {}

    def influential_citations(self, paper_id: str, limit: int = 20) -> list[dict]:
        fields = "intent,isInfluential,citingPaper.title,citingPaper.year,citingPaper.authors,citingPaper.externalIds"
        params = {"fields": fields, "limit": limit}
        r = _get(f"{self.base}/paper/{paper_id}/citations", params)
        if not r or r.status_code != 200:
            return []
        rows = r.json().get("data", [])
        return [x for x in rows if x.get("isInfluential")]


# ─────────────────────────────────────────────────────────────────────────────
# 6. arXiv — 2.4M preprints
#    http://export.arxiv.org/api/query
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class Arxiv:
    base: str = "http://export.arxiv.org/api/query"
    ns: dict = field(default_factory=lambda: {
        "atom": "http://www.w3.org/2005/Atom",
        "arxiv": "http://arxiv.org/schemas/atom",
    })

    def search(self, query: str, limit: int = 20) -> list[dict]:
        params = {"search_query": f"all:{query}", "start": 0, "max_results": limit,
                  "sortBy": "relevance"}
        r = _get(self.base, params, headers={"Accept": "application/atom+xml"})
        if not r or r.status_code != 200:
            return []
        try:
            root = ET.fromstring(r.text)
        except ET.ParseError:
            return []
        out = []
        for entry in root.findall("atom:entry", self.ns):
            out.append({
                "id": (entry.findtext("atom:id", "", self.ns) or "").strip(),
                "title": (entry.findtext("atom:title", "", self.ns) or "").strip(),
                "summary": (entry.findtext("atom:summary", "", self.ns) or "").strip(),
                "published": entry.findtext("atom:published", "", self.ns),
                "authors": [a.findtext("atom:name", "", self.ns)
                            for a in entry.findall("atom:author", self.ns)],
            })
        return out


# ─────────────────────────────────────────────────────────────────────────────
# 7. bioRxiv / medRxiv — preprints (very 2024-2026 active for thyroid/HLA)
#    https://api.biorxiv.org/details/{server}/...
#    They don't have free-text search; we use crossref_relations endpoint or
#    Europe PMC's PPR (preprint) source.
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class BioRxiv:
    """Search via Europe PMC since bioRxiv has no full-text search API."""
    epmc: EuropePMC = field(default_factory=EuropePMC)

    def search(self, query: str, server: str = "biorxiv", limit: int = 20) -> list[dict]:
        # PPR = preprint source in Europe PMC; bioRxiv DOIs start with 10.1101/
        src_filter = "PPR"
        q = f'({query}) AND SRC:{src_filter}'
        if server == "biorxiv":
            q += ' AND DOI:"10.1101/*"'
        return self.epmc.search(q, limit=limit)


# ─────────────────────────────────────────────────────────────────────────────
# 8. ICite (NIH) — RCR (relative citation ratio)
#    https://icite.od.nih.gov/api/pubs?pmids=...
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class ICite:
    base: str = "https://icite.od.nih.gov/api/pubs"

    def get(self, pmids: list[str]) -> dict[str, dict]:
        if not pmids:
            return {}
        out: dict[str, dict] = {}
        # batch in 200s
        for i in range(0, len(pmids), 200):
            batch = pmids[i:i + 200]
            r = _get(self.base, {"pmids": ",".join(batch)})
            if not r or r.status_code != 200:
                continue
            for row in r.json().get("data", []):
                out[str(row.get("pmid"))] = row
            time.sleep(0.3)
        return out


# ─────────────────────────────────────────────────────────────────────────────
# 9. AFND (Allele Frequency Net Database) — HLA allele freq
#    No JSON API; we scrape the public results page.
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class AFND:
    base: str = "http://www.allelefrequencies.net/hla6006a.asp"

    def allele_frequency(self, allele: str, population_filter: str | None = None,
                         limit: int = 50) -> list[dict]:
        """
        allele: e.g. 'DPB1*05:01', 'DRB1*04:05'
        population_filter: case-insensitive substring match on population name
                           (e.g. 'Korea' matches 'South Korea pop2'). Applied
                           client-side after AFND returns the global table for
                           the allele.
        """
        # AFND requires the full hidden-form param set to actually return
        # the results table. `hla_show=>` (URL-encoded `%3E`) opts into
        # frequencies > 0.0001 (i.e. all real entries).
        params = {
            "hla_locus_type": "Classical",
            "hla_locus": allele.split("*")[0],
            "hla_allele1": allele,
            "hla_allele2": "",
            "hla_selection": "",
            "hla_pop_selection": "",
            "hla_population": "",
            "hla_country": "",
            "hla_dataset": "",
            "hla_region": "",
            "hla_ethnic": "",
            "hla_study": "",
            "hla_order": "order_1",
            "hla_sample_size_pattern": "bigger_equal_than",
            "hla_sample_size": "",
            "hla_sample_year_pattern": "equal",
            "hla_sample_year": "",
            "hla_loci": "",
            "standard": "g",
            "hla_show": ">",
        }
        r = _get(self.base, params, headers={"Accept": "text/html"}, timeout=45)
        if not r or r.status_code != 200:
            return []
        import re
        rows: list[dict] = []
        # Result row layout (12 cells, after stripping HTML / nbsp):
        #   [0]=rank  [1]=allele  [2]=''  [3]=population  [4]=''
        #   [5]=allele_freq  [6]=''  [7]=sample_size  [8]='See' …
        for row_match in re.finditer(r"<tr[^>]*>(.*?)</tr>", r.text, re.S | re.I):
            cells = re.findall(r"<td[^>]*>(.*?)</td>", row_match.group(1), re.S | re.I)
            cells = [re.sub(r"<[^>]+>", "", c).replace("&nbsp;", " ").strip() for c in cells]
            if len(cells) < 8:
                continue
            freq_idx = None
            for i, c in enumerate(cells):
                if re.match(r"^0?\.\d{3,}$", c):
                    freq_idx = i
                    break
            if freq_idx is None or freq_idx < 3:
                continue
            pop = cells[freq_idx - 2]
            if not pop or pop.lower() in ("population", ""):
                continue
            if population_filter and population_filter.lower() not in pop.lower():
                continue
            rows.append({
                "population": pop,
                "allele_freq": cells[freq_idx],
                "sample_size": cells[freq_idx + 2] if len(cells) > freq_idx + 2 else "",
                "allele": cells[freq_idx - 4] if freq_idx >= 4 else allele,
            })
            if len(rows) >= limit:
                break
        return rows


# ─────────────────────────────────────────────────────────────────────────────
# 10. ClinicalTrials.gov v2 API
#     https://clinicaltrials.gov/api/v2/studies?query.term=...
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class ClinicalTrials:
    base: str = "https://clinicaltrials.gov/api/v2/studies"

    def search(self, query: str, limit: int = 30, status: str | None = "RECRUITING|ACTIVE_NOT_RECRUITING") -> list[dict]:
        params: dict[str, Any] = {
            "query.term": query,
            "pageSize": min(limit, 100),
            "format": "json",
        }
        if status:
            params["filter.overallStatus"] = status
        r = _get(self.base, params)
        if not r or r.status_code != 200:
            return []
        return r.json().get("studies", [])[:limit]


# ─────────────────────────────────────────────────────────────────────────────
# 11. Unpaywall — OA discovery
#     https://api.unpaywall.org/v2/{doi}?email=...
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class Unpaywall:
    base: str = "https://api.unpaywall.org/v2"

    def lookup(self, doi: str) -> dict:
        r = _get(f"{self.base}/{doi}", {"email": EMAIL})
        return r.json() if r and r.status_code == 200 else {}


# ─────────────────────────────────────────────────────────────────────────────
# 12. CORE — 200M OA papers (requires free API key in env CORE_API_KEY)
#     https://api.core.ac.uk/v3/search/works
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class CORE:
    import os
    base: str = "https://api.core.ac.uk/v3"

    def _api_key(self) -> str | None:
        import os
        return os.environ.get("CORE_API_KEY")

    def search(self, query: str, limit: int = 20) -> list[dict]:
        key = self._api_key()
        if not key:
            return []  # gracefully skip if no key
        headers = {"Authorization": f"Bearer {key}", "User-Agent": UA}
        params = {"q": query, "limit": min(limit, 100)}
        r = _get(f"{self.base}/search/works", params, headers)
        if not r or r.status_code != 200:
            return []
        return r.json().get("results", [])
