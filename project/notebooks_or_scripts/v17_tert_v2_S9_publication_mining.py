"""S9 — PubMed-driven publication mining.

Uses NCBI eUtils to find recent thyroid TERT promoter papers, fetches PMC OA
full text where available, and scans for TCGA barcodes paired with TERT
mutation status.
"""
from __future__ import annotations

import asyncio
import json
import re
from pathlib import Path

import aiohttp

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from v17_tert_v2_common import (
    PARSED,
    RAW,
    SourceResult,
    configure_logger,
    extract_tcga_barcodes,
    fetch_url,
    make_session_kwargs,
    write_attempt_log,
    write_result,
)

LOG = configure_logger("S9_pubmed")
SOURCE_ID = "S9_publication_mining"

ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
PMC_OA = "https://www.ncbi.nlm.nih.gov/pmc/articles"


QUERIES = [
    "TERT promoter AND papillary thyroid AND TCGA",
    "TERT promoter AND thyroid AND (table OR supplementary)",
    "TERT C228T AND thyroid",
    "TERT C250T AND thyroid",
    "TERT promoter mutation thyroid carcinoma 2020:2026[dp]",
    "TERT promoter thyroid review 2022:2026[dp]",
]


async def search_query(session: aiohttp.ClientSession, q: str, retmax: int = 30) -> list[str]:
    url = f"{ESEARCH}?db=pubmed&term={q.replace(' ', '+')}&retmax={retmax}&retmode=json"
    att = await fetch_url(session, url, save_path=RAW / f"{SOURCE_ID}_search_{re.sub(r'[^a-z0-9]', '_', q.lower())[:60]}.json", timeout=30)
    if not att.ok or not att.saved_to:
        return []
    try:
        payload = json.loads(Path(att.saved_to).read_text())
        return payload.get("esearchresult", {}).get("idlist", [])
    except Exception:  # noqa: BLE001
        return []


async def get_pmc_id(session: aiohttp.ClientSession, pmid: str) -> str | None:
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&db=pmc&id={pmid}&retmode=json"
    try:
        async with session.get(url, timeout=aiohttp.ClientTimeout(total=20)) as r:
            if r.status != 200:
                return None
            d = await r.json(content_type=None)
            for ls in d.get("linksets", []):
                for link in ls.get("linksetdbs", []):
                    if link.get("dbto") == "pmc":
                        ids = link.get("links", [])
                        if ids:
                            return f"PMC{ids[0]}"
    except Exception:  # noqa: BLE001
        return None
    return None


async def run() -> SourceResult:
    result = SourceResult(source_id=SOURCE_ID, label="PubMed publication mining")
    all_pmids: set[str] = set()
    async with aiohttp.ClientSession(**make_session_kwargs()) as session:
        for q in QUERIES:
            ids = await search_query(session, q)
            LOG.info("query=%r -> %d ids", q, len(ids))
            all_pmids.update(ids)
        LOG.info("Total unique pmids: %d", len(all_pmids))

        # Map to PMC ids (only OA full text accessible)
        pmid_to_pmc = {}
        sem = asyncio.Semaphore(10)

        async def map_one(pmid):
            async with sem:
                pmc = await get_pmc_id(session, pmid)
                if pmc:
                    pmid_to_pmc[pmid] = pmc

        await asyncio.gather(*[map_one(p) for p in list(all_pmids)[:200]])
        LOG.info("PMID -> PMC mapped: %d", len(pmid_to_pmc))

        # Fetch full text XML for each PMC ID (OA subset)
        async def fetch_pmc_xml(pmid, pmcid):
            url = f"{EFETCH}?db=pmc&id={pmcid.replace('PMC', '')}&rettype=full&retmode=xml"
            return await fetch_url(session, url, save_path=RAW / f"{SOURCE_ID}_{pmcid}.xml", timeout=60)

        attempts = await asyncio.gather(
            *[fetch_pmc_xml(p, pmc) for p, pmc in pmid_to_pmc.items()],
            return_exceptions=True,
        )
        for a in attempts:
            if isinstance(a, Exception):
                continue
            result.attempts.append(a)

    # Scan recovered XMLs for TCGA barcodes + TERT context
    rows = []
    paper_summaries = []
    for xml_path in RAW.glob(f"{SOURCE_ID}_PMC*.xml"):
        try:
            text = xml_path.read_text(errors="ignore")
            if "TERT" not in text:
                continue
            # Strip XML tags for keyword scan
            stripped = re.sub(r"<[^>]+>", " ", text)
            barcodes = extract_tcga_barcodes(stripped)
            promoter_hits = re.findall(r"C228T|C250T|chr5\D+1295\d{3}|-124[CG]>T|-146[CG]>T", stripped, flags=re.I)
            if not (barcodes or promoter_hits):
                continue
            paper_summaries.append({
                "pmcid": xml_path.stem.replace(f"{SOURCE_ID}_", ""),
                "n_barcodes": len(barcodes),
                "n_promoter_keyword": len(promoter_hits),
                "first_barcodes": ",".join(barcodes[:5]),
            })
            for bc in barcodes:
                # Look for TERT/promoter within 300 char window
                for m in re.finditer(re.escape(bc), stripped):
                    ctx = stripped[max(0, m.start() - 250): m.end() + 250]
                    if re.search(r"TERT|C228T|C250T|promoter", ctx, flags=re.I):
                        rows.append({
                            "pmcid": xml_path.stem.replace(f"{SOURCE_ID}_", ""),
                            "tcga_barcode": bc,
                            "context": ctx[:500].replace("\n", " "),
                            "has_C228T": bool(re.search(r"C228T", ctx, flags=re.I)),
                            "has_C250T": bool(re.search(r"C250T", ctx, flags=re.I)),
                        })
        except Exception as e:  # noqa: BLE001
            LOG.warning("xml parse %s: %s", xml_path.name, e)

    import pandas as pd
    if paper_summaries:
        out = PARSED / f"{SOURCE_ID}_paper_summaries.tsv"
        pd.DataFrame(paper_summaries).to_csv(out, sep="\t", index=False)
        result.output_files.append(str(out))
    if rows:
        out = PARSED / f"{SOURCE_ID}_tcga_matched_tert.tsv"
        pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
        result.output_files.append(str(out))
        result.n_records = len(rows)
        result.n_tcga_matched = len({r["tcga_barcode"] for r in rows})
        result.success = True

    result.notes.append(f"Total PMIDs searched: {len(all_pmids)}")
    result.notes.append("Heuristic scan; manual triage required for confidence-grading")
    write_attempt_log(SOURCE_ID, result.attempts)
    write_result(result)
    return result


if __name__ == "__main__":
    asyncio.run(run())
