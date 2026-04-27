"""S8 — bioRxiv / medRxiv preprint search for thyroid TERT TCGA papers.

Uses the public bioRxiv API and tries to extract supplementary URLs and
GitHub/Zenodo data repositories from the abstracts/landing pages.
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
    fetch_url,
    make_session_kwargs,
    write_attempt_log,
    write_result,
)

LOG = configure_logger("S8_preprint")
SOURCE_ID = "S8_preprint"

BIORXIV_DETAILS = "https://api.biorxiv.org/details/biorxiv/2024-01-01/2026-04-27/0"
MEDRXIV_DETAILS = "https://api.biorxiv.org/details/medrxiv/2024-01-01/2026-04-27/0"

QUERY_TERMS = [
    re.compile(r"\bTERT\b.*\bthyroid\b", re.I),
    re.compile(r"\bthyroid\b.*\bTERT\b", re.I),
    re.compile(r"papillary thyroid.*TERT", re.I),
    re.compile(r"TCGA.*thyroid.*TERT", re.I),
]


async def crawl(session: aiohttp.ClientSession, base_url: str, max_pages: int = 5) -> list[dict]:
    """Walk biorxiv pages until we have enough or hit cap."""
    out = []
    cursor = base_url
    for _ in range(max_pages):
        try:
            async with session.get(cursor, timeout=aiohttp.ClientTimeout(total=60)) as resp:
                if resp.status != 200:
                    break
                data = await resp.json(content_type=None)
                msgs = data.get("collection", []) or data.get("messages", []) or []
                # api returns "collection" with paper records
                papers = data.get("collection", [])
                out.extend(papers)
                # Pagination via "messages" cursor
                msg = data.get("messages", [{}])[0]
                next_cursor = msg.get("cursor")
                if not next_cursor:
                    break
                cursor = base_url.rsplit("/", 1)[0] + f"/{next_cursor}"
        except Exception as e:  # noqa: BLE001
            LOG.warning("crawl error: %s", e)
            break
    return out


async def run() -> SourceResult:
    result = SourceResult(source_id=SOURCE_ID, label="bioRxiv/medRxiv preprints")
    matched_papers = []

    async with aiohttp.ClientSession(**make_session_kwargs()) as session:
        # Use the bioRxiv search endpoint (faster than crawling all 2024-2026)
        # Try direct query API first
        for term, server in [
            ("TERT+thyroid", "biorxiv"),
            ("thyroid+TERT+promoter", "biorxiv"),
            ("TERT+thyroid", "medrxiv"),
        ]:
            url = f"https://api.biorxiv.org/details/{server}/{term}/0"
            att = await fetch_url(session, url, save_path=RAW / f"{SOURCE_ID}_{server}_{term}.json", timeout=60)
            result.attempts.append(att)
            if att.ok and att.saved_to:
                try:
                    payload = json.loads(Path(att.saved_to).read_text())
                    papers = payload.get("collection", [])
                    LOG.info("server=%s term=%s -> %d papers", server, term, len(papers))
                    for p in papers:
                        text = (
                            (p.get("title", "") + " " + p.get("abstract", "") + " " + p.get("authors", ""))
                        )
                        if any(q.search(text) for q in QUERY_TERMS):
                            matched_papers.append({**p, "matched_via": f"{server}:{term}"})
                except Exception as e:  # noqa: BLE001
                    LOG.warning("biorxiv parse %s", e)

        # Also try a direct Europe PMC search (covers preprints + published)
        epmc_url = (
            "https://www.ebi.ac.uk/europepmc/webservices/rest/search?"
            "query=%22TERT+promoter%22+AND+%22thyroid%22&format=json&pageSize=50"
        )
        att = await fetch_url(session, epmc_url, save_path=RAW / f"{SOURCE_ID}_europepmc.json", timeout=60)
        result.attempts.append(att)
        if att.ok and att.saved_to:
            try:
                payload = json.loads(Path(att.saved_to).read_text())
                hits = payload.get("resultList", {}).get("result", [])
                LOG.info("EuropePMC hits: %d", len(hits))
                for h in hits[:50]:
                    matched_papers.append(
                        {
                            "doi": h.get("doi"),
                            "title": h.get("title"),
                            "abstract": h.get("abstractText", "")[:500],
                            "year": h.get("pubYear"),
                            "pmid": h.get("pmid"),
                            "pmcid": h.get("pmcid"),
                            "matched_via": "europepmc",
                        }
                    )
            except Exception as e:  # noqa: BLE001
                LOG.warning("epmc parse %s", e)

    if matched_papers:
        import pandas as pd

        # Deduplicate by DOI/title
        seen = set()
        deduped = []
        for p in matched_papers:
            key = p.get("doi") or p.get("title") or json.dumps(p)[:100]
            if key in seen:
                continue
            seen.add(key)
            deduped.append(p)
        out = PARSED / f"{SOURCE_ID}_papers.tsv"
        pd.DataFrame(deduped).to_csv(out, sep="\t", index=False)
        result.output_files.append(str(out))
        result.n_records = len(deduped)
        result.success = True
    result.notes.append("Preprint papers only — sample-level data extraction would require manual triage")
    write_attempt_log(SOURCE_ID, result.attempts)
    write_result(result)
    return result


if __name__ == "__main__":
    asyncio.run(run())
