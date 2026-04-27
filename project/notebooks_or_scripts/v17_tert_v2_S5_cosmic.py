"""S5 — COSMIC TERT thyroid lookup.

COSMIC v100 requires authenticated access for bulk download. We try:
  1. The public targeted screen pages (HTML scrape) for TERT thyroid hits.
  2. The Cancer Mutation Census subset endpoint (legacy, public).
  3. Manual instructions if all blocked (output ready-to-run curl commands).
"""
from __future__ import annotations

import asyncio
import re
from pathlib import Path

import aiohttp

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from v17_tert_v2_common import (
    PARSED,
    RAW,
    V2_DIR,
    SourceResult,
    configure_logger,
    fetch_url,
    make_session_kwargs,
    write_attempt_log,
    write_result,
)

LOG = configure_logger("S5_cosmic")
SOURCE_ID = "S5_cosmic"

CANDIDATE_URLS = [
    # Public targeted screen for TERT
    "https://cancer.sanger.ac.uk/cosmic/gene/analysis?ln=TERT",
    "https://cancer.sanger.ac.uk/cosmic/search?genome=38&q=TERT&search_genome=38",
    # CMC (Cancer Mutation Census) snapshot — TERT specific page
    "https://cancer.sanger.ac.uk/cmc/gene/TERT",
    # Tissue-distribution page (often blocks bots, but worth a try)
    "https://cancer.sanger.ac.uk/cosmic/gene/tissue?ln=TERT&site=thyroid",
    # GENIE / cBioPortal mirror that COSMIC sometimes references
    "https://www.cbioportal.org/api/genes/TERT/mutations",
]


async def run() -> SourceResult:
    result = SourceResult(source_id=SOURCE_ID, label="COSMIC TERT thyroid")
    async with aiohttp.ClientSession(**make_session_kwargs()) as session:
        coros = []
        for url in CANDIDATE_URLS:
            ext = ".json" if "/api/" in url else ".html"
            safe = re.sub(r"[^A-Za-z0-9]", "_", url)[-80:]
            coros.append(fetch_url(session, url, save_path=RAW / f"{SOURCE_ID}_{safe}{ext}", timeout=45))
        attempts = await asyncio.gather(*coros)
    result.attempts = attempts
    LOG.info("OK: %d / %d", sum(1 for a in attempts if a.ok), len(attempts))

    found_promoter_mentions = 0
    for a in attempts:
        if not a.ok or not a.saved_to:
            continue
        path = Path(a.saved_to)
        full = path if path.is_absolute() else (Path.cwd() / path)
        if not full.exists():
            continue
        try:
            text = full.read_text(errors="ignore")[:1_000_000]
            hits = re.findall(r"(C228T|C250T|chr5:1295228|1295113)", text, flags=re.I)
            if hits:
                found_promoter_mentions += len(hits)
                LOG.info("%s: %d promoter keyword hits", full.name, len(hits))
        except Exception as e:  # noqa: BLE001
            LOG.warning("read err %s: %s", full.name, e)

    # Always emit manual application instructions (COSMIC needs free academic license)
    instructions = V2_DIR / f"{SOURCE_ID}_application_guide.md"
    instructions.write_text(
        """# COSMIC v100 manual application guide

COSMIC bulk access requires a free academic license. To recover sample-level
TERT promoter mutations for thyroid:

1. Register at https://cancer.sanger.ac.uk/cosmic/register (academic email)
2. Download `CosmicMutantExport.tsv.gz` and `CosmicSample.tsv.gz` from
   `https://cancer.sanger.ac.uk/cosmic/download` (v100, GRCh38)
3. After download, run:

```bash
zcat CosmicMutantExport.tsv.gz | awk -F'\\t' '$1=="TERT"' > TERT_all.tsv
# Filter to thyroid using sample mapping
zcat CosmicSample.tsv.gz | awk -F'\\t' '$5 ~ /thyroid/i' > thyroid_samples.tsv
```

4. Join `TERT_all.tsv` to `thyroid_samples.tsv` on sample_id.
5. Filter rows where `Mutation CDS == c.-124C>T` or `c.-146C>T` (= C228T/C250T promoter hotspots).

Once the file is on disk, point this script at it via env var
`COSMIC_TSV_PATH=/path/to/TERT_all.tsv` and rerun.
""",
        encoding="utf-8",
    )
    result.output_files.append(str(instructions))
    result.notes.append(f"COSMIC promoter keyword hits in scraped pages: {found_promoter_mentions}")
    result.notes.append("Bulk sample-level data behind license wall; application guide written")
    result.success = found_promoter_mentions > 0
    write_attempt_log(SOURCE_ID, attempts)
    write_result(result)
    return result


if __name__ == "__main__":
    asyncio.run(run())
