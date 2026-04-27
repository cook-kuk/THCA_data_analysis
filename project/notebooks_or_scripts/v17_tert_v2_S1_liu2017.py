"""S1 — Liu et al. JCO 2017 (PMID 27979994) TERT promoter Sanger supplement.

Tries 8 alternative URLs for the supplementary file. If a PDF is recovered,
attempts table extraction with pdfplumber to find TERT mutation status per
TCGA barcode.
"""
from __future__ import annotations

import asyncio
import io
import json
import re
from pathlib import Path

import aiohttp

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from v17_tert_v2_common import (
    Attempt,
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

LOG = configure_logger("S1_liu2017")
SOURCE_ID = "S1_liu2017"

CANDIDATE_URLS = [
    # Direct ASCO suppl files (multiple guesses for filename schemes)
    "https://ascopubs.org/doi/suppl/10.1200/JCO.2016.71.5654/suppl_file/CR_DS_JCO.2016.71.5654-1.pdf",
    "https://ascopubs.org/doi/suppl/10.1200/JCO.2016.71.5654/suppl_file/JCO.2016.71.5654-1.xlsx",
    "https://ascopubs.org/doi/suppl/10.1200/JCO.2016.71.5654/suppl_file/JCO.2016.71.5654-1.pdf",
    "https://ascopubs.org/doi/suppl/10.1200/JCO.2016.71.5654/suppl_file/protocol_JCO.2016.71.5654.pdf",
    "https://ascopubs.org/doi/abs/10.1200/JCO.2016.71.5654",
    # PubMed Central
    "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5455594/",
    "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5455594/bin/NIHMS876716-supplement-Supplemental_Table_1.xlsx",
    "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5455594/bin/NIHMS876716-supplement-Supplemental_Tables.xlsx",
    "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5455594/bin/NIHMS876716-supplement-1.xlsx",
    "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5455594/bin/NIHMS876716-supplement-1.pdf",
    "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5455594/bin/supplemental_table1.xlsx",
    "https://europepmc.org/article/MED/27979994",
    # ResearchGate (often blocked but try anyway)
    "https://www.researchgate.net/publication/311936147",
    # Google Scholar fallback (just landing page; won't have file but useful as audit trail)
    "https://scholar.google.com/scholar?q=Liu+2017+TERT+promoter+TCGA+thyroid+JCO",
]


async def run() -> SourceResult:
    result = SourceResult(source_id=SOURCE_ID, label="Liu 2017 JCO supplement")
    LOG.info("Attempting %d URLs", len(CANDIDATE_URLS))
    async with aiohttp.ClientSession(**make_session_kwargs()) as session:
        coros = []
        for url in CANDIDATE_URLS:
            ext = ".pdf" if ".pdf" in url else (".xlsx" if ".xlsx" in url else ".html")
            safe = re.sub(r"[^A-Za-z0-9]", "_", url)[-80:]
            save_path = RAW / f"{SOURCE_ID}_{safe}{ext}"
            coros.append(fetch_url(session, url, save_path=save_path, timeout=45))
        attempts = await asyncio.gather(*coros, return_exceptions=False)
    result.attempts = attempts
    n_ok = sum(1 for a in attempts if a.ok)
    LOG.info("OK fetches: %d / %d", n_ok, len(attempts))

    # Try to parse any recovered PDF or XLSX for TERT promoter status
    parsed_records = []
    barcode_hits = set()
    for a in attempts:
        if not a.ok or not a.saved_to:
            continue
        path = Path(a.saved_to)
        full = path if path.is_absolute() else (Path.cwd() / path)
        if not full.exists():
            continue
        if a.bytes_ < 1024:
            continue  # likely error page
        text = ""
        try:
            if full.suffix.lower() == ".pdf":
                import pdfplumber

                with pdfplumber.open(full) as pdf:
                    for page in pdf.pages:
                        t = page.extract_text() or ""
                        text += t + "\n"
            elif full.suffix.lower() == ".xlsx":
                import pandas as pd

                xl = pd.ExcelFile(full)
                for sn in xl.sheet_names:
                    df = xl.parse(sn).astype(str)
                    text += " ".join(df.values.flatten().tolist()) + "\n"
            else:
                text = full.read_text(errors="ignore")[:200_000]
        except Exception as e:  # noqa: BLE001
            LOG.warning("parse failure %s: %s", full.name, e)
            continue

        codes = extract_tcga_barcodes(text)
        if codes:
            barcode_hits.update(codes)
            for bc in codes:
                # heuristic: scan ±150 chars around barcode for TERT/C228T/C250T markers
                for m in re.finditer(re.escape(bc), text, flags=re.I):
                    ctx = text[max(0, m.start() - 200): m.end() + 200]
                    flag_promoter = bool(
                        re.search(r"C228T|C250T|chr5:1295|promoter", ctx, flags=re.I)
                    )
                    flag_tert = bool(re.search(r"\bTERT\b", ctx))
                    if flag_tert or flag_promoter:
                        parsed_records.append(
                            {
                                "tcga_barcode": bc,
                                "source_file": full.name,
                                "context": ctx[:400].replace("\n", " "),
                                "promoter_keyword": flag_promoter,
                                "tert_keyword": flag_tert,
                            }
                        )
    result.notes.append(f"Unique TCGA barcodes mentioned: {len(barcode_hits)}")
    result.n_records = len(parsed_records)
    if parsed_records:
        import pandas as pd

        out = PARSED / f"{SOURCE_ID}_tert_status.tsv"
        pd.DataFrame(parsed_records).to_csv(out, sep="\t", index=False)
        result.output_files.append(str(out))
        result.success = True

    write_attempt_log(SOURCE_ID, attempts)
    write_result(result)
    return result


if __name__ == "__main__":
    asyncio.run(run())
