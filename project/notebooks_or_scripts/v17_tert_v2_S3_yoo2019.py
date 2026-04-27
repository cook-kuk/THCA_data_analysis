"""S3 — Yoo et al. Nature Communications 2019 (PMID 31388009) Korean PTC TERT.

Korean cohort with TERT promoter Sanger. If recovered, opens door to a
cross-cohort validation that doesn't depend on TCGA-THCA WES coverage.
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
    SourceResult,
    configure_logger,
    fetch_url,
    make_session_kwargs,
    write_attempt_log,
    write_result,
)

LOG = configure_logger("S3_yoo2019")
SOURCE_ID = "S3_yoo2019"

CANDIDATE_URLS = [
    # Nature Communications article + supplements
    "https://www.nature.com/articles/s41467-019-11551-9",
    "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-019-11551-9/MediaObjects/41467_2019_11551_MOESM1_ESM.pdf",
    "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-019-11551-9/MediaObjects/41467_2019_11551_MOESM2_ESM.xlsx",
    "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-019-11551-9/MediaObjects/41467_2019_11551_MOESM3_ESM.xlsx",
    "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-019-11551-9/MediaObjects/41467_2019_11551_MOESM4_ESM.xlsx",
    "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-019-11551-9/MediaObjects/41467_2019_11551_MOESM5_ESM.xlsx",
    "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-019-11551-9/MediaObjects/41467_2019_11551_MOESM6_ESM.xlsx",
    # PMC
    "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6692329/",
    "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6692329/bin/41467_2019_11551_MOESM2_ESM.xlsx",
    "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6692329/bin/41467_2019_11551_MOESM3_ESM.xlsx",
]


async def run() -> SourceResult:
    result = SourceResult(source_id=SOURCE_ID, label="Yoo 2019 Korean cohort")
    async with aiohttp.ClientSession(**make_session_kwargs()) as session:
        coros = []
        for url in CANDIDATE_URLS:
            ext = ".xlsx" if ".xlsx" in url else (".pdf" if ".pdf" in url else ".html")
            safe = re.sub(r"[^A-Za-z0-9]", "_", url)[-80:]
            coros.append(fetch_url(session, url, save_path=RAW / f"{SOURCE_ID}_{safe}{ext}", timeout=60))
        attempts = await asyncio.gather(*coros)
    result.attempts = attempts
    LOG.info("OK: %d / %d", sum(1 for a in attempts if a.ok), len(attempts))

    rows = []
    for a in attempts:
        if not a.ok or not a.saved_to or a.bytes_ < 1024:
            continue
        path = Path(a.saved_to)
        full = path if path.is_absolute() else (Path.cwd() / path)
        if not full.exists():
            continue
        try:
            if full.suffix.lower() == ".xlsx":
                import pandas as pd

                xl = pd.ExcelFile(full)
                for sn in xl.sheet_names:
                    df = xl.parse(sn)
                    df_str = df.astype(str)
                    text = " ".join(df_str.values.flatten().tolist())
                    if "TERT" not in text.upper():
                        continue
                    LOG.info("%s sheet %s mentions TERT", full.name, sn)
                    # Save sheet for manual review
                    out = PARSED / f"{SOURCE_ID}_{full.stem}_{sn}.tsv"
                    df.to_csv(out, sep="\t", index=False)
                    result.output_files.append(str(out))
                    # Look for TERT mutation rows / annotations
                    for i, r in df.iterrows():
                        row_text = " ".join(str(v) for v in r.values)
                        if re.search(r"TERT", row_text):
                            promoter = bool(re.search(r"C228T|C250T|promoter|chr5.*1295|-124|-146", row_text, flags=re.I))
                            rows.append(
                                {
                                    "source_file": full.name,
                                    "sheet": sn,
                                    "row_idx": i,
                                    "promoter_keyword": promoter,
                                    "row_data": row_text[:400],
                                }
                            )
            elif full.suffix.lower() == ".pdf":
                import pdfplumber

                with pdfplumber.open(full) as pdf:
                    for pg_idx, page in enumerate(pdf.pages):
                        t = page.extract_text() or ""
                        if "TERT" in t:
                            LOG.info("PDF %s page %d mentions TERT", full.name, pg_idx)
        except Exception as e:  # noqa: BLE001
            LOG.warning("parse err %s: %s", full.name, e)

    if rows:
        import pandas as pd

        out = PARSED / f"{SOURCE_ID}_korean_tert_rows.tsv"
        pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
        result.output_files.append(str(out))
        result.n_records = len(rows)
        result.success = True
    result.notes.append(
        "Korean cohort patient IDs do not map to TCGA barcodes; usable as external cohort only"
    )
    write_attempt_log(SOURCE_ID, attempts)
    write_result(result)
    return result


if __name__ == "__main__":
    asyncio.run(run())
