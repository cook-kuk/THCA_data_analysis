"""S2 — Landa et al. JCI 2016 (PMID 27613709) PDTC/ATC TERT promoter Sanger.

Note: corrected paper reference is JCI 2016 ("Genomic and Transcriptomic
Hallmarks..."), not Cancer Cell. cBioPortal study id: thca_msk_2016 / pdtc_msk_2017.
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
    is_promoter_position,
    make_session_kwargs,
    write_attempt_log,
    write_result,
)

LOG = configure_logger("S2_landa2016")
SOURCE_ID = "S2_landa2016"

CANDIDATE_URLS = [
    # JCI / Cell direct
    "https://www.jci.org/articles/view/85271",
    "https://www.jci.org/articles/view/85271/sd/1",
    "https://www.jci.org/articles/view/85271/sd/2",
    "https://dm5migu4zj3pb.cloudfront.net/manuscripts/85000/85271/JCI85271sd.pdf",
    # PMC
    "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5113236/",
    "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5113236/bin/JCI85271sd1.xlsx",
    "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5113236/bin/JCI85271sd2.xlsx",
    # cBioPortal datahub (thca_msk_2016 — Landa et al. PDTC / ATC)
    "https://github.com/cBioPortal/datahub/raw/master/public/thca_mskcc_2016/data_mutations.txt",
    "https://github.com/cBioPortal/datahub/raw/master/public/thca_mskcc_2016/data_mutations_extended.txt",
    "https://github.com/cBioPortal/datahub/raw/master/public/thca_mskcc_2016/meta_mutations_extended.txt",
    "https://github.com/cBioPortal/datahub/raw/master/public/thca_mskcc_2016/data_clinical_sample.txt",
    # Alternative PDTC IDs
    "https://github.com/cBioPortal/datahub/raw/master/public/pdtc_mskcc_2017/data_mutations.txt",
    "https://github.com/cBioPortal/datahub/raw/master/public/atc_mskcc_2017/data_mutations.txt",
]


async def run() -> SourceResult:
    result = SourceResult(source_id=SOURCE_ID, label="Landa 2016 PDTC/ATC")
    async with aiohttp.ClientSession(**make_session_kwargs()) as session:
        coros = []
        for url in CANDIDATE_URLS:
            ext = (".xlsx" if ".xlsx" in url else (".pdf" if ".pdf" in url else ".txt"))
            safe = re.sub(r"[^A-Za-z0-9]", "_", url)[-80:]
            coros.append(fetch_url(session, url, save_path=RAW / f"{SOURCE_ID}_{safe}{ext}", timeout=60))
        attempts = await asyncio.gather(*coros, return_exceptions=False)
    result.attempts = attempts
    LOG.info("OK fetches: %d / %d", sum(1 for a in attempts if a.ok), len(attempts))

    # Look for TERT mutation rows in any recovered MAF
    promoter_records = []
    other_tert_records = []
    for a in attempts:
        if not a.ok or not a.saved_to or a.bytes_ < 100:
            continue
        path = Path(a.saved_to)
        full = path if path.is_absolute() else (Path.cwd() / path)
        if not full.exists():
            continue
        try:
            if full.suffix == ".txt":
                import pandas as pd

                # MAF format
                df = pd.read_csv(full, sep="\t", comment="#", low_memory=False)
                if "Hugo_Symbol" not in df.columns:
                    continue
                tert = df[df["Hugo_Symbol"].astype(str).str.upper() == "TERT"]
                LOG.info("%s TERT rows: %d", full.name, len(tert))
                for _, r in tert.iterrows():
                    chrom = r.get("Chromosome", "")
                    pos = r.get("Start_Position", r.get("Start_position", 0))
                    sample = r.get("Tumor_Sample_Barcode", "")
                    rec = {
                        "sample_barcode": sample,
                        "variant_class": r.get("Variant_Classification", ""),
                        "chrom": chrom,
                        "start": pos,
                        "ref": r.get("Reference_Allele", ""),
                        "alt": r.get("Tumor_Seq_Allele2", ""),
                        "source_file": full.name,
                    }
                    if is_promoter_position(chrom, pos):
                        promoter_records.append(rec)
                    else:
                        other_tert_records.append(rec)
            elif full.suffix == ".xlsx":
                import pandas as pd

                try:
                    xl = pd.ExcelFile(full)
                    for sn in xl.sheet_names:
                        df = xl.parse(sn).astype(str)
                        text = " ".join(df.values.flatten().tolist())
                        if "TERT" in text.upper():
                            LOG.info("xlsx %s sheet %s mentions TERT", full.name, sn)
                except Exception as e:  # noqa: BLE001
                    LOG.warning("xlsx parse err: %s", e)
        except Exception as e:  # noqa: BLE001
            LOG.warning("parse err %s: %s", full.name, e)

    if promoter_records:
        import pandas as pd

        p = PARSED / f"{SOURCE_ID}_pdtc_atc_tert_promoter.tsv"
        pd.DataFrame(promoter_records).to_csv(p, sep="\t", index=False)
        result.output_files.append(str(p))
        result.n_promoter_mutations = len(promoter_records)
        result.success = True
    if other_tert_records:
        import pandas as pd

        p2 = PARSED / f"{SOURCE_ID}_pdtc_atc_tert_other.tsv"
        pd.DataFrame(other_tert_records).to_csv(p2, sep="\t", index=False)
        result.output_files.append(str(p2))
    result.n_records = len(promoter_records) + len(other_tert_records)
    result.notes.append("PDTC/ATC samples are MSK cohort (no TCGA barcode overlap) — useful for external validation only")
    write_attempt_log(SOURCE_ID, attempts)
    write_result(result)
    return result


if __name__ == "__main__":
    asyncio.run(run())
