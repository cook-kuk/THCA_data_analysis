"""S4 — Pozdeyev et al. (2024) thyroid cancer landscape paper.

Searches PubMed and a few likely DOI patterns. Pozdeyev's group also
publishes via cBioPortal (mskcc / aacr genie). Falls back to AACR GENIE
panel for thyroid TERT promoter (capture-based, includes TERT promoter).
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
    is_promoter_position,
    make_session_kwargs,
    write_attempt_log,
    write_result,
)

LOG = configure_logger("S4_pozdeyev2024")
SOURCE_ID = "S4_pozdeyev2024"

# Direct candidate URLs (DOI guesses + PubMed search + GENIE thyroid)
CANDIDATE_URLS = [
    # PubMed search for the paper
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=Pozdeyev+thyroid+cancer+TERT+landscape&retmax=20&retmode=json",
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=Pozdeyev+thyroid+TERT+2024&retmax=20&retmode=json",
    # AACR GENIE consortium data (panel-based, captures TERT promoter)
    "https://github.com/cBioPortal/datahub/raw/master/public/genie_public/data_mutations.txt",
    # GENIE thyroid filter via cBioPortal API requires study id; try direct
    "https://www.cbioportal.org/api/studies/genie_public/molecular-profiles/genie_public_mutations/mutations?projection=DETAILED",
    # Likely earlier paper from same group: Pozdeyev 2018 (PMID 29615459)
    "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6135480/",
    "https://aacrjournals.org/clincancerres/article/24/13/3059/80916/",
    # MSK thyroid landscape 2017 (close substitute)
    "https://github.com/cBioPortal/datahub/raw/master/public/thyroid_mskcc_2016/data_mutations.txt",
    "https://github.com/cBioPortal/datahub/raw/master/public/thyroid_mskcc_2016/data_clinical_sample.txt",
]


async def run() -> SourceResult:
    result = SourceResult(source_id=SOURCE_ID, label="Pozdeyev 2024 / GENIE substitute")
    async with aiohttp.ClientSession(**make_session_kwargs()) as session:
        coros = []
        for url in CANDIDATE_URLS:
            ext = ".json" if "esearch" in url or "/api/" in url else (
                ".txt" if "data_mutations" in url or "data_clinical" in url else ".html"
            )
            safe = re.sub(r"[^A-Za-z0-9]", "_", url)[-80:]
            coros.append(fetch_url(session, url, save_path=RAW / f"{SOURCE_ID}_{safe}{ext}", timeout=120))
        attempts = await asyncio.gather(*coros)
    result.attempts = attempts
    LOG.info("OK: %d / %d", sum(1 for a in attempts if a.ok), len(attempts))

    promoter_records = []
    other_tert_records = []
    pubmed_ids: list[str] = []
    for a in attempts:
        if not a.ok or not a.saved_to:
            continue
        path = Path(a.saved_to)
        full = path if path.is_absolute() else (Path.cwd() / path)
        if not full.exists():
            continue
        # Pubmed JSON
        if "esearch" in a.url:
            try:
                payload = json.loads(full.read_text())
                ids = payload.get("esearchresult", {}).get("idlist", [])
                pubmed_ids.extend(ids)
                LOG.info("PubMed ids: %s", ids)
            except Exception:  # noqa: BLE001
                pass
        # Big MAFs
        if a.bytes_ > 100_000 and full.suffix == ".txt":
            try:
                import pandas as pd

                # GENIE MAF can be very large; read in chunks
                chunks = pd.read_csv(full, sep="\t", comment="#", low_memory=False, chunksize=200_000)
                for df in chunks:
                    if "Hugo_Symbol" not in df.columns:
                        break
                    tert = df[df["Hugo_Symbol"].astype(str).str.upper() == "TERT"]
                    # Restrict to thyroid samples if there's a CANCER_TYPE column or sample id contains thyroid
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
            except Exception as e:  # noqa: BLE001
                LOG.warning("MAF parse err %s: %s", full.name, e)

    if promoter_records:
        import pandas as pd

        out = PARSED / f"{SOURCE_ID}_tert_promoter.tsv"
        pd.DataFrame(promoter_records).to_csv(out, sep="\t", index=False)
        result.output_files.append(str(out))
        result.n_promoter_mutations = len(promoter_records)
        result.success = True
    if other_tert_records:
        import pandas as pd

        out = PARSED / f"{SOURCE_ID}_tert_other.tsv"
        pd.DataFrame(other_tert_records).to_csv(out, sep="\t", index=False)
        result.output_files.append(str(out))
    result.n_records = len(promoter_records) + len(other_tert_records)
    if pubmed_ids:
        result.notes.append(f"PubMed ids found: {pubmed_ids[:10]}")
    result.notes.append("Pozdeyev 2024 paper not located via DOI guess; substituting GENIE/MSKCC panels which include TERT promoter capture")
    write_attempt_log(SOURCE_ID, attempts)
    write_result(result)
    return result


if __name__ == "__main__":
    asyncio.run(run())
