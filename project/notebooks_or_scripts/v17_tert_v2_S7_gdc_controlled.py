"""S7 — GDC controlled-access scout for TCGA-THCA WGS / Targeted Sequencing.

The standard TCGA-THCA WES open-access MAFs do not capture TERT promoter.
This script lists which TCGA-THCA samples have WGS or targeted sequencing
data in GDC (potentially TERT-promoter-containing) and writes a dbGaP
application guide.
"""
from __future__ import annotations

import asyncio
import json
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

LOG = configure_logger("S7_gdc")
SOURCE_ID = "S7_gdc_controlled"

GDC_FILES_API = "https://api.gdc.cancer.gov/files"


def gdc_filter_payload(experimental_strategies: list[str]) -> dict:
    return {
        "filters": {
            "op": "and",
            "content": [
                {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-THCA"]}},
                {"op": "in", "content": {"field": "experimental_strategy", "value": experimental_strategies}},
            ],
        },
        "fields": "file_id,file_name,experimental_strategy,data_type,data_format,data_category,access,cases.submitter_id,cases.case_id,cases.samples.submitter_id,cases.samples.sample_type",
        "size": 5000,
        "format": "json",
    }


async def query(session: aiohttp.ClientSession, strategies: list[str], save_path: Path) -> dict:
    payload = gdc_filter_payload(strategies)
    try:
        async with session.post(GDC_FILES_API, json=payload, timeout=aiohttp.ClientTimeout(total=120)) as resp:
            data = await resp.read()
            if resp.status == 200:
                save_path.write_bytes(data)
                return {"ok": True, "saved": str(save_path), "bytes": len(data), "status": resp.status}
            return {"ok": False, "status": resp.status, "bytes": len(data)}
    except Exception as e:  # noqa: BLE001
        return {"ok": False, "error": f"{type(e).__name__}: {e}"}


async def run() -> SourceResult:
    result = SourceResult(source_id=SOURCE_ID, label="GDC controlled access scout")
    queries = {
        "wgs": ["WGS"],
        "targeted": ["Targeted Sequencing"],
        "wxs": ["WXS"],  # for completeness — already known limitation
    }
    summary = {}
    async with aiohttp.ClientSession(**make_session_kwargs()) as session:
        for k, strats in queries.items():
            save = RAW / f"{SOURCE_ID}_{k}.json"
            res = await query(session, strats, save)
            summary[k] = res
            LOG.info("query=%s -> %s", k, res)

    # Parse the JSONs
    rows = []
    for k in queries:
        p = RAW / f"{SOURCE_ID}_{k}.json"
        if not p.exists():
            continue
        try:
            payload = json.loads(p.read_text())
            hits = payload.get("data", {}).get("hits", [])
            LOG.info("%s hits: %d", k, len(hits))
            for h in hits:
                cases = h.get("cases", [{}])
                for c in cases:
                    rows.append(
                        {
                            "experimental_strategy": h.get("experimental_strategy"),
                            "data_type": h.get("data_type"),
                            "data_format": h.get("data_format"),
                            "access": h.get("access"),
                            "file_name": h.get("file_name"),
                            "file_id": h.get("file_id"),
                            "case_submitter_id": c.get("submitter_id"),
                            "case_id": c.get("case_id"),
                            "sample_types": ",".join(
                                s.get("sample_type", "") for s in c.get("samples", [])
                            ),
                            "sample_submitters": ",".join(
                                s.get("submitter_id", "") for s in c.get("samples", [])
                            ),
                        }
                    )
        except Exception as e:  # noqa: BLE001
            LOG.warning("parse %s: %s", p.name, e)

    if rows:
        import pandas as pd

        df = pd.DataFrame(rows)
        out = PARSED / f"{SOURCE_ID}_thca_seq_files.tsv"
        df.to_csv(out, sep="\t", index=False)
        result.output_files.append(str(out))
        wgs = df[df["experimental_strategy"] == "WGS"]
        targeted = df[df["experimental_strategy"] == "Targeted Sequencing"]
        result.notes.append(f"WGS files: {len(wgs)}, unique cases: {wgs['case_submitter_id'].nunique()}")
        result.notes.append(f"Targeted Seq files: {len(targeted)}, unique cases: {targeted['case_submitter_id'].nunique()}")
        controlled = df[df["access"] == "controlled"]
        result.notes.append(f"Controlled-access files: {len(controlled)}")
        result.n_records = len(df)
        result.success = True

    # Always write the application guide
    guide = V2_DIR / f"{SOURCE_ID}_application_guide.md"
    guide.write_text(
        """# dbGaP / GDC controlled-access application guide for TCGA-THCA TERT recovery

The TCGA-THCA WES MAFs in open access **do not capture TERT promoter** because
the promoter falls outside the WES exonic baits. To recover TERT promoter
status for TCGA-THCA, you need WGS BAMs (controlled access).

## Steps

1. **dbGaP application** — study accession `phs000178` (TCGA), data use limitation
   "General Research Use" or "Health/Medical/Biomedical".
2. Get an eRA Commons account if you don't have one (PI or sponsoring PI).
3. Submit a Data Access Request (DAR) on dbGaP web.
4. After approval (typically 2–6 weeks):
   - Download GDC token from the GDC portal.
   - Use `gdc-client` to download THCA WGS BAMs:
     ```bash
     gdc-client download -t gdc-user-token.txt -m thca_wgs_manifest.tsv
     ```
   - Generate a manifest from `parsed/S7_gdc_controlled_thca_seq_files.tsv` filtering to
     `experimental_strategy == "WGS"` and `data_format == "BAM"`.
5. **Variant calling at TERT promoter region**:
   - Region: chr5:1,294,500-1,295,800 (GRCh38) or chr5:1,295,000-1,296,000 (GRCh37).
   - Use `mutect2` or `varscan2` with low-VAF settings (TERT promoter is often
     subclonal, ~5-15% VAF).
   - Recommended: `gatk Mutect2 -R Homo_sapiens_assembly38.fasta -I tumor.bam -L TERT_promoter.bed`

## Realistic timeline
- DAR submission to data download: 4–8 weeks.
- Variant calling: 1–2 days for ~50 BAMs on a modest cluster.

This source is **flagged for future work**, not actionable in the current sprint.
""",
        encoding="utf-8",
    )
    result.output_files.append(str(guide))

    write_attempt_log(SOURCE_ID, result.attempts)
    write_result(result)
    return result


if __name__ == "__main__":
    asyncio.run(run())
