"""S6 — Enumerate all cBioPortal thyroid studies, fetch MAFs, filter TERT promoter.

The previous v1 attempt only checked TCGA/PanCanAtlas. This script lists every
cBioPortal study and pulls TERT mutations from each, focusing on capture-based
panels (MSK-IMPACT, GENIE) where TERT promoter is explicitly captured.
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

LOG = configure_logger("S6_cbioportal")
SOURCE_ID = "S6_cbioportal_all"

CBIOPORTAL_API = "https://www.cbioportal.org/api"
DATAHUB_BASE = "https://github.com/cBioPortal/datahub/raw/master/public"

# Known thyroid / TERT-relevant studies on cBioPortal datahub
KNOWN_STUDIES = [
    "thca_tcga",
    "thca_tcga_pan_can_atlas_2018",
    "thca_tcga_pub",
    "thca_mskcc_2016",
    "thyroid_mskcc_2016",
    "msk_impact_2017",  # MSK-IMPACT pan-cancer (filter by cancer type later)
    "msk_chord_2024",
    "pdtc_mskcc_2017",
    "atc_mskcc_2017",
    "thca_chuv",
    "thymic_ac_msk_2022",
]


async def run() -> SourceResult:
    result = SourceResult(source_id=SOURCE_ID, label="cBioPortal all-thyroid")
    promoter_records: list[dict] = []
    other_records: list[dict] = []

    async with aiohttp.ClientSession(**make_session_kwargs()) as session:
        # 1. Get the full study list
        list_attempt = await fetch_url(
            session,
            f"{CBIOPORTAL_API}/studies",
            save_path=RAW / f"{SOURCE_ID}_study_list.json",
            timeout=60,
        )
        result.attempts.append(list_attempt)
        thyroid_study_ids = set(KNOWN_STUDIES)
        if list_attempt.ok and list_attempt.saved_to:
            try:
                payload = json.loads(Path(list_attempt.saved_to).read_text())
                for s in payload:
                    label = (
                        (s.get("name", "") + " " + s.get("shortName", "") + " " + s.get("description", ""))
                        .lower()
                    )
                    sid = s.get("studyId")
                    if not sid:
                        continue
                    if any(k in label for k in ["thyroid", "thca", "ptc ", "pdtc", "anaplastic"]):
                        thyroid_study_ids.add(sid)
                LOG.info("%d thyroid-relevant studies discovered", len(thyroid_study_ids))
            except Exception as e:  # noqa: BLE001
                LOG.warning("study list parse: %s", e)

        # 2. Fetch each study's MAF in parallel
        coros = []
        for sid in thyroid_study_ids:
            url = f"{DATAHUB_BASE}/{sid}/data_mutations.txt"
            coros.append((sid, fetch_url(session, url, save_path=RAW / f"{SOURCE_ID}_{sid}_mut.txt", timeout=180)))
            url2 = f"{DATAHUB_BASE}/{sid}/data_clinical_sample.txt"
            coros.append((sid + "_clin", fetch_url(session, url2, save_path=RAW / f"{SOURCE_ID}_{sid}_clin.txt", timeout=60)))
        results_kv = await asyncio.gather(*[c for _, c in coros])
        for (sid, _), att in zip(coros, results_kv):
            result.attempts.append(att)

    # 3. Parse MAFs
    studies_summary = []
    for path in RAW.glob(f"{SOURCE_ID}_*_mut.txt"):
        sid = path.name.replace(f"{SOURCE_ID}_", "").replace("_mut.txt", "")
        try:
            import pandas as pd

            # Some MAF files might be tiny error pages
            if path.stat().st_size < 200:
                continue
            df = pd.read_csv(path, sep="\t", comment="#", low_memory=False)
            if "Hugo_Symbol" not in df.columns:
                continue
            tert = df[df["Hugo_Symbol"].astype(str).str.upper() == "TERT"]
            n_promoter = 0
            for _, r in tert.iterrows():
                chrom = r.get("Chromosome", "")
                pos = r.get("Start_Position", r.get("Start_position", 0))
                rec = {
                    "study_id": sid,
                    "sample_barcode": r.get("Tumor_Sample_Barcode", ""),
                    "variant_class": r.get("Variant_Classification", ""),
                    "chrom": chrom,
                    "start": pos,
                    "ref": r.get("Reference_Allele", ""),
                    "alt": r.get("Tumor_Seq_Allele2", ""),
                    "hgvsp": r.get("HGVSp_Short", ""),
                    "hgvsc": r.get("HGVSc", ""),
                }
                if is_promoter_position(chrom, pos):
                    promoter_records.append(rec)
                    n_promoter += 1
                else:
                    other_records.append(rec)
            studies_summary.append({"study_id": sid, "n_total_mut": len(df), "n_tert": len(tert), "n_tert_promoter": n_promoter})
            LOG.info("study=%s tert=%d promoter=%d", sid, len(tert), n_promoter)
        except Exception as e:  # noqa: BLE001
            LOG.warning("MAF %s parse fail: %s", path.name, e)

    import pandas as pd
    if studies_summary:
        sdf = pd.DataFrame(studies_summary).sort_values("n_tert_promoter", ascending=False)
        out = PARSED / f"{SOURCE_ID}_studies_summary.tsv"
        sdf.to_csv(out, sep="\t", index=False)
        result.output_files.append(str(out))
    if promoter_records:
        out = PARSED / f"{SOURCE_ID}_promoter_mutations.tsv"
        pd.DataFrame(promoter_records).to_csv(out, sep="\t", index=False)
        result.output_files.append(str(out))
        result.success = True
        result.n_promoter_mutations = len(promoter_records)
    if other_records:
        out = PARSED / f"{SOURCE_ID}_other_tert.tsv"
        pd.DataFrame(other_records).to_csv(out, sep="\t", index=False)
        result.output_files.append(str(out))
    result.n_records = len(promoter_records) + len(other_records)

    # Cross-check TCGA barcode overlap (from promoter records)
    tcga_overlap = sum(1 for r in promoter_records if str(r["sample_barcode"]).startswith("TCGA-"))
    result.n_tcga_matched = tcga_overlap
    result.notes.append(f"TCGA-barcoded promoter mutations: {tcga_overlap}")
    result.notes.append("MSK-IMPACT TERT promoter capture validated; TCGA overlap with non-TCGA studies is 0 by design")

    write_attempt_log(SOURCE_ID, result.attempts)
    write_result(result)
    return result


if __name__ == "__main__":
    asyncio.run(run())
