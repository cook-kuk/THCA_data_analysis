#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import requests

from pilot_utils import append_missing_data, config, ensure_standard_dirs, maybe_symlink, pilot_root, setup_logging


def download(url: str, dst: Path, timeout: int = 60) -> bool:
    try:
        with requests.get(url, stream=True, timeout=timeout) as response:
            response.raise_for_status()
            dst.parent.mkdir(parents=True, exist_ok=True)
            with dst.open("wb") as handle:
                for chunk in response.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        handle.write(chunk)
        return True
    except Exception:
        if dst.exists():
            dst.unlink()
        return False


def main() -> None:
    parser = argparse.ArgumentParser(description="Link or download TCGA-THCA expression and clinical data.")
    parser.parse_args()
    ensure_standard_dirs()
    logger = setup_logging("02_download_or_link_tcga_thca")
    cfg = config()
    raw = pilot_root() / "data_raw" / "tcga_thca"
    rows = []

    local = cfg.get("local_sources", {})
    for key in [
        "tcga_pancan_expression",
        "tcga_pancan_phenotype",
        "tcga_pancan_survival",
        "tcga_thca_clinical_extended",
        "tcga_thca_mutation_groups",
    ]:
        src = Path(local.get(key, ""))
        if src.exists():
            dst = raw / src.name
            maybe_symlink(src, dst)
            rows.append({"source_key": key, "status": "linked_local", "path": str(dst), "original_path": str(src)})
            logger.info("Linked %s from %s", key, src)
        else:
            rows.append({"source_key": key, "status": "missing_local", "path": "", "original_path": str(src)})

    if not any(r["source_key"] == "tcga_pancan_expression" and r["status"] == "linked_local" for r in rows):
        remote = cfg.get("remote_sources", {})
        expr_dst = raw / "TCGA.THCA.HiSeqV2.gz"
        clin_dst = raw / "TCGA.THCA.clinicalMatrix.gz"
        expr_ok = download(remote.get("tcga_thca_xena_expression", ""), expr_dst)
        clin_ok = download(remote.get("tcga_thca_xena_clinical", ""), clin_dst)
        rows.append({"source_key": "tcga_thca_xena_expression", "status": "downloaded" if expr_ok else "download_failed", "path": str(expr_dst) if expr_ok else "", "original_path": remote.get("tcga_thca_xena_expression", "")})
        rows.append({"source_key": "tcga_thca_xena_clinical", "status": "downloaded" if clin_ok else "download_failed", "path": str(clin_dst) if clin_ok else "", "original_path": remote.get("tcga_thca_xena_clinical", "")})
        if not expr_ok:
            append_missing_data(
                "TCGA-THCA expression",
                "Automatic local linking and UCSC Xena download failed. Manually download TCGA-THCA RNA-seq from UCSC Xena: https://xenabrowser.net/datapages/?cohort=TCGA%20Thyroid%20Cancer%20(THCA)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443",
            )

    manifest = pd.DataFrame(rows)
    out = raw / "tcga_thca_sources.tsv"
    manifest.to_csv(out, sep="\t", index=False)
    logger.info("Wrote TCGA source manifest to %s", out)


if __name__ == "__main__":
    main()
