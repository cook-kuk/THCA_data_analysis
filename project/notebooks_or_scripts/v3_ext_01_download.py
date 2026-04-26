#!/usr/bin/env python
"""v3_ext_01_download.py

Acquire external cohorts for v3 external validation sprint:
  A) PRJEB11591 (Yoo 2016 PLOS Genet) -- processed counts via supplementary.
     If unobtainable, mark SKIPPED (spec-approved) and continue.
  B) TCGA fusions -- cBioPortal structural-variants API; TumorFusions fallback.
  C) TCGA SCNA -- Broad firehose GISTIC2 all_data_by_genes.
  D) GSE33630 (GPL570 PTC)
  E) GSE29265 (GPL570 PTC)

5 tasks in parallel ProcessPoolExecutor. Log to logs/v3_ext_01_download.log.
"""
import os
import sys
import json
import logging
import time
import io
import gzip
import tarfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import urllib.request
import urllib.error

ROOT = Path("/opt/thyroid-dash/project")
RAW = ROOT / "data_raw" / "v3_ext"
LOGDIR = ROOT / "logs"
RAW.mkdir(parents=True, exist_ok=True)
LOGDIR.mkdir(parents=True, exist_ok=True)

LOGFILE = LOGDIR / "v3_ext_01_download.log"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s %(message)s",
    handlers=[logging.FileHandler(LOGFILE, mode="w"), logging.StreamHandler()],
)
log = logging.getLogger("v3_ext_01")

STATUS_FILE = RAW / "download_status.json"


def _download(url: str, dest: Path, timeout: int = 120) -> bool:
    try:
        log.info(f"GET {url}")
        req = urllib.request.Request(url, headers={"User-Agent": "thca-v3-ext/1.0"})
        with urllib.request.urlopen(req, timeout=timeout) as r:
            dest.write_bytes(r.read())
        log.info(f"  wrote {dest} ({dest.stat().st_size} bytes)")
        return True
    except Exception as e:
        log.warning(f"  failed: {e}")
        return False


# ---------------- A) PRJEB11591 ---------------- #
def task_prjeb11591() -> dict:
    """Attempt to fetch Yoo 2016 PLOS Genet PRJEB11591 processed counts.

    Strategy:
      1) PLOS Genet article pgen.1006448 supplementary data (counts table).
      2) ENA filereport API -> (if processed counts unavailable, skip).
    """
    out = RAW / "prjeb11591"
    out.mkdir(exist_ok=True)
    candidates = [
        # Yoo et al 2016 PLOS Genet supplementary
        "https://journals.plos.org/plosgenetics/article/file?type=supplementary&id=info:doi/10.1371/journal.pgen.1006448.s008",
        "https://journals.plos.org/plosgenetics/article/file?type=supplementary&id=info:doi/10.1371/journal.pgen.1006448.s009",
    ]
    got = False
    for url in candidates:
        name = url.split("/")[-1].replace("info:doi/", "").replace("/", "_")
        dest = out / f"{name}.bin"
        if _download(url, dest, timeout=90):
            got = True
    ena = "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB11591&result=read_run&fields=run_accession,sample_accession,fastq_ftp,sample_title&format=tsv"
    _download(ena, out / "ena_filereport.tsv", timeout=60)
    # Processed count matrix commonly unobtainable from PLOS; mark SKIPPED if nothing usable
    status = "OK" if got else "SKIPPED"
    # Even if supplementary downloaded, we can't reliably parse RNA-seq counts without inspection.
    # The spec allows SKIPPED so we tag anchor downstream as histology_proxy using ENA metadata only.
    return {"task": "prjeb11591", "status": status, "dest": str(out)}


# ---------------- B) TCGA fusions ---------------- #
def task_tcga_fusions() -> dict:
    out = RAW / "tcga_fusions"
    out.mkdir(exist_ok=True)
    ok = False
    # cBioPortal public api - structural variants for thca_tcga_pan_can_atlas_2018
    base = "https://www.cbioportal.org/api"
    candidates = [
        f"{base}/structural-variant/fetch?structuralVariantMultiGeneQueries=%5B%5D&projection=SUMMARY",
    ]
    # Simpler: fetch all SV for thca
    sv_url = "https://www.cbioportal.org/api/molecular-profiles/thca_tcga_pan_can_atlas_2018_structural_variants/structural-variants"
    if _download(sv_url, out / "cbio_structural_variants.json", timeout=90):
        ok = True
    # TumorFusions fallback: static bulk table
    tf = "https://tumorfusions.org/downloads/pancanfus.txt.gz"
    _download(tf, out / "pancanfus.txt.gz", timeout=90)
    return {"task": "tcga_fusions", "status": "OK" if ok else "SKIPPED", "dest": str(out)}


# ---------------- C) TCGA SCNA (GISTIC2) ---------------- #
def task_tcga_scna() -> dict:
    out = RAW / "tcga_scna"
    out.mkdir(exist_ok=True)
    # Broad firehose GISTIC2 -- large tar.gz. Best-effort.
    urls = [
        "https://gdac.broadinstitute.org/runs/analyses__2016_01_28/data/THCA-TP/20160128/gdac.broadinstitute.org_THCA-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0.tar.gz",
    ]
    ok = False
    for u in urls:
        dest = out / "gistic2.tar.gz"
        if _download(u, dest, timeout=180):
            ok = True
            break
    return {"task": "tcga_scna", "status": "OK" if ok else "SKIPPED", "dest": str(out)}


# ---------------- D/E) GEO download (GSE33630 / GSE29265) ---------------- #
def _geo(gse: str) -> dict:
    out = RAW / gse
    out.mkdir(exist_ok=True)
    try:
        import GEOparse  # noqa
        log.info(f"GEOparse {gse} ...")
        gp = GEOparse.get_GEO(geo=gse, destdir=str(out), silent=True)
        # Write a compact expression matrix
        gsms = gp.gsms
        rows = []
        cols = []
        for gsm_name, gsm in gsms.items():
            df = gsm.table  # usually ID_REF, VALUE
            if df is None or df.empty:
                continue
            if "ID_REF" in df.columns and "VALUE" in df.columns:
                s = df.set_index("ID_REF")["VALUE"]
                rows.append(s)
                cols.append(gsm_name)
        if rows:
            import pandas as pd
            import numpy as np
            M = pd.concat(rows, axis=1)
            M.columns = cols
            # Force numeric, drop all-nan rows
            M = M.apply(pd.to_numeric, errors="coerce").dropna(how="all")
            # GPL570 is usually log2 already or linear. Heuristic: if median > 100, log2.
            med = float(np.nanmedian(M.values))
            if med > 50:
                M = np.log2(M.clip(lower=0) + 1.0)
            M.index.name = "probe_id"
            M.to_csv(out / f"{gse}_probe_matrix_log2.tsv", sep="\t")
            # Metadata
            meta_rows = []
            for gsm_name, gsm in gsms.items():
                row = {"sample_id": gsm_name}
                for k, v in gsm.metadata.items():
                    if isinstance(v, list):
                        v = "; ".join(map(str, v))
                    row[k] = v
                meta_rows.append(row)
            import pandas as pd
            pd.DataFrame(meta_rows).to_csv(out / f"{gse}_metadata.tsv", sep="\t", index=False)
            return {"task": gse, "status": "OK", "n_samples": len(cols), "dest": str(out)}
        return {"task": gse, "status": "SKIPPED", "dest": str(out), "reason": "no_tables"}
    except Exception as e:
        log.warning(f"{gse} failed: {e}")
        return {"task": gse, "status": "SKIPPED", "dest": str(out), "reason": str(e)[:200]}


def task_gse33630():
    return _geo("GSE33630")


def task_gse29265():
    return _geo("GSE29265")


def main():
    log.info("v3_ext_01_download start")
    tasks = [task_prjeb11591, task_tcga_fusions, task_tcga_scna, task_gse33630, task_gse29265]
    results = []
    # ProcessPool would re-import; simpler to run sequentially since each task is I/O-bound
    # but we parallelise with threads for I/O.
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=5) as ex:
        futs = {ex.submit(t): t.__name__ for t in tasks}
        for f in as_completed(futs):
            name = futs[f]
            try:
                r = f.result()
                results.append(r)
                log.info(f"[done] {name}: {r.get('status')}")
            except Exception as e:
                results.append({"task": name, "status": "ERROR", "error": str(e)})
                log.error(f"[err] {name}: {e}")
    STATUS_FILE.write_text(json.dumps(results, indent=2))
    log.info(f"status -> {STATUS_FILE}")
    log.info("v3_ext_01_download done")


if __name__ == "__main__":
    main()
