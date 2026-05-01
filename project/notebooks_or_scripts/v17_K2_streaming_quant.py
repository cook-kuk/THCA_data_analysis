#!/usr/bin/env python3
"""K2 streaming quant — download R1 → kallisto SE quant → delete FASTQ → repeat.

Processes the 151+ runs not yet quanted in /data/thca/PRJEB11591_quant_se/.
After each batch of 10 runs, re-runs K2 v4 fix to update predictions.
"""
from __future__ import annotations
import os, sys, time, subprocess, logging
from pathlib import Path
import pandas as pd

ROOT = Path("/opt/thyroid-dash/project")
RES  = ROOT/"results/v17_korean"
LOG  = ROOT/"logs/v17_K2_streaming.log"
QUANT = Path("/data/thca/PRJEB11591_quant_se")
FASTQ = Path("/data/thca/PRJEB11591_fastq")
KAL = "/opt/thyroid-dash/project/.venv/lib/python3.12/site-packages/kb_python/bins/linux/kallisto/kallisto"
IDX = Path("/data/thca/reference_kallisto/gencode.v44.8gene.kallisto.idx")
WORK = Path("/data/thca/k2_work"); (WORK/"tmp").mkdir(parents=True, exist_ok=True)

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(message)s",
                    handlers=[logging.FileHandler(LOG, mode="a"), logging.StreamHandler(sys.stdout)])
log = logging.getLogger(__name__)

def get_todo_runs():
    manifest = pd.read_csv(RES/"K1A_prjeb11591_runs.tsv", sep="\t")
    all_runs = set(manifest["run_accession"].tolist())
    done = {d.name for d in QUANT.iterdir() if (d/"abundance.tsv").exists()}
    todo = sorted(all_runs - done)
    # Get FTP URL per todo run
    url_by_run = {}
    for _, row in manifest.iterrows():
        if row["run_accession"] in todo:
            urls = str(row["fastq_ftp"]).split(";")
            r1 = next((u for u in urls if "_1.fastq.gz" in u), None)
            if r1:
                url_by_run[row["run_accession"]] = r1
    return todo, url_by_run

def download_r1(run, url):
    """Download R1 FASTQ for one run; returns path or None on fail."""
    out_dir = FASTQ/run; out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir/Path(url).name
    if out_file.exists() and out_file.stat().st_size > 500_000_000:
        return out_file
    full_url = url if url.startswith(("http://","https://","ftp://")) else f"https://{url}"
    t0 = time.time()
    try:
        r = subprocess.run(["curl","-sSL","--retry","3","--retry-delay","5",
                            "--max-time","1800","-o",str(out_file), full_url],
                           capture_output=True, text=True, timeout=2000)
        if r.returncode != 0 or not out_file.exists() or out_file.stat().st_size < 100_000_000:
            log.warning(f"  ! [{run}] download fail rc={r.returncode}, size={out_file.stat().st_size if out_file.exists() else 0}")
            try: out_file.unlink()
            except: pass
            return None
        log.info(f"  ↓ [{run}] downloaded {out_file.stat().st_size/1e9:.2f} GB in {time.time()-t0:.0f}s")
        return out_file
    except Exception as e:
        log.warning(f"  ! [{run}] dl exception: {e}")
        return None

def quant_se(run, fq):
    out = QUANT/run; out.mkdir(parents=True, exist_ok=True)
    if (out/"abundance.tsv").exists():
        return out/"abundance.tsv"
    t0 = time.time()
    try:
        r = subprocess.run([KAL,"quant","-i",str(IDX),"-o",str(out),
                            "--single","-l","200","-s","30","-t","4", str(fq)],
                           cwd=str(WORK), capture_output=True, text=True, timeout=600)
        if r.returncode != 0:
            log.warning(f"  ! [{run}] quant rc={r.returncode}: {r.stderr[:200]}")
            return None
        log.info(f"  ✓ [{run}] quanted in {time.time()-t0:.0f}s")
        return out/"abundance.tsv"
    except Exception as e:
        log.warning(f"  ! [{run}] quant exception: {e}")
        return None

def cleanup_fastq(fq):
    try: fq.unlink()
    except: pass

def rerun_predictions():
    """Re-run K2 v4 fix to update predictions on full quant set."""
    log.info("=== Re-running K2 v4 fix on updated quant set ===")
    r = subprocess.run(["python3","-u","notebooks_or_scripts/v17_KOREAN_K2_v4_fix.py"],
                       cwd=str(ROOT), capture_output=True, text=True, timeout=300)
    if r.returncode == 0:
        # Print last few lines of output
        log.info("K2 v4 update OK; tail of stdout:")
        for line in r.stdout.split("\n")[-10:]:
            if line.strip(): log.info(f"  {line}")
    else:
        log.warning(f"K2 v4 update fail rc={r.returncode}: {r.stderr[:200]}")

def main():
    log.info("=== K2 streaming quant START ===")
    todo, url_by_run = get_todo_runs()
    log.info(f"TODO: {len(todo)} runs (manifest 262 minus existing quants)")

    BATCH = 10
    success = 0
    for i, run in enumerate(todo):
        if run not in url_by_run:
            log.warning(f"  [{i+1}/{len(todo)}] {run}: no FTP URL"); continue
        log.info(f"  [{i+1}/{len(todo)}] {run}: download...")
        fq = download_r1(run, url_by_run[run])
        if fq is None:
            log.warning(f"  [{i+1}/{len(todo)}] {run}: skipped"); continue
        log.info(f"  [{i+1}/{len(todo)}] {run}: quant...")
        ab = quant_se(run, fq)
        cleanup_fastq(fq)  # save disk regardless
        if ab is None: continue
        success += 1
        # Re-predict every BATCH runs
        if success % BATCH == 0:
            rerun_predictions()

    log.info(f"=== DONE: {success}/{len(todo)} new quants ===")
    rerun_predictions()  # final update

if __name__ == "__main__":
    sys.exit(main())
