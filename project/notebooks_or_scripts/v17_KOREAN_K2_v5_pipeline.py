#!/usr/bin/env python3
"""v17 KOREAN K2 v5 — Pipelined download + quant + immediate cleanup.

Architecture:
  3 download workers (I/O bound) → bounded queue → 2 quant workers (CPU, kallisto -t 4) → cleanup
  Peak disk: ~5–7 R1 FASTQs * 1.3 GB ≈ 7 GB (vs 130 GB if all downloaded first).
  Target: 100 total quanted (existing 9 already done; this script processes the next 91).

Run: .venv/bin/python notebooks_or_scripts/v17_KOREAN_K2_v5_pipeline.py [N=100]
"""
from __future__ import annotations
import json
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from threading import Lock

import pandas as pd

ROOT = Path("/opt/thyroid-dash/project/results/v17_korean")
KAL = "/opt/thyroid-dash/project/.venv/lib/python3.12/site-packages/kb_python/bins/linux/kallisto/kallisto"
DATA_REF = Path("/data/thca/reference_kallisto")
DATA_FASTQ = Path("/data/thca/PRJEB11591_fastq")
QUANT = Path("/data/thca/PRJEB11591_quant_se")
WORK = Path("/data/thca/k2_work")
DATA_FASTQ.mkdir(parents=True, exist_ok=True)
QUANT.mkdir(parents=True, exist_ok=True)
(WORK / "tmp").mkdir(parents=True, exist_ok=True)

DL_WORKERS = 3
QUANT_WORKERS = 2
KALLISTO_THREADS = 4

print_lock = Lock()


def log(msg):
    with print_lock:
        print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def download_r1(run, ftp_str):
    """Download R1 FASTQ only. Returns path or None."""
    urls = [u.strip() for u in str(ftp_str).split(";") if u.strip()]
    urls = ["http://" + u if not u.startswith("http") else u for u in urls]
    r1_url = next((u for u in urls if "_1.fastq.gz" in u), None)
    if r1_url is None:
        log(f"  ✗ {run} no R1 url")
        return None
    out_dir = DATA_FASTQ / run
    out_dir.mkdir(parents=True, exist_ok=True)
    f1 = out_dir / r1_url.split("/")[-1]
    if f1.exists() and f1.stat().st_size > 500_000_000:
        log(f"  ↻ {run} R1 cached ({f1.stat().st_size//1_000_000} MB)")
        return str(f1)
    t0 = time.time()
    try:
        r = subprocess.run(
            ["wget", "-q", "-c", "--timeout=600", "--tries=2", "-O", str(f1), r1_url],
            capture_output=True, timeout=1800,
        )
        if r.returncode != 0 or f1.stat().st_size < 200_000_000:
            log(f"  ✗ {run} download failed (rc={r.returncode}, size={f1.stat().st_size})")
            try: f1.unlink()
            except Exception: pass
            return None
        sz = f1.stat().st_size // 1_000_000
        log(f"  ⬇  {run} R1 {sz} MB in {time.time()-t0:.0f}s")
        return str(f1)
    except subprocess.TimeoutExpired:
        log(f"  ✗ {run} download timeout")
        try: f1.unlink()
        except Exception: pass
        return None


def quant_and_cleanup(run, f1, idx):
    out = QUANT / run
    out.mkdir(parents=True, exist_ok=True)
    if (out / "abundance.tsv").exists():
        log(f"  ↻ {run} quant cached")
        try: Path(f1).unlink()
        except Exception: pass
        return run
    t0 = time.time()
    try:
        r = subprocess.run(
            [KAL, "quant", "-i", str(idx), "-o", str(out),
             "--single", "-l", "200", "-s", "30", "-t", str(KALLISTO_THREADS),
             f1],
            cwd=str(WORK), capture_output=True, text=True, timeout=900,
        )
        if r.returncode != 0:
            log(f"  ✗ {run} quant rc={r.returncode}: {r.stderr[:120]}")
            return None
        log(f"  ✓ {run} quanted in {time.time()-t0:.0f}s")
    except Exception as e:
        log(f"  ✗ {run} quant exception: {e}")
        return None
    finally:
        try: Path(f1).unlink()
        except Exception: pass
    return run


def main():
    target_n = int(sys.argv[1]) if len(sys.argv) > 1 else 100

    idx = DATA_REF / "gencode.v44.8gene.kallisto.idx"
    if not idx.exists():
        print(f"[ERROR] index missing: {idx}")
        return 1

    runs_df = pd.read_csv(ROOT / "K1A_prjeb11591_runs.tsv", sep="\t").sort_values("run_accession")
    already_done = {d.name for d in QUANT.iterdir() if (d / "abundance.tsv").exists()}
    log(f"  index OK · {len(runs_df)} runs in metadata · {len(already_done)} already quanted")

    pending = runs_df[~runs_df["run_accession"].isin(already_done)].head(target_n - len(already_done))
    log(f"  target N={target_n} → {len(pending)} new samples to process")
    if len(pending) == 0:
        log("  nothing to do")
        return 0

    t_start = time.time()
    completed = []
    failed = []

    with ThreadPoolExecutor(max_workers=DL_WORKERS, thread_name_prefix="dl") as dl_pool, \
         ThreadPoolExecutor(max_workers=QUANT_WORKERS, thread_name_prefix="qt") as qt_pool:

        dl_futures = {dl_pool.submit(download_r1, row["run_accession"], row["fastq_ftp"]): row["run_accession"]
                      for _, row in pending.iterrows()}
        qt_futures = {}

        for dl_fut in as_completed(dl_futures):
            run = dl_futures[dl_fut]
            f1 = dl_fut.result()
            if f1 is None:
                failed.append(run)
                continue
            qt_fut = qt_pool.submit(quant_and_cleanup, run, f1, idx)
            qt_futures[qt_fut] = run

        for qt_fut in as_completed(qt_futures):
            run = qt_futures[qt_fut]
            res = qt_fut.result()
            if res:
                completed.append(res)
            else:
                failed.append(run)
            elapsed = time.time() - t_start
            log(f"  [{len(completed)}/{len(pending)}] {run} pipeline done "
                f"(elapsed {elapsed/60:.1f} min, ~{elapsed/max(len(completed),1):.0f}s/sample)")

    log(f"\n  ✓ {len(completed)} new samples quanted in {(time.time()-t_start)/60:.1f} min "
        f"({len(failed)} failed)")
    log(f"  total quanted: {len(already_done) + len(completed)}")

    log("  → run notebooks_or_scripts/v17_KOREAN_K2_v4_fix.py to aggregate + LogReg")
    return 0


if __name__ == "__main__":
    sys.exit(main())
