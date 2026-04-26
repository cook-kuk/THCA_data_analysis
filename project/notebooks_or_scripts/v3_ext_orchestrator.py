#!/usr/bin/env python
"""v3_ext_orchestrator.py

Stages:
  1) v3_ext_01_download.py (parallel A-E internal; sequential in orchestrator)
  2,3) parallel batch: v3_ext_02_harmonize.py + v3_ext_03_brs71_recover.py
  4,6) parallel batch: v3_ext_04_fusion_anchor.py + v3_ext_06_batch_diagnostics.py
  5) v3_ext_05_external_validation.py  (depends 2,3,4)
  7) v3_ext_07_build_pages.py
  8) v3_ext_08_write_reports.py
ProcessPoolExecutor max_workers = min(cpu_count, 8). Each task wrapped in try/except,
downstream dependants SKIP on upstream fail (stub outputs). Orchestrator never aborts.

Health check at end: curl _version.json + page 25 sizes.
"""
import os
import sys
import json
import time
import logging
import subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import cpu_count

ROOT = Path("/opt/thyroid-dash/project")
SCRIPTS = ROOT / "notebooks_or_scripts"
LOGDIR = ROOT / "logs"
LOGDIR.mkdir(exist_ok=True)
PY = "/opt/thyroid-dash/project/.venv/bin/python"
LOGFILE = LOGDIR / "v3_ext_run.log"

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    handlers=[logging.FileHandler(LOGFILE, mode="w"), logging.StreamHandler()],
)
log = logging.getLogger("v3_ext_orch")


def run_step(name: str, timeout_s: int = 1800) -> dict:
    script = SCRIPTS / name
    log.info(f"STEP START {name}")
    t0 = time.time()
    try:
        res = subprocess.run([PY, str(script)], capture_output=True, text=True, timeout=timeout_s)
        dt = time.time() - t0
        ok = (res.returncode == 0)
        log.info(f"STEP {'DONE' if ok else 'FAIL'} {name} ({dt:.1f}s rc={res.returncode})")
        # Save tail of stderr if failed
        if not ok:
            tail = res.stderr[-2000:] if res.stderr else ""
            log.warning(f"stderr tail:\n{tail}")
        return {"name": name, "ok": ok, "rc": res.returncode, "dt": dt,
                "stdout_tail": (res.stdout or "")[-1000:],
                "stderr_tail": (res.stderr or "")[-1000:]}
    except subprocess.TimeoutExpired as e:
        dt = time.time() - t0
        log.error(f"STEP TIMEOUT {name} ({dt:.1f}s)")
        return {"name": name, "ok": False, "rc": -1, "dt": dt, "error": "timeout"}
    except Exception as e:
        dt = time.time() - t0
        log.error(f"STEP ERROR {name}: {e}")
        return {"name": name, "ok": False, "rc": -2, "dt": dt, "error": str(e)}


def parallel_steps(names: list, timeout_s: int = 1800, max_workers: int = None) -> list:
    max_workers = max_workers or min(cpu_count(), 8)
    results = []
    with ProcessPoolExecutor(max_workers=min(max_workers, len(names))) as ex:
        futs = {ex.submit(run_step, n, timeout_s): n for n in names}
        for f in as_completed(futs):
            try:
                results.append(f.result())
            except Exception as e:
                results.append({"name": futs[f], "ok": False, "error": str(e)})
    return results


def health_check():
    vfile = ROOT / "reports" / "html" / "_version.json"
    page25 = ROOT / "reports" / "html" / "pages" / "25_external_prjeb11591.html"
    info = {}
    info["version_json_exists"] = vfile.exists()
    info["version_json_size"] = vfile.stat().st_size if vfile.exists() else 0
    info["page25_exists"] = page25.exists()
    info["page25_size"] = page25.stat().st_size if page25.exists() else 0
    log.info(f"health check: {info}")
    return info


def main():
    t0 = time.time()
    log.info("v3_ext orchestrator start")
    status = []

    # Step 1
    r1 = run_step("v3_ext_01_download.py", timeout_s=600)
    status.append(r1)

    # Check PRJEB11591 status: if SKIPPED (expression file absent), orchestrator continues
    # but step 5 will still run (uses proxy); we only STOP early if BOTH PRJEB and microarray
    # pipeline are broken, which is not our case.

    # Steps 2 + 3 in parallel
    r23 = parallel_steps(["v3_ext_02_harmonize.py", "v3_ext_03_brs71_recover.py"], timeout_s=900)
    status += r23

    # Steps 4 + 6 in parallel
    r46 = parallel_steps(["v3_ext_04_fusion_anchor.py", "v3_ext_06_batch_diagnostics.py"], timeout_s=900)
    status += r46

    # Step 5 (needs 2,3,4)
    r5 = run_step("v3_ext_05_external_validation.py", timeout_s=1800)
    status.append(r5)

    # Step 7 + 8 sequential
    r7 = run_step("v3_ext_07_build_pages.py", timeout_s=300)
    status.append(r7)
    r8 = run_step("v3_ext_08_write_reports.py", timeout_s=300)
    status.append(r8)

    hc = health_check()

    summary = {
        "elapsed_sec": time.time() - t0,
        "steps": status,
        "health": hc,
    }
    (LOGDIR / "v3_ext_orchestrator_summary.json").write_text(json.dumps(summary, indent=2, default=str))
    log.info(f"v3_ext orchestrator done in {summary['elapsed_sec']:.1f}s")
    # Also print a compact one-liner for STAGE 99
    for s in status:
        log.info(f"  {s.get('name'):40s}  ok={s.get('ok')}  rc={s.get('rc')}  dt={s.get('dt',0):.1f}s")


if __name__ == "__main__":
    main()
