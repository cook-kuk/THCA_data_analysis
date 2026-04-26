#!/usr/bin/env python3
"""v5 Orchestrator.

Phase 1 (parallel via ProcessPoolExecutor): Tracks 1, 2, 3 (max_workers=3).
Phase 2 (sequential): Track 4 (depends on T1+T2 outputs).
Phase 3 (sequential): build pages + reports + append banner.
Phase 4: health check curl.

Logs to logs/v5_run.log.
"""
from __future__ import annotations

import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

ROOT = Path("/opt/thyroid-dash/project")
SCRIPTS = ROOT / "notebooks_or_scripts"
LOG_PATH = ROOT / "logs" / "v5_run.log"
PY = ROOT / ".venv" / "bin" / "python"


def log(msg: str) -> None:
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    with LOG_PATH.open("a") as fh:
        fh.write(line + "\n")


def run_script(name: str, timeout: int = 2400):
    """Run a v5 track script, capture output, return dict."""
    t0 = time.time()
    script = SCRIPTS / name
    out_log = ROOT / "logs" / f"{script.stem}.log"
    try:
        with out_log.open("w") as fh:
            proc = subprocess.run(
                [str(PY), str(script)],
                stdout=fh, stderr=subprocess.STDOUT,
                timeout=timeout, cwd=str(ROOT))
        rc = proc.returncode
        status = "OK" if rc == 0 else f"FAIL(rc={rc})"
    except subprocess.TimeoutExpired:
        rc = -1
        status = "TIMEOUT"
    except Exception as e:
        rc = -2
        status = f"ERROR:{e}"
    elapsed = round(time.time() - t0, 1)
    return {"name": name, "status": status, "rc": rc,
            "elapsed_sec": elapsed, "log": str(out_log)}


def main() -> int:
    LOG_PATH.parent.mkdir(exist_ok=True)
    log("=== v5 orchestrator start ===")

    # Phase 1: parallel
    phase1 = [
        "v5_track1_nonlinear_correction.py",
        "v5_track2_adversarial.py",
        "v5_track3_cross_cancer.py",
    ]
    log(f"Phase 1: parallel run of {phase1} (max_workers=3, per-script timeout=2400s)")
    results = {}
    with ProcessPoolExecutor(max_workers=3) as ex:
        futs = {ex.submit(run_script, s, 2400): s for s in phase1}
        for fut in as_completed(futs):
            name = futs[fut]
            try:
                r = fut.result()
            except Exception as e:
                r = {"name": name, "status": f"CRASH:{e}", "rc": -3,
                     "elapsed_sec": -1, "log": ""}
            results[name] = r
            log(f"  phase1 done: {r['name']} -> {r['status']} ({r['elapsed_sec']}s)")

    # Phase 2: Track 4 (sequential, depends on T1+T2)
    log("Phase 2: sequential Track 4")
    r = run_script("v5_track4_robust_targets.py", timeout=900)
    results["v5_track4_robust_targets.py"] = r
    log(f"  phase2 done: {r['status']} ({r['elapsed_sec']}s)")

    # Phase 3: synthesis
    log("Phase 3: synthesis builder (pages + reports + banner)")
    r = run_script("v5_synth_builder.py", timeout=300)
    results["v5_synth_builder.py"] = r
    log(f"  phase3 done: {r['status']} ({r['elapsed_sec']}s)")

    # Phase 4: health check
    log("Phase 4: health check")
    base = "http://127.0.0.1:8012/reports/html/pages"
    for pg in ["34_v5_nonlinear_correction.html",
               "35_v5_adversarial_method.html",
               "36_v5_cross_cancer.html",
               "37_v5_robust_targets.html",
               "38_v5_synthesis.html"]:
        url = f"{base}/{pg}"
        try:
            proc = subprocess.run(
                ["curl", "-o", "/dev/null", "-s", "-w", "%{http_code}", url],
                capture_output=True, text=True, timeout=5)
            log(f"  curl {pg} -> {proc.stdout.strip()}")
        except Exception as e:
            log(f"  curl {pg} -> ERROR {e}")

    log("=== v5 orchestrator complete ===")

    # print a compact summary
    print()
    print("=== v5 ORCHESTRATION SUMMARY ===")
    for k, r in results.items():
        print(f"  {k}: {r['status']} ({r['elapsed_sec']}s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
