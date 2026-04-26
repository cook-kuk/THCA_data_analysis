#!/usr/bin/env python3
"""
v4 Quadruple-Track Sprint orchestrator.

Runs Tracks A-D in parallel (ProcessPoolExecutor max_workers=4), then runs
the synth builder sequentially, then a health check.

- Any single-track failure is logged but does NOT abort the remaining pipeline.
- Logs to logs/v4_run.log.
"""
from __future__ import annotations

import logging
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

ROOT = Path("/opt/thyroid-dash/project")
SCRIPTS = ROOT / "notebooks_or_scripts"
LOGS = ROOT / "logs"
LOGS.mkdir(parents=True, exist_ok=True)

PY = "/opt/thyroid-dash/project/.venv/bin/python"

TRACKS = {
    "A": SCRIPTS / "v4_trackA_combat_rescue.py",
    "B": SCRIPTS / "v4_trackB_bethesda_strengthen.py",
    "C": SCRIPTS / "v4_trackC_prjeb11591_retry.py",
    "D": SCRIPTS / "v4_trackD_honest_narrative.py",
}
SYNTH = SCRIPTS / "v4_synth_builder.py"

log_path = LOGS / "v4_run.log"
logging.basicConfig(
    filename=str(log_path),
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
logger = logging.getLogger("v4")
console = logging.StreamHandler()
console.setLevel(logging.INFO)
console.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
logger.addHandler(console)


def run_track(name: str, script: Path) -> tuple[str, int, str, float]:
    t0 = time.time()
    log_file = LOGS / f"v4_track{name}.log"
    try:
        proc = subprocess.run(
            [PY, str(script)],
            capture_output=True,
            text=True,
            timeout=3600,
        )
        out = proc.stdout + "\n---STDERR---\n" + proc.stderr
        log_file.write_text(out)
        rc = int(proc.returncode)
    except subprocess.TimeoutExpired:
        rc = -1
        out = "TIMEOUT"
        log_file.write_text(out)
    except Exception as e:
        rc = -2
        out = f"exception: {e}"
        log_file.write_text(out)
    return name, rc, str(log_file), time.time() - t0


def main() -> int:
    logger.info("v4 orchestrator start")
    logger.info(f"logs at {log_path}")

    # Track D depends on A+C to exist; allow it to still run even if empty (uses defaults)
    # and re-run after A+C if they produced data. Easiest: run all 4 in parallel, then run
    # Track D again after synth (we instead run D as part of synth flow by re-running D
    # synchronously after A/C finish, then synth).
    # Spec says run A-D in parallel, then synth. We will also re-run D after A/C completion
    # to pick up updated data, because D only reads outputs.
    results = {}

    with ProcessPoolExecutor(max_workers=4) as pool:
        futs = {pool.submit(run_track, name, path): name for name, path in TRACKS.items()}
        for f in as_completed(futs):
            name, rc, log_file, dt = f.result()
            results[name] = {"rc": rc, "log": log_file, "seconds": dt}
            logger.info(f"track {name} rc={rc} time={dt:.1f}s log={log_file}")

    # If Track D ran before A/B finished, re-run it once synchronously
    try:
        logger.info("re-running Track D to pick up latest A/B/C outputs...")
        name, rc, log_file, dt = run_track("D", TRACKS["D"])
        results["D_rerun"] = {"rc": rc, "log": log_file, "seconds": dt}
        logger.info(f"track D_rerun rc={rc} time={dt:.1f}s")
    except Exception as e:
        logger.warning(f"D rerun failed: {e}")

    # Synth
    try:
        logger.info("running synth builder...")
        log_file = LOGS / "v4_synth.log"
        t0 = time.time()
        proc = subprocess.run([PY, str(SYNTH)], capture_output=True, text=True, timeout=900)
        log_file.write_text(proc.stdout + "\n---STDERR---\n" + proc.stderr)
        rc = int(proc.returncode)
        results["SYNTH"] = {"rc": rc, "log": str(log_file), "seconds": time.time() - t0}
        logger.info(f"synth rc={rc}")
    except Exception as e:
        logger.error(f"synth failed: {e}")
        results["SYNTH"] = {"rc": -3, "log": "", "seconds": 0}

    # Health check: ensure key output files exist
    required = [
        ROOT / "results" / "ml" / "v4_trackA_summary.json",
        ROOT / "results" / "tables" / "v4_bethesda_operating_points.tsv",
        ROOT / "results" / "ml" / "v4_trackC_status.txt",
        ROOT / "results" / "ml" / "v4_trackD_audit_checklist.json",
        ROOT / "reports" / "html" / "pages" / "29_combat_rescue.html",
        ROOT / "reports" / "html" / "pages" / "30_bethesda_clinical.html",
        ROOT / "reports" / "html" / "pages" / "31_prjeb11591_status.html",
        ROOT / "reports" / "html" / "pages" / "32_honest_audit.html",
        ROOT / "reports" / "html" / "pages" / "33_v4_synthesis.html",
        ROOT / "reports" / "v4_synthesis.md",
    ]
    missing = [str(p) for p in required if not p.exists()]
    if missing:
        logger.warning(f"HEALTH: missing {len(missing)} artifact(s):")
        for m in missing:
            logger.warning(f"  - {m}")
    else:
        logger.info("HEALTH: all required artifacts present")

    logger.info(f"results: {results}")
    logger.info("v4 orchestrator done")
    # Orchestrator succeeds even if individual tracks fail, per spec.
    return 0


if __name__ == "__main__":
    sys.exit(main())
