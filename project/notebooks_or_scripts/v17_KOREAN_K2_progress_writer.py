#!/usr/bin/env python3
"""Writes /opt/thyroid-dash/project/submission/npj/k2_progress.json every 10 seconds
with current state of the v5 pipeline. Reads from /tmp/k2_v5_262.log + /data/thca/PRJEB11591_quant_se."""
from __future__ import annotations
import json
import re
import shutil
import subprocess
import time
from pathlib import Path

LOG = Path("/tmp/k2_v5_262.log")
QUANT = Path("/data/thca/PRJEB11591_quant_se")
FASTQ = Path("/data/thca/PRJEB11591_fastq")
OUT_JSON = Path("/opt/thyroid-dash/project/submission/npj/k2_progress.json")
TARGET = 262

PIPELINE_PID_HINT = 1617201  # initial PID; we re-detect each cycle


def is_pipeline_running() -> bool:
    try:
        r = subprocess.run(["pgrep", "-f", "v17_KOREAN_K2_v5_pipeline.py"],
                           capture_output=True, text=True, timeout=3)
        return bool(r.stdout.strip())
    except Exception:
        return False


def main():
    t0 = time.time()
    initial_done = sum(1 for d in QUANT.iterdir() if (d / "abundance.tsv").exists()) if QUANT.exists() else 0
    while True:
        try:
            now = time.time()
            running = is_pipeline_running()
            done = sum(1 for d in QUANT.iterdir() if (d / "abundance.tsv").exists()) if QUANT.exists() else 0
            in_flight = []
            if FASTQ.exists():
                for d in FASTQ.iterdir():
                    if not d.is_dir():
                        continue
                    for f in d.iterdir():
                        if f.name.endswith("_1.fastq.gz") and f.stat().st_size > 100_000_000:
                            in_flight.append({"run": d.name, "size_mb": f.stat().st_size // 1_000_000})
            in_flight = in_flight[:8]

            recent = []
            if LOG.exists():
                try:
                    lines = LOG.read_text(encoding="utf-8", errors="replace").splitlines()
                    recent = lines[-20:]
                except Exception:
                    recent = []

            # progress + rate
            elapsed = now - t0
            new_done = max(done - initial_done, 0)
            rate_per_sec = new_done / elapsed if elapsed > 60 and new_done > 0 else 0
            remaining = max(TARGET - done, 0)
            eta_sec = remaining / rate_per_sec if rate_per_sec > 0 else None

            # disk
            try:
                du = shutil.disk_usage("/data")
                disk = {"total_gb": du.total // (1024**3),
                        "used_gb": du.used // (1024**3),
                        "free_gb": du.free // (1024**3),
                        "pct_used": round(du.used / du.total * 100, 1)}
            except Exception:
                disk = None

            payload = {
                "ts": int(now),
                "ts_human": time.strftime("%Y-%m-%d %H:%M:%S"),
                "running": running,
                "target": TARGET,
                "done": done,
                "remaining": remaining,
                "pct": round(done / TARGET * 100, 1),
                "new_done_this_session": new_done,
                "elapsed_sec": int(elapsed),
                "rate_per_min": round(rate_per_sec * 60, 2) if rate_per_sec else None,
                "eta_sec": int(eta_sec) if eta_sec else None,
                "eta_human": (
                    f"~{int(eta_sec/3600)}h {int((eta_sec%3600)/60)}m" if eta_sec and eta_sec > 60
                    else f"~{int(eta_sec/60)}m" if eta_sec else None
                ),
                "in_flight_downloads": in_flight,
                "recent_log_lines": recent,
                "disk": disk,
            }
            OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
            OUT_JSON.write_text(json.dumps(payload, indent=2, ensure_ascii=False))
        except Exception as e:
            try:
                OUT_JSON.write_text(json.dumps({"error": str(e), "ts_human": time.strftime("%Y-%m-%d %H:%M:%S")}))
            except Exception:
                pass
        time.sleep(10)


if __name__ == "__main__":
    main()
