#!/usr/bin/env python3
"""Watch the active RunPod 6VRN MD job and auto-refresh local reports.

This script is deliberately conservative:
- while the run is active, it syncs state/PDB metadata only;
- after completion, it syncs the full trajectory and reruns the complete
  trajectory-dependent analysis stack.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path


REPO = Path(__file__).resolve().parents[3]
OUT = REPO / "project/results/cross_neo_md_audit_2026_05_10"
WATCH_DIR = OUT / "watchers"
LOCAL_POD1 = OUT / "remote_sync/pod1_thca_neo_bayesian_aux"

SSH_KEY = Path.home() / ".runpod/ssh/RunPod-Key-Go"
POD1 = "root@135.84.176.142"
POD1_PORT = "20878"
REMOTE_BASE = "/workspace/openmm_pilot_10ns_package"
REMOTE_RUN = f"{REMOTE_BASE}/prod_10ns_6VRN_1fs300K"


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def run(cmd: str, *, timeout: int | None = None, check: bool = True) -> subprocess.CompletedProcess:
    return subprocess.run(
        cmd,
        shell=True,
        cwd=REPO,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=timeout,
        check=check,
    )


def ssh(remote_cmd: str, *, timeout: int = 60, check: bool = False) -> subprocess.CompletedProcess:
    quoted = remote_cmd.replace("'", "'\"'\"'")
    cmd = (
        f"ssh -i {SSH_KEY} -p {POD1_PORT} -o StrictHostKeyChecking=no "
        f"-o ConnectTimeout=10 {POD1} '{quoted}'"
    )
    return run(cmd, timeout=timeout, check=check)


def remote_status() -> dict:
    cmd = f"""
set -e
cd {REMOTE_RUN}
tail -n 1 state.tsv || true
echo __FILES__
ls -l state.tsv trajectory.dcd final.pdb final.chk 2>/dev/null || true
echo __PROCS__
pgrep -af '^python run_openmm_pilot.py' || true
pgrep -af 'submit_hmtevvrhc_replicate_screens.sh' || true
"""
    cp = ssh(cmd, timeout=45, check=False)
    text = cp.stdout
    lines = text.splitlines()
    file_lines: list[str] = []
    proc_lines: list[str] = []
    section = "state"
    for line in lines[1:]:
        if line == "__FILES__":
            section = "files"
            continue
        if line == "__PROCS__":
            section = "procs"
            continue
        if section == "files":
            file_lines.append(line)
        elif section == "procs":
            proc_lines.append(line)
    status = {
        "checked_at_utc": now(),
        "ssh_returncode": cp.returncode,
        "remote_text": text,
        "running": any("python run_openmm_pilot.py" in line for line in proc_lines),
        "replicate_watcher_running": any("submit_hmtevvrhc_replicate_screens.sh" in line for line in proc_lines),
        "has_final_pdb": any(line.rstrip().endswith(" final.pdb") for line in file_lines),
        "has_final_chk": any(line.rstrip().endswith(" final.chk") for line in file_lines),
        "trajectory_bytes_remote": None,
        "time_ps": None,
        "temperature_k": None,
        "speed_ns_per_day": None,
        "step": None,
    }
    first = text.splitlines()[0] if text.splitlines() else ""
    parts = first.split("\t")
    if len(parts) >= 5:
        try:
            status["step"] = int(float(parts[0]))
            status["time_ps"] = float(parts[1])
            status["temperature_k"] = float(parts[3])
            status["speed_ns_per_day"] = float(parts[4])
        except Exception:
            pass
    for line in file_lines:
        if line.rstrip().endswith(" trajectory.dcd"):
            fields = line.split()
            if len(fields) >= 5:
                try:
                    status["trajectory_bytes_remote"] = int(fields[4])
                except Exception:
                    pass
    status["complete_by_state"] = bool(status["time_ps"] is not None and status["time_ps"] >= 9999.0)
    status["complete"] = (not status["running"]) and (status["complete_by_state"] or status["has_final_pdb"])
    return status


def sync_state_only() -> None:
    LOCAL_POD1.mkdir(parents=True, exist_ok=True)
    cmd = (
        f"ssh -i {SSH_KEY} -p {POD1_PORT} -o StrictHostKeyChecking=no -o ConnectTimeout=10 {POD1} "
        f"'cd {REMOTE_BASE} && tar --exclude=trajectory.dcd --exclude=final.chk -czf - prod_10ns_6VRN_1fs300K' "
        f"| tar -xzf - -C {LOCAL_POD1}"
    )
    run(cmd, timeout=120, check=True)


def sync_full_run() -> None:
    target = LOCAL_POD1 / "prod_10ns_6VRN_1fs300K"
    target.mkdir(parents=True, exist_ok=True)
    if shutil.which("rsync"):
        cmd = (
            f"rsync -avP -e 'ssh -i {SSH_KEY} -p {POD1_PORT} -o StrictHostKeyChecking=no' "
            f"{POD1}:{REMOTE_RUN}/ {target}/"
        )
        run(cmd, timeout=3600, check=True)
    else:
        cmd = (
            f"ssh -i {SSH_KEY} -p {POD1_PORT} -o StrictHostKeyChecking=no {POD1} "
            f"'cd {REMOTE_BASE} && tar -czf - prod_10ns_6VRN_1fs300K' "
            f"| tar -xzf - -C {LOCAL_POD1}"
        )
        run(cmd, timeout=3600, check=True)


def refresh_light() -> None:
    for script in [
        "00_parse_openmm_status.py",
        "11_generate_md_extra_visuals.py",
        "08_integrate_md_with_cross_neo.py",
        "12_write_md_decision_report.py",
        "10_build_md_visual_dossier.py",
    ]:
        run(f"python project/scripts/cross_neo_md/{script}", timeout=240, check=True)


def refresh_full() -> None:
    for script in [
        "00_parse_openmm_status.py",
        "01_annotate_complex.py",
        "02_md_qc.py",
        "03_analyze_pmhc_contacts.py",
        "04_analyze_tcr_contacts.py",
        "05_counterfactual_md_analysis.py",
        "06_replicate_consistency.py",
        "07_md_evidence_score.py",
        "08_integrate_md_with_cross_neo.py",
        "09_generate_md_figures.py",
        "11_generate_md_extra_visuals.py",
        "12_write_md_decision_report.py",
        "10_build_md_visual_dossier.py",
    ]:
        run(f"python project/scripts/cross_neo_md/{script}", timeout=1800, check=True)


def write_status(status: dict) -> None:
    WATCH_DIR.mkdir(parents=True, exist_ok=True)
    (WATCH_DIR / "hmtevvrhc_autowatch_latest.json").write_text(json.dumps(status, indent=2))
    line = json.dumps(status, sort_keys=True)
    with (WATCH_DIR / "hmtevvrhc_autowatch_events.jsonl").open("a") as handle:
        handle.write(line + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--interval-seconds", type=int, default=300)
    parser.add_argument("--max-hours", type=float, default=12.0)
    parser.add_argument("--once", action="store_true")
    args = parser.parse_args()

    WATCH_DIR.mkdir(parents=True, exist_ok=True)
    deadline = time.time() + args.max_hours * 3600
    completed_once = False
    while True:
        status = remote_status()
        write_status(status)
        if status["ssh_returncode"] != 0:
            pass
        elif status["complete"]:
            sync_full_run()
            refresh_full()
            completed_once = True
            status["post_completion_sync_done_at_utc"] = now()
            write_status(status)
        else:
            sync_state_only()
            refresh_light()
            status["state_sync_done_at_utc"] = now()
            write_status(status)

        if args.once or completed_once or time.time() >= deadline:
            break
        time.sleep(max(30, args.interval_seconds))


if __name__ == "__main__":
    main()
