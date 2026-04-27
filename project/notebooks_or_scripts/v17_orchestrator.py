#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path

from v17_common import FIG, LOG, RPT, TAB, VENV_PY, log_line

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/notebooks_or_scripts")

TASKS = {
    "task1": ROOT / "v17_idea1_dark_matter.py",
    "task2": ROOT / "v17_idea2_quad_group.py",
    "task3": ROOT / "v17_idea3_fusion_landscape.py",
    "task4": ROOT / "v17_idea4_dial_audit.py",
    "task5": ROOT / "v17_idea5_trajectory.py",
}


def launch(name: str) -> subprocess.Popen:
    log_line(LOG, f"launch {name}")
    return subprocess.Popen([str(VENV_PY), str(TASKS[name])], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)


def drain(name: str, proc: subprocess.Popen) -> int:
    out = proc.communicate()[0]
    if out:
        with (TAB / f"{name}.log").open("w", encoding="utf-8") as fh:
            fh.write(out)
    log_line(LOG, f"finish {name} rc={proc.returncode}")
    return int(proc.returncode)


def write_summary() -> None:
    lines = ["# v17 Summary", "", "## Quick findings", ""]
    p1 = TAB / "dark_matter_summary.json"
    if p1.exists():
        d = json.loads(p1.read_text())
        lines.append(f"- Dark Matter TCGA cohort: {d['n_samples']} samples, best K={d['best_k']}, cluster counts={d['cluster_counts']}.")
    p2 = TAB / "sample_master_v17_tert.tsv"
    if p2.exists():
        df = __import__("pandas").read_csv(p2, sep="\t")
        lines.append(f"- TERT mutated TCGA tumors: {int((df['tert_status']=='mutated').sum())} / {len(df)}.")
    p3 = TAB / "driver_landscape_v17_summary.tsv"
    if p3.exists():
        df = __import__("pandas").read_csv(p3, sep="\t")
        lines.append(f"- v17 driver landscape top groups: {', '.join((df['driver_anchor_v17'] + '=' + df['n_samples'].astype(str)).head(6).tolist())}.")
    p4 = TAB / "dial_audit_v17.tsv"
    if p4.exists():
        df = __import__("pandas").read_csv(p4, sep="\t")
        if not df.empty:
            lines.append(f"- DIAL audit rows: {len(df)}; mean DIA-AUC={df['DIA_AUC'].mean():.3f}, mean identifiability={df['identifiability'].mean():.3f}.")
    p5 = TAB / "trajectory_pseudotime.tsv"
    if p5.exists():
        df = __import__("pandas").read_csv(p5, sep="\t")
        lines.append(f"- Trajectory pseudotime computed for {len(df)} tumor samples.")
    (RPT / "v17_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    html = f"""<!doctype html><html><head><meta charset='utf-8'><title>v17 summary</title></head>
<body style="font-family:system-ui;max-width:1100px;margin:40px auto;padding:0 16px">
<h1>v17 dashboard stub</h1>
<ul>
  <li><a href="../../results/v17/figs/dark_matter_umap.html">dark matter umap</a></li>
  <li><a href="../../results/v17/figs/quad_group_km_pfi.html">quad group PFI KM</a></li>
  <li><a href="../../results/v17/figs/driver_landscape_oncoplot.html">driver landscape oncoplot</a></li>
  <li><a href="../../results/v17/figs/dial_forest_plot.html">dial forest plot</a></li>
  <li><a href="../../results/v17/figs/trajectory_umap_3d.html">trajectory umap 3d</a></li>
</ul>
<p><a href="v17_summary.md">summary markdown</a></p>
</body></html>"""
    (RPT / "index.html").write_text(html, encoding="utf-8")


def main() -> int:
    log_line(LOG, "v17 orchestrator start")
    procs = {k: launch(k) for k in ("task1", "task2", "task3")}
    rc = {}
    rc["task1"] = drain("task1", procs["task1"])
    rc["task2"] = drain("task2", procs["task2"])
    rc["task3"] = drain("task3", procs["task3"])
    if rc["task1"] == 0:
        p4 = launch("task4")
        rc["task4"] = drain("task4", p4)
    else:
        rc["task4"] = 99
        log_line(LOG, "skip task4 because task1 failed")
    if rc["task3"] == 0:
        p5 = launch("task5")
        rc["task5"] = drain("task5", p5)
    else:
        rc["task5"] = 99
        log_line(LOG, "skip task5 because task3 failed")
    write_summary()
    log_line(LOG, f"v17 orchestrator done {json.dumps(rc)}")
    return 0 if all(v == 0 for k, v in rc.items() if k != "task4" and k != "task5") else 1


if __name__ == "__main__":
    sys.exit(main())
