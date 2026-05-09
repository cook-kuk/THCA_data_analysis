#!/usr/bin/env python3
"""Run CLEAN-NeoBench + BAR-Neo end to end."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=".", help="Repository root")
    parser.add_argument("--output-root", required=True, help="Output directory")
    return parser.parse_args()


def run_step(script_dir: Path, script_name: str, repo_root: Path, output_root: Path) -> None:
    cmd = [
        sys.executable,
        str(script_dir / script_name),
        "--repo-root",
        str(repo_root),
        "--output-root",
        str(output_root),
    ]
    print("[clean-neobench-runner]", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def main() -> None:
    args = parse_args()
    repo_root = Path(args.repo_root).resolve()
    output_root = Path(args.output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    script_dir = Path(__file__).resolve().parent

    manifest_path = output_root / "run_manifest.json"
    manifest_path.write_text(
        json.dumps(
            {
                "pipeline": "CLEAN-NeoBench + BAR-Neo",
                "repo_root": str(repo_root),
                "output_root": str(output_root),
                "scripts": [],
                "warnings": [],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )

    steps = [
        "build_clean_neobench_master.py",
        "build_overlap_flags.py",
        "collect_method_scores.py",
        "build_public_tool_overlap_audit.py",
        "evaluate_clean_neobench.py",
        "train_barneo.py",
        "barneo_bma_agent.py",
        "analyze_failure_modes.py",
        "apply_failure_aware_reliability.py",
        "write_clean_neobench_reports.py",
    ]
    for step in steps:
        run_step(script_dir, step, repo_root, output_root)
        try:
            manifest = json.loads(manifest_path.read_text())
        except Exception:
            manifest = {}
        manifest.setdefault("scripts", [])
        manifest["scripts"].append(step)
        manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(f"[clean-neobench-runner] complete output_root={output_root}", flush=True)


if __name__ == "__main__":
    main()
