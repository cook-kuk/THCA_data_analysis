#!/usr/bin/env python3
"""00 — verify required libraries, log environment, ensure data symlink resolves.

Usage:  python3 scripts/00_setup_environment.py
"""
from __future__ import annotations
import importlib
import platform
import sys
from datetime import datetime
from pathlib import Path

REQUIRED = ["pandas", "numpy", "scipy", "sklearn", "matplotlib", "yaml", "requests", "GEOparse"]
OPTIONAL = ["openpyxl"]

ROOT = Path(__file__).resolve().parent.parent
LOG = ROOT / "logs" / "env.txt"


def check(pkg: str) -> tuple[str, str]:
    try:
        mod = importlib.import_module(pkg)
        ver = getattr(mod, "__version__", "unknown")
        return pkg, ver
    except Exception as e:
        return pkg, f"MISSING ({e})"


def main() -> None:
    LOG.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        f"# RAI atlas env audit — {datetime.now().isoformat(timespec='seconds')}",
        f"python:   {sys.version.splitlines()[0]}",
        f"platform: {platform.platform()}",
        f"cwd:      {ROOT}",
        "",
        "## required libraries",
    ]
    missing = []
    for pkg in REQUIRED:
        name, ver = check(pkg)
        line = f"  {name:14s} {ver}"
        lines.append(line)
        print(line)
        if "MISSING" in ver:
            missing.append(name)
    lines.append("")
    lines.append("## optional libraries")
    for pkg in OPTIONAL:
        name, ver = check(pkg)
        lines.append(f"  {name:14s} {ver}")
        print(f"  {name:14s} {ver}")
    lines.append("")
    data_link = ROOT / "data"
    lines.append(f"## data symlink")
    if data_link.is_symlink():
        target = data_link.resolve()
        ok = target.exists() and target.is_dir()
        lines.append(f"  data → {target}  [{'OK' if ok else 'BROKEN'}]")
        print(f"  data → {target}  [{'OK' if ok else 'BROKEN'}]")
    else:
        lines.append("  data is not a symlink (consider running mkdir under /data2/rai_atlas + ln -s)")
        print("  WARN: data is not a symlink")

    LOG.write_text("\n".join(lines) + "\n")
    if missing:
        print(f"\nMISSING REQUIRED PACKAGES: {missing}")
        sys.exit(1)
    print(f"\nenv log → {LOG}")


if __name__ == "__main__":
    main()
