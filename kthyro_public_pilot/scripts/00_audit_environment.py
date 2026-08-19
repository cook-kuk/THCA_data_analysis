#!/usr/bin/env python3
from __future__ import annotations

import argparse
import importlib
import importlib.metadata
import platform
import subprocess
import sys
from pathlib import Path

from pilot_utils import ensure_standard_dirs, pilot_root, setup_logging


MODULES = [
    ("pandas", "pandas"),
    ("numpy", "numpy"),
    ("scipy", "scipy"),
    ("sklearn", "scikit-learn"),
    ("statsmodels", "statsmodels"),
    ("matplotlib", "matplotlib"),
    ("seaborn", "seaborn"),
    ("scanpy", "scanpy"),
    ("anndata", "anndata"),
    ("squidpy", "squidpy"),
    ("gseapy", "gseapy"),
    ("lifelines", "lifelines"),
    ("requests", "requests"),
    ("GEOparse", "GEOparse"),
    ("pyarrow", "pyarrow"),
    ("openpyxl", "openpyxl"),
    ("tqdm", "tqdm"),
    ("yaml", "pyyaml"),
]


def check_module(import_name: str, package_name: str) -> dict:
    try:
        importlib.import_module(import_name)
        try:
            version = importlib.metadata.version(package_name)
        except importlib.metadata.PackageNotFoundError:
            version = "installed_version_unknown"
        return {"package": package_name, "import": import_name, "status": "OK", "version": version}
    except Exception as exc:
        return {"package": package_name, "import": import_name, "status": "MISSING", "version": "", "error": repr(exc)}


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit Python environment for K-Thyro public pilot.")
    parser.add_argument("--install-missing", action="store_true", help="Install missing packages with pip --user.")
    args = parser.parse_args()

    ensure_standard_dirs()
    logger = setup_logging("00_audit_environment")
    rows = [check_module(import_name, package_name) for import_name, package_name in MODULES]
    missing = [r["package"] for r in rows if r["status"] != "OK"]

    out = pilot_root() / "results" / "logs" / "environment_audit.txt"
    with out.open("w") as handle:
        handle.write(f"python_executable\t{sys.executable}\n")
        handle.write(f"python_version\t{sys.version.replace(chr(10), ' ')}\n")
        handle.write(f"platform\t{platform.platform()}\n")
        handle.write("\npackage\timport\tstatus\tversion\terror\n")
        for row in rows:
            handle.write(
                f"{row['package']}\t{row['import']}\t{row['status']}\t{row.get('version','')}\t{row.get('error','')}\n"
            )

    if missing:
        req = pilot_root() / "requirements_public_pilot.txt"
        logger.info("Missing packages: %s", ", ".join(missing))
        logger.info("Requirements file available at %s", req)
        if args.install_missing:
            cmd = [sys.executable, "-m", "pip", "install", "--user", *missing]
            logger.info("Installing missing packages: %s", " ".join(cmd))
            subprocess.run(cmd, check=False)
    else:
        logger.info("All required packages are importable.")
    logger.info("Wrote environment audit to %s", out)


if __name__ == "__main__":
    main()
