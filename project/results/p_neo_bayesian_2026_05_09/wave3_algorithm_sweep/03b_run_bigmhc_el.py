#!/usr/bin/env python3
"""Score with BigMHC EL (eluted-ligand presentation head) for completeness."""
import subprocess
import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parent
BIGMHC = Path("/data/thca/_repo_offload/wave3_tools/bigmhc")
SRC = BIGMHC / "src"
PYTHON = "/data/thca/_repo_offload/wave3_tools/venv_mhc/bin/python"

inp = ROOT / "bigmhc_input.csv"
out = ROOT / "bigmhc_el_output.csv"

cmd = [
    PYTHON, str(SRC / "predict.py"),
    f"-i={inp}",
    f"-o={out}",
    "-m=el",
    "-a=0",
    "-p=1",
    "-c=1",
    "-d=cpu",
    "-v=1",
]
res = subprocess.run(cmd, capture_output=True, text=True, cwd=str(SRC))
print(res.stdout[-1000:])
print(res.stderr[-1000:])
o = pd.read_csv(out)
print("EL output columns:", list(o.columns))
print(o.head())
