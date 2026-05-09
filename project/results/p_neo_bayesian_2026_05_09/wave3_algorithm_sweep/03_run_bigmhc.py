#!/usr/bin/env python3
"""Score ITSNdb scoring set with BigMHC IM (immunogenicity head)."""
import subprocess
import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parent
BIGMHC = Path("/data/thca/_repo_offload/wave3_tools/bigmhc")
SRC = BIGMHC / "src"
PYTHON = "/data/thca/_repo_offload/wave3_tools/venv_mhc/bin/python"

# BigMHC takes input CSV with mhc,pep[,tgt]
# Run on bigmhc_input.csv we already prepared
inp = ROOT / "bigmhc_input.csv"
out = ROOT / "bigmhc_im_output.csv"

# Header count = 1, mhc col = 0, pep col = 1
cmd = [
    PYTHON, str(SRC / "predict.py"),
    f"-i={inp}",
    f"-o={out}",
    "-m=im",
    "-a=0",
    "-p=1",
    "-c=1",
    "-d=cpu",
    "-v=1",
]
print("Running:", " ".join(cmd))
res = subprocess.run(cmd, capture_output=True, text=True, cwd=str(SRC))
print("STDOUT tail:", res.stdout[-1500:])
print("STDERR:", res.stderr[-1500:])

# Inspect output
o = pd.read_csv(out)
print("BigMHC output columns:", list(o.columns))
print(o.head())
print(f"rows: {len(o)}")

# Merge back with scoring set
df = pd.read_csv(ROOT / "scoring_set_itsndb.tsv", sep="\t")
# BigMHC output has same column names (mhc, pep, tgt) plus prediction columns
# Scoring column for IM is typically "BigMHC_IM" — check
o = o.rename(columns={c: f"bigmhc_{c}" if c not in ["mhc", "pep", "tgt"] else c for c in o.columns})
df_m = df.merge(
    o.rename(columns={"mhc": "HLA_norm", "pep": "peptide"})
     .drop(columns=["tgt"], errors="ignore"),
    on=["peptide", "HLA_norm"], how="left"
)
df_m.to_csv(ROOT / "bigmhc_scored.tsv", sep="\t", index=False)
print(f"saved {ROOT / 'bigmhc_scored.tsv'}")
