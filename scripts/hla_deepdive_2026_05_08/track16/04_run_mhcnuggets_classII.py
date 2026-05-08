#!/usr/bin/env python3
"""
Track 16 Step 4: mhcnuggets-2.4 HLA-II binding for the 14 thyroid self-antigens
× top-10 Korean HLA-DRB1/DPB1 alleles. Uses subprocess-per-allele to avoid TF/Keras
state pollution that causes mid-batch failures.
Boundary: thyroid autoimmunity peptidomics — not cancer neoantigen prediction.
"""
from __future__ import annotations
import os
import sys
import time
import subprocess
from pathlib import Path
import pandas as pd

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUTDIR = ROOT / "project/results/hla_deepdive_2026_05_08/track16_peptidomics"
PEP_PATH = OUTDIR / "tables" / "T03_peptide_library_15mer.tsv"
WORK = OUTDIR / "_classII_work"
WORK.mkdir(exist_ok=True)
PYBIN = str(ROOT / ".venv/bin/python")
WORKER = str(ROOT / "scripts/hla_deepdive_2026_05_08/track16/04b_classII_per_allele.py")

TOP_DRB1 = ["DRB1*15:01", "DRB1*01:01", "DRB1*09:01", "DRB1*08:03", "DRB1*13:02",
            "DRB1*04:05", "DRB1*12:02", "DRB1*07:01", "DRB1*04:06", "DRB1*14:54"]
TOP_DPB1 = ["DPB1*05:01", "DPB1*02:01", "DPB1*04:02", "DPB1*04:01", "DPB1*03:01",
            "DPB1*02:02", "DPB1*13:01", "DPB1*09:01", "DPB1*135:01", "DPB1*14:01"]

ALL_ALLELES = TOP_DRB1 + TOP_DPB1


def safe_filename(allele: str) -> str:
    return allele.replace("*", "_").replace(":", "_")


def main() -> int:
    pep_df = pd.read_csv(PEP_PATH, sep="\t")
    n_pep = len(pep_df)
    print(f"Loaded {n_pep:,} 15-mer peptides across {pep_df.gene.nunique()} antigens.")

    pep_file = str(WORK / "all_peptides.txt")
    with open(pep_file, "w") as f:
        for p in pep_df["peptide"]:
            f.write(p + "\n")

    rows = []
    failed = []
    t0 = time.time()
    for a in ALL_ALLELES:
        ta = time.time()
        out_csv = str(WORK / f"{safe_filename(a)}.csv")
        rank_csv = str(WORK / f"{safe_filename(a)}_ranks.csv")

        if os.path.exists(out_csv) and os.path.getsize(out_csv) > 1000:
            pass  # cached
        else:
            try:
                proc = subprocess.run(
                    [PYBIN, WORKER, a, pep_file, out_csv],
                    capture_output=True, text=True, timeout=600,
                )
                if proc.returncode != 0:
                    print(f"  FAILED {a}: rc={proc.returncode}\n  stderr-tail: {proc.stderr[-300:]}", file=sys.stderr)
                    failed.append(a)
                    continue
            except subprocess.TimeoutExpired:
                print(f"  TIMEOUT {a}", file=sys.stderr)
                failed.append(a)
                continue

        if not os.path.exists(out_csv):
            print(f"  WARNING: {a}: missing {out_csv}", file=sys.stderr)
            failed.append(a)
            continue
        df = pd.read_csv(out_csv)
        if os.path.exists(rank_csv):
            rdf = pd.read_csv(rank_csv)
            df["human_proteome_rank"] = rdf["human_proteome_rank"].values
        else:
            df["human_proteome_rank"] = float("nan")
        if len(df) != n_pep:
            print(f"  WARNING: {a} length mismatch ({len(df)} vs {n_pep})", file=sys.stderr)
            failed.append(a)
            continue

        out_df = pep_df.copy()
        out_df["ic50_nM"] = df["ic50"].values
        out_df["rank_pct"] = df["human_proteome_rank"].values
        out_df["allele"] = a
        rows.append(out_df)
        print(f"  {a}: n={len(out_df):,} preds, dt={time.time()-ta:.1f}s")

    out = pd.concat(rows, ignore_index=True)
    out["binder_class"] = pd.cut(
        out["rank_pct"],
        bins=[-0.001, 2.0, 10.0, 100.0],
        labels=["SB", "WB", "NB"],
    )
    out_path = OUTDIR / "tables" / "T07_HLA_II_binding_long.tsv.gz"
    out.to_csv(out_path, sep="\t", index=False, compression="gzip")
    print(f"wrote {out_path} ({out_path.stat().st_size/1e6:.1f} MB) total rows={len(out):,}")

    summ = (out.groupby(["allele", "binder_class"], observed=True)
              .size().unstack(fill_value=0).reset_index())
    summ.to_csv(OUTDIR / "tables" / "T08_HLA_II_binders_per_allele.tsv", sep="\t", index=False)

    by_ag = (out.groupby(["gene", "allele", "binder_class"], observed=True)
               .size().unstack(fill_value=0).reset_index())
    by_ag.to_csv(OUTDIR / "tables" / "T09_HLA_II_binders_per_antigen_allele.tsv", sep="\t", index=False)
    print(f"elapsed {time.time()-t0:.0f}s; failed alleles: {failed}")
    return 0 if not failed else 0  # don't fail run; report failures


if __name__ == "__main__":
    sys.exit(main())
