"""TESLA manual ingest — runs after a human has dropped the Cell suppl xlsx files
into /data/neoantigen_vaccine_hub/data_raw/tesla/cell_supp_tables/

Reads any *.xlsx in that dir and produces a master-schema-ready table.
Then run scripts/build_master_dataset.py to rebuild the master with TESLA included.
"""
import sys
from pathlib import Path
import pandas as pd

SRC = Path("/data/neoantigen_vaccine_hub/data_raw/tesla/cell_supp_tables")
OUT = Path("/data/neoantigen_vaccine_hub/data_raw/tesla/tesla_master_ready.tsv")


def main():
    files = list(SRC.glob("*.xlsx"))
    files = [f for f in files if f.stat().st_size > 5000]
    if not files:
        print("[TESLA] no xlsx files larger than 5 KB found.")
        print(f"        Drop Cell suppl xlsx into {SRC}")
        print(f"        Then rerun: python3 {Path(__file__).name}")
        sys.exit(1)

    rows = []
    for f in sorted(files):
        try:
            xls = pd.ExcelFile(f)
        except Exception as e:
            print(f"   {f.name}: not a real xlsx ({e})")
            continue
        print(f"   {f.name} sheets: {xls.sheet_names}")
        for sn in xls.sheet_names:
            df = pd.read_excel(f, sheet_name=sn)
            cols_lower = [c.lower() for c in df.columns]
            has_pep = any(k in c for c in cols_lower for k in
                ("peptide","epitope","neopeptide","mutpep","mt_peptide"))
            has_hla = any(k in c for c in cols_lower for k in ("hla","mhc","allele"))
            print(f"     {sn}: {len(df)} rows · cols: {list(df.columns)[:10]}")
            if has_pep and has_hla:
                df["__source_file"] = f.name
                df["__source_sheet"] = sn
                rows.append(df)

    if not rows:
        print("[TESLA] no peptide+HLA tables found in any sheet")
        sys.exit(1)
    big = pd.concat(rows, ignore_index=True, sort=False)
    big.to_csv(OUT, sep="\t", index=False)
    print(f"[TESLA] wrote {OUT} · n={len(big)} · cols={big.shape[1]}")
    print("        Now run: python3 scripts/build_master_dataset.py")


if __name__ == "__main__":
    main()
