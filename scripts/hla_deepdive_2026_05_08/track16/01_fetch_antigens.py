#!/usr/bin/env python3
"""
Track 16 Step 1: Fetch UniProt canonical sequences for 14 thyroid self-antigens.
Boundary: thyroid autoimmunity peptidomics — not cancer neoantigen prediction.
"""
from __future__ import annotations
import os
import sys
import time
import urllib.request
import urllib.error
from pathlib import Path

OUTDIR = Path("/home/seungho/personal/THCA_data_analysis/project/results/hla_deepdive_2026_05_08/track16_peptidomics")
FASTA_DIR = OUTDIR / "fasta"
FASTA_DIR.mkdir(parents=True, exist_ok=True)

ANTIGENS = [
    ("TG",      "P01266", "thyroglobulin"),
    ("TPO",     "P07202", "thyroid peroxidase"),
    ("TSHR",    "P16473", "TSH receptor"),
    ("SLC5A5",  "Q92911", "Na/I symporter (NIS)"),
    ("DUOX1",   "Q9NRD9", "dual oxidase 1"),
    ("DUOX2",   "Q9NRD8", "dual oxidase 2"),
    ("DIO1",    "P49895", "iodothyronine deiodinase 1"),
    ("DIO2",    "Q92813", "iodothyronine deiodinase 2"),
    ("IYD",     "Q6PHW0", "iodotyrosine deiodinase"),
    ("PAX8",    "Q06710", "paired box 8"),
    ("FOXE1",   "O00358", "forkhead box E1"),
    ("NKX2-1",  "P43699", "thyroid transcription factor 1"),
    ("HHEX",    "Q03014", "hematopoietically expressed homeobox"),
    ("CALCA",   "P01258", "calcitonin"),
]


def fetch_fasta(accession: str) -> str | None:
    url = f"https://www.uniprot.org/uniprotkb/{accession}.fasta"
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 track16"})
            with urllib.request.urlopen(req, timeout=30) as resp:
                return resp.read().decode("utf-8")
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as e:
            print(f"  attempt {attempt+1} failed for {accession}: {e}", file=sys.stderr)
            time.sleep(2)
    return None


def parse_seq(fasta_text: str) -> str:
    lines = [l.strip() for l in fasta_text.strip().splitlines()]
    return "".join(l for l in lines if not l.startswith(">"))


def main() -> int:
    rows = ["gene\taccession\tdescription\tlength\tfasta_path"]
    combined_path = FASTA_DIR / "all_antigens.fasta"
    combined_handle = open(combined_path, "w")
    for gene, accession, desc in ANTIGENS:
        out_path = FASTA_DIR / f"{gene}_{accession}.fasta"
        if out_path.exists() and out_path.stat().st_size > 50:
            text = out_path.read_text()
            print(f"  cached: {gene} {accession}")
        else:
            print(f"  fetching {gene} {accession} ...")
            text = fetch_fasta(accession)
            if text is None:
                print(f"  FAILED to fetch {gene} {accession}; aborting.", file=sys.stderr)
                return 2
            out_path.write_text(text)
        seq = parse_seq(text)
        rows.append(f"{gene}\t{accession}\t{desc}\t{len(seq)}\t{out_path.relative_to(OUTDIR)}")
        combined_handle.write(text if text.endswith("\n") else text + "\n")
    combined_handle.close()

    summary = OUTDIR / "tables" / "T01_antigen_panel.tsv"
    summary.parent.mkdir(parents=True, exist_ok=True)
    summary.write_text("\n".join(rows) + "\n")
    print(f"wrote {summary}")
    print(f"wrote {combined_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
